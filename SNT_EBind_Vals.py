#! /usr/bin/python

import os
import numpy as np
from sys import argv
import matplotlib.pyplot as plt
import scipy.spatial as scispat
from ortools.algorithms import pywrapknapsack_solver
import math
import scipy.optimize as scopt
import random
import numpy.linalg as npla

def SearchDirectory(path):
    files = os.listdir(path)
    pdbs  = []
    for handle in files:
        abspath = os.path.join(path, handle)
        if os.path.isdir(abspath):
            pdbs += SearchDirectory(abspath)
        elif abspath[-4:] == ".map":
            pdbs.append(abspath)
    return pdbs


def CalculateBoltz(filename):
    data     = np.genfromtxt(filename)
    theta    = data[:,1]
    energy   = data[:,2]
    sinTheta = np.sin(theta * np.pi / 180.0)
    boltz    = np.sum(energy * sinTheta * np.exp(-1.0 * energy)) / np.sum(sinTheta * np.exp(-1.0 * energy))
    return boltz

def MaxBoltz(filename):
    data = np.genfromtxt(filename)
    maxLoc = np.argmin(data[:,2])
    return (data[maxLoc])

def MeanBoltzAtAngle(filename,angle):
    data = np.genfromtxt(filename)

    uprightSet =   np.logical_and(  data[:,1] > angle-5 , data[:,1] < angle+5 )  
    theta    = data[ uprightSet,1]
    energy   = data[ uprightSet,2]
    sinTheta = 1#np.sin(theta * np.pi / 180.0)
    boltz    = np.sum(energy * sinTheta * np.exp(-1.0 * energy)) / np.sum(sinTheta * np.exp(-1.0 * energy))
    return boltz

def getTarget(data, size, protein,shape,zp):
    return data[ np.logical_and(  np.logical_and ( np.logical_and(data[:,4] == size , data[:,0]==protein), data[:,3]==shape )    ,data[:,5] == zp)  ]



def getAtomCoords(filename):
    fileIn = open(filename,"r")
    coordList = []
    for line in fileIn:
        lineData = line.split()
        if lineData[0] == "ATOM":
            coordList.append([ float(lineData[6]) ,  float(lineData[7]) , float(lineData[8])])
    fileIn.close()
 
    return np.array(coordList)
 
def rotatePDB(coords,phiVal,thetaVal):
    phiRotated = -1.0 * phiVal
    thetaRotated = np.pi - thetaVal
    rotCoords = np.copy(coords)
    rotCoords[:,0] = coords[:,0] * np.cos(phiRotated) - coords[:,1] * np.sin(phiRotated)
    rotCoords[:,1] = coords[:,0] * np.sin(phiRotated) + coords[:,1] * np.cos(phiRotated)
    finalCoords = np.copy(rotCoords)
    finalCoords[:,0] = rotCoords[:,0] * np.cos(thetaRotated) + rotCoords[:,2] * np.sin(thetaRotated) 
    finalCoords[:,2] = -1.0 * rotCoords[:,0] * np.sin(thetaRotated) + rotCoords[:,2] * np.cos(thetaRotated)
    return finalCoords

def getPDBAreaXY(coords):
    projectedConvexHull = scispat.ConvexHull( coords[:,(0,1) ] )
    return projectedConvexHull.volume

def getAreaHeight(pointCoords):
    coords=pointsToBeads(pointCoords)
    projectedConvexHull = scispat.ConvexHull( coords[:,(0,1) ] )
    heightAboveSurface = pointCoords[:,2] - np.mean(pointCoords[:,2]) #Distance from the centre (average z location) to the 
    return (projectedConvexHull.volume, heightAboveSurface)

def CalculateBoltzArea(filename,pdbFile):
    data     = np.genfromtxt(filename)
    mfptData = np.genfromtxt(filename[:-4]+"_mfpt.map")
    diffusionCoeff = 1.0e-10
    mfptVals = mfptData[:,2] * 1.0e-18 / diffusionCoeff
    koffVals = 1.0/mfptVals
    rawCoords =  getAtomCoords(pdbFile)*0.1 # convert to nm
    theta    = data[:,1]
    energy   = data[:,2]
    sinTheta = np.sin(theta * np.pi / 180.0)
    boltz    = np.sum(energy * sinTheta * np.exp(-1.0 * energy)) / np.sum(sinTheta * np.exp(-1.0 * energy))
    weightedArea = 0
    energyByArea = 0
    AreaResList =[]
    m, c = npla.lstsq( np.vstack([energy,np.ones(len(energy))]).T , np.log(koffVals),rcond=None)[0]
    for dataLine in data:
        (crossSectionalArea,heightAboveSurface) =  getAreaHeight(rotatePDB(rawCoords, dataLine[0]*np.pi/180, dataLine[1]*np.pi/180))
        weightedArea =  weightedArea + crossSectionalArea * np.sin( dataLine[1]*np.pi/180) * np.exp(-1.0 * dataLine[2] )
        energyByArea = energyByArea + ( dataLine[2]/crossSectionalArea * np.sin( dataLine[1]*np.pi/180) * np.exp(-1.0 * dataLine[2] )) #if we want <E/A> we need to do the averaging all in one go to avoid getting <E>/<A>
        AreaResList.append([crossSectionalArea,dataLine[2],heightAboveSurface])
    return (boltz,  weightedArea / np.sum(sinTheta * np.exp(-1.0 * energy))   , energyByArea/np.sum(sinTheta * np.exp(-1.0 * energy)) ,np.array(AreaResList),np.exp(c),np.exp(c)* np.exp(boltz)) 



#Returns an array of areas and binding energies
def CalculateEnergyAreaList(filename, pdbFile):
    data     = np.genfromtxt(filename)
    rawCoords =  getAtomCoords(pdbFile)*0.1 # convert to nm
    theta    = data[:,1]
    energy   = data[:,2]
    sinTheta = np.sin(theta * np.pi / 180.0)
    boltz    = np.sum(energy * sinTheta * np.exp(-1.0 * energy)) / np.sum(sinTheta * np.exp(-1.0 * energy))
    weightedArea = 0
    AreaResList = []
    for dataLine in data:
        crossSectionalArea =  getPDBAreaXY(pointsToBeads(rotatePDB(rawCoords, dataLine[0]*np.pi/180, dataLine[1]*np.pi/180)))
        AreaResList.append([crossSectionalArea,dataLine[2]])
    return np.array(AreaResList)




def knapsackGreedy(capacity, weights, values):
    valuePerWeight = values/weights
    dataArray = (np.vstack( (weights, values, valuePerWeight))).T
    sortedData = np.flip(dataArray[ dataArray[:,2].argsort()],0)
    usedCapacity = 0
    totalValue = 0
    for dataLine in sortedData:
        weight = dataLine[0]
        value = dataLine[1]
        remainingCapacity = capacity - usedCapacity 
        numObjects = int(math.floor( remainingCapacity/weight ))
        usedCapacity = usedCapacity + numObjects*weight
        totalValue = totalValue + numObjects*value
    return totalValue
       
def geneticKnapsackObjective(quants, weights,vals,capacity):
    integerQuants = np.round(quants)
    valArray = integerQuants*vals
    weightArray = integerQuants*weights
    totalValue = np.sum(valArray)
    totalWeight = np.sum(weightArray)
    overweightPenalty = 0
    if totalWeight > capacity:
        overweightPenalty = 500 + (totalWeight-capacity) - totalValue
    return (totalValue + overweightPenalty)

def geneticKnapsack(weights, vals, capacity):
    maxNumber = 100
    numberOrientations = len(weights)
    boundsList = [] 
    maxNum = 0
    for i in range(numberOrientations): 
        capacityofWeight = int(round(capacity/weights[i]))
        boundsList.append( (0,capacityofWeight) )
        if maxNum < capacityofWeight:
            maxNum = capacityofWeight

    #The genetic algorithm here computes the minimum of the objective function, i.e. the lowest energy possible.
    #the default init procedure doesn't do so well with this problem as it scatters points with equal density over all possible configurations, almost all of which are too high energy.
    #we therefore generate a more reasonable initial set as follows
    initMemberList =[]
    for j in range(30* numberOrientations):
       memberArray = np.zeros(numberOrientations)
       for i in range(int(round(maxNum*0.5))):
            randomIndex = np.random.randint(0, len(memberArray))
            memberArray[randomIndex] = memberArray[randomIndex]+1
       initMemberList.append(memberArray)
    initialPop = np.array(initMemberList)
    solver = scopt.differential_evolution( geneticKnapsackObjective, boundsList, args=(weights,vals,capacity) ,init=initialPop,recombination = 0.8)
    integerSol = np.round(solver.x)
    return np.sum(integerSol * vals)

def pointsToBeads(coords):
    if(len(coords)>2):
        projectedConvexHull = scispat.ConvexHull( coords[:,(0,1) ] )
        edgeCoords = coords[projectedConvexHull.vertices ]
    else:
        edgeCoords = coords[:,(0,1)]
    beadRadius = 0.5
    numPoints = 8
    thetaSet = np.arange(0,numPoints)/(1.0*numPoints)  * np.pi*2
    circlex = beadRadius*np.cos(thetaSet)
    circley = beadRadius*np.sin(thetaSet)
    newCoordList = []
    for baseCoords in edgeCoords[:,(0,1)]:
        newCoordList.append(  (np.stack((circlex + baseCoords[0],circley + baseCoords[1])) ).T  )
    return np.reshape(np.array(newCoordList),(-1,2))


def AppendArea(filename,pdbFile):
    data     = np.genfromtxt(filename)
    rawCoords =  getAtomCoords(pdbFile)*0.1 # convert to nm
    theta    = data[:,1]
    energy   = data[:,2]
    sinTheta = np.sin(theta * np.pi / 180.0)
    boltz    = np.sum(energy * sinTheta * np.exp(-1.0 * energy)) / np.sum(sinTheta * np.exp(-1.0 * energy))
    weightedArea = 0
    AreaResList = []
    for dataLine in data:
        crossSectionalArea =  getPDBAreaXY( pointsToBeads(  rotatePDB(rawCoords, dataLine[0]*np.pi/180, dataLine[1]*np.pi/180)[:,(0,1)]     ))
        AreaResList.append([dataLine[0],dataLine[1],dataLine[2],crossSectionalArea])
    return np.array(AreaResList)




# Main

resList = []
'''
spaWeight = 0.1
spdWeight = 0.05
dppcWeight = 1-(spaWeight+spdWeight)

rut001FaceEnergy = 1.28
rut100FaceEnergy = 1.68
rut110FaceEnergy = 2.19


 
ana101Weight = 0.96
ana100Weight = 0
ana001Weight = 0.04
rut001Weight = 0
rut100Weight = 0.44
rut110Weight = 0.56
'''

folderNamesSphere = ["results_fe2o3_sphere"]
folderNamesCylinder = []
#folderNamesSphere = ["results_rutile_alltargets_sphere","results_anatase_alltargets_sphere"]
#folderNamesCylinder = ["results_anatase_alltargets_cylinder", "results_swcnt_alltargets","results_swcntOH4_alltargets","results_swcntOH14_alltargets","results_swcntCOOH3_alltargets","results_swcntCOOH30_alltargets","results_swcntNH2-2_alltargets", "results_mwcnt_alltargets","results_mwcntOH4_alltargets","results_mwcntOH14_alltargets","results_mwcntCOOH3_alltargets","results_mwcntCOOH30_alltargets","results_mwcntNH2-2_alltargets"]
folderNamesCube = []
#folderNames = ["results_spa12d_to2kbt_anatase_001_sphere"]


resList = []
for folderName in folderNamesSphere:
#from the folder name we extract the details about what type of NP we're looking at
    folderNameTerms = folderName.split("_")
#    print folderNameTerms[3]
    pdbs = (SearchDirectory(folderName)) 
    pdbs.sort()
    for pdb in pdbs:
        name = pdb.split('/')[-1].split('.')[0]
        if len(pdb) > 9 and pdb[-9:] == "_mfpt.map":
            continue
        proteinName = name.split("_")[0]
        if proteinName == "spa2":
            continue
        boltzEnergy = CalculateBoltz(pdb)
        nameTerms = name.split("_")
        print folderNameTerms[1], "sphere", nameTerms[1], float(nameTerms[2])/1000.0, nameTerms[0],  boltzEnergy
        resList.append([ folderNameTerms[1], nameTerms[1], nameTerms[2], nameTerms[0], boltzEnergy])


for folderName in folderNamesCylinder:
#from the folder name we extract the details about what type of NP we're looking at
    folderNameTerms = folderName.split("_")
#    print folderNameTerms[3]
    pdbs = (SearchDirectory(folderName))
    pdbs.sort()
    for pdb in pdbs:
        name = pdb.split('/')[-1].split('.')[0]
        if len(pdb) > 9 and pdb[-9:] == "_mfpt.map":
            continue
        proteinName = name.split("_")[0]
        if proteinName == "spa2":
            continue
        boltzEnergy = CalculateBoltz(pdb)
        nameTerms = name.split("_")
        print folderNameTerms[1], "cylinder", nameTerms[1], float(nameTerms[2])/1000.0, nameTerms[0],  boltzEnergy
        resList.append([ folderNameTerms[1], nameTerms[1], nameTerms[2], nameTerms[0], boltzEnergy])



for folderName in folderNamesCube:
#from the folder name we extract the details about what type of NP we're looking at
    folderNameTerms = folderName.split("_")
#    print folderNameTerms[3]
    pdbs = (SearchDirectory(folderName))
    pdbs.sort()
    for pdb in pdbs:
        name = pdb.split('/')[-1].split('.')[0]
        if len(pdb) > 9 and pdb[-9:] == "_mfpt.map":
            continue
        proteinName = name.split("_")[0]
        if proteinName == "spa2":
            continue
        boltzEnergy = CalculateBoltz(pdb)
        nameTerms = name.split("_")
        print folderNameTerms[1], "cube", nameTerms[1], float(nameTerms[2])/1000.0, nameTerms[0],  boltzEnergy
        resList.append([ folderNameTerms[1], nameTerms[1], nameTerms[2], nameTerms[0], boltzEnergy])


'''
#return data[ np.logical_and(  np.logical_and ( np.logical_and(data[:,4] == size , data[:,0]==protein), data[:,3]==shape )    ,data[:,5] == zp)  ]
resultArray =np.array(resList)
ana001 =  resultArray[   np.logical_and(resultArray[:,2]=="anatase"  ,  resultArray[:,1]=="001"  )  ]
ana100 =  resultArray[ np.logical_and( resultArray[:,2]=="anatase" ,  resultArray[:,1]=="100" )  ]
ana101 =  resultArray[ np.logical_and( resultArray[:,2]=="anatase" ,  resultArray[:,1]=="101" )  ]

for npType in ["sphere","cylinder","cube"]:
    sizeList = ana001[np.logical_and(np.logical_and(  ana001[:,0 ]=="3dbz" , ana001[:,3]==npType), ana001[:,5]=="-7468"),4]
    zpList = ana001[np.logical_and(np.logical_and(  ana001[:,0 ]=="3dbz" , ana001[:,3]==npType), ana001[:,4]=="3"),5]
    for proteinString in ["3dbz","spa1","DPPCHead","DPPCTail", "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HID", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "EST","CHL","PHO"]:
        for sizeString in sizeList:
            for zp in zpList:
                ana001Data = getTarget(ana001, sizeString, proteinString, npType,zp)[0]
                ana100Data = getTarget(ana100, sizeString, proteinString, npType,zp)[0]
                ana101Data = getTarget(ana101, sizeString, proteinString, npType,zp)[0]
                ana001Energy = ana001Data[6].astype(np.float)
                ana100Energy = ana100Data[6].astype(np.float)
                ana101Energy = ana101Data[6].astype(np.float)
                ana001EByContact = ana001Data[8].astype(np.float)
                ana100EByContact = ana100Data[8].astype(np.float)
                ana101EByContact = ana101Data[8].astype(np.float)
                ana001EByContactKS = ana001Data[9].astype(np.float)
                ana100EByContactKS = ana100Data[9].astype(np.float)
                ana101EByContactKS = ana101Data[9].astype(np.float)
                ana001KOn = ana001Data[10].astype(np.float)
                ana100KOn = ana100Data[10].astype(np.float)
                ana101KOn = ana101Data[10].astype(np.float)
                ana001KOff = ana001Data[11].astype(np.float)
                ana100KOff = ana100Data[11].astype(np.float)
                ana101KOff = ana101Data[11].astype(np.float)
                print proteinString, "anatase", npType, sizeString,zp, ana001Energy*ana001Weight + ana100Energy*ana100Weight + ana101Energy*ana101Weight, ana001EByContact*ana001Weight + ana100EByContact*ana100Weight + ana101EByContact*ana101Weight, ana001EByContactKS*ana001Weight + ana100EByContactKS*ana100Weight + ana101EByContactKS*ana101Weight, ana001KOn*ana001Weight + ana100KOn*ana100Weight + ana101KOn*ana101Weight, ana001KOff*ana001Weight + ana100KOff*ana100Weight + ana101KOff*ana101Weight


rut001 =  resultArray[ np.logical_and( resultArray[:,2]=="rutile" ,  resultArray[:,1]=="001" )  ]
rut100 =  resultArray[ np.logical_and( resultArray[:,2]=="rutile" ,  resultArray[:,1]=="100" )  ]
rut110 =  resultArray[ np.logical_and( resultArray[:,2]=="rutile" ,  resultArray[:,1]=="110" )  ]
 
for npType in ["sphere","cylinder","cube"]:
    sizeList = rut001[np.logical_and(  np.logical_and(  rut001[:,0 ]=="3dbz" , rut001[:,3]==npType)  , rut001[:,5]=="-7468") ,4]
    zpList = rut001[np.logical_and(  np.logical_and(  rut001[:,0 ]=="3dbz" , rut001[:,3]==npType)  , rut001[:,4]=="3") ,5]
    for proteinString in ["3dbz","spa1","DPPCHead","DPPCTail", "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HID", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "EST","CHL","PHO"]:
        for sizeString in sizeList:
            for zp in zpList:
                rut001Data = getTarget(rut001, sizeString, proteinString, npType,zp)[0]
                rut100Data = getTarget(rut100, sizeString, proteinString, npType,zp)[0]
                rut110Data = getTarget(rut110, sizeString, proteinString, npType,zp)[0]
                rut001Energy = rut001Data[6].astype(np.float)
                rut100Energy = rut100Data[6].astype(np.float)
                rut110Energy = rut110Data[6].astype(np.float)
                rut001EByContact = rut001Data[8].astype(np.float)
                rut100EByContact = rut100Data[8].astype(np.float)
                rut110EByContact = rut110Data[8].astype(np.float)
                rut001EByContactKS = rut001Data[9].astype(np.float)
                rut100EByContactKS = rut100Data[9].astype(np.float)
                rut110EByContactKS = rut110Data[9].astype(np.float)
                rut001KOn = rut001Data[10].astype(np.float)
                rut100KOn = rut100Data[10].astype(np.float)
                rut110KOn = rut110Data[10].astype(np.float)
                rut001KOff = rut001Data[11].astype(np.float)
                rut100KOff = rut100Data[11].astype(np.float)
                rut110KOff = rut110Data[11].astype(np.float)
                print proteinString, "rutile", npType, sizeString, zp, rut001Energy*rut001Weight + rut100Energy*rut100Weight + rut110Energy*rut110Weight, rut001EByContact*rut001Weight + rut100EByContact*rut100Weight +rut110EByContact*rut110Weight, rut001EByContactKS*rut001Weight + rut100EByContactKS*rut100Weight +rut110EByContactKS*rut110Weight, rut001KOn*rut001Weight + rut100KOn*rut100Weight +rut110KOn*rut110Weight, rut001KOff*rut001Weight + rut100KOff*rut100Weight +rut110KOff*rut110Weight
'''

