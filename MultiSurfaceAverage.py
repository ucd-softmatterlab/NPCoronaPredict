#! /usr/bin/python

import os
import numpy as np
from sys import argv
import matplotlib.pyplot as plt
import scipy.spatial as scispat
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
        elif abspath[-4:] == ".uam":
            pdbs.append(abspath)
    return pdbs

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


def averageEnergy( surfaceSet, uaMapFile, pdbFile="", weightArea = False, areaEpsilon = 0.1):
    ''' Calculate simple and canonical averages for a set of protein binding states.  '''
    allStates = []
    gotAreaWeights = False
    if weightArea == True:
        try:
            pdbStructure = getAtomCoords(pdbFile)*0.1
        except:
            print("Failed to find PDB file for current target, attempted ", pdbFile)
            weightArea = False
        areaWeightList = []
    for surface in surfaceSet:
        data = np.genfromtxt( surface[0]+"/"+uaMapFile )
        weights = np.sin( ( data[:,1]+2.5) * np.pi / 180.0 ) * surface[1] 
        if weightArea == True:
            if gotAreaWeights == False:
                areaWeightList = []
                for dataLine in data:
                    crossSectionalArea =  getPDBAreaXY( pointsToBeads(  rotatePDB(pdbStructure, dataLine[0]*np.pi/180, dataLine[1]*np.pi/180)[:,(0,1)]     ))
                    areaWeightList.append(crossSectionalArea + areaEpsilon)
                gotAreaWeights = True
                areaWeights = 1.0/np.array(areaWeightList)
            weights = weights * areaWeights
        allStates.append( np.stack( [data[:,2] , weights ] ,axis=-1))
    allStateArr = np.concatenate( allStates,axis=0) 
    simpleAverage = np.sum( allStateArr[:,0] * allStateArr[:,1] )/np.sum(  allStateArr[:,1] )
    boltzAverage = np.sum( allStateArr[:,0] * allStateArr[:,1] * np.exp(-1.0 * allStateArr[:,0]) )/np.sum(  allStateArr[:,1] * np.exp(-1.0 * allStateArr[:,0]))
    #freeEAverage = -1.0 * np.log(  np.sum(  allStateArr[:,1] *  np.exp(-1.0 * allStateArr[:,0])    )/np.sum(allStateArr[:,1] )   )
    return (simpleAverage, boltzAverage, np.amin( allStateArr[:,0]))
         

#Define the set of surfaces to be averaged within. The first entry is the results folder containing UA output and the second is the weight (e.g. Wulff weight) for the surface. These are normalised such that they sum to unity in case some surfaces are omitted. 
goldSurfaceSet = [
["results_testproject-anatase-oct/np1R_5_ZP_0", 1.0],
["results_testproject-gold-oct/np1R_5_ZP_0", 0.001]
]

allSurfaceSets = [goldSurfaceSet]

#For each set of surfaces, find all UAM files common to all surfaces and average over these.
for surfaceSet in allSurfaceSets:
    totalSurfaceWeight = 0
    for surface in surfaceSet:
        totalSurfaceWeight += surface[1]
    for i in range(len(surfaceSet)):
        surfaceSet[i][1] = surfaceSet[i][1] / totalSurfaceWeight   
    print(surfaceSet)
    allFoundTargets = []
    for surface in surfaceSet:
        targetCandidates = []
        targetCandidatePaths = SearchDirectory(surface[0])
        for target in targetCandidatePaths:
            uamFileName = target.split("/")[-1]
            targetCandidates.append(uamFileName)
        allFoundTargets.append(targetCandidates)
    uniqueTargets = list( set.intersection(*map(set,list(allFoundTargets))) )
    uniqueTargets.sort()
    for target in uniqueTargets:
        pdbNameGuess = "all_proteins/"+target.split("-")[0]+".pdb"
        simple,boltz,maxEnergy = averageEnergy(surfaceSet, target,pdbNameGuess,weightArea=True)
        print(target, simple, boltz,maxEnergy)
# Main

