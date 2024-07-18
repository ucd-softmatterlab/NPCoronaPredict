'''
Given a list of protein (or other molecule) concentrations and binding energies runs an MC model to determine the coverage in the steady-state.

The state is stored internally as a list of lists:

[newProteinID,newPhi,newC2 ,collidingNP ]
where proteinID associates an instance of a protein with its dataset entry, phi,C2 are the co-ordinates (C2=theta for sphere, z for cylinder) and NP is the ID number of the NP to which the protein is bound

'''
import os
import numpy as np
import scipy as sp
import scipy.special as scspec
import random
import argparse
import matplotlib.pyplot as plt

def analyticSol(t,conc0,kon,koff,numBindingSites):
    return conc0*kon*numBindingSites/(koff+ conc0*kon) - conc0*kon*numBindingSites * np.exp(-t *(koff+conc0 * kon))/(koff+conc0*kon)

def analyticSolNoAreaLimit(t,conc0,kon,koff,numBindingSites):
    return conc0*kon*numBindingSites/(koff) - conc0*kon*  np.exp(-t *(koff ))/(koff )

#check to see if protein in position entryID is overlapping with any other proteins
def collisionDetect(state, entryID):
    if(len(state))<2:
        return 0
    entrySize = proteinData[  int(state[int( entryID),0]),5]
    for i in range(len(state)):
        if i== entryID:
            continue
        proteinID = int(state[int(i),0])
        dist = np.sqrt( (state[int(i),1] -  state[entryID,1] )**2 + (state[int(i),2] -  state[entryID,2] )**2 + (state[int(i),3] -  state[entryID,3] )**2)
        if dist < proteinData[ proteinID,5] + entrySize:
            return 1 
    return 0
    

#Routines for testing if the proposed protein would overlap with any existing ones

enableMFSpaceTest = 1 #enables the second test in the meanfield model to see if this would bring the surface coverage above 1

#This function returns a variable denoting if any collision was found and a list of all proteins which would overlap with the new one
def adsorbCollisionDetect(state,newType,newC1,newC2):
    collisionDetected = 0
    detectedOverlaps = []
    if meanFieldApprox == 1:
        if np.random.random() < 1 - surfaceCoverage :
            if enableMFSpaceTest == 1:
                newCoverage = surfaceCoverage + 1.0/proteinBindingSites[newType]
                if newCoverage < 1:
                    collisionDetected = 0
                else:
                    collisionDetected = 1
            else:
                collisionDetected = 0
        else:
            return 1
    if npShape == 1:
        collisionDetected, detectedOverlaps =  SphereCollisionDetect(state,newType,newC1,newC2)
    elif npShape == 2:
        if hasPBC == True:
            collisionDetected, detectedOverlaps = CylinderCollisionDetect(state,newType,newC1,newC2) 
            collisionDetected2, detectedOverlaps2 = CylinderCollisionDetect(state,newType,newC1,newC2,2*cylinderHalfLength) 
            collisionDetected3, detectedOverlaps3 =  CylinderCollisionDetect(state,newType,newC1,newC2,-2*cylinderHalfLength)
            collisionDetected = max(collisionDetected , collisionDetected2 , collisionDetected3)
            detectedOverlaps  = list(set( detectedOverlaps  + detectedOverlaps2 + detectedOverlaps3 ))
        else:
            collisionDetected, detectedOverlaps = CylinderCollisionDetect(state,newType,newC1,newC2)  
    elif npShape == 3:
        if hasPBC == True:
            PlaneCollisionDetected = 0
            for xo in [-2*planeHalfLength, 0, 2*planeHalfLength]:
                for yo in [-2*planeHalfLength, 0, 2*planeHalfLength]:
                    newcollisionDetected, newdetectedOverlaps = PlaneCollisionDetect(state,newType,newC1,newC2,xo,yo) 
                    collisionDetected = max(collisionDetected , newcollisionDetected)
                    detectedOverlaps = detectedOverlaps + newdetectedOverlaps
            detectedOverlaps = list(set(detectedOverlaps))
        else:
            collisionDetected, detectedOverlaps =  PlaneCollisionDetect(state,newType,newC1,newC2) 
    else:
        collisionDetected, detectedOverlaps =  SphereCollisionDetect(state,newType,newC1,newC2)
    return collisionDetected, detectedOverlaps

def SphereCollisionDetect(state, newType,newPhi, newTheta):
    collisionDetected = 0
    detectedOverlaps = []
    if(len(state))<1:
        return collisionDetected, detectedOverlaps
    newr = proteinData[ newType,1 ]
    radiusArray = proteinData[ state[:,0].astype(int), 1] 
    if hardSphereMode == 1:
        allowedDists = radiusArray + newr
        distSq = ((npRadius+radiusArray)*np.cos(state[:,1])*np.sin(state[:,2]) - (npRadius+newr)*np.cos(newPhi)*np.sin(newTheta) )**2 + ((npRadius+radiusArray)*np.sin(state[:,1])*np.sin(state[:,2]) - (npRadius+newr)*np.sin(newPhi)*np.sin(newTheta) )**2 + ((npRadius+radiusArray)*np.cos(state[:,2]) - (npRadius+newr)*np.cos(newTheta) )**2
        if np.any( np.sqrt(distSq)  < allowedDists):
            collisionDetected =  1
            detectedOverlaps =  np.nonzero( np.sqrt(distSq)  < allowedDists )[0].tolist()
        return collisionDetected, detectedOverlaps
            
    #second test: surface projection detection
    minimumAngle = np.arcsin(radiusArray/(radiusArray+npRadius))  + np.arcsin(newr/(newr+npRadius))
    angleDist = np.arccos( np.cos(newTheta)*np.cos(state[:,2]) + np.sin(newTheta)*np.sin(state[:,2])*np.cos( np.abs(newPhi - state[:,1]) )    )
    if np.any(angleDist < minimumAngle):
        collisionDetected =  1
        detectedOverlaps =    np.nonzero(angleDist < minimumAngle)[0].tolist()
    return collisionDetected, detectedOverlaps


def shuffleCollisionDetect(state, shiftedNum ,newPhi, newTheta):
    if(len(state))<2:
        return 0
    proteinType = int(state[shiftedNum,0])
    newr = proteinData[ proteinType,1 ]
    otherStateMask = np.ones((state[:,0]).shape,bool)
    otherStateMask[shiftedNum] = False
    otherState = state[otherStateMask]
    radiusArray = proteinData[ otherState[:,0].astype(int), 1]
    heightArray = np.where(radiusArray  > newr, radiusArray , newr   ) + npRadius
    dists = heightArray * np.sqrt((  np.cos(otherState[:,1]) * np.sin(otherState[:,2]) -np.cos(newPhi) * np.sin(newTheta)   )**2 + (  np.sin(otherState[:,1]) * np.sin(otherState[:,2]) -np.sin(newPhi) * np.sin(newTheta))**2  + ( np.cos(otherState[:,2]) - np.cos(newTheta)  )**2 )
    allowedDists = proteinData[  otherState[:,0].astype(int) ,1] + newr
    if np.any(allowedDists > dists):
        return 1
    return 0


def CylinderCollisionDetect(state, newType,newPhi, newZ,zoffset=0):
    collisionDetected = 0
    detectedOverlaps = []
    if(len(state))<1:
        return collisionDetected,detectedOverlaps
    newr = proteinData[ newType,1 ]
    radiusArray = proteinData[ state[:,0].astype(int), 1]
    #first pass: detect physical overlap
    allowedDists = radiusArray + newr
    overlapCondition = np.sqrt( ( (npRadius+radiusArray)*np.cos(state[:,1]) - (newr+npRadius)*np.cos(newPhi)  )**2  +  ( (npRadius+radiusArray)*np.sin(state[:,1])  - (newr+npRadius)*np.sin(newPhi)  )**2   +  (  state[:,2] + zoffset  - newZ  )**2      )     < allowedDists
    if np.any( overlapCondition):
        collisionDetected = 1
        detectedOverlaps = np.nonzero(overlapCondition)[0].tolist()
    #At this point in the code, no hard-sphere overlaps between proteins have been detected. so if we're enabling this mode we can now return 0
    if hardSphereMode == 1:
        return collisionDetected,detectedOverlaps
    #second pass: project all sphere-pairs up to the same radial distance such that the larger is still touching the cylinder and check again for overlap
    minRD = np.where( radiusArray > newr, radiusArray, newr)
    overlapCondition =  np.sqrt( ( (npRadius+minRD)*np.cos(state[:,1]) - (minRD+npRadius)*np.cos(newPhi)  )**2  +  ( (npRadius+minRD)*np.sin(state[:,1])  - (minRD+npRadius)*np.sin(newPhi)  )**2    +  (  state[:,2] + zoffset  - newZ  )**2      )     < allowedDists
    if np.any(overlapCondition):
        detectedOverlaps = np.nonzero(overlapCondition)[0].tolist()
        collisionDetected  = 1
    return collisionDetected,detectedOverlaps


def PlaneCollisionDetect(state, newType,newC1, newC2,xoffset=0,yoffset=0):
    collisionDetected = 0
    detectedOverlaps = []
    if(len(state))<1:
        return collisionDetected,detectedOverlaps
    newr = proteinData[ newType,1 ]
    radiusArray = proteinData[ state[:,0].astype(int), 1]
    #first pass: detect physical overlap
    allowedDists = radiusArray + newr
    overlapCondition = np.sqrt(    (  state[:,1] + xoffset  -  newC1   )**2  +  (  state[:,2]  +yoffset-  newC2   )**2   +  (  radiusArray - newr )**2      )     < allowedDists
    if np.any( overlapCondition ):
        detectedOverlaps = np.nonzero(overlapCondition)[0].tolist()
        collisionDetected  = 1
    #At this point in the code, no hard-sphere overlaps between proteins have been detected. so if we're enabling this mode we can now return 0
    if hardSphereMode == 1:
        return collisionDetected,detectedOverlaps
    #second pass: as before, except ignoring z such that all proteins have their centers on the same plane
    minRD = np.where( radiusArray > newr, radiusArray, newr)
    overlapCondition = np.sqrt(    (  state[:,1] + xoffset -  newC1   )**2  +  (  state[:,2]  +yoffset-  newC2   )**2       )     < allowedDists
    if np.any(overlapCondition ):
        detectedOverlaps = np.nonzero(overlapCondition)[0].tolist()
        collisionDetected  = 1
    return collisionDetected,detectedOverlaps


def bindingArea(rnp,ri):
    if npShape == 1:
        area = bindingAreaSphere(rnp,ri)
    elif npShape == 2:
        area = bindingAreaCylinder(rnp,ri)
    elif npShape == 3:
        area = bindingAreaSphere(1000,ri)
    else:
        area = bindingAreaSphere(rnp,ri)
    return area

def bindingAreaSphere(rnp, ri):
    return 2*np.pi * (rnp**2) * (  1 - np.sqrt(rnp*(2*ri+rnp))/(ri+rnp) ) 

def bindingAreaCylinder(rnp,ri):
    return ri*rnp* 4 * np.sqrt(  rnp*(2 + rnp/ri)/ri    ) * (  scspec.ellipe(-1.0/( 2*rnp/ri + rnp*rnp/(ri*ri) )) - scspec.ellipk(-1.0/( 2*rnp/ri + rnp*rnp/(ri*ri))  )  )



def outputStateAll():
    if len(state) > 0:
        stateArr = (np.array(state)[:,0]).astype(int)
    else:
        stateArr = np.array([[]])
    resEntry = [lastUpdate]
    totalProteins = 0
    totalCoverage = 0
    print(lastUpdate,end=' ' )
    for id in proteinIDList:
        numProteins = len(stateArr[stateArr == id])
        if coarseGrainAtEnd ==0:
            print(numProteins,end=' ' )
        resEntry.append(numProteins)
        totalProteins += numProteins
        totalCoverage +=  numProteins/proteinBindingSites[id]
    if doAnalytic == 1:
        analyticState =  steadyStateA -np.matmul( sp.linalg.expm(t * aCoeffMatrix), steadyStateA)
        for id in proteinIDList:
            print(analyticState[id],end=' ' )
            resEntry.append(analyticState[id])
    print(totalProteins, " ", totalCoverage,end=' ' )
    print("")
    resList.append(resEntry)


#totalCoverage = 0

def outputState():
    if len(state) > 0:
        stateArr = (np.array(state)[:,0]).astype(int)
    else:
        stateArr = np.array([[]])
    resEntry = [lastUpdate]
    totalProteins = 0
    totalCoverage = 0
    uniqueProteinNums = np.zeros(len(uniqueProteins))
    uniqueProteinCoverage = np.zeros(len(uniqueProteins))
    lastUpdateRounded = round(lastUpdate/updateInterval) * updateInterval

    outString = str(lastUpdateRounded)
    outStringCoverage = str(lastUpdateRounded)
    print(lastUpdate,end=' ' )
    #old method: check every possible protein and see if any are bound
    #new method: only look at bound proteins - this will be quick if nbound < ntypes, which will very usually be true
    scanBound = True
    if scanBound == True and len(state) > 0:
        #print("debug starts here")
        #print(stateArr)
        for id in stateArr:
            #print(id)
            proteinName = proteinNames[id] 
            upIndex = np.nonzero(uniqueProteins == proteinName)[0][0]
            uniqueProteinNums[upIndex] += 1
            uniqueProteinCoverage[upIndex] += 1.0/proteinBindingSites[id]
            totalProteins += 1
            totalCoverage += 1.0/proteinBindingSites[id]
    else:
        for id in proteinIDList: #scan through all protein-orientations and sum up the total number of each protein
            proteinName = proteinNames[id]
            upIndex = np.nonzero(uniqueProteins == proteinName)[0][0]
            #print upIndex, uniqueProteins[upIndex]
            numProteins = len(stateArr[stateArr == id])
            uniqueProteinNums[upIndex] += numProteins
            uniqueProteinCoverage[upIndex] += numProteins/proteinBindingSites[id]
            totalProteins += numProteins
            totalCoverage +=  numProteins/proteinBindingSites[id]
    


    for upIndex in range(len(uniqueProteins)):
        if coarseGrainAtEnd ==0:
            print(uniqueProteinNums[upIndex]/numNPs,end=' ' )
            outString = outString+" "+str(uniqueProteinNums[upIndex]/numNPs)
            outStringCoverage = outStringCoverage+" "+str(uniqueProteinCoverage[upIndex]/numNPs)
        resEntry.append(uniqueProteinNums[upIndex]/numNPs )
    if doAnalytic == 1:
        analyticState =  steadyStateA -np.matmul( sp.linalg.expm(t * aCoeffMatrix), steadyStateA)
        for id in proteinIDList:
            print(analyticState[id],end=' ' )
            resEntry.append(analyticState[id])
    if coarseGrainAtEnd != 0:
        deltaGValCoverage,deltaGValNumberAverage = estimateDeltaG()
    else:
        deltaGValCoverage = ""
        deltaGValNumberAverage = ""
        
    #print( (np.abs( np.sum(proteinAdsorptionEvents,axis=1) )/dybeckNE )[:3] , end = ' ' )   
    #print( (dybeckNF*2*dybeckrs/( dybeckRateForwards + dybeckRateBackwards)), end = ' ' )  
    #print( quasiEquil[:3],end = ' ')
    #print( sufficientlyExecuted[:3],end = ' ')
    #print( empiricalAcceptance[:3],end = ' ')
    #print( lastSBEscape, end = ' ')
    #print(dybeckAlpha)
    #print(sufficientlyExecuted)
    print(totalProteins/numNPs, " ", totalCoverage/numNPs, deltaGValCoverage, deltaGValNumberAverage,end=' ' )
    print("")
    outString = outString+" "+str(totalProteins/numNPs)+" "+str(totalCoverage/numNPs)+"\n"
    outStringCoverage = outStringCoverage+" "+str(totalProteins/numNPs)+" "+str(totalCoverage/numNPs)+"\n"
    runningFile.write(outString)
    runningFileCoverage.write(outStringCoverage)
    resList.append(resEntry)
    #print( proteinData[:,0] * proteinData[:,2] * proteinBindingSites ) 
    #print(proteinCollisionEvents)
    #print(empiricalAcceptance)
def estimateDeltaG():
    #if coarseGrainAtEnd != 0:
    if len(state) > 0:
        stateArr = (np.array(state)[:,0]).astype(int)
    else:
        stateArr = np.array([[]])
        return 0,0 
    totalProteins = 0
    totalCoverage = 0
    energyAverage = 0
    for id in proteinIDList:
        numProteins = len(stateArr[stateArr == id])
        totalProteins += numProteins
        totalCoverage +=  numProteins/proteinBindingSites[id]
        energyAverage += numProteins * proteinData[id,4]
    effectiveRadius = 2*npRadius*(-2*totalCoverage**2 + 2*totalCoverage*totalProteins + totalProteins*np.sqrt(totalCoverage*(totalProteins-totalCoverage)))/( (totalProteins - 2 * totalCoverage)**2   )
    #print totalProteins, " ", totalCoverage
    if meanFieldApprox != 1:
        effectiveKeqConc = totalCoverage/(1 - totalCoverage) * np.exp(  3*totalCoverage/(1- totalCoverage) + totalCoverage**2 / (1-totalCoverage)**2 )
    else:
        effectiveKeqConc = totalCoverage/(1-totalCoverage)
    deltaGEst =  -np.log(effectiveKeqConc/ np.sum(proteinData[:,0]) )
    return deltaGEst, energyAverage/totalProteins



#how often the routine should print out information to the screen and store the number of adsorbed proteins for final output. 
#updateInterval =0.05



#these parameters control the surface diffusion - as long as the binding is reversible this can usually be disabled by not setting the -d option, as surface restructing occurs through deadsorption and readsorption.
maxShuffleTrials = 4 #attempt to diffuse each particle this many times, then give up
surfaceDiffusionCoeff = 1e5  #in nm^2/s



#if this value is set to 1 then the mean-field model is evaluated and reported at each time step.
#the mean-field model will predict different results to the hard-sphere model simulated here.
doAnalytic = 0



#defines the command line arguments. to keep things simple, only 3 can be set here: a string used to label the output file,  the radius of the NP and the file containing protein definitions
parser = argparse.ArgumentParser(description = "Parameters for the KMC simulation")
parser.add_argument('-f','--fileid',help="String to append to the filename",default="0")
parser.add_argument('-r','--radius',type=float,help="Radius of the NP",default=35.0)
parser.add_argument('-p','--proteins',help="Protein file set",default="")
parser.add_argument('-d','--diffuse',help="Enable surface diffusion",default=0)
parser.add_argument('-c', '--coarse',help="Treat input as orientations of single protein, report only total numbers",default=0)
parser.add_argument('-s','--shape',help="Shape of the NP, 1 = sphere, 2 = cylinder, 3 = plane, 4 = truncated sphere",default=1)
parser.add_argument('-m','--meanfield',help="Enable mean field approximation",default=0,type=int)
parser.add_argument('-t','--time',help="Number of hours of simulated time",default=1.0,type=float)
parser.add_argument('-n','--numnp',help="Number of NPs to simulate simultaneously", default = 1, type=int)
parser.add_argument('-x','--npconc',help="Concentration of NPs", default = 0, type=float)
parser.add_argument('-H','--hardsphere',help="Enable true hard sphere modelling", default = 0, type=int)
parser.add_argument('--demo',help="Enable the live demo mode", default = 0, type = int)
parser.add_argument('--timedelta',help="Time step [s] between showing updates", default = 1e-5, type=float)
parser.add_argument('-l','--loadfile',help="KMC file for previous run (precoating)", default="")
parser.add_argument('-P','--projectname',help="Name for project", default="testproject")
parser.add_argument('-R','--runningfile',help="Save output snapshots", default=0, type=int)
parser.add_argument('-b','--boundary',help="Boundary type: 0 = vacuum, 1 = periodic", default=1, type=int)
parser.add_argument('-D','--displace',help="Allow incoming protein to displace bound protein, nonzero = yes", default = 0, type = int)
parser.add_argument('-A','--accelerate',help="Experimental feature for quasiequilibriation scaling, nonzero = yes", default = 0, type = int)
parser.add_argument('-S','--steady', help="Fix the adsorption rate to a constant and adjust desorption rate for same steady-state behaviour", action="store_true")



doMovie = False
args = parser.parse_args()
endTime = args.time*3600

displaceSwitchOff = 0

hardSphereMode = args.hardsphere
if hardSphereMode == 1:
    print("Enabling actual hard spheres")

allowDisplace = False
displaceWater = False
if args.displace != 0:
    allowDisplace = True
    print("Enabling displacement, please make sure all binding energies are correct")
    displaceWater = True
    displaceSwitchOff = endTime + 100

useDybeckAcceleration = False
useDybeckUltra = False
if args.accelerate > 0:
    useDybeckAcceleration = True
    print("Using Dybeck acceleration")
    if args.accelerate > 1:
        useDybeckUltra = True
        print("And then accelerating further with scaling freezing")
    
forceSteadyState = False
if args.steady==True:
    print("Will attempt to find steady state and not time-resolved")
    #useDybeckAcceleration = True
    if allowDisplace == True:
        print("WARNING: Displacement mode will be automatically disabled after pre-seeding is complete and does not need to be manually passed")
    allowDisplace = True #enables displacement features if not already enabled
    displaceWater = True
    forceSteadyState = True
    displaceSwitchOff =  5 #switches off displacement mode after pre-seeding the NP 


updateInterval = args.timedelta
    

#if args.steady==True:
#    updateInterval = 0.1

npShape = int(args.shape)

if args.demo==1:
    print("Enabling demonstration mode")
    plt.ion()
    fig = plt.figure(figsize=(12.8,4.8))
    ax1 = plt.subplot(121)
    ax2 = plt.subplot(122,projection='3d')
    #ax3 = plt.subplot(133,projection='3d')
    #plt.tight_layout()
    #fig, (ax1,ax2) = plt.subplots(nrows=1,ncols=2)

cylinderHalfLength = 10 #end-to-centre length of the cylinder

hasPBC = False
#If we're using periodic boundary conditions, set a finite size (to avoid issues with large spheres). Else use the size of the NP.
if args.boundary == 1:
    planeHalfLength = 40 
    hasPBC = True
else:
    planeHalfLength = args.radius
#define the molar concentration of NPs and the number to explicitly simulate
#if you're using a vanishing concentration of NPs its best to simulate one at a time and post-average to get statistics
#but simulating multiple NPs for non-zero concentrations helps suppress fluctuations in concentration.
npConc = args.npconc #1.41e-9
numNPs = args.numnp #10


if npShape == 1:
    print("Spherical NP")
    npSurfaceArea = 4*np.pi * args.radius**2
    thetaMax = np.pi
elif npShape == 2:
    print("Cylindrical NP")
    npSurfaceArea = (2*cylinderHalfLength)*2*np.pi * args.radius
elif npShape == 3:
    print("Planar NP")
    npSurfaceArea =    4 * (planeHalfLength**2)
elif npShape == 4:
    print("Truncated sphere")
    surfaceFraction = 0.1 #must be between 0 and 1
    thetaMax =  np.arccos( 1 - 2* surfaceFraction) # np.pi/8
    vmax = 1
    vlower = 0.5*(1 + np.cos(thetaMax) )
    npSurfaceArea =  2 * np.pi * args.radius**2 * (1 - np.cos(thetaMax) )
else:
    print("Shape not recognised, defaulting to spherical")
    npShape = 1
    npSurfaceArea = 4*np.pi * args.radius**2


meanFieldApprox = args.meanfield
if meanFieldApprox != 0:
    print("Mean-field enabled")


#if this flag is set to 1 then all proteins in the input are assumed to be different orientations of the same protein. the script runs as normal, then at the end prints out an average equilibrium constant and radius
#note that this radius and equilibrium constant are dependent on the size of the NP
coarseGrainAtEnd = args.coarse

updateNum=0

#define some constants
kbtVal = 1

baseTemperature = 300 #reference value for binding energies from UA or other sources - it is assumed these are at 300K
simulationTemp = 300
tempScaleVal = simulationTemp/baseTemperature #apply a rescaling factor to energies, only used for steady-state mode and only as part of pseudo-simulated-annealing


doShuffle = args.diffuse

#define the NP area
npRadius = args.radius
proteinInput =  args.proteins

#Protein data with [0] = concentration/M, [1] =  radius/nm, [2] = kon/ M^{-1} s^{-1}, [3] = koff/s^{-1}. [4] is the binding energy, nominally -np.log(Kon/KOff)  [5] is currently unused but may be the binding area for historical reasons
#this defines the trial set of HSA, HDL and Fib using parameters from the literature. 
proteinDataOriginal = np.array([

["HDL","1.5e-5","5","3e4","3e-5", "-20.7233",str(np.pi*5**2)],
["HSA","6e-4","4","2.4e3","2e-3","-16.3004",str(np.pi*4**2)],
["Fib","8.8e-6","8.3","2e3","2e-3","-16.1181",str(np.pi*8.3**2)]

])


#kmcFileOut.write("#Name,Conc,Size,KOn,Koff,EAds,Area,C1,C2,NP\n")

precoatNames = []
precoatData = []
precoatState = []
precoatIDs = []
precoatIDLast = 0
if args.loadfile != "":
    print("Loading pre-coating file", args.loadfile)
    precoatFile = open(args.loadfile, "r")
    for line in precoatFile:
        if line[0]=="#":
            continue
        lineTerms = line.strip().split(",")
        precoatName = "PC-"+lineTerms[0]
        if precoatName not in precoatNames:
            precoatNames.append(precoatName)
            precoatData.append( [precoatName, lineTerms[1],lineTerms[2],lineTerms[3],lineTerms[4],lineTerms[5],lineTerms[6] ])
            precoatIDs.append(precoatIDLast)
            precoatIDLast += 1
        precoatID = precoatNames.index(  precoatName)
        precoatState.append(  [  precoatName, float(lineTerms[7]),float(lineTerms[8]),float(lineTerms[9])    ] )
    precoatFile.close()

if proteinInput == "":
    proteinDataInput = proteinDataOriginal
else:
    print("loading from file ", proteinInput)
    proteinDataInput = np.genfromtxt(proteinInput, dtype=str)
    #print(proteinDataInput)
#proteinData = proteinDataInput[:,1:].astype(float)

if len(precoatData)>0:
    proteinDataInput = np.concatenate(( np.array(precoatData), proteinDataInput))
    #print(proteinDataInput)


if proteinDataInput.ndim == 1:
    proteinDataInput = np.array([proteinDataInput])
proteinData = proteinDataInput[:,1:].astype(float)

#If displacement is allowed then we need to apply a correction to the desorption rates to maintain the same equilibrium constant, since kads is effectively decreased
if allowDisplace == True and displaceWater == True:
    proteinData[:,3] = proteinData[:,3] * scspec.expit(  -(proteinData[:,4]) )

proteinNamesAll = proteinDataInput[:,0]
#print(proteinNamesAll)
print("Producing name list")
proteinNameList = []
for proteinNameOrientation in proteinNamesAll:
    nameTerms = proteinNameOrientation.split(":")
    proteinNameList.append(nameTerms[0])
proteinNames = np.array(proteinNameList,dtype=str)

uniqueProteinList,uniqueProteinIDs = np.unique(proteinNames,return_index=True)
uniqueProteins = proteinNames[np.sort(uniqueProteinIDs)]
#print uniqueProteinIDs
proteinBindingSites = npSurfaceArea /( bindingArea(npRadius, proteinData[:,1]) )

largestProteinRadius = np.amax(proteinData[:,1])

#Rescale the adsorption rate to be uniform and rescale the desorption rate so that it behaves correctly in the steady-state
if forceSteadyState == True:
    print(proteinData[:,2]/proteinData[:,3] )

    slowestDesorbIndex = np.argmin( proteinData[:,3] )
   
    print("Original min desorption rate: ", proteinNamesAll[ slowestDesorbIndex ] , proteinData[ slowestDesorbIndex]  , (1.0/proteinData[ slowestDesorbIndex,3])/(24*3600.0) )
    collisionRateConst = 1.0  #/len(proteinData)
    desorbRateInitial = np.copy( proteinData[:,3] )
    #collideRateInitial = (proteinData[:,0] * proteinData[:,2] * proteinBindingSites * numNPs)
    #desorbRateAdjusted = collisionRateConst * desorbRateInitial/collideRateInitial
    
    bindingAreaSet = bindingArea(npRadius, proteinData[:,1] )
    #print(proteinData[:,4] - np.log(proteinData[:,0]) , bindingAreaSet )
    #energyDensity = (proteinData[:,4] - np.log(proteinData[:,0])  ) * proteinBindingSites #since this is inversely proportionate to area
    #print(energyDensity)
    #zeroKelvinProtein = np.argmin( energyDensity)
    totalEnergyFull =  (proteinData[:,4] - np.log(proteinData[:,0])  )* proteinBindingSites
    #print(totalEnergyFull)
    zeroKelvinProtein = np.argmin(totalEnergyFull)
    #print(np.sqrt(bindingAreaSet/np.pi))
    #quit()
    strongBindArea = bindingAreaSet[zeroKelvinProtein]
    #print(zeroKelvinProtein)
    steadyStateCoverage = 0.5
    oneMCov = 1 - steadyStateCoverage
    areaFactorSq = bindingAreaSet * steadyStateCoverage**2 / strongBindArea**2
    term1 = bindingAreaSet * steadyStateCoverage/strongBindArea
    term2 = 2 * np.sqrt(bindingAreaSet) * steadyStateCoverage / np.sqrt(strongBindArea)
    steadyStateInsertionProb = oneMCov *  np.exp( -1.0*(  areaFactorSq   /  (oneMCov**2)      )    - 1.0*( (term1+term2) /  (oneMCov)   )   ) 
    print(steadyStateInsertionProb)
    steadyStateScale = steadyStateInsertionProb/np.amax(steadyStateInsertionProb)
    print(steadyStateScale)
    #quit()
    #quit()


    kaInitial = np.copy(proteinData[:,2] )
    print(collisionRateConst)
    print(proteinData[:,0])
    print(proteinBindingSites)
    kaAdjustFactor = collisionRateConst/(   proteinData[:,0] * proteinBindingSites  * proteinData[:,2] ) 
    #kaAdjustFactor =  kaAdjustFactor /  steadyStateScale 
    #kaAdjustFactor = collisionRateConst/( proteinData[:,0] )
    print("ka adjust: ", kaAdjustFactor)
    rescaleAdsorption = True
    if rescaleAdsorption == True:
        #print("Scaled rates: ", proteinData[:,2] * proteinData[:,0] * proteinBindingSites )
        #rescale KA to produce approx. equivalent collision rates for all species
        proteinData[:,2] = proteinData[:,2] * kaAdjustFactor
        proteinData[:,3] = proteinData[:,3] * kaAdjustFactor
        #totalIncomingRate = np.sum( proteinData[:,2] * proteinData[:,0] * proteinBindingSites )
        #extraFactor = totalIncomingRate
        #print("Rescaled ka: ", proteinData[:,2])
        #print("Rescaled rates: ", proteinData[:,2] * proteinData[:,0] * proteinBindingSites ) 
        #proteinData[:,2] = 1e3 *  proteinData[:,2] / totalIncomingRate
        #proteinData[:,3] = 1e3 *  proteinData[:,3] / totalIncomingRate
    else:
        #second option: rescale KD to be constant for all proteins 
        kdRescale = 1.0/proteinData[:,3] 
        proteinData[:,2] =  proteinData[:,2] * kdRescale
        proteinData[:,3] =  proteinData[:,3] * kdRescale
        
    totalIncomingRate = np.sum( proteinData[:,2] * proteinData[:,0] * proteinBindingSites )
    #print("total rate: ", totalIncomingRate)
    proteinData[:,2] = 1e3 *  proteinData[:,2] / totalIncomingRate
    proteinData[:,3] = 1e3 *  proteinData[:,3] / totalIncomingRate


    #print(proteinData[:,2])
    #print( proteinData[:,0] * proteinData[:,2] * proteinBindingSites ) 
    #quit()

    #quit()
    #then adjust again to account for the protein size, since larger proteins are exponentially suppressed. The leading term in the acceptance probability is exp(-Ri^2/Rj^2) and we set Rj=2 as a baseline
    #largestRadius = np.amax( proteinData[:,1] )
    #proteinRescaleFromRadius = np.exp( 1.0 * (proteinData[ :,1]**2) / largestRadius**2)
    #alternate trial: rescale based on reduction in acceptance at half-coverage for smallest
    #smallestRadius = max(1.0, np.amin( proteinData[:,1]) )
    #print(smallestRadius)
    #proteinRescaleFromRadius = np.exp( 1.0 *   ( 2 * proteinData[:,1] * (smallestRadius + proteinData[:,1])  )/ (smallestRadius**2)    )
    #proteinRescaleFromRadius = proteinRescaleFromRadius/np.amin(proteinRescaleFromRadius)
    #proteinRescaleFromRadius = 1.0

    #proteinData[:,2] = proteinData[:,2] * proteinRescaleFromRadius
    #proteinData[:,3] = proteinData[:,3] * proteinRescaleFromRadius
    print(proteinData[:,2])
    print(proteinData[:,3])
    print("Equil const:" , proteinData[:,2]/proteinData[:,3])
    finalIncomingRate = np.sum( proteinData[:,2] * proteinData[:,0] * proteinBindingSites ) 
    print("Final incoming rate: ", finalIncomingRate )
    print("Recommended timestep per event: ", 1.0/finalIncomingRate )
    print("Recommended timestep: slowest desorption: ", 1.0/np.amin( proteinData[:,3] ) )
    print("Recommended timestep: fastest desorption: ", 1.0/np.amax( proteinData[:,3] ) )
    #quit()
    #print(proteinRescaleFromRadius) 
    #quit()

    print(proteinData[:,2]/proteinData[:,3] )
    print("Final collision rates: ", proteinData[:,0] * proteinData[:,2] * proteinBindingSites)
    print("Final desorption rates: ", proteinData[:,3] )
    print("Min. desorption rate: ", np.amin( proteinData[:,3] ) )
    #quit()

    #(proteinData[:,0] * proteinData[:,2]) * proteinBindingSites * numNPs == collisionRateConst
    #proteinData[:,3]/collisionRateConst = desorbRateInitial/(proteinData[:,0] * proteinData[:,2] * proteinBindingSites * numNPs)

#print(proteinBindingSites)
if doAnalytic!=0:
    aCoeffMatrix = np.zeros( (len(proteinData), len(proteinData)))
    bVector = proteinData[:,2] * proteinData[:,0] * proteinBindingSites
    for i in range(len(proteinData)):
        aCoeffMatrix[i] =  - proteinData[i,2] * proteinData[i,0] * proteinBindingSites[i]/proteinBindingSites[:]
        aCoeffMatrix[i,i] -= proteinData[i,3] 
    steadyStateA  = proteinBindingSites* proteinData[:,0]*proteinData[:,2]/proteinData[:,3]/(1 + np.sum((proteinData[:,0]*proteinData[:,2]/proteinData[:,3])) )


#calculate the collision rates used to determine which protein to add = concentration * k_on * NP surface area / protein cross sectional area
collisionRates = (proteinData[:,0] * proteinData[:,2]) * proteinBindingSites * numNPs

print(collisionRates)
#quit()
#To find the steady-state we make everything arrive with the same probability
#if forceSteadyState == True:
#    collisionRateFactor = np.copy(collisionRates)
#    collisionRates = np.zeros_like(collisionRates) + 1.0 

additionProb = np.cumsum(collisionRates)

slowestAdsorber = min(10,np.amin(collisionRates) )

state = [] #list containing the protein ID and xyz coordinates



proteinIDList = np.arange(len(proteinData))



#orientationFactors = np.zeros_like(proteinData[:,0])

#for id in proteinIDList:
proteinTotalConcs = np.zeros_like(proteinData[:,0] )



dybeckDelta  = 1.0 #factor for the error tolerance for determining quasi-equilibriation



#estimate the length of the event


dybeckNE = 40#np.amax(np.round(0.25*proteinBindingSites).astype(int)) #number of events to record for determining QE

dybeckExecutionNumber = 20

#if dybeckNE < 20:
#    dybeckNE = 20
print("Event queue size: ", dybeckNE)

dybeckNF = 5  #target number of events per superbasin

if forceSteadyState == True:
    dybeckNF = 20

dybeckAutoNF = True

dybeckMinAlpha = 0.01


timeInSuperbasin = 1e-20
proteinAdsorptionEvents = np.zeros( (len(proteinData[:,0]), dybeckNE)) #Nprot x dybeckNE matrix for tracking events, 1 = adsorb, -1 = desorb
proteinCollisionEvents = np.zeros( (len(proteinData[:,0]), 2)) #Nprot x 2 matrix for tracking collisions, c1 = successful c2 = total
proteinCollisionEvents[:,:] = 0.0
proteinCollisionMemory = 0.999999


#proteinCollisionMemory = proteinData[:,2] 



empiricalAcceptance = 1.0  + proteinCollisionEvents[:,0]/(1e-6 + proteinCollisionEvents[:,1])  #for an empty corona there's a unit chance of acceptance

proteinSBExecutions= np.zeros_like(proteinData[:,0])
proteinSBExecutionsFormer= np.zeros_like(proteinData[:,0])
lastSBEscape = ""


dybeckRateForwards = np.zeros_like(proteinData[:,0]) + 1e-10
dybeckRateBackwards = np.zeros_like(proteinData[:,0]) + 1e-10
dybeckRateBackwardsTerm = np.zeros_like(proteinData[:,0])
#dybeckDelta = np.zeros_like(proteinData[:,0]) + dybeckDeltaMin

dybeckAlpha = np.ones_like(proteinData[:,0])
dybeckAlphaFrozen = np.ones_like(proteinData[:,0])
dybeckIsFrozen = False

quasiEquil = np.logical_and(  np.abs( np.sum(proteinAdsorptionEvents,axis=1) )/dybeckNE < dybeckDelta   , np.sum(np.abs(proteinAdsorptionEvents),axis=1) >= dybeckNE ) 


sufficientlyExecuted = np.logical_and(  proteinSBExecutions >=dybeckExecutionNumber, quasiEquil)

totalIDs = len(proteinIDList)
print("Preparing total concentrations")


#reverse-engineer the orientation factor applied to each protein so the concentrations are properly adjusted if depletion is enabled
'''
numIDsDone = 0
for id in proteinIDList:
    proteinName = proteinNames[id]
    proteinTotalConcs[ proteinNames == proteinName   ] += proteinData[id,0]
    if numIDsDone % 1000 == 0:
        print( numIDsDone,"/",totalIDs)
    numIDsDone += 1
''' 
for uniqueName in uniqueProteins:
    print(uniqueName)
    nameMask = proteinNames == uniqueName
    totalNameConc = np.sum(  proteinData[  nameMask     , 0]  )
    proteinTotalConcs[nameMask] = totalNameConc
    print("Set all ", uniqueName , " to base total concentration", totalNameConc)
orientationFactors = proteinData[:,0]/(1e-15+proteinTotalConcs)
#print(orientationFactors)


resList =[]
t = 0
lastUpdate =0
surfaceCoverage = 0 

boundProteinAll = np.zeros( len(proteinData[:,0]) )





#Add in pre-existing proteins from the save file, if any

for i in range(len(precoatState)):
    existingProteinName = precoatState[i][0]
    #print(existingProteinName)
    #print(proteinNamesAll)
    newProteinID = np.where( proteinNamesAll == existingProteinName)[0][0]
    #print(newProteinID)
    #quit()
    newPhi = precoatState[i][1]
    newC2 = precoatState[i][2]
    collidingNP = precoatState[i][3]
    state.append([newProteinID,newPhi,newC2 ,collidingNP ])
    surfaceCoverage += 1.0/(  numNPs*  proteinBindingSites[newProteinID])
    boundProteinAll[  proteinNames == proteinNames[newProteinID]   ]  += 1






if proteinInput == "":
    outputTag = "hsahdlfib"
else:
    outputTag = proteinInput.split("/")[-1]

if meanFieldApprox == 1:
    mfTag = "mf"
else:
    mfTag = "hs"

coronaSaveDir = "CoronaPredictionProjects/"+args.projectname+"/coronakmc_output"

os.makedirs(coronaSaveDir,exist_ok=True)

if doMovie==True:
    os.makedirs(coronaSaveDir+"/movie",exist_ok=True)

finalName = coronaSaveDir+"/kmc_"+outputTag+"_"+str(npRadius)+"_s"+str(doShuffle)+"_"+mfTag+"_"+args.fileid+".txt"
runningName = coronaSaveDir+"/kmc_running_"+outputTag+"_"+str(npRadius)+"_s"+str(doShuffle)+"_"+mfTag+"_"+args.fileid+".txt"
runningNameCoverage = coronaSaveDir+"/kmc_running_"+outputTag+"_"+str(npRadius)+"_s"+str(doShuffle)+"_"+mfTag+"_"+args.fileid+"_coverage.txt"

runningFile = open(runningName, "w")
runningFileHeader = ",".join( [a for a in uniqueProteins])
runningFile.write("t/s,"+runningFileHeader+",total,total_coverage\n")

runningFileCoverage = open(runningNameCoverage, "w")
runningFileCoverageHeader = ",".join( ["s_"+a for a in uniqueProteins])
runningFileCoverage.write("t/s,"+runningFileCoverageHeader+",total,total_coverage\n")


print("t/s",end=' ')
for proteinName in uniqueProteins:
    print(proteinName,end=' ')
print("total", "total_coverage",end=' ')
if coarseGrainAtEnd != 0:
    print("DeltaG")
else:
    print("")

localRateRescale = np.ones_like(proteinData[:,0])
disabledDisplace = False
#Kinetic Monte Carlo approach
while t < endTime:

    if disabledDisplace == False and t > displaceSwitchOff:
        allowDisplace = False
        displaceWater = False
        #print(proteinData[:,3])
        #print(proteinData[:,4], scspec.expit(  -(proteinData[:,4]))
        proteinData[:,3] = proteinData[:,3] / scspec.expit(  -(proteinData[:,4]) ) #remove the water-displacement-adjustment to KD
        #print(proteinData[:,3] )
        disabledDisplace = True
        print("Displace off, estimated time for settling:" + str(1.0/averageDesorbRate))

    empiricalProbs =  np.where(  proteinCollisionEvents[:,1] < 1, 1.0,  empiricalAcceptance)
    #localRateRescale = 1.0/np.clip( empiricalProbs, 1e-8, 1.0) #rescale adsorption and desorption rates: proteins being frequently rejected (eP->0) are boosted
    #localRateRescale = np.clip(empiricalProbs,0.1,1.0)  # frequent rejections are slowed - selected less often but removed less often

    #localRateRescale = 1.0 /( empiricalProbs )  #or rescale so unlikely ones occur more but fall off quicker to try to accelerate convergence
    #print(empiricalProbs)
    #print(empiricalProbs * proteinData[:,0] * (proteinData[:,2] / proteinData[:,3] )* proteinBindingSites ) 

    #print(empiricalProbs,localRateRescale)

    #dybeck acceleration https://doi.org/10.1021/acs.jctc.6b00859
    resetSuperbasin = False
    '''
    if forceSteadyState == True:
        timeScaleConst = 0.2
        newTemp  = 300.0 + np.exp( -t/timeScaleConst) * 100.0
        tempScaleVal = newTemp/300.0
    '''
    
        #print("rescaling by ", localRateRescale)

    #get the quasi-equilibriated states for all protein-orientations: those which are current approximately equal in terms of adsorption and desorption and which have had a sufficient number of events
    #quasiEquil = np.logical_and(  np.abs( np.sum(proteinAdsorptionEvents,axis=1) )/dybeckNE < dybeckDelta   , np.sum(np.abs(proteinAdsorptionEvents),axis=1) >= dybeckNE ) 
    quasiEquil = np.logical_and(  np.abs( np.sum(proteinAdsorptionEvents,axis=1) )  < dybeckDelta*np.sqrt(dybeckNE)   , np.sum(np.abs(proteinAdsorptionEvents),axis=1) >= dybeckExecutionNumber ) 
    
    
    #print(proteinAdsorptionEvents[2])
    #select those which are "sufficiently executed" to be rescaled
    sufficientlyExecuted = np.logical_and(  proteinSBExecutions >= dybeckExecutionNumber, quasiEquil)
    
    dybeckrs = 0.5*np.sum(   (dybeckRateForwards + dybeckRateBackwards) * np.where(sufficientlyExecuted, 0, 1) )
    #get the scaling factor for processes to decelerate fast ones
    if timeInSuperbasin > 1e-19 and dybeckIsFrozen == False:
        #dybeckAlphaUnscaled = 2*dybeckrs/( dybeckRateForwards + dybeckRateBackwards)
        dybeckNFVal = dybeckNF
        dybeckAlpha = dybeckNFVal*2*dybeckrs/( dybeckRateForwards + dybeckRateBackwards)
        
        #dybeckAlpha = dybeckAlpha*dybeckMinAlpha/np.amin(dybeckAlpha)
        dybeckAlpha = np.clip(dybeckAlpha,dybeckMinAlpha,1.0)
    elif timeInSuperbasin > 1e-19 and dybeckIsFrozen == True:
        dybeckAlpha = dybeckAlphaFrozen
    else:
        dybeckAlpha = np.ones_like( proteinData[:,0])
        #dybeckAlphaUnscaled = np.ones_like( proteinData[:,0])
    #print(dybeckAlpha)
    #main code
    leavingRates = []
    dybeckRateBackwardsTerm[:] = 0
    
    #print(t)
    stateArr = np.array(state)
    '''
    if len(state) > 0:
        stateIDs = stateArr[:,0].astype(int)
        desorbRates = proteinData[stateIDs ,3]
        dybeckRateBackwardsTerm[ stateIDs ] += desorbRates
        if useDybeckAcceleration == True:
            desorbRates = desorbRates * dybeckAlpha[ stateIDs ]
        
    else:
        desorbRates = []
    '''
    averageDesorbRate = 1e-12
    for i in range(len(state)):
        currentProtein = state[i]
        deltaG = proteinData[ currentProtein[0], 4] 
        desorbRate = proteinData[ currentProtein[0],3]
        desorbRate0 = desorbRate
        if forceSteadyState == True: #for simulated annealing convert to current desorption rate
            desorbRate = desorbRate * np.exp( deltaG*(1.0  - tempScaleVal)/tempScaleVal )
            desorbRate = desorbRate * localRateRescale[ currentProtein[0] ] 
            #print("Adjusted", desorbRate0 , "to", desorbRate, "for", deltaG, " at temp factor", tempScaleVal)
            #quit()
        if useDybeckAcceleration == True:
            desorbRate = desorbRate * dybeckAlpha[ currentProtein[0] ] 
        leavingRates.append( desorbRate  )
        averageDesorbRate += desorbRate
        dybeckRateBackwardsTerm[ currentProtein[0] ] += desorbRate
    #print(leavingRates, desorbRates)
    numAds = len(state)
    if numAds > 0:
        averageDesorbRate = averageDesorbRate/numAds
    #leavingRates = desorbRates
    if len(stateArr) < 1:
        #if forceSteadyState == True:
        #    collisionrates = proteinData[:,2]
        #else:
        collisionRates =(proteinData[:,0]  ) * proteinData[:,2] * proteinBindingSites * numNPs
        if forceSteadyState == True:
            collisionRates = collisionRates * localRateRescale


        #originalCollisionRates =(proteinData[:,0]  ) * proteinData[:,2] * proteinBindingSites * numNPs
    else:
        #boundProteinAll = np.zeros(len(proteinData[:,0]))
        #allBoundForOrientation = np.zeros( len(proteinIDList))
        #for id in proteinIDList: #scan through all protein-orientations and sum up the total number of each protein
        #    proteinName = proteinNames[id]
        #    numProteins = len(stateArr[stateArr == id])
        #    boundProteinAll[  proteinNames == proteinName   ] += numProteins
        adjustedConc = (proteinTotalConcs - npConc * boundProteinAll / numNPs)*orientationFactors
        adjustedConc = np.where(adjustedConc > 0, adjustedConc, 0)
        if forceSteadyState == True:
            collisionRates =  proteinData[:,0]  *   proteinData[:,2] * proteinBindingSites * numNPs * localRateRescale
            #print(collisionRates)
        else:
            collisionRates = adjustedConc * proteinData[:,2] * proteinBindingSites * numNPs
        #originalCollisionRates = adjustedConc * proteinData[:,2] * proteinBindingSites * numNPs
    #print boundProteinAll[::600]
    #print(collisionRates / np.array(leavingRates))
    '''
    uniqueProteinNums = np.zeros(len(uniqueProteins))

    for id in proteinIDList: #scan through all protein-orientations and sum up the total number of each protein
        proteinName = proteinNames[id]
        upIndex = np.nonzero(uniqueProteins == proteinName)[0][0]
        numProteins = len(stateArr[stateArr == id])
        uniqueProteinNums[upIndex] += numProteins


    '''
    
    
    if useDybeckAcceleration == True:
        collisionRates = collisionRates * dybeckAlpha
    
    
    #base concentration - estimated bound concentration
    #for a given protein, each orientation has concentration C0 * factor , where C0 is sum of all orientation-specific concentrations
    # C[theta,phi] = C0*f[ theta,phi] => f[theta,phi] = C[theta,phi]/C0
    # the time-dependent concentration is the total concentration - npConc numProteinsBound/numNPs
    # C[t, theta, phi] = (C0 - npConc*proteinsPerNP[t]) * f[theta,phi]
    #                  = max(C0 - npConc*proteinsPerNP[t],0) * C[t=0,theta,phi ]/C0

    #collisionRates = (proteinData[:,0] - npConc*numProteinsBound/numNPs   ) * proteinData[:,2] * proteinBindingSites * numNPs 
    #additionProb = np.cumsum(collisionRates)
    allProcesses = np.concatenate((collisionRates, leavingRates))
    processCS = np.cumsum(allProcesses)
    chosenProcess= np.argmax( processCS > processCS[-1] * np.random.random() )
    


    
    deltat= 1.0/(processCS[-1]) * np.log( 1.0 / np.random.random() )
    deltatTotal = deltat
    # print deltatTotal
    #print out the state of the system at the specified updating interval, using an empirical correction for the current acceptance rate of each protein type so that Dybeck rates are adsorption and not collision
    #when forcing steady-state these are also used to rescale rate constants to favour proteins which can actually adsorb
    empiricalAcceptance = proteinCollisionEvents[:,0]/(1e-6 + proteinCollisionEvents[:,1])
    dybeckRateForwards += ( (collisionRates*deltatTotal * empiricalAcceptance)/dybeckAlpha) 
    dybeckRateBackwards += ((dybeckRateBackwardsTerm*deltatTotal)/dybeckAlpha) 
    
    while deltat > updateInterval:
        deltat = deltat - updateInterval
        t += updateInterval
        if t > lastUpdate:
            outputState()
            lastUpdate += updateInterval
        print 
        
    #we then (if required) print out the state at the specified update time before making the most recent change, since this change by definition takes place after this update time.
    t += deltat
    timeInSuperbasin += deltatTotal
    
    
    
    
    if t > lastUpdate:
        outputState()
        #print( np.nonzero(sufficientlyExecuted) , np.nonzero(quasiEquil)  )
        #if len(np.nonzero(quasiEquil) ) > 0:
        #    print( proteinSBExecutions[ np.nonzero(quasiEquil)  ] ,  dybeckAlpha[ np.nonzero(quasiEquil)  ] , dybeckAlphaUnscaled[ np.nonzero(quasiEquil)  ]  )
        #print state
        lastUpdate += updateInterval
        if args.demo==1:
            shapeClass = "sphere"
            if npShape == 2:
                shapeClass = "cylinder"
            if npShape == 3:
                shapeClass = "plane"
            resArray = np.array(resList)
            ax1.clear()
            for i in range(len(uniqueProteins)):
                ax1.plot( resArray[:,0]/60, resArray[:,i+1] ,label=uniqueProteins[i],color="C"+str(i))
            ax1.set_xlabel("Time [min]")
            ax1.set_ylabel("Num. absorbed")
            ax1.legend()
            ax2.clear()
            if shapeClass == "sphere":
                u = np.linspace(0, 2 * np.pi, 20)
                v = np.linspace(0, thetaMax, 20)
                npx = npRadius * np.outer(np.cos(u), np.sin(v))
                npy = npRadius * np.outer(np.sin(u), np.sin(v))
                npz = npRadius * np.outer(np.ones(np.size(u)), np.cos(v))
                plotBound = npRadius + largestProteinRadius
            elif shapeClass == "cylinder":
                u = np.linspace(0,2*np.pi,20)
                v = np.linspace(-1,1,20,endpoint=True)
                npx = npRadius * np.outer(np.cos(u), np.ones(np.size(v)))
                npy = npRadius * np.outer(np.sin(u), np.ones(np.size(v)))
                #print("NP surface element coords")
                #print(x)
                #print(y)
                npz = cylinderHalfLength* np.outer(np.ones(np.size(u)), v)
                plotBound = npRadius + largestProteinRadius
                #print(z)
            elif shapeClass == "plane":
                u = np.linspace(-1,1,20,endpoint=True)
                v = np.linspace(-1,1,20,endpoint=True)
                npx = planeHalfLength * np.outer(u, np.ones(np.size(v)))
                npy = planeHalfLength * np.outer(np.ones(np.size(u)), v)
                #print("NP surface element coords")
                #print(x)
                #print(y)
                npz =  0.0 * np.outer(np.ones(np.size(u)), v)
                plotBound = planeHalfLength
            #ax2.plot_surface(x, y, z,color='gray')            
            #ax2.scatter(npx,npy,npz,color='gray')
            #print("setting plotBound: ", plotBound, " from ", npRadius, " and large protein: ", largestProteinRadius)
            ax2.set_xlim3d( -plotBound ,plotBound)
            ax2.set_ylim3d( -plotBound , plotBound)
            if shapeClass == "plane":
                ax2.set_zlim3d(0 , plotBound)
            else:
                ax2.set_zlim3d( -plotBound , plotBound)
            #ax3.plot_surface(x, y, z,color='gray')
            #ax3.set_xlim3d( -plotBound ,plotBound)
            #ax3.set_ylim3d( -plotBound , plotBound)
            #ax3.set_zlim3d( -plotBound , plotBound)
            proteinOffset = 0.3
            for i in range(len(state)):
                currentProtein = state[i]
                proteinName = proteinNames[currentProtein[0]]
                upIndex = np.nonzero(uniqueProteins == proteinName)[0][0]
                proteinRadius = proteinData[ currentProtein[0] ,1]
                if shapeClass == "sphere":
                    proteinX = (proteinRadius + npRadius +proteinOffset)*np.cos( currentProtein[1] ) * np.sin(currentProtein[2])
                    proteinY = (proteinRadius + npRadius+proteinOffset)*np.sin( currentProtein[1] ) * np.sin(currentProtein[2])
                    proteinZ = (proteinRadius + npRadius+proteinOffset)*np.cos(currentProtein[2] )
                    #print("Protein center: ", proteinX, proteinY, proteinZ, currentProtein[1] * 180.0/np.pi, currentProtein[2] * 180.0/np.pi)
                elif shapeClass == "plane":
                    proteinX = currentProtein[1]
                    proteinY = currentProtein[2]
                    proteinZ = (proteinRadius + proteinOffset)
                else:
                    proteinX = (proteinRadius + npRadius+proteinOffset)*np.cos( currentProtein[1] ) 
                    proteinY = (proteinRadius + npRadius+proteinOffset)*np.sin( currentProtein[1] ) 
                    proteinZ = (currentProtein[2] )
                u = np.linspace(0, 2 * np.pi, 20)
                v = np.linspace(0, np.pi, 20)
                x = proteinX + proteinRadius * np.outer(np.cos(u), np.sin(v))
                y = proteinY + proteinRadius * np.outer(np.sin(u), np.sin(v))
                z = proteinZ + proteinRadius * np.outer(np.ones(np.size(u)), np.cos(v))
                ax2.plot_surface(x, y, z,color="C"+str(upIndex))
                #thetaProject =np.arctan2(  np.sqrt(x**2 + y**2) ,z  )
                #phiProject = np.arctan2(  y, x )
                #shadowOffset = 0.1
                #pointRadius = np.sqrt(x**2+y**2+z**2)
                #xp = (shadowOffset+npRadius)*x/(pointRadius)
                #yp = (shadowOffset+npRadius )* y/(pointRadius)
                #zp = (shadowOffset+npRadius) * z /(pointRadius)
                #ax3.plot_surface(xp,yp,zp, color="C"+str(upIndex)) 
            #ax2.scatter(npx,npy,npz,color='gray')
            plt.pause(0.05)
            if doMovie == True:
                plt.savefig(coronaSaveDir+"/movie/frame_"+str(updateNum)+".png")
            updateNum+=1
    #next, diffuse proteins around the surface if enabled
    if doShuffle != 0 and npShape==1:
        #diffuse proteins around the surface
        shuffleOrder = range(len(state))
        random.shuffle(shuffleOrder)
        diffusionSigma = np.sqrt((surfaceDiffusionCoeff*2*deltatTotal))
        for i in shuffleOrder:
            numAttempt = 0
            while numAttempt < maxShuffleTrials:
                newtheta = state[i][2] + np.random.normal(0, np.sqrt(surfaceDiffusionCoeff*2*deltatTotal / (npRadius**2)) ) #using result from PRE 96 022606
                phiVariance =  (1.0/(state[i][2]**2) + 1.0/((state[i][2] - np.pi)**2)   )* (2*surfaceDiffusionCoeff*deltatTotal)/(npRadius**2)
                newphi = state[i][1] + np.random.normal(0,np.sqrt(phiVariance))
                #bring phi to interval [0,2pi]
                while newphi > 2*np.pi:
                    newphi -= 2*np.pi
                while newphi < 0:
                    newphi += 2*np.pi
               #bring theta to interval [0,2pi] then to interval [0,pi] by reflection
                while newtheta > 2*np.pi:
                    newtheta -= 2*np.pi
                while newtheta < 0:
                    newtheta += 2*np.pi
                if newtheta > np.pi:
                    newtheta = 2*np.pi - newtheta
                stateArr = np.array(state)
                currentNP =  state[i][3]
                isCollision = shuffleCollisionDetect(stateArr[ stateArr[:,3] == currentNP], i ,newphi, newtheta) #check to see if the proposed step would lead to a collision
                if isCollision == 0:
                    state[i][1] = newphi
                    state[i][2] = newtheta
                    numAttempt = maxShuffleTrials
                else:
                    numAttempt+=1


    #this finally brings the system up to the correct point in time for the protein to adsorb or desorb, which is then done.
    #print chosenProcess
    if chosenProcess < len(collisionRates):
        #attempt to add a protein but reject if this would cause an overlap (i.e. the incoming protein sticks with prob 1 to bare NP, reflected otherwise)
        newProteinID = chosenProcess
        proteinCollisionEvents[newProteinID,1] = 1 + proteinCollisionMemory*proteinCollisionEvents[newProteinID,1]
        newPhi = 2*np.pi * np.random.random()
        collidingNP = np.random.randint(0,numNPs)
        if npShape == 1:
            newC2 = np.arccos( 2*np.random.random() - 1) #coordinate 2 is theta for a sphere, z for a cylinder
        elif npShape == 2:
            newC2 =  2*(np.random.random()-0.5)*(cylinderHalfLength )
        elif npShape == 3:
            newPhi =  2*(np.random.random()-0.5)*(planeHalfLength ) #for a plane we remap C1 to x and C2 to y. 2* (random - 0.5) gives -1 to +1 , so multiplying by the half length gives the full range
            newC2 =  2*(np.random.random()-0.5)*(planeHalfLength )
        elif npShape == 4:
            v = (vmax - vlower)*np.random.random() + vlower
            newC2 = np.arccos(2 * v - 1)
        else:
            newC2 = np.arccos( 2*np.random.random() - 1) #coordinate 2 is theta for a sphere, z for a cylinder
        stateArr = np.array(state)
        #stateArr[ stateArr[:,3] == collidingNP  ]
        
        if len(state) < 1:
            isOvercrowded = 0
        else:
            isOvercrowded, blockingProteinStateIDs =  adsorbCollisionDetect(stateArr[ stateArr[:,3] == collidingNP  ], newProteinID, newPhi, newC2)

        directAccept = False
        
        if isOvercrowded == 0:
            if displaceWater == False:   #if displacement of water isn't necessary, an unblocked adsorbate is automatically accepted
                directAccept = True
            else: #
                enew = proteinData[newProteinID,4] # test to see if the adsorbate can displace water
                #acceptanceProbability = np.exp( -(enew ))/( np.exp( -(enew)) + 1.0)
                acceptanceProbability = scspec.expit( -enew/tempScaleVal )
                if np.random.random() < acceptanceProbability:
                    directAccept = True
        
        if directAccept == True:
            state.append([newProteinID,newPhi,newC2 ,collidingNP ])
            surfaceCoverage += 1.0/(  numNPs*  proteinBindingSites[newProteinID])
            boundProteinAll[  proteinNames == proteinNames[newProteinID]   ]  += 1
            proteinAdsorptionEvents[newProteinID] = np.roll( proteinAdsorptionEvents[newProteinID],-1 )
            proteinAdsorptionEvents[newProteinID,-1] = 1
            proteinCollisionEvents[newProteinID,0] = 1 + proteinCollisionMemory*proteinCollisionEvents[newProteinID,0]
            proteinSBExecutions[newProteinID] +=1
            if quasiEquil[newProteinID] == False:
                resetSuperbasin = True
            #print "accepted protein ", newProteinID
        elif allowDisplace == True and isOvercrowded == True: #if the adsorbate is blocked but displacement is possible test for this
            #print("Proteins in the way state IDs: ", blockingProteinStateIDs)
            blockingProteinTypes =   stateArr[ blockingProteinStateIDs ,0 ].astype(int)  
            #print("Blocking protein types:", blockingProteinTypes)
            blockingProteinEnergies = proteinData[ blockingProteinTypes, 4]
            #print("Blocking protein energies: ", blockingProteinEnergies)
            erep = np.sum(blockingProteinEnergies)
            enew = proteinData[newProteinID,4]
            #acceptanceProbability = np.exp( - (proteinData[ newProteinID, 4] - np.sum(blockingProteinEnergies))) 
            #acceptanceProbability = np.exp( -(enew - erep))/( np.exp( -(enew-erep)) + 1.0) #calculate the acceptance probability such that the two states are canonically distributed
            acceptanceProbability = scspec.expit(  -(enew-erep)/tempScaleVal  )
            #print("Total energy:", erep, " new protein energy ", enew , "acceptance probability", acceptanceProbability)
            if   np.random.random() < acceptanceProbability:
                #print("Replacement occured")
                blockingProteinStateIDs.sort()
                blockingProteinStateIDs.reverse()
                for bpid in blockingProteinStateIDs:
                    removedProtein = state.pop( bpid)
                    surfaceCoverage -= 1.0/(numNPs*proteinBindingSites[ removedProtein[0]])
                    boundProteinAll[  proteinNames == proteinNames[ removedProtein[0]]   ]  -= 1
                    proteinAdsorptionEvents[removedProtein[0] ] = np.roll( proteinAdsorptionEvents[removedProtein[0]],-1 )
                    proteinAdsorptionEvents[removedProtein[0],-1] = -1
                    proteinSBExecutions[removedProtein[0]] +=1
                    #if quasiEquil[removedProtein[0] ] == False:
                    #    resetSuperbasin = True
                state.append([newProteinID,newPhi,newC2 ,collidingNP ])
                surfaceCoverage += 1.0/(  numNPs*  proteinBindingSites[newProteinID])
                boundProteinAll[  proteinNames == proteinNames[newProteinID]   ]  += 1
                proteinAdsorptionEvents[newProteinID] = np.roll( proteinAdsorptionEvents[newProteinID],-1 )
                proteinAdsorptionEvents[newProteinID,-1] = 1
                proteinCollisionEvents[newProteinID,0] = 1 + proteinCollisionMemory*proteinCollisionEvents[newProteinID,0]
                proteinSBExecutions[newProteinID] +=1
                if quasiEquil[newProteinID] == False:
                    resetSuperbasin = True
        else:
            proteinCollisionEvents[newProteinID,0] = proteinCollisionMemory*proteinCollisionEvents[newProteinID,0]
    else:
        #remove the protein 
        removedProtein = state.pop( chosenProcess - len(collisionRates)   )
        surfaceCoverage -= 1.0/(numNPs*proteinBindingSites[ removedProtein[0]])
        boundProteinAll[  proteinNames == proteinNames[ removedProtein[0]]   ]  -= 1
        proteinAdsorptionEvents[removedProtein[0] ] = np.roll( proteinAdsorptionEvents[removedProtein[0]],-1 )
        proteinAdsorptionEvents[removedProtein[0],-1] = -1
        proteinSBExecutions[removedProtein[0]] +=1
        #if quasiEquil[removedProtein[0] ] == False:
        #    resetSuperbasin = True
        #print "removed protein"
    if resetSuperbasin == True:
       if useDybeckUltra == True and lastSBEscape == proteinNames[newProteinID] :  #freeze scaling rates if that option is on and the same protein escapes twice in a row
           proteinSBExecutionsFormer = np.copy(proteinSBExecutions)
           dybeckAlphaFrozen = np.copy(dybeckAlpha)
           dybeckIsFrozen = True
       proteinSBExecutions[:] = 0
       dybeckRateForwards = np.zeros_like(proteinData[:,0]) + 1e-10
       dybeckRateBackwards = np.zeros_like(proteinData[:,0]) + 1e-10
       timeInSuperbasin = 1e-20
       lastSBEscape = proteinNames[newProteinID] 
       #print("Resetting superbasin due to: ", proteinNames[newProteinID])
    if useDybeckUltra == True and dybeckIsFrozen == True:
        previousMask = proteinSBExecutionsFormer > dybeckExecutionNumber #get all that were previously sufficiently executed, test if these are now sufficiently executed
        if np.all( proteinSBExecutions[previousMask] > dybeckExecutionNumber):
            dybeckIsFrozen = False

stateArray = np.array(state)
radiusArray = np.array([proteinData[ stateArray[:,0].astype(int)   , 1] ])


if npShape == 1:
    heightAboveSurface = (np.array([proteinData[ stateArray[:,0].astype(int)   , 1] ]) + npRadius) 
    outputArray =  np.vstack( (stateArray[:,0],  heightAboveSurface * np.cos(stateArray[:,1]) * np.sin(stateArray[:,2])   , heightAboveSurface * np.sin(stateArray[:,1]) * np.sin(stateArray[:,2]) ,heightAboveSurface * np.cos(stateArray[:,2]) ,radiusArray ))
elif npShape == 2:
    heightAboveSurface = (np.array([proteinData[ stateArray[:,0].astype(int)   , 1] ]) + npRadius)
    outputArray =  np.vstack( (stateArray[:,0],  heightAboveSurface * np.cos(stateArray[:,1])   , heightAboveSurface * np.sin(stateArray[:,1])  ,stateArray[:,2] ,radiusArray ))
elif npShape == 3:
    heightAboveSurface = (np.array([proteinData[ stateArray[:,0].astype(int)   , 1] ]) + 0)
    outputArray =  np.vstack( (stateArray[:,0],  stateArray[:,1]   , stateArray[:,2]  ,np.zeros_like(stateArray[:,2]) ,radiusArray ))
elif npShape == 4:
    heightAboveSurface = (np.array([proteinData[ stateArray[:,0].astype(int)   , 1] ]) + npRadius) 
    outputArray =  np.vstack( (stateArray[:,0],  heightAboveSurface * np.cos(stateArray[:,1]) * np.sin(stateArray[:,2])   , heightAboveSurface * np.sin(stateArray[:,1]) * np.sin(stateArray[:,2]) ,heightAboveSurface * np.cos(stateArray[:,2]) ,radiusArray ))
else:
    heightAboveSurface = (np.array([proteinData[ stateArray[:,0].astype(int)   , 1] ]) + npRadius) 
    outputArray =  np.vstack( (stateArray[:,0],  heightAboveSurface * np.cos(stateArray[:,1]) * np.sin(stateArray[:,2])   , heightAboveSurface * np.sin(stateArray[:,1]) * np.sin(stateArray[:,2]) ,heightAboveSurface * np.cos(stateArray[:,2]) ,radiusArray ))

outputTranspose = outputArray.T
proteinNames[ stateArray[:,0].astype(int) ]
outputTranspose 


if coarseGrainAtEnd != 0:
    if len(state) > 0:
        stateArr = (np.array(state)[:,0]).astype(int)
    else:
        stateArr = np.array([[]])
    totalProteins = 0
    totalCoverage = 0
    for id in proteinIDList:
        numProteins = len(stateArr[stateArr == id])
        totalProteins += numProteins
        totalCoverage +=  numProteins/proteinBindingSites[id]
    effectiveRadius = 2*npRadius*(-2*totalCoverage**2 + 2*totalCoverage*totalProteins + totalProteins*np.sqrt(totalCoverage*(totalProteins-totalCoverage)))/( (totalProteins - 2 * totalCoverage)**2   )
    print(totalProteins, " ", totalCoverage)
    effectiveKeqConc = totalCoverage/(1 - totalCoverage) * np.exp(  3*totalCoverage/(1- totalCoverage) + totalCoverage**2 / (1-totalCoverage)**2 )
    print("Ri:", effectiveRadius, "KEq*conc",  effectiveKeqConc,  " delta G", -np.log(effectiveKeqConc/ np.sum(proteinData[:,0]) ))





np.savetxt(finalName,np.array(resList))
#np.savetxt("corona_results_testing/"+outputTag+"_coords_"+str(npRadius)+"_s"+str(doShuffle)+".txt", outputArray.T)

print( "Number adsorbed: " , len(stateArray))
print( outputTranspose.shape )


coordFileOut = open(coronaSaveDir+"/"+outputTag+"_finalcoords_"+str(npRadius)+"_s"+str(doShuffle)+".txt", "w")
coordFileOut.write("#Protein type, x, y, z \n")
for i in range(len(stateArray)):
    coordData = outputTranspose[i] 
    outputData =[  proteinNamesAll[   stateArray[i,0].astype(int) ] , str(coordData[1]), str(coordData[2]), str(coordData[3] ) ]
    coordFileOut.write( ",".join(outputData) + "\n" )
coordFileOut.close()
#proteinDataOneLarge = np.array([
#
#["HDL","1.5e-5","20","3e4","3e-5", "-20.7233",str(np.pi*20**2)]
#])

# proteinData:
# concentration, size, kon, koff, eads, area
#state:
#ID, C1, C2, NP

kmcFileOut=open(coronaSaveDir+"/"+outputTag+"_"+str(npRadius)+".kmc","w")
kmcFileOut.write("#Name,Conc,Size,KOn,Koff,EAds,Area,C1,C2,NP\n")
for i in range(len(stateArray)):
    proteinID = stateArray[i,0].astype(int)
    proteinName = proteinNamesAll[proteinID]
    proteinDataLine = proteinData[proteinID]
    #print(proteinID,proteinName,proteinDataLine, stateArray[i])
    outputSet = [ proteinName, 0, proteinDataLine[1], proteinDataLine[2], proteinDataLine[3],proteinDataLine[4], proteinDataLine[5] ,stateArray[i,1], stateArray[i,2], stateArray[i,3]]
    outputLine = ",".join( [ str(a) for a in outputSet])
    kmcFileOut.write(outputLine+"\n")
kmcFileOut.close()
runningFile.close()
runningFileCoverage.close()
