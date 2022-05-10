'''
Given a list of protein (or other molecule) concentrations and binding energies runs an MC model to determine the coverage in the steady-state.

'''

import numpy as np
import scipy as sp
import scipy.special as scspec
import random
import argparse

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

#This function returns 1 if an overlap is detected and 0 else. 
def adsorbCollisionDetect(state,newType,newC1,newC2):
    if meanFieldApprox == 1:
        if np.random.random() < 1 - surfaceCoverage :
            if enableMFSpaceTest == 1:
                newCoverage = surfaceCoverage + 1.0/proteinBindingSites[newType]
                if newCoverage < 1:
                    return 0
                else:
                    return 1
            else:
                return 0
        else:
            return 1
    if npShape == 1:
        return SphereCollisionDetect(state,newType,newC1,newC2)
    elif npShape == 2:
        return CylinderCollisionDetect(state,newType,newC1,newC2) + CylinderCollisionDetect(state,newType,newC1,newC2,2*cylinderHalfLength) + CylinderCollisionDetect(state,newType,newC1,newC2,-2*cylinderHalfLength)
    else:
        return SphereCollisionDetect(state,newType,newC1,newC2)

def SphereCollisionDetect(state, newType,newPhi, newTheta):
    if(len(state))<1:
        return 0
    newr = proteinData[ newType,1 ]
    radiusArray = proteinData[ state[:,0].astype(int), 1] 
    if hardSphereMode == 1:
        allowedDists = radiusArray + newr
        distSq = ((npRadius+radiusArray)*np.cos(state[:,1])*np.sin(state[:,2]) - (npRadius+newr)*np.cos(newPhi)*np.sin(newTheta) )**2 + ((npRadius+radiusArray)*np.sin(state[:,1])*np.sin(state[:,2]) - (npRadius+newr)*np.sin(newPhi)*np.sin(newTheta) )**2 + ((npRadius+radiusArray)*np.cos(state[:,2]) - (npRadius+newr)*np.cos(newTheta) )**2
        if np.any( np.sqrt(distSq)  < allowedDists):
            return 1
        else:
            return 0
    #second test: surface projection detection
    minimumAngle = np.arcsin(radiusArray/(radiusArray+npRadius))  + np.arcsin(newr/(newr+npRadius))
    angleDist = np.arccos( np.cos(newTheta)*np.cos(state[:,2]) + np.sin(newTheta)*np.sin(state[:,2])*np.cos( np.abs(newPhi - state[:,1]) )    )
    if np.any(angleDist < minimumAngle):
        return 1
    return 0


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

'''
def hardSphereCollisionDetect(state, newType,newPhi, newTheta):
    if(len(state))<1:
        return 0
    newr = proteinData[ newType,1 ]
    radiusArray = proteinData[ state[:,0].astype(int), 1]
    xsq = (radiusArray+npRadius) * np.cos(state[:,1]) * np.sin(state[:,2]) 
    minimumAngle = np.arcsin(radiusArray/(radiusArray+npRadius))  + np.arcsin(newr/(newr+npRadius))
    angleDist = np.arccos( np.cos(newTheta)*np.cos(state[:,2]) + np.sin(newTheta)*np.sin(state[:,2])*np.cos( np.abs(newPhi - state[:,1]) )   $
    if np.any(angleDist < minimumAngle):
        return 1
    return 0
'''

'''

old version
def CylinderCollisionDetect(state, newType,newPhi, newZ):
    if(len(state))<1:
        return 0
    newr = proteinData[ newType,5 ]
    radiusArray = proteinData[ state[:,0].astype(int), 5]
    #first pass: detect physical overlap
    allowedDists = radiusArray + newr
    if np.any( np.sqrt( ( (npRadius+radiusArray)*np.cos(state[:,1]) - (newr+npRadius)*np.cos(newPhi)  )**2     +  ( (npRadius+radiusArray)*np.sin(state[:,1])  - (newr+npRadius)*np.sin(newPhi)  )**2    +  (  state[:,2]  - newZ  )**2     
        return 1
    #second pass: project all sphere-pairs up to the same radial distance such that the larger is still touching the cylinder and check again for overlap
    minRD = np.where( radiusArray > newr, radiusArray, newr)
    if np.any( np.sqrt( ( (npRadius+minRD)*np.cos(state[:,1]) - (minRD+npRadius)*np.cos(newPhi)  )**2     +  ( (npRadius+minRD)*np.sin(state[:,1])  - (minRD+npRadius)*np.sin(newPhi)  )**2    +  (  state[:,2]  - newZ  )**2      )     < a$
        return 1
    return 0
'''

def CylinderCollisionDetect(state, newType,newPhi, newZ,zoffset=0):
    if(len(state))<1:
        return 0
    newr = proteinData[ newType,1 ]
    radiusArray = proteinData[ state[:,0].astype(int), 1]
    #first pass: detect physical overlap
    allowedDists = radiusArray + newr
    if np.any( np.sqrt( ( (npRadius+radiusArray)*np.cos(state[:,1]) - (newr+npRadius)*np.cos(newPhi)  )**2  +  ( (npRadius+radiusArray)*np.sin(state[:,1])  - (newr+npRadius)*np.sin(newPhi)  )**2   +  (  state[:,2] + zoffset  - newZ  )**2      )     < allowedDists):
        return 1
    #At this point in the code, no hard-sphere overlaps between proteins have been detected. so if we're enabling this mode we can now return 0
    if hardSphereMode == 1:
        return 0
    #second pass: project all sphere-pairs up to the same radial distance such that the larger is still touching the cylinder and check again for overlap
    minRD = np.where( radiusArray > newr, radiusArray, newr)
    if np.any( np.sqrt( ( (npRadius+minRD)*np.cos(state[:,1]) - (minRD+npRadius)*np.cos(newPhi)  )**2  +  ( (npRadius+minRD)*np.sin(state[:,1])  - (minRD+npRadius)*np.sin(newPhi)  )**2    +  (  state[:,2] + zoffset  - newZ  )**2      )     < allowedDists):
        return 1
    return 0


def bindingArea(rnp,ri):
    if npShape == 1:
        area = bindingAreaSphere(rnp,ri)
    elif npShape == 2:
        area = bindingAreaCylinder(rnp,ri)
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




def outputState():
    if len(state) > 0:
        stateArr = (np.array(state)[:,0]).astype(int)
    else:
        stateArr = np.array([[]])
    resEntry = [lastUpdate]
    totalProteins = 0
    totalCoverage = 0
    uniqueProteinNums = np.zeros(len(uniqueProteins))
    outString = str(lastUpdate)
    print(lastUpdate,end=' ' )
    for id in proteinIDList: #scan through all protein-orientations and sum up the total number of each protein
        proteinName = proteinNames[id]
        upIndex = np.nonzero(uniqueProteins == proteinName)[0][0]
        #print upIndex, uniqueProteins[upIndex]
        numProteins = len(stateArr[stateArr == id])
        uniqueProteinNums[upIndex] += numProteins
        totalProteins += numProteins
        totalCoverage +=  numProteins/proteinBindingSites[id]
    for upIndex in range(len(uniqueProteins)):
        if coarseGrainAtEnd ==0:
            print(uniqueProteinNums[upIndex]/numNPs,end=' ' )
            outString = outString+" "+str(uniqueProteinNums[upIndex]/numNPs)
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
    print(totalProteins/numNPs, " ", totalCoverage/numNPs, deltaGValCoverage, deltaGValNumberAverage,end=' ' )
    print("")
    outString = outString+" "+str(totalProteins/numNPs)+" "+str(totalCoverage/numNPs)+"\n"
    runningFile.write(outString)
    resList.append(resEntry)


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
updateInterval =0.05



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
parser.add_argument('-s','--shape',help="Shape of the NP, 1 = sphere, 2 = cylinder",default=1)
parser.add_argument('-m','--meanfield',help="Enable mean field approximation",default=0,type=int)
parser.add_argument('-t','--time',help="Number of hours of simulated time",default=1.0,type=float)
parser.add_argument('-n','--numnp',help="Number of NPs to simulate simultaneously", default = 1, type=int)
parser.add_argument('-x','--npconc',help="Concentration of NPs", default = 0, type=float)
parser.add_argument('-H','--hardsphere',help="Enable true hard sphere modelling", default = 0, type=int)



args = parser.parse_args()
endTime = args.time*3600

hardSphereMode = args.hardsphere
if hardSphereMode == 1:
    print("Enabling actual hard spheres")

npShape = int(args.shape)

cylinderHalfLength = 50 #end-to-centre length of the cylinder


#define the molar concentration of NPs and the number to explicitly simulate
#if you're using a vanishing concentration of NPs its best to simulate one at a time and post-average to get statistics
#but simulating multiple NPs for non-zero concentrations helps suppress fluctuations in concentration.
npConc = args.npconc #1.41e-9
numNPs = args.numnp #10


if npShape == 1:
    print("Spherical NP")
    npSurfaceArea = 4*np.pi * args.radius**2
elif npShape == 2:
    print("Cylindrical NP")
    npSurfaceArea = (2*cylinderHalfLength)*2*np.pi * args.radius
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



#define some constants
kbtVal = 1
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


proteinDataSlow = np.array([

["HDL","1.5e-5","5","3e4","3e-6", "-20.7233",str(np.pi*5**2)],
["HSA","6e-4","4","2.4e3","2e-4","-16.3004",str(np.pi*4**2)],
["Fib","8.8e-6","8.3","2e3","2e-4","-16.1181",str(np.pi*8.3**2)]

])

proteinDataExtra = np.array([

["HDL","1.5e-5","5","3e4","3e-5", "-20.7233",str(np.pi*5**2)],
["HSA","6e-4","4","2.4e3","2e-3","-16.3004",str(np.pi*4**2)],
["Fib","8.8e-6","8.3","2e3","2e-3","-16.1181",str(np.pi*8.3**2)],
["FP1","6e-4","8.3","2e3","2e-3","-16.1181",str(np.pi*8.3**2)],
["FP2","1.5e-5","5.5","3e4","3e-5", "-20.7233",str(np.pi*5.5**2)]


])

proteinDataFastA = np.array([

["HDL","1.5e-5","5","3e6","3e-5", "-20.7233",str(np.pi*5**2)],
["HSA","6e-4","4","2.4e5","2e-3","-16.3004",str(np.pi*4**2)],
["Fib","8.8e-6","8.3","2e5","2e-3","-16.1181",str(np.pi*8.3**2)]

])


proteinDataBigSmall  = np.array([

["Big","1.5e-5","5","3e6","3e-5", "-20.7233",str(np.pi*5**2)],
["Small","6e-4","1","2.4e5","2e-3","-16.3004",str(np.pi*1**2)]

])

proteinDataOneLarge = np.array([

["HDL","1.5e-5","20","3e4","3e-5", "-20.7233",str(np.pi*20**2)]

])



if proteinInput == "":
    proteinDataInput = proteinDataOriginal
else:
    print("loading from file ", proteinInput)
    proteinDataInput = np.genfromtxt(proteinInput, dtype=str)
    print(proteinDataInput)
proteinData = proteinDataInput[:,1:].astype(float)
proteinNames = proteinDataInput[:,0]
uniqueProteinList,uniqueProteinIDs = np.unique(proteinNames,return_index=True)
uniqueProteins = proteinNames[np.sort(uniqueProteinIDs)]
#print uniqueProteinIDs
proteinBindingSites = npSurfaceArea /( bindingArea(npRadius, proteinData[:,1]) )

print(proteinBindingSites)
if doAnalytic!=0:
    aCoeffMatrix = np.zeros( (len(proteinData), len(proteinData)))
    bVector = proteinData[:,2] * proteinData[:,0] * proteinBindingSites
    for i in range(len(proteinData)):
        aCoeffMatrix[i] =  - proteinData[i,2] * proteinData[i,0] * proteinBindingSites[i]/proteinBindingSites[:]
        aCoeffMatrix[i,i] -= proteinData[i,3] 
    steadyStateA  = proteinBindingSites* proteinData[:,0]*proteinData[:,2]/proteinData[:,3]/(1 + np.sum((proteinData[:,0]*proteinData[:,2]/proteinData[:,3])) )


#calculate the collision rates used to determine which protein to add = concentration * k_on * NP surface area / protein cross sectional area
collisionRates = (proteinData[:,0] * proteinData[:,2]) * proteinBindingSites * numNPs
additionProb = np.cumsum(collisionRates)

state = [] #list containing the protein ID and xyz coordinates



proteinIDList = np.arange(len(proteinData))

#orientationFactors = np.zeros_like(proteinData[:,0])

#for id in proteinIDList:
proteinTotalConcs = np.zeros_like(proteinData[:,0] )

for id in proteinIDList:
    proteinName = proteinNames[id]
    #upIndex = np.nonzero(uniqueProteins == proteinName)[0][0]
    proteinTotalConcs[ proteinNames == proteinName   ] += proteinData[id,0]

orientationFactors = proteinData[:,0]/proteinTotalConcs
print(orientationFactors)


resList =[]
t = 0
lastUpdate =0
surfaceCoverage = 0 

boundProteinAll = np.zeros( len(proteinData[:,0]) )




if proteinInput == "":
    outputTag = "hsahdlfib"
else:
    outputTag = proteinInput.split("/")[-1]

if meanFieldApprox == 1:
    mfTag = "mf"
else:
    mfTag = "hs"
finalName = "corona_results_testing/kmc_"+outputTag+"_"+str(npRadius)+"_s"+str(doShuffle)+"_"+mfTag+"_"+args.fileid+".txt"
runningName ="corona_results_testing/kmc_running_"+outputTag+"_"+str(npRadius)+"_s"+str(doShuffle)+"_"+mfTag+"_"+args.fileid+".txt"

runningFile = open(runningName, "w")


print("t/s",end=' ')
for proteinName in uniqueProteins:
    print(proteinName,end=' ')
print("total", "total_coverage",end=' ')
if coarseGrainAtEnd != 0:
    print("DeltaG")
else:
    print("")
#Kinetic Monte Carlo approach
while t < endTime:
    leavingRates = []
    for i in range(len(state)):
        currentProtein = state[i]
        leavingRates.append( proteinData[ currentProtein[0] ,3]  )

    stateArr = np.array(state)
    if len(stateArr) < 1:
        collisionRates =(proteinData[:,0]  ) * proteinData[:,2] * proteinBindingSites * numNPs
    else:
        #boundProteinAll = np.zeros(len(proteinData[:,0]))
        #allBoundForOrientation = np.zeros( len(proteinIDList))
        #for id in proteinIDList: #scan through all protein-orientations and sum up the total number of each protein
        #    proteinName = proteinNames[id]
        #    numProteins = len(stateArr[stateArr == id])
        #    boundProteinAll[  proteinNames == proteinName   ] += numProteins
        adjustedConc = (proteinTotalConcs - npConc * boundProteinAll / numNPs)*orientationFactors
        adjustedConc = np.where(adjustedConc > 0, adjustedConc, 0)
        collisionRates = adjustedConc * proteinData[:,2] * proteinBindingSites * numNPs
    #print boundProteinAll[::600]
    '''
    uniqueProteinNums = np.zeros(len(uniqueProteins))
    for id in proteinIDList: #scan through all protein-orientations and sum up the total number of each protein
        proteinName = proteinNames[id]
        upIndex = np.nonzero(uniqueProteins == proteinName)[0][0]
        numProteins = len(stateArr[stateArr == id])
        uniqueProteinNums[upIndex] += numProteins


    '''
    #base concentration - estimated bound concentration
    #for a given protein, each orientation has concentration C0 * factor , where C0 is sum of all orientation-specific concentrations
    # C[theta,phi] = C0*f[ theta,phi] => f[theta,phi] = C[theta,phi]/C0
    # the time-dependent concentration is the total concentration - npConc numProteinsBound/numNPs
    # C[t, theta, phi] = (C0 - npConc*proteinsPerNP[t]) * f[theta,phi]
    #                  = max(C0 - npConc*proteinsPerNP[t],0) * C[t=0,theta,phi ]/C0

    #collisionRates = (proteinData[:,0] - npConc*numProteinsBound/numNPs   ) * proteinData[:,2] * proteinBindingSites * numNPs 
    additionProb = np.cumsum(collisionRates)
    allProcesses = np.concatenate((collisionRates, leavingRates))
    processCS = np.cumsum(allProcesses)
    chosenProcess= np.argmax( processCS > processCS[-1] * np.random.random() )
    deltat= 1.0/(processCS[-1]) * np.log( 1.0 / np.random.random() )
    deltatTotal = deltat
    # print deltatTotal
    #print out the state of the system at the specified updating interval.

    while deltat > updateInterval:
        deltat = deltat - updateInterval
        t += updateInterval
        if t > lastUpdate:
            outputState()
            lastUpdate += updateInterval
    #we then (if required) print out the state at the specified update time before making the most recent change, since this change by definition takes place after this update time.
    t += deltat
    if t > lastUpdate:
        outputState()
        #print state
        lastUpdate += updateInterval
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
        newPhi = 2*np.pi * np.random.random()
        collidingNP = np.random.randint(0,numNPs)
        if npShape == 1:
            newC2 = np.arccos( 2*np.random.random() - 1) #coordinate 2 is theta for a sphere, z for a cylinder
        else:
            newC2 =  2*(np.random.random()-0.5)*(cylinderHalfLength )
        stateArr = np.array(state)
        #stateArr[ stateArr[:,3] == collidingNP  ]
        if len(state) < 1:
            isOvercrowded = 0
        else:
            isOvercrowded =  adsorbCollisionDetect(stateArr[ stateArr[:,3] == collidingNP  ], newProteinID, newPhi, newC2)
        if isOvercrowded == 0:
            state.append([newProteinID,newPhi,newC2 ,collidingNP ])
            surfaceCoverage += 1.0/(  numNPs*  proteinBindingSites[newProteinID])
            boundProteinAll[  proteinNames == proteinNames[newProteinID]   ]  += 1
            #print "accepted protein ", newProteinID
    else:
        #remove the protein 
        removedProtein = state.pop( chosenProcess - len(collisionRates)   )
        surfaceCoverage -= 1.0/(numNPs*proteinBindingSites[ removedProtein[0]])
        boundProteinAll[  proteinNames == proteinNames[ removedProtein[0]]   ]  -= 1
        #print "removed protein"


stateArray = np.array(state)
radiusArray = np.array([proteinData[ stateArray[:,0].astype(int)   , 1] ])


if npShape == 1:
    heightAboveSurface = (np.array([proteinData[ stateArray[:,0].astype(int)   , 1] ]) + npRadius) 
    outputArray =  np.vstack( (stateArray[:,0],  heightAboveSurface * np.cos(stateArray[:,1]) * np.sin(stateArray[:,2])   , heightAboveSurface * np.sin(stateArray[:,1]) * np.sin(stateArray[:,2]) ,heightAboveSurface * np.cos(stateArray[:,2]) ,radiusArray ))
elif npShape == 2:
    heightAboveSurface = (np.array([proteinData[ stateArray[:,0].astype(int)   , 1] ]) + npRadius)
    outputArray =  np.vstack( (stateArray[:,0],  heightAboveSurface * np.cos(stateArray[:,1])   , heightAboveSurface * np.sin(stateArray[:,1])  ,stateArray[:,2] ,radiusArray ))
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


coordFileOut = open("corona_results_testing/"+outputTag+"_coords_"+str(npRadius)+"_s"+str(doShuffle)+".txt", "w")
coordFileOut.write("#Protein type, x, y, z \n")
for i in range(len(stateArray)):
    coordData = outputTranspose[i] 
    outputData =[  proteinNames[   stateArray[i,0].astype(int) ] , str(coordData[1]), str(coordData[2]), str(coordData[3] ) ]
    coordFileOut.write( ",".join(outputData) + "\n" )

coordFileOut.close()

runningFile.close()

