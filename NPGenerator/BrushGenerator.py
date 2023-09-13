import numpy as np
import scipy.special as scspec
import math
import scipy.interpolate as scint
import scipy.optimize as scopt
import matplotlib.pyplot as plt
def interpolateRoot(x,func,var):
    return func(x) - var

def probRound(x):
    xround = np.round(x)
    if np.random.random() < x-xround:
        xround+=1
    return xround.astype(int)

def updatem(m,numPoints):
    return (m*np.pi*numPoints*(2*scspec.ellipe(-m**2) - scspec.ellipk(-m**2)))/(numPoints*np.pi*scspec.ellipe(-m**2) - numPoints*np.pi*scspec.ellipk(-m**2) + m * scspec.ellipe(-m**2)**2)

def updateTheta(theta,m,j):
    return theta + ((2*j-1)*np.pi - m*((2*j-1)*np.pi/m))/(m * np.sqrt(  (1 + m**2*np.sin(theta)**2 ) ) )

def getPhiThetaSet(dist,beadRadius):
    totalPoints = math.floor( 4*np.pi * dist**2 / (np.pi*1.1*beadRadius**2)) - 1
    mVal = np.sqrt(totalPoints*np.pi)
    for i in range(10):
        mVal = updatem(mVal,totalPoints)
    thetaVals = np.arccos( 1-  ( (1+2*  np.arange(totalPoints) )/totalPoints   ))
    for i in range(21):
        thetaVals = updateTheta(thetaVals,mVal,1)
    phiVals = thetaVals*mVal
    return (phiVals % 2*np.pi,thetaVals)

#print(getPhiThetaSet(2,0.5))

#beadsPerLayer = [20,15,10,9,7,6,2]
#npRadius = 5
#beadRadius = 0.5


def generateBeadSet(beadRadius,npRadius,normedDensity):
    (phiSet,thetaSet) = getPhiThetaSet(npRadius+beadRadius,beadRadius)
    totalChainsPossible = len(phiSet)
    #print("Chain sites available: ", totalChainsPossible)
    densityFunc = np.copy(normedDensity)
    densityFunc[:,1] = densityFunc[:,1] * 4 * np.pi* ( densityFunc[:,0] + npRadius)**2
    maxBeadsPerArea = 1/(( 2*beadRadius)**2)
    occupationProb = densityFunc[0,1]/totalChainsPossible
    #print(occupationProb)
    if(occupationProb > 1):
        occupationProb = 1
        print("warning: occupation probability > 1, check density and bead radius")
    #lengthCDF =np.transpose( np.array(  [densityFunc[:,0], 1 - densityFunc[:,1]/densityFunc[0,1] ]))
    lengthCDF =np.transpose( np.array(  [densityFunc[:,0], 1 - densityFunc[:,1]/densityFunc[0,1]  ]))
    cdfInterpolate = scint.interp1d(lengthCDF[:,0],lengthCDF[:,1],bounds_error = False,fill_value=1)
    beadList = []
    for i in range(len(phiSet)):
        phiVal = phiSet[i]
        thetaVal = thetaSet[i]
        if occupationProb > np.random.random():
            randomVar = np.random.random()
            optres = scopt.root_scalar(interpolateRoot,  args=(cdfInterpolate,randomVar),bracket = (np.amin(densityFunc[:,0]),np.amax(densityFunc[:,0])), x0=(np.amax(densityFunc[:,0]-np.amin(densityFunc[:,0]))/2.0))
            numBeads = probRound(optres.root / (2*beadRadius))
            for l in range(numBeads):
                beadList.append([phiVal,thetaVal,l*2*beadRadius+beadRadius+npRadius])
    beadArray = np.array(beadList)
    beadArrayC = np.zeros_like(beadArray)
    beadArrayC[:,0] = beadArray[:,2] * np.cos(beadArray[:,0]) * np.sin(beadArray[:,1])
    beadArrayC[:,1] = beadArray[:,2] * np.sin(beadArray[:,0]) * np.sin(beadArray[:,1])
    beadArrayC[:,2] = beadArray[:,2] *   np.cos(beadArray[:,1])
    return beadArrayC
