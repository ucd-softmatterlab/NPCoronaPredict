import numpy as np
import scipy.optimize as scopt
import scipy.spatial as scispat
import scipy.integrate as scint

import argparse

def SphereProjectedRadius(RNP,radius):
    return np.sqrt( - 2 * RNP**2 * (-1 + np.sqrt(RNP*(2*radius+RNP))/(radius+RNP)) )

#probability to insert a protein into the corona of a spherical NP as a function of the current coverage and (effective) protein radius.
def InsertionProbSphere(coverage, radius):
    area = np.pi * radius**2
    perimeter = 2*np.pi * radius
    freeArea = 1 - np.sum(coverage)
    if freeArea < 0:
        freeAreaScale = 0
    else:
        freeAreaScale = 1
    return  freeArea * np.exp( - area * (np.sum(perimeter*coverage/area)**2)/(4.0*np.pi*freeArea**2) -( area * np.sum(coverage/area) + perimeter*np.sum(perimeter*coverage/(2 * np.pi * area)  ))/freeArea )

#defines the rate/kdadsorption. if an element of this is equal to zero then this protein is in (meta)stable equilibrium
def DCoverageDT(coverage, radius, konconc,koff):
    return konconc*InsertionProbSphere(coverage,radius) - coverage*koff 

#as above, but with the time included as a variable
def DCoverageDTTime(t,coverage,radius,konconc,koff):
    return konconc*InsertionProbSphere(coverage,radius) - coverage*koff


#scint.odeint(DCoverageDTTime, [0,0,0], tvals,args=(radiusArr,affinityArr), tfirst=True)

#to solve this set of equations we define a wrapper function that takes a set of variables x in [-inf,inf] and rescales them to [0,1] by s = x**2/(1+x**2).
def CoverageWrapperFunc(x, radius,konconc,koff):
    return DCoverageDT(x**2/(1+x**2), radius,konconc ,koff)

def CoverageWrapperFuncOptX(x, radius,konconc,koff):
    areaCheck = np.sum( x**2/(1+x**2))
    if areaCheck > 1:
        penaltyTerm = 0
    else:
        penaltyTerm = areaCheck**2
    return np.sum( DCoverageDT(x**2/(1+x**2), radius,konconc,koff )**2) #+ penaltyTerm


def CoverageWrapperFuncOpt(x, radius,konconc,koff):
    return np.sum( DCoverageDT(x, radius,konconc,koff )**2)


def CalculateBoltz(filename,pdbFile):
    data     = np.genfromtxt(filename)
    theta    = data[:,1]
    energy   = data[:,2]
    sinTheta = np.sin(theta * np.pi / 180.0)
    boltz    = np.sum(energy * sinTheta * np.exp(-1.0 * energy)) / np.sum(sinTheta * np.exp(-1.0 * energy))
    return boltz

def getAtomCoords(filename):
    fileIn = open(filename,"r")
    coordList = []
    for line in fileIn:
        lineData = line.split()
        if lineData[0] == "ATOM":
            coordList.append([ float(line[30:38]) , float(line[38:46]) , float(line[46:54])])
    fileIn.close()
    return np.array(coordList)

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

def getAreaHeight(pointCoords):
    coords=pointsToBeads(pointCoords)
    projectedConvexHull = scispat.ConvexHull( coords[:,(0,1) ] )
    heightAboveSurface = pointCoords[:,2] - np.mean(pointCoords[:,2]) #Distance from the centre (average z location) to the 
    return (projectedConvexHull.volume, heightAboveSurface)

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

def CalculateBoltzArea(filename,pdbFile):
    data     = np.genfromtxt(filename)
    rawCoords =  getAtomCoords(pdbFile)*0.1 # convert to nm
    theta    = data[:,1]
    energy   = data[:,2]
    sinTheta = np.sin(theta * np.pi / 180.0)
    boltz    = np.sum(energy * sinTheta * np.exp(-1.0 * energy)) / np.sum(sinTheta * np.exp(-1.0 * energy))
    weightedArea = 0
    energyByArea = 0
    AreaResList =[]
    for dataLine in data:
        (crossSectionalArea,heightAboveSurface) =  getAreaHeight(rotatePDB(rawCoords, dataLine[0]*np.pi/180, dataLine[1]*np.pi/180))
        weightedArea =  weightedArea + crossSectionalArea * np.sin( dataLine[1]*np.pi/180) * np.exp(-1.0 * dataLine[2] )
        energyByArea = energyByArea + ( dataLine[2]/crossSectionalArea * np.sin( dataLine[1]*np.pi/180) * np.exp(-1.0 * dataLine[2] )) #if we want <E/A> we need to do the averaging all in one go to avoid getting <E>/<A>
        AreaResList.append([crossSectionalArea,dataLine[2],heightAboveSurface])
    return (boltz,  weightedArea / np.sum(sinTheta * np.exp(-1.0 * energy))  ) 



def projectOntoSphere(coords,npRadius,beadRadius = 0.5):
    coords[:,0] = coords[:,0] - np.mean(coords[:,0])
    coords[:,1] = coords[:,1] - np.mean(coords[:,1])
    coords[:,2] = coords[:,2] - np.amin(coords[:,2]) + npRadius + beadRadius #after this transformation the proteins centre is at x=0,y=0 and z such that the protein is just touching the NP
    beadDist = np.sqrt( coords[:,0]**2 + coords[:,1]**2 + coords[:,2]**2)
    phiProjections = np.arctan2( coords[:,1], coords[:,0]) 
    thetaProjections = np.arctan2( np.sqrt(coords[:,1]**2 + coords[:,0]**2), coords[:,2]) + np.arcsin( beadRadius/beadDist) #projected theta, corrected for non-zero bead radius
    if len(thetaProjections) < 5: #if the protein consists of only a few beads we just calculate the area from the maximum value of theta. for a single bead, this is consistent with the value for a sphere.
        return 2*np.pi*npRadius**2*(1 - np.cos( np.amax(thetaProjections) ))
    #if not, then we employ the convex-hull routine to find the shape of the projected protein and get its area.  we map the points onto x = theta cos phi , y = theta sin phi
    projectionHull = scispat.ConvexHull(np.column_stack(( thetaProjections * np.cos(phiProjections), thetaProjections * np.sin(phiProjections))))
    phiThetaVals = np.column_stack((phiProjections,thetaProjections))
    phiThetaOutlinePoints =  phiThetaVals[projectionHull.vertices]
    pTOSorted = phiThetaOutlinePoints[np.lexsort( (phiThetaOutlinePoints[:,1], phiThetaOutlinePoints[:,0]))]
    PTOCyclic =  np.append(pTOSorted,[pTOSorted[0] + [2*np.pi,0]],0)
    areaTerms = (PTOCyclic[:-1,0] - PTOCyclic[1:,0])/(PTOCyclic[:-1,1] - PTOCyclic[1:,1]) *(   PTOCyclic[1:,1] - PTOCyclic[:-1,1] - np.sin(PTOCyclic[1:,1]) + np.sin(PTOCyclic[:-1,1]))
    return np.sum(npRadius**2 * areaTerms)

parser = argparse.ArgumentParser(description = "Parameters for corona calculation")
parser.add_argument('-r','--radius',type=float,help="Radius of the NP",default=19)
parser.add_argument('-z','--zeta',type=float,help="Zeta potential of the NP",default=0)
parser.add_argument('-a','--average',type=float,help="average over orientations (does nothing right now)",default=0)
parser.add_argument('-f','--folder',type=str,help="folder containing UA heatmaps",default="results_anatase_alltargets_sphere")

args = parser.parse_args()


proteinDataOriginal = np.array([

[1.5e-5,5,3e4,3e-5, 0,np.pi*5**2],
[6e-4,4,2.4e3,2e-3,0,np.pi*4**2],
[8.8e-6,8.3,2e3,2e-3,0,np.pi*8.3**2],
[1.5e-5,7,3e4,3e-5, 0,np.pi*7**2],
[8.8e-6,2,2e3,2e-3,0,np.pi*2**2]
])





#defines the list of proteins to include and their concentrations. here it assumes a comma-seperated list with PDBID,conc on each line.
#the PDB ID is used to find the corresponding UA output file.

concentrationData = np.genfromtxt("ProteinConcs.csv",delimiter=",",dtype=np.str,skip_header=1)


#define where to look for the required input: a pdb file and the output from united atom
pdbFolder = "pdbs/All"
#energyMapFolder = "results_swcnt_alltargets"
energyMapFolder = args.folder
#parameters for calculating rate constants
avogadroNumber = 6.022e23
kbT = 300 * 1.381e-23 
eta = 8.9e-4 #viscosity of the medium

npRadius = args.radius
npZp = args.zeta
orientationAverage = args.average

proteinIDList = []
proteinDataList = []

diffusionCoeff = 1e-11
outputSet = []

for proteinData in concentrationData:
    filename=proteinData[0]+"_"+str(int(npRadius))+"_"+str(int(npZp))+".map"
    #load in the file as before, but calculate the projected area,kon and koff separately for each orientation
    rawCoords =  getAtomCoords( pdbFolder+"/"+proteinData[0]+".pdb")*0.1
    data     = np.genfromtxt(energyMapFolder+"/"+filename)
    thetaOffset = (np.amax(data[1:,1] - data[0:-1,1]))/2.0 #UA output gives the left-edge of the bin containing the angular orientations, this corrects for that
    phiOffset = (np.amax(data[1:,0] - data[0:-1,0]))/2.0 #UA output gives the left-edge of the bin containing the angular orientations, this corrects for that
    #print thetaOffset
    data[:,1] +=thetaOffset
    data[:,0] += phiOffset
    sinThetaNorm = 1.0/np.sum( np.sin(data[:,1] *np.pi/180.0   ))
    affinityList = []
    radiusList = []
    sinThetaAll = np.sin(data[:,1] * np.pi / 180.0)
    boltz    = np.sum(data[:,2] * sinThetaAll * np.exp(-1.0 * data[:,2])) / np.sum(sinThetaAll * np.exp(-1.0 * data[:,2]))
    konList = []
    koffList = []
    for orientation in data:
        theta    = orientation[1]
        phi      = orientation[0]
        energy   = orientation[2]
        sinTheta = np.sin(theta * np.pi / 180.0)
        rotatedCoords = rotatePDB(rawCoords,phi*np.pi/180.0,theta*np.pi/180.0)
        projectedArea = projectOntoSphere(rotatedCoords, npRadius) #gets the projected area of the protein
        effectiveRadius = np.sqrt(projectedArea/np.pi) #the equivalent radius of a circle with the same area as the projection
        effectiveRadius3D = ( -npRadius* effectiveRadius**4 + 2*(npRadius**3) *effectiveRadius*(2*effectiveRadius + np.sqrt(4*npRadius**2 - effectiveRadius**2))    )/( (-2*npRadius**2 + effectiveRadius**2 )**2   ) #the radius of the spherical protein which has a shadow with radius effectiveRadius
        #print effectiveRadius, effectiveRadius3D,SphereProjectedRadius(npRadius,effectiveRadius3D)
        konApprox = projectedArea/( npRadius**2) * avogadroNumber * (npRadius + effectiveRadius3D) * kbT/(6*np.pi*eta) * (1.0/npRadius + 1.0/effectiveRadius3D)
        koffApprox = konApprox * np.exp(energy)
        outputSet.append([proteinData[0], float(proteinData[1])*sinTheta * sinThetaNorm, effectiveRadius3D, konApprox, koffApprox, energy, projectedArea])
        print proteinData[0], float(proteinData[1])*sinTheta * sinThetaNorm, effectiveRadius3D, konApprox, koffApprox, energy, projectedArea
        affinityList.append(float(proteinData[1])*sinTheta * sinThetaNorm * konApprox/koffApprox)
        radiusList.append(effectiveRadius)
        konList.append(float(proteinData[1])*sinTheta * sinThetaNorm*konApprox)
        koffList.append(koffList)
    #np.savetxt("cg_corona_data/"+proteinData[0]+"_"+str(npRadius)+".csv", outputSet)
    #concScaleFactor = 1
    #timeScaleFactor = 1e-9 #amount to rescale kon and koff by 
    #mediumAffinity = np.array(affinityList)*concScaleFactor
    #radiusArr = np.array(radiusList)
    #MFCoverage = mediumAffinity/(1+np.sum(mediumAffinity))
    #totalMFCoverage = np.sum(MFCoverage)
    #if totalMFCoverage > 0.5:
    #    MFCoverage2 = 0.5 * MFCoverage /totalMFCoverage
    #    HSCoverage0 = MFCoverage2
    #else:
    #    HSCoverage0 = MFCoverage
    #print HSCoverage0, np.sum(HSCoverage0), np.amax(HSCoverage0)
    #peakBinding = np.argmax(HSCoverage0)
    #print data[peakBinding]
    #x0Vals = np.sqrt(HSCoverage0)/np.sqrt(1 - HSCoverage0)
    #tvals = np.arange(5)*1e-9
    #initialVals = HSCoverage0
    #oneHourCorona = scint.odeint(DCoverageDTTime, initialVals, tvals,args=( np.array(radiusList) ,np.array(konList),np.array(koffList)), tfirst=True)
    #print oneHourCorona[-1]
    #print np.sum(oneHourCorona[-1])
    #print konApprox
    #HSCoverageX =  scopt.minimize( CoverageWrapperFuncOptX, x0Vals, args=(np.array(radiusList),np.array(konList),np.array(koffList)) )
    #HSCoverage= HSCoverageX.x**2/(1+HSCoverageX.x**2)
    #HSCoverage = HSCoverageX.x
    #print HSCoverage
    #numbindingSites = (4*np.pi*npRadius**2)/(np.pi*radiusArr**2)
    #print proteinData[0], "MF: ", totalMFCoverage, np.sum(MFCoverage*numbindingSites)   ,  "HS: ", np.sum(HSCoverage), np.sum(HSCoverage*numbindingSites)
    #coverage = np.sum(HSCoverage)
    #cgKEqConc = coverage * np.exp(  3*coverage/(1-coverage) + coverage**2/((1-coverage)**2)     )/(1-coverage)
    #print cgKEqConc/(concScaleFactor*float(proteinData[1])), np.exp(-boltz)


np.savetxt("cg_corona_data/"+energyMapFolder+"_"+str(npRadius)+".csv", np.array(outputSet) , fmt="%s")


'''
if orientationAverage!=0:
    for proteinData in concentrationData:
        filename=proteinData[0]+"_"+str(int(npRadius))+"_"+str(int(npZp))+".map"
        (eads,area) = CalculateBoltzArea(energyMapFolder+"/"+filename,  pdbFolder+"/"+proteinData[0]+".pdb")
        proteinRadius = np.sqrt(area/np.pi)
        #calculate kon from smulchowski collision rate
        #kon exp(eads) =   koff 
        konApprox = 466874 *  SphereProjectedRadius(npRadius,proteinRadius)**2 * (proteinRadius+npRadius)**2/( proteinRadius * npRadius**3)   
        koffApprox = konApprox * np.exp(eads)
        print [float(proteinData[1]), proteinRadius, konApprox, koffApprox, eads, np.pi*proteinRadius**2]
        proteinDataList.append([float(proteinData[1]), proteinRadius, konApprox, koffApprox, eads, np.pi*proteinRadius**2])
else:
    for proteinData in concentrationData:
        filename=proteinData[0]+"_"+str(int(npRadius))+"_"+str(int(npZp))+".map"
        #load in the file as before, but calculate the projected area,kon and koff separately for each orientation
        rawCoords =  getAtomCoords( pdbFolder+"/"+proteinData[0]+".pdb")*0.1
        data     = np.genfromtxt(energyMapFolder+"/"+filename)
        sinThetaNorm = 1.0/np.sum( np.sin(data[:,1] *np.pi/180.0   ))
        for orientation in data:
            theta    = orientation[1]
            phi      = orientation[0]
            energy   = orientation[2]
            sinTheta = np.sin(theta * np.pi / 180.0)
            rotatedCoords = rotatePDB(rawCoords,phi*np.pi/180.0,theta*np.pi/180.0)
            projectedArea = projectOntoSphere(rotatedCoords, npRadius)
            effectiveRadius = np.sqrt(projectedArea/np.pi)
            konApprox = projectedArea/( npRadius**2) * avogadroNumber * (npRadius + effectiveRadius) * kbT/(6*np.pi*eta) * (1.0/npRadius + 1.0/effectiveRadius)
            koffApprox = konApprox * np.exp(energy)
            print proteinData[0], float(proteinData[1])*sinTheta * sinThetaNorm, effectiveRadius, konApprox, koffApprox, energy, projectedArea
'''


#print np.array(proteinDataList)

'''
proteinDataArray = np.array(proteinDataList)
#print proteinDataArray
radiusList = proteinDataArray[:,1]
numbindingSites = (4*np.pi*npRadius**2)/(np.pi*radiusList**2)


#medium affinity is the medium-specific affinity with which a protein binds to an NP, given by concentration*kads/kdads = concentration*exp(-ebind)
mediumAffinity = proteinDataArray[:,0]



#the coverages in the mean-field approximation are given by the analytical expression below.
MFCoverage = mediumAffinity/(1+np.sum(mediumAffinity))
x0Vals = np.sqrt(MFCoverage)/(np.sqrt(1 - MFCoverage))
'''
'''
root-finding doesn't really succeed because the scipy routine isn't great.

#non-linear root finding to obtain the HS coverages, using the rescaled variables x such that s = x**2/(1+x**2) and the MF coverages as the initial estimate.
HSCoverageX = scopt.root(  CoverageWrapperFunc, x0Vals, args=(proteinDataOriginal[:,1],mediumAffinity) )
print HSCoverageX
print HSCoverageX.x**2/(1+HSCoverageX.x**2)

thus, we instead employ a minimisation routine. we minimise the sum-of-squares of the rates calculated using the unbounded variable x where s=x**2/(1+x**2).
this isn't bijective but works well enough.

HSCoverageX =  scopt.minimize( CoverageWrapperFuncOpt, x0Vals, args=(radiusList,mediumAffinity) ,tol=1e-12)
HSCoverage= HSCoverageX.x**2/(1+HSCoverageX.x**2)

HSNumBound = HSCoverage*numbindingSites
#array ranking from stack overflow
HSorder = np.flip(HSNumBound.argsort())
HSNRanks = HSorder.argsort()


reslist = []
print "protein-ID, R/nm, ebind/kbT, conc,  MF-coverage, MF-num, HS-coverage, HS-num"
for i in range(len(proteinDataArray)):
    print proteinIDList[i], radiusList[i], -np.log(mediumAffinity[i]/float(concentrationData[i,1])) , float(concentrationData[i,1]),MFCoverage[i], MFCoverage[i]*numbindingSites[i], HSCoverage[i], HSCoverage[i]*numbindingSites[i],HSNRanks[i]+1
    #reslist.append(" ".join( map( str,   [proteinIDList[i],radiusList[i], -np.log(mediumAffinity[i]/float(concentrationData[i,1])) , float(concentrationData[i,1]),MFCoverage[i], MFCoverage[i]*numbindingSites[i], HSCoverage[i], HSCoverage[i]*numbindingSites[i]]    )  )     )
    reslist.append(   [proteinIDList[i],radiusList[i], -np.log(mediumAffinity[i]/float(concentrationData[i,1])) , float(concentrationData[i,1]),MFCoverage[i], MFCoverage[i]*numbindingSites[i], HSCoverage[i], HSCoverage[i]*numbindingSites[i],HSNRanks[i]+1]     )

print np.array(reslist, dtype=np.str)
np.savetxt("rutile_corona_2020/corona_R"+str(int(npRadius))+"_ZP_"+str(int(npZp))+".txt", np.array(reslist, dtype=np.str),  header="protein-ID R/nm ebind/kbT conc  MF-coverage MF-num HS-coverage HS-num HS-num-rank", fmt="%s")
'''

