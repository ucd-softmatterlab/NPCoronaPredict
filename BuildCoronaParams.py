import numpy as np
import scipy.optimize as scopt
import scipy.spatial as scispat
import scipy.integrate as scint
import scipy.special as scspec
import os

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

def readLigandFile(filepath):
    ligandSet = {}
    fileIn = open(filepath,"r")
    seenLigands = []
    lastLigand = ""
    for line in fileIn:
        if line[0] == "#": 
            continue
        lineTerms = line.strip().split(",")
        if len(lineTerms) < 2:
            continue
        ligandID = lineTerms[0]
        ligandType = lineTerms[0].split("-")[0]
        targetAA = lineTerms[1]
        if ligandType != lastLigand:
            #ligand type has changed so register the last one
            if ligandType not in seenLigands:
                seenLigands.append(lastLigand)
            lastLigand = ligandType
        if ligandType not in seenLigands: 
            #this ligand type has not been completed so we can append a new version
            if ligandID not in ligandSet.keys():
                ligandSet[ligandID] = targetAA
            else:
                print("mapping for ", ligandID, " already set, ignoring")
        else:
            print("ligand type", ligandType, "already completed, ignoring")
    fileIn.close()
    print("Ligand overrides: ", ligandSet)
    return ligandSet

#ligandData=readLigandFile("ligand-data/testliganddef2.csv")
#print(ligandData)
def getAtomCoords(filename, caOnly=False,ligandFile=""):
    fileIn = open(filename,"r")
    coordList = []
    moleculeBeadLabels = []
    moleculeBeadRadiusSet = []
    print("looking for ligand file", ligandFile)
    if ligandFile != "":
        ligandLookup = readLigandFile(ligandFile)
    else:
        ligandLookup = {}
    #ligandLookup = {"NAG-C2":"DGL" }
    for line in fileIn:
        lineData = line.split()
        if lineData[0] == "ATOM" and (  lineData[13:15]=="CA" or caOnly == False  ):
            coordList.append([ float(line[30:38]) , float(line[38:46]) , float(line[46:54])])   #x.emplace_back(0.1 * std::stod(line.substr(30, 8)));
            moleculeBeadLabelIn =  line[17:20]
            moleculeBeadLabels.append( moleculeBeadLabelIn ) # 17, 3
            moleculeBeadRadiusSet.append( beadRadiusSet.get(  moleculeBeadLabelIn, 0.5 ) ) 
        elif line[:6] == "HETATM":
            trialLigandID = line[17:20].strip()+"-"+line[12:16].strip()
            #print(trialLigandID)
            if trialLigandID in ligandLookup.keys():
                print("mapping ", trialLigandID, "to", ligandLookup[trialLigandID])
                ligandBeadName = ligandLookup[trialLigandID]
                coordList.append([ float(line[30:38]) , float(line[38:46]) , float(line[46:54])])   #x.emplace_back(0.1 * std::stod(line.substr(30, 8)));
                moleculeBeadLabelIn =  ligandBeadName
                moleculeBeadLabels.append( moleculeBeadLabelIn ) # 17, 3
                moleculeBeadRadiusSet.append( beadRadiusSet.get(  moleculeBeadLabelIn, 0.5 ) )

    fileIn.close()
    # return np.array(coordList)
    coordData = np.atleast_2d( np.array(coordList) )
    print(filename)
    print(coordData)
    coordData[:,0] = coordData[:,0] - np.mean(coordData[:,0])
    coordData[:,1] = coordData[:,1] - np.mean(coordData[:,1])
    coordData[:,2] = coordData[:,2] - np.mean(coordData[:,2])
    return coordData, moleculeBeadLabels, np.array(moleculeBeadRadiusSet)

def pointsToBeads(coords,beadRadius):
    if(len(coords)>2):
        projectedConvexHull = scispat.ConvexHull( coords[:,(0,1) ] )
        edgeCoords = coords[projectedConvexHull.vertices ]
    else:
        edgeCoords = coords[:,(0,1)]
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
    phiRotated =  -phiVal
    thetaRotated = np.pi - thetaVal
    rotCoords = np.copy(coords)
    rotCoords[:,0] = coords[:,0] * np.cos(phiRotated) - coords[:,1] * np.sin(phiRotated)
    rotCoords[:,1] = coords[:,0] * np.sin(phiRotated) + coords[:,1] * np.cos(phiRotated)
    finalCoords = np.copy(rotCoords)
    finalCoords[:,0] = rotCoords[:,0] * np.cos(thetaRotated) + rotCoords[:,2] * np.sin(thetaRotated) 
    finalCoords[:,2] = -1.0 * rotCoords[:,0] * np.sin(thetaRotated) + rotCoords[:,2] * np.cos(thetaRotated)
    return finalCoords


def rotatePDB3Old(coords,phiVal,thetaVal,omegaVal):
    phiRotated =  -phiVal
    thetaRotated = np.pi - thetaVal
    omegaRotated = omegaVal
    rotCoords = np.copy(coords)
    rotCoords[:,0] = coords[:,0] * np.cos(phiRotated) - coords[:,1] * np.sin(phiRotated)
    rotCoords[:,1] = coords[:,0] * np.sin(phiRotated) + coords[:,1] * np.cos(phiRotated)
    finalCoords = np.copy(rotCoords)
    finalCoords[:,0] = rotCoords[:,0] * np.cos(thetaRotated) + rotCoords[:,2] * np.sin(thetaRotated)
    finalCoords[:,2] = -1.0 * rotCoords[:,0] * np.sin(thetaRotated) + rotCoords[:,2] * np.cos(thetaRotated)

    finalCoords2 = np.copy(finalCoords)
    finalCoords2[:,0] = finalCoords[:,0] * np.cos(omegaRotated) - finalCoords[:,1] * np.sin(omegaRotated)
    finalCoords2[:,1] = finalCoords[:,0] * np.sin(omegaRotated) + finalCoords[:,1] * np.cos(omegaRotated)


    return finalCoords2

def rotatePDB3(coords,phiVal,thetaVal,omegaVal):
    phi =  -phiVal
    theta = np.pi - thetaVal
    omega = omegaVal
    rotCoords = np.copy(coords)
    #rotCoords[:,0] = coords[:,0] * np.cos(phiRotated) - coords[:,1] * np.sin(phiRotated)
    #rotCoords[:,1] = coords[:,0] * np.sin(phiRotated) + coords[:,1] * np.cos(phiRotated)
    #finalCoords = np.copy(rotCoords)
    #finalCoords[:,0] = rotCoords[:,0] * np.cos(thetaRotated) + rotCoords[:,2] * np.sin(thetaRotated)
    #finalCoords[:,2] = -1.0 * rotCoords[:,0] * np.sin(thetaRotated) + rotCoords[:,2] * np.cos(thetaRotated)

    #finalCoords2 = np.copy(finalCoords)
    #finalCoords2[:,0] = finalCoords[:,0] * np.cos(omegaRotated) - finalCoords[:,1] * np.sin(phiRotated)
    #finalCoords2[:,1] = finalCoords[:,0] * np.sin(phiRotated) + finalCoords[:,1] * np.cos(phiRotated)
    #rxx = np.cos(theta)*np.sin(phi)
    rxx = np.cos(theta) * np.cos(phi) * np.cos(omega) - np.sin(omega) * np.sin(phi);
    rxy = -1.0 * np.cos(theta) * np.sin(phi) * np.cos(omega) - np.cos(phi)*np.sin(omega);
    rxz = np.sin(theta)*np.cos(omega);
    ryx = np.sin(phi)*np.cos(omega) + np.cos(phi) * np.cos(theta)*np.sin(omega);
    ryy = np.cos(phi)*np.cos(omega) - np.cos(theta)*np.sin(omega)*np.sin(phi);
    ryz = np.sin(omega) * np.sin(theta);
    rzx = -1.0 * np.sin(theta) * np.cos(phi);
    rzy = np.sin(theta) * np.sin(phi);
    rzz = np.cos(theta);
    rotCoords[:,0] = rxx*coords[:,0] + rxy*coords[:,1] + rxz*coords[:,2] 
    rotCoords[:,1] = ryx*coords[:,0] + ryy*coords[:,1] + ryz*coords[:,2]
    rotCoords[:,2] = rzx*coords[:,0] + rzy*coords[:,1] + rzz*coords[:,2]

    return rotCoords



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



def projectOntoSphere(coords,npRadius,beadRadius, ccd=-1):
    coords[:,0] = coords[:,0] - np.mean(coords[:,0])
    coords[:,1] = coords[:,1] - np.mean(coords[:,1])
    zOffset = - np.amin(coords[:,2]) + npRadius + beadRadius[  np.argmin(coords[:,2])  ]
    if ccd > -0.1:
        zOffset = - np.mean(coords[:,2] ) + ccd
    coords[:,2] = coords[:,2] + zOffset #after this transformation the proteins centre is at x=0,y=0 and z such that the protein is just touching the NP
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

def projectOntoCylinder(coords,npRadius,beadRadius = 0.5,ccd=-1):
    coords[:,0] = coords[:,0] - np.mean(coords[:,0])
    coords[:,1] = coords[:,1] - np.mean(coords[:,1])
    zOffset = - np.amin(coords[:,2]) + npRadius + beadRadius[  np.argmin(coords[:,2])  ]
    if ccd > -0.1:
        zOffset = - np.mean(coords[:,2] ) + ccd
    #print(zOffset)
    coords[:,2] = coords[:,2] + zOffset 
    angleCoord = np.arctan2(coords[:,2] , coords[:,1])
    #print( beadRadius/ np.sqrt(1e-6 +  coords[:,1]**2 + coords[:,2]**2) )
    beadOffsetAngle = np.arcsin( beadRadius/ np.sqrt( beadRadius**2 + coords[:,1]**2 + coords[:,2]**2))
    #print(beadOffsetAngle)
    if np.any( np.isnan(beadOffsetAngle) ):
        print( beadOffsetAngle )
        for i in range(len(coords)):
            print(ccd, zOffset, coords[i] , 1.0/np.sqrt( beadRadius**2 + coords[i,1]**2 + coords[i,2]**2 ) , np.arcsin( beadRadius[i]/ np.sqrt( beadRadius**2 + coords[i,1]**2 + coords[i,2]**2)) )
        quit
    beadSet1 =   np.column_stack(( coords[:,0] + beadRadius, angleCoord + beadOffsetAngle   ))
    beadSet2 =   np.column_stack( (coords[:,0] - beadRadius, angleCoord + beadOffsetAngle   ))
    beadSet3 =   np.column_stack( (coords[:,0] + beadRadius, angleCoord - beadOffsetAngle   ))
    beadSet4 =   np.column_stack( (coords[:,0] - beadRadius, angleCoord - beadOffsetAngle   ))
    allBeads = np.row_stack( (beadSet1,beadSet2,beadSet3,beadSet4))
    if len(allBeads) < 20:
        return bindingAreaCylinder(npRadius, np.amax(  np.abs( allBeads[:,0])) ) #for a very small protein we just look at the most extreme point on the x-axis and calculate the area as if it were spherical
    projectionHull = scispat.ConvexHull(allBeads)
    return npRadius * projectionHull.volume

    
    
def bindingAreaCylinder(rnp,ri):
    return ri*rnp* 4 * np.sqrt(  rnp*(2 + rnp/ri)/ri    ) * (  scspec.ellipe(-1.0/( 2*rnp/ri + rnp*rnp/(ri*ri) )) - scspec.ellipk(-1.0/( 2*rnp/ri + rnp*rnp/(ri*ri))  )  )

def cylinderRadiusWrapperFunc(ri,area,rnp):
    return (area - bindingAreaCylinder(rnp,np.sqrt(ri**2)))**2
    
    
    
parser = argparse.ArgumentParser(description = "Parameters for corona calculation")
parser.add_argument('-r','--radius',type=float,help="Radius of the NP [nm]",default=19)
parser.add_argument('-z','--zeta',type=float,help="Zeta potential of the NP [mV]",default=0)
parser.add_argument('-a','--average',type=float,help="average over orientations (does nothing right now)",default=0)
parser.add_argument('-f','--folder',type=str,help="folder containing UA heatmaps",default="results_anatase_alltargets_sphere")
parser.add_argument('-s','--shape',type=int,help="NP shape: 1 = sphere 2 = cylinder 3 = plane", default=1)
parser.add_argument('-p','--proteins',type=str,help="protein definition file",default="ProteinConcs.csv")
parser.add_argument('-c','--coordfolder',type=str,help="location of PDB files",default="pdbs/All")
parser.add_argument('-v','--verbose',type=int,help="verbose",default=0)
parser.add_argument('-n','--nullfile',type=int,default=0,help="Write out the null corona input (planar projection, zero binding energy)")
parser.add_argument('-b','--beadset',type=str,help="Bead definition file", default="beadsets/StandardAABeadSet.csv" )
parser.add_argument('-o','--outputname',type=str,help="Name for the output file", default="buildcoronaparams-demo.csv")
parser.add_argument('-I','--inneroverride',help="For custom .np files with a manually specified inner bound, set this value to the inner bound so BCP can find the correct .uam files, else do not use", default=-1,type=int)
parser.add_argument("-L","--ligand-file", type=str, default = "", help = "Path to a UA ligand override file, leave blank to skip")


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


concentrationData = np.genfromtxt(args.proteins,delimiter=",",dtype=str,skip_header=1)
#print concentrationData
#if len(concentrationData) < 2:
#    concentrationData =np.array([concentrationData])
if concentrationData.ndim == 1:
    concentrationData = np.array([concentrationData])


#print concentrationData
#define where to look for the required input: a pdb file and the output from united atom
pdbFolder = args.coordfolder
#energyMapFolder = "results_swcnt_alltargets"
energyMapFolder = args.folder
#parameters for calculating rate constants
avogadroNumber = 6.022e23
kbT = 300 * 1.381e-23 
eta = 8.9e-4 #viscosity of the medium

npRadius = args.radius
npZp = args.zeta
orientationAverage = args.average

npShape = args.shape
if npShape !=1 and npShape !=2 and npShape !=3:
    npShape = 1

doNullCorona = False
if args.nullfile != 0:
    doNullCorona = True
    outputSetNull = []
    nullCoronaDenom = 1
proteinIDList = []
proteinDataList = []

diffusionCoeff = 1e-11
outputSet = []
averageValSet = []

#load in radii for single-bead models
beadRadiusSet = {}
try:
    beadRadiusFile = open(args.beadset,"r")
    for line in beadRadiusFile:
        if line[0]=="#":
            continue
        lineTerms = line.strip().split(",")
        beadRadiusSet[ lineTerms[0] ] = float(lineTerms[2])
    beadRadiusFile.close()
except:
    print("Bead definition file could not be found, using default radius for all beads")



filenameRadius = str(int(npRadius))
if args.inneroverride > 0:
    filenameRadius = str(int( args.inneroverride ))
    print("Looking for .uam files with radius ", filenameRadius, " in the name.")


for proteinData in concentrationData:
    #print(proteinData[0] )
    if npShape==1 or npShape == 3:
        filename=proteinData[0]+"_"+filenameRadius+"_"+str(int(npZp))+".uam"
        print("looking for: ", filename)
        #load in the file as before, but calculate the projected area,kon and koff separately for each orientation
        #rawCoords =  getAtomCoords( pdbFolder+"/"+proteinData[0]+".pdb")*0.1
        data     = np.genfromtxt(energyMapFolder+"/"+filename)
        omegaSet = np.reshape( np.zeros_like(data[:,0]), (-1,1) )
        data = np.concatenate( (data,omegaSet) , axis=1)
    else:
        alldatafiles = []
        for omega in [0,45,90,135]:
            filename = proteinData[0]+"_"+filenameRadius+"_"+str(int(npZp))+"_"+str(omega)+".uam"
            print("Looking for: ", filename)
            newdata = np.genfromtxt(energyMapFolder+"/"+filename)
            omegaSet = np.reshape( np.zeros_like(newdata[:,0]) + omega, (-1,1))
            #print(omegaSet)
            newdata = np.concatenate( (newdata, omegaSet) ,axis=1)
            alldatafiles.append(newdata)
        data = np.concatenate(alldatafiles,axis=0)
    #print(data[0:2])
 
    singleBeadMode = False
    print("input ligand file", args.ligand_file)
    rawCoords, beadLabels, moleculeBeadRadiusSet =  getAtomCoords( pdbFolder+"/"+proteinData[0]+".pdb", ligandFile=args.ligand_file)
    rawCoords = rawCoords*0.1




    thetaOffset = (np.amax(data[1:,1] - data[0:-1,1]))/2.0 #UA output gives the left-edge of the bin containing the angular orientations, this corrects for that
    phiOffset = (np.amax(data[1:,0] - data[0:-1,0]))/2.0 #UA output gives the left-edge of the bin containing the angular orientations, this corrects for that
    #print thetaOffset
    data[:,1] +=thetaOffset
    data[:,0] += phiOffset
    
    
    if len(rawCoords)  == 1:
        print("Using single bead mode")
        singleBeadMode = True
        data = np.array( [  data[0] ] ) #a single bead has no orientation dependency, so discard everything but the first entry.
        data[0,1] = 90.0 * np.pi / 180.0 #manually set the single bead theta value to 90 degrees to avoid the inevitable confusing division by zero
    sinThetaNorm = 1.0/np.sum( np.sin(data[:,1] *np.pi/180.0   ))
    affinityList = []
    radiusList = []
    sinThetaAll = np.sin(data[:,1] * np.pi / 180.0)
    boltz    = np.sum(data[:,2] * sinThetaAll * np.exp(-1.0 * data[:,2])) / np.sum(sinThetaAll * np.exp(-1.0 * data[:,2]))
    konList = []
    koffList = []
    #nameList = []
    for orientation in data:
        theta    = orientation[1]
        phi      = orientation[0]
        energy   = orientation[2]
        omega = orientation[-1]
        ccd = -1.0
        if len(orientation) > 10:
            ccd = orientation[10]
        #nameList.append( proteinData[0]+":"+str(theta)+"-"+str(phi) )
        sinTheta = np.sin(theta * np.pi / 180.0)
        rotatedCoords = rotatePDB3(rawCoords,phi*np.pi/180.0,theta*np.pi/180.0, omega*np.pi/180.0)
        if npShape == 1:
            projectedArea = projectOntoSphere(rotatedCoords, npRadius,moleculeBeadRadiusSet,ccd) #gets the projected area of the protein
            effectiveRadius = np.sqrt(projectedArea/np.pi) #the equivalent radius of a circle with the same area as the projection
            effectiveRadius3D = ( -npRadius* effectiveRadius**4 + 2*(npRadius**3) *effectiveRadius*(2*effectiveRadius + np.sqrt(4*npRadius**2 - effectiveRadius**2))    )/( (-2*npRadius**2 + effectiveRadius**2 )**2   ) #the radius of the sph$
        elif npShape == 2:
            #print(theta,phi,omega)
            projectedArea = projectOntoCylinder(rotatedCoords,npRadius,moleculeBeadRadiusSet,ccd)
            radiusApprox = np.sqrt(projectedArea/np.pi)
            #print scopt.root(  cylinderRadiusWrapperFunc, radiusApprox, args=(projectedArea,npRadius) )
            optRes= scopt.minimize(  cylinderRadiusWrapperFunc, np.sqrt( projectedArea/np.pi), args=(projectedArea, npRadius),tol=0.01 )
            effectiveRadius3D= np.abs( (optRes.x)[0])
            if effectiveRadius3D > 5*radiusApprox:
                effectiveRadius3D = 5*radiusApprox #cap the effective radius to the value computed from the projected area 
            projectedAreaFromCalc = bindingAreaCylinder(npRadius,effectiveRadius3D)
            #print("First approx: ", radiusApprox, " second approx: ", effectiveRadius3D, "projected area:" , projectedArea, " re-projected area: ", projectedAreaFromCalc," num binding", (4*np.pi*npRadius**2)/projectedArea)
            #effectiveRadius = np.sqrt(projectedArea/np.pi) #figure out how to do this for the cylindrical projection!
        elif npShape == 3:
            projectedArea = projectOntoSphere(rotatedCoords,1000,moleculeBeadRadiusSet,ccd)
            effectiveRadius = np.sqrt(projectedArea/np.pi)
            effectiveRadius3D = effectiveRadius 
        else:
            projectedArea = projectOntoSphere(rotatedCoords,npRadius,moleculeBeadRadiusSet,ccd)
            effectiveRadius = np.sqrt(projectedArea/np.pi) #the equivalent radius of a circle with the same area as the projection
            effectiveRadius3D = ( -npRadius* effectiveRadius**4 + 2*(npRadius**3) *effectiveRadius*(2*effectiveRadius + np.sqrt(4*npRadius**2 - effectiveRadius**2))    )/( (-2*npRadius**2 + effectiveRadius**2 )**2   ) #the radius of the spherical protein which has a shadow with radius effectiveRadius
        #print effectiveRadius, effectiveRadius3D,SphereProjectedRadius(npRadius,effectiveRadius3D)
        pairDiffusionCoeff = kbT/(6*np.pi*eta) * (1.0/ (1e-9*npRadius) + 1.0/(1e-9*effectiveRadius3D))
        #print "DCoeff: ", pairDiffusionCoeff, "pair radius", (npRadius + effectiveRadius3D), "[nm]"
        pairCollisionRate = 4 * np.pi *(npRadius*1e-9 + effectiveRadius3D*1e-9) * pairDiffusionCoeff * avogadroNumber
        numBindingSites = (4*np.pi*npRadius**2)/projectedArea
        #print "Pair collision rate (total per NP):", pairCollisionRate
        #print "Num binding sites: ", numBindingSites
        konApprox =1000* pairCollisionRate/numBindingSites #prefactor is to go from m^3 / mol to 1/M
        #konApprox =  projectedArea/( npRadius**2) * avogadroNumber * (npRadius + effectiveRadius3D) * kbT/(6*np.pi*eta) * (1.0/npRadius + 1.0/effectiveRadius3D) #in dm^3/mol / s 
        koffApprox = konApprox * np.exp(energy)
        outputSet.append([proteinData[0]+":"+str(theta)+"-"+str(phi), float(proteinData[1])*sinTheta * sinThetaNorm, effectiveRadius3D, konApprox, koffApprox, energy, projectedArea])
        if args.verbose == 1:
            print(proteinData[0]+":"+str(theta)+"-"+str(phi), float(proteinData[1])*sinTheta * sinThetaNorm, effectiveRadius3D, konApprox, koffApprox, energy, projectedArea)
        if doNullCorona == True:
            projectedAreaPlane = projectOntoSphere(rotatedCoords,npRadius+500.0)
            effectiveRadiusPlane = np.sqrt(projectedAreaPlane/np.pi)
            outputSetNull.append([proteinData[0], proteinData[0]+":"+str(theta)+"-"+str(phi), float(proteinData[1])*sinTheta * sinThetaNorm, effectiveRadiusPlane, konApprox, konApprox, 0, projectedAreaPlane])
            
            
outputPathFull = args.outputname
outputBasePath = "/".join(  outputPathFull.split("/")[:-1]  )
os.makedirs(outputBasePath,exist_ok=True)
#"/".join(energyMapFolder.split("/")[:-1])
np.savetxt(outputPathFull , np.array(outputSet) , fmt="%s")
if doNullCorona == True:
    nullConcArray = np.array(outputSetNull)
    np.savetxt("cg_corona_data/"+energyMapFolder+"_nullcorona.csv", nullConcArray , fmt="%s")
    allConcs = nullConcArray[:,2].astype(float)
    mfDenom = 1.0 + np.sum(allConcs)
    numBindingSites = 1.0 / ( nullConcArray[:,7].astype(float)) #num binding sites per nm^2 implicitly
    mfProteinNums = allConcs*numBindingSites/mfDenom
    uniqueProteins = np.unique(nullConcArray[:,0])
    proteinNumberDict = {}
    for up in uniqueProteins:
        proteinNumberDict[up] = 0
    for i in range(len(nullConcArray[:,0])):
        proteinNumberDict[nullConcArray[i,0] ] += mfProteinNums[i]
    print(proteinNumberDict)
    
    #write out the steady-state mean-field corona prediction 
#np.savetxt("cg_corona_data/"+energyMapFolder+"_"+str(int(npRadius))+"_"+str(int(npZp))+".csv", np.array(outputSet) , fmt="%s")
