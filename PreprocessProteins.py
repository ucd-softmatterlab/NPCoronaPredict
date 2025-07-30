from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.Residue  import DisorderedResidue
import Bio.PDB.Superimposer as Superimposer
import warnings
from Bio import BiopythonWarning


import numpy as np
import copy
import os
import numpy.random as npr
import numpy.linalg as npla
import argparse



warnings.simplefilter('ignore',BiopythonWarning)


def SearchDirectory(path):
    files = os.listdir(path)
    pdbs  = []
    for handle in files:
        abspath = os.path.join(path, handle)
        if os.path.isdir(abspath):
            pdbs += SearchDirectory(abspath)
        elif abspath[-4:] == ".pdb":
            pdbs.append(abspath)
    return pdbs

def listDiff(mainSet,subSet):
    subSet = set(subSet)
    return [item for item in mainSet if item not in subSet]

def zeroCOM(coordArr):
    coords = np.copy(coordArr)
    coords[:,0] = coords[:,0] - np.mean(coords[:,0])
    coords[:,1] = coords[:,1] - np.mean(coords[:,1])
    coords[:,2] = coords[:,2] - np.mean(coords[:,2])
    return coords

def getTransformMatrix(coordArr):
    ixx = np.sum(coordArr[:,1]**2 + coordArr[:,2]**2)
    ixy = np.sum(- coordArr[:,0]*coordArr[:,1])
    ixz = np.sum(-coordArr[:,0]*coordArr[:,2])
    iyy = np.sum(coordArr[:,0]**2 + coordArr[:,2]**2)
    iyx = ixy
    iyz = np.sum(-coordArr[:,1]*coordArr[:,2])
    izz = np.sum(coordArr[:,0]**2 + coordArr[:,1]**2)
    izx = ixz
    izy = iyz
    inertialArray = np.array([ [ixx, ixy,ixz],[iyx,iyy,iyz],[izx,izy,izz] ])
    eigvals,eigvecs = npla.eig(inertialArray)
    #invarr = npla.inv(eigvecs)
    sortIndex = eigvals.argsort()[::-1]
    invarr = npla.inv(eigvecs[:,sortIndex])
    return invarr



parser = PDBParser(PERMISSIVE=1)
io = PDBIO()
sup = Superimposer()


argparser = argparse.ArgumentParser(description="Parameters for PreprocessProteins")
argparser.add_argument("-p","--ph", type=float,default=7,help="Target pH")

#argparser.add_argument("-P","--propka", type=int,default=1,help="If non-zero run propka")
#argparser.add_argument("-r","--rotate", type=int, default=1,help="Rotate proteins to canonical form if non-zero (default yes)")
argparser.add_argument("-P", "--propka", action="store_true", help="Enable running PROPKA")
argparser.add_argument("-r", "--rotate", action="store_true", help="Rotate proteins to canonical form")


argparser.add_argument("-f","--folder", type=str, default="all_proteins",help="Target protein folder for conversion")
argparser.add_argument("-o","--outputfolder",type=str,default="",help="Output folder, if left blank will be auto-generated")
#argparser.add_argument("-i","--interpolation",type=int,default=1,help="If non-zero residues are interpolated, if zero set to most likely form")
argparser.add_argument("-L","--likelyonly", action="store_true", help="Force residues to only take the most likely form")

argparser.add_argument("--hid", type=float,default = 0.2, help="Fraction of HID relative to HIE")
argparser.add_argument("--aaset", help="path to a .csv file containing standard/protonated/deprotonated bead names to override default", default="" )
args = argparser.parse_args()


makeCanonical = True
canonicalString = "-canonical"
if args.rotate == False:
    makeCanonical = False
    canonicalString = ""

targetPH = args.ph

if args.outputfolder=="":
    outputFolder = args.folder+"_pH"+str(targetPH)+canonicalString
else:
    outputFolder = args.outputfolder
    
    
runPropka = True
if args.propka == False:
    runPropka = False
    
os.makedirs(outputFolder,exist_ok=True)
os.makedirs(outputFolder+"/propkaoutput",exist_ok=True)
#if this is set to 0 then residues are set to their most likely form, i.e., protonated if pH < pKa and deprotonated HIS set to HIE. 
#if it is non-zero we include both types of residue and use the occupancy flag to set their relative weight, UA then weights
#the input potentials according to this weight. This results in larger files but smoother behaviour, especially for pH \approx pKa

residueInterpolation = 1
if args.likelyonly == True:
    residueInterpolation = 0

hidFraction = args.hid #HIE (epsilon-protonated) is favoured over HID (delta-protonated) by a ratio of 4:1 in solution





#names for the charge variants
#this table is defined as: standard name, deprotonated name, protonated name
#HIS is a special case, here the deprotonated is HIE and protonated is HIP. HID is added in using the hidFraction later.
aaData = [
["ALA", "ALA", "ALA"],
["ARG", "ARG", "ARG"],
["ASN", "ASN", "ASN"],
["ASP", "ASP", "ASP"],
["CYS", "CYM", "CYS"],
["GLN", "GLN", "GLN"],
["GLU", "GLU", "GAN"],
["HIS", "HIE", "HIP"],
["ILE", "ILE", "ILE"],
["LEU", "LEU", "LEU"],
["LYS", "LYS", "LYS"],
["MET", "MET", "MET"],
["PHE", "PHE", "PHE"],
["SER", "SER", "SER"],
["THR", "THR", "THR"],
["TRP", "TRP", "TRP"],
["TYR", "TYR", "TYR"],
["VAL", "VAL", "VAL"]
]


aaDataExtra =  [
["ALA", "ALA", "ALA"],
["ARG", "ARN", "ARG"],
["ASN", "ASN", "ASN"],
["ASP", "ASP", "AAN"],
["CYS", "CYM", "CYS"],
["GLN", "GLN", "GLN"],
["GLU", "GLU", "GAN"],
["HIS", "HIE", "HIP"],
["ILE", "ILE", "ILE"],
["LEU", "LEU", "LEU"],
["LYS", "LYN", "LYS"],
["MET", "MET", "MET"],
["PHE", "PHE", "PHE"],
["SER", "SER", "SER"],
["THR", "THR", "THR"],
["TRP", "TRP", "TRQ"],
["TYR", "TYM", "TYR"],
["VAL", "VAL", "VAL"]
]

if args.aaset == "":
    print("Using default AA states: GLU/GAN protonation and HIS -> HIE/HID/HIP only")
elif args.aaset == "extended":
    print("Using extended AA states: GLU/GAN, HIS -> HIE/HID/HIP, ASP/AAN, CYS/CYM, TRP/TRQ, TYR/TYM")
    aaData = aaDataExtra
else:
    trialFile = args.aaset
    aaData = []
    if os.path.exists(trialFile):
        fileIn = open(trialFile, "r")
        for line in fileIn:
            if line[0] == "#":
                continue
            lineTerms = line.strip().split(",")
            aaData.append( [lineTerms[0].strip() , lineTerms[1].strip(), lineTerms[2].strip() ] )
        fileIn.close()
        print("Finished reading custom definition from ", trialFile)
    else:
        print("Could not find custom AA definition, quitting")
        quit()

print(aaData)
#quit()

aaDataIndexDict = {}
for i in range(len(aaData)):
    aaDataIndexDict[aaData[i][0]] = i





#define the per-AA properties used to generate property-dipoles for defining the handedness of proteins
partialChargeDict={ "ARG":1, "HIS":1, "LYS":1, "ASP":-1, "GLU":-1 }
hydrophobicityDict = {"ALA":0}

allTargets0  = []
targetsFound = SearchDirectory(args.folder)
for targetFound in targetsFound:
    allTargets0.append(targetFound)
allTargets = sorted(  allTargets0   )


xrot180 = np.array( [ [1,0,0],[0,-1,0],[0,0,-1] ])
yrot180 = np.array( [ [-1,0,0],[0,1,0],[0,0,-1] ])
zrot180 = np.array( [ [-1,0,0], [0,-1,0],[0,0,1]])


missingStates = []


for pdbPath in allTargets:
    pathTerms = pdbPath.split("/")
    fileStem = pathTerms[-1][:-4]
    print("Processing" , pdbPath, "at pH ", targetPH, " file stem: ", fileStem)
    
    
    try:
        structureData = parser.get_structure(fileStem+"-struct" , pdbPath)
    except:
        print("Structure data for ", fileStem, " not found or not readable")
        continue
    if makeCanonical == True:
        coordSet = []
        chargeList = []
        
        
        
        for atom in structureData.get_atoms():
            if atom.get_name() != "CA":
                continue
            residueName = (atom.get_parent()).get_resname()
            coords = atom.get_coord()
            chargeList.append( partialChargeDict.get(residueName,0.0) )
            coordSet.append( coords )
        if len(chargeList) >= 2:
            chargeArr = np.array(chargeList)
            coordArr = zeroCOM(np.array(coordSet))
            transformMatrix = getTransformMatrix(coordArr)
            if npla.det(transformMatrix) < 0: #sometimes this generates an inversion, cancel this by mirroring in the xy plane to ensure we're generating a rotation
                #print(transformMatrix)
                hhXY = np.array( [ [1,0,0],[0,1,0],[0,0,-1] ])
                transformMatrix = np.matmul( hhXY,transformMatrix)
                #print(transformMatrix)
            #the matrix generated so far is a rotation onto the principle axes, but the co-ords are not uniquely defined
            #we use the convention that at most the dipole moment on x can have a negative sign and the other two are positive
            dipoleX = np.dot( chargeArr , coordArr[:,0])
            dipoleY = np.dot(chargeArr, coordArr[:,1])
            dipoleZ = np.dot(chargeArr, coordArr[:,2])
            dipoleArr =  np.array( [  [dipoleX,dipoleY,dipoleZ]           ])
            dipoleArrTransformed = (np.matmul( transformMatrix, dipoleArr.T)).T
            #print(dipoleX,dipoleY,dipoleZ)
            #print(dipoleArr,dipoleArrTransformed)
            if dipoleArrTransformed[0,2] >= 0:
                if dipoleArrTransformed[0,1] < 0:
                    transformMatrix = np.matmul( zrot180, transformMatrix )
            else:
                if dipoleArrTransformed[0,1] <0:
                    transformMatrix = np.matmul( xrot180, transformMatrix)
                else:
                    transformMatrix = np.matmul( yrot180, transformMatrix)
            for atom in structureData.get_atoms():
                oldCoords = np.reshape(atom.get_coord() ,(1,3)   )
                newCoords = np.matmul( transformMatrix, oldCoords.T).T
                atom.set_coord(newCoords[0])        
        
        
    else:
        print("Skipping canonical rotation")
    if runPropka == True:
        if not os.path.isfile(outputFolder+"/propkaoutput/"+fileStem+".pka"):
            print("PROPKA file not found for ", fileStem, "attempting to run propka")
            os.system("propka3 "+pdbPath)
            os.system("mv "+fileStem+".pka "+outputFolder+"/propkaoutput/"+fileStem+".pka")
        try:
            propkaFile = open(outputFolder+"/propkaoutput/"+fileStem+".pka","r")
            propkaFileData = propkaFile.readlines()
            propkaFile.close()
        except:
            print("PROPKA file could not be read for ", fileStem)
            continue
        summaryStart = 0
        newIDVal = 0
        for line in propkaFileData:
            if line == "SUMMARY OF THIS PREDICTION\n":
                 summaryStart = 1
                 continue
            if summaryStart == 1 and line == "--------------------------------------------------------------------------------------------------------\n":
                 break
            if summaryStart == 1:
                 if (line.strip())[:5] != "Group":
                    rowData =  line.strip().split()
                    #rowData: 0 = AA tag, 1 = residue number, 2 = chain ID, 3 = adjusted  pKa val
                    #if pH < pKa the side chain is protonated
                    residuepKa = float(rowData[3])
                    if rowData[0] not in aaDataIndexDict:
                        try:
                            print("Unrecognised residue ", rowData[0])
                        except:
                            print("Other failure: ", rowData)
                        continue
                    for model in structureData:
                         for residue in model.get_residues():
                            if newIDVal < residue.get_id()[1]:
                                newIDVal = residue.get_id()[1]
                         newIDVal = newIDVal + 1
                         deprotonatedFrac = np.exp(targetPH)/( np.exp(targetPH) + np.exp(  residuepKa ))
                         singleResidue = 0
                         if deprotonatedFrac < 0.05 or deprotonatedFrac > 0.95: #to reduce UA run time, if less than 5% is protonated/deprotonated then only use one residue. 
                             singleResidue = 1
                         #if the protonated form = deprotonated form but propka thinks there's an important difference issue a warning - this may mean new PMFs need generation.
                         if aaData[   aaDataIndexDict[rowData[0]]][2] == aaData[   aaDataIndexDict[rowData[0]]][1]:
                             if aaData[   aaDataIndexDict[rowData[0]]][2] not in missingStates:
                                 print("Propka thinks this residue should have multiple forms but no options are set. Please update rowData if PMFs are avaliable")
                                 print(aaData[aaDataIndexDict[rowData[0]]])
                                 missingStates.append(  aaData[   aaDataIndexDict[rowData[0]]][2] )
                             singleResidue = 1
                         if residueInterpolation == 0 or singleResidue == 1:
                              if targetPH <  residuepKa:
                                   newName  = aaData[   aaDataIndexDict[rowData[0]]][2]
                              else:
                                   newName  = aaData[   aaDataIndexDict[rowData[0]]][1]
                              model[rowData[2]][int(rowData[1])].resname = newName
                         else:
                              deprotonatedFrac = np.exp(targetPH)/( np.exp(targetPH) + np.exp(  residuepKa ))
                              protonatedFrac = 1 - deprotonatedFrac
                              targetResidue = model[rowData[2]][int(rowData[1])]
                              hisFactor = 1
                              if rowData[0] == "HIS": #add in the HID form of histidine using the weighting defined earlier
                                   hisFactor = 1-hidFraction
                                   print("HIE: ", deprotonatedFrac*hisFactor, " HID ", deprotonatedFrac*hidFraction, "HIP" , protonatedFrac)
                                   hidCopy = copy.deepcopy(targetResidue)
                                   hidCopy.resname = "HID"
                                   hidCopy.id =  ( hidCopy.get_id()[0],  newIDVal,hidCopy.get_id()[2])
                                   newIDVal += 1
                                   for atom in hidCopy.get_atoms():
                                       atom.set_occupancy(atom.get_occupancy() * deprotonatedFrac * hidFraction   )
                                   model[rowData[2]].add(hidCopy)
                              #set the initially present residue to the deprotonated form
                              #targetResidue = model[rowData[2]][int(rowData[1])]
                              protonatedResidue = copy.deepcopy(targetResidue)
                              targetResidue.resname = aaData[   aaDataIndexDict[rowData[0]]][1]
                              for atom in targetResidue.get_atoms():
                                  atom.set_occupancy(atom.get_occupancy() * deprotonatedFrac *hisFactor  )
                              protonatedResidue.resname = aaData[   aaDataIndexDict[rowData[0]]][2]
                              protonatedResidue.id = ( protonatedResidue.get_id()[0],  newIDVal, protonatedResidue.get_id()[2])
                              newIDVal += 1
                              for atom in protonatedResidue.get_atoms():
                                  atom.set_occupancy(atom.get_occupancy() * protonatedFrac   )
                              model[rowData[2]].add(protonatedResidue)
    io.set_structure(structureData)
    io.save(outputFolder+"/"+fileStem+"-pH"+str(targetPH).replace(".","p")+"model"+str(residueInterpolation)+canonicalString+".pdb")
print("Summary of missing states:")
for missingRes in missingStates:
    print(missingRes)
