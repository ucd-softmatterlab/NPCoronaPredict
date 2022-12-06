'''
Rotate a protein from the initial co-ordinates to principal-axis-aligned form, with the longest axis along Z and the second longest along Y
The final rotation convention is set such that the x-component of the dipole moment is negative
This has the benefit that the UA rotation angles of theta = 90 aligns the longest part of the protein with the surface of the NP
'''



from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO
import numpy.random as npr
import Bio.PDB.Superimposer as Superimposer
import warnings
from Bio import BiopythonWarning


import os
import numpy as np
import numpy.linalg as npla

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
sup = Superimposer()
io = PDBIO()


targetFolders = [ "pdb"]

#targetFolders = ["pdbs/ExtraStructures"]
allTargets0  = []


#define the per-AA properties used to generate property-dipoles for defining the handedness of proteins
partialChargeDict={ "ARG":1, "HIS":1, "LYS":1, "ASP":-1, "GLU":-1 }
hydrophobicityDict = {"ALA":0}


for targetFolder in targetFolders:
    targetsFound = SearchDirectory(targetFolder)
    for targetFound in targetsFound:
        allTargets0.append(targetFound)
allTargets = sorted(  allTargets0   )

results = []


xrot180 = np.array( [ [1,0,0],[0,-1,0],[0,0,-1] ])
yrot180 = np.array( [ [-1,0,0],[0,1,0],[0,0,-1] ])
zrot180 = np.array( [ [-1,0,0], [0,-1,0],[0,0,1]])


for targetPath in allTargets:
    coordSet = []
    pathTerms = targetPath.split("/")
    target = pathTerms[-1]
    try:
        inputStructure = parser.get_structure("inputstructure",targetPath)
    except:
        print("Error loading ",target)
        continue
    chargeList = []
    for atom in inputStructure.get_atoms():
        if atom.get_name() != "CA":
            continue
        residueName = (atom.get_parent()).get_resname()
        coords = atom.get_coord()
        chargeList.append( partialChargeDict.get(residueName,0.0) )
        coordSet.append( coords )
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
    for atom in inputStructure.get_atoms():
        oldCoords = np.reshape(atom.get_coord() ,(1,3)   )
        newCoords = np.matmul( transformMatrix, oldCoords.T).T
        atom.set_coord(newCoords[0])
    io.set_structure(inputStructure)
    outputPath = "rotated_pdbs/"+"/".join( pathTerms[:-1])
    os.makedirs(outputPath,exist_ok = True)
    io.save(outputPath+"/"+target)
    print("Saved to: ", outputPath+"/"+target)
