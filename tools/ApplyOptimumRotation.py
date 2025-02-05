import numpy as np


from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO
import Bio.PDB.Superimposer as Superimposer
import warnings
from Bio import BiopythonWarning

import os
import argparse

warnings.simplefilter('ignore',BiopythonWarning)



#The rotation convention here is the standard such that it roates by phi around the z axis, rotates by y around theta.
#note that the UA orientation definition is defined such that the orientation labelled "phi, theta" must be rotated by "-phi, pi - theta" to produce the optimal binding orientation relative to an NP "under" the protein.
#if you're ever unsure check using the Dipole.pdb or Tetra.pdb models as these are designed to be very easy to interpret on NPs of large zeta potentials.
def rotatePDB(coords,phi,theta):
    rxx = np.cos(theta)*np.cos(phi)
    rxy = -1 * np.cos(theta)*np.sin(phi)
    rxz = np.sin(theta)
    ryx = np.sin(phi)
    ryy = np.cos(phi)
    ryz = 0
    rzx = -1*np.sin(theta)*np.cos(phi)
    rzy = np.sin(theta)*np.sin(phi)
    rzz = np.cos(theta)
    finalCoords = np.zeros_like(coords)
    finalCoords[:,0] = coords[:,0] * rxx + coords[:,1] * rxy + coords[:,2]*rxz 
    finalCoords[:,1] = coords[:,0] * ryx + coords[:,1] * ryy + coords[:,2]*ryz 
    finalCoords[:,2] = coords[:,0] * rzx + coords[:,1] * rzy + coords[:,2]*rzz 
    return finalCoords

def getRotationMatrix(phi,theta):
    rxx = np.cos(theta)*np.cos(phi)
    rxy = -1 * np.cos(theta)*np.sin(phi)
    rxz = np.sin(theta)
    ryx = np.sin(phi)
    ryy = np.cos(phi)
    ryz = 0
    rzx = -1*np.sin(theta)*np.cos(phi)
    rzy = np.sin(theta)*np.sin(phi)
    rzz = np.cos(theta)
    return [  [rxx,rxy,rxz] , [ryx,ryy,ryz], [rzx,rzy,rzz] ]

parser = argparse.ArgumentParser(description = "Settings for rotating PDB")
parser.add_argument("-p", "--proteinfile", help="Path to initial protein .pdb structure")
parser.add_argument("-u", "--uafile", help="Path to UA heatmap in .uam format")
parser.add_argument("-z", "--zoffset", help="Additional offset to apply along z axis in nm",default=0)
parser.add_argument("-o", "--outputfolder", help="Folder to save output",default="rotated_pdbs") 
parser.add_argument("-R", "--radius", help="NP radius in nm", default=20)
parser.add_argument("-s", "--shiftradius", help="If non-zero, shift protein by NP radius", default = 0)
args = parser.parse_args()


npRadius = 0 #an additional offset to apply
#Tetra.pdb is a test "protein" using the CHL and PHO residues arranged in a tetrahedron with 3 CHL and 1 PHO. CHL carries a +1 charge and PHO a -1 charge, so we expect the orientation with PHO facing away from the NP to occur on NPs with -ve zeta potential and PHO in contact with the NP on NPs of +ve zeta potential
proteinStructurePath = args.proteinfile
proteinHeatmapPath = args.uafile
npRadius = float(args.radius) * 10
zoffset = float(args.zoffset) * 10
outputFolder = args.outputfolder
shiftRadius = 0
if int(args.shiftradius) != 0:
    shiftRadius = 1

parser = PDBParser(PERMISSIVE=1)
sup = Superimposer()
io = PDBIO()


proteinHeatmapData     = np.genfromtxt(proteinHeatmapPath)
proteinHeatmapData[:,0] = proteinHeatmapData[:,0] + 2.5 #UA orientation output is left-hand bin edges, this pushes it back to centre-bin edges
proteinHeatmapData[:,1] = proteinHeatmapData[:,1] + 2.5
optimalBindingLine = proteinHeatmapData[ np.argmin( proteinHeatmapData[:,2]   ) ]

sinthetaVals = np.sin( proteinHeatmapData[:,1] * np.pi / 180.0)
simpleAverage = np.sum( proteinHeatmapData[:,2] * sinthetaVals )/np.sum(   sinthetaVals  )
boltzAverage = np.sum( proteinHeatmapData[:,2] * sinthetaVals *np.exp(-1.0 * proteinHeatmapData[:,2]))/np.sum(   sinthetaVals *np.exp(-1.0 * proteinHeatmapData[:,2]))


print("Protein binding energy stats: ")
print("Minimum energy: ", optimalBindingLine[2])
print("Simple average: ", simpleAverage)
print("Boltzmann average: ", boltzAverage)
print(optimalBindingLine)


minTheta = optimalBindingLine[1]
minPhi = optimalBindingLine[0]

rotationMatrix = np.array( getRotationMatrix( -minPhi * np.pi/180.0 , np.pi - minTheta*np.pi/180.0    ) )

#load in the original structure
proteinStructure = parser.get_structure("originalstructure",proteinStructurePath)
proteinNPID = (proteinHeatmapPath.split("/")[-1] )[:-4]
xcentre=0
ycentre=0
zcentre=0
numCA = 0

#find centre-of-CA mass
for atom in proteinStructure.get_atoms():
    if atom.get_name() != "CA":
        continue
    coords = atom.get_coord()
    xcentre += coords[0]
    ycentre += coords[1]
    zcentre += coords[2]
    numCA += 1
xcentre = xcentre/numCA
ycentre = ycentre/numCA
zcentre = zcentre/numCA

zmin  = 0
#rotate to optimum configuration
for atom in proteinStructure.get_atoms():
    coords0 = atom.get_coord()
    xshift = coords0[0] - xcentre
    yshift = coords0[1] - ycentre
    zshift = coords0[2] - zcentre
    xrot =  rotationMatrix[0,0] * xshift + rotationMatrix[0,1]*yshift + rotationMatrix[0,2] * zshift
    yrot =  rotationMatrix[1,0] * xshift + rotationMatrix[1,1]*yshift + rotationMatrix[1,2] * zshift
    zrot =  rotationMatrix[2,0] * xshift + rotationMatrix[2,1]*yshift + rotationMatrix[2,2] * zshift 
    if zrot < zmin:
        zmin = zrot
    coordsFinal = [xrot,yrot,zrot]
    atom.set_coord(coordsFinal)


beadNPDistanceSet = []
#finally shift along z to the target location
for atom in proteinStructure.get_atoms():
   coords0 = atom.get_coord()
   #print(coords0)
   zshift = coords0[2] - zmin + optimalBindingLine[4]*10  + zoffset + shiftRadius*npRadius
   zshiftNP = coords0[2] - zmin + optimalBindingLine[4]*10  + zoffset +  npRadius
   atom.set_coord( [ coords0[0], coords0[1], zshift])
   if atom.get_name() == "CA":
       atomOffsetDist = np.sqrt(  coords0[0]**2 + coords0[1]**2 + zshiftNP**2  ) - npRadius
       residueID = atom.get_parent()
       
       #print( residueID.__repr__(),residueID.get_resname(), atomOffsetDist)
       beadNPDistanceSet.append( [residueID.get_resname(), 0.1*atomOffsetDist])
beadNPDistanceSet.sort( key = lambda x: float(x[1]) )
print(beadNPDistanceSet)

#save the output
io.set_structure(proteinStructure)
os.makedirs(outputFolder,exist_ok = True)
print("Saving to: ", outputFolder+"/"+proteinNPID+"_opt.pdb")
io.save(outputFolder+"/"+proteinNPID+"_opt.pdb")
outFile = open( outputFolder+"/"+proteinNPID+"_beaddists.csv","w")
outFile.write("#Type,distance[nm]\n")
for bead in beadNPDistanceSet:
    outFile.write( bead[0]+","+str(bead[1]) +"\n")
outFile.close()

#(rawCoords,resNames) =getAtomCoordsNames(proteinStructurePath)
#mean-center protein, convert to NM
#rawCoords[:,0] = (rawCoords[:,0] - np.mean(rawCoords[:,0]))*0.1
#rawCoords[:,1] = (rawCoords[:,1] - np.mean(rawCoords[:,1]))*0.1
#rawCoords[:,2] = (rawCoords[:,2] - np.mean(rawCoords[:,2]))*0.1

#rotated the PDB into the preferred orientation
#rotCoords =  rotatePDB(rawCoords ,-1*minPhi*np.pi/180, np.pi - minTheta*np.pi/180)
#displace the protein to the optimal location
#rotCoords[:,2] = rotCoords[:,2] - np.amin(rotCoords[:,2]) + optimalBindingLine[4] + npRadius

