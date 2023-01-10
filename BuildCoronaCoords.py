import numpy as np
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB import Structure
import Bio.PDB.Superimposer as Superimposer
import warnings
from Bio import BiopythonWarning
import os

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
    

def rotateToPoint(phi,theta):
    rxx = np.cos(phi)*np.cos(theta)
    rxy = - np.sin(phi)
    rxz = np.cos(phi)*np.sin(theta)
    ryx = np.cos(theta)*np.sin(phi)
    ryy = np.cos(phi)
    ryz = np.sin(phi)*np.sin(theta)
    rzx = -np.sin(theta)
    rzy = 0
    rzz = np.cos(theta)
    return  [  [rxx,rxy,rxz] , [ryx,ryy,ryz], [rzx,rzy,rzz] ]

    
inputFilename = "corona_results_testing/sio2_pristine_birch_h2o_ph47_40_-29.csv_coords_5.0_s0.txt"
proteinDir = "pdbs/birch_pdbs_phmod"

outputFolder = "viewcoronadata-salzburg-5nm-1s-v2"
os.makedirs(outputFolder,exist_ok=True)


coronaData = []
inputFile = open(inputFilename,"r")
coronaLines = inputFile.readlines()
inputFile.close()

for line in coronaLines:
    if line[0]!="#":
        lineTerms = line.strip().split(",")
        proteinName,thetaphi = lineTerms[0].split(":")
        #print(thetaphi)
        theta,phi =  thetaphi.split("-") 
        x = float(lineTerms[1])*10
        y=float(lineTerms[2])*10
        z=float(lineTerms[3])*10
        coronaData.append([proteinName, float(theta), float(phi), x, y, z])
        
print(coronaData)


parser = PDBParser(PERMISSIVE=1)
sup = Superimposer()
io = PDBIO()


proteinNum = 0

completeStructure = Structure.Structure("complete")
csI = 0
coronaOutFile = outputFolder+"/corona_all.pdb"


for proteinLine in coronaData:
    proteinName = proteinLine[0]
    proteinStructurePath = proteinDir+"/"+proteinName+".pdb"
    minPhi = proteinLine[2]
    minTheta = proteinLine[1]
    xc = proteinLine[3] 
    yc = proteinLine[4]
    zc = proteinLine[5]
    rc = np.sqrt(xc**2+yc**2+zc**2)
    comPhi = np.arctan2( yc, xc)
    comTheta = np.arccos( zc/rc)
    print("protein COM at ", xc, yc, zc)
    
    proteinStructure = parser.get_structure("originalstructure",proteinStructurePath)
    #procedure: place the protein in the optimum configuration on top of the NP, then rotate to the correct location
    #find centre-of-CA mass
    
    xcentre = 0
    ycentre=0
    zcentre=0
    numCA=0
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
    print( minPhi, minTheta)
    rotationMatrix = np.array( getRotationMatrix( -minPhi * np.pi/180.0 , np.pi - minTheta*np.pi/180.0    ) )
    print(rotationMatrix)
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
        #print(coords0, coordsFinal)
        atom.set_coord(coordsFinal)
    #next shift along z so that the centre of mass is at the correct location
    for atom in proteinStructure.get_atoms():
       coords0 = atom.get_coord()
       zshift = coords0[2] +rc
       atom.set_coord( [ coords0[0], coords0[1], zshift])
    io.set_structure(proteinStructure)
    io.save(outputFolder+"/"+proteinName+"_"+str(proteinNum)+"_nptop.pdb")
       
    rotationMatrix2 = np.array( rotateToPoint(comPhi,comTheta) )
    print("Rotating COM to ", comPhi, comTheta) 
    for atom in proteinStructure.get_atoms():
        coords0 = atom.get_coord()
        xshift = coords0[0]
        yshift = coords0[1] 
        zshift = coords0[2] 
        xrot =  rotationMatrix2[0,0] * xshift + rotationMatrix2[0,1]*yshift + rotationMatrix2[0,2] * zshift
        yrot =  rotationMatrix2[1,0] * xshift + rotationMatrix2[1,1]*yshift + rotationMatrix2[1,2] * zshift
        zrot =  rotationMatrix2[2,0] * xshift + rotationMatrix2[2,1]*yshift + rotationMatrix2[2,2] * zshift 
        coordsFinal = [xrot,yrot,zrot]
        atom.set_coord(coordsFinal)
        
    #recheck the COM
    xcentre = 0
    ycentre = 0
    zcentre = 0
    for atom in proteinStructure.get_atoms():
        if atom.get_name() != "CA":
            continue
        coords = atom.get_coord()
        xcentre += coords[0]
        ycentre += coords[1]
        zcentre += coords[2]
    xcentre = xcentre/numCA
    ycentre = ycentre/numCA
    zcentre = zcentre/numCA
    print("Positioned COM at ", xcentre, ycentre, zcentre)
    io.set_structure(proteinStructure)
    print("Saving to: ", outputFolder+"/"+proteinName+"_"+str(proteinNum)+"_final.pdb")
    io.save(outputFolder+"/"+proteinName+"_"+str(proteinNum)+"_final.pdb")
    for model in list(proteinStructure):
        modelOut = model.copy()
        modelOut.id = csI
        modelOut.serial_num = csI+1
        csI = csI+1
        completeStructure.add(modelOut)
    proteinNum = proteinNum + 1
io.set_structure(completeStructure)
#io.save(outputFolder+"/completeCorona.pdb")
# "foreach pdb [glob *.pdb ] { mol new $pdb }
