import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial as scispat
from mpl_toolkits.mplot3d import Axes3D

def getAtomCoords(filename):
    fileIn = open(filename,"r")
    coordList = []
    for line in fileIn:
        lineData = line.split()
        if lineData[0] == "ATOM" and lineData[2] == "CA":
            #coordList.append([ float(lineData[6]) ,  float(lineData[7]) , float(lineData[8])])
            coordList.append([ float(line[30:38]) ,  float(line[38:46]) , float(line[46:54])])
    fileIn.close()
 
    return np.array(coordList)


def getAtomCoordsNames(filename):
    fileIn = open(filename,"r")
    coordList = []
    nameList = []
    for line in fileIn:
        lineData = line.split()
        if lineData[0] == "ATOM":
            #coordList.append([ float(lineData[6]) ,  float(lineData[7]) , float(lineData[8])])
            coordList.append([ float(line[30:38]) ,  float(line[38:46]) , float(line[46:54])])
            nameList.append([line[17:20],line[12:16], line[22:26]   ])
    fileIn.close()
 
    return (np.array(coordList),np.array(nameList))


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

def getAreaHeight(coords):
    projectedConvexHull = scispat.ConvexHull( coords[:,(0,1) ] )
    heightAboveSurface = coords[:,2] - np.mean(coords[:,2]) #Distance from the centre (average z location) to the 
    return (projectedConvexHull.volume, heightAboveSurface)



npRadius = 40

#Tetra.pdb is a test "protein" using the CHL and PHO residues arranged in a tetrahedron with 3 CHL and 1 PHO. CHL carries a +1 charge and PHO a -1 charge, so we expect the orientation with PHO facing away from the NP to occur on NPs with -ve zeta potential and PHO in contact with the NP on NPs of +ve zeta potential
proteinStructurePath = "../pdbs/4A88.pdb"
proteinHeatmapPath = "P02662_40_-50.uam"
'''
proteinHeatmapData     = np.genfromtxt(proteinHeatmapPath)
proteinHeatmapData[:,0] = proteinHeatmapData[:,0] + 2.5 #UA orientation output is left-hand bin edges, this pushes it back to centre-bin edges
proteinHeatmapData[:,1] = proteinHeatmapData[:,1] + 2.5
optimalBindingLine = proteinHeatmapData[ np.argmin( proteinHeatmapData[:,2]   ) ]
#print(optimalBindingLine)
minTheta = optimalBindingLine[1]
minPhi = optimalBindingLine[0]
'''

minTheta = 82.5 
minPhi = 92.5 
(rawCoords,resNames) =getAtomCoordsNames(proteinStructurePath)
#mean-center protein, convert to NM
rawCoords[:,0] = (rawCoords[:,0] - np.mean(rawCoords[:,0]))*0.1
rawCoords[:,1] = (rawCoords[:,1] - np.mean(rawCoords[:,1]))*0.1
rawCoords[:,2] = (rawCoords[:,2] - np.mean(rawCoords[:,2]))*0.1


#print(rawCoords[:,0])
#rotated the PDB into the preferred orientation
rotCoords =  rotatePDB(rawCoords ,-1*minPhi*np.pi/180, np.pi - minTheta*np.pi/180)
#print(rotCoords[:,0])
'''
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(rotCoords[:,0],rotCoords[:,1],rotCoords[:,2])
ax.scatter(rawCoords[:,0],rawCoords[:,1],rawCoords[:,2])
'''


#displace the protein to the optimal location

optDist = 0 #optimalBindingLine[4]
rotCoords[:,2] = rotCoords[:,2] - np.amin(rotCoords[:,2]) #+ optDist + npRadius

resNameArray = np.array(resNames)


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(rotCoords[:,0],rotCoords[:,1],rotCoords[:,2])
#plot an NP
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = npRadius * np.outer(np.cos(u), np.sin(v))
y = npRadius * np.outer(np.sin(u), np.sin(v))
z = npRadius * np.outer(np.ones(np.size(u)), np.cos(v))
ax.plot_surface(x, y, z, color='r')

plotRadius = npRadius+5
ax.set_xlim3d(-plotRadius,plotRadius)
ax.set_ylim3d(-plotRadius,plotRadius)
ax.set_zlim3d(-plotRadius,plotRadius)
plt.show()

atomNum = 1


print("COMPND NP_CORONA")
print("AUTHOR VNPS")
print("CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P1          1")
print("MODEL 0")
for k in range(len(rawCoords)):
    pdbLine= "ATOM  {:5} {} {} A{:4}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}           C".format(atomNum,resNameArray[k,1],resNameArray[k,0],resNameArray[k,2],rotCoords[k,0]*10,rotCoords[k,1]*10,rotCoords[k,2]*10,1,0)
    print(pdbLine)
    #print("ATOM      "+str(atomNum)+"  "+resNameArray[k,0]+"  "+resNameArray[k,1]+" A 106      "+str(round(rotCoords[k,0]*10,2))+"  "+str(round(rotCoords[k,1]*10,2))+"   "+str(round(rotCoords[k,2]*10,2))+"  1.00  0.00           C")
    atomNum+=1

print("MASTER        0    0    0    0    0    0    0    0   12    0   12    0")
print("END")

