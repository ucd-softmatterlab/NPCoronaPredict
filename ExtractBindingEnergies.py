
import os
import numpy as np

#Scan a target directory and all sub-directories for .uam files
def SearchDirectory(path):
    files = os.listdir(path)
    pdbs  = []
    for handle in files:
        abspath = os.path.join(path, handle)
        if os.path.isdir(abspath):
            pdbs += SearchDirectory(abspath)
        elif abspath[-4:] == ".uam":
            pdbs.append(abspath)
    return pdbs


#given a .uam file, calculate the simple and Boltzmann-weighted average binding energies. 
def CalculateEnergies(filename):
    data     = np.genfromtxt(filename)
    theta    = data[:,1] + 2.5 #left-edge correction
    energy   = data[:,2]
    sinTheta = np.sin(theta * np.pi / 180.0)
    simple    = np.sum(energy * sinTheta ) / np.sum(sinTheta  )
    boltz    = np.sum(energy * sinTheta * np.exp(-1.0 * energy)) / np.sum(sinTheta * np.exp(-1.0 * energy))
    return simple, boltz




#This has to be a list of target folders to look at, the script will scan sub-folders automatically and get all the .uam files it can find
targetFolderNames = ["results"] 


#You can change this filename if needed too
resultsFile = open("bindingEnergyResults.csv","w")



print("FilePath NPName ProteinName SimpleEnergy[kbT] BoltzWeightedEnergy[kbT]")
resultsFile.write("FilePath,NPName,ProteinName,SimpleEnergy[kbT],BoltzWeightedEnergy[kbT]\n")
resList = []
for folderName in targetFolderNames:
#from the folder name we extract the details about what type of NP we're looking at
    folderNameTerms = folderName.split("_")
    uamFiles = (SearchDirectory(folderName)) 
    uamFiles.sort()
    for uam in uamFiles:
        name = uam.split('/')[-1].split('.')[0]
        simpleEnergy,boltzEnergy = CalculateEnergies(uam)
        nameTerms = name.split("_")
        proteinName = nameTerms[0]
        npname = uam.split("/")[-2]
        print(uam.split(".")[0], npname, proteinName, simpleEnergy, boltzEnergy)
        resultsFile.write(uam.split(".")[0]+","+npname+","+ proteinName+","+str(simpleEnergy)+","+str(boltzEnergy)+"\n")
resultsFile.close()



