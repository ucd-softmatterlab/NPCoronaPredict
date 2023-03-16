#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#This notebook is designed to simplify the procedure for setting up and running a set of CoronaKMC simulations
#Because both the UnitedAtom and CoronaKMC steps can take a significant amount of time, these are not run directly
#Instead it prints out the commands you should run for CoronaKMC 
#For help please ask Ian Rouse, ian.rouse@ucd.ie


# In[14]:


import os
import argparse


parser = argparse.ArgumentParser(description = "Parameters for corona calculation")
parser.add_argument('-r','--radius',type=float,help="Radius of the NP [nm]",default=5)
parser.add_argument('-z','--zeta',type=float,help="Zeta potential of the NP [mV]",default=0)
parser.add_argument('-m','--material',type=str,help="Material",default="anatase101")
parser.add_argument('-p','--projectname', type=str,help="Name of project", default="test_project")
parser.add_argument('-o','--otherproteins',type=str,help="File containing a set of other proteins to include in format ID,concentration" , default="")
parser.add_argument('-a','--autorun',type=int,help="Auto-run scripts (0 to disable)", default=1)
args = parser.parse_args()


#Define some naming conventions used throughout the setup.

NPRadius = int(args.radius) #in nm
NPZeta = int(args.zeta) #in mV
NPMaterial = args.material
CGBeadFile = "beadsets/StandardAABeadSet.csv"
CoronaSimTime = 1e-5
isCylinder = False

availableMaterials = []
materialDefs = {}
materialFile = open("MaterialSet.csv","r") 
for line in materialFile:
    if line[0] == "#":
        continue
    lineTerms = line.split(",")
    if len(lineTerms)<4:
        print("Problem reading material line: ", line)
        continue
    availableMaterials.append(lineTerms[0])
    materialDefs[ lineTerms[0] ] = lineTerms
    #materialSet[ lineTerms[0]] = [lineTerms[1],lineTerms[2],int(lineTerms[3])]
#print(materialSet)
print(availableMaterials)

#availableMaterials = ["silicaquartz","silicaamorph","anatase100","anatase101","rutile110","rutile100","fe2o3","CdSe","gold","carbonblack"]
if NPMaterial not in availableMaterials:
    print("Could not find material ", NPMaterial, " check the spelling. ")
    raise ValueError("End")

npShape=int(materialDefs[NPMaterial][3])
if npShape == 2 or  npShape == 4 or npShape == 5:
    isCylinder = True
    print("Enabling cylinder mode")
ProjectName = args.projectname
ProteinStorageFolder = "all_proteins"
ProteinWorkingFolder = "proteins_"+ProjectName
UAResultsFolderBase = "results_"+ProjectName
UAResultsFolder = UAResultsFolderBase+"/np1R_"+str(round(NPRadius))+"_ZP_"+str(round(NPZeta))

allFolders = [ProteinStorageFolder,ProteinWorkingFolder,UAResultsFolderBase,UAResultsFolder,"cg_corona_data"]
for folderName in allFolders:
    if not os.path.isdir(folderName):
        print("Making folder ",folderName)
        os.mkdir(folderName)


    
#Defines the set of proteins that calculations should be run for. 
#In all cases, if you have no proteins in this category leave the list empty
#E.g. if you have no proteins with structures from the PDB then leave this as
#  PDBProteinSet = []



#These have the format [UniprotID, number concentration,label] and should be separated by commas
#Structures for these proteins are found from AlphaFold
UniProtProteinSet = [
    ["Q9BYF1", 1e-3,"ACE2Human"],
    ["A4GE70", 1e-3,"CoffeeEnzyme"],
    ["P02769", 1e-3,"BSA"]
]


#Structures for this set are instead found from the RSC PDB repository
PDBProteinSet = [
    ["1AX8", 1e-3,"Leptin"]
]


#Structures for these are taken directly from the storage folder
OtherProteinSet =[
]


proteinItNum = 0
loadSerumFile = "" 
if args.otherproteins != "":
    loadSerumFile = args.otherproteins # "daphnia_serum/daphnia_serum_mass1.0_fixedmassfrac.csv"
if loadSerumFile!= "":
    serumFile = open(loadSerumFile,"r")
    for line in serumFile:
        if line[0]=="#":
            continue
        lineTerms = line.strip().split(",")
        proteinItNum+=1
        OtherProteinSet.append([ lineTerms[0],lineTerms[1],"P"+str(proteinItNum) ])



AllProteins = []
for protein in UniProtProteinSet:
    AllProteins.append( [protein[0],protein[1], "AF", protein[2]])
for protein in PDBProteinSet:
    AllProteins.append( [protein[0],protein[1], "PDB",protein[2]])
for protein in OtherProteinSet:
    AllProteins.append( [protein[0],protein[1], "Other",protein[2]])
    
    
#For all the proteins specified in the list
#    Check to see if the working folder already has the structure
#    Check to see if the protein already has a structure available in the storage folder
#    If so, copy this to the working folder.
#    If not, download from AlphaFold/PDB assuming
#        If this succeeds, copy to the working folder
#        If this fails, print an error and move on to the next protein

successfulProteins = []
for proteinLine in AllProteins:
    foundProtein = 0
    proteinID = proteinLine[0]
    proteinSource = proteinLine[2]
    if os.path.exists( ProteinWorkingFolder+"/"+proteinID+".pdb"  ):
        successfulProteins.append( proteinLine)
        print("Found "+proteinID+" in working folder")
        continue
    print("Looking for ", ProteinStorageFolder+"/"+proteinID+".pdb")
    if os.path.exists( ProteinStorageFolder+"/"+proteinID+".pdb"):
        os.system('cp '+ProteinStorageFolder+"/"+proteinID+".pdb "+ProteinWorkingFolder+"/"+proteinID+"-"+proteinLine[3]+".pdb")
        print("Found "+proteinID+" in storage folder, copied to working")
        foundProtein = 1
    else:
        if proteinSource=="AF":
            #download from AlphaFold
            try:
                os.system('wget  https://alphafold.ebi.ac.uk/files/AF-'+proteinID+'-F1-model_v2.pdb -P '+ProteinStorageFolder+' -O '+ProteinStorageFolder+"/"+proteinID+'.pdb')
                os.system('cp '+ProteinStorageFolder+"/"+proteinID+".pdb "+ProteinWorkingFolder+"/"+proteinID+"-"+proteinLine[3] +".pdb")
                foundProtein = 1
            except:
                print("AlphaFold download failed, please try manually")
                foundProtein = 0
        elif proteinSource=="PDB":
            #download from PDB
            try:
                os.system('wget  https://files.rcsb.org/download/'+proteinID+'.pdb -P '+ProteinStorageFolder)
                os.system('cp '+ProteinStorageFolder+"/"+proteinID+".pdb "+ProteinWorkingFolder+"/"+proteinID+"-"+proteinLine[3]+   ".pdb")
                foundProtein = 1
            except:
                print("PDB download failed, please try manually")
                foundProtein = 0
        else:
            print("Could not find a structure for protein "+proteinID)
    if foundProtein == 1:
        successfulProteins.append(proteinLine)

if len(successfulProteins) == 0:
     raise ValueError("No proteins were found. Please check again")
serumOutputFile = open(ProjectName+"_serum.csv","w")
serumOutputFile.write("#ProteinID, NumberConcentration\n")
for protein in successfulProteins:
    serumOutputFile.write(protein[0]+"-"+protein[3]+","+str(protein[1])+"\n")
serumOutputFile.close()

print("Now run UnitedAtom with pdbs set to "+ProteinWorkingFolder)
print("Suggested autorun command: ")

UACommandString = "python3 RunUA.py -r "+str(round(NPRadius))+" -z "+str(NPZeta/1000.0)+" -p "+ProteinWorkingFolder+ " -o "+UAResultsFolderBase+ " -m "+NPMaterial+" --operation-type=pdb-folder -b "+CGBeadFile
print(UACommandString)

BCPCommandString = "python3 BuildCoronaParams-P3.py -r "+str(round(NPRadius))+" -z "+str(int(NPZeta))+" -f "+UAResultsFolder+" -p "+ProjectName+"_serum.csv -c "+ProteinWorkingFolder
KMCCommandString = "python3 CoronaKMC-P3.py -r "+str(round(NPRadius))+" -f 0 -p cg_corona_data/"+UAResultsFolder+".csv -t "+str(CoronaSimTime)+" --demo 1 --timedelta 0.00001 -P "+ProjectName

if isCylinder == True:
    print("Adding cylinder argument")
    BCPCommandString = BCPCommandString+" -s 2"
    KMCCommandString = KMCCommandString+" -s 2"


print(BCPCommandString)
print(KMCCommandString)
print("If python3 isn't installed, use python BuildCoronaParams.py, python CoronaKMC.py instead")


# In[15]:

if args.autorun != 0:
    os.system(UACommandString)
    os.system(BCPCommandString)
    os.system(KMCCommandString)


# In[ ]:




