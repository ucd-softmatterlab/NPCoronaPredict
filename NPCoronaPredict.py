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




availableMaterials = []
materialDefs = {}
#targetMaterialFile = "MaterialSet.csv"


for targetMaterialFile in ["MaterialSet.csv", "surface-pmfp/MaterialSetPMFP.csv"]:
    if not os.path.exists(targetMaterialFile):
        continue
    materialFile = open(targetMaterialFile ,"r") 
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
#print(availableMaterials)

parser = argparse.ArgumentParser(description = "Parameters for corona calculation")
parser.add_argument('-r','--radius',type=float,help="Radius of the NP [nm]",default=5)
parser.add_argument('-z','--zeta',type=float,help="Zeta potential of the NP [mV]",default=0)
parser.add_argument('-m','--material',type=str,help="Material",default="")
parser.add_argument('-p','--projectname', type=str,help="Name of project", default="test_project")
parser.add_argument('-o','--otherproteins',type=str,help="File containing a set of all proteins to include in format ID,concentration, disables auto-finding structures (all proteins wanted must be included)" , default="")
parser.add_argument('-a','--autorun',type=int,help="Auto-run scripts (0 to disable, 1 = UA only, 2 = UA+BCP, 3 = UA+BCP+CKMC (default)", default=3)
parser.add_argument('-d','--demonstration',type=int,help="If non-zero show the live footage in CoronaKMC", default = 0)
parser.add_argument('-t','--time',type=float,help="Corona simulation run-time in seconds", default = 5e-4)
parser.add_argument('-D','--displace',help="Allow  incoming protein to displace bound protein, nonzero = yes", default = 0, type = int)
parser.add_argument('-A','--accelerate',help="Experimental feature for quasiequilibriation scaling, nonzero = yes", default = 0, type = int)
parser.add_argument('-H','--hamaker',help="Enable Hamaker interaction in UA, default = nonzero = yes, 0 = no", default = 1, type = int)
parser.add_argument('-N','--npfile',help="Sets a target NP file for use in UA", default = "", type=str)
parser.add_argument('-I','--inneroverride',help="For custom .np files with a manually specified inner bound, set this value to the inner bound so BCP can find the correct .uam files, else do not use", default=-1,type=int)

args = parser.parse_args()


#Define some naming conventions used throughout the setup.

NPRadius = int(args.radius) #in nm
NPZeta = int(args.zeta) #in mV
NPMaterial = args.material
CGBeadFile = "beadsets/StandardAABeadSet.csv"
CoronaSimTime = args.time/3600 #1 second in hours
isCylinder = False
isPlane = False
boundaryType = 1
planarRadiusCutoff = 500 #if a spherical NP radius is set larger than this, then the corona tools approximate the NP as a plane. Note that UA will still treat the NP as whatever the MaterialSet shape is, e.g. a 501 nm sphere is a sphere for UA but a plane for BCP/CKMC.

enableHamaker = True
if args.hamaker == 0:
    enableHamaker = False

setupFailed = False
predefNP = False
npName = "autonp"

filenameRadius = NPRadius
if args.inneroverride > 0:
    filenameRadius = args.inneroverride
    print("BCP will be told to look for .uam files with radius ", filenameRadius, " in the name.")

if args.npfile != "":
    print("Using custom NP file for UA input. Note that the radius given as input will be used for corona simulations, please check this is physically meaningful for your system. Passing default material anatase101.")
    predefNP = True
    #npName = args.npfile[:-3]
    npName = (args.npfile.split("/")[-1])[:-3]
    NPMaterial = "anatase101"
else:
    print("NP file will be auto generated")
    if NPMaterial == "":
        print("Please set a material using the -m flag. Available materials are: ")
        print(availableMaterials)
        setupFailed = True
        #raise ValueError("End")
    if NPMaterial not in availableMaterials :
        print("Could not find material ", NPMaterial, " check the spelling. Available materials are: ")
        print(availableMaterials)
        setupFailed = True
        #raise ValueError("End")

if setupFailed == True:
    raise ValueError("End")

if predefNP == False and NPMaterial[-5:] == "-pmfp":
    print("Auto-detected PMFPredictor output material, changing to extended bead set")
    CGBeadFile = "pmfp-beadsetdef/PMFP-BeadSet.csv"

if predefNP == False:
    npShape=int(materialDefs[NPMaterial][3])
else:
    print("Predefined NP, assuming spherical. For cylindrical NPs you will need to run commands individually.")
    npShape = 1
    

if npShape == 2 or  npShape == 4 or npShape == 5:
    isCylinder = True
    print("Enabling cylinder mode")
    boundaryType = 1
if npShape== 3:
    isPlane = True
    print("Enabling planar mode (UA+Corona, vacuum boundary)")
    boundaryType = 0
if NPRadius >= planarRadiusCutoff and npShape == 1:
    isPlane = True
    print("Enabling large-sphere planar mode (Corona only, periodic boundary)")
    boundaryType = 1
    
BaseStorageFolder = "CoronaPredictionProjects"
ProjectName = args.projectname
ProteinStorageFolder = "all_proteins"
ProteinWorkingFolder = BaseStorageFolder+"/"+ProjectName+"/proteins"
UAResultsFolderBase =  BaseStorageFolder+"/"+ProjectName+"/results"
if predefNP == False:
    UAResultsFolder = UAResultsFolderBase+"/np1R_"+str(round(NPRadius))+"_ZP_"+str(round(NPZeta))
else:
    UAResultsFolder = UAResultsFolderBase+"/"+ npName

allFolders = [BaseStorageFolder,ProteinStorageFolder,ProteinWorkingFolder,UAResultsFolderBase,UAResultsFolder,"cg_corona_data"]
for folderName in allFolders:
    if not os.path.isdir(folderName):
        print("Making folder ",folderName)
        os.makedirs(folderName,exist_ok=True)


    
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
    PDBProteinSet = [  ]
    UniProtProteinSet = [ ]
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
        if proteinSource=="AF" :
            #download from AlphaFold
            try:
                os.system('wget  https://alphafold.ebi.ac.uk/files/AF-'+proteinID+'-F1-model_v4.pdb -P '+ProteinStorageFolder+' -O '+ProteinStorageFolder+"/"+proteinID+'.pdb')
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
serumFileLocation = BaseStorageFolder+"/"+ProjectName+"/serum.csv"


serumOutputFile = open(serumFileLocation,"w")
serumOutputFile.write("#ProteinID, NumberConcentration\n")
for protein in successfulProteins:
    serumOutputFile.write(protein[0]+"-"+protein[3]+","+str(protein[1])+"\n")
serumOutputFile.close()

print("Now run UnitedAtom with pdbs set to "+ProteinWorkingFolder)
print("Suggested autorun command: ")

UACommandString = "python3 RunUA.py -r "+str(round(NPRadius))+" -z "+str(NPZeta/1000.0)+" -p "+ProteinWorkingFolder+ " -o "+UAResultsFolderBase+ " -m "+NPMaterial+" --operation-type=pdb-folder -b "+CGBeadFile+" -c "+(BaseStorageFolder+"/"+ProjectName)+" -n "+ProjectName+" -H "+str(args.hamaker)

if predefNP == True:
    #npName = args.npfile[:-3]
    UACommandString += " -N " + args.npfile

print(UACommandString)
kmcFileLocation = BaseStorageFolder+"/"+ProjectName+"/coronakmcinput.csv"
BCPCommandString = "python3 BuildCoronaParams.py -r "+str(round(NPRadius))+" -z "+str(int(NPZeta))+" -f "+UAResultsFolder+" -p "+ serumFileLocation+" -c "+ProteinWorkingFolder+" -b "+CGBeadFile+" -o "+kmcFileLocation
if args.inneroverride > 0:
    BCPCommandString += " -I "+str(filenameRadius)
KMCCommandString = "python3 CoronaKMC.py -r "+str(round(NPRadius))+" -f 0 -p "+kmcFileLocation+" -t "+str(CoronaSimTime)+" --timedelta 0.0000001 -P "+ProjectName+" --demo "+str(args.demonstration)+" -b "+str(boundaryType)+" -D "+str(args.displace)+" -A "+str(args.accelerate)

if isCylinder == True:
    print("Adding cylinder argument")
    BCPCommandString = BCPCommandString+" -s 2"
    KMCCommandString = KMCCommandString+" -s 2"
elif isPlane == True:
    print("Adding plane argument")
    BCPCommandString = BCPCommandString+" -s 3"
    KMCCommandString = KMCCommandString+" -s 3"
    

print(BCPCommandString)
print(KMCCommandString)
print("If python3 isn't installed, use python BuildCoronaParams-p2.py, python CoronaKMC-p2.py instead")


# In[15]:

if args.autorun > 0:
    os.system(UACommandString)
if args.autorun > 1:
    os.system(BCPCommandString)
if args.autorun > 2:
    os.system(KMCCommandString)


# In[ ]:




