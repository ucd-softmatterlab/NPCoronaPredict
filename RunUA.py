import os
import argparse
import numpy as np

#Load in all known materials
materialSet = {}

firstMaterial = ""
for targetMaterialFile in ["MaterialSet.csv", "surface-pmfp/MaterialSetPMFP.csv"]:
    if not os.path.exists(targetMaterialFile):
        continue
    materialFile = open(targetMaterialFile,"r")
    for line in materialFile:
        if line[0] == "#":
            continue
        lineTerms = line.strip().split(",")
        if len(lineTerms)<4:
            print("Problem reading material line: ", line)
            continue
        #print(lineTerms)
        materialSet[ lineTerms[0]] = [lineTerms[1],lineTerms[2],int(lineTerms[3]),float(lineTerms[4])]
        if firstMaterial == "":
            firstMaterial = lineTerms[0]
if len( materialSet.keys() ) == 0:
    print("No materials found. Please try again")
    quit()        
        

#print(materialSet)
parser = argparse.ArgumentParser(description="Parameters for UA Config File Generation")
parser.add_argument('--operation-type', default = "pdb-folder", choices = ['pdb' , 'pdb-folder'], type = str, help = 'Currently only \'pdb\' or pdb-folder are valid')
parser.add_argument("-p","--input-file",  required=True, help="Path to protein PDB file")
parser.add_argument("-r","--radius", type=float,help="NP radius [nm]", default = 5)
parser.add_argument("-z","--zeta", type=float, help="NP zeta potential [V]", default = 0)
parser.add_argument("-o","--outputfolder", default="UAOutput", help="Working folder for results")
parser.add_argument("-m","--material",  choices =  materialSet.keys() ,help="Chosen material", default="none")
parser.add_argument("-P", "--postprocess",default = 1, help="Post-process results")
parser.add_argument("-b","--beadset",default="", help="Bead parameter file, leave blank to auto-chose")
parser.add_argument("-c","--configloc",default="", help="Location to save the generated configuration file")
parser.add_argument("-T","--temperature", type=float, default=300.0, help="Nominal temperature")
parser.add_argument("-i","--ionicstrength",type=float,default=0.15,help="Ionic strength in Mol (one-half * sum:conc*chargeSquared)")
parser.add_argument("-n","--name",type=str,default="uaautorun",help="Output file name")
parser.add_argument("-H","--hamaker",type=int,default=1,help="Enable Hamaker interaction (0 to disable, enabled by default)")
parser.add_argument("-N","--nps",type=str,default="",help="NP target [file/folder], leave blank for automatic generation from radius/zeta. If enabled this will override radius, zeta, material.")
parser.add_argument("-j","--jitter",type=float, help= "S. dev. of random noise to apply to each CG bead position per-axis [nm]", default=0.0)
parser.add_argument("-B", "--boltzmode",type=int,default=0, help="If >0 enables Boltzmann local averaging in UA")
parser.add_argument("-S", "--shapeoverride", type=int, default=-1 , help="If  > 0 overrides the default shape for a given material to the specified shape number")
parser.add_argument("-L","--ligand-file", type=str, default = "", help = "Path to a UA ligand override file, leave blank to skip")
parser.add_argument("-f","--flexres", type=float, default = -0.01, help="Resolution for flexible beads if > 0.005, else disabled")
parser.add_argument("--relaxsteps", type=int, default = 0, help="Number of relaxation steps")

args = parser.parse_args()



flexOn = False
flexRes = 0
if args.flexres > 0.005:
    flexOn = True
    flexRes = args.flexres
    print("enabling flexibility")


relaxOn = False
relaxSteps = 0
if args.relaxsteps > 0:
    relaxSteps = int(args.relaxsteps)
    print("enabling relaxation (fewer samples will be used in UA)")
    relaxOn = True
canRun = 0
checkMaterial =""

isNPFile = False
if args.material in materialSet:
    pmfFolder,hamakerFile,shape,pmfLJCutoff = materialSet[ args.material]
    canRun = 1
    checkMaterial = args.material
elif args.nps != "":
    print("No material found but NP folder set, using default material which may not have sufficient PMFs/Hamaker constants. ")
    canRun = 1
    checkMaterial = firstMaterial
    pmfFolder, hamakerFile, shape,pmfLJCutoff = materialSet[ firstMaterial ]


if args.nps != "":
    isNPFile = True

if canRun == 0:
    print("An error has occured, cancelling run. Please check material name again.")
    print("Input material: ", args.material)
    print("Known materials: ", materialSet.keys() )
    quit()

useNPFolder = False
if args.nps != "":
    useNPFolder = True
    npTargetFolder = args.nps

enableHamaker = True
if args.hamaker == 0:
    enableHamaker = False

enableBoltz = 0
if args.boltzmode > 0:
    enableBoltz = 1

jitterMag = 0
if args.jitter > 0.001:
    jitterMag = args.jitter


tubeCorrectionShapes = [4,5]
planeCorrectionShapes = [1,2,3]

overrideShape = False
if args.shapeoverride > 0:
    overrideShape = True
    newShape = args.shapeoverride
    if shape in tubeCorrectionShapes and newShape in planeCorrectionShapes:
        print("Material is a tube PMF format but you have selected a planar NP shape. Please be aware this may cause errors")
    if newShape in tubeCorrectionShapes and shape in planeCorrectionShapes:
        print("Material is a plane PMF format but you have selected a tube NP  shape. Please be aware this may cause errors")
    shape = newShape


'''
if predefNP == False and NPMaterial[-5:] == "-pmfp":
    print("Auto-detected PMFPredictor output material, changing to extended bead set")
    CGBeadFile = "pmfp-beadsetdef/PMFP-BeadSet.csv"
'''

chosenBeadSet =  "beadsets/StandardAABeadSet.csv"

if args.material[-5:] == "-pmfp":
    chosenBeadSet = "pmfp-beadsetdef/PMFP-BeadSet.csv"

if args.beadset != "":
    chosenBeadSet = args.beadset



beadSetFile = open(chosenBeadSet,"r")
beadNames = []
beadCharges = []
beadRadii = []
useDefaultBeadSet = 0

allPMFFolders = []
allHamakerFiles = []

def SearchNPDirectory(path):
    files = os.listdir(path)
    pdbs  = []
    for handle in files:
        abspath = os.path.join(path, handle)
        if os.path.isdir(abspath):
            pdbs += SearchDirectory(abspath)
        elif abspath[-3:] == ".np":
            pdbs.append(abspath)
    return pdbs


def getNPPMFs( npFilePath):
    npMaterials = []
    npFileIn = open(npFilePath,"r")
    for line in npFileIn:
        #print(line)
        if line[0]=="#":
            continue
        lineTerms = line.split(",")
        pmfDir = lineTerms[9]
        #print(pmfDir)
        if pmfDir not in npMaterials:
            npMaterials.append(pmfDir)
            #print("registering", pmfDir)
    return list(set(npMaterials))

def getNPHamakers( npFilePath):
    npMaterials = []
    npFileIn = open(npFilePath,"r")
    for line in npFileIn:
        #print(line)
        if line[0]=="#":
            continue
        lineTerms = line.split(",")
        pmfDir = lineTerms[8]
        #print(pmfDir)
        if pmfDir not in npMaterials:
            npMaterials.append(pmfDir)
            #print("registering", pmfDir)
    return list(set(npMaterials))


if isNPFile == True:
    allHamakerFiles =[]
    #scan over NP file and get all the PMF folders to include, including resetting some of the default parameters to try to keep things functioning
    if args.nps[-3:] == ".np":
        print("checking NP file")
        allMaterialPMFs = getNPPMFs( args.nps )
        allPMFFolders = allPMFFolders + allMaterialPMFs
        allHamakerFiles = allHamakerFiles + getNPHamakers(args.nps)
    else:
        allNPFiles = SearchNPDirectory( args.nps ) 
        for npFile in allNPFiles:
            allMaterialPMFs  = getNPPMFs( npFile )
            allPMFFolders = list(set(   allPMFFolders + allMaterialPMFs) )
            allHamakerFiles = list(set( allHamakerFiles + getNPHamakers(npFile) )) 
    #pmfFolder, hamakerFile, shape,pmfLJCutoff = materialSet[ firstMaterial ]
    pmfFolder = allPMFFolders[0]
    hamakerFile = allHamakerFiles[0] 
    print("setting defaults to ", pmfFolder, hamakerFile )
else:
    allPMFFolders.append(pmfFolder) 
    
print(allPMFFolders)
#quit()

for line in beadSetFile:
    if line[0]=="#":
        continue
    lineTerms = line.strip().split(",")
    if len(lineTerms) < 3:
        print("Could not read line", line)
        continue
    #test if bead exists for this material
    beadName = lineTerms[0]
    materialNotFound = False
    for pmfFolderCheck in allPMFFolders:
        if not os.path.exists( pmfFolderCheck+"/"+ beadName +".dat" ):
            print("Bead type ", lineTerms[0], " could not be located for the selected material", pmfFolderCheck)
            print("This will be omitted from the UA input file and may cause issues with ligands or non-AA beads")
            materialNotFound = True
    if materialNotFound == True:
        continue
    else:
        #print("bead defined for all materials")
        beadNames.append(lineTerms[0])
        beadCharges.append(lineTerms[1])
        beadRadii.append(lineTerms[2])


if len(beadNames) == 0:
    print("No AA beads found in file, attempting with default parameters")
    useDefaultBeadSet = 1
if len(beadNames)!=len(beadCharges) or len(beadNames)!=len(beadRadii):
    print("Mismatch between number of bead parameters, check input. Attempting with default parameters")
    useDefaultBeadSet = 1
beadNameString = ", ".join(beadNames)
beadChargeString = ", ".join( [ str(round(float(a),3)) for a in beadCharges])
beadRadiiString = ", ".join( [ str( round(float(a),3)) for a in beadRadii])

configOutputName = args.name+"_configautogen.config"

inputTemp = args.temperature
inputTempC = inputTemp - 273.15
if inputTempC > 100.0:
    print("Warning: temperature greater than 100 C, medium is likely to evaporate")
if inputTempC < 0.0:
    print("Warning: temperature less than 0 C, medium is likely to freeze")

chargeUnit = 1.6e-19
epsRel = 87.740 - 0.4008*inputTempC + 9.398*(10**(-4))*(inputTempC**2)  - 1.410*(10**(-6))*(inputTempC**3)#malmberg and maryott, J. Res. National Bureau of Standards, 56, 1956, "Dielectric constant of water from 0 deg to 100 deg C
print("Estimated eps rel at temperature T=", inputTempC, ": ", epsRel)
epsZero = 8.854e-12
nAvo = 6.02214076e23
kb = 1.380649e-23;
ionicStrengthNumber = args.ionicstrength * nAvo * 1000
#calculate the Debye (useful) and Bjerrum (not useful but mandatory) lengths
debyeLength = 1e9 * np.sqrt( epsRel * epsZero * kb * inputTemp/ (2 * (chargeUnit**2) * ionicStrengthNumber) )
print("Calculated debye length:" ,  str(round(debyeLength,3)))
bjerrumLength = 1e9*chargeUnit**2/(4 * np.pi * epsZero * epsRel * kb * inputTemp )
print("Calculated Bjerrum length:" ,  str(round(bjerrumLength ,3)))

def writeConfigFile(configOutputLoc):
    outputConfigFile = open(configOutputLoc,"w")
    outputConfigFile.write("#Autogenerated UA Config file\n")
    outputConfigFile.write("output-directory = "+args.outputfolder+"\n")
    outputConfigFile.write("pdb-target = "+args.input_file+"\n")
    if useNPFolder == True:
        outputConfigFile.write("np-target = "+npTargetFolder+"\n")
    outputConfigFile.write("nanoparticle-radius = [" + str(args.radius)+"]\n")
    outputConfigFile.write("np-type = " + str(shape)+"\n")
    outputConfigFile.write("pmf-directory = " + pmfFolder + "\n")
    outputConfigFile.write("hamaker-file = " + hamakerFile + "\n")
    outputConfigFile.write("enable-surface \n")
    if enableHamaker == True:
        outputConfigFile.write("enable-core \n")
    outputConfigFile.write("enable-electrostatic \nsimulation-steps = 2000 \npotential-cutoff=5.0 \npotential-size = 1000 \nangle-delta = 5.0 \nbjerum-length="+str(round(bjerrumLength ,3))+" \ndebye-length="+str(round(debyeLength,3))+" \n")
    outputConfigFile.write("temperature = "+str(round(inputTemp,2))+"\n")
    outputConfigFile.write("zeta-potential = [" + str(args.zeta) + "] \n")
    outputConfigFile.write("pdb-jitter-magnitude = "+str(jitterMag)+" \n")
    outputConfigFile.write("pmf-cutoff="+str( pmfLJCutoff)+"\n")

    if flexOn == True:
        outputConfigFile.write("bead-flexibility-method=3\n")
        outputConfigFile.write("flex-sdev=3\n")
        outputConfigFile.write("flex-resolution="+str(flexRes)+"\n")
    if relaxOn == True:
        outputConfigFile.write("num-random-samples = 16\n")
        outputConfigFile.write("enable-pdb-relax\n")
        outputConfigFile.write("enable-local-boltz\n")
        outputConfigFile.write("relax-steps="+str(relaxSteps)+"\n")
        outputConfigFile.write("relax-gradient=0.0\n")
    if args.ligand_file != "":
        outputConfigFile.write("ligand-file = "+args.ligand_file+"\n")

    if enableBoltz > 0:
        outputConfigFile.write("enable-local-boltz \n")
    if useDefaultBeadSet == 1:
        outputConfigFile.write("amino-acids         = [ ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY, HIS, ILE, LEU, LYS, MET, PHE, PRO, SER, THR, TRP, TYR, VAL] \n")
        outputConfigFile.write("amino-acid-charges  = [ 0.0, 1.0, 0.0, -1.0, 0.0, 0.0, -1.0, 0.0, 0.5, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ] \n")
        outputConfigFile.write("amino-acid-radii    = [ 0.323, 0.429, 0.362, 0.356, 0.352, 0.386, 0.376, 0.285, 0.302, 0.401, 0.401, 0.405, 0.402, 0.421, 0.362, 0.328, 0.357, 0.449, 0.425, 0.377] \n")
    else:
        outputConfigFile.write("amino-acids         = [ "+beadNameString+"] \n")
        outputConfigFile.write("amino-acid-charges  = [ "+beadChargeString+"] \n")
        outputConfigFile.write("amino-acid-radii    = [ "+beadRadiiString+"] \n")
    outputConfigFile.close()
writeConfigFile(configOutputName)
if args.configloc != "":
    writeConfigFile(args.configloc+"/"+configOutputName)


print("Generated config file, running UA")
os.system("./UnitedAtom --config-file="+configOutputName)
print("UA run complete")

NPRadius = int(args.radius) #in nm
NPZetaMV = int(1000*args.zeta) #in mV

finalOutputFolder = args.outputfolder+"/np1R_"+str(round(NPRadius))+"_ZP_"+str(round(NPZetaMV))
if args.postprocess == 1 and args.operation_type=="pdb":
    proteinTargetName = (args.input_file.split("/")[-1])[:-4]
    if shape==1:
        uaResultFileNameOriginal = proteinTargetName+"_"+str(int(args.radius))+"_"+str(int(1000*args.zeta))+".uam"
    else:
        uaResultFileNameOriginal = proteinTargetName+"_"+str(int(args.radius))+"_"+str(int(1000*args.zeta))+"_0.uam"
    uaResultFileName = proteinTargetName+"_"+args.material+"_"+str(args.radius)+"_"+str(args.zeta)+".uam"
    os.system("mv "+finalOutputFolder+"/"+uaResultFileNameOriginal+" "+ finalOutputFolder+"/"+uaResultFileName)
    print("Generating heatmap")
    os.system("python3 tools/plotmap "+ finalOutputFolder+"/"+uaResultFileName)
    print("Generating rotated PDB")
    os.system("python3 tools/ApplyOptimumRotation.py -p " + args.input_file + " -u " + finalOutputFolder+"/"+uaResultFileName + " -o "+finalOutputFolder)



