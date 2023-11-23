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
        lineTerms = line.split(",")
        if len(lineTerms)<4:
            print("Problem reading material line: ", line)
            continue
        materialSet[ lineTerms[0]] = [lineTerms[1],lineTerms[2],int(lineTerms[3])]
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
parser.add_argument("-b","--beadset",default="beadsets/StandardAABeadSet.csv", help="Bead parameter file")
parser.add_argument("-c","--configloc",default="", help="Location to save the generated configuration file")
parser.add_argument("-T","--temperature", type=float, default=300.0, help="Nominal temperature")
parser.add_argument("-i","--ionicstrength",type=float,default=0.15,help="Ionic strength in Mol (one-half * sum:conc*chargeSquared)")
parser.add_argument("-n","--name",type=str,default="uaautorun",help="Output file name")
parser.add_argument("-H","--hamaker",type=int,default=1,help="Enable Hamaker interaction (0 to disable, enabled by default)")
parser.add_argument("-N","--nps",type=str,default="",help="NP target [file/folder], leave blank for automatic generation from radius/zeta. If enabled this will override radius, zeta, material.")
args = parser.parse_args()




canRun = 0


if args.material in materialSet:
    pmfFolder,hamakerFile,shape = materialSet[ args.material]
    canRun = 1
elif args.nps != "":
    print("No material found but NP folder set, using default material which may not have sufficient PMFs/Hamaker constants. ")
    canRun = 1
    pmfFolder, hamakerFile, shape = materialSet[ firstMaterial ]
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

beadSetFile = open(args.beadset,"r")
beadNames = []
beadCharges = []
beadRadii = []
useDefaultBeadSet = 0
for line in beadSetFile:
    if line[0]=="#":
        continue
    lineTerms = line.strip().split(",")
    if len(lineTerms) < 3:
        print("Could not read line", line)
        continue
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



