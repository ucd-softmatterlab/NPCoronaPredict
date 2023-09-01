'''Decompose a target molecule defined by a SMILES code into fragments using the BRICS algorithm, then matches these fragments to pre-existing fragments where possible'''

from rdkit import Chem
from rdkit.Chem import BRICS
from rdkit.Chem import AllChem
from rdkit import RDLogger
import numpy as np
import datetime
import os


#Define the molecule you want to process - name can be anything 
targetChemName = "TestChemical"
targetMolSMILES = "CC(NCC(O)COC1=CC=CC2=CC=CC=C21)C"


#Preparation
RDLogger.DisableLog('rdApp.*') 
dummyAtom = Chem.MolFromSmiles("*")
baseDir = "CGMolecules"
os.makedirs(baseDir,exist_ok=True)
dateString =  datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
shortDate = datetime.datetime.now().strftime("%d-%b-%y")


#load in a library of pre-defined fragments with the form beadID, smiles code,  UA 3-letter code
canonSmilesLookup = {}
libraryFile=open("FragmentLookup.csv","r")
for line in libraryFile:
   if line[0] == "#":
       continue
   lineTerms = line.split(",")
   entrySMILES = lineTerms[2].replace("<HASH>","#").replace("<COMMA>",",")
   canonSMILES = Chem.CanonSmiles(entrySMILES)
   if canonSMILES not in canonSmilesLookup:
       canonSmilesLookup[canonSMILES] = [ lineTerms[0], lineTerms[1] ]
libraryFile.close()
#print(canonSmilesLookup) 


def removeDummySMILES(smilesIn):
    cin = Chem.MolFromSmiles(smilesIn)
    return Chem.MolToSmiles( Chem.RemoveHs(  Chem.AllChem.ReplaceSubstructs(cin, dummyAtom, Chem.MolFromSmiles("[H]"), True)[0] ))


def getChargeFromSmiles(smilesIn):
    return Chem.GetFormalCharge( Chem.MolFromSmiles(smilesIn) )

def getFragmentSMILES(frag):
    return Chem.CanonSmiles(Chem.MolToSmiles( Chem.RemoveHs(  Chem.AllChem.ReplaceSubstructs(frag, dummyAtom, Chem.MolFromSmiles("[H]"), True)[0] )) )


def beadsetToPDB(chemName, chemSmiles, beadset):
    '''Writes out a set of beads [beadID, CodeID, x, y, z, smiles, longname] to a pdb-like file for use in UnitedAtom'''
    beadOut = open(baseDir+"/"+chemName+".pdb","w")
    beadOut.write("HEADER "+ '{:40.40}'.format(chemName)+shortDate+"\n")
    #beadOut.write("REMARK Bead ID:" + uaTargets[chem]+"\n")
    beadOut.write("TITLE   "+chemName+"\n")
    beadOut.write("REMARK SMILES: "+chemSmiles+"\n")
    beadOut.write("REMARK Generated: "+dateString+"\n")
    for atom in beadset:
        beadOut.write("REMARK "+str(atom[0])+":"+atom[6]+":"+atom[1]+":"+atom[5]+"\n")
    for atom in beadset:
        beadOut.write("ATOM  "+str(atom[0]).rjust(5)+" "+" CA "+" "+atom[1]+" A   1    {:8.3}{:8.3}{:8.3}  1.00  0.00\n".format(atom[2],atom[3],atom[4]))
    beadOut.write("END\n")
    beadOut.close()




targetMol  = Chem.MolFromSmiles(targetMolSMILES)
for atom in targetMol.GetAtoms():
    atom.SetIntProp("index0", atom.GetIdx())

targetMolH = Chem.AddHs( targetMol)
AllChem.EmbedMolecule(targetMolH)
positionList = np.array(targetMolH.GetConformer().GetPositions() )

molBroken = Chem.BRICS.BreakBRICSBonds(targetMol)
molFragments = Chem.GetMolFrags(molBroken, asMols=True)
numMissingFragments = 0
beadID = 0
pmfpredOutList = []
finalBeadSet = []
for frag in molFragments:
    fragSMILES = getFragmentSMILES(frag)
    print("Using fragment",   fragSMILES )
    beadCode = "ZZZ"
    beadName = "NewFragment"+str(numMissingFragments)
    beadLoc = np.array([0,0,0])
    totalBeadMass = 0.00001
    for atom in frag.GetAtoms():
        #print(atom.GetIdx())
        try: 
            atomBaseIndex = atom.GetIntProp("index0")
            #print("Original index: ", atomBaseIndex, "location: ",  positionList[   atom.GetIntProp("index0") ] )
            beadLoc = beadLoc +  atom.GetMass() * positionList[   atom.GetIntProp("index0") ]
            totalBeadMass += atom.GetMass()
        except:
            #print("Dummy atom, skipping")
            skipped=1
    beadLoc = beadLoc/totalBeadMass
    #print(beadLoc)
    if fragSMILES in canonSmilesLookup:
        #print( "Found", fragSMILES, "mapped to",  canonSmilesLookup[fragSMILES] )
        beadCode = canonSmilesLookup[fragSMILES][1]
        beadName = canonSmilesLookup[fragSMILES][0]
    else:
        #print("Need to generate: ", fragSMILES)
        charge = getChargeFromSmiles(fragSMILES)
        description="Auto-generated fragment"
        pmfpredOutList.append( beadName+","+fragSMILES+","+str(charge)+","+description )
        numMissingFragments += 1
    finalBeadSet.append( [beadID, beadCode, beadLoc[0], beadLoc[1], beadLoc[2]  , fragSMILES, beadName ] )
    beadID += 1
beadsetToPDB(targetChemName,targetMolSMILES, finalBeadSet)
if len(pmfpredOutList) > 0:
    print("-----")
    print("Fragments requiring PMF generation (Name,SMILES,FormalCharge,Description). Please assign names and descriptions, generate PMFs, update this script and rerun")        
    for newFrag in pmfpredOutList:
        print(newFrag)
else:
    print("All fragments successfully mapped")
