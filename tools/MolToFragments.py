'''Decompose a target molecule defined by a SMILES code into fragments using the BRICS algorithm, then matches these fragments to pre-existing fragments where possible'''

from rdkit import Chem
from rdkit.Chem import BRICS, Draw
from rdkit.Chem import AllChem
from rdkit import RDLogger
import numpy as np
import datetime
import os
import argparse
import itertools


parser = argparse.ArgumentParser(description="Parameters for PMFPredictor")
parser.add_argument("-m","--molecule", type=str,default="NewMolecule", help="Name for new molecle")
parser.add_argument("-s","--smiles", type=str,default="CCCC",help="SMILES string for new molecule")
parser.add_argument("-M","--method",type=str,default="EqualParts",help="Method used to decompose molecule. BRICS: use BRICS only. PMFP: Pre-match to existing fragments where possible. EqualParts: Recursive bisection by breakable bonds, matching to existing PMFs where possible ")
args = parser.parse_args()


targetChemName = args.molecule
targetMolSMILES = Chem.CanonSmiles(args.smiles)
method = args.method


#Preparation
RDLogger.DisableLog('rdApp.*') 
dummyAtom = Chem.MolFromSmiles("*")
baseDir = "CGMolecules"
os.makedirs(baseDir,exist_ok=True)
os.makedirs(baseDir+"_images",exist_ok=True)
dateString =  datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
shortDate = datetime.datetime.now().strftime("%d-%b-%y")

missingBeadString = "!!!"

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


def beadsetToPDB(chemName, chemSmiles, beadset,outputFolder=baseDir):
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
    print("Written output file to ", outputFolder+"/"+chemName+".pdb")


def fragmentsToBeads(targetMol, molFragments):
    targetMolH = Chem.AddHs( targetMol)
    AllChem.EmbedMolecule(targetMolH)
    positionList = np.array(targetMolH.GetConformer().GetPositions() )
    numMissingFragments = 0
    beadID = 0
    pmfpredOutList = []
    finalBeadSet = []
    reverseBeadLookup = {}
    for frag in molFragments:
        fragSMILES = getFragmentSMILES(frag)
        print("Using fragment",   fragSMILES )
        beadCode = missingBeadString
        beadName = targetChemName+"New"+str(numMissingFragments)
        beadLoc = np.array([0,0,0])
        totalBeadMass = 0.00001
        for atom in frag.GetAtoms():
            #print(atom.GetIdx())
            try: 
                atomBaseIndex = atom.GetIntProp("index0")
                reverseBeadLookup[atomBaseIndex] = beadID
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
            canonSmilesLookup[fragSMILES] = [ beadName, missingBeadString] 
        finalBeadSet.append( [beadID, beadCode, beadLoc[0], beadLoc[1], beadLoc[2]  , fragSMILES, beadName ] )
        beadID += 1
    beadsetToPDB(targetChemName,targetMolSMILES, finalBeadSet)
    
    for atom in targetMol.GetAtoms():
        try:
            cgBeadID = reverseBeadLookup[atom.GetIntProp("index0")]
            atom.SetProp("atomNote", str(cgBeadID) )
        except:
            skipped = 1
    Chem.Draw.MolToFile( targetMol, baseDir+"_images/"+targetChemName+".png", size=(600,600)  )
    
    if len(pmfpredOutList) > 0:
        print("-----")
        print("Fragments requiring PMF generation (Name,SMILES,FormalCharge,Description). Please assign names and descriptions, generate PMFs, update this script and rerun")        
        extraBeadFile = open("SuggestedNewBeads.csv","a")
        for newFrag in pmfpredOutList:
            print(newFrag)
            extraBeadFile.write( newFrag+"\n" )
        extraBeadFile.close()
    else:
        print("All fragments successfully mapped")
        
        

def fragmentCostFunction(frag):
    oversizeCost = 5
    targetSize = 5 #fragments larger than this are penalised by oversizeCost per extra  breakable bond
    totalCost = 0
    fragSize = frag.GetNumAtoms()
    bisFrag = frag.GetSubstructMatches(Chem.MolFromSmarts('[!R][!R]')) + frag.GetSubstructMatches(Chem.MolFromSmarts('[R][!R]'))
    breakableBondsFrag =[frag.GetBondBetweenAtoms(x,y).GetIdx() for x,y in bisFrag]
    totalCost += oversizeCost*(fragSize - targetSize)**2
    totalCost += oversizeCost*(  max(0, len(breakableBondsFrag) - 3 )  )**2
    return totalCost
        
        
def getBreakableBonds(frag): #get all bonds which aren't between two ring atoms
    bis = frag.GetSubstructMatches(Chem.MolFromSmarts('[!R][!R]')) + frag.GetSubstructMatches(Chem.MolFromSmarts('[R][!R]'))
    #breakableBonds =   [frag.GetBondBetweenAtoms(x,y).GetIdx() for x,y in bis]        
    nonringBonds = [frag.GetBondBetweenAtoms(x,y)  for x,y in bis] 
    breakableBonds = []
    for bond in nonringBonds:
        if bond.GetBondTypeAsDouble() == 1.0:
            breakableBonds.append( bond.GetIdx() )
    return breakableBonds
    
def bisectFrag(frag): #break a fragment along one bond such that the two remaining fragments have roughly equal numbers of bonds remaining
    breakableBonds = getBreakableBonds(frag)
    bestSplitBond = breakableBonds[0]
    bestScore = 10 + len(breakableBonds)*np.log(0.001 + len(breakableBonds) )   
    bestSizeDisparity = frag.GetNumAtoms()
    #print(breakableBonds)
    for trialBond in breakableBonds:
        trialFrags = Chem.GetMolFrags(  Chem.FragmentOnBonds(frag ,[trialBond], addDummies=False) , asMols=True) #split along the bond and count the number of bonds remaining in each
        l1=len(getBreakableBonds(trialFrags[0]))
        l2=len(getBreakableBonds(trialFrags[1]))
        score = l1 * np.log(0.0001+l1) + l2*np.log(0.0001+l2)
        #print(trialBond,score)
        if score < bestScore:
            bestSplitBond  = trialBond
            bestScore = score
            bestSizeDisparity = np.abs( trialFrags[0].GetNumAtoms() - trialFrags[1].GetNumAtoms() )
        if score  < bestScore + 0.1: #tie-breaker - select the split which produces the most equally-sized fragments
            newSizeDisparity = np.abs( trialFrags[0].GetNumAtoms() - trialFrags[1].GetNumAtoms() )
            if newSizeDisparity < bestSizeDisparity :
                bestSplitBond  = trialBond
                bestScore = score
                bestSizeDisparity = newSizeDisparity
    #print(bestScore, bestSplitBond)
    bestFragments = Chem.GetMolFrags(  Chem.FragmentOnBonds(frag, [bestSplitBond], addDummies=False) , asMols=True) 
    #print( [ Chem.CanonSmiles(Chem.MolToSmiles(frag)) for frag in bestFragments ])
    return bestFragments

def isValidFrag(frag, newCanonSmiles, allowLoneAtoms = False, maxFragSize = 5, maxPathLength = 3):
    '''Defines whether a fragment is accepted or not. It is accepeted if a fragment with the same SMILES code has already been accepted, if it can't be broken further, if breaking it in the next step would produce a lone atom,  if it has fewer atoms than maxFragSize or if all non-ring paths are less than maxPathLength'''
    fragSmiles = Chem.CanonSmiles( Chem.MolToSmiles(frag) )
    if fragSmiles in canonSmilesLookup  or fragSmiles in newCanonSmiles:
        #print("Known fragment: ", Chem.CanonSmiles( Chem.MolToSmiles(frag) ) )
        return True
    if frag.GetNumAtoms() <= maxFragSize:
        #print("Fewer than five atoms remaining: ", Chem.CanonSmiles( Chem.MolToSmiles(frag) ) )
        return True
    if len( getBreakableBonds(frag) )  == 0:
        print("No more breakable bonds left: ", Chem.CanonSmiles( Chem.MolToSmiles(frag) ) )
        return True
    if allowLoneAtoms == False:
        trialBisect = bisectFrag(frag)
        if trialBisect[0].GetNumAtoms() < 2 or trialBisect[1].GetNumAtoms()<2:
            #print("Next move would produce a lone atom", Chem.CanonSmiles( Chem.MolToSmiles(frag) ) ) 
            return True
    
    '''
    #break into ring and non-ring components        
    bis =   frag.GetSubstructMatches(Chem.MolFromSmarts('[R][!R]'))
    ringtononringBonds = [frag.GetBondBetweenAtoms(x,y)  for x,y in bis]
    longestPath = 0
    if len(ringtononringBonds) > 0:   
        ringAndSideChains = Chem.GetMolFrags(  Chem.FragmentOnBonds(frag , ringtononringBonds , addDummies=False) , asMols=True)
        for chain in ringAndSideChains:
            if isRing(chain):
                continue
            if pathLength(chain) > longestPath:
                longestPath = pathLength(chain)
    elif isRing(frag) == False:
        longestPath = pathLength(frag)
    if longestPath < maxPathLength: #the longest path is either shorter than the cut-off or the structure is a ring with short side chains
        return True
    '''
    return False

def multipleBisect(targetMol):
    processingSet = [ targetMol ]
    acceptedNewSMILES = []
    completeSet = []
    while len(processingSet) > 0:
        currentMol = processingSet.pop(0)
        trialSMILES = Chem.CanonSmiles( Chem.MolToSmiles(currentMol) ) 
        if isValidFrag(currentMol, acceptedNewSMILES)   :
            completeSet.append(currentMol)
            if trialSMILES not in acceptedNewSMILES:
                acceptedNewSMILES.append(trialSMILES)
        else:
            newFrags = bisectFrag(currentMol)
            print("Split ", trialSMILES, "to", Chem.CanonSmiles( Chem.MolToSmiles(newFrags[0])) , Chem.CanonSmiles( Chem.MolToSmiles(newFrags[1])  ))
            processingSet.append(newFrags[0])
            processingSet.append(newFrags[1])
    return completeSet

targetMol  = Chem.MolFromSmiles(targetMolSMILES)
for atom in targetMol.GetAtoms():
    atom.SetIntProp("index0", atom.GetIdx())
    

if targetMolSMILES in canonSmilesLookup:
    print(targetMolSMILES, " is already defined ", canonSmilesLookup[targetMolSMILES]  )
    fragmentsToBeads(targetMol, [targetMol])  
elif method=="BRICS":
    molBroken = Chem.BRICS.BreakBRICSBonds(targetMol)
    molFragments = Chem.GetMolFrags(molBroken, asMols=True)
    fragmentsToBeads(targetMol, molFragments)   
elif method=="PMFP":
    #find all non-ring bonds - these can be broken
    bis = targetMol.GetSubstructMatches(Chem.MolFromSmarts('[!R][!R]')) + targetMol.GetSubstructMatches(Chem.MolFromSmarts('[R][!R]'))
    breakableBonds =   [targetMol.GetBondBetweenAtoms(x,y).GetIdx() for x,y in bis]
    #print(bis)
    #print(breakableBonds)
    if len(breakableBonds) == 0:
        print("No non-ring bonds exist, structure cannot be decomposed, using single-bead model")
        fragmentsToBeads(targetMol, [targetMol])
    quit()
    #define a set of costs - the best mapping is the one that minimises the total cost.
    #newBeadTypeCost should be less than perBeadCost to encourage the use of larger beads where possible

    
    perBeadCost = 5
    newBeadTypeCost = 4
    newPMFCost = 10
    ringAtomWeight = 1
    loneAtomPenalty = 5 #extra penalty for lone C atoms
    
    #initial best is just the entire molecule - cost for one bead, plus cost for one type of new bead, + new pmf cost + potential oversize penalty
    bestCost =   perBeadCost  + newBeadTypeCost + newPMFCost
    bestFragments = [targetMol]
    bestCost += fragmentCostFunction(targetMol)
    
    maxBreaks = 10
    numBreaks = min(maxBreaks, len(breakableBonds))
    currentBeadNum = 1
    breakableBondPermutations = itertools.chain(*(list(itertools.combinations(breakableBonds, i + 1)) for i in range(numBreaks)))
    improvedWithExtra = False
    for bondSet in breakableBondPermutations:
        totalCost = 0
        molFragments = Chem.GetMolFrags(  Chem.FragmentOnBonds(targetMol, bondSet, addDummies=False) , asMols=True)
        allFragmentSMILES = [ Chem.CanonSmiles(Chem.MolToSmiles(frag)) for frag in molFragments ]
        numBeads = len(molFragments)
        if numBeads > currentBeadNum:
            currentBeadNum = numBeads
            print("Starting sets of ", currentBeadNum, "beads. Best so far: ",  [ Chem.CanonSmiles(Chem.MolToSmiles(frag)) for frag in bestFragments ]  , "with total cost ", bestCost)
        totalCost += perBeadCost * numBeads
        numBeadTypes = 0
        numNewBeadTypes = 0
        usedSMILES = []
        #print(molFragments)
        for frag in molFragments:
            fragCS = Chem.CanonSmiles(Chem.MolToSmiles(frag))
            totalCost += fragmentCostFunction(frag)
            if fragCS not in usedSMILES: #if the SMILES code is not already used for this molecule, apply a penalty. if it needs a new PMF generated, apply an extra penalty
                totalCost += newBeadTypeCost
                usedSMILES.append(fragCS)
                #print("Adding ", fragCS, "to SMILES set for this target")
                if fragCS not in canonSmilesLookup:
                     totalNewPMFCost = newPMFCost #set the base cost for generating a new PMF 
                     #if a PMF has a lot of breakable bonds its likely to be flexible and not as useful a building block
                     #print("New target PMF has: ", len(breakableBondsPMF), "breakable bonds", fragCS)
                     #surplusPMFBonds = max(0 , len(breakableBondsPMF) - 4 )
                     totalCost += totalNewPMFCost #+ 2 * surplusPMFBonds
                     #print("Adding ", fragCS, "to PMF generation set for this target")
        if totalCost < bestCost:
            bestCost = totalCost
            bestFragments = molFragments
    print("Best obtained set:",  [ Chem.CanonSmiles(Chem.MolToSmiles(frag)) for frag in bestFragments ]  , "with total cost ", bestCost)
    fragmentsToBeads(targetMol, bestFragments)
elif method=="EqualParts":
    bestFragments = multipleBisect(targetMol)
    print("Best obtained set:",  [ Chem.CanonSmiles(Chem.MolToSmiles(frag)) for frag in bestFragments ] )
    fragmentsToBeads(targetMol, bestFragments)
else:
    print("Molecule is not a direct match and no known fragmentation method provided. Check input and try again")
