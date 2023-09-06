'''Decompose a target molecule defined by a SMILES code into fragments using the BRICS algorithm (or others), then matches these fragments to pre-existing fragments where possible
Currently implemented algorithms (chose one using the -M option):

BRICS: Uses the BRICS method to split molecule into fragments corresponding to reasonable chemical reactions. Fast, reliable, but likely to need generation of many PMFs and produces stray atoms.
Exhaustive: Considers all possible bond breakages and finds the global minimum according to a defined score function (currently, fragments of a certain size without too many non-ring bonds). This can be very slow and potentially run out of memory, use with caution. 
EqualParts: Recursive bisection by breakable bonds, until each fragment is under a certain size, cannot be decomposed further, or corresponds to a known PMF. This tends to produce fewer lone atoms but may need more PMFs generated.
ForwardsMatching: For each bead template, gets all possible matches, identifies the atom with the fewest options and selects the largest applicable template (moving to smaller templates if this setp fails), removes the corresponding section from the molecule and iterates over the new fragments until everything is mapped to a known fragment. Not guaranteed to work, especially if there are complex aromatic sections.
'''

from rdkit import Chem
from rdkit.Chem import BRICS, Draw
from rdkit.Chem import AllChem
from rdkit import RDLogger
import numpy as np
import datetime
import os
import argparse
import itertools



#Preparation
RDLogger.DisableLog('rdApp.*') 
dummyAtom = Chem.MolFromSmiles("*")

dateString =  datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
shortDate = datetime.datetime.now().strftime("%d-%b-%y")
missingBeadString = "!!!"





def removeDummySMILES(smilesIn):
    cin = Chem.MolFromSmiles(smilesIn)
    return Chem.MolToSmiles( Chem.RemoveHs(  Chem.AllChem.ReplaceSubstructs(cin, dummyAtom, Chem.MolFromSmiles("[H]"), True)[0] ))


def getChargeFromSmiles(smilesIn):
    return Chem.GetFormalCharge( Chem.MolFromSmiles(smilesIn) )

def getFragmentSMILES(frag):
    return Chem.CanonSmiles(Chem.MolToSmiles( Chem.RemoveHs(  Chem.AllChem.ReplaceSubstructs(frag, dummyAtom, Chem.MolFromSmiles("[H]"), True)[0] )) )


def beadsetToPDB(chemName, chemSmiles, beadset,outputFolder,atomLabel="CA", extraRemarks = []):
    '''Writes out a set of beads [beadID, CodeID, x, y, z, smiles, longname, elementlabel] to a pdb-like file for use in UnitedAtom'''
    beadOut = open(outputFolder+"/"+chemName+".pdb","w")
    beadOut.write("HEADER "+ '{:40.40}'.format(chemName)+shortDate+"\n")
    #beadOut.write("REMARK Bead ID:" + uaTargets[chem]+"\n")
    beadOut.write("TITLE   "+chemName+"\n")
    beadOut.write("REMARK SMILES: "+chemSmiles+"\n")
    beadOut.write("REMARK Generated: "+dateString+"\n")
    for remark in extraRemarks:
        beadOut.write("REMARK "+remark+"\n")
    seenBeadIDs = []
    for atom in beadset:
        if atom[0] not in seenBeadIDs:  
            beadOut.write("REMARK "+str(atom[0])+":"+atom[6]+":"+atom[1]+":"+atom[5]+"\n")
            seenBeadIDs.append(atom[0])
    atomNum = 1
    for atom in beadset:
        beadOut.write("ATOM  "+str(atom[0]).rjust(5)+" "+" "+atomLabel+" "+" "+atom[1]+" A"+str(atomNum).rjust(4)+"    {:8.3}{:8.3}{:8.3}  1.00  0.00".format(atom[2],atom[3],atom[4])+atom[7].rjust(12)+"\n")
        atomNum += 1
    beadOut.write("END\n")
    beadOut.close()
    print("Written output file to ", outputFolder+"/"+chemName+".pdb")


def fragmentsToBeads(targetMol, molFragments,outputDir,canonSmilesLookup , usedMethod="DirectMatch"):
    targetMolH = Chem.AddHs( targetMol)
    embedRes = AllChem.EmbedMolecule(targetMolH)
    if embedRes == -1:
        AllChem.EmbedMolecule(targetMolH, maxAttempts=5000)
    positionList = np.array(targetMolH.GetConformer().GetPositions() )
    numMissingFragments = 0
    beadID = 0
    numAtomsFound = 0
    pmfpredOutList = []
    finalBeadSet = [] #set of CG beads
    fullBeadSet = [] #set of atoms with CG bead IDs
    reverseBeadLookup = {}
    hit_bonds = []
    highlightAtomSet = []
    for frag in molFragments:
        fragSMILES = getFragmentSMILES(frag)
        print("Using fragment",   fragSMILES )
        beadCode = missingBeadString
        beadName = targetChemName+"New"+str(numMissingFragments)
        beadLoc = np.array([0,0,0])
        totalBeadMass = 0.00001
        hit_ats = []
        
        for atom in frag.GetAtoms():
            #print(atom.GetIdx())
            try: 
                atomBaseIndex = atom.GetIntProp("index0")
                numAtomsFound +=1
                reverseBeadLookup[atomBaseIndex] = beadID
                hit_ats.append(atomBaseIndex)
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
        finalBeadSet.append( [beadID, beadCode, beadLoc[0], beadLoc[1], beadLoc[2]  , fragSMILES, beadName, "C" ] )
        for atom in frag.GetAtoms():
            try:
                atomLoc = positionList[   atom.GetIntProp("index0") ]
                fullBeadSet.append(   [beadID, beadCode, atomLoc[0], atomLoc[1], atomLoc[2]  , fragSMILES, beadName,atom.GetSymbol() ] )
            except:
                skipped = 1
        beadID += 1
        #hit_ats = list(targetMol.GetSubstructMatch(frag))
        #hit_ats
        
        #print( hit_ats )
        if len(hit_ats) == 1:
            highlightAtomSet.append( hit_ats[0]  )
        for bond in frag.GetBonds():
            aid1 = hit_ats[bond.GetBeginAtomIdx()]
            aid2 = hit_ats[bond.GetEndAtomIdx()]
            hit_bonds.append(targetMol.GetBondBetweenAtoms(aid1,aid2).GetIdx())
    beadsetToPDB(targetChemName,targetMolSMILES, finalBeadSet,outputDir, extraRemarks = ["Method: " + usedMethod])
    
    for atom in targetMol.GetAtoms():
        try:
            cgBeadID = reverseBeadLookup[atom.GetIntProp("index0")]
            atom.SetProp("atomNote", str(cgBeadID) )
        except:
            skipped = 1
            

    print(highlightAtomSet)
    Chem.Draw.MolToFile( targetMol, outputDir+"_images/"+targetChemName+".png", size=(1200,1200), highlightAtoms=highlightAtomSet ,highlightBonds=hit_bonds )
    
    for atomID in range(numAtomsFound,len(positionList)):
        beadLocAll = positionList[atomID]
        fullBeadSet.append( [beadID, "HHH", beadLocAll[0], beadLocAll[1], beadLocAll[2]  , "H", "H" ,"H"]  )
    beadsetToPDB(targetChemName+"-full", targetMolSMILES, fullBeadSet, outputDir+"_mappings",atomLabel="CX")
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

def isValidFrag(frag, newCanonSmiles,canonSmilesLookup, minFragSize=2, maxFragSize = 5):
    '''Defines whether a fragment is accepted or not. It is accepted if a fragment with the same SMILES code has already been accepted, if it can't be broken further, if breaking it in the next step would produce a fragment of size under minFragSize or  if the fragment has fewer atoms than maxFragSize'''
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
    testMinSize = True
    if testMinSize == True:
        trialBisect = bisectFrag(frag)
        if trialBisect[0].GetNumAtoms() < minFragSize or trialBisect[1].GetNumAtoms()< minFragSize:
            #print("Next move would produce a lone atom", Chem.CanonSmiles( Chem.MolToSmiles(frag) ) ) 
            return True
    return False

def multipleBisect(targetMol,canonSmilesLookup):
    processingSet = [ targetMol ]
    acceptedNewSMILES = []
    completeSet = []
    while len(processingSet) > 0:
        currentMol = processingSet.pop(0)
        trialSMILES = Chem.CanonSmiles( Chem.MolToSmiles(currentMol) ) 
        if isValidFrag(currentMol, acceptedNewSMILES,canonSmilesLookup)   :
            completeSet.append(currentMol)
            if trialSMILES not in acceptedNewSMILES:
                acceptedNewSMILES.append(trialSMILES)
        else:
            newFrags = bisectFrag(currentMol)
            print("Split ", trialSMILES, "to", Chem.CanonSmiles( Chem.MolToSmiles(newFrags[0])) , Chem.CanonSmiles( Chem.MolToSmiles(newFrags[1])  ))
            processingSet.append(newFrags[0])
            processingSet.append(newFrags[1])
    return completeSet

def coarseGrainMolecule(targetChemName, targetMolSMILESIn, method,libraryFilename="FragmentLookup.csv",baseDir = "CGMolecules", allowNewPMFs=True):
    os.makedirs(baseDir,exist_ok=True)
    os.makedirs(baseDir+"_images",exist_ok=True)
    os.makedirs(baseDir+"_mappings",exist_ok=True)
    targetMolSMILES = Chem.CanonSmiles(targetMolSMILESIn)
    #load in a library of pre-defined fragments with the form beadID, smiles code,  UA 3-letter code
    canonSmilesLookup = {}
    libraryFile=open(libraryFilename,"r")
    for line in libraryFile:
       if line[0] == "#":
           continue
       lineTerms = line.split(",")
       entrySMILES = lineTerms[2].replace("<HASH>","#").replace("<COMMA>",",")
       canonSMILES = Chem.CanonSmiles(entrySMILES)
       if canonSMILES not in canonSmilesLookup:
           canonSmilesLookup[canonSMILES] = [ lineTerms[0], lineTerms[1] ]
    libraryFile.close()
    knownMethods = ["BRICS", "Exhaustive", "EqualParts", "ForwardsMatching"]
    
    targetMol  = Chem.MolFromSmiles(targetMolSMILES)
    for atom in targetMol.GetAtoms():
        atom.SetIntProp("index0", atom.GetIdx())
    if targetMolSMILES in canonSmilesLookup:
        print(targetMolSMILES, " is already defined ", canonSmilesLookup[targetMolSMILES]  )
        fragmentsToBeads(targetMol, [targetMol],baseDir,canonSmilesLookup, usedMethod="DirectMatch")  
        return True
    elif method=="BRICS":
        print("Warning: BRICS is not compatible with requirement that fragments match pre-existing PMFs")
        molBroken = Chem.BRICS.BreakBRICSBonds(targetMol)
        molFragments = Chem.GetMolFrags(molBroken, asMols=True)
        fragmentsToBeads(targetMol, molFragments,baseDir,canonSmilesLookup, usedMethod=method)   
        return True
    elif method=="Exhaustive":
        #find all non-ring bonds - these can be broken
        bis = targetMol.GetSubstructMatches(Chem.MolFromSmarts('[!R][!R]')) + targetMol.GetSubstructMatches(Chem.MolFromSmarts('[R][!R]'))
        breakableBonds =   [targetMol.GetBondBetweenAtoms(x,y).GetIdx() for x,y in bis]
        #print(bis)
        #print(breakableBonds)
        if len(breakableBonds) == 0:
            print("No non-ring bonds exist, structure cannot be decomposed, using single-bead model")
            fragmentsToBeads(targetMol, [targetMol],baseDir,canonSmilesLookup, usedMethod=method)
            return True
        #define a set of costs - the best mapping is the one that minimises the total cost.
        #newBeadTypeCost should be less than perBeadCost to encourage the use of larger beads where possible
        perBeadCost = 5
        newBeadTypeCost = 4
        newPMFCost = 10
        ringAtomWeight = 1
        loneAtomPenalty = 5 #extra penalty for lone C atoms
        if allowNewPMFs == False:
            newPMFCost = 1000000
            loneAtomPenalty = 1
        #initial best is just the entire molecule - cost for one bead, plus cost for one type of new bead, + new pmf cost + potential oversize penalty
        bestCost =   perBeadCost  + newBeadTypeCost + newPMFCost
        bestFragments = [targetMol]
        bestCost += fragmentCostFunction(targetMol)
    
        maxBreaks = 10
        numBreaks = min(maxBreaks, len(breakableBonds))
        currentBeadNum = 1
        if len(breakableBonds) > 10:
            print("Warning: this structure is too large for exhaustive mode to be used automatically. Partial fragmentation will be performed, the results will require further fragmentation.")
            numBreaks = 2
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
                if method=="KnownFrags":
                    totalCost += 0
                else:
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
        fragmentsToBeads(targetMol, bestFragments,baseDir,canonSmilesLookup, usedMethod=method)
        return True
    elif method=="EqualParts":
        bestFragments = multipleBisect(targetMol,canonSmilesLookup)
        print("Best obtained set:",  [ Chem.CanonSmiles(Chem.MolToSmiles(frag)) for frag in bestFragments ] )
        fragmentsToBeads(targetMol, bestFragments,baseDir,canonSmilesLookup, usedMethod=method)
        return True
    elif method=="ForwardsMatching":
        cslSort= sorted( canonSmilesLookup.keys(), key=len)
        cslSort.reverse()
        
        completedFragments = []
        processingFragments = [targetMol]
        failedFragments = []
        while len(processingFragments) > 0 :
            currentFrag = processingFragments.pop(0)
            print("Beginning fragment: ", Chem.MolToSmiles(currentFrag) )
            atomMatches = []
            for atom in currentFrag.GetAtoms():
                atomMatches.append( [] )
            for trialFragSmiles in cslSort:
                matchSet = currentFrag.GetSubstructMatches(Chem.MolFromSmiles(trialFragSmiles)) 
                for matchAtoms in matchSet:
                    #print(matchAtoms)
                    for matchAtomIndex in matchAtoms:
                        atomMatches[ matchAtomIndex ].append([ trialFragSmiles, matchAtoms])
            currentTarget = 0
            currentMatches = len( cslSort)
            possibleMatchSet = []
            for i in range(len(atomMatches)):
                atom = atomMatches[i]
                
                #print(atom)
                if len(atom) < currentMatches:
                    currentMatches = len(atom)
                    currentTarget = i
                    possibleMatchSet = atom 
                    #print("Switched to atom ", i, "with ", atom)
                if len(atom) == 0:
                    print("An atom exists with no potential matches for the current set of PMFs.")
                    return false
            currentFragSet = []
            successfulFragmentation = False
            for currentSmiles in possibleMatchSet:
                print("removing based on pattern",currentSmiles )
                removeAtoms = list( currentSmiles[1] )
                bondsToBreak = []
                for atomIdx in removeAtoms:
                    atom = currentFrag.GetAtomWithIdx(atomIdx)
                    for bond in atom.GetBonds():
                        otherAtom = bond.GetOtherAtomIdx(atomIdx)
                        #print(atomIdx, otherAtom)
                        if otherAtom not in removeAtoms:
                            bondsToBreak.append(bond.GetIdx())
                            #print("removing bond", bond, bond.GetIdx(), atomIdx, otherAtom)
                #print(bondsToBreak)
                try:
                    currentFragSet = Chem.GetMolFrags(  Chem.FragmentOnBonds(currentFrag, bondsToBreak, addDummies=False) , asMols=True)
                    successfulFragmentation = True
                    break
                except:
                    print("Current fragment", Chem.MolToSmiles(currentFrag) , "could not be separated using pattern", currentSmiles)
                    continue
                    #print(currentFragSet)
                    #print( [ Chem.CanonSmiles(Chem.MolToSmiles(frag)) for frag in currentFragSet ] )
            if successfulFragmentation == False:
                print("Could not process fragment")
                failedFragments.append(currentFrag)
            for frag in currentFragSet:
                    if Chem.CanonSmiles(Chem.MolToSmiles(frag)) in cslSort:
                       completedFragments.append(frag)
                       print(  Chem.CanonSmiles(Chem.MolToSmiles(frag)), "recognised")
                    else:
                       processingFragments.append(frag)
               # except:
               #     print("Current fragment", Chem.MolToSmiles(currentFrag) , "could not be separated")
               #     failedFragments.append(currentFrag)
        print("Successfully mapped:")
        for frag in completedFragments:
            print( Chem.MolToSmiles(frag) )
        if len(failedFragments) > 0:
            print("Not mapped:")
            for frag in failedFragments:
                print( Chem.MolToSmiles(frag) )      
        else:
            fragmentsToBeads(targetMol,completedFragments,baseDir,canonSmilesLookup, usedMethod=method) 
    else:
        print("Molecule is not a direct match and no known fragmentation method provided. Check input and try again")
        print("Allowed methods are: ", knownMethods)
        return False


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parameters for PMFPredictor")
    parser.add_argument("-m","--molecule", type=str,default="NewMolecule", help="Name for new molecle")
    parser.add_argument("-s","--smiles", type=str,default="CCCC",help="SMILES string for new molecule")
    parser.add_argument("-M","--method",type=str,default="EqualParts",help="Method used to decompose molecule. BRICS: use BRICS only.  Exhaustive: Considers all possible bond breakages and finds the global minimum by scoring (slow). EqualParts: Recursive bisection by breakable bonds, matching to existing PMFs where possible. ForwardsMatching: Maps atoms to pre-existing PMFs, not guaranteed to work. ")
    parser.add_argument("-k","--knownonly",type=int,default=0,help = "If non-zero, attempts to generate structures using only pre-existing PMFs. Only implemented for Exhaustive")
    args = parser.parse_args()
    targetChemName = args.molecule
    targetMolSMILES = Chem.CanonSmiles(args.smiles)
    method = args.method
    allowNewPMFs = True
    if args.knownonly != 0:
        allowNewPMFs = False
    coarseGrainMolecule(targetChemName, targetMolSMILES, method , allowNewPMFs = allowNewPMFs)

    


