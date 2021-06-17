from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.Residue  import DisorderedResidue
import numpy as np
import copy
import os
parser = PDBParser(PERMISSIVE=1)
io = PDBIO()


#targetPH = 7
#phRange = [1, 2, 3, 3.79,3.81, 4.49,4.51, 6.49,6.51, 7,  8.99,9.01, 10.49,10.51, 14]
phRange = [2,3]
proteinSet = [ "1bkz"]

outputFolder = "proteinpHStructuresBaseVals"

#if this is set to 0 then residues are set to their most likely form, i.e., protonated if pH < pKa and deprotonated HIS set to HIE. 
#if it is non-zero we include both types of residue and use the occupancy flag to set their relative weight, UA then weights
#the input potentials according to this weight. This results in larger files but smoother behaviour, especially for pH \approx pKa
residueInterpolation = 1
hidFraction = 0.2 #HIE (epsilon-protonated) is favoured over HID (delta-protonated) by a ratio of 4:1 in solution

#names for the charge variants
#this table is defined as: standard name, deprotonated name, protonated name
#HIS is a special case, here the deprotonated is HIE and protonated is HIP. HID is added in using the hidFraction later.
aaData = [
["ALA", "ALA", "ALA"],
["ARG", "ARG", "ARG"],
["ASN", "ASN", "ASN"],
["ASP", "ASP", "ASP"],
["CYS", "CYM", "CYS"],
["GLN", "GLN", "GLN"],
["GLU", "GLU", "GAN"],
["HIS", "HIE", "HIP"],
["ILE", "ILE", "ILE"],
["LEU", "LEU", "LEU"],
["LYS", "LYS", "LYS"],
["MET", "MET", "MET"],
["PHE", "PHE", "PHE"],
["SER", "SER", "SER"],
["THR", "THR", "THR"],
["TRP", "TRP", "TRP"],
["TYR", "TYR", "TYR"],
["VAL", "VAL", "VAL"]
]
aaDataIndexDict = {}
for i in range(len(aaData)):
    aaDataIndexDict[aaData[i][0]] = i
#print aaDataIndexDict

#fileStem = "1bkz"

#proteinSet = [ "1bkz" , "prot2"]

#print [ [ [ph, fileStem] for fileStem in proteinSet    ]    for ph in phRange     ]
inputSet = []
for ph in phRange:
    for protein in proteinSet:
        inputSet.append( (ph,protein))


for targetPH,fileStem in inputSet:
    print "Processing" , fileStem, "at pH ", targetPH
    if not os.path.isfile(fileStem+".pka"):
        print "PROPKA file not found for ", fileStem, "attempting to run propka"
        os.system("propka3 "+fileStem+".pdb")
    try:
        propkaFile = open(fileStem+".pka","r")
        propkaFileData = propkaFile.readlines()
        propkaFile.close()
    except:
        print "PROPKA file could not be read for ", fileStem
        #os.system("propka3 "+fileStem+".pdb")
        continue
    try:
        structureData = parser.get_structure(fileStem+"-struct" , fileStem+".pdb")
    except:
        print "Structure data for ", fileStem, " not found or not readable"
        continue
    summaryStart = 0
    newIDVal = 0
    for line in propkaFileData:
        if line == "SUMMARY OF THIS PREDICTION\n":
             summaryStart = 1
             continue
        if summaryStart == 1 and line == "--------------------------------------------------------------------------------------------------------\n":
             break
        if summaryStart == 1:
             if (line.strip())[:5] != "Group":
                rowData =  line.strip().split()
                #rowData: 0 = AA tag, 1 = residue number, 2 = chain ID, 3 = adjusted  pKa val
                #if pH < pKa the side chain is protonated
                residuepKa = float(rowData[3])
                if rowData[0] not in aaDataIndexDict:
                    print "Unrecognised residue ", rowData[0]
                    continue
                #print rowData
                #print "Deprotonated fraction" ,  np.exp(targetPH)/( np.exp(targetPH) + np.exp( float(rowData[3])))
                for model in structureData:
                     for residue in model.get_residues():
                        if newIDVal < residue.get_id()[1]:
                            newIDVal = residue.get_id()[1]
                     newIDVal = newIDVal + 1
                     if residueInterpolation == 0:
                          if targetPH <  residuepKa:
                              #print rowData[0], " at " , rowData[1] , " chain " , rowData[2] , "is protonated, new name ", aaData[   aaDataIndexDict[rowData[0]]][2]
                               newName  = aaData[   aaDataIndexDict[rowData[0]]][2]
                          else:
                               #print rowData[0], " at " , rowData[1] , " chain " , rowData[2] , "is deprotonated, new name ", aaData[   aaDataIndexDict[rowData[0]]][1]
                               newName  = aaData[   aaDataIndexDict[rowData[0]]][1]
                          model[rowData[2]][int(rowData[1])].resname = newName
                     else:
                          deprotonatedFrac = np.exp(targetPH)/( np.exp(targetPH) + np.exp(  residuepKa ))
                          protonatedFrac = 1 - deprotonatedFrac
                          targetResidue = model[rowData[2]][int(rowData[1])]
                          hisFactor = 1
                          if rowData[0] == "HIS": #add in the HID form of histidine using the weighting defined earlier
                               hisFactor = 1-hidFraction
                               print "HIE: ", deprotonatedFrac*hisFactor, " HID ", deprotonatedFrac*hidFraction, "HIP" , protonatedFrac
                               hidCopy = copy.deepcopy(targetResidue)
                               hidCopy.resname = "HID"
                               hidCopy.id =  ( hidCopy.get_id()[0],  newIDVal,hidCopy.get_id()[2])
                               newIDVal += 1
                               for atom in hidCopy.get_atoms():
                                   atom.set_occupancy(atom.get_occupancy() * deprotonatedFrac * hidFraction   )
                               model[rowData[2]].add(hidCopy)
                          #set the initially present residue to the deprotonated form
                          #targetResidue = model[rowData[2]][int(rowData[1])]
                          protonatedResidue = copy.deepcopy(targetResidue)
                          targetResidue.resname = aaData[   aaDataIndexDict[rowData[0]]][1]
                          for atom in targetResidue.get_atoms():
                              atom.set_occupancy(atom.get_occupancy() * deprotonatedFrac *hisFactor  )
                          protonatedResidue.resname = aaData[   aaDataIndexDict[rowData[0]]][2]
                          protonatedResidue.id = ( protonatedResidue.get_id()[0],  newIDVal, protonatedResidue.get_id()[2])
                          newIDVal += 1
                          for atom in protonatedResidue.get_atoms():
                              atom.set_occupancy(atom.get_occupancy() * protonatedFrac   )
                          model[rowData[2]].add(protonatedResidue)
    io.set_structure(structureData)
    io.save(outputFolder+"/"+fileStem+"-pH"+str(targetPH).replace(".","p")+"model"+str(residueInterpolation)+".pdb")
