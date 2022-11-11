''' 
Helper script for concatenating corona input files together to allow for multi-surface NPs (e.g. wulff crystals, mixed-density brushes)
'''
import numpy as np


#Define the surface components. Each patch here corresponds to the output of a BuildCoronaParams run (first entry in list) 
#and a weight (second entry). The normalised weight is used to scale kon so that small patches have smaller targets. 
#the weights should reflect the surface area of each patch type
'''
npPatches = [
["sio2_pristine_birch_h2o_ph47_40_-29.csv", 4]
["anatase_patch.csv",1]
]
'''
#would generate an NP with 4/5ths surface area silica and 1/5th anatase.
#note that you will get weird results if you have inconsistent protein inputs!
#if a protein entry is missing for a given patch it will behave as if it is completely nonbinding to that patch



outputFileName = "cg_corona_data/mergedpatches.csv"
outputFileHandle = open(outputFileName,"w")

npPatches = [
["cg_corona_data/sio2_pristine_birch_h2o_ph47_40_-29.csv", 4],
["cg_corona_data/anatase_patch.csv",1]
]

#Normalise
weightNorm = 0
for i in range(len(npPatches)):
    weightNorm += npPatches[i][1]
for i in range(len(npPatches)):
    patch = npPatches[i]
    patchWeight = patch[1] / weightNorm
    patchFile = open(patch[0],"r")
    for line in patchFile:
        lineTerms = line.strip().split()
        lineTerms[0] = lineTerms[0]+":patch"+str(i)
        lineTerms[3] = str( float(lineTerms[3])*patchWeight )
        outputFileHandle.write(" ".join(lineTerms)+"\n" )
    patchFile.close()
outputFileHandle.close()
