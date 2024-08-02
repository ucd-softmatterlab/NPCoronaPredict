import numpy as np
import os

def SearchDirectory(path):
    files = os.listdir(path)
    pdbs  = []
    for handle in files:
        abspath = os.path.join(path, handle)
        if os.path.isdir(abspath):
            pdbs += SearchDirectory(abspath)
        elif abspath[-4:] == ".dat":
            pdbs.append(abspath)
    return pdbs


pmfFolder = "AlFCC111"

pmfSet  = SearchDirectory(pmfFolder)
for pmf in pmfSet:
    print(pmf)
    try:
        pmfIn = np.genfromtxt(pmf)
        #print(pmfIn)
        np.savetxt(pmf, pmfIn, delimiter=",", fmt='%10.6f', header="#h[nm], U[kJ/mol]")
    except:
        print("Failed for ", pmf)
