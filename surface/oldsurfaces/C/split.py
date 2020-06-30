#! /usr/bin/python

import os
from sys import argv

# Read File
handle = open(argv[1])
lines = handle.read().splitlines()
handle.close()

# Split File
Materials = {}
for line in lines:
    if line[0] == '#':
        continue
    args = filter(None, line.split('\t'))
    material  = args[0]
    aminoAcid = args[1]
    distance  = args[2]
    energy    = args[3]
    if material not in Materials:
        Materials[material] = {}
    if aminoAcid not in Materials[material]:
        Materials[material][aminoAcid] = []
    Materials[material][aminoAcid].append((float(distance), float(energy)))

# Write New Files
for material, AminoAcids in Materials.iteritems():
    if not os.path.isdir(material):
        os.mkdir(material)
    for aminoAcid, data in AminoAcids.iteritems():
        handle = open(material + "/" + aminoAcid + ".dat", 'w')
        handle.write("# z (nm)   (kJ/mol)\n")
        for distance, energy in data:
            handle.write("{:8.4f}{:11.5f}\n".format(distance, energy))
        handle.close()
