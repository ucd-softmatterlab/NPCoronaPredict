#! /bin/bash

#OUT=results-3dbz
OUT=results-1n5u

RADIUS=[2.5,5,10,20,40,60,80,120,160,240,320,500,640,10000]
ZETA=[-0.05,0,0.05]
PDB=pdbs/crystal/1n5u.pdb

#RADIUS=50
#ZETA=0.0

./UnitedAtom --hamaker-file=hamaker/Au.dat --pmf-directory=surface/Au/FCC/100/sca --output-directory=${OUT} --pdb-target=${PDB} --nanoparticle-radius=${RADIUS} --zeta-potential=${ZETA} $@



