#! /bin/bash

OUT=results-3dbz
ZETA=[-20,0,20]

./UnitedAtom --hamaker-file=hamaker/Au.dat --pmf-directory=surface/Au/FCC/100/sca --output-directory=${OUT} --pdb-target=pdbs/crystal/3dbz.pdb --nanoparticle-radius=[5,10,20,30,40,50,60,70,80,90,100,200,10000] --zeta-potential=${ZETA}


