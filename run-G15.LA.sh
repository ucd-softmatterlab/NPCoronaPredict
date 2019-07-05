#! /bin/bash

OUT=results-G15.LA
RADIUS=7.95
ZETA=-0.02648

# Antithrombin-3
./UnitedAtom --hamaker-file=hamaker/Au.dat --pmf-directory=surface/Au/FCC/100/sca --output-directory=results-G15.LA --pdb-target=pdbs/G15.LA/P01008.pdb --nanoparticle-radius=${RADIUS} --zeta-potential=${ZETA}

# Prothrombin
./UnitedAtom --hamaker-file=hamaker/Au.dat --pmf-directory=surface/Au/FCC/100/sca --output-directory=results-G15.LA --pdb-target=pdbs/G15.LA/P00734.pdb --nanoparticle-radius=${RADIUS} --zeta-potential=${ZETA}

# Vitronectin
./UnitedAtom --hamaker-file=hamaker/Au.dat --pmf-directory=surface/Au/FCC/100/sca --output-directory=results-G15.LA --pdb-target=pdbs/G15.LA/P04004.pdb --nanoparticle-radius=${RADIUS} --zeta-potential=${ZETA}

# Apolipopretein E
./UnitedAtom --hamaker-file=hamaker/Au.dat --pmf-directory=surface/Au/FCC/100/sca --output-directory=results-G15.LA --pdb-target=pdbs/G15.LA/P02649.pdb --nanoparticle-radius=${RADIUS} --zeta-potential=${ZETA}

# Histidine-rich Glycoprotein
./UnitedAtom --hamaker-file=hamaker/Au.dat --pmf-directory=surface/Au/FCC/100/sca --output-directory=results-G15.LA --pdb-target=pdbs/G15.LA/P04196.pdb --nanoparticle-radius=${RADIUS} --zeta-potential=${ZETA}

