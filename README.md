####NPCoronaPredict####

####Quick installation####

A typical installation runs like:


mkdir NPCoronaPredict; 
cd NPCoronaPredict; 
git clone https://github.com/ucd-softmatterlab/NPCoronaPredict.git;
cd NPCoronaPredict;
sudo apt install libboost-all-dev libboost-filesystem-dev wget perl parallel;
sudo apt install build-essential libsqlite3-dev libpng-dev libjpeg-dev;
sudo apt install python3 python3-pip;
pip3 install numpy scipy matplotlib argparse Bio;
make clean;
make;
sudo apt install qtchooser qtbase5-dev qt5-qmake;
cd UAQuickRun;
qmake -qt=qt5;
make;
cd ../NPDesigner;
qmake -qt=qt5;
make;
cd ../

This will build the required executables including the grahpical interfaces, and ensure the majority of the Python libraries needed are installed.

####Contents####
Two main utilities are included in the NPCoronaPredict repository:

UnitedAtom is used to compute the adsorption affinity of large biomolecules to nanoparticles via a multiscale coarse-graining model.

CoronaKMC is used to predict corona abundances of a set of adsorbates given adsorption and desorption rate constants and serum abundances. 

The NPCoronaPredict.py script provides a unified pipeline for processing a set of potential adsorbates and running UnitedAtom and CoronaKMC together. In particular this uses the supplied BuildCoronaParams Python script to map UA adsorption energies to CoronaKMC rate constants.


####Credits ####


Contributors (Code and model)
David Power (UCD)
Ian Rouse (UCD)
Hender Lopez (UCD(
Vladimir Lobaskin (UCD)


Contributors (PMF input files)
Erik Brandt (Stockholm University)
Marzieh Saeedimasine (Stockholm University)
Alexander Lyubartsev (Stockholm University)
Julia Subbotina (UCD)
Parinaz Mosadeggi (UCD)
