####NPCoronaPredict####
NPCoronaPredict is a package to allow the prediction of the corona content of a nanoparticle immersed in a mixture of biomolecules via a multiscale coarse-grained approach.

Please refer to the manual included in the /docs folder for detailed information about how to run the programs.

For works published using results generated using scripts in this repository, please cite at least:

    NPCoronaPredict: A Computational Pipeline for the Prediction of the Nanoparticle–Biomolecule Corona
    Ian Rouse, David Power, Julia Subbotina, and Vladimir Lobaskin
    Journal of Chemical Information and Modeling Article ASAP
    DOI: 10.1021/acs.jcim.4c00434 


For UnitedAtom the original publication for the current methodology is:
    
    David Power et al 2019 Modelling Simul. Mater. Sci. Eng. 27 084003, DOI: 10.1088/1361-651X/ab3b6e

For CoronaKMC cite:

    A hard-sphere model of protein corona formation on spherical and cylindrical nanoparticles
    Rouse, Ian et al.
    Biophysical Journal, Volume 120, Issue 20, 4457 - 4471


####Quick installation####

A typical installation runs like:


mkdir NPCoronaPredict; 
cd NPCoronaPredict; 
git clone https://github.com/ucd-softmatterlab/NPCoronaPredict.git;
cd NPCoronaPredict;
sudo apt install libboost-all-dev libboost-filesystem-dev wget perl parallel;
sudo apt install build-essential libsqlite3-dev libpng-dev libjpeg-dev;
sudo apt install python3 python3-pip;
pip3 install numpy===1.26.3 scipy matplotlib argparse Bio;
make clean;
make;
sudo apt install qtchooser qtbase5-dev qt5-qmake;
cd NPCoronaPredict-GUI;
qmake -qt=qt5;
make;
cd ../NPDesigner;
qmake -qt=qt5;
make;
cd ../

This will build the required executables including the graphical interfaces, and ensure the majority of the Python libraries needed are installed.

If you are not an administrator, you won't be able to run anything requiring sudo so you will need to make sure your system admin has installed Python, Boost and QT. In this case, you can just run:

mkdir NPCoronaPredict;
cd NPCoronaPredict;
git clone https://github.com/ucd-softmatterlab/NPCoronaPredict.git;
cd NPCoronaPredict;
pip3 install numpy===1.26.3 scipy matplotlib argparse Bio;
make clean;
make;
cd NPCoronaPredict-GUI;
qmake -qt=qt5;
make;
cd ../NPDesigner;
qmake -qt=qt5;
make;
cd ../



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
