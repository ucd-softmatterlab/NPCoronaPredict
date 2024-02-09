Please note: Some of the information below is outdated, for more up to date instructions please refer to the manuals included in the docs folder. A typical installation runs like:


mkdir UnitedAtom
git clone https://github.com/ucd-softmatterlab/NPCoronaPredict.git
cd UnitedAtom
sudo apt install libboost-all-dev libboost-filesystem-dev wget perl parallel  
sudo apt install build-essential libsqlite3-dev libpng-dev libjpeg-dev
sudo apt install python3 python3-pip
pip3 install numpy scipy matplotlib argparse Bio
make clean
make
sudo apt install qtchooser qtbase5-dev qt5-qmake 
cd UAQuickRun
qmake -qt=qt5
make
cd ../NPDesigner
qmake -qt=qt5
make
cd ../

This will build the required executables including the grahpical interfaces, and ensure the majority of the Python libraries needed are installed.

####UNITEDATOM######

How to build:
Install GCC and Boost with development headers (more detailed instructions to come). Also install python3, numpy, scipy and matplotlib as some of the tools use these.

Once everything is installed run
make clean; make
This will bring up a number of warnings, primarily due to unused parameters which are left as hooks for future versions of the code. These are generally ok to ignore. 
Boost errors usually indicate that it's not finding the correct version of Boost, which can happen if you have multiple versions installed. Try editing the Makefile to add the
-Lpath/to/recent/boost
option.

Please note that if you have a slightly older version of boost there may be a compilation error:
src/main.cpp:821:110: error: ‘boost::filesystem::copy_options’ has not been declared
caused by the line:
boost::filesystem::copy_file(configFileIn, config.m_outputDirectory+"/"+configFileIn, boost::filesystem::copy_options::overwrite_existing);
To fix this either upgrade boost, comment this line out, or change it to :
boost::filesystem::copy_file(configFileIn, config.m_outputDirectory+"/"+configFileIn, boost::filesystem::copy_option::overwrite_if_exists);


How to run:
./UnitedAtom --config-file=name_of_config_file.config

The expected output is a set of .map files in the results folder named in the config file. Each map file is named:
PDBID_(NP radius in nm)_(zeta potential in microvolts)_(angle of rotation of NP in degrees).map and contains data in the form
phi theta energy error

A sample config file is provided and most options are self-explanatory. Some require more explanation as provided below:

The geometry of the NP is defined by the radius and the np-type variable. The options for this are given by:
1 = sphere (Planar PMF)
2 = cylinder  (Planar PMF) 
3 = cube (Planar PMF)
4 = SWCNT ( Cylindrical PMF of diameter 1.5nm , finite thickness tube for vdW)
5 = MWCNT (Cylindrical PMF of diameter 1.5nm, full cylinder for vdW)
For a spherical NP, no rotation of the NP is performed and the output names do not have the final _ANGLE, e.g. they are 1AX8_5_0.map

If the recalculate-zp option is set to 1 then the input zeta potential is taken to be the zeta potential for a particle of radius 1 nm in a solution with debye length 1 nm and Bjerrum length 1 nm. The zeta potential applied to
a given particle is then adjusted to ensure that the surface charge density of the NP is a constant regardless of the radius. 

Setting the calculate-mfpt flag to 1 makes the code estimate the mean-first passage time for each orientation and save this to _mfpt.map files. This is very slow and is mostly of use if you want to be more rigorous for corona studies. 
This output must also be corrected for the diffusion coefficient of the protein in question.

Temperature: This can be set in the config file using the temperature option. Output energies are given in units of kb * T=300 K for backwards compatability and also in kJ/mol for ease of comparison.

save-potentials : if set to 1, enables saving the protein-NP potential for each orientation. This takes up a lot of space and is buggy due to multithreading, its mostly useful for debugging


If np-target is set to a folder or specific .np file this enables the multi-component NP mode for NP complexes. This is relatively stable for core-shell models but still has some issues with brushes
Multi-component NPs are defined by a .np file in which any line starting with a # is a header and all other lines define component NPs in the complex.



NPs with non-spherical structures require a bounding radius to be defined. By default, this is estimated by drawing a sphere around the entire complex and taking the minimum radius
which encapsulates all the NPs. For low-density brushes this can be too generous and cause the protein to try to bind to empty space. Manually setting the bounding-radius parameter might lead to better results but may also lead to overlaps and invalid energies.


The NP file format is as follows:

x,y,z,radius,zeta,coreScaleFactor,surfaceScaleFactor,shape,HamakerFile,PMFDirectory,PMFCutoff,Correction Type

where x,y,z are the co-ordinates in nm of the component centre, radius in NM is the radius, zeta in V is the zeta potential. Core scale and Surface scale enable scaling of these components by arbitrary constants.This allows for subtraction of a small NP from a large one to produce a shell.
HamakerFile and PMF directory have the same form as in the regular config file.
PMFCutoff is given in nm and is used for mapping input PMFs and should match the LJ cutoff distance (not the final distance in the PMF). Typically, this is 1.0 for GAFF PMFs (Stockholm) and 1.2 for CHARMM PMFs (UCD).
Correction type defines how the PMF is re-mapped to the target geometry. Correction type 1 is plane-to-sphere, correction type 2 is 1.5nm diameter tube to cylinder, type 3 is no correction, type 4 and 5 is 1.5nm diameter tube to tube
ANy other value or 0 leaves applies zero correction - this can be used for beads with fixed radius
The PMFCutoff value is also used for calculating the Hamaker lens potential so make sure this is set to be a sensible value, else beads overlap and divergences happen.







###BuildCoronaParams.py###
Converts UA output to a file suitable for loading into CoronaKMC.

How to run: (requires python 3)
python BuildCoronaParams-P3.py -r [RADIUS] -z [ZETA potential] -f [folder containing UA heatmaps] -s [1 or 2 for sphere or cylinder] -p [protein definition file] -c [coordinate (in PDB form) folder)

This requires a file named ProteinConcs.csv (or whatever you pass using the -p argument)  of the format
protein1name,protein1conc
protein2name,protein2conc

where proteinXname is the name of the protein used by UnitedAtom (e.g. 1AX8)  and proteinXconc is the molar concentration of that protein. The script scans the target folder for all proteins named in the ProteinConcs file for the target radius and zeta potential, then creates a file containing the binding parameters for each orientation.

The script also requires the PDB files used by UA for each protein, which it gets from the folder given by the -c argument. 

For example, if the protein concentration file is named ProteinConcTest.csv and contains
PROTA,0.0001
and the input command is
python BuildCoronaParams.py -r 5 -z 0 -f results -s 1 -p ProteinConcTest.csv -f pdbs
then the following happens:
1) The script looks in the folder "results" for a file called "PROTA_5_0.map" and loads this in.
2) The script looks in the folder "pdbs" for a file called "PROTA.pdb" and loads this in.
3) For each orientation in PROTA_5_0.map, the script calculates a set of binding coefficients based on the UA heatmap and geometry of the protein.

The shape specification should be used together with the correct folder for UA heatmaps for spherical or cylindrical geometries, e.g.
python BuildCoronaParams.py -r 5 -z 0 -f results_anatase_sphere -s 1
or
python BuildCoronaParams.py -r 5 -z 0 -f results_anatase_cylinder -s 2
as there is presently no way to distinguish between heatmaps for different geometeries. The -s argument mostly exists to tell the script how to calculate binding areas.

The expected output from this script is a file containing an entry for every orientation of every protein included in ProteinConcs.csv . As the UA default is 3600 orientations this can be quite a long file.
The output is of the form
proteinname concentration radius kon koff ebind area

with one entry for each orientation of each protein. The concentrations for orientations of a given protein are weighted by \sin \theta and normalised such that the sum over all orientations is equal to the total concentration set in the input file
The radius is calculated by projecting the protein onto the surface of the NP, determining the area of this projection, and calculating the radius of a sphere required to produce the same area.
kon is estimated from diffusion theory for sphere-sphere collisions.
koff is estimated from kon/koff = exp(-ebind)
ebind is taken from the UA heatmap for that orientation.
area is the projected area.  

####CoronaKMC####
Coarse-grained Kinetic Monte Carlo simulation of corona evolution.

How to run: (requires python 3)
python CoronaKMC-P3.py  -r [RADIUS] -p [protein file] -f [File ID - number to add to filename of output]

The input specified by the -p argument is a file with one entry per line of the form:
proteinname concentration radius kon koff ebind area

This output is produced by BuildCoronaParams.py. The radius given as input must be set to match whatever is specified when calculating the protein input file.
There are a number of further optional parameters:
-s: Defaults to 1 for a spherical NP. If set to 2 this produces a cylinder instead.
-d: If set to 1 then proteins on the surface of the NP are allowed to diffuse. This is very slow, but necessary for irreversible adsorption to ensure surface restructuring.
-c: If set to 1 then all the output is suppressed apart from the total surface coverage and total number of proteins. At the end these are used to calculate a protein with radius and binding energy that would produce these values.
-m: If set to 1 then the meanfield approximation is enabled. 

The output is the number of adsorbed proteins at each timestep. If multiple proteins in the input file have the same name then they are counted together in this output.
This enables the use of an input file containing multiple orientations of the same protein without producing extremely complex output, while still allowing for a difference in binding between different orientations.

It should be noted that if the meanfield approximation is enabled then there is an analytical solution available through matrix methods to check the results of output, but this is unlikely to be viable for large systems.
The steady-state in the meanfield approximation is given by:
N_i = (Area_NP/Binding_area )* (    conc_i * k_on_i/k_off_i)/( 1 + sum_j conc_j k_on_j/k_off_j)

####Credits ###


Contributors (Code and model)
David Power (UCD)
Ian Rouse (UCD)
Hender Lopez (UCD(
Vladimir Lobaskin (UCD)


Contributors (PMF input files)
Erik Brandt (Stockholm University)
Marzieh Saeedimasine (Stockholm University)
Alexander Lyubartsev (Stockholm University)
