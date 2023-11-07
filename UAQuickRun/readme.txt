###UAQuickRun README file, last updated 19/10/2023 by Ian Rouse, ian.rouse@ucd.ie



----------------------
Introduction
----------------------
UAQuickRun is a simple graphical interface for running UnitedAtom for a single protein - NP pair and visualising the results. Its functionality is very limited by design compared to the full command-line program to make it easy to use for beginners without overwhelming them with too many options. 



----------------------
Compilation
----------------------
Dependencies:
UAQuickRun needs a working qt5 installation. The quick and easy way to get this on Linux is:
sudo apt install qtchooser qtbase5-dev qt5-qmake 

Compiling UAQuickRun:
The source code is available as part of the NPCoronaPredict package and requires gcc and Qt5 for installation. To compile the program, go into the UAQuickRun folder in a terminal window and run the commands:
qmake -qt=qt5
make

If all goes well there will be a few compiler warnings about unused variables, and then a UAQuickRun executable will be generated in the UAQuickRun folder. You should be able to run this like any other program. Depending on what other packages you have installed, you might need to do further work to get a working QT runtime - check the usual sources (google, stackoverflow, etc). 



----------------------
Basics
----------------------

When you   run the program, it automatically checks the directory one level up for a UA executable and a MaterialSet.csv file - this is so that it knows its in the right place and knows what types of nanomaterial are currently available. If it can't find UA, click the "Find UA Folder" button on the "UA Setup" tab and navigate to the folder containing UA. You might then also need to click the "Load Material Set" button to find a MaterialSet.csv file. 

Once it all works, there are two main tabs, UA Setup and View Results.


----------------------
UA Setup
----------------------
This tab allows you to do a basic run of UA. To do this, it needs the following information:
1) A UA executable - this is in the folder set by "Find UA Folder"
2) A PDB file for input for UA (click the PDB target button and find one).
3) An NP definition. If the auto NP button is checked, an NP will be made for you based on the options set: you can set the radius of the NP, the value of the electrostatic potential at the surface in mV and choose a material from all the ones it found in the MaterialSet.csv . You can load in extra materials with the Load Material Set button if you know what you're doing and have another prepared - "surface-pmfp/MaterialSetPMFP.csv" would be a typical choice here to load in the materials generated via machine learning. More advanced users can uncheck the auto NP button and use the NP file path box and button to find a target NP.
4) Somewhere to store the results - for quick testing the default location is fine, but if you want to change it click the "Result folder" button.

Once all that's done, you can click the "Run UA" button. This will start printing information to the big text box at the bottom of this tab and bring up an alert when it's finished. You can then move to the "View Results" tab.


----------------------
View Results
----------------------
This tab is designed to let you look at the results from a UA run quickly and easily - there are better tools for publication/presentation graphics so use them for final figures. If you've just ran something using the UA Setup tab some of the options here will be pre-set to speed things up, but this can also be used stand-alone to look at other data. Typically, you'll want to do the following steps:
1) Load in a .uam result file ("Load .uam" button) - this may be preset to the result folder from UA Setup to speed things up, but you still need to select the file you want to look at. Remember UA output usually has the format PROTEINNAME_RADIUS_ZETA.uam . Once this is loaded, you should see a heatmap plot appear in the left-hand pane.
2) Load in a .pdb file - this may be preset to the one used on UA Setup but you need to press the button. If all goes well, you'll see the biomolecule and a representation of the NP in the right hand pane.
3) Make sure the NP radius is set to the correct value
4) (Optional) Click "Colour by Boltz. distance". This assigns colours to each biomolecule bead based on their average distance to the surface of the NP - red ones are in close contact more often. This is computed based on a boltzmann-weighted average over all orientations. The "Set to Boltz minimum" button reorients the protein to get as many beads as close to their average distance as possible. 
5) Click "set to min energy" to orient the protein to its most energetically favourable conformation. You can also click directly on the heatmap to set the orientation to the corresponding values of theta, phi, or manually adjust these using the spin boxes phi, theta. The omega dial lets you spin the protein around the z-axis - this doesn't change binding energies due to the assumption of an isotropic NP, but can help you see why a given conformation is favourable.
 6) The simple and Boltzmann weighted energies are displayed automatically, as is the energy for the current orientation.


At the moment the graphics view assumes a spherical NP and shows the xz plane, more options will be available later.

