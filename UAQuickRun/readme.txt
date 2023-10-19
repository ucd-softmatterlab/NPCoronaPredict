###UAQuickRun README file, last updated 19/10/2023 by Ian Rouse, ian.rouse@ucd.ie

UAQuickRun is a simple graphical interface for running UnitedAtom for a single protein - NP pair and visualising the results. Its functionality is very limited by design compared to the full command-line program to make it easy to use for beginners without overwhelming them with too many options. 

Dependencies:
UAQuickRun needs a working qt5 installation. The quick and easy way to get this on Linux is:
sudo apt install qtchooser qtbase5-dev qt5-qmake 

Compiling UAQuickRun:
The source code is available as part of the NPCoronaPredict package and requires gcc and Qt5 for installation. To compile the program, go into the UAQuickRun folder in a terminal window and run the commands:
qmake -qt=qt5
make

If all goes well there will be a few compiler warnings about unused variables, and then a UAQuickRun executable will be generated in the UAQuickRun folder. You should be able to run this like any other program. Depending on what other packages you have installed, you might need to do further work to get a working QT runtime - check the usual sources (google, stackoverflow, etc). 


