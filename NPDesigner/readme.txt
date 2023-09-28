###NPDesigner README file, last updated 28/09/2023 by Ian Rouse, ian.rouse@ucd.ie

NPDesigner is a tool to make it easy to build spherical NPs for use in the NPCoronaPredict set of tools, particularly UnitedAtom. 

Dependencies:
NPDesigner needs a working qt5 installation. The quick and easy way to get this on Linux is:
sudo apt install qtchooser qtbase5-dev qt5-qmake 

Compiling NPDesigner:
The source code is available as part of the NPCoronaPredict package and requires gcc and Qt5 for installation. To compile the program, go into the NPDesigner folder in a terminal window and run the commands:
qmake -qt=qt5
make

If all goes well there will be a few compiler warnings about unused variables, and then an NPDesigner executable will be generated in the NPDesigner folder. You should be able to run this like any other program.

Using NPDesigner:
A typical workflow would look like this:

1) Set the UnitedAtom base directory - this is usually one level up from NPDesigner and should be the folder containing the UA executable. This is necessary so that path names are generated properly when you're saving output.
2) Add the first NP bead type using the "Add Bead Type" button. In this dialogue you can set the Hamaker file and PMF directory, a surface electrostatic potential and a radius for the bead. The other parameters might need adjusting depending on your exact needs. If you've made any mistakes you can double-click the entry in the table and make changes, then press "Update from tables" to save changes.
3) Add a bead of this type using "Add NP Bead". This lets you specify the co-ordinates and bead-type for a new bead. Again, if you make mistakes you can edit them directly in the table. Be careful editing the bead-type for specific NP beads.
4) Add a shell - this is a continuum-level model for a solid layer of a different material. This dialogue lets you specify the surface/hamaker files, the inner and outer radius of the shell, and the zeta potential at once, and produces a pair of bead types and the corresponding NP beads all in one go to make a "large NP" and a "negative small NP". They are centered at coordinates 0,0,0 but can be moved later as needed.
5) Add some brush - these are single-beads attached to the outside of the NP. These must be a pre-defined bead type, so if you want them to be a different material or size add this first. The brush dialogue lets you choose the bead type to place, the radius at which the bead centers should be placed (with a suggested value underneath to attach them to the previous layer), and set the occupancy of this layer, i.e., how many potential sites should be filled with beads. If the "allow floating beads" option is checked then beads can be free-standing and will be generated at any valid location up to the allowed occupancy, else they will only be generated if they are in direct contact with another bead (which in turn must also be in contact with another and so on). This algorithm iterates until either no more beads can be placed or the target number has been met. The bottom two fields show the expecte density and number of beads for this layer - if you have a target value for either of these then alter occupancy until its at the correct value.
6) Keep doing this until you've built your NP. The window in the upper-left will show a very rough image of what it looks like. Remember that if you edit the tables manually you need to press "Update from tables" to make the changes take place.
7) When you're happy with the NP it can be saved for output. If you want to set 0,0,0 to the center of mass (all beads are assigned a mass based on their volume) then press "Re-center". If a lot of the mass of the NP is excluded from the blue dashed circle in the output window you should uncheck "Auto-set binding radii" and increase R0 until it covers what can reasonably be thought of as the main part of the NP. Make sure R1 is equal to or greater than this number. The red dashed line shows R1 + 2 nm, since this is the limit used in UnitedAtom - integration of protein energies takes place over the region between the blue dashed line to the red dashed line. 
8) If you want multiple random orientations of the NP, set this appropriately. If its left to 1 only the orientation you've been building is saved out. If its >1 , your manual orientation is saved and X random orientations are generated and saved too. 
9) Press Save NP(s) and save the NP somewhere. This will produce myname.np with co-ordinates matching those in the table. If you've set orientations to a number X > 1, it will also produce myname_1.np, myname_2.np .... myname_X.np , each of which have random orientations. 
Note that NPs generated using this tool use the new NP file convention where bead types and beads are separated. An example output will look like this:


#NPDesigner NP file 
#NP-parameters 
INNERBOUND,1
OUTERBOUND,6
#Bead definitions 
#Radius,zeta,corefactor,surfacefactor,shape,hamakerFile,surfaceDir,pmfCutoff,correctionType 
TYPE,1,0,1,1,1,hamaker/MyHamaker.dat,surface/mymaterial,1,0
#Beads
#type,x,y,z
BEAD,0,-4.95083,0.64961,0.259408
BEAD,0,-3.94877,2.40231,-1.90685
BEAD,0,-4.1395,-0.79659,-2.68886
BEAD,0,-4.51572,-1.78342,-1.19488
BEAD,0,-4.24174,1.59553,2.11232
BEAD,0,-3.67238,3.13976,1.28667
BEAD,0,-2.64058,3.71633,-2.05335
BEAD,0,-3.08075,-2.42032,-3.10661
BEAD,0,-3.69276,-1.12893,3.17632
BEAD,0,-2.85582,2.38099,3.34293
BEAD,0,-2.26697,3.80569,2.31896
BEAD,0,-1.67733,4.64549,0.778475
BEAD,0,-0.892335,-0.662153,-4.87497
BEAD,0,-1.9713,-4.40273,-1.31528
BEAD,0,-1.66725,0.887298,4.62958
BEAD,0,0.56513,4.95875,-0.302327
BEAD,0,0.952797,4.4621,-2.04496
BEAD,0,1.20903,1.79218,-4.50847
BEAD,0,1.07202,-0.0173105,-4.8837
BEAD,0,0.431408,-3.38505,-3.65449
BEAD,0,0.0457123,-4.48346,-2.21281
BEAD,0,0.244542,0.937528,4.90523
BEAD,0,1.40093,3.782,2.95531
BEAD,0,1.9733,4.40429,1.30703
BEAD,0,2.4483,4.33284,-0.481984
BEAD,0,2.90375,2.22723,-3.40701
BEAD,0,2.61705,-1.33415,-4.04612
BEAD,0,2.27352,-2.99261,-3.29779
BEAD,0,1.88716,-4.20018,-1.94862
BEAD,0,1.79077,-0.351002,4.6551
BEAD,0,2.31149,1.36895,4.21699
BEAD,0,2.90783,2.68131,3.05861
BEAD,0,3.47917,3.30124,1.41322
BEAD,0,3.92427,3.07665,-0.366539
BEAD,0,4.16384,2.03658,-1.87476
BEAD,0,3.9396,-1.42719,-2.72812
BEAD,0,3.57874,-2.98635,-1.80952
BEAD,0,3.20716,-3.83002,-0.212286
BEAD,0,2.96898,-3.68707,1.60955
BEAD,0,3.28679,-0.861201,3.66815
BEAD,0,3.82104,0.809249,3.12166
BEAD,0,4.40973,1.70095,1.63127
BEAD,0,4.81971,1.32586,-0.111882
BEAD,0,4.57077,-1.91834,-0.654229

10) If you have a pre-built NP you can use the File->Open interface to load this in. This supports old-style NPs too.
