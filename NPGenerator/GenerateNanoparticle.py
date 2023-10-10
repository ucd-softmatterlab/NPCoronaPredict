#A script for generating a set of nanoparticles for use with UnitedAtom
#This tool assumes a core-shell-brush structure is required. 
#For aggregates or raspberry models with no well-defined structure use something else.
import numpy as np
import BrushGenerator as brushgen
import RaspberryGenerator as raspgen
import os


def getArvoMatrix():
    piVal = np.pi;
    u1= np.random.random();
    u2= np.random.random();
    u3= np.random.random();
    #print(u1,u2,u3)
    
    rxx = -(np.cos(2*piVal*u1)*(1 - 2*u3*np.power(np.cos(2*piVal*u2),2))) - 2*u3*np.cos(2*piVal*u2)*np.sin(2*piVal*u1)*np.sin(2*piVal*u2)
    rxy = -((1 - 2*u3*np.power(np.cos(2*piVal*u2),2))*np.sin(2*piVal*u1)) +    2*u3*np.cos(2*piVal*u1)*np.cos(2*piVal*u2)*np.sin(2*piVal*u2)
    rxz = 2*np.sqrt(1 - u3)*np.sqrt(u3)*np.cos(2*piVal*u2)
    ryx = 2*u3*np.cos(2*piVal*u1)*np.cos(2*piVal*u2)*np.sin(2*piVal*u2) + np.sin(2*piVal*u1)*(1 - 2*u3*np.power(np.sin(2*piVal*u2),2))
    ryy = 2*u3*np.cos(2*piVal*u2)*np.sin(2*piVal*u1)*np.sin(2*piVal*u2) - np.cos(2*piVal*u1)*(1 - 2*u3*np.power(np.sin(2*piVal*u2),2))
    ryz = 2*np.sqrt(1 - u3)*np.sqrt(u3)*np.sin(2*piVal*u2)
    rzx = 2*np.sqrt(1 - u3)*np.sqrt(u3)*np.cos(2*piVal*u1)*np.cos(2*piVal*u2) -    2*np.sqrt(1 - u3)*np.sqrt(u3)*np.sin(2*piVal*u1)*np.sin(2*piVal*u2)
    rzy = 2*np.sqrt(1 - u3)*np.sqrt(u3)*np.cos(2*piVal*u2)*np.sin(2*piVal*u1) +    2*np.sqrt(1 - u3)*np.sqrt(u3)*np.cos(2*piVal*u1)*np.sin(2*piVal*u2)
    rzz = -1 + 2*(1 - u3)
    rotateMatrix = np.array( [[rxx,rxy,rxz],[ryx,ryy,ryz],[rzx,rzy,rzz]])
    #print(rotateMatrix)
    return rotateMatrix

def beadMatch( bead1, bead2):
    for i in range( len(bead1)):
        if bead1[i] != bead2[i]:
            return False
    return True

def saveNP(outputFile,npComponents, doRotate = False, legacyMode = False):
    f = open(outputFile,"w")
    f.write("#"+npName+"\n")
    beadTypeSet = []
    beadLocSet   = []
    innerBoundAuto = -1
    outerBoundAuto = -1
    if doRotate == True:
        rotateMatrix = getArvoMatrix()
    else:
        rotateMatrix = np.array([[1,0,0],[0,1,0],[0,0,1]] )
    for npComponent in npComponents:
        xp = rotateMatrix[0,0]*npComponent[0] + rotateMatrix[0,1]*npComponent[1]  + rotateMatrix[0,2]*npComponent[2] 
        yp = rotateMatrix[1,0]*npComponent[0] + rotateMatrix[1,1]*npComponent[1]  + rotateMatrix[1,2]*npComponent[2] 
        zp = rotateMatrix[2,0]*npComponent[0] + rotateMatrix[2,1]*npComponent[1]  + rotateMatrix[2,2]*npComponent[2] 
        coordString = ",".join( str(w) for w in [xp,yp,zp ] )
        if legacyMode == True:
            f.write( coordString+","  + (",".join(str(w) for w in npComponent[3:]   ))+"\n"  )
        else:
            matchIndex = - 1
            for beadTypeI in range(len(beadTypeSet)):
                if beadMatch( npComponent[3:] , beadTypeSet[beadTypeI] ):
                    matchIndex = beadTypeI + 1
            if matchIndex == -1:
                beadTypeSet.append( npComponent[3:] )
                matchIndex = len(beadTypeSet) 
            #update autobounds
            innerBoundAuto = max( innerBoundAuto, npComponent[3] ) 
            outerBoundAuto = max( outerBoundAuto, npComponent[3] + np.sqrt( npComponent[0]**2 + npComponent[1]**2 + npComponent[2]**2) )
            beadLocSet.append( [ matchIndex-1, xp, yp, zp ])
    if legacyMode == False:
        f.write("INNERBOUND,"+str(innerBoundAuto)+"\n")
        f.write("OUTERBOUND,"+str(outerBoundAuto)+"\n")
        f.write("#Bead definitions \n")
        f.write("#Radius,zeta,corefactor,surfacefactor,shape,hamakerFile,surfaceDir,pmfCutoff,correctionType\n")
        for beadType in beadTypeSet:
            f.write("TYPE," + ",".join( str(w) for w in beadType) +"\n")
        f.write("#Beads\n")
        f.write("#type,x,y,z\n")
        for beadLoc in beadLocSet:
            f.write("BEAD," + ",".join( str(w) for w in beadLoc) +"\n")
    f.close()

targetDir = "../nps/NPGeneratorOutput-test"
os.makedirs(targetDir, exist_ok=True)

#parameters used to control the overall structure

outerShellSurfaceOnly = True #If this parameter is set only the outermost shell (and not the core or inner layers) will have the surface PMF term included.
defaultCutoff = 0.001
defaultSurfaceFactor = 0
finalCutoff = 1 #vdW cutoff of PMF for the outermost layer
#the above parameters are deprecated, these are now set per-material during definitions below

numRandomRotations = 5 #if greater than 0, random orientations of the initial configuration will be generated
saveLegacyNP = False #if true, NPs are saved in the old legacy format.


if outerShellSurfaceOnly == False:
    defaultCutoff = 1
    defaultSurfaceFactor  = 1

#switch on raspberry-style NP generation. core and shell terms are then constructed from sets of beads rather than one large bead.
coreRaspberry = False
raspberryCoreBeadRadius = 0.5
shellRaspberry = False
raspberryShellBeadRadius = 0.4

#A list of the radii for the core NP
coreRadiusRange = [5, 10, 50,100]
#define the hamaker file, surface PMF directory, and the name of the material, LJ cutoff, PMF correction type
#the name is used only for reference and generating file names.
coreMaterials = [
["hamaker/TiO2_Anatase.dat","surface/TiO2-ana-101","Anatase101",1.0,1],
["hamaker/TiO2_Rutile.dat","surface/TiO2-rut-110","Rutile110",1.0,1]
]

#coreHamaker = "hamaker/TiO2_Anatase.dat"
#corePMF = "surface/TiO2-ana-101"


#By default only the outermost shell (and potentially brush elements) have a zeta potential.
#In the Debye-Huckel model the potential only needs to be specified at one radius (and at infinity) so we just set a value for the outer shell
outermostZPRange = [ -0.01, 0, 0.01]

#Define the layers to apply to each NP. Each set can contain multiple layers as shown in the examples below. 
#here we  require the radius, hamaker file, surface directory and a short name for the material, the LJ cutoff used, and the surface correction type (generally 1 for plane PMFs, 0 for beads)
#as an example we define one coating which is a layer of gold and a thin layer of carbon black

shellSetType1 = [
[0.25,"hamaker/Au.dat", "surface/Au/FCC/100/sca","Au",1.0,1],
[0.1,"hamaker/CarbonBlack.dat","surface/CarbonBlack","CB",1.0,1]
]

#this produces a thicker layer of gold.
shellSetType2  = [
[0.5,"hamaker/Au.dat", "surface/Au/FCC/100/sca","Au",1.0,1]
]

#define a null-coating . this causes bare NPs to be generated and is required to make sure the shellVariants has the correct dimensionality.
uncoatedType = []

#set all the shells to loop over. Note that this doesn't combine layers from different sets, it produces e.g. an NP with shell type 1 and another NP with shell type 2
#if you don't want any extra shells then leave only the uncoatedType option
#shellVariants = [uncoatedType]
shellVariants = [shellSetType1, shellSetType2, uncoatedType]

#Definitions for the brush. set applyBrush = False to disable brush entirely, set True to enable the generation of brushes. 
applyBrush = True
#radius, Hamaker, surface, short name, vdw cutoff, zeta
brushBeadDefinition = [0.5,"hamaker/CarbonBlack.dat","surface/CarbonBlack","CB",1.0,0]

#define the brush density as a function of the distance from the surface of the NP as a list of pairs, [d,rho(d)]
#the density here is normalised by the surface area at d- i.e. you count the number of beads at a distance d, then divide by tne surface area 4 Pi (d+RNP)**2
#the un-normalised density is required to be non-increasing as a function of d for the target geometry so take care when generating multiple NP radii
#in general, up-scaling is safe (i.e. inputting the density calculated for a small NP and using this to generate larger NPs) because a given value of d will correspond to a larger area
#conversely, it is not guaranteed that this procedure works for normalised densities generated for a planar NP and down-scaled to smaller spherical NPs so be careful
brushDensityTable =np.array( [ [1, 0.995],[2, 0.03537],[3,0.01989],[4,0.0095493],[5,0.00442],[6,0.002446],[7,0.00124],[8,0] ])
#you can also load this in from a suitably formatted file e.g.
#brushDensityTable = np.genfromtxt("brush_density.csv")







#finally we loop over all shell variants, outermost ZPs and core radii to produce the full set of NPs.
#this will produce quite a few as it makes all combinations - remove those you don't want to keep.
#e.g. for the above example of three coating versions, 3 zeta potentials and four core radii we get 36 NPs out - 12 bare anatase, 12 coated with gold and 12 coated with gold & carbon black

#x,y,z,radius,zeta,coreFactor,surfaceFactor,shape,hamaker,PMF-directory,pmf-vdw-cutoff, correction override
for coreMaterial in coreMaterials:
    for shellSet in shellVariants:
        for outerZP in outermostZPRange:
            for coreRadius in coreRadiusRange:
                npRunningRadius = 0
                npLayers = []
                #print "starting NP"
                if coreRaspberry == True:
                    raspberryBeads = raspgen.buildRaspberry( 0, coreRadius, raspberryCoreBeadRadius   )
                    beadSet = []
                    for bead in raspberryBeads:
                        beadSet.append([bead[0],bead[1],bead[2],raspberryCoreBeadRadius,0,1,defaultSurfaceFactor,1,coreMaterial[0],coreMaterial[1],coreMaterial[3] , coreMaterial[4]])
                    npName = coreMaterial[2]+"-rasp-"+str(coreRadius)
                    if len(beadSet)>0:
                        npLayers.append(beadSet)
                    else:
                        print("Warning: zero beads generated for core raspberry")
                else:
                    npLayers.append([[0,0,0,coreRadius,0,1,defaultSurfaceFactor,1,coreMaterial[0],coreMaterial[1],coreMaterial[3] , coreMaterial[4]]])
                    npName = coreMaterial[2]+"-"+str(coreRadius)
                npRunningRadius += coreRadius
                #npName = coreMaterial[2]+"-"+str(coreRadius)
                for shell in shellSet:
                   if shellRaspberry == True:
                       beadSet = []
                       raspberryBeads = raspgen.buildRaspberry( npRunningRadius, npRunningRadius+shell[0], raspberryShellBeadRadius)
                       for bead in raspberryBeads:
                           beadSet.append([bead[0],bead[1],bead[2], raspberryShellBeadRadius,0,1,defaultSurfaceFactor,1,shell[1],shell[2],defaultCutoff])
                       if len(beadSet)>0:
                           npLayers.append(beadSet)
                       else:
                           print("Warning: zero beads generated for shell raspberry")
                       npName = npName+"_"+shell[3]+"-rasp-"+str(shell[0])
                   else:
                       #first add the anti-NP to subtract the unneeded potential, then the actual NP. together this produces a shell NP.
                       npLayers.append([[0,0,0,npRunningRadius,0,-1, -defaultSurfaceFactor, 1,  shell[1], shell[2], shell[4], shell[5]] ])
                       npLayers.append([[0,0,0,npRunningRadius +shell[0],0,1,defaultSurfaceFactor,1, shell[1], shell[2], shell[4], shell[5]  ]])
                       #npRunningRadius += shell[0]
                       npName = npName+"_"+shell[3]+"-"+str(shell[0])
                   npRunningRadius += shell[0]
                #Next we switch back on the surface PMFs for the outermost material. If shells are present in a non-raspberry model we have to set the anti-NP layer to have a negative surface term to correctly subtract this
                if outerShellSurfaceOnly==True:
                    for index,value in enumerate(npLayers[-1]):
                        npLayers[-1][index][6] = 1
                        npLayers[-1][index][10] = finalCutoff 
                    if len(npLayers) > 1 and shellRaspberry == False:
                        for index,value in enumerate(npLayers[-2]):
                            npLayers[-2][index][6] = -1
                            npLayers[-2][index][10] = finalCutoff
                #next we set the ZP of the outermost layer to be equal to the target value
                for index,value in enumerate(npLayers[-1]):
                    npLayers[-1][index][4] = outerZP
                npName = npName+"_zp"+str(outerZP)
                #generate the brush if needed
                #print npName
                if applyBrush == True:
                    beadRadius = brushBeadDefinition[0]
                    beadZP = brushBeadDefinition[5]
                    brushBeadSet = brushgen.generateBeadSet(beadRadius,npRunningRadius,brushDensityTable)
                    materialBeads = []
                    for bead in brushBeadSet:
                        materialBeads.append([ bead[0], bead[1], bead[2], beadRadius, beadZP, 1,1,1,brushBeadDefinition[1],brushBeadDefinition[2],brushBeadDefinition[4],brushBeadDefinition[5]])
                    npLayers.append(materialBeads)
                    npName = npName+"_brush-"+brushBeadDefinition[3]
                    
                allComponents = []
                for npLayer in npLayers:
                    for npComponent in npLayer:
                        allComponents.append( npComponent)
                print("Generated: ", npName)
                saveNP(targetDir+"/"+npName+".np",allComponents, doRotate=False,legacyMode=saveLegacyNP)
                #saveNP(targetDir+"/"+npName+".np",allComponents, doRotate=False,legacyMode=True)
                for n in range(numRandomRotations):
                    saveNP(targetDir+"/"+npName+"_"+str(n+1)+".np",allComponents, doRotate = True, legacyMode = saveLegacyNP)



