import numpy as np
import scipy.special as scspec

'''
Given a surface potential (approximately the zeta potential) for a NP of a specified size and type in a solution defined by relative permittivty epsilon_1 and Debye parameter kappa1, this script returns the surface potentials for nanoparticles of a range of radii in a solution defined by epsilon_2,kappa2.


Three geometries are considered - sphere, infinitly long cylinder, cube - each defined in terms of a size parameter R. For the sphere and cylinder R gives the radius. For the cube, R is equal to the half-length of a side such that a concentric sphere of radius R touches each face at the centre.

 
The relation between surface charge and R,kappa for the sphere and cylinder can be found analytically by integrating the known potential over all space outside the NP to find the total charge in solution and equating this to the charge on the NP. For a cube, this cannot be performed analytically and instead an approximate power series is obtained by performing this integration numerically over the numerically obtained potential.
'''

def sigmaByKappaCubeFunc(Rkappa):
    return (1 + 0.128311/(Rkappa**2  )+ 0.297472/(Rkappa ))

def sigmaToPsiSphere(R, kappa, eps, sigma):
    return sigma*R/(eps* (1+ kappa*R))

def sigmaToPsiCylinder(R, kappa, eps, sigma):
    return sigma* scspec.kn(0,kappa*R)/scspec.kn(1,kappa*R) / (kappa*eps)



def sigmaToPsiCube(R, kappa, eps, sigma):
    return sigma/(kappa*eps*sigmaByKappaCubeFunc(R*kappa)) 


def sigmaToPsiPlane(R, kappa, eps, sigma):
    return sigma/(kappa*eps)

def psiToSigmaSphere(R, kappa, eps, psi):
    return kappa*eps*psi*(1+kappa*R)/(kappa*R)  

def psiToSigmaCylinder(R, kappa, eps, psi):
    return kappa*eps*psi * scspec.kn(1,kappa*R)/scspec.kn(0,kappa*R)

def psiToSigmaCube(R, kappa, eps, psi):
    return kappa*eps*psi *sigmaByKappaCubeFunc(R*kappa)
        

def psiToSigmaPlane(R, kappa, eps, psi):
    return kappa*eps*psi 

#These define the measurement conditions for the known zeta potential
psiIn = -0.02 #Measured zeta potential - output values will be in the same units

eps1 = 80 #Relative permittivity of solution of the measured NP
kappa1 = 1.0#Inverse debye length of solution of the measured NP

radiusIn = 5#Radius of the measured NP. Must have same units as 1/kappa1
inputType = "sphere" #Morphology of measured NP

#These define the parameters for the target solution
kappa2 = 1.5
eps2 = 70
targetRadii = [1,2,3,4,5,10,15,20,25,50,100] #size parameters to calculate zeta potentials for

#Calculate the surface charge density for the measured NP
if inputType == "sphere":
    sigmaEst = psiToSigmaSphere(radiusIn, kappa1, eps1, psiIn)
elif inputType == "cylinder":
    sigmaEst = psiToSigmaCylinder(radiusIn, kappa1, eps1, psiIn)
elif inputType == "cube":
    sigmaEst = psiToSigmaCube(radiusIn, kappa1, eps1, psiIn)
else:
    print "Unrecognised input type (sphere/cylinder/cube), defaulting to sphere"
    sigmaEst = psiToSigmaSphere(radiusIn, kappa1, eps1, psiIn)


#The surface charge density is a property of the material, and so can be used to calculate zeta potentials for a range of NPs of different shapes and sizes in a new solution


print "Sphere zeta potentials:"
for radius in targetRadii:
    print radius, sigmaToPsiSphere(radius, kappa2, eps2, sigmaEst)

print "Cylinder zeta potentials:"
for radius in targetRadii:
    print radius, sigmaToPsiCylinder(radius, kappa2, eps2, sigmaEst)


print "Cube zeta potentials:"
for radius in targetRadii:
    print radius, sigmaToPsiCube(radius, kappa2, eps2, sigmaEst)
