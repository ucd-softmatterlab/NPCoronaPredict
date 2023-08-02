import numpy as np


def getColor(colorID, npColorID):
    if colorID < npColorID+1:
        return colorID - 1
    if colorID >= npColorID+1:
        return colorID 

coronaFilePath = "CoronaPredictionProjects/AuFCC100UCD-pmfp-daphnia-1ms/coronakmc_output/coronakmcinput.csv_5.0.kmc"
coronaData = open(coronaFilePath,"r")
seenProteins = {"NPDummy":0}
npRadius = 5
npColorID = 3

print("mol new")
print("graphics 0 color "+str(getColor(proteinIDNum,npColorID)))
print("graphics 0 material BrushedMetal")
print("graphics 0 sphere {"+ str(round(0))+" "+str(round(0))+" "+str(round(0))+"} radius "+str(round(npRadius*10))+" resolution 80")
for line in coronaData:
    if line[0] == "#":
        continue
    lineTerms = line.split(",")
    proteinIDSet = lineTerms[0].split(":")
    if proteinIDSet[0] in seenProteins:
        proteinIDNum = seenProteins[proteinIDSet[0]]
    else:
        proteinIDNum = len( seenProteins)
        seenProteins[proteinIDSet[0] ] = proteinIDNum
        print("mol new")
        print("graphics "+str(proteinIDNum)+" color "+str(getColor(proteinIDNum,npColorID)))
        print("graphics "+str(proteinIDNum)+" material Diffuse")
    #print(lineTerms[0], proteinIDNum )
    beadRadius = float(lineTerms[2])
    #print("bead radius: ", beadRadius  )
    beadPhi = float(lineTerms[7])
    beadTheta = float(lineTerms[8])
    beadX = np.cos(beadPhi)*np.sin(beadTheta)*(npRadius+beadRadius)
    beadY = np.sin(beadPhi)*np.sin(beadTheta)*(npRadius+beadRadius)
    beadZ = np.cos(beadTheta)*(npRadius+beadRadius)
    print("graphics "+str(proteinIDNum)+" sphere {"+ str(round(beadX*10,3))+" "+str(round(beadY*10,3))+" "+str(round(beadZ*10,3))+"} radius "+str(round(beadRadius*10))+" resolution 80")



