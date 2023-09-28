import numpy as np
import pandas as pd
import sklearn.cluster as skcluster
import sklearn.decomposition as skdecomp
from scipy.cluster.vq import vq
import warnings
import os
import shutil


warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)
proteinSet = pd.read_csv("DaphniaProteinStats.csv",skipinitialspace=True)
#print(proteinSet.keys().tolist() )



serumProteinConcs = []
coronaSerumProteinFile = "daphnia_serum/daphnia_serum_mass1.0_fixedmassfrac.csv"

if coronaSerumProteinFile!= "":
    serumFile = open(coronaSerumProteinFile,"r")
    serumFile.readline()
    for line in serumFile:
        if line[0]=="#":
            continue
        lineTerms = line.strip().split(",")
        serumProteinConcs.append(float(lineTerms[1]))
serumProteinConcs = np.array(serumProteinConcs)


targetDescriptors = ["Mass (Da)", "Charge", 'AlaNumber',  'CysNumber',  'AspNumber', 'GluNumber',  'PheNumber',  'GlyNumber',  'HisNumber', 'IleNumber', 'LysNumber',  'LeuNumber','MetNumber', 'AsnNumber',  'ProNumber',  'GlnNumber', 'ArgNumber',  'SerNumber', 'ThrNumber',  'ValNumber',  'TrpNumber',  'TyrNumber',  'TinyNumber',  'SmallNumber',  'AliphaticNumber', 'AromaticNumber', 'NonPolarNumber', 'PolarNumber','ChargedNumber',  'BasicNumber','AcidicNumber']

print("Material," + ",".join([ str(a) for a in targetDescriptors])  +"," +   ",".join([ "enrich_"+str(a) for a in targetDescriptors])     +"," +   ",".join([ "massenrich_"+str(a) for a in targetDescriptors])   )

materialSet = ["silicaquartz","silicaamorph","anatase101","rutile110","rutile100","fe2o3","CdSe","carbonblack","grapheneoxide","redgrapheneoxide","graphene","gold","cnt"]
fingerprintOut = open("fingerprint_output.txt","w")
fingerprintOut.write("Material," + ",".join([ str(a) for a in targetDescriptors])  +"," +   ",".join([ "enrich_"+str(a) for a in targetDescriptors])     +"," +   ",".join([ "massenrich_"+str(a) for a in targetDescriptors])+"\n"  )

for targetMaterial in materialSet:
    coronaStatsFile = open("corona_results/daphnia-fingerprint-"+targetMaterial+"/kmc_running_np1R_5_ZP_0.csv_5.0_s0_hs_0.txt","r")
    headerLine = coronaStatsFile.readline()
    headerTerms = headerLine.strip().split(",")
    proteinSetTemp = headerTerms[1:-2]
    proteinSetNames = []
    totalSerumMass = 0
    massWeightedDescriptors = np.zeros( len(targetDescriptors) )
    proteinMassList = []
    proteinDescriptorVals = []
    for protein in proteinSetTemp:
        proteinTerms = protein.split("-P")
        #print(proteinTerms)
        proteinSetNames.append(proteinTerms[0])
        proteinMass = proteinSet.loc[   proteinSet["Name"] == proteinTerms[0], "Mass (Da)"].to_numpy()[0]
        proteinMassList.append(proteinMass)
        proteinDescriptors = proteinSet.loc[   proteinSet["Name"] == proteinTerms[0], targetDescriptors].to_numpy()[0]
        #print(proteinTerms[0], proteinDescriptors)
        proteinDescriptorVals.append(proteinDescriptors)
        massWeightedDescriptors += proteinDescriptors*proteinMass
        totalSerumMass +=  proteinMass
    #print(proteinSetNames)
    allDescriptorVals = np.stack( proteinDescriptorVals)
    #print(allDescriptorVals)
    proteinMassArr = np.array(proteinMassList)
    massDotConc = np.dot(proteinMassArr, serumProteinConcs)
    descriptorsDotConc = np.dot( np.transpose(allDescriptorVals), serumProteinConcs)
    unitVector = np.ones_like(serumProteinConcs)
    totalConc = np.sum(serumProteinConcs)
    massDotUnity = np.dot(proteinMassArr, unitVector)
    descriptorsDotUnity = np.dot( np.transpose(allDescriptorVals), unitVector)
    #
    for line in coronaStatsFile:
        lineTerms = line.strip().split()
        #print(lineTerms[0])
        proteinNums = np.array( [float(a) for a in lineTerms[1:-2]] )
        #print(proteinNums)
        descriptorNums =  np.dot( np.transpose(allDescriptorVals), proteinNums)
        #print( str(lineTerms[0]) +"," + )
        totalAdsorbedMass = np.dot( proteinNums, proteinMassArr)
        totalAdsorbedNumber = np.sum(proteinNums)
        rawString = ",".join([ str(a) for a in descriptorNums]) 
        lastLine=targetMaterial +"," + rawString 
    print(lastLine)
    fingerprintOut.write(lastLine+"\n")
fingerprintOut.close()
    
