#!/usr/bin/env python3 


from __future__ import print_function
import os
import re
import argparse
from Bio.PDB import *
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

__author__ = "Yana Bodnar"
__email__ = "yb143526@uni-greifswald.de"
__status__ = "development"


def average(List):
    if not List:
        x = ""
    else:
        x = sum(List)/len(List)
    return x

parser = argparse.ArgumentParser(
    description='This script is made for the analysing of the ouput of yb_database_analysis.py script.\n' + 
                'It averages the distances between Cysteins in oxidised and reduced structures of one protein (which was required for proteins with several PDB structures)\n'+
                'and saves the results of the calcultions into the outputfiles: "difference.txt" and "biggest_difference.txt".\n' +
                'In the latter only proteins with the distance differences above 2Angstrom are printed out.\n'+ 
                'Furthermore, this script creates histograms (in *svg and *png formats) containing descributions \n'+
                'of occurence of the distances between S-S of the Cysteins and distribution of the distance differences between ox and red proteins.\n')
parser.add_argument('-f', '--file', required=True, type=str,
                    help="Specify the input file for the calculation.")
args = parser.parse_args()


ProteinDict = {}
SSRedList = []
SSRedDict = {}
SSOxList = []
SSOxDict = {}
DSRedDict = {}
DSOxDict = {}

OX = "OX"
RED = "RED"
with open(args.file, "r") as file:
    for line in file:
        if line.startswith("Protein"):
            continue
        else:
            CysResList = []
            CysResDict = {}
            line = line.split("\t")
            Protein = line[0].strip(" ")
            PDB_OX = line[1]
            PDB_RED = line[2]
            CysRes1 = line[3]
            CysRes2 = line[4]
            SS_OX = line[5]
            SS_RED = re.sub(r"\n", "", line[6])
            CysResList.append(CysRes1)
            CysResList.append(CysRes2)
            CysResStr = str(CysRes1) + "," + str(CysRes2)
            OX = "OX"
            RED = "RED"
            if Protein not in ProteinDict:
                RedOxDict = {}
                ProteinDict[Protein] = {}
                #print("1")
                CysResDict[CysResStr] = {}
                RedOxDict[OX] = []
                RedOxDict[RED] = []
                RedOxDict["PDB_OX"] = []
                RedOxDict["PDB_RED"] = []
                CysResDict[CysResStr] = RedOxDict
                ProteinDict[Protein] = CysResDict
                if not SS_OX == "":
                    if not SS_OX.startswith("Could not"):
                        RedOxDict[OX].append(float(SS_OX))
                        SSOxList.append(float(SS_OX))
                if not SS_RED == "":
                    if not SS_RED.startswith("Could not"):
                        RedOxDict[RED].append(float(SS_RED))
                        SSRedList.append(float(SS_RED))
                if not PDB_OX == "" and not SS_OX.startswith("Could not"):
                    RedOxDict["PDB_OX"].append(PDB_OX)
                if not PDB_RED == "" and not SS_RED.startswith("Could not"):
                    RedOxDict["PDB_RED"].append(PDB_RED)
            else:
                if CysResStr not in ProteinDict[Protein].keys():
                    RedOxDict = {}
                    #print("2")
                    RedOxDict[OX] = []
                    RedOxDict[RED] = []
                    RedOxDict["PDB_OX"] = []
                    RedOxDict["PDB_RED"] = []
                    if not SS_OX == "":
                        if not SS_OX.startswith("Could not"):
                            RedOxDict[OX].append(float(SS_OX))
                            SSOxList.append(float(SS_OX))
                    if not SS_RED == "":
                        if not SS_RED.startswith("Could not"):
                            RedOxDict[RED].append(float(SS_RED))
                            SSRedList.append(float(SS_RED))
                    if not PDB_OX == "" and not SS_OX.startswith("Could not"):
                        RedOxDict["PDB_OX"].append(PDB_OX)
                    if not PDB_RED == "" and not SS_RED.startswith("Could not"):
                        RedOxDict["PDB_RED"].append(PDB_RED)
                    ProteinDict[Protein][CysResStr] = RedOxDict
                else:
                    #print("3")
                    if not SS_OX == "":
                        if not SS_OX.startswith("Could not"):
                            ProteinDict[Protein][CysResStr][OX].append(float(SS_OX))
                            SSOxList.append(float(SS_OX))
                    if not SS_RED == "":
                        if not SS_RED.startswith("Could not"):
                            ProteinDict[Protein][CysResStr][RED].append(float(SS_RED))
                            SSRedList.append(float(SS_RED))
                    if not PDB_OX == "" and not SS_OX.startswith("Could not"):
                        ProteinDict[Protein][CysResStr]["PDB_OX"].append(PDB_OX)
                    if not PDB_RED == "" and not SS_RED.startswith("Could not"):
                        ProteinDict[Protein][CysResStr]["PDB_RED"].append(PDB_RED)

difference_file = "difference.txt"
with open(difference_file, "w") as file:
            file.write("#Protein \t CysPair \t Ox [A] \t Red[A] \t difference red-ox [A] \n"  )
with open("biggest_difference.txt", "w") as file:
            file.write("#Protein \t CysPair \t Ox [A] \t Red[A] \t difference red-ox [A] \t PDB OX \t PDB RED\n"  )
for k1, v1 in ProteinDict.items():  #Protein 
    for k2, v2 in v1.items(): #CysRes
        for k3, v3 in v2.items(): #OX/RED
            if not k3.startswith("PDB"):
                ProteinDict[k1][k2][k3] = average(v3)
            else:
                ProteinDict[k1][k2][k3] = v3

RedOxdiffList =  []
CutOff = 2
for k1, v1 in ProteinDict.items():  #Protein
    for k2, v2 in v1.items(): #CysRes
        try:
            RedOxdiff = v2["RED"] - v2["OX"]
        except:
            RedOxdiff = "Could not calculate"
        RedOxdiffList.append(RedOxdiff)
        if type(RedOxdiff) != str:
            if (RedOxdiff >= CutOff):
                with open("biggest_difference.txt", "a") as file:
                    file.write(str(k1) + "\t" + str(k2) + "\t" + str(v2["OX"])
                                + "\t" + str(v2["RED"]) + "\t" + str(RedOxdiff) 
                                +  "\t" + str(v2["PDB_OX"])+  "\t" + str(v2["PDB_RED"])+" \t \n")
            else:
                print(":C")
        with open(difference_file, "a") as file:
            file.write(str(k1) + "\t" + str(k2) + "\t" + str(v2["OX"])+ "\t" + str(v2["RED"]) + "\t" + str(RedOxdiff) +  "\t" + str(v2["PDB_OX"])+  "\t" + str(v2["PDB_RED"])+ " \t \n" )
RedOxDiffListClean = []
for x in RedOxdiffList:
    if type(x) == str:
        continue
    else:
        RedOxDiffListClean.append(x)   
fig, ax = plt.subplots()
plt.figure(figsize=(5,5), dpi= 250) 
plt.xlabel("Distance in Å")
plt.ylabel("Probability")
plt.title("Distribution of Red-Ox differences")
plt.axvline(2 + int(CutOff), color="grey")
plt.hist(RedOxDiffListClean, alpha=0.8, color="dodgerblue",  bins=200, density=True, stacked = True)
plt.tick_params(axis="both", direction="in", length=10, which="major" )
plt.minorticks_on()
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
plt.tick_params(axis="both", direction="in", length=5, which="minor" )
plt.savefig("disulfides_distributionDIFF.svg")
plt.savefig("disulfides_distributionDIFF.png")
plt.close()


plt.figure(figsize=(5,5), dpi= 250)
plt.title("Distribution of S-S-distance")
plt.xlabel("Distance in Å")
plt.ylabel("Probability")
plt.tick_params(axis="both", direction="in", length=10, which="major" )
plt.minorticks_on()
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
plt.tick_params(axis="both", direction="in", length=5, which="minor" )
plt.hist(SSRedList, alpha=0.8, color="dodgerblue", bins=200, density=True, stacked = True)
plt.savefig("disulfides_distributionRED2.svg")
plt.savefig("disulfides_distributionRED2.png")
plt.hist(SSOxList, alpha=0.8, color="orange", bins = 200, density=True, stacked = True)
plt.xlim(1,15)
#sns.displot(SSRedList, color="dodgerblue", label="reduced", kind="kde")
plt.savefig("disulfides_distribution2.svg")
plt.savefig("disulfides_distribution2.png")