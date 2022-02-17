#!/usr/bin/env python3 


from __future__ import print_function
import os
import re
import argparse
#import matplotlib.pyplot as plt
from Bio.PDB import *
import pandas as pd 
import numpy as np 

__author__ = "Yana Bodnar"
__email__ = "yana.bodnar@stud.uni-greifswald.de"
__status__ = "development"


parser = argparse.ArgumentParser(
    description='Script for analysis of a database excel file, which calculates the distance between two cysteines in reduced and oxidised form.\n' 
                +'The script downloads the desired structures in two formats: *.ent and *.pdb.\n'
                + 'The script requires an input file in excel format xlsx, which would contain following columns: \n'  
                +'\"PDBs with Bond\", \"PDB missing Bond\", \"PDB Cys1\", \"PDB Cys2\", \"Compound\", \"Organism\".\n'  
                +'The cells containing PDB-IDs of the analysed proteins can contain several PDB-IDs, if they are separated by comma.\n'
                +'This script uses following non-standard python modules: Bio.PDB, pandas, numpy. Make sure those are installed on your system!')
parser.add_argument('-f', '--file', required=True, type=str,
                    help="Specify the input excel file which contains the database information.")
parser.add_argument('-s', '--sheet', required=False, type=str,
                    help="Specify the name of the sheet in the input file which contains the database [default = Sheet1] ")
parser.add_argument('-o', '--output', required=False, type=str,
                    help = "Specify the name of the output file. [default = results.txt]")
parser.add_argument('--atom', required=False, type=str, help = "Specify the name of the atoms you want the distance to be measured between. Important: the atom name should be as specified in the PDB file. [default = SG ]")
args = parser.parse_args()

def download_calculate(pdbID, CysRes1, CysRes2, atom_name): # fucntion for downloading the given PDB structure (PDB-ID) and calculating the Distance between the given atoms
    if len(pdbID) > 4:
        chainID = pdbID[4]
        pdbID = str(pdbID[0]) + str(pdbID[1]) + str(pdbID[2]) + str(pdbID[3])
    try:
        pdbl.retrieve_pdb_file(pdbID, pdir = "PDB-database", file_format = "pdb")
        p = PDBParser()
        pdbID_file = "PDB-database/pdb" + str(pdbID.lower()) + ".ent"
        structure = p.get_structure ("X", pdbID_file)
        io = PDBIO()
        io.set_structure(structure)
        if not os.path.exists("PDB-database/PDB/"):
            os.makedirs("PDB-database/PDB/")
            print("created a folder for PDBs")
        io.save("PDB-database/PDB/"+ pdbID + ".pdb")
        r = structure.get_residues()
        try:
            atom1 = structure[0][str(chainID)][int(CysRes1)][atom_name]
            atom2 = structure[0][str(chainID)][int(CysRes2)][atom_name]
            X1, Y1, Z1 = atom1.get_coord()
            X2, Y2, Z2 = atom2.get_coord()
            SSdistance = ((X1 - X2)**2 + (Y1 - Y2)**2 + (Z1 - Z2)**2)**0.5
        except:
            print("Could not calculate the distance!")
            SSdistance = "Could not download the pdb-file: " + str(pdbID)
    except:
        print("Could not download the pdb-file: " + str(pdbID))
        SSdistance = "Could not download the pdb-file: " + str(pdbID)
    return SSdistance


##READING THE USER PARSER INPUT OF SETTING TO DEFAULT

filename = args.file

if args.sheet:
    sheet = args.sheet
else:
    sheet = "Sheet1"
if args.atom:
    atomname = args.atom
else:
    atomname = "SG"
if args.output:
    Results_file = args.output
else:
    Results_file = "results.txt"
##___________________________________________


pdbl = PDBList() #loads Bio.Python

###READING THE DATABASE FILE
database = pd.read_excel(io=filename, sheet_name=sheet, header=0, index_col=10, engine='openpyxl')
df = pd.DataFrame(database, columns = ["PDBs with Bond", "PDB missing Bond", "PDB Cys1", "PDB Cys2", "Compound", "Organism"])
##___________________




with open(Results_file, "w") as file:
    file.write("Protein"+ "\t" + "pdbID-OX"+ "\t" + "pdbID_RED"+ "\t" + "CysRes1"+ "\t" + "CysRes1"+ "\t" + "SSdistance-OX"+ "\t" + "SSdistance-RED" + "\n")
for index, row in df.iterrows():
    PDBListOX = []
    PDBListRED = []
    ####FILLINF THE PDB LISTS with the PDB-ids
    print("check it out: " + str(row["PDBs with Bond"]))
    if row["PDBs with Bond"]!="nan":
        if re.search(r"[\s\[\]\']", row["PDBs with Bond"]): # if the column contains spaces or []-brackets -> remove them
            #print(row["PDBs with Bond"])
            row["PDBs with Bond"] = re.sub(r"[\s\[\]\']", "", row["PDBs with Bond"])
        #print(row["PDBs with Bond"])
            if re.search(r",", row["PDBs with Bond"]): # if the column contains commas, split the Cell by the commas
                PDBListOX = row["PDBs with Bond"].split(",")
            else:
                PDBListOX.append(row["PDBs with Bond"])
    if row["PDB missing Bond"]!="nan":
        if re.search(r"[\s\[\]\']", row["PDB missing Bond"]): # if the column contains spaces or []-brackets -> remove them
            row["PDB missing Bond"] = re.sub(r"[\s\[\]\']", "", row["PDB missing Bond"])

        if re.search(r",", row["PDB missing Bond"]):
            PDBListRED = row["PDB missing Bond"].split(",")
        else:
            PDBListRED.append(row["PDB missing Bond"])
    ####Download the pdbs based on the pdb-ids from the lists, calculate the distance and save out the results in outputfile
    CysRes1 = row["PDB Cys1"]
    CysRes2 = row["PDB Cys2"]
    ProteinName = row["Compound"]
    OrganismName = row["Organism"]
    #print(PDBListOX)
    for protein in PDBListOX:
        print("Working on protein: " + str(ProteinName)+ " " + str(protein))
        SSdistanceOX = download_calculate(protein, CysRes1, CysRes2, atomname)
        with open(Results_file, "a") as file:
            file.write(ProteinName + "\t"+ str(protein) + "\t" 
                                + "\t" + str(CysRes1)+ "\t" + str(CysRes2) 
                                + "\t" + str(SSdistanceOX) + "\t"+"\t"  +"\n")
    for protein in PDBListRED:
        print("Working on protein: " + str(ProteinName) + " " + str(protein))
        SSdistanceRED = download_calculate(protein, CysRes1, CysRes2, atomname)
        with open(Results_file, "a") as file:
            file.write(ProteinName + "\t" + "\t" 
                                + str(protein) + "\t" + str(CysRes1)+ "\t" + str(CysRes2) 
                                + "\t"  + "\t"+ str(SSdistanceRED)  +"\n")

