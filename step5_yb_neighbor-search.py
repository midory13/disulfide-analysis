#!/usr/bin/env python3 

import os
from string import ascii_uppercase
from Bio.PDB import *
import argparse
import os.path
import numpy as np
import sys
import re

__author__ = "Yana Bodnar"
__email__ = "yb143526@uni-greifswald.de"
__status__ = "development"





parser = argparse.ArgumentParser(
    description='This script is there for you to find all residues that are neighbouring the sulfur of the target cystein (or methionine) in specified radius.'
                +'For each entry, this script create a file with information about contacts found in the specified radius.'
                +'Since it creates a single file for each entry, it is advisible to create a separate folder for these files [can be specified under --folder]'
                +'The so created contact files can be further evaluated via yb_neighbor-search_part2.py script'
                +'Important note: the atom types of the sulfur that are supported are: SG, S, SG1, SG2, SD; other types will not be recognised')
parser.add_argument('-i', '--input', required=True, type=str,
                    help="Specify the input database. The file should contain three columns:" 
                    +"1st: name of the protein, 2nd: pdb-id, 3rd: residue number" 
                    +"Please note that this version of the script can only work with single entries, e.g. if your protein has several pdb structures, all of them should be entered as separate entries")
parser.add_argument('-d', '--distance', required=False, type=float,
                    help="Specify the radius for the neighbours searching in Angstrom, [default = 3.5A]")
parser.add_argument('-e', '--error', required=False, type=str,
                    help="Specify the name  of the error file, which would contain a list of entries unprocessed by the script [default: Error.dat]")
parser.add_argument('-f', '--folder', required=False, type=str,
                    help="path for saving out the output data. Note, that the folder has to exist prior to start of the calculation [default folder name: contacts]")
parser.add_argument('-p', '--pdbpath', required=False, type=str,
                    help="Specify path to the folder containing pdb structures. If none specified, the script will search for PDB-Database/PDB folder, and create one if ot does not exists"
                    +"If the PDB structure could not be found, the script will automatically download it from PDB")
args = parser.parse_args()

database_file = args.input #'Cys_database_SS_bigger.dat'

if args.folder:
    folder = args.folder
else:
    folder = "contacts"

if args.error:
    error_file = args.error
else:
    error_file = "Error.dat"

if args.distance:
    radius = args.distance
else:
    radius = float(3.5)

if args.pdbpath:
    pdb_path = args.pdbpath
else:
    pdb_path = "PDB-database/PDB/"

pdb_database_path = "PDB-database/"

pdbl = PDBList()
p = PDBParser()
with open(database_file, "r") as file:
    for line in file:
        if not line.startswith("#"):
            line = line.strip()
            line = line.split("\t")
            name = line[0]
            if re.search(r"\W", name):
                name = re.sub("\W", "_", name)
            pdb = line[1]
            residue = line[2]
            print(residue)
            residue_original = line[2]
            if len(pdb) > 4:
                chainID = pdb[4]
                pdbID = str(pdb[0]) + str(pdb[1]) + str(pdb[2]) +str(pdb[3])
                chain_manual = False
            else: 
                pdbID = pdb
                chainID = 'A'
                chain_manual = True
            if not os.path.exists(folder):
                os.makedirs(folder)
            results_file = str(folder)+"/"+str(name)+ "_" + str(pdb) + "_"+ str(residue)+ "_"+ str(radius) + "_contacts.dat"
            pdb_file = str(pdb_path)+ pdbID + '.pdb'
            if not os.path.isfile(pdb_file):
                print("Downloading the pdb file for " + str(name))
                pdbl.retrieve_pdb_file(pdbID, pdir = pdb_database_path, file_format = "pdb") 
                print("Download done!")
                pdbID_file = str(pdb_database_path) + "pdb" + str(pdbID.lower()) + ".ent"
                structure = p.get_structure ("X", pdbID_file)
            else:
                structure = p.get_structure ("X", pdb_file)
            structure_1 = structure[0]
            if chain_manual == True:
                for letter in ascii_uppercase:
                    #print(letter)
                    if structure_1.has_id(letter):
                        chainID = letter
                        break
            center_chain = structure[0][chainID]
            if center_chain.has_id(int(residue)):
                print("X")
                center_residue = structure[0][chainID][int(residue)]
                if center_residue.has_id("SG"):
                    target_atom = structure[0][chainID][int(residue)]["SG"]
                elif center_residue.has_id("S"):
                    target_atom = structure[0][chainID][int(residue)]["S"]
                elif center_residue.has_id("SG1"):
                    target_atom = structure[0][chainID][int(residue)]["SG1"]
                elif center_residue.has_id("SG2"):
                    target_atom = structure[0][chainID][int(residue)]["SG2"]
                elif center_residue.has_id("SD"):
                    target_atom = structure[0][chainID][int(residue)]["SD"]
            else:
                print("Could not find residue " + str(residue) + " in prtein " + str(pdb) + "\t"+str(name))
                continue
            #try:
            center_atom = Selection.unfold_entities(center_residue, "A")
            ns = NeighborSearch(center_atom)
            residues = Selection.unfold_entities(structure, "R") 
            atoms = Selection.unfold_entities(structure, "A")#returns a list of all atoms in structure, A stays for Atoms, can be either R residues C chains
            ns = NeighborSearch(atoms)
            close_atoms = ns.search(target_atom.coord, radius)
            close_residues = Selection.unfold_entities(close_atoms, "R")
            for ca in close_atoms:
                cr = Selection.unfold_entities(ca, "R")
                distance = target_atom - ca
                #print(ca, cr, distance)
                with open(results_file, "a") as file:
                    #file.write(str(name)+ "\t" + str(pdb) + "\t"+ str(residue) + "\t" + str(cr) + "\t" +str(ca) + "\t" +str(distance)+ "\n")
                    file.write(str(name)+ "\t;" + str(pdb) + "\t;"+ str(residue) + "\t;" + str(cr) + "\t;" +str(ca) + "\t;" +str(distance)+ "\n")
            #except:
            #    print("Could not find neighbors")
            #    with open(error_file, "a") as file:
            #        file.write(str(name)+ "\t;" + str(pdb) + "\t;"+ str(residue) + "\n")