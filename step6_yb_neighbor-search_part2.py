#!/usr/bin/env python3 
import os
from Bio.PDB import *

import os.path
import numpy as np
import sys
import re
import argparse



parser = argparse.ArgumentParser(
    description='This script is there for you to analyse the output contact data from yb_neighbor-search.py.'
                +'The script counts the number of each aminoacids occuring in the contact data (1 count = 1 aa)'
                +'Furthermore, the script creates a file for each AA, containing details about all occuring neigbour atoms (Protein, pdb, residue, atom name)')
parser.add_argument('-i', '--input', required=True, type=str,
                    help="Specify the input file for the calculation:"+
                    "1st column: protein name, 2nd column: pdb-id, 3rd column: residue id" + 
                    "if the list mode was switched on, the pdb-ids should be each separated by commas")
parser.add_argument('-m', '--mode', required=False, type=bool,
                    help="Specify the list mode of action, if true[yes, 1] another type of hierarchy of the input file is required" 
                    +"1st column: name of the protein; 2nd column: list of PDB-ids separated with commas" 
                    +"[default: False]")
parser.add_argument('-d', '--distance', required=False, type=float,
                    help="Specify the radius for the neighbours searching in Angstrom, [default = 3.5A]")
parser.add_argument('-f', '--folder', required=False, type=str,
                    help="path for saving out the output data. Note, that the folder has to exist prior to start of the calculation [default folder name: contacts]")
parser.add_argument('-o', '--output', required=False, type=str,
                    help = "Specify the name of the output file. [default = contacts.dat]")
args = parser.parse_args()


Ala_count =0
Arg_count =0
Asn_count =0
Asp_count =0
Cys_count =0
Glu_count =0
Gln_count =0
Gly_count =0
His_count =0
Ile_count =0
Leu_count =0
Lys_count =0
Met_count =0
Phe_count =0
Pro_count =0
Ser_count =0
Thr_count =0
Trp_count =0
Tyr_count =0
Val_count =0

def count_neigbours(neighbours_file):
    global Ala_count, Arg_count,Asn_count,Asp_count,Cys_count,Glu_count,Gln_count,Gly_count,His_count,Ile_count,Leu_count,Lys_count,Met_count,Phe_count,Pro_count,Ser_count,Thr_count,Trp_count,Tyr_count,Val_count
    if os.path.isfile(neighbours_file):
        print(neighbours_file)
        with open(neighbours_file, "r") as f:
            Ala_lst = []
            Cys_lst = []
            Arg_lst = []
            Asn_lst = []
            Asp_lst = []
            Cys_lst = []
            Glu_lst = []
            Gln_lst = []
            Gly_lst = []
            His_lst = []
            Ile_lst = []
            Leu_lst = []
            Lys_lst = []
            Met_lst = []
            Phe_lst = []
            Pro_lst = []
            Ser_lst = []
            Thr_lst = []
            Trp_lst = []
            Tyr_lst = []
            Val_lst = []
            for line in f:
                line = line.strip()
                line = line.split("\t;")
                cr = line[3]
                match_string = "[<Residue CYS het=\w*  resseq=(\d+)" + str(residue)
                if re.search("CYS", line[3]):
                    searching = re.search(r"\[<Residue CYS het=\w*\s* resseq=(\d+)", line[3])
                    cys_res = searching.group(1)
                    with open(Cys_file, "a") as f1:
                            f1.write(str(line) + "\n")
                    if (float(line[5])>0) and (cr not in Cys_lst) and cys_res != residue:
                        Cys_count += 1
                        Cys_lst.append(cr)
                elif re.search('ALA', line[3]):
                    with open(Ala_file, "a") as f1:
                            f1.write(str(line) + "\n")
                    if cr not in Ala_lst:
                        Ala_lst.append(cr)
                        Ala_count += 1
                elif re.search('ARG', line[3]):
                    with open(Arg_file, "a") as f1:
                            f1.write(str(line) + "\n")
                    if cr not in Arg_lst:
                        Arg_lst.append(cr)
                        Arg_count += 1
                elif re.search('ASN', line[3]):
                    with open(Asn_file, "a") as f1:
                            f1.write(str(line) + "\n")
                    if cr not in Asn_lst:
                        Asn_lst.append(cr)
                        Asn_count += 1
                elif re.search('ASP', line[3]):
                    with open(Asp_file, "a") as f1:
                            f1.write(str(line) + "\n")
                    if cr not in Asp_lst:
                        Asp_lst.append(cr)
                        Asp_count += 1
                elif re.search('GLU', line[3]):
                    with open(Glu_file, "a") as f1:
                            f1.write(str(line) + "\n")
                    if cr not in Glu_lst:
                        Glu_lst.append(cr)
                        Glu_count += 1
                elif re.search('GLN', line[3]):
                    with open(Gln_file, "a") as f1:
                            f1.write(str(line) + "\n")
                    if cr not in Gln_lst:
                        Gln_count += 1
                        Gln_lst.append(cr)
                elif re.search('GLY', line[3]):
                    with open(Gly_file, "a") as f1:
                            f1.write(str(line) + "\n")
                    if cr not in Gly_lst:
                        Gly_count += 1
                        Gly_lst.append(cr)
                elif re.search('HIS', line[3]):
                    with open(His_file, "a") as f1:
                            f1.write(str(line) + "\n")
                    if cr not in His_lst:
                        His_count += 1
                        His_lst.append(cr)
                elif re.search('ILE', line[3]):
                    with open(Ile_file, "a") as f1:
                            f1.write(str(line) + "\n")
                    if cr not in Ile_lst:
                        Ile_count += 1
                        Ile_lst.append(cr)
                elif re.search('LEU', line[3]):
                    with open(Leu_file, "a") as f1:
                            f1.write(str(line) + "\n")
                    if cr not in Leu_lst:
                        Leu_count += 1
                        Leu_lst.append(cr)
                elif re.search('LYS', line[3]):
                    with open(Lys_file, "a") as f1:
                            f1.write(str(line) + "\n")
                    if cr not in Lys_lst:
                        Lys_count += 1
                        Lys_lst.append(cr)
                elif re.search('MET', line[3]):
                    with open(Met_file, "a") as f1:
                            f1.write(str(line) + "\n")
                    searching_MET = re.search(r"\[<Residue MET het=\w*\s* resseq=(\d+)", line[3])
                    met_res = searching_MET.group(1)
                    if (float(line[5])>0) and (cr not in Met_lst) and met_res != residue:
                        if cr not in Met_lst:
                            Met_count += 1
                            Met_lst.append(cr)
                elif re.search('PHE', line[3]):
                    with open(Phe_file, "a") as f1:
                            f1.write(str(line) + "\n")
                    if cr not in Phe_lst:
                        Phe_count += 1
                        Phe_lst.append(cr)
                elif re.search('PRO', line[3]):
                    with open(Pro_file, "a") as f1:
                            f1.write(str(line) + "\n")
                    if cr not in Pro_lst:
                        Pro_count += 1
                        Pro_lst.append(cr)
                elif re.search('SER', line[3]):
                    with open(Ser_file, "a") as f1:
                            f1.write(str(line) + "\n")
                    if cr not in Ser_lst:
                        Ser_count += 1
                        Ser_lst.append(cr)
                elif re.search('THR', line[3]):
                    with open(Thr_file, "a") as f1:
                            f1.write(str(line) + "\n")
                    if cr not in Thr_lst:
                        Thr_count += 1
                        Thr_lst.append(cr)
                elif re.search('TRP', line[3]):
                    with open(Trp_file, "a") as f1:
                            f1.write(str(line) + "\n")
                    if cr not in Trp_lst:
                        Trp_count += 1
                        Trp_lst.append(cr)
                elif re.search('TYR', line[3]):
                    with open(Tyr_file, "a") as f1:
                            f1.write(str(line) + "\n")
                    if cr not in Tyr_lst:
                        Tyr_count += 1
                        Tyr_lst.append(cr)
                elif re.search('VAL', line[3]):
                    with open(Val_file, "a") as f1:
                            f1.write(str(line) + "\n")
                    if cr not in Val_lst:
                        Val_count += 1
                        Val_lst.append(cr)
        with open(cys_file, "a") as f2:
            f2.write(str(name)+ "_" + str(pdb) + "_"+ str(residue) +"\n")
            for x in Cys_lst:
                #print(x)
                match = re.search(r"\[<Residue CYS het=  resseq=(\d+)", x)
                if match:
                    f2.write("\t"+ str(match.group(1)) +"\n")
            f2.write("-------------------\n")
    else:
        print("Could not open file " + str(neighbours_file))




cys_file = "The_Cys-List.dat"

if args.input:
    database_file = args.input
else:
    database_file = 'Cys_database.dat'
if args.output:
    contact_file = args.output
else:
    contact_file = "contacts.dat"

if args.folder:
    folder = args.folder
else:
    folder = "contacts"

if args.distance:
    radius = args.distance
else:
    radius = float(3.5)


Ala_file="Ala_collection_" + str(radius) + "_contacts.dat"
Arg_file="Arg_collection_" + str(radius) + "_contacts.dat"
Asn_file="Asn_collection_" + str(radius) + "_contacts.dat"
Asp_file="Asp_collection_" + str(radius) + "_contacts.dat"
Cys_file="Cys_collection_" + str(radius) + "_contacts.dat"
Glu_file="Glu_collection_" + str(radius) + "_contacts.dat"
Gln_file="Gln_collection_" + str(radius) + "_contacts.dat"
Gly_file="Gly_collection_" + str(radius) + "_contacts.dat"
His_file="His_collection_" + str(radius) + "_contacts.dat"
Ile_file="Ile_collection_" + str(radius) + "_contacts.dat"
Leu_file="Leu_collection_" + str(radius) + "_contacts.dat"
Lys_file="Lys_collection_" + str(radius) + "_contacts.dat"
Met_file="Met_collection_" + str(radius) + "_contacts.dat"
Phe_file="Phe_collection_" + str(radius) + "_contacts.dat"
Pro_file="Pro_collection_" + str(radius) + "_contacts.dat"
Ser_file="Ser_collection_" + str(radius) + "_contacts.dat"
Thr_file="Thr_collection_" + str(radius) + "_contacts.dat"
Trp_file="Trp_collection_" + str(radius) + "_contacts.dat"
Tyr_file="Tyr_collection_" + str(radius) + "_contacts.dat"
Val_file="Val_collection_" + str(radius) + "_contacts.dat"



with open(database_file, "r") as file:
    for line in file:
        if not line.startswith("#"):
            line = line.strip()
            line = line.split("\t")
            if args.mode ==True:
                name = line[0]
                if re.search(r"\W", name):
                    name = re.sub("\W", "_", name)
                residue = line[2]
                if line[1] != "-":
                    for i in line[1].split():
                        match = re.search(r"\'*(\w*)\'*(,*)", i)
                        pdb = match.group(1)
                        neighbours_file = str(folder)+"/"+str(name)+ "_" + str(pdb) + "_"+ str(residue)+ "_"+ str(radius) + "_contacts.dat"
                        count_neigbours(neighbours_file)
            else:
                name = line[0]
                if re.search(r"\W", name):
                    name = re.sub("\W", "_", name)
                pdb = line[1]
                residue = line[2]
                residue_original = line[2]
                neighbours_file = str(folder)+"/"+str(name)+ "_" + str(pdb) + "_"+ str(residue)+ "_"+ str(radius) + "_contacts.dat"
                count_neigbours(neighbours_file)


with open(contact_file, "w") as file:
    file.write('Ala_count ' + str(Ala_count) + '\n')
    file.write('Arg_count ' + str(Arg_count) + '\n')
    file.write('Asn_count ' + str(Asn_count) + '\n')
    file.write('Asp_count ' + str(Asp_count) + '\n')
    file.write('Cys_count ' + str(Cys_count) + '\n')
    file.write('Glu_count ' + str(Glu_count) + '\n')
    file.write('Gln_count ' + str(Gln_count) + '\n')
    file.write('Gly_count ' + str(Gly_count) + '\n')
    file.write('His_count ' + str(His_count) + '\n')
    file.write('Ile_count ' + str(Ile_count) + '\n')
    file.write('Leu_count ' + str(Leu_count) + '\n')
    file.write('Lys_count ' + str(Lys_count) + '\n')
    file.write('Met_count ' + str(Met_count) + '\n')
    file.write('Phe_count ' + str(Phe_count) + '\n')
    file.write('Pro_count ' + str(Pro_count) + '\n')
    file.write('Ser_count ' + str(Ser_count) + '\n')
    file.write('Thr_count ' + str(Thr_count) + '\n')
    file.write('Trp_count ' + str(Trp_count) + '\n')
    file.write('Tyr_count ' + str(Tyr_count) + '\n')
    file.write('Val_count ' + str(Val_count) + '\n')