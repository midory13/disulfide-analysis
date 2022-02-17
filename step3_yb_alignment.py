#!/usr/bin/env python3 


from __future__ import print_function
import os
import re
import argparse

import Bio.PDB 
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt


__author__ = "Yana Bodnar"
__email__ = "yana.bodnar@stud.uni-greifswald.de"
__status__ = "development"


parser = argparse.ArgumentParser(
    description='Script for alignment and RMSD calculation between structures specified on the database file (--file)'  
                +'This script only alignes two structures at one time: if the database entry for the protein contains several PDB-ids, '
                +'one will be selected as the reference structure and all other pdb structures will be compared to this one reference structure'
                +'Caution: Current version of the script does not save the aligned structures, but only writes out the RMSD.'
                +'Please, be aware that this script requires a directory containing *ent files created by BioPython under following path (relative to your current path): "PDB-database/" '
                +'You can refer to the yb_database_analysis.py to create this database directory')
parser.add_argument('-f', '--file', required=False, type=str,
                    help="Specify the input text file containing the PDB-ids. [default = difference.txt] ")
parser.add_argument('-o', '--output', required=False, type=str,
                    help = "Specify the name of the output file. [default = all_local_alignment_results.txt]")
parser.add_argument('-p', '--path', required=False, type=str,
                    help="")
args = parser.parse_args()




parser = Bio.PDB.PDBParser()

def overall_alignment(Protein,ref_PDB, ref_chain, mov_PDB, mov_chain,redox): #Function for the overall alignment of the two pdb structures
    mov_atom = []
    ref_atom = []
    mov_atom_2 = []
    ref_atom_2 = []
    ref_residues_2 = []
    mov_residues_2 = []
    ref_pdbID_file = str(path_to_pdb) +"pdb" + str(ref_PDB.lower()) + ".ent"
    mov_pdbID_file = str(path_to_pdb) +"pdb" + str(mov_PDB.lower()) + ".ent"
    ref_structure = parser.get_structure(Protein, ref_pdbID_file)
    mov_structure = parser.get_structure(Protein, mov_pdbID_file)
    ref_model = ref_structure[0]
    mov_model = mov_structure[0]
    ref_residues = ref_model[ref_chain]
    mov_residues = mov_model[mov_chain]
    print("--")
    for r in ref_residues:
        if r.get_id()[0] == " ":
            #print(r.get_id()[1])
            if r.get_id()[1] in mov_residues:
                if not ref_residues_2:
                    ref_residues_2 = r
                else:
                    ref_residues_2.add(r)
            else:
                continue
        else:
            continue
    for m in mov_residues:
        if m.get_id()[0] == " ":
            #print(m.get_id()[1])
            if m.get_id()[1] in ref_residues:
                if not mov_residues_2:
                    mov_residues_2 = m
                else:
                    mov_residues_2.add(m)
            else:
                continue
        else:
            continue
    for r in ref_residues_2:
            try:
                    ref_a = r["CA"]
                    ref_atom.append(ref_a)
            except:
                    continue
    for r in mov_residues_2:
            try:
                    mov_a = r["CA"]
                    mov_atom.append(mov_a)
            except:
                    continue
    if len(mov_atom) > len(ref_atom):
        del mov_atom[-1]
    elif len(mov_atom) < len(ref_atom):
        del ref_atom[-1]
    print("(((----)))")
    try:
        super_imposer.set_atoms(ref_atom, mov_atom)
        super_imposer.apply(mov_model.get_atoms())
        print(super_imposer.rms)
        print(Protein)
        Protein = Protein.split("  ")
        io = Bio.PDB.PDBIO()
        return(super_imposer.rms, mov_structure, ref_structure)
    except:
        print("could not align the pdbs: " + str(ref_PDB) + ", " + str(mov_PDB))
        rmsd_error = "NO"
        return("Error1", "Error1", "Error1")

def local_alignment(Protein,ref_PDB, ref_chain, mov_PDB, mov_chain,redox, CysRes1, CysRes2): #Function for the local alignment of the two pdb structures +/-5 residues around two given Cys
    mov_atom = []
    ref_atom = []
    mov_atom_2 = []
    ref_atom_2 = []
    ref_residues_2 = []
    mov_residues_2 = []
    min1 = CysRes1-5
    max1 = CysRes1+5
    min2 = CysRes2-5
    max2 = CysRes2+5
    ref_pdbID_file = str(path_to_pdb) +"pdb" + str(ref_PDB.lower()) + ".ent"
    mov_pdbID_file = str(path_to_pdb) +"pdb" + str(mov_PDB.lower()) + ".ent"
    ref_structure = parser.get_structure(Protein, ref_pdbID_file)
    mov_structure = parser.get_structure(Protein, mov_pdbID_file)
    ref_model = ref_structure[0]
    mov_model = mov_structure[0]
    ref_residues = ref_model[ref_chain]
    mov_residues = mov_model[mov_chain]
    print("--")
    for r in ref_residues:
        if r.get_id()[0] == " ":
            if r.get_id()[1] in mov_residues and (r.get_id()[1] in range(min1, max1+1) or r.get_id()[1] in range(min2, max2+1)):
                if not ref_residues_2:
                    ref_residues_2 = r
                else:
                    ref_residues_2.add(r)
            else:
                continue
        else:
            continue
    for m in mov_residues:
        if m.get_id()[0] == " ":
            #print(m.get_id()[1])
            if m.get_id()[1] in ref_residues and (m.get_id()[1] in range(min1, max1+1) or m.get_id()[1] in range(min2, max2+1)):
                if not mov_residues_2:
                    mov_residues_2 = m
                else:
                    mov_residues_2.add(m)
            else:
                continue
        else:
            continue
    for r in ref_residues_2:
            try:
                    ref_a = r["CA"]
                    ref_atom.append(ref_a)
            except:
                    continue
    for r in mov_residues_2:
            try:
                    mov_a = r["CA"]
                    mov_atom.append(mov_a)
            except:
                    continue
    for i in range(0,1):
        if len(mov_atom) > len(ref_atom):
            del mov_atom[-1]
        elif len(mov_atom) < len(ref_atom):
            del ref_atom[-1]
        else:
            continue
        i+=1

    try:
        super_imposer.set_atoms(ref_atom, mov_atom)
        super_imposer.apply(mov_model.get_atoms())
        print(super_imposer.rms)
        print(Protein)
        Protein = Protein.split("  ")
        io = Bio.PDB.PDBIO()
        return(super_imposer.rms, mov_structure, ref_structure)
    except:
        print("could not align the pdbs: " + str(ref_PDB) + ", " + str(mov_PDB))
        rmsd_error = "NO"
        return("Error1", "Error1", "Error1")



def rmsd_calc(Protein, ref_PDB, ref_chain, mov_PDB, mov_chain, CysRes1, CysRes2): #Fuction for the RMSD calculation between two given structures
        mov_atoms = []
        ref_atoms = []
        ref_pdbID_file = str(path_to_pdb) + "pdb" + str(ref_PDB.lower()) + ".ent"
        mov_pdbID_file = str(path_to_pdb) + "pdb" + str(mov_PDB.lower()) + ".ent"
        ref_structure = parser.get_structure(Protein, ref_pdbID_file)
        mov_structure = parser.get_structure(Protein, mov_pdbID_file)
        r = ref_structure[0][ref_chain]
        m = mov_structure[0][mov_chain]
        min1 = CysRes1-5
        max1 = CysRes1+5
        min2 = CysRes2-5
        max2 = CysRes2+5
        print(ref_PDB, mov_PDB)
        for i in range(min1, max1+1):
            #print(i)
            try:
                ref_atom = ()
                X1,Y1,Z1 = r[i]["CA"].get_coord()
                ref_atom = X1,Y1,Z1
                ref_atoms.append(ref_atom)
            except:
                continue
            try:
                mov_atom = ()
                X2,Y2,Z2 = m[i]["CA"].get_coord()
                mov_atom = X2,Y2,Z2
                mov_atoms.append(mov_atom)
            except:
                continue
        #        print(m_a)
        #print('#####')
        for i in range(min2, max2+1):
            try:
                ref_atom = ()
                X1,Y1,Z1 = r[i]["CA"].get_coord()
                ref_atom = X1,Y1,Z1
                ref_atoms.append(ref_atom)
            except:
                continue
            try:
                mov_atom = ()
                X2,Y2,Z2 = m[i]["CA"].get_coord()
                mov_atom = X2,Y2,Z2
                mov_atoms.append(mov_atom)
            except:
                continue
        print(len(ref_atoms))
        print(len(mov_atoms))
        try:
            if len(ref_atoms) >= len(mov_atoms):
                range_X = len(mov_atoms)
                print("ref is bigger")
            elif len(ref_atoms) < len(mov_atoms):
                range_X = len(ref_atoms)
                print("mov is bigger")
            for x in range(1,range_X):
                X1,Y1,Z1 = ref_atoms[x]
                X2,Y2,Z2 = mov_atoms[x]
                sd = (X1-X2)**2 + (Y1-Y2)**2 + (Z1-Z2)**2
                sd += sd
                distance = ((X1-X2)**2 + (Y1-Y2)**2 + (Z1-Z2)**2)**0.5
                    #print(sd)
            rmsd = ( (1/len(ref_atoms) ) * (sd) )**(0.5)
            return(rmsd)
        except:
            print("could not align the pdbs: " + str(ref_PDB) + ", " + str(mov_PDB))
            rmsd_error = "NO"
            return("Error1")
            #print(distance)


#pdbl = PDBList()
super_imposer = Bio.PDB.Superimposer()
io = Bio.PDB.PDBIO()
PDBListOX = []
PDBListRED = []
ProteinDict = {}
if args.file:
    PDB_list_file = args.file
else:
    PDB_list_file = "difference.txt"

if args.output:
    output_file = args.output
else: 
    output_file = "all_local_alignment_results.txt"

if args.path:
    path_to_pdb = args.path
else:
    path_to_pdb = "PDB-database/"


with open(output_file, "w") as file:
    file.write("#Protein Name \t OX/RED? \treference \t target \t Cys1 \t Cys2 \t total RMSD \t local RMSD (+/- 5 residues from Cys) \n" )


with open(PDB_list_file, "r") as file:
    for line in file:
        CysRes_lst = []
        if line.startswith("#"):
            continue
        else:
            PDBListOX = []
            PDBListRED = []
            line = line.split("\t")
            ProteinName = line[0]
            CysRes = line[1]
            print(CysRes)
            CysRes_lst = CysRes.split(",")
            print(CysRes_lst, CysRes_lst[0], type(CysRes_lst[0]))
            if re.search(r"[\s\[\]\']", line[5]):
                line[5] = re.sub(r"[\s\[\]\']", "", line[5])
            #print(line[5])
            if re.search(r"[\s\[\]\']", line[6]):
                line[6] = re.sub(r"[\s\[\]\']", "", line[6])
            if re.search(r",", line[5]):
                PDBListOX = line[5].split(",")
            else:
                PDBListOX.append(line[5])
            if re.search(r",", line[6]):
                PDBListRED = line[6].split(",")
            else:
                PDBListRED.append(line[6])
            first_model = PDBListOX[0]
            match = re.match(r"(\w{4})(\w)", first_model)
            if match:
                ref_model = match[1]
                ref_chain = match[2]
            else:
                ref_model = first_model
                ref_chain = "A"
            mov_OX = PDBListOX[1:]
            mov_red = PDBListRED
            print(PDBListOX)
            print(ProteinName, mov_OX, mov_red)
            print(ref_model)
            print(len(mov_OX))
            if len(mov_OX) != 0:
                print(mov_OX)
                for ox_PDB in mov_OX:
                    match = re.match(r"(\w{4})(\w)", ox_PDB)
                    if match:
                        mov = match[1]
                        mov_chain = match[2]
                    rmsd, mov_str, ref_str = overall_alignment(ProteinName, ref_model, ref_chain, mov, mov_chain, "ox")
                    rmsd_local, mov_str_local, ref_str_local = local_alignment(ProteinName, ref_model, ref_chain, mov, mov_chain, "ox", int(CysRes_lst[0]), int(CysRes_lst[1]))
                    RMSD = rmsd_calc(ProteinName, ref_model, ref_chain, mov, mov_chain, int(CysRes_lst[0]), int(CysRes_lst[1]))
                    print("local RMSD: "+ str(RMSD))
                    with open(output_file, "a") as file:
                        file.write(ProteinName + "\t OX \t" +  ref_model + ref_chain + "\t" + mov + mov_chain + "\t" + str(CysRes_lst[0]) + "\t" + str(CysRes_lst[1]) + "\t" + str(rmsd) + "\t" + str(rmsd_local) + "\n" )
            if len(mov_red) != 0 and len(mov_OX)!=0:
                for red_PDB in mov_red:
                    match = re.match(r"(\w{4})(\w)", red_PDB)
                    if match:
                        mov = match[1]
                        mov_chain = match[2]
                    rmsd, mov_str, ref_str = overall_alignment(ProteinName, ref_model, ref_chain, mov, mov_chain, "red")
                    rmsd_local, mov_str_local, ref_str_local = local_alignment(ProteinName, ref_model, ref_chain, mov, mov_chain, "red", int(CysRes_lst[0]), int(CysRes_lst[1]))
                    RMSD = rmsd_calc(ProteinName, ref_model, ref_chain, mov, mov_chain, int(CysRes_lst[0]), int(CysRes_lst[1]))
                    print("local RMSD: "+ str(RMSD))
                    with open(output_file, "a") as file:
                        file.write(ProteinName + "\t RED \t" +  ref_model + ref_chain + "\t" + mov + mov_chain + "\t" + str(CysRes_lst[0]) + "\t" + str(CysRes_lst[1]) + "\t" + str(rmsd) + "\t" + str(rmsd_local) + "\n" )
