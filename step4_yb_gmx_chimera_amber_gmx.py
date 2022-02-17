#!/usr/bin/env python3 

from asyncio import subprocess
import os
import gmxapi as gmx
from Bio.PDB import *
import re
import os.path
import numpy as np
import sys
import argparse

__author__ = "Yana Bodnar"
__email__ = "yb143526@uni-greifswald.de"
__status__ = "development"





parser = argparse.ArgumentParser(
    description='Script for calculating the solvent acessible surface area of the given Cys residue (total and SASA of the sulfur) and total protein.'+
                'This script uses following software: USCF Chimera (non-gui mode), Gromacs 2020.1 and AmberTools21'+
                'The script is able to mutate Selenomethionine and pyroglutamic acid to methionine and glutamic acid, respectivelly.' + 
                'Only residues standard to Amber99SB-ILDN force field are supported, make sure your protein does not have non-std amino acids' + 
                'All additives from crystallographic experiments (inlc. solvent and ions) will not be considered and do not desturbe the calculations'+
                'While the calculation, the script creates multiple files for every protein entry, including: ')
parser.add_argument('-i', '--input', required=False, type=str,
                    help="Specify the input file. [default = database.dat] The file should contain three columns separated with tab (all following column will be ignored):"+
                    "The first column of the file should contain the name of the protein, second -- PDB-id, third -- residue number." +
                    "If your PDB structure contains several chains, you can select a specific chain by adding its id to the PDB-id: e.g. Chain A of 6h2v -> 6h2vA."+
                    "If you do not specify the chain, the script will only search for the residue in the first chain occuring in the pdb file"+
                    "Please, make sure the protein names do not contain spaces!!")
parser.add_argument('-o', '--output', required=False, type=str,
                    help="Specify the name of the output file [default = SASA_results.dat]")
parser.add_argument('-m', '--mdp', required=False, type=str,
                    help="Specify the mdp file for gromacs calculations. When no mdp file specified, the script will produce the standard mdp \"ions.mdp\"")
parser.add_argument('-op', '--outputpath', required=False, type=str,
                    help="Specify the name of folder for the intermidiate files [default= SASA_calc_output]")
parser.add_argument('-p', '--pdbpath', required=False, type=str,
                    help="Specify path to the folder containing PDB structures, if non existing, the script is able to download the structures from PDB")
args = parser.parse_args()


if args.output:
    sasa_results_file = args.output
else:
    sasa_results_file = "SASA_results.dat"

if args.input:
    database_file = args.input
else:
    database_file = 'database.dat'

if args.mdp:
    mdp_file = args.mdp
else:
    if not os.path.isfile("ions.mdp"):
        with open("ions.mdp", "w") as f1:
            f1.write("integrator  = steep  \n"+
                    "emtol       = 1000.0   \n"+
                    "emstep      = 0.01 \n"+
                    "nsteps      = 50000         \n"+
                    "nstlist         = 1     \n"+
                    "cutoff-scheme	= Verlet  \n"+
                    "ns_type         = grid \n"+
                    "coulombtype     = cutoff \n"+
                    "rcoulomb        = 1.0 \n"+
                    "rvdw            = 1.0 \n"+
                    "pbc             = xyz \n"+
                    ";distance restraints\n"+
                    "disre = simple\n"+
                    "disre-fc = 1000\n")
    mdp_file = "ions.mdp"




if args.outputpath:
    output_path = args.outputpath
else:
    output_path = "SASA_calc_output"

if not os.path.exists(output_path):
    os.makedirs(output_path)

if args.pdbpath:
    pdb_path = args.pdbpath
else:
    pdb_path = "PDB-database/PDB/"

pdb_database_path = "PDB-database/"


error_log = "error.log"
error_short = "error.dat"
missing_log = "missing_proteins.dat"

pdbl = PDBList()
p = PDBParser()
if os.path.isfile(sasa_results_file)!=True:
    with open(sasa_results_file, "w") as file:
        file.write("#protein_name \t pdb_id \t residue \t SASA Cys \t SASA S \t SASA total protein \n")
with open(database_file, "r") as file:
    for line in file:
        if line.startswith("#"):
            continue
        else:
            line = line.strip()
            line = line.split("\t")
            print(line)
            name = line[0]
            if re.search(r"\W", name):
                name = re.sub("\W", "_", name)
            pdb = line[1]
            residue = line[2]
            residue_original = line[2]
            if len(pdb) > 4:
                chainID = pdb[4]
                pdbID = str(pdb[0]) + str(pdb[1]) + str(pdb[2]) +str(pdb[3])
                chain_important = True
            else:
                chain_important = False
                pdbID = pdb
            pdb_file = pdb_path+ pdbID + '.pdb'
            if os.path.isfile(pdb_file)!=True:
                print("Downloading the pdb file for " + str(name))
                pdbl.retrieve_pdb_file(pdbID, pdir = str(pdb_path), file_format = "pdb") 
                print("Download done!")
                pdbID_file = str(pdb_path)+"pdb" + str(pdbID.lower()) + ".ent"
                structure = p.get_structure ("X", pdbID_file)
                io = PDBIO()
                io.set_structure(structure)
                database_path = str(pdb_database_path)+"/PDB/"
                if not os.path.exists(database_path):
                    os.makedirs(database_path)
                    print("created a folder for PDBs")
                io.save(str(database_path)+ pdbID + ".pdb")
                print("working on the " +str(name) + " protein, residue #" + str(residue))
            try:
                coordinate = str(output_path) + "/" + str(name) + "_" + str(pdb) + ".gro"
                topology = str(output_path) + "/" + str(name) + "_" + str(pdb) + ".top"
                tpr = str(output_path) + "/" + str(name) + "_" + str(pdb) + ".tpr"
                posre = str(output_path) + "/" + str(name) + "_" + str(pdb) + "_posre.itp"
                amber_output = str(output_path) + "/" + str(name)+ "_" + str(pdbID) + "_amber.pdb"
                #print(amber)
                ####SEARCH FOR SELENOMETHIONINE(SME) AND PYROGLUTAMIC ACID(PCA)
                MSE_found = False
                PCA_found = False
                with open(pdb_file, "r") as file:
                    for line in file:
                        if re.search(r"MSE", line):
                            MSE_found = True
                        if re.search(r"PCA", line):
                            PCA_found = True
                ###MUTATION OF SME TO MET AND PCA TO GLN
                if MSE_found == True or PCA_found == True:
                    try:
                        chimera_file = str(output_path) + "/" + str(name) + "_" + str(pdb) + "_chimera-script.dat"
                        chimera_output = str(output_path) + "/" + str(name) + "_" + str(pdb) + "_std-aa-only.pdb"
                        with open(chimera_file, "w") as file:
                            if MSE_found == True:
                                file.write("swapaa MET ::MSE \nwrite #0 " + str(chimera_output))
                            if PCA_found == True:
                                file.write("swapaa GLN ::PCA \nwrite #0 " + str(chimera_output))
                        chimera_command = "chimera --nogui " + str(pdb_file) + " 0< " + str(chimera_file) 
                        print("Swaping.....")
                        os.system(chimera_command)
                        print("Done swaping!")
                        pdb_file = chimera_output
                    except:
                        print("Chimera didnt work")
                try:
                    if chain_important == True and os.path.isfile(amber_output) !=True:
                        input_ndx_chain = 'chain ' + chainID +'\n q \n'
                        gmx_ndx_chain = str(output_path) + "/" + str(name) + "_" + str(pdb) + "_chain"+str(chainID) + ".ndx"
                        ndx_chain= gmx.commandline_operation(
                            'gmx', ['make_ndx'],
                            input_files={
                                '-f': pdb_file},
                            output_files={'-o': gmx_ndx_chain},
                            stdin=input_ndx_chain
                        )
                        input_editconf_chain = "ch" + chainID + "\n"
                        gmx_editconf_chain = str(output_path) + "/" + str(name) + "_" + str(pdb) + "_chain"+str(chainID) + ".pdb"
                        editconf_chain= gmx.commandline_operation(
                            'gmx', ['editconf'],
                            input_files={
                                '-f': pdb_file,
                                '-n':gmx_ndx_chain},
                            output_files={'-o': gmx_editconf_chain},
                            stdin=input_editconf_chain)
                        ndx_chain.run()  #runs the command
                        editconf_chain.run()
                        amber_command = " pdb4amber -i " +str(gmx_editconf_chain) +" -o "+ str(amber_output) +" --add-missing-atoms -p "
                        os.system(amber_command)
                    else:
                        if os.path.isfile(amber_output) !=True:
                            amber_command = " pdb4amber -i " +str(pdb_file) +" -o "+ str(amber_output) +" --add-missing-atoms -p "
                            os.system(amber_command)
                except:
                    print("Amber didnt work!")
                    with open(error_log, "a") as f:
                        f.write("Amber didnt work!")
                if os.path.isfile(coordinate)==False or os.path.isfile(topology)==False:
                    pdb2gmx = gmx.commandline_operation(
                        'gmx', ['pdb2gmx'],
                        input_files={
                            '-f': amber_output,
                            '-water': 'tip3p',
                            '-ignh': '1',
                            '-ff':'amber99sb-ildn'},
                        output_files={'-o': coordinate, '-p': topology, '-i': posre},
                    )
                    if pdb2gmx.output.returncode.result() !=0: #checks wheather there was an error
                    #print(editconf.output.erroroutput.result())
                        print(pdb2gmx.output.erroroutput.result()) 
                if os.path.isfile(tpr) == False:
                    grompp = gmx.commandline_operation(
                        'gmx', ['grompp'],
                        input_files={
                            '-c': coordinate,
                            '-f': mdp_file,
                            '-p': topology},
                        output_files={'-o': tpr})
                    if grompp.output.returncode.result() !=0:
                        print(grompp.output.erroroutput.result())
                renumbering = str(output_path) + "/" + str(name)+ "_" + str(pdbID) + "_amber_renum.txt"
                if os.path.isfile(renumbering):
                    found_new_number= False
                    with open(renumbering, "r") as file_r:
                        for line_r in file_r:
                            line_r = line_r.strip()
                            line_r = line_r.split()
                            if line_r[1] == residue and found_new_number == False:
                                residue = line_r[3]
                                found_new_number = True
                #-------------------------------------------------------------------------------------------
                #-------------------------------------------------------------------------------------------
                #-------------------------------------------------------------------------------------------
                #### MAKE NDX for residue selection #########
                ##-------------------------------------------------------------------------------------------
                input_ndx = 'r' + residue +'\n' + 'r'+ residue + '&tS* \n q \n'
                gmx_ndx_cys = str(output_path) + "/" + str(name) + "_" + str(pdb) + "_Cys"+str(residue_original) + ".ndx"
                ndx_cys= gmx.commandline_operation(
                    'gmx', ['make_ndx'],
                    input_files={
                        '-f': tpr},
                    output_files={'-o': gmx_ndx_cys},
                    stdin=input_ndx
                )
                ndx_cys.run()  #runs the command

                if ndx_cys.output.returncode.result() !=0: #checks wheather there was an error
                    print(ndx_cys.output.erroroutput.result())
                #-------------------------------------------------------------------------------------------
                #-------------------------------------------------------------------------------------------
                cystein = 'r_' + str(residue)
                sulfur = '\"r_'+str(residue)+'_&_S*\"'
                cystein_out = str(output_path) + "/" + str(name) + "_" + str(pdb) + "_Cys"+str(residue_original) + "_total.xvg"
                sulfur_out = str(output_path) + "/" + str(name) + "_" + str(pdb) + "_Cys"+str(residue_original) + "_S.xvg"
                sasa_cys= gmx.commandline_operation(
                    'gmx', ['sasa'],
                    input_files={
                        '-f': coordinate,
                        '-n': ndx_cys.output.file['-o'],
                        '-s': tpr,
                        '-surface': 'Protein',
                        '-output': cystein},
                    output_files={'-o': cystein_out}
                )

                sasa_s= gmx.commandline_operation(
                    'gmx', ['sasa'],
                    input_files={
                        '-f': coordinate,
                        '-n': ndx_cys.output.file['-o'],
                        '-s': tpr,
                        '-surface': 'Protein',
                        '-output': sulfur
                    },
                    output_files={'-o': sulfur_out}
                )

                if sasa_cys.output.returncode.result() !=0 or sasa_s.output.returncode.result() !=0:
                    print(sasa_cys.output.erroroutput.result())
                    print(sasa_s.output.erroroutput.result())



                ###--------------------------------READ-XVG-FILE-------------------------------

                xvg_file_cystein = np.loadtxt(cystein_out, comments="@", skiprows=14)
                xvg_file_sulfur = np.loadtxt(sulfur_out, comments="@", skiprows=14)
                with open(sasa_results_file, "a") as file:
                    file.write(str(name) + '\t' + str(pdb) + '\t'+ str(residue_original)+ '\t' + str(xvg_file_cystein[2]) + '\t'+ str(xvg_file_sulfur[2]) +'\t' + str(xvg_file_cystein[1]) + '\n')
                ####----------------------------------------------------------------------------  
            except:
                with open(missing_log, "a") as file:
                    file.write(str(name) + '\t' + str(pdb) + '\t'+ str(residue_original)+ '\n')
                with open(error_short, "a") as file:
                    file.write("I had problems with protein: " + str(name)+ " pdbID: " + str(pdb) + " residue " + str(residue)+"\n")
                    file.write("\n __________________________________________________________\n")
                    print("I had problems with protein: " + str(name)+ " pdbID: " + str(pdb) + " residue " + str(residue)+"\n")
                    print("\n __________________________________________________________\n")
                with open(error_log, "a") as file:
                    file.write("I had problems with protein: " + str(name)+ " pdbID: " + str(pdb) + " residue " + str(residue)+"\n")
                    try:
                        file.write(pdb2gmx.output.erroroutput.result()) 
                        file.write(grompp.output.erroroutput.result())
                        file.write(ndx_cys.output.erroroutput.result())
                        file.write(sasa_cys.output.erroroutput.result())
                        file.write(sasa_s.output.erroroutput.result())
                    except:
                        continue
                    file.write("\n __________________________________________________________\n")

