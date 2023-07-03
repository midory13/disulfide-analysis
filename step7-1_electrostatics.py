#! /usr/bin/ python3


from __future__ import print_function
import pandas as pd
import numpy as np
import shutil
import os
import sys
import logging
#import Biopython
import subprocess
from Bio.PDB import *
from datetime import datetime
import re

__author__ = "Yana Bodnar"
__email__ = "yana.bodnar@uni-greifswald.de"
__status__ = "development"


""" This script requires two following installations for functioning: 1) MutComp https://github.com/WillyBruhn/MutComp and 2)ComparingProteins https://github.com/BerensF/ComparingProteins 3) UCSF chimera https://www.cgl.ucsf.edu/chimera/"""


file_handler = logging.FileHandler("program.log")
stream_handler = logging.StreamHandler(sys.stdout)
parser = PDBParser()
io = PDBIO()

#####################GLOBAL VARIABLES######################
path_R ="/home/sysgen/server/user/YB/Experimental/electrostatics_/ComparingProteins_FB/ComparingProteins-master/EMDandClustering/AllLowerB_EMD_Clust.R"
path_cpp = "/home/sysgen/server/user/YB/Experimental/electrostatics_/ComparingProteins_FB/ComparingProteins-master/LowerBounds/FirstLowerBound/main"
fixing_on = True
start_time = datetime.now()
input_file_name = "results"
chimera = "/home/sysgen/.local/UCSF-Chimera64-1.16/bin/chimera"

input_file_txt = input_file_name + ".txt"
input_file_pck = input_file_name + ".pkl"
M = 5
N = 10
path_db = "/home/sysgen/server/user/YB/Experimental/MM/Hogg_database/PDB-database/PDB/"
path_out = "/home/sysgen/server/user/YB/Experimental/MM/Hogg_database/Electrostatics/All-proteins/"
##################################################################


###############################FUNCTIONS############################################
########################################################################################
########################################################################################
########################################################################################
def read_database(csv):
    df = pd.read_csv(csv, sep="\t",index_col=False)
    print(df.head())
    protein = df.Protein
    pdb_ox = df["pdbID-OX"]
    pdb_red = df["pdbID_RED"]
    ca_ox = df["Ca-Ca-distance-OX"]
    ca_red = df["Ca-Ca-distance-RED"]
    df.Protein = protein.str.replace('\s', '',regex=True)
    df.Protein = df.Protein.str.replace('[^a-zA-Z0-9]', '_',regex=True)
    df.drop(df[ca_ox.str.startswith("Could not",na=False)].index, inplace = True)
    df.drop(df[ca_red.str.startswith("Could not",na=False)].index, inplace = True)
    df.drop("CysRes1", inplace=True, axis=1)
    df.drop("CysRes1.1", inplace=True, axis=1)
    df.drop("Ca-Ca-distance-OX", inplace=True, axis=1)
    df.drop("Ca-Ca-distance-RED", inplace=True, axis=1)
    df["PDB"] = pdb_ox.combine_first(pdb_red)
    pdb = df["PDB"]
    protein = df.Protein
    return df, protein, pdb_ox, pdb_red, pdb


def check_MSE(pdb_file):
    MSE_found = False
    PCA_found = False
    with open(pdb_file, "r") as file:
        for line in file:
            if re.search(r"MSE", line):
                MSE_found = True
            if re.search(r"PCA", line):
                PCA_found = True
    if MSE_found == True or PCA_found == True:
        return True
    else:
        return False
    ###MUTATION OF SME TO MET AND PCA TO GLN
def mutate_MSE(pdb_file, chimera_output, output_path):
    if check_MSE(pdb_file) ==True:
        try:
            pdb_file = path_db + pdbid + ".pdb"
            pdb = str(pdb_file)
            chimera_file = pdb[:4] + "_chimera-script.dat"
            with open(chimera_file, "w") as file:
                if MSE_found == True:
                    file.write("swapaa MET ::MSE \nwrite #0 " + str(chimera_output))
                if PCA_found == True:
                    file.write("swapaa GLN ::PCA \nwrite #0 " + str(chimera_output))
            chimera_command = str(chimera) + " --nogui " + str(pdb_file) + " 0< " + str(chimera_file) 
            print("Swaping.....")
            subprocess.run(chimera_command,shell=True)
            print("Done swaping!")
        except:
            print("Chimera didnt work")

def fix_pdb(pdb_file):
    new_pdb_file = pdb_file[:-4] + "_amber.pdb"
    amber_command = " pdb4amber -i " +str(pdb_file) +" -o "+ str(new_pdb_file) +" --add-missing-atoms -p "
    subprocess.run(amber_command, shell=True)
    subprocess.run(["mv", new_pdb_file, pdb_file])

def search_pdb(protein, pdbid, df):
    #search for the PDB-ID in the path, copy it to a folder named protein, if such does not exists
    index = df[(df.Protein == protein) & (df.PDB == pdbid)].index
    if len(pdbid) > 4:
        pdbid = pdbid[:4]
    pdb_file = pdbid + ".pdb"
    pdb_db = path_db + pdbid + ".pdb"
    print(pdb_db)
    if os.path.exists(pdb_db)==False:
        print("Could not find the pdb in the database...")
        df.at[index, "PDB_exists"] = 0 
    else:
        df.at[index, "PDB_exists"] = int(1) 
    df.to_pickle(input_file_pck)    
    df.to_csv("test_output_df.txt") 
    return df

def copy_pdb(protein, pdbid):
    protein_dir = os.path.join(path_out, protein)
    input_dir = protein + "/Input/"
    pdb_file = pdbid + ".pdb"
    index = df[(df.Protein == protein) & (df.PDB == pdbid)].index
    directory_out = os.path.join(path_out, input_dir)
    if os.path.exists(protein_dir)==False:
        os.mkdir(protein_dir)
    if os.path.exists(directory_out)==False:
        os.mkdir(directory_out)
    if len(pdbid) > 4:
        #split with BioPython
        dst = os.path.join(directory_out, pdb_file)
        chain = pdbid[4]
        pdb_file_short = pdbid[:4] + ".pdb"
        pdb_file_short_path = os.path.join(path_db, pdb_file_short)
        structure=parser.get_structure("X", pdb_file_short_path)
        correct_chain = structure[0][chain]
        io.set_structure(correct_chain)
        io.save(dst)

        print("Saving PDB chain to the Input folder.... ")
    else:
        src = os.path.join(path_db, pdb_file)
        dst = os.path.join(directory_out, pdb_file)
        mutate_MSE(src, dst, directory_out)
        if os.path.exists(dst) == False:
            shutil.copy(src, dst)
    mutate_MSE(dst, dst, directory_out)
    if fixing_on == True:
        fix_pdb(dst)    
    df.at[index, "Copied_files"] = int(1) 
    df.to_pickle(input_file_pck)    
    df.to_csv("test_output_df.txt") 
    return directory_out

def mutcomp(path_to_protein, protein, pdbid):
    #runs mutcomp for the protein of interest, pdbs found under path/protein/Input
    print("Working on the electrostatics calculations of " + protein + " " + pdbid + "...............")
    output_dir = str(protein + "/")
    directory_out = str(os.path.join(path_out, output_dir))
    index = df[(df.Protein == protein)].index
    cmd = "Rscript /home/sysgen/server/user/Z_Final_Backups_Alumni/20191217_final_Backup_WLB/PredictingProteinInteractions/Classification/./predictingProteinInteractions.R --mode SingleDistance --doClustering TRUE  --pdb_folder " + path_to_protein + " --doMutComp TRUE --q_val 1 --labels_train /home/sysgen/server/user/Z_Final_Backups_Alumni/20191217_final_Backup_WLB/PredictingProteinInteractions/data/106Test/labels.txt --numberOfPoints 4 --rounds 10 --MutCompOutPut " + directory_out
    subprocess.run([cmd], shell=True)
    df.at[index, "MutComp_done"] = int(1)     
    print("Done with electrostatics of "+ protein)
    print("Paths:")
    print(directory_out)
    print(path_to_protein)
    directory_out = str(os.path.join(directory_out, "Output/"))
    df.to_pickle(input_file_pck)
    df.to_csv("test_output_df.txt")
    return directory_out

def mutcomp_output_dir(protein):
    output_dir = str(protein + "/Output/")
    directory_out = str(os.path.join(path_out, output_dir))
    return directory_out

def compare(path, protein, pdbid):
    #creates a parameter file for the CompareProteins
    #runs CompareProteins
    index = df[(df.Protein == protein) & (df.PDB == pdbid)].index.tolist()
    #print(index)
    if df.at[index[0], "CompareProteins_done"] == 0:
        parameter_file = os.path.join(path, "parameter.txt")
        path_in = path
        print("-----")
        print("Working on the protein ", protein)
        print("-----")
        outputpath = "OutputCompareProteins"
        path_out = os.path.join(path[:-8], outputpath)
        if os.path.exists(path_out)==False:
            os.mkdir(path_out)
        with open(parameter_file, "w") as file:
            file.write("# Parameter file for the clustering of protein by there isosurfaces \n# comments are marked with an '#' \n# '=' is a seperator between paramtername and parameter \n# \n# Full Path to the protein files \nPathToData= " + path_in +" \n# Full path to the directory where the output should be stored \nPathToOutput=" + path_out + "\n# Number of points to select, normally 100 \nn = " + str(N) +" \n# Numbe of of rounds, normaly 500 \nm = " + str(M) + " \n############################################################################################ \n# Path to the program, where the C++ and R code is stored, only change if the file is moved \nPathToRProgram= " + path_R + " \nPathToCPPProgram= "+ path_cpp + " \n################################ \n")
        subprocess.run(["Rscript /home/sysgen/server/user/YB/Experimental/electrostatics_/ComparingProteins_FB/ComparingProteins-master/./CompareIsosurfaces.R " + parameter_file], shell=True)
        index = df[(df.Protein == protein)].index
        print(index)
        df.at[index, "CompareProteins_done"] = int(1)   
        df.to_pickle(input_file_pck)
        df.to_csv("test_output_df_1.txt")
        mutcomp_output_dir(protein)
        print("-----")


def postprocessing(pdb_ox, pdb_red, protein, output_compare):
    positive = "ListOfEMD_positive_" + n + ".csv"
    negative = "ListOfEMD_negative_" + n + ".csv"
    pos_path = os.path.join(output_compare, positive)
    neg_path = os.path.join(output_compare, negative)
    
##################################################################
##################################################################
##################################################################




###########################################MAIN############################################
df, protein, pdb_ox, pdb_red, pdb = read_database(input_file_txt)
if os.path.exists(input_file_pck)== False:
    df["MutComp_done"] = int(0)
    df["PDB_exists"] = int(0)
    df["Copied_files"] = int(0)
    df["CompareProteins_done"] = int(0)
    df.to_pickle(input_file_pck)
columns = ['Protein', 'PDB', 'pdbID-OX', 'pdbID_RED', 'MutComp_done', 'PDB_exists', 'Copied_files','CompareProteins_done', 'Input_dir']
df_pck = pd.read_csv("test_output_df.txt", index_col=0, usecols=columns)
#df_pck = pd.read_pickle(input_file_pck)
df_txt = df[["Protein", "PDB", "pdbID-OX","pdbID_RED"]]
print(df_txt.head())
df = pd.merge(df_txt, df_pck, how="left", on = ["PDB", "Protein", "pdbID-OX","pdbID_RED"], suffixes=("", ""))
df[["MutComp_done", "PDB_exists", "Copied_files", "CompareProteins_done"]] = df[["MutComp_done", "PDB_exists", "Copied_files", "CompareProteins_done"]].fillna(int(0))
df.to_pickle(input_file_pck)
df.to_csv("test_output_df.txt")
print(df.head())
df.apply(lambda x: search_pdb(x["Protein"], x["PDB"], df) if (x["PDB_exists"] == 0 and x["Copied_files"] == 0) else None, axis=1)
print(type(df))
df["Input_dir"] = df.apply(lambda x: copy_pdb(x["Protein"], x["PDB"]) if (x["PDB_exists"] == 1 and x["Copied_files"] == 0) else x["Input_dir"], axis=1)
end_time1 = datetime.now()
print('Copying Duration: {}'.format(end_time1 - start_time))
print(df.head())

print("__________________________________________________________________________________________________________________________________")
df["Output_dir"] = df.apply(lambda x: mutcomp(x["Input_dir"], x["Protein"], x["PDB"]) if (x["MutComp_done"] == 0) else mutcomp_output_dir(x["Protein"]), axis=1)
df.to_pickle(input_file_pck)
df.to_csv("test_output_df.txt")
print(df.head())
end_time2 = datetime.now()
print('MutComp Duration: {}'.format(end_time2 - end_time1))
print("__________________________________________________________________________________________________________________________________")
df.apply(lambda x: compare(x["Output_dir"], x["Protein"], x["PDB"]) if (x["MutComp_done"]==1 and x["CompareProteins_done"] == 0) else None, axis=1)

end_time3 = datetime.now()
print('Compare Duration: {}'.format(end_time3 - end_time2))
print('Total Duration: {}'.format(end_time3 - start_time))
