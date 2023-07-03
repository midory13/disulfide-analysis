

import re
import pandas as pd
import os
import numpy as np


__author__ = "Yana Bodnar"
__email__ = "yana.bodnar@stud.uni-greifswald.de"
__status__ = "development"










path_out ="/home/sysgen/server/user/YB/Experimental/MM/Hogg_database/Electrostatics/All-proteins/"
#path_out = "/home/yana/Server/user/YB/Experimental/MM/Hogg_database/Electrostatics/All-proteins"
ox_EMD, red_EMD, ox_red_EMD, error_EMD = [], [], [], []

def sort_EMD(protein, A,B, EMD, ox_lst,red_lst, ox_EMD, red_EMD, ox_red_EMD, error_EMD):
    if A in ox_lst and B in ox_lst:
        ox_EMD.append(EMD)
    elif A in red_lst and B in red_lst:
        red_EMD.append(EMD)
    elif (A in ox_lst and B in red_lst) or (A in red_lst and B in ox_lst):
        ox_red_EMD.append(EMD)
    else:
        error_EMD.append(EMD)
    #df2 = pd.DataFrame(data=[[protein, ox_EMD, red_EMD, ox_red_EMD, error_EMD]], columns=["Protein", "OX", "RED", "OX-RED", "ERR"])
    lst = protein, ox_EMD, red_EMD, ox_red_EMD, error_EMD
    #print(i)
    return lst


def create_df(csv_start,path):
    for filename in os.listdir(path):
        if filename.startswith(csv_start):
            filename = os.path.join(path, filename)
            return pd.read_csv(filename)

def sort_pos_neg(protein, pos_or_neg):
    path_to_EMD = path_out +"/" +  protein + "/OutputCompareProteins/"
    if pos_or_neg=="+":
        return  create_df("ListOfEMD_positive_100", path_to_EMD)
    elif pos_or_neg =="-":
        return  create_df("ListOfEMD_negative_100", path_to_EMD)
done_dict = {}
def EMD_df(protein):
    if protein not in done_dict:
        done_dict[protein] = 0 
        print("Working on " + protein)
    if done_dict[protein] == 0:
        ox_EMD, red_EMD, ox_red_EMD, error_EMD = [], [], [], []
        ox_lst = df.loc[df["Protein"]==protein, "pdbID-OX"].tolist()
        red_lst = df.loc[df["Protein"]==protein, "pdbID_RED"].tolist()
        try:
            df_pos = sort_pos_neg(protein, "+")
            print(df_pos)
            print(protein)
            df_pos["Protein"] = protein
            i=0
            POS = df_pos.apply(lambda x: sort_EMD(x["Protein"], x["NameforEMDA"], x["NameforEMDB"], x["EMD"], ox_lst, red_lst,ox_EMD, red_EMD, ox_red_EMD, error_EMD), axis=1,result_type="expand")
            POS = POS.loc[0].tolist()
        except:
            print("Could not find positive")
            POS = [protein, ["","","",""]]
            with open("EMD_error_2023.txt", "a") as fe:
                fe.write(protein + "\tpos\n")
        
        try:
            df_neg = sort_pos_neg(protein, "-")
            df_neg["Protein"] = protein
            ox_EMD, red_EMD, ox_red_EMD, error_EMD = [], [], [], []
            NEG = df_neg.apply(lambda x: sort_EMD(x["Protein"], x["NameforEMDA"], x["NameforEMDB"], x["EMD"], ox_lst, red_lst, ox_EMD, red_EMD, ox_red_EMD, error_EMD ), axis=1,result_type="expand")
            NEG = NEG.loc[0].tolist()
        except:
            print("Could not find negative")
            NEG = [protein, ["","","",""]]
            with open("EMD_error_2023.txt", "a") as fe:
                fe.write(protein + "\tneg\n")
        #POS = [x for x in POS if not re.search(r'amber', str(x))]
        #NEG = [x for x in NEG if not re.search(r'amber', str(x))]
        done_dict[protein] = 1
        return POS, NEG

def mean(lst):
    if len(lst)>0:
        return sum(lst)/len(lst)
    else:
        return "nan"


#df = pd.read_csv("results_el_overview.txt")

df = pd.read_csv("redoing_compare.txt")

df_pos_out = pd.DataFrame(columns = ["Protein", "OX", "RED", "OX-RED", "ERR"])
df_neg_out = pd.DataFrame(columns = ["Protein", "OX", "RED", "OX-RED", "ERR"])

p = df.apply(lambda x: EMD_df(x["Protein"]), axis=1).dropna()
print(p)
p = p.apply(pd.Series)


df_pos_out = pd.DataFrame(p[0].tolist(), columns = ["Protein", "OX", "RED", "OX-RED", "ERR"]).dropna()
df_neg_out = pd.DataFrame(p[1].tolist(), columns = ["Protein", "OX", "RED", "OX-RED", "ERR"]).dropna()


df_pos_out.to_csv("EMD_positive_2023.csv")
df_neg_out.to_csv("EMD_negative_2023.csv")
col_lst = ["OX", "RED", "OX-RED", "ERR"]
col_lst_out = ["OX_mean", "RED_mean", "OX-RED_mean", "ERR_mean"]
col_lst_out_std = ["OX_mean_std", "RED_mean_std", "OX-RED_mean_std", "ERR_mean_std"]
print(df_pos_out)
df_pos_out[col_lst_out] = pd.DataFrame((df_pos_out.apply(lambda x: [np.mean(x[col]) for col in col_lst], axis=1)).tolist())
df_pos_out[col_lst_out_std] = pd.DataFrame((df_pos_out.apply(lambda x: [np.std(x[col]) for col in col_lst], axis=1)).tolist())
df_neg_out[col_lst_out] = pd.DataFrame((df_neg_out.apply(lambda x: [np.mean(x[col]) for col in col_lst], axis=1)).tolist())
df_neg_out[col_lst_out_std] = pd.DataFrame((df_neg_out.apply(lambda x: [np.std(x[col]) for col in col_lst], axis=1)).tolist())

left = df_pos_out[["Protein","OX_mean","OX_mean_std","RED_mean","RED_mean_std","OX-RED_mean","OX-RED_mean_std","ERR_mean","ERR_mean_std"]]
right= df_neg_out[["Protein","OX_mean","OX_mean_std","RED_mean","RED_mean_std","OX-RED_mean","OX-RED_mean_std","ERR_mean","ERR_mean_std"]]

df_total = pd.merge(left,right, on="Protein", suffixes=("_positive","_negative"))

df_total.to_csv("EMD_grouped-by-protein_mean_output_2023.csv")
pos_lst = [sum(df_pos_out[x].tolist(), []) for x in col_lst]
neg_lst = [sum(df_neg_out[x].tolist(), []) for x in col_lst]
df_group_1 = pd.DataFrame(pos_lst).T
df_group_2 = pd.DataFrame(neg_lst).T
df_group = pd.concat([df_group_1, df_group_2], axis=1)
df_group.columns = "OX_+", "RED_+", "OX-RED_+", "ERR_+", "OX_-", "RED_-", "OX-RED_-", "ERR_-"
df_overall = pd.DataFrame()
df_overall["Mean"] = df_group.mean(axis=0)
df_overall["Std"] = df_group.std(axis=0)
print(df_overall)
df_overall.to_csv("EMD_total_mean_output_2023.csv")
