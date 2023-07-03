 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy


__author__ = "Yana Bodnar"
__email__ = "yana.bodnar@uni-greifswald.de"
__status__ = "development"


df = pd.read_csv("EMD_grouped-by-protein_mean_output_20230213.csv", index_col = 0)
df_lst =['Protein', 'OX_mean_positive',  'RED_mean_positive',  'OX-RED_mean_positive',  'OX_mean_negative',  'RED_mean_negative',  'OX-RED_mean_negative']

#
df = df[df_lst]

print(df.columns.values.tolist())
with open("statistics/20230224_EMD_summary-stat.dat", "w") as f1:
    f1.write(str(df.agg(["count", "min", "max", "median", "mean", "skew"])))
    print(df.agg(["count", "min", "max", "median", "mean", "skew"]))



with open("statistics/20230224_EMD_stat_two-sided.dat", "w") as f2:
    f2.write(str("POS OX vs RED\n"))
    print(scipy.stats.mannwhitneyu(x=df["OX_mean_positive"].dropna(), y = df['RED_mean_positive'].dropna(), alternative='two-sided'))
    f2.write(str(scipy.stats.mannwhitneyu(x=df["OX_mean_positive"].dropna(), y = df['RED_mean_positive'].dropna(), alternative='two-sided')))

    f2.write(str("\nPOS OX vs OX-RED\n"))
    print(scipy.stats.mannwhitneyu(x=df["OX_mean_positive"].dropna(), y = df['OX-RED_mean_positive'].dropna(), alternative='two-sided'))
    f2.write(str(scipy.stats.mannwhitneyu(x=df["OX_mean_positive"].dropna(), y = df['OX-RED_mean_positive'].dropna(), alternative='two-sided')))

    f2.write(str("\nPOS RED vs OX-RED\n"))
    print(scipy.stats.mannwhitneyu(x=df["RED_mean_positive"].dropna(), y = df['OX-RED_mean_positive'].dropna(), alternative='two-sided'))
    f2.write(str(scipy.stats.mannwhitneyu(x=df["RED_mean_positive"].dropna(), y = df['OX-RED_mean_positive'].dropna(), alternative='two-sided')))

    f2.write(str("\nPOS OX vs RED\n"))
    print(scipy.stats.mannwhitneyu(x=df["OX_mean_negative"].dropna(), y = df['RED_mean_negative'].dropna(), alternative='two-sided'))
    f2.write(str(scipy.stats.mannwhitneyu(x=df["OX_mean_negative"].dropna(), y = df['RED_mean_negative'].dropna(), alternative='two-sided')))

    f2.write(str("\nPOS OX vs OX-RED\n"))
    print(scipy.stats.mannwhitneyu(x=df["OX_mean_negative"].dropna(), y = df['OX-RED_mean_negative'].dropna(), alternative='two-sided'))
    f2.write(str(scipy.stats.mannwhitneyu(x=df["OX_mean_negative"].dropna(), y = df['OX-RED_mean_negative'].dropna(), alternative='two-sided')))

    f2.write(str("\nPOS RED vs OX-RED\n"))
    print(scipy.stats.mannwhitneyu(x=df["RED_mean_negative"].dropna(), y = df['OX-RED_mean_negative'].dropna(), alternative='two-sided'))
    f2.write(str(scipy.stats.mannwhitneyu(x=df["RED_mean_negative"].dropna(), y = df['OX-RED_mean_negative'].dropna(), alternative='two-sided')))
