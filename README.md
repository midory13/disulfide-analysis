# disulfide-analysis
Collection of Python3 scripts for analysing and comparison of disulfides in oxidised and reduced proteins.
## Software requirements:
- python3 modules:
  - BioPDB (Biopython)
  - pandas
  - NumPy
  - Matplotlib
- [Gromacs (gmxapi)](https://manual.gromacs.org/documentation/2020.1/gmxapi/index.html)
- [AmberTools](https://ambermd.org/AmberTools.php)
- [UCSF Chimera](https://www.cgl.ucsf.edu/chimera/)


## How to...
### Step 1: calculation of the Cys-Cys distances 
In this step, using the `step1_yb_database_analysis.py` script, it is possible to calculate the distance between two cysteins in reduced and oxidised crystallographic structures from Protein Data Bank (PDB).
For calculation, the script requires an excel file (for now only .xlsx extension is supported), which would contain following columns with the exact titles (columns with other headers will be ignored):
|PDBs with Bond | PDB missing Bond| PDB Cys1 | PDB Cys2 | Compound | Organism |
| --- | --- | --- | --- | --- | --- |

The example of such file can be found in the `example_input` file, which represents a part of the table S2 from *Pijning et al 2018 R. Soc. open sci.*

Example of usage (will calculate the distance between CA atoms of Cys specified in Sheet1 of the database file):
```python3
step1_yb_database_analysis.py -f Step1_database_for_analysis.xlsx -s Sheet1 -o step1_results.txt --atom CA
```
### Step 2: averaging the Cys-Cys distances over multiple PDB-structures
In this step, using the `step2_yb_SS_analysis.py` script, we can further analyse the results of the calculations from Step 1. It averages the distances between Cysteins in oxidised and reduced structures of one protein (which was required for proteins with several PDB structures) and saves out the results of the calcultions into two output files: `difference.txt` and `biggest_difference.txt`. In the latter, only proteins with the distance differences above 2 Angstrom are printed out.
Furthermore, this script creates histograms (both in .svg and .png formats) containing descributions of occurence of the distances between S-S of the Cysteins and distribution of the distance differences between ox and red proteins. The examples of all these output files can be found under example_output/step2

### Step 3: structural alignment of ox and red proteins and RMSD calculation
### Step 4: calculation of solvent acessible surface area of Cys
### Step 5: neighbours search


If you use this script, please do not forget to cite us:
Bodnar Y, Lillig CH, 2022, http://rcocf.de/
