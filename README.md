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
### Input files:
The examples of needed input filed can be found under example_input.
The example of the database file for step1 represents a part of the table S2 from *Pijning et al 2018 R. Soc. open sci.*
The inputfiles for steps 4 to 6 are identicall.

### Step 1: calculation of the Cys-Cys distances 
In this step, using the `step1_yb_database_analysis.py` script, it is possible to calculate the distance between two cysteins in reduced and oxidised crystallographic structures from Protein Data Bank (PDB).
For calculation, the script requires an excel file (for now only .xlsx extension is supported), which would contain following columns with the exact titles (columns with other headers will be ignored):
|PDBs with Bond | PDB missing Bond| PDB Cys1 | PDB Cys2 | Compound | Organism |
| --- | --- | --- | --- | --- | --- |

Example of usage (will calculate the distance between CA atoms of Cys specified in Sheet1 of the database file):
```bash
step1_yb_database_analysis.py -f Step1_database_for_analysis.xlsx -s Sheet1 -o step1_results.txt --atom CA
```
### Step 2: averaging the Cys-Cys distances over multiple PDB-structures
In this step, using the `step2_yb_SS_analysis.py` script, we can further analyse the results of the calculations from Step 1. It averages the distances between Cysteins in oxidised and reduced structures of one protein (which was required for proteins with several PDB structures) and saves out the results of the calcultions into two output files: `difference.txt` and `biggest_difference.txt`. In the latter, only proteins with the distance differences above 2 Angstrom are printed out.
Furthermore, this script creates histograms (both in .svg and .png formats) containing descributions of occurence of the distances between S-S of the Cysteins and distribution of the distance differences between ox and red proteins. The examples of all these output files can be found under example_output/Step2

Example of usage (will calculate the distance between CA atoms of Cys specified in Sheet1 of the database file):
```bash
step2_yb_SS_analysis.py -f step1_results.txt
```

### Step 3: structural alignment of ox and red proteins and calculation of root-mean-square-diviation (RMSD)
Using the `difference.txt` as an input file for `step3_yb_alignment.py`, the structures of one protein in oxidised and reduced forms can be aligned. This script does two types of alignment: 
1) overall (or total): based on the backbone of the whole protein (output: total RMSD values) 
2) local: based on +/-5 residues from Cysteines specified in the input file (output: local RMSD values)

Important to note that the current version of the sctipr:
- only alignes two structures at one time: if the database entry for the protein contains several PDB-ids, one oxidised structure will be selected as the reference and all other pdb structures will be compared to it
- does not save the aligned structures, but only writes out the RMSD
- needs the PDB-database created in previous steps containing .ent files

### Step 4: calculation of solvent acessible surface area of Cys
`step4_yb_gmx_chimera_amber_gmx.py` is a script for calculating the solvent acessible surface area of the given Cys residue (total and SASA of the sulfur) and total protein. This script uses following software: [USCF Chimera](https://www.cgl.ucsf.edu/chimera/)(non-gui mode), [Gromacs 2020.1](https://manual.gromacs.org/documentation/2020.1/gmxapi/index.html) and [AmberTools21](https://ambermd.org/AmberTools.php). 

The script is able to mutate Selenomethionine and pyroglutamic acid to methionine and glutamic acid, respectivelly. 
Only residues standard to Amber99SB-ILDN force field are supported, make sure your protein does not have non-std amino acids. All additives from crystallographic experiments (inlc. solvent and ions) will not be considered and do not disturbe the calculations. 

While the calculation, the script creates large amount of intermediate files for every protein entry. For those files, the script conveniatelly creates a folder, name of which can be specified under `--outputpath`. **Note**, that the output file containing the calculation overview is saved in the working directory!

### Steps 5 and 6: neighbours search
The `step5_yb_neighbor-search.py` script is there for you to find all residues that are neighbouring the sulfur of the target cystein (or methionine) in specified radius (default: 3.5A)
For each entry, this script create a file with information about contacts found in that radius.
Since it creates a single file for each entry, it is advisible to create a separate folder for these files (can be specified under --folder)
The so created contact files can be further evaluated via `step6_yb_neighbor-search_part2.py` script
Important note: the atom types of the sulfur that are supported are: SG, S, SG1, SG2, SD; other types will not be recognised

If you use this script, please do not forget to cite us:
Bodnar Y, Lillig CH, 2022, http://rcocf.de/
