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
```python
yb_database_analysis.py -f input.txt
```
### Step 2: averaging the Cys-Cys distances over multiple PDB-structures
### Step 3: structural alignment of ox and red proteins and RMSD calculation
### Step 4: calculation of solvent acessible surface area of Cys
### Step 5: neighbours search


If you use this script, please do not forget to cite us:
Bodnar Y, Lillig CH, 2022, http://rcocf.de/
