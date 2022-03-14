# Phylovar
Phylovar is a likelihood-based method for joint inference of SNVs and cell lineages from SCS datasets consisting of a large number of loci. It is implemented fully in Python and benefits from the vectorized operations in NumPy to scale up to millions of loci. 

## How to install required packages 
To run Phylovar, you need to install Python (version >2.7). In the recent versions of Python, Numpy, and Scipy are installed as built-in libraries. However, if you do not have them installed, you can install them using the following commands for either pip or conda installation:
### Pip installation
```
pip install numpy
pip install scipy
```
### Conda installation
```
conda install -c anaconda numpy
conda install -c anaconda scipy
```
In addition to basic Python libraries, Phylvar requires Ete3 and Dendropy to be installed on your system. 
### Pip installation
```
pip install ete3
pip install DendroPy
```
### Conda installation
``` conda install -c etetoolkit ete3 ```
``` conda install -c bioconda dendropy```
## Usage
After you downloaded the repository, go to the scripts, and use the following command to make the code executable:
```
chmod +x ./phylovar.py
```
To test the code, in command-line, enter:
```
./phylovar.py --help
```
This will show the description of the arguments. Below you can see an example command for running Phylovar.
```
python main.py -names cellnames.txt -out ./output/ -indir ./input/ 
-infile mpileup.mpileup -mdthr 1 -niter 100000 -stoch yes -M 1 -c 1 -verbose yes -errN 100
```

## Definition of arguments
* -indir : path to the input directory
* -out : path to the output directory
* -names : name of the file containing the cell names (in the input directory)
* -infile : name of the input mpileup file (in the input directory)
* -mdthr : missing data threshold (default value is 1)
* -niter : number of iterations
* -stoch : whether the search is stochastic or not (yes, no)
* -w : number of iterations to wait before termination (for regular hill-climbing search, when stoch=no)
* -M : number of CPUs to use 
* -c : number of hill-climbing chains with different seeds
* -verbose : whether to print the iteration info or not (yes, no)
* -errN : number of iterations per which to sample new error rates


## Reproducibility
Here are the instructions for reproducing the results of the TNBC data containing 32 single-cells from a triple-negative breast cancer patient that is presented in the Phylovar paper.
1. ### Filtering the non-informative sites
   - #### Install SCIPhi with the modified scripts
     To run Phylovar on an example dataset whose results are presented in the paper, you would need to install SCIPhi following the instructions at its Github repository https://github.com/cbg-ethz/SCIPhI. Before installing SCIPhi, copy the modified scripts named `findBesttrees.cpp` and `readData.h` from the directory `sciphi_modified_scripts`. Then go to the `src` directory of SCIPhi and replace the files having the same names in there. These two modified files make SCIPhi to run only the initial statistic test for filtering the non-informative genomic loci without running the entire SCIPhi algorithm.
   - #### Run `global_index_conversion.py` on the mpileup file
     Since the output of SCIPhi's filtering algorithm rewrites the actual genomic positions in the original mpileup files into a *global* indexing starting from 1 to *N*  (the total number of positions), we wrote a simple script to make a copy of the original mpileup file with global indices. Go to `indexing_scripts` directory, change the `mpileup_path` (path to the original mpileup file) and `out_path` (path to the new mpileup file with global indices). Then, run the code using:
     ```
     python global_index_conversion.py
     ```
     This will output a new mpileup (e.g. `tnbc_global_idx.mpileup`)
   - #### Run `indexer.py` on the mpileup file
     To preserve the original positions, we stored the global index, chromosome index, and the original local index of the positions into a file named `index.csv`. To create this file (which will be used later for recovering the local indices), go to `indexing_scripts`, change the `mpileup_path` (path to the original mpileup file) and `index_path` (path to the csv file containing all the indices), and run the following:
     ```
     python indexer.py
     ```
     This will output the csv file (e.g. `index.csv`)
   - #### Prepare the list of cell names
     Along with the mpileup file, SCIPhi takes as input a list of cell names with their labels. According to the instructions in the Github repository of SCIPhi:
     > run SCIPhI using the cell names provided in cellNames.txt (same order as in the mpileup file). Note that cellNames.txt is a tab delimited file with the cell name in the first column and a cell type identifier in the second column. The cell type can be either CT (tumor cell), CN (control normal cell), or BN (control bulk normal)
     The cell names of the TNBC data is provided in `data` directory of this repository named `cellNames.txt`. 
   - #### Run SCIPhi's statistic test algorithm to filter the non-informative sites
     
2. ### Running Phylovar on TNBC data
## Contact
If you have any questions, please contact edrisi@rice.edu or edrisi.rice@gmail.com
