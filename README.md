# Phylovar
Phylovar is a likelihood-based method for joint inference of SNVs and cell lineages from SCS datasets consisting of a large number of loci. It is implemented fully in Python and benefits from the vectorized operations in NumPy to scale up to millions of loci. You can read the manuscript [Here].(https://www.biorxiv.org/content/biorxiv/early/2022/01/18/2022.01.16.476509.full.pdf) 

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
     Since the output of SCIPhi's filtering algorithm rewrites the actual genomic positions in the original mpileup files into a *global* indexing starting from 1 to *N*  (the total number of positions), we wrote a simple script to make a copy of the original mpileup file with global indices. Go to `indexing_scripts` directory, change the `mpileup_path` (path to the original mpileup file) and `out_path` (path to the new mpileup file with global indices): https://github.com/NakhlehLab/Phylovar/blob/e4dc7b32767d0aa7ac7354d3695a3f0c959530f6/indexing_scripts/global_index_conversion.py#L2-L3
     Then, run the code using:
     ```
     python global_index_conversion.py
     ```
     This will output a new mpileup (e.g. `tnbc_global_idx.mpileup`)
   - #### Run `indexer.py` on the mpileup file
     To preserve the original positions, we stored the global index, chromosome index, and the original local index of the positions into a file named `index.csv`. To create this file (which will be used later for recovering the local indices), go to `indexing_scripts`, change the `mpileup_path` (path to the original mpileup file) and `index_path` (path to the csv file containing all the indices): https://github.com/NakhlehLab/Phylovar/blob/8f707b6f552a5058402310c784b1b08f7f928b4b/indexing_scripts/indexer.py#L2-L3
     Then, run the following:
     ```
     python indexer.py
     ```
     This will output the csv file (e.g. `index.csv`)
   - #### Prepare the list of cell names
     Along with the mpileup file, SCIPhi takes as input a list of cell names with their labels. According to the instructions in the Github repository of SCIPhi:
     > run SCIPhI using the cell names provided in cellNames.txt (same order as in the mpileup file). Note that cellNames.txt is a tab delimited file with the cell name in the first column and a cell type identifier in the second column. The cell type can be either CT (tumor cell), CN (control normal cell), or BN (control bulk normal)
     
     The cell names of the TNBC data is provided in `data` directory of this repository named `cellNames.txt`. The diploid cells are tagged with *n* and used as control normal samples, so we labeled them with *CN* 
   - #### Run SCIPhi's statistic test algorithm to filter the non-informative sites
     Given the mpileup file with global indices (say, `tnbc_global_idx.mpileup`) and the list of cell names (say, `cellNames.txt`), go to SCIPhi `build` directory and run SCIPhi using the following command:
     ```
     cat <path to tnbc_global_idx.mpileup> | ./sciphi -o results --in <path to cellNames.txt>
     ```
     This command will generate a csv file named `genotype_matrix.csv` at the same location (`build`) that contains the selected genomic loci after the statistic test.
   - #### Run `local_index_recovery.py`
     To retrieve the actual positions of the genomic sites, you need to run `local_index_recovery.py` given the copy of the original mpileup file with global indices (`tnbc_global_idx.mpileup`), the csv file containing all the indices (`index.csv`), and the output of SCIPhi filtering (`genotype_matrix.csv`). Change the paths to these files plus the path to the output: https://github.com/NakhlehLab/Phylovar/blob/faf19554c148ededcbe05f89ddd1f074b0bbdee2/indexing_scripts/local_index_recovery.py#L2-L5 
     Then, run the following:
     ```
     python local_index_recovery.py
     ```
2. ### Running Phylovar on TNBC data
   After creating the files described in [Filtering the non-informative sites](https://github.com/NakhlehLab/Phylovar/blob/main/README.md#filtering-the-non-informative-sites), we have all the inputs ready for Phylovar. Put the list of cell names (e.g. `cellNames.txt`) and the pre-processed mpileup (e.g. `tnbc_local_idx.mpileup`) into the same folder (like `data` in this repository). Then run the following command to run Phylovar:
   ```
   phylovar.py -names cellNames.txt -out ./data/output/ -indir ./data/ -infile tnbc_local_idx.mpileup -mdthr 1 -mode 0 -niter 100000 -stoch 1 -M 5 -c 10 -verbose yes -errN 100
   ```
   The above command runs Phylovar on `tnbc_local_idx.mpileup` with 10 hill-climbing chains (`-c 10`), on a pool of 5 cores (`-M 5`) each for 100000 iterations (stochastic hill-climbing). After each 100 iterations, the false-positive and false-negative error rates are sampled (`-errN 100`).
3. ### Description of outputs
   
## Contact
If you have any questions, please contact edrisi@rice.edu or edrisi.rice@gmail.com
