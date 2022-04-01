# Phylovar
## Towards scalable phylogeny-aware inference of single-nucleotide variations from single-cell DNA sequencing data
Phylovar is a likelihood-based method for joint inference of SNVs and cell lineages from SCS datasets consisting of a large number of loci. It is implemented fully in Python and benefits from the vectorized operations in NumPy to scale up to millions of loci. You can read the manuscript [Here](https://www.biorxiv.org/content/biorxiv/early/2022/01/18/2022.01.16.476509.full.pdf).

## Description of the directories
   - `src`: contains all the source codes of Phylovar
   - `data`: contains an example dataset (referred to as TNBC data in the paper). It includes the necessary files for reproducing the result of Phylovar on this data as presented in the paper (the list of cell names and pre-processed mpileup) and a zipped file containing the outputs of Phylovar.
   - `indexing_scripts`: contains the three scripts needed for processing the mpileup files before and after SCIPhi's initial filtering.
   - `sciphi_modified_scripts`: contains two scripts we modified in order to utilize SCIPhi's filtering algorithm in our work.

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
   Here are the instructions for filtering the non-informative sites from a given mpileup file, before passing it as input to Phylovar. Not all genomic loci contain mutations, moreover, Phylovar has a limit on the number of loci so we need to filter out the unnecessary loci from the original mpileup. Note, filtering out non-informative sites can be done in other ways too and it is up to the user how to prepare the mpileup file for Phylovar. Here, we are just describing the technique we used in our original analysis. 
   - #### Install SCIPhi with the modified scripts
     To run Phylovar on an example dataset whose results are presented in the paper, you would need to install SCIPhi following the instructions at its Github repository https://github.com/cbg-ethz/SCIPhI. Before installing SCIPhi, copy the modified scripts named `findBesttrees.cpp` and `readData.h` from the directory `sciphi_modified_scripts`. Then go to the `src` directory of SCIPhi and replace the files having the same names in there. These two modified files make SCIPhi to run only the initial statistic test for filtering the non-informative genomic loci without running the entire SCIPhi algorithm.
   - #### Run `global_index_conversion.py` on the original mpileup file
     Since the output of SCIPhi's filtering algorithm rewrites the actual genomic positions in the original mpileup files into a *global* indexing starting from 1 to *N*  (the total number of positions), we wrote a simple script to make a copy of the original mpileup file with global indices. Go to `indexing_scripts` directory, and run `global_index_conversion.py` by specifying two arguments, namely, the `-mpileup` which is path to the original mpileup file and `-out` which is path to the new mpileup file with global indices. For example:
     ```
     python global_index_conversion.py -mpileup ./tnbc.mpileup -out ./tnbc_global_idx.mpileup
     ```
     This will output a new mpileup (e.g. `tnbc_global_idx.mpileup`)
   - #### Prepare the list of cell names
     Along with the mpileup file, SCIPhi takes as input a list of cell names with their labels. According to the instructions in the Github repository of SCIPhi:
     > run SCIPhI using the cell names provided in cellNames.txt (same order as in the mpileup file). Note that cellNames.txt is a tab delimited file with the cell name in the first column and a cell type identifier in the second column. The cell type can be either CT (tumor cell), CN (control normal cell), or BN (control bulk normal)
     
     The cell names of the TNBC data is provided in `data` directory of this repository named `cellNames.txt`. The tumor cells are tagged with *CT* and the diploid cells with prefix *n* are used as control normal samples, so we labeled them with *CN*. Note that we you are going to run Phylovar on healthy cells, you would need to tag them with *CT*, otherwise, SCIPhi considers all of them as control normals and outputs an empty file.
   - #### Run SCIPhi's statistic test algorithm to filter the non-informative sites
     Given the mpileup file with global indices (say, `tnbc_global_idx.mpileup`) and the list of cell names (say, `cellNames.txt`), go to SCIPhi `build` directory and run SCIPhi using the following command:
     ```
     cat <path to tnbc_global_idx.mpileup> | ./sciphi -o results --in ./cellNames.txt
     ```
     This command will generate a csv file named `genotype_matrix.csv` at the same location (`build`) that contains the selected genomic loci after the statistic test. After this step, run the following command on `genotype_matrix.csv` to remove the duplicated lines in the file:
     ```
     cut -f1-17 -d',' genotype_matrix.csv | uniq > genotype_matrix_new.csv
     ```
     Since the number of cancer cells in the TNBC data is 16, we used 17 (16+1) as the number of fields in the above command (`NF==17{print}{}`). The general form of the above command is:
     ```
     cut -f1-<number of cancer cells + 1> -d',' genotype_matrix.csv | uniq > genotype_matrix_new.csv
     ```
   - #### Run `local_index_recovery.py`
     To retrieve the actual positions of the genomic sites, you need to run `local_index_recovery.py`. `local_index_recovery.py` works with three arguments, `-mpileup` which is the original mpileup file with actual indices (`tnbc.mpileup`), `-sciphi` which is the output of SCIPhi filtering (`genotype_matrix_new.csv` after removing the duplocated lines), and `-out` which is the path to the output. The following is an example command to run this code:
     ```
     python local_index_recovery.py -mpileup ./tnbc.mpileup -sciphi ./genotype_matrix_new.csv -out ./tnbc_local_idx.mpileup
     ```
     This command will produce `tnbc_local_idx.mpileup` which is the mpileup that only contains the loci selected by SCIPhi. This is the input to Phylovar.
2. ### Running Phylovar on TNBC data
   After creating the files described in [Filtering the non-informative sites](https://github.com/NakhlehLab/Phylovar/blob/main/README.md#filtering-the-non-informative-sites), we have all the inputs ready for Phylovar. Put the list of cell names (e.g. `cellNames.txt`) and the pre-processed mpileup (e.g. `tnbc_local_idx.mpileup`) into the same folder (like `data` in this repository). Then run the following command to run Phylovar:
   ```
   phylovar.py -names cellNames.txt -out ./data/output/ -indir ./data/ -infile tnbc_local_idx.mpileup -mdthr 1 -mode 0 -niter 100000 -stoch 1 -M 5 -c 10 -verbose yes -errN 100
   ```
   The above command runs Phylovar on `tnbc_local_idx.mpileup` with 10 hill-climbing chains (`-c 10`), on a pool of 5 cores (`-M 5`) each for 100000 iterations (`-niter 100000`) using stochastic hill-climbing (`-stoch 1`). After each 100 iterations, the false-positive and false-negative error rates are sampled (`-errN 100`).
   A zip file containing the outputs of Phylovar is provided in `data` directory, named `output.zip`. We ran Phylovar with Python 2.7, on a pool of five CPU’s, each with 48 cores (AMD EPYC 7642) on a node with 192 GB RAM.
3. ### Description of the outputs
   The output folder contains:
   - `phylovar_stats.txt`: summary of running time, estimated false-negative and false-positive rates, and the log-likelihood value from the best hill-climbing chain.
   - `info.csv`: summary table of Robinson-Foulds distances between the trees found by different hill-climbing chains.
   - `snv.vcf`: the VCF file containing the inferred genotypes from the best hill-climbing chain.
   - `outputs`: a folder containing the thread-specific results including the VCF files from each hill-climbing chain. The names of the subdirectories are in the form of `p<number of the thread>`. In each subdirectory, the corresponding VCF file is named as `snv.vcf`.
   
## Contact
If you have any questions, please contact edrisi@rice.edu or edrisi.rice@gmail.com
