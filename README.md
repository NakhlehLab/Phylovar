# Phylovar
Phylovar is a likelihood-based method for joint inference of SNVs and cell lineages from SCS datasets consisting of a large number of loci. It is implemented fully in Python and benefits from the vectorized operations in NumPy to scale up to millions of loci. 
[alg.pdf](https://github.com/mae6/Phylovar/files/8043612/alg.pdf)

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
This will show the description of the arguments.
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

## Contact
If you have any questions, please contact edrisi@rice.edu
