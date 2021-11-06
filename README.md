# Phylovar
Phylovar is a likelihood-based method for joint inference of SNVs and cell lineages from SCS datasets consisting of a large number of loci. It is implemented fully in Python and benefits from the vectorized operations in NumPy to scale up to millions of loci. 
To run Phylovar, you need to install Python (version 2.7>). In the recent versions of Python, Numpy, and Scipy are installed as built-in libraries. However, if you do not have them installed, you can install them using the following commands for either pip or conda installation:
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
## How to install required packages 
### Pip installation
```
pip install ete3
pip install DendroPy
```
### Conda installation
``` conda install -c etetoolkit ete3 ```
``` conda install -c bioconda dendropy```
## Usage
```
python main.py -names cellnames.txt -out ./output/ -indir ./input/ 
-infile mpileup.mpileup -mdthr 1 -niter 100000 -stoch yes -M 1 -c 1 -verbose yes -errN 100
```

## option-definition
* -names : path to the file containing cell names
* -out : path to the output directory
* -mthdr : missing data threshold (default value is 1)
* -fn : false negative rate (default value is 0.1)
* -fp : false positive rate (default value is 1e-8)
