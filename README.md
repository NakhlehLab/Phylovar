# Phylovar 
## How to install required packages 
Phylovar requires ete3 and Dendropy python libraries
## Usage
```
python main.py -names cellnames.txt -out ./output/ -indir ./input/ 
-infile mpileup.mpileup -mdthr 1 -niter 100000 -stoch yes -M 1 -c 1 -verbose yes -errN 100
```

## option-definition
* -in : path to the mpileup file
* -names : path to the file containing cell names
* -out : path to the output directory
* -mthdr : missing data threshold (default value is 1)
* -fn : false negative rate (default value is 0.1)
* -fp : false positive rate (default value is 1e-8)
