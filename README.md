# PhyloVar 
## Usage
```
python NJ_and_search.py -in loci.mpileup -names cellNames.txt -out ./results -mthdr 1 -fn 0.2 
```

## option-definition
* -in : path to the mpileup file
* -names : path to the file containing cell names
* -out : path to the output directory
* -mthdr : missing data threshold (default value is 1)
* -fn : false negative rate (default value is 0.1)
* -fp : false positive rate (default value is 1e-8)
