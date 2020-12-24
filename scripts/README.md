# Description of the codes
## Neighbor-joining and tree search
Neighbor-joining algorithm is used for initializing the tree topology. Given a topology, the best mutation placement is found under infinite-sites assumption (ISA). Next, using the tree rearrangement techniques, all the possible new topologies out of the currently best topologies are generated. If there are new topologies with better likelihood score than the previous one(s), they are accepted and search continues to the next iteration. 
`NJ_and_search.py` performs this search.
## Nearest Neighbor Interchange (NNI)
## Subtree Pruning and Regrafting (SPR)
