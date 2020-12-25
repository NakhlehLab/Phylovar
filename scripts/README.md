# Description of the codes
## Neighbor-joining and tree search
Neighbor-joining algorithm is used for initializing the tree topology. Given a topology, the best mutation placement is found under infinite-sites assumption (ISA). Next, using the tree rearrangement techniques, all the possible new topologies out of the currently best topologies are generated. If there are new topologies with better likelihood score than the previous one(s), they are accepted and search continues to the next iteration.  
`NJ_and_search.py` performs this search and terminates when no new topology has a better score. 
## Nearest Neighbor Interchange (NNI)
One of the techniques to search for new topologies given a topology. NNI swaps one of the subtrees under a certain internal node with the sister subtree of that internal node. `NNI.py` performs NNI on a given tree and returns all the possible topologies. NNI emulates a local search in the space of tree topologies. 
## Subtree Pruning and Regrafting (SPR)
SPR refers to a technique where a subtree is pruned from the main tree and is reattached/regrafted to another internal branch of the main tree. `SPR.py` performs this technique on all the subtrees and returns all the possible new topologies. 
