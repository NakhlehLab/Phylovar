# =============================================================================
# File name: main.py
# Created By  : Mohammadamin Edrisi
# Created Date: Sun April 11 2021
# =============================================================================
"""The module is the main code in the implementation of Phylovar algorithm 
 for SNV detection. It detects the mutations from the candidate loci while 
 contrsucting a phylogeny though a hill climbing search """
# =============================================================================

# ========================================== #
# import the required packages and libraries #
# ========================================== #
import numpy as np
import argparse
import copy
import sys
import os
import csv
import random
import math
import scipy
from scipy.stats import *
from scipy.special import comb
from scipy.special import gammaln
import time
from ete3 import Tree
import dendropy as dnd
import mpileup_parser
import VCF_writer
from VCF_writer import gen_VCF
import VCF_parser
from VCF_parser import parse_vcf
import NNI
import SPR
import SNL
import matplotlib
import multiprocessing as mp
import uuid
import pandas as pd

DATA_TYPE = np.float64

def ThereAreDuplicates(listOfElems):
    ''' Check if given list contains any duplicates '''
    if len(listOfElems) == len(set(listOfElems)):
        return False
    else:
        return True

def gen_unique_ids(num):
    ids = []
    while True:
        ids = []
        for i in range(num):
            ids.append(str(uuid.uuid4()))
        if not ThereAreDuplicates(ids):
            break
    return ids

def unique_rows(a):
    a = np.ascontiguousarray(a)
    unique_a = np.unique(a.view([('', a.dtype)]*a.shape[1]))
    return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))
#================================================================== #
# Binomial_pmf fucntion calculates the pmf of binomial distribution #
#================================================================== #

def Binomial_pmf_(k, n, p):
    tmp = (gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)) + (
        k * np.log(p, dtype=DATA_TYPE) + (n - k) * np.log(1 - p, dtype=DATA_TYPE))
    return np.exp(tmp, dtype=DATA_TYPE)+np.nextafter(0,1)

#===================================================================================== #
# compute_HammingDistance is used for calculating pairwise distance matrix for NJ tree #
#===================================================================================== #
def compute_HammingDistance(X):
    ### calculate the pairwise distance matrix
    ### the missing data do not contribute in the distances
    ### the distance is defined as the number of (1,0) or (0,1) pairs between the two vectors
    pdist = np.zeros((X.shape[0], X.shape[0]))
    for p in range(X.shape[0]):
        q = p + 1
        while q < X.shape[0]:
            x = len(np.argwhere(((X[p] - X[q] == 1) | (X[q] - X[p] == 1)) & (X[p] != 0.5) & (X[q] != 0.5)))
            pdist[p][q] = x
            pdist[q][p] = x
            q += 1
    return pdist

#============================================================== #
# write_csv writes the pairwise distance matrix into a txt file #
#============================================================== #
def write_csv(dist_matrix, names_, dir_):
    outf = open(dir_ + "pdist_matrix.csv", "w")
    st = ""
    for name in names_:
        st += ("," + name)
    outf.write(st + "\n")
    for i in range(len(dist_matrix)):
        st = names_[i]
        for val in dist_matrix[i]:
            st += ("," + str(val))
        outf.write(st + "\n")
    outf.close()
    return

#================================================================================================================== #
# vectorized_mutation_placement computes the maximum likelihood given a topology and L arrays in vectorized fashion #
#================================================================================================================== #
def vectorized_mutation_placement(splits, inverted_splits, ones, zeros):
    s_time = time.time()
    tmp = np.matmul(splits, ones) + np.matmul(inverted_splits, zeros)
    # print("matrix multiplication took %s seconds" %(time.time()-s_time))
    s_time = time.time()
    split_indexes = np.argmax(tmp, axis=0)
    # print("calling argmax took %s seconds" %(time.time()-s_time))
    genotypes = splits[split_indexes]
    s_time = time.time()
    L = np.sum(np.max(tmp, axis=0))
    # print("summing over the genotype likelihoods took %s seconds" %(time.time()-s_time))
    return (L, genotypes)

#================================================================================================== #
# vectorized_mutation_placement computes the maximum likelihood and returns the corresponding split #
#================================================================================================== #
def vectorized_mutation_placement_with_splits(splits, inverted_splits, ones, zeros, d_):
    tmp = np.matmul(splits, ones) + np.matmul(inverted_splits, zeros)
    split_indexes = np.argmax(tmp, axis=0)
    selected_splits = splits[split_indexes]
    corr_names = [d_[str(selected_splits[x])] for x in range(selected_splits.shape[0])]

    L = np.sum(np.max(tmp, axis=0))
    return (L, selected_splits, corr_names)

#=========================================================================================== #
# ML is the main function that computes the maximum likelihood given a topology and L arrays #
#=========================================================================================== #
def ML(tree, names, mu0, mu1, MD, one_Ls_, zero_Ls_, data_type):
    s_time = time.time()
    cells = copy.copy(names)
    leaf_dict = dict((c, i) for i, c in enumerate(cells))
    split_list = []

    for node in tree.traverse():
        leaves = node.get_leaf_names(is_leaf_fn=None)
        z = np.zeros(len(cells))
        for l in leaves:
            z[leaf_dict[l]] = 1
        split_list.append(np.array(z, np.bool))

    split_list.append(np.zeros(len(cells), np.bool))
    split_list = np.array(split_list)
    split_counts = [np.count_nonzero(split_list[i]) for i in range(len(split_list))]
    indexes = np.argsort(split_counts)
    splits_sorted = split_list[indexes, :]
    splits_inverted = np.invert(splits_sorted)
    # print("computing the split matrix took %s seconds" %(time.time()-s_time))
    s_time = time.time()
    new_genotypes = np.zeros((one_Ls_.shape[1], one_Ls_.shape[0]), np.bool)
    L_wg = data_type(0)
    # print("preparing the genotypes matrix took %s seconds" %(time.time()-s_time))
    s_time = time.time()
    (L_wg, new_genotypes) = vectorized_mutation_placement(splits=splits_sorted, inverted_splits=splits_inverted,ones=one_Ls_, zeros=zero_Ls_)
    # print("calling vectorized mutation placement took %s seconds" %(time.time()-s_time))
    return (new_genotypes.T, L_wg)

#================================================================ #
# best_split is another version of ML that returns the best split #
#================================================================ #
def best_split(tree, names, mu0, mu1, MD, one_Ls_, zero_Ls_, data_type):
    if tree.get_tree_root().name == "":
        tree.get_tree_root().name = "artificial_root"
    cells = copy.copy(names)
    leaf_dict = dict((c, i) for i, c in enumerate(cells))
    split_list = []
    d = {}
    n = []
    for node in tree.traverse():
        leaves = node.get_leaf_names(is_leaf_fn=None)
        z = np.zeros(len(cells))
        for l in leaves:
            z[leaf_dict[l]] = 1
        split_list.append(np.array(z, np.bool))
        n.append(node.name)

    n.append("*")
    split_list.append(np.zeros(len(cells), np.bool))
    for i in range(len(split_list)):
        d[str(split_list[i])] = n[i]
    split_list = np.array(split_list)
    split_counts = [np.count_nonzero(split_list[i]) for i in range(len(split_list))]
    indexes = np.argsort(split_counts)
    splits_sorted = split_list[indexes, :]
    splits_inverted = np.invert(splits_sorted)

    new_genotypes = np.zeros((one_Ls_.shape[1], one_Ls_.shape[0]), np.bool)
    L_wg = data_type(0)
    (L_wg, new_genotypes, corr_names_) = vectorized_mutation_placement_with_splits(splits=splits_sorted,
        inverted_splits=splits_inverted,
        ones=one_Ls_, zeros=zero_Ls_, d_=d)
    return (new_genotypes.T, L_wg, corr_names_)


#===================================================================== #
# sort_mat performs hierarchical clustering on a given genotype matrix #
#===================================================================== #
def sort_mat(mat_unsorted):
    tmp_array = copy.copy(mat_unsorted)
    R1 = tmp_array.T
    distMatrix = dist.pdist(R1)
    distSquareMatrix = dist.squareform(distMatrix)
    linkageMatrix = hier.linkage(distMatrix, method='ward')
    dendro = hier.dendrogram(linkageMatrix)
    leaves1 = dendro['leaves']
    transformedData = R1[leaves1, :]
    ## leaves1 for the positions

    R2 = tmp_array
    distMatrix = dist.pdist(R2)
    distSquareMatrix = dist.squareform(distMatrix)
    linkageMatrix = hier.linkage(distMatrix, method='ward')
    dendro = hier.dendrogram(linkageMatrix)
    leaves2 = dendro['leaves']
    transformedData = transformedData[:, leaves2]

    return (leaves1, leaves2)

#============================================================================ #
# get_tree returns the neighbor-joining tree given a pairwise distance matrix #
#============================================================================ #
def get_tree(matrx, nams, out_path):
    names = copy.copy(nams)
    names.append("normal")
    matrx = np.array(matrx)
    matrix_ = copy.copy(matrx)
    matrix_ = np.vstack((matrix_, [0 for i in range(matrix_.shape[1])]))
    write_csv(dist_matrix=compute_HammingDistance(X=matrix_), names_=names, dir_=out_path)
    # write_csv(dist_matrix=compute_L1Distance(X=matrix_), names_=names, dir_=out_path)
    pdm = dnd.PhylogeneticDistanceMatrix.from_csv(src=open(out_path + "pdist_matrix.csv"))
    nj_tree = pdm.nj_tree()
    nj_newick = nj_tree.as_string("newick", suppress_edge_lengths=True)
    nj_newick = nj_newick.replace("[&U] ", "")
    ete_nj = Tree(nj_newick, format=8)
    ete_nj.set_outgroup(ete_nj & "normal")
    (ete_nj & "normal").detach()
    names.remove("normal")
    newroot = ete_nj.get_tree_root().get_children()[0]
    ete_nj = newroot.detach()
    matrix_ = np.delete(matrix_, -1, axis=0)
    ete_nj.get_tree_root().name = "artificial_root"
    new_names = gen_unique_ids(num=len(ete_nj)-2)
    counter = 0
    for node in ete_nj.traverse():
        if (not node.is_leaf()) and (not node.is_root()):
            node.name = new_names[counter]
            counter+=1
    return (ete_nj, nj_tree)

#================================================================================================ #
# Search performs the search for the best topology given a false negative and positive error rate #
#================================================================================================ #
def Search_sequential(arguments):


    runtimes = []
    [args, names, fns, fps, thread_num, read_counts, B_mat, A_mat] = arguments
    sample_params = None
    sample_params = True
    # sample_rate = 100
    sample_rate = args.errN

    a = read_counts[:, :, 0]
    b = read_counts[:, :, 1]
    # find the indices of non-missing entries
    idx = np.argwhere(a + b >= args.mdthr)


    random.seed(thread_num+1)
    # load the numpy array of initial matrix
    p_ = args.out + "outputs/"+"p"+str(thread_num)+"/"
    if not os.path.isdir(p_):
        os.makedirs(p_)
    s_time = time.time()
    zeros = []
    ones = []
    for fp in fps:
        zero_path = args.out + "zeros/"+"fp"+format(fp)+".npy"
        zero_mat = np.load(zero_path)
        zero_mat = zero_mat.astype(DATA_TYPE)
        zeros.append(zero_mat)

    for fn in fns:
        one_path = args.out + "ones/"+"fn"+format(fn)+".npy"
        one_mat = np.load(one_path)
        one_mat = one_mat.astype(DATA_TYPE)
        ones.append(one_mat)

    initialization = np.load(args.out + "init_noparams.npy")
    # initialize the matrices of zeros and ones
    zero_path = args.out + "zero_noparams.npy"
    one_path = args.out + "one_noparams.npy"
    current_one_mat = np.load(one_path)
    current_zero_mat = np.load(zero_path)
    current_one_mat = one_mat.astype(DATA_TYPE)
    current_zero_mat = zero_mat.astype(DATA_TYPE)
    ete_nj_init = None
    verbose = args.verbose


    # ======================================== #
    # constructing or reading the initial tree #
    # ======================================== #
    if args.intree == None:
        print("constructing the NJ tree in thread "+str(thread_num))
        tree_time = time.time()
        (ete_nj_init, nj_tree_init) = get_tree(matrx=initialization, nams=names,
                                               out_path=p_)
        print("NJ tree in thread "+str(thread_num)+ " was computed in %s seconds" % (
                    time.time() - tree_time))
        ete_nj_init.write(format=8, outfile=p_+ "init_tree.nw")

    else:
        print("reading the initial topology from an existing newick file ...")
        ete_nj_init = Tree(args.intree, format=8)
        ete_nj_init.get_tree_root().name = "artificial_root"

    # ======================================================================= #
    # calculate the best likelihood (mutation placement) for the initial tree #
    # ======================================================================= #
    st = time.time()
    (new_mat, Likelihood) = ML(tree=ete_nj_init, names=names, mu0=1e-3, mu1=0.5, MD=args.mdthr, one_Ls_=current_one_mat,zero_Ls_=current_zero_mat, data_type=DATA_TYPE)

    print("computing the best likelihood took %s seconds" %(time.time()-st))
    # the best likelihood value so far
    best_L = Likelihood
    # the current likelihood for tree search
    current_L = Likelihood
    # save the best tree in one variable for hill climbing search
    best_tree = ete_nj_init
    # save the best alpha in one variable for hill climbing search
    best_alpha = fps[0]
    # save the best beta in one variable for hill climbing search
    best_beta = fns[0]
    # save the current tree for topology search
    current_tree = ete_nj_init
    # save the current false-positive rate for the search
    current_alpha = fps[0]
    # save the current false-negative rate for the search
    current_beta = fns[0]
    # number of iterations
    n_iterations = args.niter
    # total number of trees processed during the search
    n_samples = 1
    # list of best likelihood values out of each iteration
    Ls = [Likelihood]
    # list of best alpha values out of each iteration
    alphas = [fps[0]]
    # list of best beta values out of each iteration
    betas = [fns[0]]
    # counter for iterations 
    iter_counter = 0
    # waiting time for terminating the search
    waiting_iter = args.w
    # counter for waiting time
    waiting_counter = 0
    # whether the search is stochastic or not
    stochastic = args.stoch
    # list of likelihoods for the steepest descent search
    tmp_ls = []
    # list of trees for the steepest descent search
    tmp_trs = []
    # list of best scores
    best_Ls = [Likelihood]
    # list of best false positive samples
    best_alphas = [current_alpha]
    # list of best false negative samples
    best_betas = [current_beta]



    # =================== #
    # search for topology #
    # =================== #
    print("searching thread # "+str(thread_num))
    while (waiting_counter<waiting_iter and stochastic==0) or (iter_counter<n_iterations and stochastic==1):

        ### temporary
        iteration_time = time.time()


        found = False
        accepted = False
        if verbose and (iter_counter % 100 == 0):
            print("iteration # %s in thread %s, best likelihood %s best alpha %s best beta %s" % (iter_counter+1, thread_num, best_L, best_alpha, best_beta))

        t = current_tree
        t_one = current_one_mat
        t_zero = current_zero_mat
        t_alpha = current_alpha
        t_beta = current_beta


        if (iter_counter%sample_rate!=0 and sample_params) or (not sample_params) or (iter_counter==0 and sample_params):
            # ========================== #
            # performing NNI on the tree #
            # ========================== #
            # list of the trees returned from NNI and SPR moves on the selected topology from bag of trees
            # choose a move randomly
            moves = [0, 1, 2]
            move = random.choice(moves)
            if move == 0:
                # print("NNI")
                a = NNI.Main(in_tree=copy.copy(current_tree), N=1)
                t = random.choice(a)
            elif move == 1:
                # print("SNL")
                t = SNL.Main(in_tree=copy.copy(current_tree), N=1, sample=True, all_=False)[0]
            else:
                # print("SPR")
                a = SPR.Main(in_tree=copy.copy(current_tree), N=1, N_dest=1)
                while len(a) == 0:
                    # print("empty list of topology by SPR")
                    a = SPR.Main(in_tree=copy.copy(current_tree), N=1, N_dest=1)
                t = a[0]

        elif (iter_counter%sample_rate==0 and sample_params):

            (mat_, Likelihood_) = ML(tree=t, names=names, mu0=1e-3, mu1=0.5, MD=args.mdthr, one_Ls_=t_one,
               zero_Ls_=t_zero, data_type=DATA_TYPE)
            fn_idx = np.argwhere((initialization == 0) & (mat_ == 1))
            non_missing_idx = np.argwhere((mat_ == 1) & (initialization != 0.5))
            t_beta = float(len(fn_idx)) / float(len(non_missing_idx)+0.0000000001)

            fp_idx = np.argwhere((initialization == 1) & (mat_ == 0))
            non_missing_idx = np.argwhere((mat_ == 0) & (initialization != 0.5))
            t_alpha = float(len(fp_idx)) / (float(len(non_missing_idx)+0.0000000001))

            t_one = np.zeros((read_counts.shape[0], read_counts.shape[1]), DATA_TYPE)
            t_zero = np.zeros((read_counts.shape[0], read_counts.shape[1]), DATA_TYPE)

            # s_time = time.time()
            tmp_mat = (np.float64(1) - np.float64(t_alpha)) * A_mat + np.float64(t_alpha) * B_mat
            t_zero[idx[:,0],idx[:,1]] = np.log(tmp_mat[idx[:,0],idx[:,1]], dtype=DATA_TYPE)
            tmp_mat = np.float64(t_beta) * A_mat + (np.float64(1) - np.float64(t_beta)) * B_mat
            t_one[idx[:,0],idx[:,1]] = np.log(tmp_mat[idx[:,0],idx[:,1]], dtype=DATA_TYPE)
            
            # print("runtime of numpy calculation %s" %(time.time()-s_time))


        n_samples += 1
        (mat_, Likelihood) = ML(tree=t, names=names, mu0=1e-3, mu1=0.5, MD=args.mdthr, one_Ls_=t_one,
               zero_Ls_=t_zero, data_type=DATA_TYPE)

        mat_ = None
        # when doing hill climbing search, once a better topology is observed, save it and move to the next iteration
        if Likelihood >= current_L:
            # found a better topology using hill climbing search
            if Likelihood > best_L:

                best_L = Likelihood
                # update the parameters depending on the iteration count:

                if (iter_counter%sample_rate!=0 and sample_params) or (not sample_params) or (iter_counter==0 and sample_params):
                    best_tree = copy.copy(t)

                elif (iter_counter%sample_rate==0 and sample_params):
                    best_alpha = t_alpha
                    best_beta = t_beta


                found = True
            # do not change bag while being inside the inner loop
            # save the best tree in another variable
            # update the current parameters depending on iteration count:

            if (iter_counter%sample_rate!=0 and sample_params) or (not sample_params) or (iter_counter==0 and sample_params):
                current_tree = copy.copy(t)

            elif (iter_counter%sample_rate==0 and sample_params):
                current_alpha = t_alpha
                current_beta = t_beta
                current_zero_mat = t_zero
                current_one_mat = t_one


            current_L = Likelihood

            Ls.append(Likelihood)
            best_Ls.append(best_L)

            alphas.append(current_alpha)
            best_alphas.append(best_alpha)

            betas.append(current_beta)
            best_betas.append(best_beta)

            accepted = True
        else:
            p_accept = None
            ## add stochasticity
            if stochastic == 1:
                p_accept = math.exp(Likelihood - current_L)
            else:
                p_accept = 0
            if random.random() < p_accept:

                current_L = Likelihood
                # update the current parameters depending on the iteration count:
                if (iter_counter%sample_rate!=0 and sample_params) or (not sample_params) or (iter_counter==0 and sample_params):
                    current_tree = copy.copy(t)

                elif (iter_counter%sample_rate==0 and sample_params):
                    current_alpha = t_alpha
                    current_zero_mat = t_zero
                    current_beta = t_beta
                    current_one_mat = t_one


                accepted = True

            Ls.append(current_L)
            best_Ls.append(best_L)

            alphas.append(current_alpha)
            best_alphas.append(best_alpha)

            betas.append(current_beta)
            best_betas.append(best_beta)

        # out of the inner loop, check for each optimization mode
        # if a better tree is found and opt is hill climbing, replace the old tree with the new one and go to next iteration with new bag
        if found:
            # found a better topology by hill climbing search
            waiting_counter = 0
            iter_counter += 1
            continue
        # if no better was found and opt is hill climbing, terminate the search
        if not found:
            # hill climbing did not find a better tree, terminating the search
            waiting_counter += 1
            iter_counter += 1

        ### temporary
        runtimes.append(time.time()-iteration_time)


    # =========================================================================== #
    # save the list of likelihoods and plot them for the current fn and fp values #
    # =========================================================================== #

    best_Ls = np.array(best_Ls, dtype=DATA_TYPE)
    Ls = np.array(Ls, dtype=DATA_TYPE)
    best_alphas = np.array(best_alphas)
    best_betas = np.array(best_betas)
    alphas = np.array(alphas)
    betas = np.array(betas)

    tmp_mat = (np.float64(1) - np.float64(best_alpha)) * A_mat + np.float64(best_alpha) * B_mat
    t_zero[idx[:,0],idx[:,1]] = np.log(tmp_mat[idx[:,0],idx[:,1]], dtype=DATA_TYPE)
    tmp_mat = np.float64(best_beta) * A_mat + (np.float64(1) - np.float64(best_beta)) * B_mat
    t_one[idx[:,0],idx[:,1]] = np.log(tmp_mat[idx[:,0],idx[:,1]], dtype=DATA_TYPE)
    p = args.out + "zeros/"+"fp"+format(best_alpha)+".npy"
    np.save(p, t_zero)
    p = args.out + "ones/"+"fn"+format(best_beta)+".npy"
    np.save(p, t_one)


    if stochastic:

        np.save(p_ + "best_Ls.npy", best_Ls)
        np.save(p_ + "Ls.npy", Ls)
        np.save(p_ + "best_alphas.npy", best_alphas)
        np.save(p_ + "alphas.npy", alphas)
        np.save(p_ + "best_betas.npy", best_betas)
        np.save(p_ + "betas.npy", betas)
    else:
        np.save(p_ + "Ls.npy", Ls)
        np.save(p_ + "alphas.npy", alphas)
        np.save(p_ + "betas.npy", betas)

    best_tree.write(format=8, outfile=p_ + "best_tree.nw")
    print("processing thread %s was done in %s seconds: best beta %s, best alpha %s, best likelihood %s" %(thread_num, time.time()-s_time, best_beta, best_alpha, best_L))
    print("average runtime per iteration %s" %(float(sum(runtimes))/float(len(runtimes))))
    return (best_L, best_tree, n_samples, iter_counter)


if __name__ == "__main__":
    # =================== #
    # parse the arguments #
    # =================== #
    ap = argparse.ArgumentParser()
    ap.add_argument("-names", "--names", required=True,
                    help="file containing the cell names in the same directory as the output directory")
    ap.add_argument("-tm", "--tm", required=False,
                    help="tab-separated file containing the true mutations (for simulation purposes)")
    ap.add_argument("-out", "--out", required=False, help="path to the output directory")
    ap.add_argument("-indir", "--indir", required=False, help="path to the input directory")
    ap.add_argument("-stats", "--stats", required=False,
                    help="name of the output file reporting estimated and input values including best tree topology, given/inferred error rates, etc.",
                    default="phylovar_stats.txt")
    ap.add_argument("-vcf", "--vcf", required=False, help="name of the output VCF file containing the SNV calls",
                    default="snv.vcf")
    ap.add_argument("-intree", "--intree", required=False,
                    help="path to the initial tree for topology search (in newick format). If not given, a NJ tree is created as initial tree")
    ap.add_argument("-infile", "--infile", required=True,
                    help="name of the input mpileup file in the same directory as the output directory")
    ap.add_argument("-mdthr", "--mdthr", required=False, help="minimum coverage for each ref-var pair, default value 1",
                    default=1, type=int)
    ap.add_argument("-mode", "--mode", required=False,
                    help="boolean argument indicating whether the code is used for simulation study (1) or real data analysis (0)",
                    default=0, type=int, choices=[0, 1])
    ap.add_argument("-niter", "--niter", required=False, help="number of iterations", type=int, default=100)
    ap.add_argument("-stoch", "--stoch", required=False, help="whether the search is stochastic or not", type=lambda x: (str(x).lower() in ['true','1', 'yes']),
                    default=True)
    ap.add_argument("-w", "--w", required=False, help="number of iterations to wait before termination", type=int,
                    default=1000)
    ap.add_argument("-M", "--M", required=False, help="number of threads to use", type=int,
                    default=1)
    ap.add_argument("-c", "--c", required=False, help="number of chains for tree search", type=int,
                    default=1)
    ap.add_argument('-verbose','--verbose', default=False, type=lambda x: (str(x).lower() in ['true','1', 'yes']), help="print the iterations")
    ap.add_argument('-errN', '--errN', default=100, type=int, help="number of iterations per which to sample new error rates")
    args = ap.parse_args()
    print(args)

    if not args.out.endswith("/"):
        args.out += "/"
    if not args.indir.endswith("/"):
        args.indir += "/"
    if not os.path.isdir(args.out):
        os.makedirs(args.out)

    # ===================================== #
    # parse the mpileup file from the input #
    # ===================================== #
    parse_time = time.time()
    (read_counts, alts, refs, chroms, positions, names, depths, tags) = mpileup_parser.Parse(
        cell_names_file=args.indir + args.names, mpileup_file=args.indir + args.infile)
    cov_ = np.sum(read_counts, axis=2)
    md_count = np.argwhere(cov_ == 0)
    md_per_total = float(len(md_count)) / float(read_counts.shape[0] * read_counts.shape[1])
    print("percentage of missing data in the dataset %s" % md_per_total)
    print("parsing the mpileup file was done in %s seconds" % (time.time() - parse_time))
    n_cells = read_counts.shape[0]
    n_positions = read_counts.shape[1]
    print("number of cells %s" % n_cells)
    print("number of candidate loci %s" % n_positions)

    # ====================================================== #
    # calculating or reading the initial likelihood matrices #
    # ====================================================== #
    init_time = time.time()
    uL = None
    fns_ = None
    fps_ = None
    best_mat = None
    Likelihood = None
    bag = None
    n_samples_ = 0
    fn_best = 0
    iters = 0
    one_mat = None
    zero_mat = None
    best_L = None
    fn_est = None
    fp_est = None
    node_names = None
    one_Ls_array = None
    zero_Ls_array = None
    zero_Ls_noparams = None
    one_Ls_noparams = None
    L_dict = None
    a_dict = None
    b_dict = None
    best_index = None
    iters = None
    max_ = None
    B_mat = None
    A_mat = None

    print("estimating the false error rates by NJ tree")
    n = read_counts.shape[0]
    l = read_counts.shape[1]
    # all the init matrices with non-zero false positive and negative rates
    init_noparams = np.full((n, l), 0.5, dtype=DATA_TYPE)
    # init_noparams = np.zeros((n, l), dtype=np.float64)

    a = read_counts[:, :, 0]
    b = read_counts[:, :, 1]
    # find the indices of non-missing entries
    idx = np.argwhere(a + b >= args.mdthr)
    print(idx.shape)
    mu0 = 1e-3
    mu1 = 0.5
    zero_Ls_noparams = np.zeros((read_counts.shape[0], read_counts.shape[1]), DATA_TYPE)
    one_Ls_noparams = np.zeros((read_counts.shape[0], read_counts.shape[1]), DATA_TYPE)
    B_mat = np.zeros((read_counts.shape[0], read_counts.shape[1]), DATA_TYPE)
    A_mat = np.zeros((read_counts.shape[0], read_counts.shape[1]), DATA_TYPE)
    print("calculating likelihood scores for non-missing entries ...")
    c_time = time.time()
    r_mat = read_counts[idx[:,0],idx[:,1],0]
    v_mat = read_counts[idx[:,0],idx[:,1],1]
    n_mat = r_mat + v_mat
    tmp_mu0 = gammaln(n_mat + 1) - gammaln(v_mat + 1) - gammaln(n_mat - v_mat + 1) + (v_mat * np.log(mu0, dtype=DATA_TYPE) + (n_mat - v_mat) * np.log(1 - mu0, dtype=DATA_TYPE))
    tmp_mu1 = gammaln(n_mat + 1) - gammaln(v_mat + 1) - gammaln(n_mat - v_mat + 1) + (v_mat * np.log(mu1, dtype=DATA_TYPE) + (n_mat - v_mat) * np.log(1 - mu1, dtype=DATA_TYPE))
    A_mat[idx[:,0],idx[:,1]] = np.exp(tmp_mu0, dtype=DATA_TYPE)+np.nextafter(0,1)
    B_mat[idx[:,0],idx[:,1]] = np.exp(tmp_mu1, dtype=DATA_TYPE)+np.nextafter(0,1)


    zero_Ls_noparams[idx[:,0],idx[:,1]] = np.log(A_mat[idx[:,0],idx[:,1]], dtype=DATA_TYPE)
    one_Ls_noparams[idx[:,0],idx[:,1]] = np.log(B_mat[idx[:,0],idx[:,1]], dtype=DATA_TYPE)

    print("the calculation was done in %s seconds" % (time.time() - c_time))
    print("matrix calculations ...")
    ones_idx_noparams = np.argwhere(one_Ls_noparams > zero_Ls_noparams)
    zero_idx_noparams = np.argwhere(one_Ls_noparams < zero_Ls_noparams)
    init_noparams[ones_idx_noparams[:, 0], ones_idx_noparams[:, 1]] = 1.0
    init_noparams[zero_idx_noparams[:, 0], zero_idx_noparams[:, 1]] = 0.0

    np.save(args.out + "init_noparams.npy", init_noparams)
    np.save(args.out + "zero_noparams.npy", zero_Ls_noparams)
    np.save(args.out + "one_noparams.npy", one_Ls_noparams)

    if not os.path.isdir(args.out + "outputs/"):
        os.makedirs(args.out + "outputs/")
    p_ = args.out + "outputs/est/"
    if not os.path.isdir(p_):
        os.makedirs(p_)

    (ete_nj_init, nj_tree_init) = get_tree(matrx=init_noparams, nams=names,
                                           out_path=p_)
    ete_nj_init.write(format=8, outfile=p_+ "init_tree.nw")

    (new_mat, Likelihood) = ML(tree=ete_nj_init, names=names, mu0=1e-3, mu1=0.5, MD=args.mdthr, one_Ls_=one_Ls_noparams,zero_Ls_=zero_Ls_noparams, data_type=DATA_TYPE)

    fn_idx = np.argwhere((init_noparams == 0) & (new_mat == 1))
    non_missing_idx = np.argwhere((new_mat == 1) & (init_noparams != 0.5))
    fn_est = float(len(fn_idx)) / float(len(non_missing_idx)+0.0000000001)
    print("estimated false-negative rate %s" % fn_est)

    fp_idx = np.argwhere((init_noparams == 1) & (new_mat == 0))
    non_missing_idx = np.argwhere((new_mat == 0) & (init_noparams != 0.5))
    fp_est = float(len(fp_idx)) / (float(len(non_missing_idx)+0.0000000001))
    print("estimated false-positive rate %s" % fp_est)
    fns_ = [fn_est]
    fps_ =[fp_est]

    one_Ls_array = [np.zeros((read_counts.shape[0], read_counts.shape[1]), DATA_TYPE) for y in fns_]
    zero_Ls_array = [np.zeros((read_counts.shape[0], read_counts.shape[1]), DATA_TYPE) for u in fps_]

    for h in range(len(fps_)):
        tmp_mat = (np.float64(1) - np.float64(fps_[h])) * A_mat + np.float64(fps_[h]) * B_mat
        zero_Ls_array[h][idx[:,0],idx[:,1]] = np.log(tmp_mat[idx[:,0],idx[:,1]],dtype=DATA_TYPE)
    for i_ in range(len(fns_)):
        tmp_mat = np.float64(fns_[i_]) * A_mat + (np.float64(1) - np.float64(fns_[i_])) * B_mat
        one_Ls_array[i_][idx[:,0],idx[:,1]] = np.log(tmp_mat[idx[:,0],idx[:,1]],dtype=DATA_TYPE)


    if not os.path.isdir(args.out + "zeros/"):
        os.makedirs(args.out + "zeros/")
    if not os.path.isdir(args.out + "ones/"):
        os.makedirs(args.out + "ones/")
    for i in range(len(fps_)):
        p = "zeros/"+"fp"+format(fps_[i])+".npy"
        np.save(args.out + p, zero_Ls_array[i])
    for j in range(len(fns_)):
        p = "ones/"+"fn"+format(fns_[j])+".npy"
        np.save(args.out + p, one_Ls_array[j])

    if not os.path.isdir(args.out + "outputs/"):
        os.makedirs(args.out + "outputs/")

    pool_args = []
    for t_ in range(args.c):
        # ======================================== # 
        # make arguments for each of the processes #
        # ======================================== #
        pool_args.append([args, names, fns_, fps_, t_+1, read_counts, B_mat, A_mat])
    pool = mp.Pool(args.M)
    pool.map(Search_sequential, pool_args)
    pool.close()
    pool.join()

    # =================================================== #
    # collect the lists of likelihoods from all processes #
    # =================================================== #
    L_dict = {}
    a_dict = {}
    b_dict = {}
    tree_dict = {}

    if not args.stoch:
        for t_ in range(args.c):
            p = args.out + "outputs/"+"p"+str(t_+1)+"/"
            L_dict[t_] = np.load(p+"Ls.npy")
            a_dict[t_] = np.load(p+"alphas.npy")
            b_dict[t_] = np.load(p+"betas.npy")
            tree_dict[t_] = Tree(p+"best_tree.nw", format=8)
    else:
        for t_ in range(args.c):
            p = args.out + "outputs/"+"p"+str(t_+1)+"/"
            L_dict[t_] = np.load(p+"best_Ls.npy")
            a_dict[t_] = np.load(p+"best_alphas.npy")
            b_dict[t_] = np.load(p+"best_betas.npy")
            tree_dict[t_] = Tree(p+"best_tree.nw", format=8)

    # ========================================================= #
    # print the best likelihood and the corresponding fn and fp #
    # ========================================================= #
    max_ = L_dict[0][-1]
    best_index = (a_dict[0][-1],b_dict[0][-1],0)
    for t_ in range(args.c):
        if L_dict[t_][-1]>max_:
            best_index = (a_dict[t_][-1],b_dict[t_][-1],t_)
            max_ = L_dict[t_][-1]
    print("the best likelihood is %s for false-negative rate %s and false-positive rate %s" %(max_, best_index[1], best_index[0]))
    print("table of reports from threads\n-----------------------------------------------------------")
    for t_ in range(args.c):
        print("best likelihood from thread %s: %s, %s, %s, %s" %(t_+1,L_dict[t_][-1],a_dict[t_][-1],b_dict[t_][-1], tree_dict[t_].get_topology_id()))
    print("-----------------------------------------------------------")
    iters = len(L_dict[best_index[2]])

    # ================================================= #
    # build a dataframe to store the pairwise distances #
    # ================================================= #
    cols = ["RF_to_p"+str(i+1) for i in range(args.c)]+["log_likelihood", "topology_id"]
    rows = ["p"+str(i+1) for i in range(args.c)]
    chains_info = [[None for i in range(len(cols))] for j in range(args.c)]
    for q1 in range(args.c):
        chains_info[q1][q1] = 0
        chains_info[q1][-2] = L_dict[q1][-1]
        chains_info[q1][-1] = tree_dict[q1].get_topology_id()
        for q2 in range(args.c):
            if q2>q1:
                rf = (tree_dict[q1]).robinson_foulds(tree_dict[q2])
                chains_info[q1][q2] = float(rf[0])/float(rf[1])
            elif q2<q1:
                chains_info[q1][q2] = chains_info[q2][q1]
            else:
                pass
    df = pd.DataFrame(chains_info, columns=cols, index=rows)
    print("pairwise RF distances between the trees\n========================================================")
    print(df)
    df.to_csv(args.out+"info.csv")
    # ========================================================== #
    # read the best tree corresponding the best fn and fp values #
    # ========================================================== #
    best_path = args.out + "outputs/"+ "p"+str(best_index[2]+1)+"/"
    
    best_path_one = args.out +"ones/"+ "fn" + format(best_index[1])+".npy"
    best_path_zero = args.out +"zeros/" + "fp" + format(best_index[0])+".npy"
    selected_tree = Tree(best_path + "/best_tree.nw", format=8)
    init_noparams = np.load(args.out + "init_noparams.npy")
    init_noparams = init_noparams.astype(DATA_TYPE)
    one_path = best_path_one
    zero_path = best_path_zero
    one_mat = np.load(one_path)
    zero_mat = np.load(zero_path)

    # ================================================================ #
    # reproduce the genotype matrix corresponding to the other results #
    # ================================================================ #
    print("Writing the VCF files of all chains")
    for chn in range(args.c):
        if chn!=best_index[2]:

            chn_path = chn_path = args.out + "outputs/"+ "p"+str(chn+1)+"/"
            
            chn_path_one = args.out +"ones/"+ "fn" + format(b_dict[chn][-1])+".npy"
            chn_path_zero = args.out +"zeros/" + "fp" + format(a_dict[chn][-1])+".npy"
            one_path = chn_path_one
            zero_path = chn_path_zero
            one_mat = np.load(one_path)
            zero_mat = np.load(zero_path)

            (tree_dict[chn]).get_tree_root().name = "artificial_root"
            (chn_mat, Likelihood, node_names) = best_split(tree=tree_dict[chn], names=names, mu0=1e-3, mu1=0.5, MD=args.mdthr,
                   one_Ls_=one_mat, zero_Ls_=zero_mat, data_type=DATA_TYPE)
            gen_VCF(out_=chn_path + args.vcf, genotype_mat=chn_mat, read_count_mat_=read_counts, chrs=chroms,
                    posits=positions, alt_counts=alts, rfs=refs, ids=names, dps=depths, n_=node_names, tags=tags)
    print("Writing the VCF files is done!")

    # ============================================================== #
    # reproduce the genotype matrix corresponding to the best result #
    # ============================================================== #
    selected_tree.get_tree_root().name = "artificial_root"
    (best_mat, Likelihood, node_names) = best_split(tree=selected_tree, names=names, mu0=1e-3, mu1=0.5, MD=args.mdthr,
           one_Ls_=one_mat, zero_Ls_=zero_mat, data_type=DATA_TYPE)


    # =============================================================================
    # calculate the accuracy of inference for when having the ground truth
    # =============================================================================
    precision_ = 0
    recall_ = 0
    F1_ = 0
    if args.mode == 1 and args.tm != None:
        # parse the ground truth genotype matrix
        acc_time = time.time()
        true_array = []
        true_positions = []
        with open(args.indir + args.tm, "r") as tfile:
            for line in tfile:
                if "#position" in line:
                    continue
                else:
                    true_array.append([])
                    tmp = line.strip().split("\t")
                    true_positions.append(int(tmp[0]))
                    for mut in tmp[1:]:
                        if int(mut) != 0 and int(mut) != 4:
                            true_array[len(true_array) - 1].append(1)
                        else:
                            true_array[len(true_array) - 1].append(0)
        true_array = np.array(true_array)
        true_array = true_array.T
        True_list = copy.copy(true_array)
        # calculate the number of false calls and accuracy measurements
        false_positive_ = 0
        false_negative_ = 0
        true_positive_ = 0
        a = set(positions)
        b = set(true_positions)
        # find the overlap between the candidate loci in the ground truth and the inferred matrix

        intersect_pos = list(a.intersection(b))
        # find the indices of intersection loci in the inferred matrix
        indx_arr = [positions.index(intersect_pos[j]) for j in range(len(intersect_pos))]
        best_submat = best_mat[:, indx_arr]
        # find the indices of the intersection loci in the ground truth matrix
        indx_arr_true = [true_positions.index(intersect_pos[j]) for j in range(len(intersect_pos))]
        # extact the submatrix containing the calls for the intersection loci
        true_submat = True_list[:, indx_arr_true]
        # calculate the true positive, false positive, and false negative calls in the intersection
        true_positive_ = sum(sum(best_submat + true_submat == 2))
        false_negative_ = sum(sum(true_submat - best_submat == 1))
        false_positive_ = sum(sum(best_submat - true_submat == 1))

        # extract the inferred submatrix containing the different loci
        f = [i for i in range(n_positions)]
        residual_indices = list(set(f) - set(indx_arr))
        residual_inference = best_mat[:, residual_indices]
        # the 1's in the residual submatrix are the falsely detected mutations
        false_positive_ += residual_inference.sum()

        # extract the ground truth submatrix contaning the different loci
        t = [i for i in range(len(true_positions))]
        residual_indices_t = list(set(t) - set(indx_arr_true))
        residual_true = True_list[:, residual_indices_t]
        # the 1's in the ground truth residual submatrix are the mutations that are not detected by the method
        false_negative_ += residual_true.sum()

        # calculate and print the accuracy measurements
        precision_ = float(true_positive_) / float(true_positive_ + false_positive_)
        recall_ = float(true_positive_) / float(true_positive_ + false_negative_)
        F1_ = 2 * float(precision_ * recall_) / float(precision_ + recall_)
        print("# true positives %s, # false negatives %s, and false positives %s" % (
        true_positive_, false_negative_, false_positive_))
        print("accuracy meansurements: precision: %s, recall: %s, F1 score: %s " % (precision_, recall_, F1_))
        print("accuracy measurements were calculated in %s seconds" % (time.time() - acc_time))
    elif args.mode == 1 and args.tm == None:
        print("please provide the name of the ground truth genotype matrix")

    # ================================ #
    # write the info into output files #
    # ================================ #
    # calculate the best mutation placement for the best topology/topologies
    run_time = time.time() - init_time
    print("total runtime was %s" % run_time)
    # write the best tree topology/topologies into a newick format file
    # save the vcf files corresponding to the best topologies
    gen_VCF(out_=args.out + args.vcf, genotype_mat=best_mat, read_count_mat_=read_counts, chrs=chroms,
    	posits=positions, alt_counts=alts, rfs=refs, ids=names, dps=depths, n_=node_names, tags=tags)
    # write the remaining information including the actuall runtime, list of best likelihood values and so on
    o_stats = open(args.out + args.stats, "w")
    o_stats.write("runtime: " + str(run_time) + " seconds" + "\n")
    o_stats.write("fn: " + format(best_index[1]) + "\n")
    o_stats.write("fp: " + format(best_index[0]) + "\n")
    o_stats.write("number of iterations: " + str(iters) + "\n")
    o_stats.write("best likelihood: %s \n" % max_)
    o_stats.write("estimated false-negative rate: %s \n" % best_index[1])
    o_stats.write("estimated false-positive rate: %s \n" % best_index[0])
    o_stats.write("missing data percentage: %s \n" % (md_per_total * 100))

    if args.mode == 1 and args.tm != None:
        o_stats.write("precision: " + str(precision_) + "\n")
        o_stats.write("recall: " + str(recall_) + "\n")
        o_stats.write("F1 score: " + str(F1_))
    o_stats.close()
