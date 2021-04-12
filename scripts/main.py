# =============================================================================
# File name: main.py
# Created By  : Mohammadamin Edrisi
# Created Date: Sun April 11 2021
# Python Version: 2.7
# =============================================================================
"""The module is the main code in the implementation of Phylovar algorithm 
 for SNV detection. It detects the mutations from the candidate loci while 
 contrsucting a phylogeny though a hill climbing search """
# =============================================================================
# imports
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
from cdecimal import Decimal
from ete3 import Tree
import dendropy as dnd
import mpileup_parser
import VCF_writer
from VCF_writer import gen_VCF
import VCF_parser
from VCF_parser import parse_vcf
import NNI
import SPR


def iterable(obj):
    try:
        iter(obj)
    except Exception:
        return False
    else:
        return True

def write_List(lst, outfile):
	with open(outfile, "w") as outf:
		outf.writelines('\t'.join(str(j) for j in i) + '\n' for i in lst)
	outf.close()
	return 
def read_List(infile):
	arr = []
	with open(infile, "r") as inf:
		for line in inf:
			arr.append([])
			splt = line.strip().split("\t")
			for s in splt:
				arr[len(arr)-1].append(int(float(s)))
	inf.close()
	return arr

def check_InCompatibility(Matrix):

	""" Count the number of character pairs which violate
	 the infinite-sites assumption """
	num_incompatibles = 0
	pairs = []
	Matrix = np.array(Matrix)
	n = Matrix.shape[0]
	l = Matrix.shape[1]
	B_prime = [[[0 for k in range(4)] for j in range(l)] for i in range(l)]
	for p in range(l):
		q = p+1
		while q<l:
			count01=0
			count10=0
			count11=0
			for cell in range(n):
				if count01+count10+count11==3:
					break  
				if Matrix[cell][p]==0 and Matrix[cell][q]==1:
					B_prime[p][q][1]=1
					count01=1
				elif Matrix[cell][p]==1 and Matrix[cell][q]==0:
					B_prime[p][q][2]=1
					count10=1
				elif Matrix[cell][p]==1 and Matrix[cell][q]==1:
					B_prime[p][q][3]=1
					count11=1
			q+=1
	for p in range(l):
		q=p+1
		while q<l:
			s = sum(B_prime[p][q])
			if s==3:
				num_incompatibles+=1
				pairs.append((p,q))
			q+=1
	print(pairs)
	return num_incompatibles

def Binomial_pmf(k,n,p):
	""" calculates the pmf of binomial distribution """
	k_decimal = Decimal(k)
	n_decimal = Decimal(n)
	p_decimal = Decimal(p)
	tmp = Decimal(gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1))+Decimal(k_decimal*p_decimal.ln()+(n_decimal-k_decimal)*Decimal(1-p_decimal).ln())
	return tmp.exp()

def compute_HammingDistance(X):
	return (2 * np.inner(X-0.5, 0.5-X) + X.shape[1] / 2)

def dist_matrix(mat):
	pdist_mat = np.zeros((mat.shape[0],mat.shape[0]))
	for p in range(pdist_mat.shape[0]):
		q=p+1
		while q<pdist_mat.shape[0]:
			a = np.array(mat[p])
			b = np.array(mat[q])
			pdist_mat[p][q] = np.count_nonzero(a!=b)
			pdist_mat[q][p] = pdist_mat[p][q]
			q+=1
	return pdist_mat

def write_csv(dist_matrix, names_, dir_):
	outf = open(dir_+"pdist_matrix_4dendropy.csv","w")
	st = ""
	for name in names_:
		st+=(","+name)
	outf.write(st+"\n")	
	for i in range(len(dist_matrix)):
		st=names_[i]
		for val in dist_matrix[i]:
			st+=(","+str(val))
		outf.write(st+"\n")
	outf.close()
	return 

def vectorized_mutation_placement(splits, inverted_splits, ones, zeros):
	tmp = np.matmul(splits,ones)+np.matmul(inverted_splits,zeros)
	split_indexes = np.argmax(tmp,axis=0)
	genotypes = splits[split_indexes]

	L = np.sum(np.max(tmp, axis=0))
	return (L, genotypes)

def ML(tree, rc, names, mu0, mu1, fn, fp, MD, one_Ls, zero_Ls):
	fp_decimal = Decimal(fp)
	fn_decimal = Decimal(fn)
	cells = copy.copy(names)
	cells = set(cells)
	split_list = []
	split_counts = []
	for node in tree.traverse("preorder"):
		### leaf names containing the mutation
		if node.is_leaf():
			mutated_names = set([node.name])
		else:
			mutated_names = set(node.get_leaf_names(is_leaf_fn=None))
		not_mutated_names = cells - mutated_names
		tmp_arr = []
		for name in names:
			if name in list(mutated_names):
				tmp_arr.append(1)
			else:
				tmp_arr.append(0)
		split_list.append(tmp_arr)
		split_counts.append(tmp_arr.count(1))


	mutated_names = set([])
	not_mutated_names = cells - mutated_names
	tmp_arr = []
	for name in names:
		if name in list(mutated_names):
			tmp_arr.append(1)
		else:
			tmp_arr.append(0)
	split_list.append(tmp_arr)
	split_counts.append(tmp_arr.count(1))
	split_list = np.array(split_list)
	indexes = np.argsort(split_counts)
	splits_sorted = split_list[indexes,:]
	splits_inverted = 1 - splits_sorted

	new_genotypes = np.zeros((rc.shape[1], rc.shape[0]))
	L_wg = np.float128(0)
	(L_wg, new_genotypes) = vectorized_mutation_placement(splits=splits_sorted, inverted_splits=splits_inverted, ones=one_Ls, zeros=zero_Ls)
	return (new_genotypes.T,L_wg)

def ML_initialization(rc, names, mu0, mu1, fn, fp, MD, one_Ls, zero_Ls):
	init_L = np.float128(0)
	fp_decimal = Decimal(fp)
	fn_decimal = Decimal(fn)
	n = rc.shape[0]
	l = rc.shape[1]
	init_mat = np.zeros((n,l))
	for i in range(n):
		print("ML initialization, cell # "+str(i+1))
		for j in range(l):
			if rc[i][j][0]+rc[i][j][1]>=MD:
					r = int(rc[i][j][0])
					v = int(rc[i][j][1])
					A = Binomial_pmf(v,r+v,mu0)
					B = Binomial_pmf(v,r+v,mu1)
					one_L = np.float128(Decimal((fn_decimal)*A+(1-fn_decimal)*B).ln())
					zero_L = np.float128(Decimal((1-fp_decimal)*A+(fp_decimal)*B).ln())
					if one_L > zero_L:
						init_mat[i][j]=1
						init_L+=one_L
					else:
						init_L+=zero_L
					one_Ls[i][j] = one_L
					zero_Ls[i][j] = zero_L
			else:
				pass
	return (init_mat, init_L, one_Ls, zero_Ls)

def sort_mat(mat_unsorted):
	tmp_array = copy.copy(mat_unsorted)
	R1 = tmp_array.T
	distMatrix = dist.pdist(R1)
	distSquareMatrix = dist.squareform(distMatrix)
	linkageMatrix = hier.linkage(distMatrix,method='ward')
	dendro = hier.dendrogram(linkageMatrix)
	leaves1 = dendro['leaves']
	transformedData = R1[leaves1,:]
	## leaves1 for the positions

	R2 = tmp_array
	distMatrix = dist.pdist(R2)
	distSquareMatrix = dist.squareform(distMatrix)
	linkageMatrix = hier.linkage(distMatrix,method='ward')
	dendro = hier.dendrogram(linkageMatrix)
	leaves2 = dendro['leaves']
	transformedData = transformedData[:,leaves2]

	return (leaves1, leaves2)

def get_tree(matrx, nams, out_path):
	names = copy.copy(nams)
	names.append("normal")
	matrx = np.array(matrx)
	matrix_ = copy.copy(matrx)
	matrix_ = np.vstack((matrix_,[0 for i in range(matrix_.shape[1])]))
	write_csv(dist_matrix=compute_HammingDistance(X=matrix_), names_=names, dir_=out_path)
	pdm = dnd.PhylogeneticDistanceMatrix.from_csv(src=open(out_path+"pdist_matrix_4dendropy.csv"))
	nj_tree = pdm.nj_tree()
	nj_newick = nj_tree.as_string("newick", suppress_edge_lengths=True)
	nj_newick = nj_newick.replace("[&U] ", "")
	ete_nj = Tree(nj_newick, format=8)
	ete_nj.set_outgroup(ete_nj&"normal")
	(ete_nj&"normal").detach()
	names.remove("normal")
	newroot = ete_nj.get_tree_root().get_children()[0]
	ete_nj = newroot.detach()
	matrix_ = np.delete(matrix_,-1,axis=0)
	return (ete_nj, nj_tree)

def tree_move(tr, k, move):
	new_tr = None
	if move==0:
		tr_list = NNI.Main(in_tree=tr, N=k)
	else:
		return tr
	return tr_list

if __name__=="__main__":
	# =============================================================================
	# parse the arguments
	# =============================================================================
	ap = argparse.ArgumentParser()
	ap.add_argument("-names", "--names", required=True, help="file containing the cell names in the same directory as the output directory")
	ap.add_argument("-tm", "--tm", required=False, help="tab-separated file containing the true mutations (for simulation purposes)")
	ap.add_argument("-out","--out",required=False, help="path to the output directory")
	ap.add_argument("-stats", "--stats", required=False, help="name of the output file reporting estimated and input values including best tree topology, given/inferred error rates, etc.", 
		default="./phylovar_stats")
	ap.add_argument("-vcf", "--vcf", required=False, help="name of the output VCF file containing the SNV calls", default="./snv.vcf")
	ap.add_argument("-intree", "--intree", required=False, help="name of the initial topology file for starting tree search (in newick format in the output directory). If not given, a NJ tree is created as initial tree")
	ap.add_argument("-infile", "--infile",required=True, help="name of the input mpileup file in the same directory as the output directory")
	ap.add_argument("-mdthr", "--mdthr",required=False, help="minimum coverage for each ref-var pair, default value 5", default=1, type=int)
	ap.add_argument("-fp", "--fp",required=False, help="false positive error rate, default value 1e-08", default=0, type=float)
	ap.add_argument("-fn", "--fn",required=False, help="false negative error rate, default value 0.1", default=0, type=float)
	ap.add_argument("-opt", "--opt", required=False, help="optimization algorithm: there are two options, 0 is the steepest descent, and 1 is hill climbing search which is recommended for large datasets more than around 30 single cells", 
		choices = [0,1], type=int)
	ap.add_argument("-mode", "--mode", required=False, 
		help="boolean argument indicating whether the code is used for simulation study (1) or real data analysis (0)", 
		default=0, type=int, choices=[0,1])
	ap.add_argument("-pre", "--pre", required=False, 
		help="boolean argument that is 1 when the user wants to utilize the output matrices from the previous runs, they are expected to be in the same directory designated as output directory", 
		default=0, type=int, choices=[0,1])
	ap.add_argument("-niter","--niter", required=False, help="number of iterations", type=int, default=100)
	args = ap.parse_args()
	print(args)

	if not args.out.endswith("/"):
		args.out+="/"
	if args.fn!=None and args.mode==1:
		args.fn = args.fn*(random.uniform(0.5,2))

	# =============================================================================
	# parse the mpileup file from the input 
	# =============================================================================
	print("parsing the mpileup file ...")
	parse_time = time.time()
	(read_counts, alts, refs, chroms, positions, names, depths) = mpileup_parser.Parse(args.out+args.names, args.out+args.infile)
	print("parsing the mpileup file was done in %s seconds" % (time.time() - parse_time))
	n_cells = read_counts.shape[0]
	n_positions = read_counts.shape[1]
	print("number of cells %s" %n_cells)
	print("number of candidate loci %s" %n_positions)

	# =============================================================================
	# calculating or reading the initial likelihood matrices
	# =============================================================================
	init_time = time.time()
	if args.pre==1:
		print("reading the matrices from existing files ...")
		one_Ls = np.load(args.out+"one.npy")
		zero_Ls = np.load(args.out+"zero.npy")
		initialization = np.load(args.out+"init.npy")
		uL = np.load(args.out+"l.npy")
		print("numpy matrices are loaded")
	if args.pre==0:
		print("initialization ...")
		zero_Ls = np.zeros((read_counts.shape[0],read_counts.shape[1]))
		one_Ls = np.zeros((read_counts.shape[0],read_counts.shape[1]))
		(initialization, uL, one_Ls, zero_Ls) = ML_initialization(rc=read_counts, names=names, mu0=1e-3, mu1=0.5, fn=args.fn, fp=args.fp, MD=args.mdthr, one_Ls=one_Ls, zero_Ls=zero_Ls)
		np.save(args.out+"zero.npy", zero_Ls)
		np.save(args.out+"one.npy", one_Ls)
		np.save(args.out+"init.npy", initialization)
		np.save(args.out+"l.npy", uL)
		print("initialization was done in %s seconds" %(time.time() - init_time))

	print("dimension of the inital matrix: ", initialization.shape)
	print("dimension of the zeros matrix: ", zero_Ls.shape)
	print("dimension of the ones matrix: ", one_Ls.shape)
	print("dimension of the likelihood matrix: ", uL.shape)
	# =============================================================================
	# constructing or reading the initial tree 
	# =============================================================================
	if args.intree==None:
		print("constructing the NJ tree ...")
		tree_time = time.time()
		(ete_nj_init, nj_tree_init) = get_tree(matrx=initialization, nams=names, out_path=args.out)
		print("NJ tree was computed in %s seconds" %(time.time()-tree_time))
	else:
		print("reading the initial topology from an existing newick file ...")
		ete_nj_init = Tree(args.out+args.intree, format=8)

	# =============================================================================
	# calculate the best likelihood (mutation placement) for the initial tree
	# =============================================================================

	(new_mat, Likelihood) = ML(tree=ete_nj_init, rc=read_counts, names=names, mu0=1e-3, mu1=0.5, fn=args.fn, fp=args.fp, MD=args.mdthr, one_Ls=one_Ls, zero_Ls=zero_Ls)
	# the best likelihood value so far
	best_L = Likelihood
	# save the best tree in one variable for hill climbing search
	best_tree = ete_nj_init
	# number of iterations
	n_iterations = args.niter
	# total number of trees processed during the search
	n_samples = 0
	# list of best likelihood values out of each iteration
	best_Ls = [best_L]
	# bag of best topologies to perform tree rearrangement moves on
	bag  = [ete_nj_init]
	# set of topology ids processed so far
	top_ids = set()
	# add the initial tree to the set of topology ids
	top_ids.add(ete_nj_init.get_topology_id())
	# save the len of bag for further investigation 
	num_bests = [1]

	# =============================================================================
	# search for topology
	# =============================================================================
	print("searching for better topologies ...")
	for it in range(n_iterations):
		print(" iteration # %s" %(it+1))
		
		Ls = []
		ts = []
		found = False

		for item_ in bag:
			# =============================================================================
			# performing NNI on the tree 
			# =============================================================================
			print("performing NNI and SPR ...")
			# list of the trees returned from NNI and SPr moves on the selected topology from bag of trees 
			proposed_trees = NNI.Main(in_tree=item_, N=n_cells)
			proposed_trees.extend(SPR.Main(in_tree=item_, N=n_cells, N_dest=n_cells))
			# shuffle the list of trees for random sampling 
			random.shuffle(proposed_trees)
			print("# of proposed trees by NNI and SPR %s" %(len(proposed_trees)))
			for tree_ in proposed_trees:
				# only the trees not in the topology set are analyzed
				if tree_.get_topology_id() not in top_ids:
					n_samples+=1
					# add the new tree to the topology set
					top_ids.add(tree_.get_topology_id())
					(mat_, Likelihood) = ML(tree=tree_, rc=read_counts, names=names, mu0=1e-3, mu1=0.5, fn=args.fn, fp=args.fp, MD=args.mdthr, one_Ls=one_Ls, zero_Ls=zero_Ls)
					# when doing hill climbing search, once a better topology is observed, save it and move to the next iteration
					if args.opt==1:
						if Likelihood>best_L:
							# found a better topology using hill climbing search
							best_L = Likelihood
							best_Ls.append(Likelihood)
							# do not change bag while being inside the inner loop
							# save the best tree in another variable
							best_tree = tree_
							found = True
							break
						else:
							continue
					# if the optimization mode is steepest descent, save all the observations into a list for further analysis
					else:
						if Likelihood>=best_L:
							print("found better or equally good tree by steepest descent")
							Ls.append(Likelihood)
							ts.append(tree_)
							found = True

		# out of the inner loop, check for each optimization mode
		# if a better tree is found and opt is hill climbing, replace the old tree with the new one and go to next iteration with new bag
		if args.opt==1 and found:
			print("found a better topology by hill climbing search")
			bag = [best_tree]
			num_bests.append(len(bag))
			continue
		# if no better was found and opt is hill climbing, terminate the search
		elif args.opt==1 and not found:
			print("hill climbing did not find a better tree, terminating the search ...")
			bag = [best_tree]
			num_bests.append(len(bag))
			break
		# opt is steepest descent and better tree was found, update the bag 
		elif args.opt==0 and found:
			max_ = max(Ls)
			if max_ > best_L:
				print("found better tree(s)")
				best_L = max_
				best_Ls.append(best_L)
				bag = []
				for i in range(len(ts)):
					if Ls[i]==max_:
						bag.append(ts[i])
			else:
				print("found equally good tree(s)")
				best_Ls.append(best_L)
				for i in range(len(ts)):
					if Ls[i]==max_:
						bag.append(ts[i])
			num_bests.append(len(bag))
		# opt is steepest descent, and no better tree was found, so terminate the search
		else:
			num_bests.append(len(bag))
			print("no better tree was found by steepest descent, terminating the search ...")

	# =============================================================================
	# calculate the accuracy of inference for when having the ground truth 
	# =============================================================================
	precision_ = 0
	recall_ = 0
	F1_ = 0
	if args.mode==1 and args.tm!=None:
		# parse the ground truth genotype matrix
		acc_time = time.time()
		true_array = []
		true_positions = []
		with open(args.out+args.tm,"r") as tfile:
			for line in tfile:
				if "#position" in line:
					continue
				else:
					true_array.append([])
					tmp = line.strip().split("\t")
					true_positions.append(int(tmp[0]))
					for mut in tmp[1:]:
						if int(mut)!=0 and int(mut)!=4:
							true_array[len(true_array)-1].append(1)
						else:
							true_array[len(true_array)-1].append(0)
		true_array=np.array(true_array)
		true_array=true_array.T
		True_list = copy.copy(true_array)
		# calculate the number of false calls and accuracy measurements 
		false_positive_=0
		false_negative_=0
		true_positive_=0
		a = set(positions)
		b = set(true_positions)
		# find the overlap between the candidate loci in the ground truth and the inferred matrix

		intersect_pos = list(a.intersection(b))
		# find the indices of intersection loci in the inferred matrix
		indx_arr = [positions.index(intersect_pos[j]) for j in range(len(intersect_pos))]
		# extract the submatrix containing the calls for the intersection loci
		best_submat = best_mat[:,indx_arr]
		# find the indices of the intersection loci in the ground truth matrix 
		indx_arr_true = [true_positions.index(intersect_pos[j]) for j in range(len(intersect_pos))]
		# extact the submatrix containing the calls for the intersection loci
		true_submat = True_list[:,indx_arr_true]
		# calculate the true positive, false positive, and false negative calls in the intersection 
		true_positive_ = sum(sum(best_submat + true_submat==2))
		false_negative_ = sum(sum(true_submat - best_submat==1))
		false_positive_ = sum(sum(best_submat - true_submat==1))


		# extract the inferred submatrix containing the different loci
		f = [i for i in range(n_positions)]
		residual_indices = list(set(f)-set(indx_arr))
		residual_inference = best_mat[:,residual_indices]
		# the 1's in the residual submatrix are the falsely detected mutations
		false_positive_+= residual_inference.sum()

		# extract the ground truth submatrix contaning the different loci 
		t = [i for i in range(len(true_positions))]
		residual_indices_t = list(set(t)-set(indx_arr_true))
		residual_true = True_list[:,residual_indices_t]
		# the 1's in the ground truth residual submatrix are the mutations that are not detected by the method
		false_negative_+= residual_true.sum()

		# calculate and print the accuracy measurements
		precision_ = float(true_positive_)/float(true_positive_+false_positive_)
		recall_ =  float(true_positive_)/float(true_positive_+false_negative_)
		F1_ = 2*float(precision_*recall_)/float(precision_+recall_)
		print("# true positives %s, # false negatives %s, and false positives %s" % (true_positive_, false_negative_, false_positive_))
		print("accuracy meansurements: precision: %s, recall: %s, F1 score: %s " % (precision_, recall_, F1_))
		print("accuracy measurements were calculated in %s seconds" %(time.time()-acc_time))
	elif args.mode==1 and args.tm==None:
		print("please provide the name of the ground truth genotype matrix")

	# =============================================================================
	# write the info into output files
	# =============================================================================
	# calculate the best mutation placement for the best topology
	(best_mat, Likelihood) = ML(tree=bag[0], rc=read_counts, names=names, mu0=1e-3, mu1=0.5, fn=args.fn, fp=args.fp, MD=args.mdthr, one_Ls=one_Ls, zero_Ls=zero_Ls)
	run_time = time.time()-init_time
	print("total runtime was %s" %run_time)
	# write the vcf file 
	gen_VCF(out_=args.out+args.vcf, genotype_mat=mat_, read_count_mat_=read_counts, chrs=chroms, posits=positions, alt_counts=alts, rfs=refs, ids=names, dps=depths)
	# write the best tree topology/topologies into a newick format file
	if args.opt==0:
		for i in range(len(bag)):
			bag[i].write(format=8, outfile=args.out+"tree"+str(i+1)+"nw")
	# write the remaining information including the actuall runtime, list of best likelihood values and so on
	o_stats = open(args.out+args.stats, "w")
	o_stats.write("runtime: "+str(run_time)+" seconds"+"\n")
	o_stats.write("fn: "+str(args.fn)+"\n")
	o_stats.write("fp: "+str(args.fp)+"\n")
	o_stats.write("total number of sampled trees: "+str(n_samples)+"\n")
	o_stats.write("likelihoods:\n")
	o_stats.write(','.join([str(x) for x in best_Ls])+"\n")
	o_stats.write(','.join([str(y) for y in num_bests])+"\n")
	if args.mode==1 and args.tm!=None:
		o_stats.write("precision: "+str(precision_)+"\n")
		o_stats.write("recall: "+str(recall_)+"\n")
		o_stats.write("F1 score: "+str(F1_))
	o_stats.close()


	






