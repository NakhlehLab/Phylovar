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
import parser
import VCF
from VCF import gen_VCF
import scipy.cluster.hierarchy as hier
import scipy.spatial.distance as dist
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

	'''Count the number of character pairs which violate
	 the infinite-sites assumption '''
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
				# if Matrix[cell][p]==0 and Matrix[cell][q]==0:
				# 	B_prime[p][q][0]=1
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
	''' calculates the pmf of binomial distribution '''
	k_decimal = Decimal(k)
	n_decimal = Decimal(n)
	p_decimal = Decimal(p)
	tmp = Decimal(gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1))+Decimal(k_decimal*p_decimal.ln()+(n_decimal-k_decimal)*Decimal(1-p_decimal).ln())
	return tmp.exp()

def compute_HammingDistance(X):
	return (2 * np.inner(X-0.5, 0.5-X) + X.shape[1] / 2)
    # return (X[:, None, :] != X).sum(2)

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

def chunk_ML(splits, inverted_splits, ones, zeros):
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
		# if node.is_root():
			# mutated_names = set([])
			# pass
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
	s = time.time()
	(L_wg, new_genotypes) = chunk_ML(splits=splits_sorted, inverted_splits=splits_inverted, ones=one_Ls, zeros=zero_Ls)
	print("total time of maximizing likelihood in seconds: ", str(time.time()-s))
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
					if np.float128(Decimal((fn_decimal)*A+(1-fn_decimal)*B).ln()) > np.float128(Decimal((1-fp_decimal)*A+(fp_decimal)*B).ln()):
						init_mat[i][j]=1
						init_L+=np.float128(Decimal((fn_decimal)*A+(1-fn_decimal)*B).ln()) 
					else:
						init_L+=np.float128(Decimal((1-fp_decimal)*A+(fp_decimal)*B).ln())
					one_Ls[i][j] = np.float128(Decimal((fn_decimal)*A+(1-fn_decimal)*B).ln())
					zero_Ls[i][j] = np.float128(Decimal((1-fp_decimal)*A+(fp_decimal)*B).ln())
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
	# write_csv(dist_matrix=dist_matrix(mat=matrix_), names_=names, dir_=out_path)
	write_csv(dist_matrix=compute_HammingDistance(X=matrix_), names_=names, dir_=out_path)
	pdm = dnd.PhylogeneticDistanceMatrix.from_csv(src=open(out_path+"pdist_matrix_4dendropy.csv"))
	# # tns = dnd.TaxonNamespace()
	nj_tree = pdm.nj_tree()
	# # nj_tree.print_plot()
	nj_newick = nj_tree.as_string("newick", suppress_edge_lengths=True)
	# print nj_newick
	nj_newick = nj_newick.replace("[&U] ", "")
	ete_nj = Tree(nj_newick, format=8)
	# print ete_nj
	ete_nj.set_outgroup(ete_nj&"normal")
	# print ete_nj
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
	print("current recursion limit: ",str(sys.getrecursionlimit()))
	Limit = 500000
	print("set the recursion limit to: ", str(Limit))
	sys.setrecursionlimit(Limit)
	#########################################################################
	fn_given = 0.1
	fp_given = 1e-8
	missing_data_threshold = 5
	out_path = "./"
	data_path = ""
	cell_names_path = ""
	true_matrix_path = ""
	ap = argparse.ArgumentParser()
	ap.add_argument("-names","--cell names", required=True, help="file containing the cell names")
	ap.add_argument("-tm","--true matrix", required=False, help="tab-separated file containing the true mutations")
	ap.add_argument("-out","--output directory",required=False, help="path to the output directory")
	ap.add_argument("-in","--input mpileup file",required=True, help="path to the input file")
	ap.add_argument("-mdthr","--missing data threshold",required=False, help="minimum coverage for each ref-var pair, default value 10")
	ap.add_argument("-fp","--false positive rate",required=False, help="false positive error rate, default value 1e-08")
	ap.add_argument("-fn","--false negative rate",required=False, help="false negative error rate, default value 0.1")
	ap.add_argument("-vio","--maximum number of violations",required=False, help="maximum number of violations of infinite-sites assumption, default value 0")
	args = vars(ap.parse_args())

	if args['cell names']!=None:
		cell_names_path = args['cell names']
	else:
		print("Please enter the path to the cell names \nUsage: python scVILP_main.py -in <path to the mpileup file> -names <path to the list of cell names>")
		sys.exit()

	if args['output directory']!=None:
		out_path = args['output directory']
	if not out_path.endswith("/"):
		out_path+="/"
	if args['input mpileup file']!=None:
		data_path = args['input mpileup file']
	else:
		print("Please enter the path to the mpileup file\nUsage: python scVILP_main.py -in <path to the mpileup file> -names <path to the list of cell names>")
		sys.exit()
	if args['true matrix']!=None:
		true_matrix_path = args['true matrix']

	if args['missing data threshold']!=None:
		missing_data_threshold = float(args['missing data threshold'])
	if args['false positive rate']!=None:
		fp_given = float(args['false positive rate'])
	if args['false negative rate']!=None:
		fn_given = float(args['false negative rate'])
	if args['maximum number of violations']!=None:
		K_ = int(args['maximum number of violations'])


	##############################################################################
	########################### Parse the mpileup file ###########################
	parse_time = time.time()
	(read_counts, alts, refs, chroms, positions, names, depths) = parser.Parse(cell_names_path, data_path)
	print(time.time()-parse_time, " seconds for parsing")

	n_cells = read_counts.shape[0]
	n_positions = read_counts.shape[1]
	print(n_positions)
	print(n_cells)


	initialization = None

	init_time = time.time()


	zero_Ls = np.zeros((read_counts.shape[0],read_counts.shape[1]))
	one_Ls = np.zeros((read_counts.shape[0],read_counts.shape[1]))


	(initialization, uL, one_Ls, zero_Ls) = ML_initialization(rc=read_counts, names=names, mu0=1e-3, mu1=0.5, fn=fn_given, fp=fp_given, MD=missing_data_threshold, one_Ls=one_Ls, zero_Ls=zero_Ls)
	np.save(out_path+"zero.npy", zero_Ls)
	np.save(out_path+"one.npy", one_Ls)
	np.save(out_path+"init.npy", initialization)
	np.save(out_path+"l.npy", uL)


	# one_Ls = None
	# zero_Ls = None
	# uL = None
	# initialization = None
	# one_Ls = np.load(out_path+"one.npy")
	# zero_Ls = np.load(out_path+"zero.npy")
	# initialization = np.load(out_path+"init.npy")
	# uL = np.load(out_path+"l.npy")

	
	print(time.time()-init_time, " seconds for the initalization")
	tree_time = time.time()
	(ete_nj_init, nj_tree_init) = get_tree(matrx=initialization, nams=names, out_path=out_path)
	print(time.time()-tree_time, " seconds for the NJ tree")
	print(ete_nj_init)


	##### reconstruct the NJ tree given the initial results
	best_mat = None
	best_L = float("-inf")
	
	start = time.time()
	(new_mat, Likelihood) = ML(tree=ete_nj_init, rc=read_counts, names=names, mu0=1e-3, mu1=0.5, fn=fn_given, fp=fp_given, MD=missing_data_threshold, one_Ls=one_Ls, zero_Ls=zero_Ls)
	# print Likelihood
	print(time.time()-start, " seconds for the final ML")

	if Likelihood>best_L:
		best_L = Likelihood

	n_iterations = 2000
	best_Ls = [best_L]
	num_bests = [1]
	stack  = [ete_nj_init]
	top_ids = set()
	top_ids.add(ete_nj_init.get_topology_id())
	for it in range(n_iterations):
		print("------------------ iteration # "+str(it+1)+" ------------------")
		Ls = []
		ts = []
		# ms = []
		for item_ in stack:
			print("NNI")
			tree_list = NNI.Main(in_tree=item_, N=n_cells)
			for tree_ in tree_list:
				if tree_.get_topology_id() not in top_ids:
					top_ids.add(tree_.get_topology_id())
					(mat_, Likelihood) = ML(tree=tree_, rc=read_counts, names=names, mu0=1e-3, mu1=0.5, fn=fn_given, fp=fp_given, MD=missing_data_threshold, one_Ls=one_Ls, zero_Ls=zero_Ls)
					Ls.append(Likelihood)
					ts.append(tree_)
			print("SPR")
			tree_list = SPR.Main(in_tree=item_, N=n_cells, N_dest=n_cells)
			for tree_ in tree_list:
				if tree_.get_topology_id() not in top_ids:
					top_ids.add(tree_.get_topology_id())
					(mat_, Likelihood) = ML(tree=tree_, rc=read_counts, names=names, mu0=1e-3, mu1=0.5, fn=fn_given, fp=fp_given, MD=missing_data_threshold, one_Ls=one_Ls, zero_Ls=zero_Ls)
					Ls.append(Likelihood)
					ts.append(tree_)
		print(len(Ls))
		max_ = max(Ls)
		if max_ > best_L:
			print("found better tree(s)")
			best_L = max_
			best_Ls.append(best_L)
			stack = []
			for i in range(len(ts)):
				if Ls[i]==max_:
					stack.append(ts[i])
			num_bests.append(len(stack))
		elif max_ == best_L:
			print("found equally good tree(s)")
			best_Ls.append(best_L)
			for i in range(len(ts)):
				if Ls[i]==max_:
					stack.append(ts[i])
			num_bests.append(len(stack))
		else:
			print("no more better proposed trees")
			print("terminating the search...")
			break


	print(best_L)
	print(best_Ls)
	print(num_bests)
	print((stack[0]).write(format=9))
	(best_mat, Likelihood) = ML(tree=stack[0], rc=read_counts, names=names, mu0=1e-3, mu1=0.5, fn=fn_given, fp=fp_given, MD=missing_data_threshold, one_Ls=one_Ls, zero_Ls=zero_Ls)
	### write the vcf file
	gen_VCF(out_dir=out_path, genotype_mat=mat_, read_count_mat_=read_counts, chrs=chroms, posits=positions, alt_counts=alts, rfs=refs, ids=names, dps=depths)