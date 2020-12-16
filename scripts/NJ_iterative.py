import numpy as np 
import argparse
# from gurobipy import *
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
from sklearn.cluster import KMeans
import dendropy as dnd
from dendropy.calculate import treecompare
import parser
import VCF
from VCF import gen_VCF
import parse_VCF
from parse_VCF import parse_vcf
import multiprocessing
import threading 
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import seaborn as sns
import scipy.cluster.hierarchy as hier
import scipy.spatial.distance as dist
import NNI
import SPR

### global variables 
# new_genotypes = None
# L_wg = 0

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

def one_split(L_arr, split_index, split_arr, site_rc, mu0, mu1, fn, fp, MD):
	fp_decimal = Decimal(fp)
	fn_decimal = Decimal(fn)
	tmp_L = 0
	for c in range(len(site_rc)):
		if site_rc[c][0]+site_rc[c][1]>=MD:
			r = int(site_rc[c][0])
			v = int(site_rc[c][1])
			A = Binomial_pmf(v,r+v,mu0)
			B = Binomial_pmf(v,r+v,mu1)
			if split_arr[split_index][c]==1:
				tmp_L += np.float128(Decimal((fn_decimal)*A+(1-fn_decimal)*B).ln())
			else:
				tmp_L += np.float128(Decimal((1-fp_decimal)*A+(fp_decimal)*B).ln())
			# tmp_L += (split_arr[i][c])*np.float128(Decimal((fn_decimal)*A+(1-fn_decimal)*B).ln()) + (1-split_arr[i][c])*np.float128(Decimal((1-fp_decimal)*A+(fp_decimal)*B).ln())
		else:
			pass
	L_arr[split_index] = tmp_L


def one_site_ML_fast(split_arr, site_rc, mu0, mu1, fn, fp, MD):
	L_array = multiprocessing.Array('d', len(split_arr))
	p_arr = []
	max_L = float('-inf')
	num_mutateds = len(site_rc)+1
	genotypes = None
	for splt_index in range(len(split_arr)):
		p_arr.append(multiprocessing.Process(target=one_split, args=(L_array, splt_index, split_arr, site_rc, mu0, mu1, fn, fp, MD)))
	for p in p_arr:
		p.start()
	for p_ in p_arr:
		p_.join()
	Ls = L_array[:]
	for i in range(len(split_arr)):
		tmp_L = Ls[i]
		if (tmp_L==max_L and num_mutateds>=split_arr[i].count(1)) or tmp_L>max_L:
			genotypes = copy.copy(split_arr[i])
			max_L = tmp_L
			num_mutateds = split_arr[i].count(1)

	return (genotypes, max_L)

def one_site_dot(split_arr, ones, zeros):
	split_arr_ = copy.copy(split_arr)
	split_arr_ = np.array(split_arr_)
	inverted_split_arr_ = 1 - split_arr_
	max_L = float('-inf')
	num_mutateds = len(ones)+1
	genotypes = None
	for splt_index in range(len(split_arr_)):
		L_ones = np.matmul(split_arr_[splt_index],ones)
		L_zeros = np.matmul(inverted_split_arr_[splt_index],zeros)
		tmp_L = L_ones + L_zeros
		if (tmp_L==max_L and num_mutateds>=list(split_arr_[splt_index]).count(1)) or tmp_L>max_L:
			genotypes = copy.copy(split_arr_[splt_index])			
			max_L = tmp_L
			num_mutateds = list(split_arr_[splt_index]).count(1)
	return (genotypes, max_L)

def one_site_mat(splits, ones, zeros, counts):
	splits_ = copy.copy(splits)
	splits_ = np.array(splits_)
	inverted_splits_ = 1 - splits_
	Ls = np.matmul(splits_,ones) + np.matmul(inverted_splits_,zeros)
	L_max = max(Ls)
	indices = [i for i, j in enumerate(Ls) if j == L_max]
	num_mutateds = counts[indices[0]]
	selected_index = indices[0]
	for index in indices:
		if counts[index]<=num_mutateds:
			selected_index = index
	genotypes = copy.copy(splits[selected_index])
	return (genotypes, L_max)

def one_site_ML(tree, site_rc, cell_names, mu0, mu1, fn, fp, MD):
	fp_decimal = Decimal(fp)
	fn_decimal = Decimal(fn)
	# print site_rc.shape
	cells = copy.copy(cell_names)
	cells = set(cells)
	genotypes = [0 for i in range(len(cell_names))]
	max_L = float('-inf')
	num_mutateds = len(cell_names)+1
	# L_arr = []
	for node in tree.iter_descendants("preorder"):
		### leaf names containing the mutation
		# if node.is_root():
			# mutated_names = set([])
			# pass
		if node.is_leaf():
			mutated_names = set([node.name])
		else:
			mutated_names = set(node.get_leaf_names(is_leaf_fn=None))
		not_mutated_names = cells - mutated_names
		tmp_L = 0
		for mutated_name in list(mutated_names):
			idx = cell_names.index(mutated_name)
			if site_rc[idx][0]+site_rc[idx][1]>=MD:
				r = int(site_rc[idx][0])
				v = int(site_rc[idx][1])
				A = Binomial_pmf(v,r+v,mu0)
				B = Binomial_pmf(v,r+v,mu1)
				tmp_L += np.float128(Decimal((fn_decimal)*A+(1-fn_decimal)*B).ln())
			else:
				### missing data do not contribute in the likelihood
				pass
		for not_mutated_name in list(not_mutated_names):
			idx = cell_names.index(not_mutated_name)
			if site_rc[idx][0]+site_rc[idx][1]>=MD:
				r = int(site_rc[idx][0])
				v = int(site_rc[idx][1])
				A = Binomial_pmf(v,r+v,mu0)
				B = Binomial_pmf(v,r+v,mu1)
				tmp_L += np.float128(Decimal((1-fp_decimal)*A+(fp_decimal)*B).ln())
			else:
				### missing data do not contribute in the likelihood
				pass
		# L_arr.append(tmp_L)
		if (tmp_L==max_L and num_mutateds>=len(list(mutated_names))) or tmp_L>max_L:
		# if tmp_L>=max_L:
			max_L = tmp_L
			num_mutateds = len(list(mutated_names))
			genotypes = [0 for i in range(len(cell_names))]
			for i in range(len(cell_names)):
				if cell_names[i] in list(mutated_names):
					genotypes[i]=1
				else:
					genotypes[i]=0
	### check the case where there is no mutation
	mutated_names = set([])
	not_mutated_names = cells - mutated_names
	tmp_L = 0
	for mutated_name in list(mutated_names):
		idx = cell_names.index(mutated_name)
		if site_rc[idx][0]+site_rc[idx][1]>=MD:
			r = int(site_rc[idx][0])
			v = int(site_rc[idx][1])
			A = Binomial_pmf(v,r+v,mu0)
			B = Binomial_pmf(v,r+v,mu1)
			tmp_L += np.float128(Decimal((fn_decimal)*A+(1-fn_decimal)*B).ln())
		else:
			### missing data do not contribute in the likelihood
			pass
	for not_mutated_name in list(not_mutated_names):
		idx = cell_names.index(not_mutated_name)
		if site_rc[idx][0]+site_rc[idx][1]>=MD:
			r = int(site_rc[idx][0])
			v = int(site_rc[idx][1])
			A = Binomial_pmf(v,r+v,mu0)
			B = Binomial_pmf(v,r+v,mu1)
			tmp_L += np.float128(Decimal((1-fp_decimal)*A+(fp_decimal)*B).ln())
		else:
			### missing data do not contribute in the likelihood
			pass
	if (tmp_L==max_L and num_mutateds>=len(list(mutated_names))) or tmp_L>max_L:
	# if tmp_L>=max_L:
		max_L = tmp_L
		num_mutateds = len(list(mutated_names))
		genotypes = [0 for i in range(len(cell_names))]
		for i in range(len(cell_names)):
			if cell_names[i] in list(mutated_names):
				genotypes[i]=1
			else:
				genotypes[i]=0

	return (genotypes, max_L)

def chunk_ML(splits, inverted_splits, ones, zeros):
	# global L_wg
	# global new_genotypes
	tmp = np.matmul(splits,ones)+np.matmul(inverted_splits,zeros)
	split_indexes = np.argmax(tmp,axis=0)
	genotypes = splits[split_indexes]
	# (a, b) = genotypes.shape

	# lock.acquire()
	L = np.sum(np.max(tmp, axis=0))
	# for i in range(a):
		# for j in range(b):
	# new_genotypes[start_index:end_index,:] = genotypes
	# lock.release()
	# print L_wg
	# return_dict[index] = (L, genotypes)
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

	# for i in range(2):

	# shared_Array = [multiprocessing.Array('d', rc.shape[0])] * rc.shape[1]
	# shared_Likelihood = multiprocessing.Value('d', 0)
	# p1 = multiprocessing.Process(target=chunk_ML, args=(shared_Likelihood, shared_Array , 0, splits_sorted, splits_inverted, one_Ls[:,0:20], zero_Ls[:,0:20]))
	# p1.start()
	# p2 = multiprocessing.Process(target=chunk_ML, args=(shared_Likelihood, shared_Array, 20, splits_sorted, splits_inverted, one_Ls[:,20:], zero_Ls[:,20:]))
	# p2.start()
	# p1.join()
	# p2.join()

	# new_genotypes = np.zeros((rc.shape[1],rc.shape[0]))
	# for i in range(rc.shape[1]):
	# 	for j in range(rc.shape[0]):
	# 		new_genotypes[i][j] = shared_Array[i][j]
	# L_wg = shared_Likelihood.value

	# Ls_ = np.matmul(splits_sorted,one_Ls)+np.matmul(1 - splits_sorted,zero_Ls)
	# split_indexes = np.argmax(Ls_,axis=0)
	# L_wg = np.sum(np.max(Ls_, axis=0))
	# new_genotypes = splits_sorted[split_indexes]

	# global new_genotypes
	# global L_wg
	new_genotypes = np.zeros((rc.shape[1], rc.shape[0]))
	L_wg = np.float128(0)
	s = time.time()
	(L_wg, new_genotypes) = chunk_ML(splits=splits_sorted, inverted_splits=splits_inverted, ones=one_Ls, zeros=zero_Ls)
	print("total time of maximizing likelihood in seconds: ", str(time.time()-s))
	# N_jobs = 1
	# lock = threading.Lock()
	# manager = multiprocessing.Manager()
	# return_dict = manager.dict()
	# chunk_ML(start_index=0, splits=splits_sorted, inverted_splits=splits_inverted, ones=one_Ls, zeros=zero_Ls)
	# jobs = []
	# chunk_size = (rc.shape[1])/N_jobs
	# for i in range(N_jobs):
		# ts.append(threading.Thread(target=chunk_ML, args=(lock, i*chunk_size, (i+1)*chunk_size, splits_sorted, splits_inverted, one_Ls[:,i*chunk_size:(i+1)*chunk_size], zero_Ls[:,i*chunk_size:(i+1)*chunk_size])))
		# p = multiprocessing.Process(target=chunk_ML, args=(return_dict, i, splits_sorted, splits_inverted, one_Ls[:,i*chunk_size:(i+1)*chunk_size], zero_Ls[:,i*chunk_size:(i+1)*chunk_size]))
		# p.start()

	# for jb in jobs:
		# jb.join()
	# print return_dict.values()
	# for key in return_dict:
		# print (return_dict[key][1]).shape
	# t1 = threading.Thread(target=chunk_ML, args=(0, splits_sorted, splits_inverted, one_Ls[:,0:100000], zero_Ls[:,0:100000]))
	# t2 = threading.Thread(target=chunk_ML, args=(100000, splits_sorted, splits_inverted, one_Ls[:,100000:], zero_Ls[:,100000:]))

	# t1.start()
	# t2.start()

	# t1.join()
	# t2.join()


	# L_wg = 0
	# (new_genotypes, L) = one_site_mat(splits=split_list, ones=one_Ls[:,0], zeros=zero_Ls[:,0], counts=split_counts)
	# # (new_genotypes, L) = one_site_dot(split_arr=split_list, ones=one_Ls[:,0], zeros=zero_Ls[:,0])
	# # (new_genotypes, L) = one_site_ML_fast(split_arr=split_list, site_rc=rc[:,0], mu0=mu0, mu1=mu1, fn=fn, fp=fp, MD=MD)
	# # (new_genotypes, L) = one_site_ML(tree=tree, site_rc=rc[:,0], cell_names=names, mu0=mu0, mu1=mu1, fn=fn, fp=fp, MD=MD)
	# L_wg+=L
	# for p in range(1,rc.shape[1]):
	# 	# print p
	# 	# if ((p+1)%10000)==0:
	# 		# print "Finished up to genomic position "+str(p+1)
	# 	(tmp, L) = one_site_mat(splits=split_list, ones=one_Ls[:,p], zeros=zero_Ls[:,p], counts=split_counts)
	# 	# (tmp, L) = one_site_dot(split_arr=split_list, ones=one_Ls[:,p], zeros=zero_Ls[:,p])
	# 	# (tmp, L) = one_site_ML_fast(split_arr=split_list, site_rc=rc[:,p], mu0=mu0, mu1=mu1, fn=fn, fp=fp, MD=MD)
	# 	# (tmp, L) = one_site_ML(tree=tree, site_rc=rc[:,p], cell_names=names, mu0=mu0, mu1=mu1, fn=fn, fp=fp, MD=MD)
	# 	L_wg+=L
	# 	new_genotypes = np.vstack((new_genotypes,tmp))
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

# def single_column_opt(index, imputed_mat, read_count_col, fp, fn, missing_data_thr, K_vios, mu0, mu1, ones_, zeros_):
# 	fp_decimal = Decimal(fp)
# 	fn_decimal = Decimal(fn)
# 	n = imputed_mat.shape[0]
# 	l = imputed_mat.shape[1]
# 	missing_data_threshold = missing_data_thr

# 	model = Model("model")
# 	B = {}
# 	Y = []
# 	obj = LinExpr()
# 	print("Add variables to the model")
# 	for i in range(len(read_count_col)):
# 		Y.append(model.addVar(vtype=GRB.BINARY, name="Y[%d]" % i))
# 	for p in range(l):
# 		if p!=index:
# 			if p>index:
# 				B["("+str(index)+","+str(p)+","+str(1)+")"]=model.addVar(vtype=GRB.BINARY)
# 				B["("+str(index)+","+str(p)+","+str(2)+")"]=model.addVar(vtype=GRB.BINARY)
# 				B["("+str(index)+","+str(p)+","+str(3)+")"]=model.addVar(vtype=GRB.BINARY)
# 				# model.addConstr(2>=B["("+str(rem_col)+","+str(p)+","+str(1)+")"]+B["("+str(rem_col)+","+str(p)+","+str(2)+")"]+B["("+str(rem_col)+","+str(p)+","+str(3)+")"])
# 			else:
# 				B["("+str(p)+","+str(index)+","+str(1)+")"]=model.addVar(vtype=GRB.BINARY)
# 				B["("+str(p)+","+str(index)+","+str(2)+")"]=model.addVar(vtype=GRB.BINARY)
# 				B["("+str(p)+","+str(index)+","+str(3)+")"]=model.addVar(vtype=GRB.BINARY)
# 	print("Add constraints to the model")
# 	for p in range(l):
# 		if p!=index:
# 			if p>index:
# 				model.addConstr(2>=B["("+str(index)+","+str(p)+","+str(1)+")"]+B["("+str(index)+","+str(p)+","+str(2)+")"]+B["("+str(index)+","+str(p)+","+str(3)+")"])
# 				for cell in range(n):
# 					if imputed_mat[cell][p]==1:
# 						model.addConstr(B["("+str(index)+","+str(p)+","+str(1)+")"]>=1-Y[cell])
# 						model.addConstr(B["("+str(index)+","+str(p)+","+str(3)+")"]>=Y[cell])
# 					else:
# 						model.addConstr(B["("+str(index)+","+str(p)+","+str(2)+")"]>=Y[cell])
# 			else:
# 				model.addConstr(2>=B["("+str(p)+","+str(index)+","+str(1)+")"]+B["("+str(p)+","+str(index)+","+str(2)+")"]+B["("+str(p)+","+str(index)+","+str(3)+")"])
# 				for cell in range(n):
# 					if imputed_mat[cell][p]==1:
# 						model.addConstr(B["("+str(p)+","+str(index)+","+str(2)+")"]>=1-Y[cell])
# 						model.addConstr(B["("+str(p)+","+str(index)+","+str(3)+")"]>=Y[cell])
# 					else:
# 						model.addConstr(B["("+str(p)+","+str(index)+","+str(1)+")"]>=Y[cell])
# 	print("Build the objective function")
# 	for i in range(n):
# 		############# This line accounts for the missing data ##############
# 		if read_count_col[i][0]+read_count_col[i][1]>=missing_data_threshold: 
# 			r = int(read_count_col[i][0])
# 			v = int(read_count_col[i][1])
# 			AA = Binomial_pmf(v,r+v,mu0)
# 			BB = Binomial_pmf(v,r+v,mu1)
# 			#obj -= (Y[i][j])*np.float128(Decimal((fn_decimal/2)*AA+(1-fn_decimal/2)*BB).ln())
# 			obj -= (Y[i])*np.float128(Decimal((fn_decimal)*AA+(1-fn_decimal)*BB).ln())
# 			obj -= (1-Y[i])*np.float128(Decimal((1-fp_decimal)*AA+(fp_decimal)*BB).ln())
# 		else:
# 			pass
# 	model.update()
# 	print("Assign the objective function")
# 	model.setObjective(obj, GRB.MINIMIZE)
# 	#####################################################
# 	######## Set the parameters of the model ############
# 	#####################################################
# 	# model.Params.timeLimit = 1000
# 	# model.Params.FeasibilityTol = 0.001
# 	# model.Params.method=3
# 	#model.Params.Threads = 31
# 	#model.Params.ConcurrentMIP = 2
# 	model.Params.MIPGap=0.0
# 	# print(model.Params.FeasibilityTol)
# 	# print(model.Params.IntFeasTol)
# 	# print(model.Params.OptimalityTol)
# 	print("Optimize the model")
# 	model.optimize()
# 	print('IsMIP: %d' % model.IsMIP)
# 	if model.status == GRB.Status.INFEASIBLE:
# 		print("The model is infeasible")
# 	print("Solved with MIPFocus: %d" % model.Params.MIPFocus)
# 	print('Obj: %g' % model.objVal)
# 	print('MIP Gap: %g' % model.Params.MIPGap)
# 	for i in range(n):
# 		imputed_mat[i][index]=int(round(Y[i].x))
# 	gc.collect()
# 	L = 0
# 	for i in range(n):
# 		for j in range(l):
# 			if imputed_mat[i][j]==1:
# 				L+=ones_[i][j]
# 			else:
# 				L+=zeros_[i][j]
# 	return (imputed_mat, L)


# def optimize(read_count_mat, fp, fn, missing_data_thr, K_vios, mu0, mu1, constrained):
    
#     #########################################################################################################################
#     ############ The arguments include the error rates, the read count matrix, and the threshold for missing data ###########
#     #########################################################################################################################
#     fp_decimal = Decimal(fp)
#     fn_decimal = Decimal(fn)
#     n = read_count_mat.shape[0]
#     l = read_count_mat.shape[1]
#     R = [[0 for i in range(l)] for j in range(n)]
#     missing_data_threshold = missing_data_thr
#     ######################################################################
#     ######################### Build the model ############################
#     ######################################################################
#     model = Model("model")
#     B = {}
#     Y = []
#     V = {}
#     ########################################
#     ### Add the variables to the model #####
#     ########################################
#     obj = LinExpr()
#     vios = LinExpr()
#     print("Add variables to the model")
#     for i in range(n):
#     	Y.append([])
#     	for j in range(l):
#     		Y[i].append(model.addVar(vtype=GRB.BINARY, name="Y[%d,%d]" % (i,j)))
#     for p in range(l):
#     	# V[p]=model.addVar(vtype=GRB.BINARY) ####
#     	# vios+=V[p] ####
#     	q=p+1
#     	while q<l:
#     		V["("+str(p)+","+str(q)+")"]=model.addVar(vtype=GRB.BINARY)
#     		vios+=V["("+str(p)+","+str(q)+")"]
#     		for k in range(3):
#     			B["("+str(p)+","+str(q)+","+str(k+1)+")"]=model.addVar(vtype=GRB.BINARY, name="B["+str(p)+","+str(q)+","+str(k+1)+"]")
#     		q+=1
#     model.update()
#     ######################################
#     ### Add constraints to the model #####
#     ######################################
#     if constrained:
# 	    print("Add constraints to the model")
# 	    for p in range(l):
# 	    	q=p+1
# 	    	while q<l:
	    		
# 	    		model.addConstr(V["("+str(p)+","+str(q)+")"]>=B["("+str(p)+","+str(q)+","+str(1)+")"]+B["("+str(p)+","+str(q)+","+str(2)+")"]+B["("+str(p)+","+str(q)+","+str(3)+")"]-2)
# 	    		# model.addConstr(V[p]+V[q]>=B["("+str(p)+","+str(q)+","+str(1)+")"]+B["("+str(p)+","+str(q)+","+str(2)+")"]+B["("+str(p)+","+str(q)+","+str(3)+")"]-2) ###
# 	    		for taxon in range(n):
# 	    			####### The constraints which control the B variables #######
# 	    			model.addConstr(B["("+str(p)+","+str(q)+","+str(1)+")"]>=Y[taxon][q]-Y[taxon][p])
# 	    			model.addConstr(B["("+str(p)+","+str(q)+","+str(2)+")"]>=Y[taxon][p]-Y[taxon][q])
# 	    			model.addConstr(B["("+str(p)+","+str(q)+","+str(3)+")"]>=Y[taxon][p]+Y[taxon][q]-1)
# 	    		q=q+1
# 	    model.addConstr(vios<=K_vios)
#     # mu0=1e-3
#     # mu1=0.5
#     #################################################################
#     ################ Build the objective function ###################
#     #################################################################
#     print("Build the objective function")
#     for i in range(n):
#     	for j in range(l):
#     		############# This line accounts for the missing data ##############
#     		if read_count_mat[i][j][0]+read_count_mat[i][j][1]>=missing_data_threshold: 
#     			r = int(read_count_mat[i][j][0])
#     			v = int(read_count_mat[i][j][1])
#     			AA = Binomial_pmf(v,r+v,mu0)
#     			BB = Binomial_pmf(v,r+v,mu1)
#     			#obj -= (Y[i][j])*np.float128(Decimal((fn_decimal/2)*AA+(1-fn_decimal/2)*BB).ln())
#     			obj -= (Y[i][j])*np.float128(Decimal((fn_decimal)*AA+(1-fn_decimal)*BB).ln())
#     			obj -= (1-Y[i][j])*np.float128(Decimal((1-fp_decimal)*AA+(fp_decimal)*BB).ln())
#     		else:
#     			pass
#     model.update()
#     ##################################################
#     ########## Assign the objective function #########
#     ##################################################
#     print("Assign the objective function")
#     model.setObjective(obj, GRB.MINIMIZE)
#     #####################################################
#     ######## Set the parameters of the model ############
#     #####################################################
#     # model.Params.timeLimit = 1000
#     # model.Params.FeasibilityTol = 0.001
#     # model.Params.method=3
#     #model.Params.Threads = 31
#     #model.Params.ConcurrentMIP = 2
#     # model.Params.MIPGap=0.0
#     print(model.Params.FeasibilityTol)
#     print(model.Params.IntFeasTol)
#     print(model.Params.OptimalityTol)
#     ########################################################
#     ######### Optimize the model and report it #############
#     ########################################################
#     print("Optimize the model")
#     model.optimize()
#     print('IsMIP: %d' % model.IsMIP)
#     if model.status == GRB.Status.INFEASIBLE:
#     	print("The model is infeasible")
#     print("Solved with MIPFocus: %d" % model.Params.MIPFocus)
#     print("The noisy model has been optimized")
#     print('Obj: %g' % model.objVal)
#     print('MIP Gap: %g' % model.Params.MIPGap)
#     # if model.Params.MIPGap > 0.01:
#     # 	model.Params.timeLimit = 2000
#     # 	model.Params.MIPGap=0.01
#     # 	model.optimize()
#     #########################################################
#     ##### Save the final array given by the ILP solver ######
#     #########################################################
#     for i in range(n):
#     	for j in range(l):
#     		R[i][j] = int(round(Y[i][j].x))
#     gc.collect()
#     R = np.array(R)
#     return R

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
	# # monovar_mat = np.delete(monovar_mat,-1,axis=1)
	matrix_ = np.delete(matrix_,-1,axis=0)
	# print ete_nj
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
	missing_data_threshold = 10
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
	##############################################################################
	######################## Initialization from Monovar #########################
	# (monovar_mat, positions, chroms, names) = parse_vcf(file_address="/Users/edrisi/Documents/scalable_scVILP/data/16_neurons_WGS/output_monovar_chr1.vcf")
	##### fill the missing entries with zero
	# positions = positions[0:10000]
	# monovar_mat = monovar_mat[:,0:1000]
	# read_counts = read_counts[:,0:608250,:]
	# read_counts = read_counts[:,0:10000,:]
	n_cells = read_counts.shape[0]
	n_positions = read_counts.shape[1]
	print(n_positions)
	print(n_cells)
	# n_positions = (n_positions/50)*50
	# read_counts = read_counts[:,0:n_positions,:]
	# ratios = np.zeros((n_cells,n_positions))
	# for i in range(n_cells):
	# 	for j in range(n_positions):
	# 		if read_counts[i][j][0]+read_counts[i][j][1]!=0:
	# 			ratios[i][j] = float(read_counts[i][j][1])/float(read_counts[i][j][0]+read_counts[i][j][1])


	initialization = None
	# for chunk in range((n_cells/50)+1):
	# 	print("---------------- Inferring Chunk # "+str(chunk+1)+" by scVILP ----------------")
	# 	if chunk==0:
	# 		initialization = optimize(read_count_mat=read_counts[chunk*50:(chunk+1)*50,:], fp=fp_given, fn=fn_given, missing_data_thr=missing_data_threshold, K_vios=0, mu0=1e-3, mu1=0.5, constrained=True)
	# 	elif chunk==(n_cells/50):
	# 		initialization = np.concatenate((initialization,optimize(read_count_mat=read_counts[chunk*50:,:], fp=fp_given, fn=fn_given, missing_data_thr=missing_data_threshold, K_vios=0, mu0=1e-3, mu1=0.5, constrained=True)), axis=0)
	# 	else:
	# 		initialization = np.concatenate((initialization,optimize(read_count_mat=read_counts[chunk*50:(chunk+1)*50,:], fp=fp_given, fn=fn_given, missing_data_thr=missing_data_threshold, K_vios=0, mu0=1e-3, mu1=0.5, constrained=True)), axis=0)

	init_time = time.time()


	# zero_Ls = np.zeros((read_counts.shape[0],read_counts.shape[1]))
	# one_Ls = np.zeros((read_counts.shape[0],read_counts.shape[1]))


	# (initialization, uL, one_Ls, zero_Ls) = ML_initialization(rc=read_counts, names=names, mu0=1e-3, mu1=0.5, fn=fn_given, fp=fp_given, MD=missing_data_threshold, one_Ls=one_Ls, zero_Ls=zero_Ls)
	# np.save(out_path+"zero.npy", zero_Ls)
	# np.save(out_path+"one.npy", one_Ls)
	# np.save(out_path+"init.npy", initialization)
	# np.save(out_path+"l.npy", uL)


	one_Ls = None
	zero_Ls = None
	uL = None
	initialization = None
	one_Ls = np.load(out_path+"one.npy")
	zero_Ls = np.load(out_path+"zero.npy")
	initialization = np.load(out_path+"init.npy")
	uL = np.load(out_path+"l.npy")

	
	# initialization = read_List(infile="./scVILP_chunks.csv")
	print(time.time()-init_time, " seconds for the initalization")
	# exact = optimize(read_count_mat=read_counts, fp=fp_given, fn=fn_given, missing_data_thr=missing_data_threshold, K_vios=0, mu0=1e-3, mu1=0.5, constrained=True)
	tree_time = time.time()
	(ete_nj_init, nj_tree_init) = get_tree(matrx=initialization, nams=names, out_path=out_path)
	print(time.time()-tree_time, " seconds for the NJ tree")
	print(ete_nj_init)
	# nj_tree_init.write(path="NJ.nex",schema="nexus")
	# write_List(lst=exact, outfile="./exact_ovarian.csv")
	# exact = read_List(infile="./exact_ovarian.csv")
	# exact = np.array(exact)

	# initialization = np.random.randint(2,size=(n_cells,n_positions))
	# for i in range(n_cells):
	# 	for j in range(n_positions):
	# 		if int(monovar_mat[i][j])==-1:
	# 			print "-1"
	# 			monovar_mat[i][j]==0
	# write_List(lst=initialization, outfile="./initialization_matrix.csv")
	# initialization = np.array(initialization)
	# initialization = initialization[:,0:n_positions]
	# print initialization.shape


	##### reconstruct the NJ tree given the initial results
	# best_mat = None
	# best_rc = None
	# best_names = None
	best_L = float("-inf")
	# Ls = []
	# current_tree = ete_nj_init
	# # monovar_mat = np.vstack((monovar_mat,[0 for i in range(n_positions)]))
	# print n_positions
	# n_positions = (n_positions/50)*50
	# read_counts = read_counts[:,0:n_positions,:]


	# initialization = np.vstack((initialization,[0 for i in range(n_positions)]))
	


	# # print len(ete_nj.get_tree_root().get_children())
	# # print ete_nj.get_tree_root().get_leaf_names(is_leaf_fn=None)
	# # one_site_ML(tree=ete_nj, site_rc=read_counts[:,0], cell_names=names, mu0=1e-3, mu1=0.5, fn=fn_given, fp=fp_given, MD=missing_data_threshold)
	# exact = optimize(read_count_mat=read_counts, fp=fp_given, fn=fn_given, missing_data_thr=missing_data_threshold, K_vios=0, mu0=1e-3, mu1=0.5, constrained=True)
	# (ete_nj, nj_tree_exact) = get_tree(matrx=exact, nams=names, out_path=out_path)
	start = time.time()
	(new_mat, Likelihood) = ML(tree=ete_nj_init, rc=read_counts, names=names, mu0=1e-3, mu1=0.5, fn=fn_given, fp=fp_given, MD=missing_data_threshold, one_Ls=one_Ls, zero_Ls=zero_Ls)
	# print Likelihood
	print(time.time()-start, " seconds for the final ML")
	# nj_tree_exact.write(path="exact.nex",schema="nexus")
	# # # print new_mat.shape
	# write_List(lst=new_mat, outfile="./NJ_supertree_repeat.csv")

	# # ### next iteration 
	# # check_InCompatibility(Matrix=new_mat)

	# ##################### Plotting 

	# (leaves1, leaves2) = sort_mat(mat_unsorted=initialization)
	# # for i in range(n_cells):
	# # 	for j in range(n_positions):
	# # 		if read_counts[i][j][0]+read_counts[i][j][1]<missing_data_threshold:
	# # 			initialization[i][j] = -1
	# # (mat1, uL) = ML_initialization(rc=read_counts, names=names, mu0=1e-3, mu1=0.5, fn=fn_given, fp=fp_given, MD=missing_data_threshold)
	# # print uL
	# mat1 = copy.copy(initialization)
	# # mat1 = copy.copy(exact)
	# mat2 = copy.copy(new_mat)

	# mat1 = mat1[:,leaves1]
	# mat1 = mat1[leaves2,:]
	# mat2 = mat2[:,leaves1]
	# mat2 = mat2[leaves2,:]

	# names_ = copy.copy(names)
	# names_ = [names_[i] for i in leaves2]

	# fig, (ax1, ax2) = plt.subplots(1,2)
	# fig.set_size_inches(11, 6)
	# sns.heatmap(mat1, ax=ax1, cmap='Blues')
	# sns.heatmap(mat2, ax=ax2, cmap='Blues')
	# ax1.set_ylabel("Cells")
	# ax1.set_xlabel("Genomic Positions")
	# ax2.set_ylabel("Cells")
	# ax2.set_xlabel("Genomic Positions")
	# # ax1.set_title("scVILP")
	# # ax2.set_title("NJ with MLE")
	# plt.savefig("heatmap_iteration0.png", dpi=300)

	# ##################### Plotting

	# # random.shuffle(leaves1)
	# # random.shuffle(leaves2)

	# # new_mat = new_mat[:,leaves1]
	# # new_mat = new_mat[leaves2,:]
	# # read_counts = read_counts[:,leaves1,:]
	# # read_counts = read_counts[leaves2,:,:]

	# # names = [names[i] for i in leaves2]
	
	# mats = [copy.copy(new_mat)]
	# # rcs = [copy.copy(read_counts)]
	# # names_s = [copy.copy(names)]

	# Ls.append(Likelihood)

	if Likelihood>best_L:
	# # 	best_mat = copy.copy(new_mat)
		best_L = Likelihood
	# # 	best_rc = copy.copy(read_counts)
	# # 	best_names = copy.copy(names)



	# new_list = tree_move(tr=ete_nj_init, k=16, move=0)
	# index = 0
	# counter = 0
	# for tree_ in new_list:
	# 	(mat_, Likelihood) = ML(tree=tree_, rc=read_counts, names=names, mu0=1e-3, mu1=0.5, fn=fn_given, fp=fp_given, MD=missing_data_threshold, one_Ls=one_Ls, zero_Ls=zero_Ls)
	# 	Ls.append(Likelihood)
	# 	if Likelihood>best_L:
	# 		best_L = Likelihood
	# 		index = counter
	# 	counter+=1
	# print(new_list[index])

	n_iterations = 2000
	best_Ls = [best_L]
	num_bests = [1]
	stack  = [ete_nj_init]
	top_ids = set()
	top_ids.add(ete_nj_init.get_topology_id())
	# best_mat = new_mat
	# best_tree = ete_nj_init
	# mats_ = [new_mat]
	for it in range(n_iterations):
		print("------------------ iteration # "+str(it+1)+" ------------------")
		# fig = plt.figure(figsize=(3,3))
		# ax = fig.add_subplot(111)
		# sns.heatmap(mats[it].T, ax=ax, cmap='Blues')
		# plt.savefig("./mat0.png",dpi=300)
		# tmp = ML_initialization(rc=read_counts, names=names, mu0=1e-3, mu1=0.5, fn=fn_given, fp=fp_given, MD=missing_data_threshold)
		# for chunk in range(n_positions/50):
		# 	print("---------------- Inferring Chunk # "+str(chunk+1)+" by scVILP ----------------")
		# 	if chunk==0:
		# 		tmp = optimize(read_count_mat=read_counts[:,chunk*50:(chunk+1)*50], fp=fp_given, fn=fn_given, missing_data_thr=missing_data_threshold, K_vios=0, mu0=1e-3, mu1=0.5, constrained=True)
		# 	else:
		# 		tmp = np.concatenate((tmp,optimize(read_count_mat=read_counts[:,chunk*50:(chunk+1)*50], fp=fp_given, fn=fn_given, missing_data_thr=missing_data_threshold, K_vios=0, mu0=1e-3, mu1=0.5, constrained=True)), axis=1)

		# fig = plt.figure(figsize=(3,3))
		# ax = fig.add_subplot(111)
		# sns.heatmap(tmp.T, ax=ax, cmap='Blues')
		# plt.savefig("./scVILP_mat0.png",dpi=300)
		# # print mats[it].shape
		# columns = [i for i in range(n_positions)]
		# columns_wr = [random.choice(columns) for _ in range(n_positions)]
		# print columns_wr

		# (ete_nj, nj_tree) = get_tree(matrx=mats[it], nams=names, out_path=out_path)
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
					# ms.append(mat_)
			print("SPR")
			tree_list = SPR.Main(in_tree=item_, N=n_cells, N_dest=n_cells)
			for tree_ in tree_list:
				if tree_.get_topology_id() not in top_ids:
					top_ids.add(tree_.get_topology_id())
					(mat_, Likelihood) = ML(tree=tree_, rc=read_counts, names=names, mu0=1e-3, mu1=0.5, fn=fn_given, fp=fp_given, MD=missing_data_threshold, one_Ls=one_Ls, zero_Ls=zero_Ls)
					Ls.append(Likelihood)
					ts.append(tree_)
					# ms.append(mat_)
		# max_ = float("-inf")
		# if len(Ls)!=0:
			# max_ = max(Ls)
		print(len(Ls))
		max_ = max(Ls)
		# else:
		# 	print("NNI did not propose a better tree")
		if max_ > best_L:
			print("found better tree(s)")
			best_L = max_
			best_Ls.append(best_L)
			stack = []
			# mats_ = []
			for i in range(len(ts)):
				if Ls[i]==max_:
					# for node in ts[i].traverse():
					# 	if len(node.get_children())==1:
					# 		print("non binary tree")
					# 		break
					stack.append(ts[i])
					# mats_.append(ms[i])
			num_bests.append(len(stack))
		elif max_ == best_L:
			print("found equally good tree(s)")
			best_Ls.append(best_L)
			for i in range(len(ts)):
				if Ls[i]==max_:
					# for node in ts[i].traverse():
					# 	if len(node.get_children())==1:
					# 		print("non binary tree")
					# 		break
					stack.append(ts[i])
					# mats_.append(ms[i])
			num_bests.append(len(stack))
		else:
			print("no more better proposed trees")
			print("terminating the search...")
			break
			# Ls = []
			# ts = []
			# for tre_ in stack:
			# 	tree_list = SPR.Main(in_tree=tre_, N=n_cells, N_dest=n_cells)
			# 	for tree__ in tree_list:
			# 		if tree__.get_topology_id() not in top_ids:
			# 			top_ids.add(tree_.get_topology_id())
			# 			(mat_, Likelihood) = ML(tree=tree__, rc=read_counts, names=names, mu0=1e-3, mu1=0.5, fn=fn_given, fp=fp_given, MD=missing_data_threshold, one_Ls=one_Ls, zero_Ls=zero_Ls)
			# 			Ls.append(Likelihood)
			# 			ts.append(tree__)
			# max_ = float("-inf")
			# if len(Ls)!=0:
			# 	max_ = max(Ls)
			# else:
			# 	print("SPR did not propose a better tree")
			# 	break
			# if max_ > best_L:
			# 	print("SPR could find better tree(s)")
			# 	best_L = max_
			# 	best_Ls.append(best_L)
			# 	stack = []
			# 	for i in range(len(ts)):
			# 		if Ls[i]==max_:
			# 			stack.append(ts[i])
			# 	num_bests.append(len(stack))
			# elif max_ == best_L:
			# 	best_Ls.append(best_L)
			# 	for i in range(len(ts)):
			# 		if Ls[i]==max_:
			# 			stack.append(ts[i])
			# 	num_bests.append(len(stack))
			# else:
			# 	print("One time SPR rearrangement is not enough")
			# 	break 
		# print ete_nj
		
		# fig = plt.figure(figsize=(3,3))
		# ax = fig.add_subplot(111)
		# sns.heatmap(mat_.T, ax=ax, cmap='Blues')
		# plt.savefig("./NJ_mat0.png",dpi=300)
		
		# # print mat_.shape
		# (mat_, Likelihood) = single_column_opt(index=it%n_positions, imputed_mat=mats[it], read_count_col=read_counts[:,it%n_positions,:], fp=fp_given, fn=fn_given, missing_data_thr=missing_data_threshold, K_vios=0, mu0=1e-3, mu1=0.5, ones_=one_Ls, zeros_=zero_Ls)
		# print Likelihood
		# Ls.append(Likelihood)
		# mats.append(copy.copy(mat_))
		# # rcs.append(copy.copy(read_counts))
		# # names_s.append(copy.copy(names))
		# if Likelihood>best_L:
		# 	# best_mat = copy.copy(mat_)
		# 	best_L = Likelihood
		# 	best_rc = copy.copy(read_counts)
		# 	best_names = copy.copy(names)

		# mat3 = ML_initialization(rc=read_counts, names=names, mu0=1e-3, mu1=0.5, fn=fn_given, fp=fp_given, MD=missing_data_threshold)
		# fig = plt.figure(figsize=(3,3))
		# ax = fig.add_subplot(111)
		# sns.heatmap(mat1.T, ax=ax, cmap='Blues')
		# plt.savefig("./raw_mat0.png",dpi=300)


		# if it==n_iterations-1:

		# 	# (leaves1, leaves2) = sort_mat(mat_unsorted=mat_)
		# 	# mat1 = copy.copy(initialization)
		# 	print Ls.index(max(Ls))
		# 	mat2 = copy.copy(mats[Ls.index(max(Ls))])
		# 	# names_ = copy.copy(names)

		# 	# mat1 = mat1[:,leaves1]
		# 	# mat1 = mat1[leaves2,:]
		# 	mat2 = mat2[:,leaves1]
		# 	mat2 = mat2[leaves2,:]
		# 	# mat3 = mat3[:,leaves1]
		# 	# mat3 = mat3[leaves2,:]

		# 	# names_ = [names_[i] for i in leaves2]

		# 	fig, (ax1, ax2) = plt.subplots(1,2)
		# 	fig.set_size_inches(11, 6)
		# 	sns.heatmap(mat1, ax=ax1, cmap='Blues')
		# 	sns.heatmap(mat2, ax=ax2, cmap='Blues')
		# 	ax1.set_ylabel("Cells")
		# 	ax1.set_xlabel("Genomic Positions")
		# 	ax2.set_ylabel("Cells")
		# 	ax2.set_xlabel("Genomic Positions")
		# 	# sns.heatmap(mat3.T, ax=ax3, cmap='Blues', xticklabels=names)
		# 	plt.savefig("heatmap_iteration"+str(it+1)+".png", dpi=300)

		# mat_ = mat_[:,leaves1]
		# mat_ = mat_[leaves2,:]
		# read_counts = read_counts[:,leaves1,:]
		# read_counts = read_counts[leaves2,:,:]

	print(best_L)
	print(best_Ls)
	print(num_bests)
	print((stack[0]).write(format=9))
	(mat_, Likelihood) = ML(tree=stack[0], rc=read_counts, names=names, mu0=1e-3, mu1=0.5, fn=fn_given, fp=fp_given, MD=missing_data_threshold, one_Ls=one_Ls, zero_Ls=zero_Ls)
	### write the vcf file
	#### read_counts, alts, refs, chroms, positions, names, depths
	gen_VCF(out_dir=out_path, genotype_mat=mat_, read_count_mat_=read_counts, chrs=chroms, posits=positions, alt_counts=alts, rfs=refs, ids=names, dps=depths)
	# print("Best Likelihood: "+str(best_L))
	# print Ls
	# print("The exact likelihood: "+str(-5449274.9516161932))



	# acc_time = time.time()
	# true_array = []
	# true_positions = []
	# with open(true_matrix_path,"r") as tfile:
	# 	for line in tfile:
	# 		if "#position" in line:
	# 			continue
	# 		else:
	# 			true_array.append([])
	# 			tmp = line.strip().split("\t")
	# 			true_positions.append(int(tmp[0]))
	# 			for mut in tmp[1:]:
	# 				if int(mut)!=0 and int(mut)!=4:
	# 					true_array[len(true_array)-1].append(1)
	# 				else:
	# 					true_array[len(true_array)-1].append(0)
	# true_array=np.array(true_array)
	# true_array=true_array.T
	# True_list = copy.copy(true_array)
	# #### Check the final result ######
	# #### Check the false positive and false negative rate in the intferred matrix ######
	# false_positive_=0
	# false_negative_=0
	# true_positive_=0
	# # print(positions)
	# a = set(positions)
	# b = set(true_positions)
	
	# # notin_inferred_in_true_pos = list(b - a)


	# intersect_pos = list(a.intersection(b))
	# indx_arr = [positions.index(intersect_pos[j]) for j in range(len(intersect_pos))]
	# best_submat = best_mat[:,indx_arr]
	# indx_arr_true = [true_positions.index(intersect_pos[j]) for j in range(len(intersect_pos))]
	# true_submat = True_list[:,indx_arr_true]
	# true_positive_ = sum(sum(best_submat + true_submat==2))
	# false_negative_ = sum(sum(true_submat - best_submat==1))
	# false_positive_ = sum(sum(best_submat - true_submat==1))


	# # in_inferred_notin_true_pos = list(a - b)
	# f = [i for i in range(n_positions)]
	# residual_indices = list(set(f)-set(indx_arr))
	# residual_inference = best_mat[:,residual_indices]
	# false_positive_+= residual_inference.sum()

	# t = [i for i in range(len(true_positions))]
	# residual_indices_t = list(set(t)-set(indx_arr_true))
	# residual_true = True_list[:,residual_indices_t]
	# false_negative_+= residual_true.sum()


	# # indices_to_remove = []
	# # for j in range(len(true_positions)):
	# # 	index_ = positions.index(true_positions[j])
	# # 	indices_to_remove.append(index_)
	# # 	print j
	# # 	for i in range(n_cells):
	# # 		if best_mat[i][index_]==1 and True_list[i][j]==0:
	# # 			false_positive_+=1
	# # 		elif best_mat[i][index_]==0 and True_list[i][j]==1:
	# # 			false_negative_+=1
	# # 		elif best_mat[i][index_]==1 and True_list[i][j]==1:
	# # 			true_positive_+=1
	# # best_mat = np.delete(best_mat, indices_to_remove, axis=1)
	# # n_positions = best_mat.shape[1]
	# # for j in range(n_positions):
	# # 	print j
	# # 	for i in range(n_cells):
	# # 		if best_mat[i][j]==1:
	# # 			false_positive_+=1


	# precision_ = float(true_positive_)/float(true_positive_+false_positive_)
	# recall_ =  float(true_positive_)/float(true_positive_+false_negative_)
	# F1_ = 2*float(precision_*recall_)/float(precision_+recall_)
	# print(precision_, recall_, F1_)
	# print(time.time()-acc_time, " seconds for accuracy measurements")


	






