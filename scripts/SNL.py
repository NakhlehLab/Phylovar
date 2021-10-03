import numpy as np 
import ete3 
from ete3 import Tree
import argparse
import random
import sys
import copy
import uuid

def Swap(tree, swap_tuple):
	nodeA_name = swap_tuple[0]
	nodeB_name = swap_tuple[1]
	nodeA = tree&nodeA_name 
	nodeB = tree&nodeB_name

	nodeA.name = nodeB_name
	nodeB.name = nodeA_name

	return tree


def Main(in_tree, N, sample=True, all_=False):

	tree_list = []
	# smallest tree should have 3 nodes 
	if len(in_tree)<=2:
		print("No valid SNL for the given tree")
		return tree_list
	# a list of leaf names 
	leaf_names = in_tree.get_leaf_names()
	# a list of tuples containing all possible and valid swaps
	swap_options = []
	# sample N pairs of nodes that are not sisters
	for p in range(len(leaf_names)):
		sister_name = None
		if ((in_tree&leaf_names[p]).get_sisters()[0].is_leaf()):
			sister_name = (in_tree&leaf_names[p]).get_sisters()[0].name
		q = p+1
		while q<len(leaf_names):
			if leaf_names[q]==sister_name:
				pass
			else:
				swap_options.append((leaf_names[p], leaf_names[q]))
			q+=1
	if sample:
		if N>len(swap_options):
			N = len(swap_options)
		selected_tuples = random.sample(swap_options, k=N)
		for tp in selected_tuples:
			# print(tp)
			t__ = copy.deepcopy(in_tree)
			tree_list.append(Swap(t__, tp))
	if all_:
		for tp in swap_options:
			t__ = copy.deepcopy(in_tree)
			tree_list.append(Swap(t__, tp))

	return tree_list

if __name__=="__main__":
	### Toy example
	### When using this module, please call the Main function
	t = Tree(name="root")

	Y = t.add_child(name="Y")
	Z = t.add_child(name="F")

	X = Y.add_child(name="X")
	V = Y.add_child(name="V")


	A = V.add_child(name="A")
	B = V.add_child(name="B")

	E = X.add_child(name="E")

	W = X.add_child(name="W")
	D = W.add_child(name="D")
	C = W.add_child(name="C")
	l = t.get_leaf_names()
	# print(l)
	# print (t&"C").get_sisters()

	# print(t)
	# t_arr = Main(in_tree=t, N=300, sample=False, all_=True)
	# for tr in t_arr:
	# 	print(tr)
	# print(len(t_arr))