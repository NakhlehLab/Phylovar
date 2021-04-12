import numpy as np 
import ete3 
from ete3 import Tree
import argparse
import random
import sys
import copy
import uuid 

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


def perform_SPR(in_tr, selection, dest):
	returned_trees = []
	selected_node = in_tr&selection

	to_remove = selected_node.up
	sister = selected_node.get_sisters()[0]
	subtree = selected_node.detach()
	to_remove.delete()
	new_internals_names = []

	for node in in_tr.traverse():
		if (not node.is_root()) and (node!=sister):
			new_internals_names.append(node.name)
	if len(new_internals_names)==0:
		print("No valid SPR rearrangement for this branch")
		returned_trees =  None
	else:
		if dest>len(new_internals_names):
			dest = len(new_internals_names)
		samples = random.sample(range(len(new_internals_names)), dest)
		for sample in samples:
			tmp_in_tr = copy.deepcopy(in_tr)
			tmp_subtree = copy.deepcopy(subtree)
			s = tmp_in_tr&new_internals_names[sample]
			parent = s.up
			new_node = parent.add_child()
			new_node.add_child(tmp_subtree)
			s.detach()
			new_node.add_child(s)
			returned_trees.append(tmp_in_tr)

	return returned_trees

def Main(in_tree, N, N_dest):
	''' Since it is not possible to deal with all the SPR rearrangements 
	given a particular topology, we randomly sample N new topologies in each call'''
	tree_list = []
	tr = in_tree
	nodes = []
	internal_names = gen_unique_ids(num=len(in_tree)-1)
	selected_names = []

	if len(in_tree)<=3:
		print("No valid SPR rearrangement for the given tree")
		return tree_list
	else:
		counter = 0
		new_t = Tree(name="artificial_root")
		new_t.add_child(copy.deepcopy(in_tree))
		(new_t.get_children()[0]).name = internal_names[counter]
		counter+=1
		for node in new_t.traverse():
			if not node.is_root():
				if not (node.up).is_root():
					if node.is_leaf():
						nodes.append(node)
						selected_names.append(node.name)
					else:
						node.name = internal_names[counter]
						nodes.append(node)
						selected_names.append(node.name)
						counter+=1
		if N>len(nodes):
			N = len(nodes)
		samples = random.sample(range(len(nodes)), k=N)
		for i in samples:
			tmp_cpy = copy.deepcopy(new_t)
			tt = perform_SPR(in_tr=tmp_cpy, selection=selected_names[i], dest=N_dest)
			if tt!=None:

				for x in tt:
					if len(x.get_tree_root().get_children())==1:
						(x.get_tree_root().get_children()[0]).delete()
						tree_list.append(x)
					else:
						tree_list.append(x)
	return tree_list

if __name__=="__main__":
	#### NOTE: each of the internal nodes must have a name for this rearrangement
	### Toy example
	### SPR must return 64 different topologies for this toy example
	t = Tree(name="root")
	Z = t.add_child(name="Z")
	Y = t.add_child(name="Y")
	# F = t.add_child(name="F")

	X = Y.add_child(name="X")
	V = Y.add_child(name="V")

	A = V.add_child(name="A")
	B = V.add_child(name="B")

	E = X.add_child(name="E")

	W = X.add_child(name="W")
	D = W.add_child(name="D")
	C = W.add_child(name="C")

	print(t)
	lst = Main(in_tree=t, N=100, N_dest=300)
	for tr_ in lst:
		print(tr_)
	print(len(lst))
