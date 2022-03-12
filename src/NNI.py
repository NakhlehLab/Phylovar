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

def perform_NNI(tree, selection, choice_):
	######## the move but it still remains as the root of the tree after the move 
	selected_node = None
	other_child = None
	D_component = None
	A_component = None

	selected_node = tree&selection
	ancestor = selected_node.up
	if choice_==1:
	### the first type of exchange
		B_component = selected_node.get_children()[0]
		C_component = selected_node.get_children()[1]
	else:
		B_component = selected_node.get_children()[1]
		C_component = selected_node.get_children()[0]

	A_component = selected_node.get_sisters()[0]
	B_component.detach()
	A_component.detach()
	selected_node.add_child(A_component)
	ancestor.add_child(B_component)


	return tree
	
def Main(in_tree, N):
	''' Since it is not possible to deal with all the 
	NNI rearrangements given a topology, we randomly sample 
	2N trees/topologies from all the possible rearrangements in each NNI.py call'''
	tree_list = []
	nodes = []
	internal_names = gen_unique_ids(num=len(in_tree)-2)
	c = 0
	for node in in_tree.traverse():
		if (not node.is_root()) and (not node.is_leaf()):
			c+=1
	selected_names = []
	if len(in_tree)<=2:
		print("No valid NNI rearrangement for the given tree")
		return tree_list
	else:
		counter = 0
		for node in in_tree.traverse():
			if (not node.is_root()) and (not node.is_leaf()):
				node.name = internal_names[counter]
				nodes.append(node)
				selected_names.append(node.name)
				counter+=1
		if N>len(nodes):
			N = len(nodes)
		samples = random.sample(range(len(nodes)), k=N)
		for i in samples:
			tmp_cpy_1 = copy.deepcopy(in_tree)
			tmp_cpy_2 = copy.deepcopy(in_tree)

			tree_list.append(perform_NNI(tree=tmp_cpy_1, selection=selected_names[i], choice_=1))
			tree_list.append(perform_NNI(tree=tmp_cpy_2, selection=selected_names[i], choice_=2))

	return tree_list

if __name__=="__main__":
	### Toy example
	### When using this module, please call the Main function
	### NNI must return 8 different topologies for this toy example
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

	print(t)
	t_arr = Main(in_tree=t, N=300)
	for tr in t_arr:
		print(tr)
	print(len(t_arr))
