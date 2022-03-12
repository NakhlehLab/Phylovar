import numpy as np
import sys
import os
import gc
import random
import re
import copy

def chr_extract(chr_):
	if chr_[-1]=="X" or chr_[-1]=="x":
		return 23
	elif chr_[-1]=="Y" or chr_[-1]=="y":
		return 24
	else:
		chrString = ""
		for i in chr_:
			if i.isdigit():
				chrString+=i
		return int(chrString)

''' These two functions are used for reading the mpileup files '''
def match(str_):
	tmp = 0
	tmp+=str_.count(",")
	tmp+=str_.count(".")
	return tmp
def mismatch(str_):
	''' first, find all the occurrences of the regular expressions 
	corresponding to the insertions and deletions'''
	insertions = re.findall(r'\+[0-9]+[ACGTNacgtn]+', str_)
	deletions = re.findall(r'\-[0-9]+[ACGTNacgtn]+', str_)
	s = copy.deepcopy(str_)
	for e in insertions:
		s = s.replace(e,"")
	for f in deletions:
		s = s.replace(f,"")
	### [As,Cs,Gs,Ts,Ns]
	alternates = [s.count("A")+s.count("a"),s.count("C")+s.count("c"),s.count("G")+s.count("g"),s.count("T")+s.count("t"), s.count("N")+s.count("n")]
	tmp=sum(alternates)

	return (tmp, alternates)

def Parse(cell_names_file, mpileup_file):
	#### inclusion map is a sorted tuple including the chromosome index and location (chr, location)
	names_file= open(cell_names_file,"r")
	names_ = []
	tags = []
	arr_names = names_file.readlines()
	names_file.close()
	for line in arr_names:
		if len(line.strip())!=0:
			tmp = line.strip().split("\t")
			names_.append(tmp[0])
			tags.append(tmp[1])

	pileup = open(mpileup_file,"r")

	read_counts_ = []
	chroms_ = []
	positions_ = []
	refs_ = []
	alts_ = []
	depth_ = []
	num_cells = len(names_)
	for i in range(num_cells):
		read_counts_.append([])
		depth_.append([])
		alts_.append([])
	print("parsing the mpileup file ...")
	with open(mpileup_file, "r") as infile:
		for line in infile:
			depth_arr = []
			rc_arr = []
			alt_arr = []
			nmc_count = 0
			sections = line.strip().split('\t')
			current_pos = int(sections[1])
			positions_.append(current_pos)
			chroms_.append(chr_extract(sections[0]))
			refs_.append(sections[2])
			split = sections[3:]
			for i in range(num_cells):
				if int(split[3*i])==0:
					rc_arr.append([0,0])
					depth_arr.append(0)
					alt_arr.append([0,0,0,0,0])
				else:
					(mis_val, alt) = mismatch(split[3*i+1])
					rc_arr.append([match(split[3*i+1]),mis_val])
					depth_arr.append(int(split[3*i]))
					alt_arr.append(alt)

			for i in range(num_cells):
				read_counts_[i].append(rc_arr[i])
				depth_[i].append(depth_arr[i])
				alts_[i].append(alt_arr[i])

	read_counts_=np.array(read_counts_)
	print("parsing the mpileup file is done")

	gc.collect()

	return (read_counts_, alts_, refs_, chroms_, positions_, names_, depth_, tags)
