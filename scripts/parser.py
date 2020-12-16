import numpy as np
import sys
import os
import gc

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
	alternates = ""
	tmp=0
	tmp+=str_.count("A")
	tmp+=str_.count("C")
	tmp+=str_.count("G")
	tmp+=str_.count("T")
	tmp+=str_.count("N")
	tmp+=str_.count("a")
	tmp+=str_.count("c")
	tmp+=str_.count("g")
	tmp+=str_.count("t")
	tmp+=str_.count("n")
	if str_.count("A")!=0 or str_.count("a")!=0:
		alternates+="A"
	if str_.count("C")!=0 or str_.count("c")!=0:
		alternates+="C"
	if str_.count("T")!=0 or str_.count("t")!=0:
		alternates+="T"
	if str_.count("G")!=0 or str_.count("g")!=0:
		alternates+="G"
	if str_.count("N")!=0 or str_.count("n")!=0:
		alternates+="N"
	if alternates=="":
		alternates="*"
	return (tmp, alternates)

def Parse(cell_names_file, mpileup_file):
	names_file= open(cell_names_file,"r")
	names_ = []
	arr_names = names_file.readlines()
	names_file.close()
	for line in arr_names:
		if len(line.strip())!=0:
			tmp = line.strip().split("\t")
			names_.append(tmp[0])

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
		alts_.append([])
		depth_.append([])

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
					alt_arr.append("*")
					depth_arr.append(0)
				else:
					(mis_val, alt_str) = mismatch(split[3*i+1])
					rc_arr.append([match(split[3*i+1]),mis_val])
					alt_arr.append(alt_str)
					depth_arr.append(int(split[3*i]))

			for i in range(num_cells):
				read_counts_[i].append(rc_arr[i])
				alts_[i].append(alt_arr[i])
				depth_[i].append(depth_arr[i])

	read_counts_=np.array(read_counts_)
	gc.collect()
	return (read_counts_, alts_, refs_, chroms_, positions_, names_, depth_)