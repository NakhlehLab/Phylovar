import numpy as np
import sys
import os
import gc
import re

def get_info(info_str):
	genes = []
	effects = []
	info_str = info_str.replace("ANN=", "")
	anns = info_str.strip().split(",")
	for ann in anns:
		secs = ann.split("|")
		if secs[2]=="HIGH" or secs[2]=="MODERATE":
			genes.append(secs[3])
			effects.append(secs[1])

	return (effects,genes)

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

def genotype(str_, mchar):
	p1 = r'^0([\/\|]0)*$'
	p2 = r'^\.([\/\|]\.)*$'
	p3 = r'^[0-9]+([\/\|][0-9]+)*$'
	if re.match(p1,str_):
		return 0
	elif re.match(p2,str_):
		return mchar
	elif re.match(p3, str_):
		return 1
	else:
		print("Unidentified genotype: "+str_)
		return None

def phred(ph):
	if ph == "./.":
		return 0.5
	else:
		ph = ph.split(",")
		g0 = int(ph[0])
		g01 = int(ph[1])
		g11 = int(ph[2])

		g0prob = pow(10,-1*g0/10)
		g01prob = pow(10,-1*g01/10)
		g11prob = pow(10,-1*g11/10)

		mutation_prob = g01prob + g11prob
		sum_prob = g0prob + g01prob + g11prob
		return mutation_prob/sum_prob

def parse_vcf(file_address, mchar, tail=True, phred_s=False, annotations=False):
	chroms = []
	positions = []
	names = []
	genotypes = []
	phred_scores = []
	genes_ = []
	effects_ = []
	with open(file_address, "r") as infile:
		for line in infile:
			if "#CHROM" in line:
				if "CT" in line:
					line = line.replace("\tCT", "")
				if "CN" in line:
					line = line.replace("\tCN", "")
				arr = line.strip().split('\t')[9:]
				names = arr
			if not "#" in line:
				genotypes.append([])
				phred_scores.append([])
				arr = line.strip().split('\t')
				current_chrom = chr_extract(arr[0])
				chroms.append(current_chrom)
				positions.append(int(arr[1]))
				if tail:
					for raw_gt in arr[9:-1]:
						gt_arr = raw_gt.split(":")
						genotypes[len(genotypes)-1].append(genotype(gt_arr[0],mchar))
						if phred_s:
							phred_scores[len(genotypes)-1].append(phred(gt_arr[-1]))
				else:
					for raw_gt in arr[9:]:
						gt_arr = raw_gt.split(":")
						genotypes[len(genotypes)-1].append(genotype(gt_arr[0],mchar))
						if phred_s:
							phred_scores[len(genotypes)-1].append(phred(gt_arr[-1]))
				if annotations:
					effs = []
					gens = []
					w = arr[7].split(";")
					if len(w)>1:
						(effs, gens) = get_info(w[1])
					else:
						pass
					effects_.append(effs)
					genes_.append(gens)

	genotypes = np.array(genotypes)
	if phred_s and not annotations:
		phred_scores = np.array(phred_scores)
		return (genotypes.T, positions, chroms, names, phred_scores.T)
	elif not phred_s and not annotations:
		return (genotypes.T, positions, chroms, names)
	elif not phred_s and annotations:
		return (genotypes.T, positions, chroms, names, genes_, effects_)
	else:
		return (genotypes.T, positions, chroms, names, phred_scores.T, genes_, effects_)

if __name__=="__main__":
	parse_vcf(file_address="/Users/edrisi/Documents/scalable_snv_calling/data/16_neurons/HC_turing/phylo_wg.vcf", mchar = -10)


