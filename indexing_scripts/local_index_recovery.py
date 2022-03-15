if __name__=="__main__":
	mpileup_path = "/shared/mae6/snv_calling_data/test/tnbc_global_idx.mpileup"
	index_path = "/shared/mae6/snv_calling_data/test/index.csv"
	out_path = "/shared/mae6/snv_calling_data/test/tnbc_local_idx.mpileup"
	genotype_path = "/shared/mae6/snv_calling_data/test/genotype_matrix.csv"
	
	genotypes = []
	line_lens = set()
	with open(genotype_path, "r") as genotypef:
		for line in genotypef:
			if "pos" in line:
				pass
			else:
				s = line.strip().split(",")
				line_lens.add(len(s))
				if len(genotypes)==0 or (s[0]!=genotypes[-1]):
					genotypes.append(s[0])
	genotypef.close()
	print(len(genotypes))
	if len(genotypes)==len(set(genotypes)):
		print("genotype matrix file does not have duplicates")
	else:
		print("genotype matrix file does have duplicates!!!")
	print("lengths of the lines in the genotype matrix file")
	print(line_lens)

	idx_dict = {}
	with open(index_path, "r") as idxf:
		for line in idxf:
			s = line.strip().split(",")
			idx_dict[s[0]] = s[2]
	idxf.close()
	print(len(idx_dict))

	line_counter = 0
	outf = open(out_path, "w")
	with open(mpileup_path, "r") as infile:
		for line in infile:
			if len(genotypes)==0:
				break
			elif str(line_counter)==genotypes[0]:
				s = line.split("\t")
				s[1] = idx_dict[str(line_counter)]
				outf.write("\t".join(s))
				genotypes.pop(0)
			line_counter+=1
	outf.close()



