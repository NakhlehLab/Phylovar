if __name__=="__main__":
	mpileup_path = "/shared/mae6/snv_calling_data/test/tnbc.mpileup"
	index_path = "/shared/mae6/snv_calling_data/test/index.csv"

	line_counter = 0
	outf = open(index_path, "w")
	with open(mpileup_path, "r") as infile:
		for line in infile:
			if line_counter!=0:
				outf.write("\n")
			s = line.strip().split("\t")
			outf.write(str(line_counter)+","+s[0]+","+s[1])
			line_counter+=1
	outf.close()
