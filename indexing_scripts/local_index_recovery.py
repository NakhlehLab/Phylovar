import argparse

if __name__=="__main__":
	ap = argparse.ArgumentParser()
	ap.add_argument("-mpileup", "--mpileup", required=True,
		help="path to the original mpileup")
	ap.add_argument("-sciphi", "--sciphi", required=True,
		help="path to the output of SCIPhi after filtering non-informative sites named genotype_matrix.csv")
	ap.add_argument("-out", "--out", required=True,
		help="path to the output mpileup")
	args = ap.parse_args()
	mpileup_path = args.mpileup
	out_path = args.out
	genotype_path = args.sciphi

	outf = open(out_path, "w")
	genotype_f = open(genotype_path, "r")
	# the first line is the headers
	genotype_f.readline()
	original_mpileup = open(mpileup_path, "r")
	line_counter = 0
	original_line = original_mpileup.readline()

	for line in genotype_f:
		global_index = int(line.strip().split(",")[0])
		#### read the original mpileup until reaching the corresponding line
		while line_counter!=global_index:
			line_counter+=1
			original_line = original_mpileup.readline()
		outf.write(original_line)
	outf.close()