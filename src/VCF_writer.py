import numpy as np
from datetime import date

def gen_VCF(out_, genotype_mat, read_count_mat_, chrs, posits, alt_counts, rfs, ids, dps, n_, tags):
	# print dps
	print("dimension of the read count matrix: ",read_count_mat_.shape)
	print("dimension of the genotype matrix: ",genotype_mat.shape)
	d = {0:"A", 1:"C", 2:"G", 3:"T", 4:"N"}
	read_count_mat_ = np.array(read_count_mat_)
	n = read_count_mat_.shape[0]
	l = read_count_mat_.shape[1]
	now = date.today().strftime("%Y%m%d")

	out_f = open(out_, "w")
	header = "##fileformat=VCF\n##fileDate="+now+"\n##FILTER=<ID=LowQual,Description=\"Low quality\">\n##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth; some reads may have been filtered\">\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n##FORMAT=<ID=RO,Number=1,Type=Integer,Description=\"Reference allele observed count\">\n##FORMAT=<ID=AO,Number=1,Type=Integer,Description=\"Alternate allele observed count\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
	for cell in ids:
		header+="\t"+cell
	out_f.write(header+"\n")
	for pos in range(l):
		if np.count_nonzero(genotype_mat[:,pos])!=0:
			total_depth = 0
			tmp_str = str(chrs[pos])+"\t"+str(posits[pos])+"\t"+n_[pos]+"\t"+rfs[pos]+"\t"
			sum_ = [0,0,0,0,0]
			for id_indx in range(n):
				total_depth=total_depth+dps[id_indx][pos]
				if genotype_mat[id_indx][pos]==1:
					sum_ = [a+b for a,b in zip(sum_,alt_counts[id_indx][pos])]

			major_alt = ""
			if sum(sum_) == 0:
				major_alt = "*"
			if sum(sum_)!=0:
				major_alt = d[sum_.index(max(sum_))]

			tmp_str+=major_alt+"\t"
			out_f.write(tmp_str)
			out_f.write(".\t")
			out_f.write("PASS\t")
			out_f.write("DP:"+str(total_depth)+"\t")
			out_f.write("GT:DP:RO:AO\t")
			sum_string = "<"
			for id_indx in range(n):
				if genotype_mat[id_indx][pos]==1:
					out_f.write("0/1:")
					sum_string+="1"
				else:
					out_f.write("0/0:")
					sum_string+="0"
				out_f.write(str(dps[id_indx][pos])+":")
				out_f.write(str(read_count_mat_[id_indx][pos][0])+":"+str(read_count_mat_[id_indx][pos][1]))
				out_f.write("\t")
			sum_string+=">"
			out_f.write(sum_string+"\n")

	out_f.close()
	return 
