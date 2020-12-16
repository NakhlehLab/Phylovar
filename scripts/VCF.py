import numpy as np
from datetime import date

def gen_VCF(out_dir, genotype_mat, read_count_mat_, chrs, posits, alt_counts, rfs, ids, dps):
	# print dps
	read_count_mat_ = np.array(read_count_mat_)
	n = read_count_mat_.shape[0]
	l = read_count_mat_.shape[1]
	now = date.today().strftime("%Y%m%d")

	out_f = open(out_dir+"snv.vcf", "w")
	header = "##fileformat=VCF\n##fileDate="+now+"\n##FILTER=<ID=LowQual,Description=\"Low quality\">\n##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth; some reads may have been filtered\">\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n##FORMAT=<ID=DPR,Number=1,Type=Integer,Description=\"Number of observation for each allele\">\n##FORMAT=<ID=RO,Number=1,Type=Integer,Description=\"Reference allele observed count\">\n##FORMAT=<ID=AO,Number=1,Type=Integer,Description=\"Alternate allele observed count\">\n#CHROM\tPOS\tID\tREF\tALT\tFILTER\tINFO\tFORMAT"
	for cell in ids:
		header+="\t"+cell
	out_f.write(header+"\n")
	for pos in range(l):
		total_depth = 0
		tmp_str = str(chrs[pos])+"\t"+str(posits[pos])+"\t*\t"+rfs[pos]+"\t"
		for id_indx in range(n):
			total_depth=total_depth+dps[id_indx][pos]
			tmp_set = []
			alt_string = alt_counts[id_indx][pos]
			for ch in alt_string:
				tmp_set.append(ch)
		alt_arr_ = list(set(tmp_set))
		if len(alt_arr_)!=0:
			for c in range(len(alt_arr_)):
				if c==len(alt_string)-1:
					tmp_str+=alt_arr_[c]
				else:
					tmp_str+=alt_string[c]+","
		else:
			tmp_str+="*"
		out_f.write(tmp_str+"\t")
		out_f.write("PASS\t")
		out_f.write("DP:"+str(total_depth)+"\t")
		out_f.write("GT:DP:DPR:RO:AO\t")
		for id_indx in range(n):
			if genotype_mat[id_indx][pos]==1:
				out_f.write("0/1:")
			else:
				out_f.write("0/0:")
			out_f.write(str(dps[id_indx][pos])+":")
			out_f.write(str(read_count_mat_[id_indx][pos][0]+read_count_mat_[id_indx][pos][1])+","+str(read_count_mat_[id_indx][pos][1])+":"+str(read_count_mat_[id_indx][pos][0])+":"+str(read_count_mat_[id_indx][pos][1]))
			if id_indx==n-1:
				out_f.write("\n")
			else:
				out_f.write("\t")



	out_f.close()
	return 