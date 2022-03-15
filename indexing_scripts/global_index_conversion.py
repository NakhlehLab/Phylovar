if __name__=="__main__":
        mpileup_path = "/shared/mae6/snv_calling_data/test/tnbc.mpileup"
        out_path = "/shared/mae6/snv_calling_data/test/tnbc_global_idx.mpileup"
        line_counter = 0
        outf = open(out_path, "w")
        with open(mpileup_path, "r") as infile:
                for line in infile:
                        s = line.split("\t")
                        s[1] = str(line_counter)
                        outf.write("\t".join(s))
                        line_counter+=1
        outf.close()
