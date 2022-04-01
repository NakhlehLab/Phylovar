import argparse

if __name__=="__main__":
        ap = argparse.ArgumentParser()
        ap.add_argument("-mpileup", "--mpileup", required=True,
                help="path to the original mpileup")
        ap.add_argument("-out", "--out", required=True,
                help="path to the output mpileup")
        args = ap.parse_args()
        mpileup_path = args.mpileup
        out_path = args.out
        line_counter = 0
        outf = open(out_path, "w")
        with open(mpileup_path, "r") as infile:
                for line in infile:
                        s = line.split("\t")
                        s[1] = str(line_counter)
                        outf.write("\t".join(s))
                        line_counter+=1
        outf.close()
