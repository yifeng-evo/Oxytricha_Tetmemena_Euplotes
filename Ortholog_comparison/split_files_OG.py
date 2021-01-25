#f1=open("Tet_exon_junction_labeled_in_CDS_pos.txt","r")
import sys
f2=open(sys.argv[1],"r")

def split_file_OG(folder,position_on_cds_file):
	for line in position_on_cds_file:
		OG=line.split("\t")[0]
		f=open(folder+OG+"."+sys.argv[3]+"_pos_on_cds.txt","a")
		f.write(line)
		f.close()

split_file_OG("./"+sys.argv[2]+"_pos_on_cds/",f2)
#split_file_OG("./Ewoo_pos_on_cds/",f2)


