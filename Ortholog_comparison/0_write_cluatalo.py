ls=open("Single_Copy_Orthologue_Sequences_ls","r")
cluscript=open("clustal_omega.sh","w")
for line in ls:
	OG=line.split(".")[0]
	cluscript.write("/N/u/yif/Karst/clustalo --in "+line[:-1]+" --out "+OG+".aln.faa --threads 8 --force --outfmt=clu --iter=3\n")

ls.close()
cluscript.close()
