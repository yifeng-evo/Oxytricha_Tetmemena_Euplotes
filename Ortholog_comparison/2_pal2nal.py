OGls=open("/N/u/yif/Karst/scratch/Euplotes_ortholog_comparison/orthofinder_Oxy_Ewoo_longestisoform/oxy_tet_ewoo/OrthoFinder/Results_Nov19/Single_Copy_Orthologue_Sequences_ls","r")
pal2nalscript=open("pal2nal.sh","w")
for line in OGls:
	OG=line.split(".")[0]
	pal2nalscript.write("/N/u/yif/Carbonate/scratch/Euplotes_ortholog_comparison/orthofinder_Oxy_Ewoo_longestisoform/pal2nal/pal2nal.pl "+OG+".aln.faa "+OG+".cds.fna -codontable 10,6,6 1>"+OG+".pal2nal.out 2>"+OG+".pal2nal.err\n")

OGls.close() 
pal2nalscript.close()
