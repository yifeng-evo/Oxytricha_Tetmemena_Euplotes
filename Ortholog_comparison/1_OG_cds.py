from Bio import SeqIO
Oxy_cds=open("oxytri_mac.CDS.fasta","r")
Ewoo_cds=open("EVM.all.CDS.fasta","r")
Tet_cds=open("tetpusstrmac_2019.CDS.fasta","r")
OGls=open("Single_Copy_Orthologue_Sequences_ls","r")

Oxyd={}
Ewood={}
Tetd={}
for record in SeqIO.parse(Oxy_cds,"fasta"):
	Oxyd[record.id]=record

for record in SeqIO.parse(Ewoo_cds,"fasta"):
	Ewood[record.id]=record

for record in SeqIO.parse(Tet_cds,"fasta"):
	Tetd[record.id]=record


for line in OGls:
	if line[0]=="O":
		OG=line[0:9]
		f=open("./Single_Copy_Orthologue_Sequences/"+OG+".cds.fna","a")
		OGfile="./Single_Copy_Orthologue_Sequences/"+line[:-1]
		for record in SeqIO.parse(OGfile,"fasta"):
			if "evm" in record.id:
				ewooid=record.id.split("_")[0]
				SeqIO.write(Ewood[ewooid],f,"fasta")
			elif "TET" in record.id:
				tetid=record.id
				SeqIO.write(Tetd[tetid],f,"fasta")
				f.close()
			else:
				oxyid=record.id.split("_")[0]
				SeqIO.write(Oxyd[oxyid],f,"fasta")				

Oxy_cds.close()
Ewoo_cds.close()
Tet_cds.close()
OGls.close()



