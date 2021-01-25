from Bio import SeqIO
import os.path
from os import path
workingdir="./Tet_gff3/"
f=open("./Orthogroups/Orthogroups_SingleCopyOrthologues_3species.txt","r")
ewoo_gff3=open("/N/u/yif/Carbonate/yif/3species_ms/MAC_genomes/Gene_annotation/tetpusstrmac_2019_telo_longestisoform_intron.gff3","r")
Ewoo_cds_fna=open("/N/u/yif/Carbonate/yif/3species_ms/MAC_genomes/Gene_annotation/tetpusstrmac_2019_telo_longestisoform.CDS.fna","r")
eijunction_cds=open("Tet_exon_junction_labeled_in_CDS_pos.txt","w")

ewoopd={}
for line in f:
	ewooP=line.split("\t")[3][:-1]
	OG=line.split("\t")[0]
	ewoopd[ewooP]=OG #with underscore

ewoop2c={}
for line in ewoo_gff3:
	id_s=line.find("Parent=")
	#id_e=line.rfind("\n")
	ewooPid=line[id_s+7:-1]#including.1,but not _1
	if ewooPid+"_1" in ewoopd:
		if line.split("\t")[2]=="intron":
			f1=open(workingdir+ewooPid+".intron.gff3","a")
			f1.write(line)
			f1.close()
		
		# f2=open(workingdir+ewooPid+".CDS.gff3","a")
		# f2.write(line)
		# f2.close()
		ewooC=line.split("\t")[0]
		ewoop2c[ewooPid]=ewooC


ewoocds_len={}
for record in SeqIO.parse(Ewoo_cds_fna,"fasta"):
	if record.id in ewoop2c:
		ewoocds_len[record.id]=len(str(record.seq))

def chromosome2cds(c_pos,p_id,length_cds): #oxy_c_pos is 0. This is good for range calculation. single nucleotide canbe represented as :[5,6) oxy_gff only contains cds lines
	cds_block=0
	cds_dict={}
	end_block=-1
	cds_pos=-1

	gff=open(workingdir+p_id+".CDS.gff3","r")

	for cds in gff:
		start=int(cds.split("\t")[3])
		end=int(cds.split("\t")[4])
		cds_dict[cds_block]=[start,end]
		strand=cds.split("\t")[6]
		if start<=c_pos<=end: 
			end_block=cds_block
			tail=c_pos-start
			i=0
			head=0
			while i < end_block:
				head=cds_dict[i][1]-cds_dict[i][0]+1+head
				i=i+1
			cds_pos=head+tail+1
			if strand =="-":
				cds_pos=length_cds-cds_pos+1
			break
		cds_block=cds_block+1

	gff.close()

	return cds_pos


# pointer_d={}
# for line in pointer_junction:
# 	new_MAC=line.split("\t")[0]
# 	if new_MAC in ewoop2c.values():
# 		f1=open(workingdir+new_MAC+".pointer.gff3","a")
# 		f1.write(line)
# 		f1.close()
	


for ewoopid in ewoop2c:
	ewooC=ewoop2c[ewoopid]
	if path.exists(workingdir+ewoopid+".intron.gff3"):
		f3=open(workingdir+ewoopid+".intron.gff3","r")
		for line in f3:
			intron_s=int(line.split("\t")[3])-1 #the junction is defined as the first in exon or last in exon
			intron_e=int(line.split("\t")[4])+1
			cds_s=chromosome2cds(intron_s,ewoopid,ewoocds_len[ewoopid])
#			eijunction_cds.write(ewoopd[ewoopid+"_1"]+"\t"+ewoopid+"\t"+ewoop2c[ewoopid][0]+"\t"+str(ewoop2c[ewoopid][1])+"\t"+str(cds_s)+"\n")
			cds_e=chromosome2cds(intron_e,ewoopid,ewoocds_len[ewoopid])
			eijunction_cds.write(ewoopd[ewoopid+"_1"]+"\t"+ewoopid+"\t"+ewoop2c[ewoopid]+"\t"+str(ewoocds_len[ewoopid])+"\t"+str(cds_s)+"\t"+str(cds_e)+"\t"+str(intron_s)+"\t"+str(intron_e)+"\n")

	


f.close()
ewoo_gff3.close()
Ewoo_cds_fna.close()

eijunction_cds.close()





