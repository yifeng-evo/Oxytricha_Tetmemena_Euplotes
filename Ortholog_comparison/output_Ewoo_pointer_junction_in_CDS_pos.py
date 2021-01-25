from Bio import SeqIO
import os.path
from os import path
workingdir="./Ewoo_gff3/"
f=open("./Orthogroups/Orthogroups_SingleCopyOrthologues_3species.txt","r")
ewoo_gff3=open("/N/u/yif/Carbonate/yif/3species_ms/MAC_genomes/Gene_annotation/sorted_Ewoo_EVM_telov4.4.CDS.gff3","r")
Ewoo_cds_fna=open("/N/u/yif/Carbonate/yif/3species_ms/MAC_genomes/Gene_annotation/Ewoo_EVM_telov4.4.CDS.fna","r")
eijunction_cds=open("Ewoo_pointer_labeled_in_CDS_pos.txt","w")
pointer_junction=open("concise_below30master_Scramble_sdrap_ewoo_11032020_pid95_add90_prod_pointers.bed","r")

ewoopd={}
for line in f:
	ewooP=line.split("\t")[1] #ewoo
	OG=line.split("\t")[0]
	ewoopd[ewooP]=OG #with underscore

ewoop2c={}
for line in ewoo_gff3:
	id_s=line.find("Parent=")
#	id_e=line.rfind(";")
	ewooPid=line[id_s+7:-1]#including.1,but not _1
	if ewooPid+"_1" in ewoopd:
		# if line.split("\t")[2]=="intron":
		# 	f1=open(workingdir+ewooPid+".intron.gff3","a")
		# 	f1.write(line)
		# 	f1.close()
		
		f2=open(workingdir+ewooPid+".CDS.gff3","a")
		f2.write(line)
		f2.close()
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


pointer_d={}
for line in pointer_junction:
	new_MAC=line.split("\t")[0]
	if new_MAC in ewoop2c.values():
		f1=open(workingdir+new_MAC+".pointer.gff3","a")
		f1.write(line)
		f1.close()
	


for ewoopid in ewoop2c:
	ewooC=ewoop2c[ewoopid]
	if path.exists(workingdir+ewooC+".pointer.gff3"):
		f3=open(workingdir+ewooC+".pointer.gff3","r")
		for line in f3:
			pointer_1=int(line.split("\t")[1])
			pointer_2=int(line.split("\t")[2])
			cds_p1=chromosome2cds(pointer_1,ewoopid,ewoocds_len[ewoopid])
			cds_p2=chromosome2cds(pointer_2,ewoopid,ewoocds_len[ewoopid])
			eijunction_cds.write(ewoopd[ewoopid+"_1"]+"\t"+ewoopid+"\t"+ewoop2c[ewoopid]+"\t"+str(ewoocds_len[ewoopid])+"\t"+str(cds_p1)+"\t"+str(cds_p2)+"\t"+str(pointer_1)+"\t"+str(pointer_2)+"\t"+line.split("\t")[3])# improved to include nonscrambled/scrambled 


f.close()
ewoo_gff3.close()
Ewoo_cds_fna.close()

eijunction_cds.close()





