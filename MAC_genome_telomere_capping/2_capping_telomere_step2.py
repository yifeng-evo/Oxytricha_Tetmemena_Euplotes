from Bio import SeqIO
assembly=open("/ifs/scratch/biochem/ll_lab/yf2342/dedupe_98_Euplotes_woodruffi_trinity_spades_cap3_1023.fasta","r")
cappingCA=open("/ifs/scratch/biochem/ll_lab/yf2342/CAcapping_1023assembly.txt","r")
cappingGT=open("/ifs/scratch/biochem/ll_lab/yf2342/GTcapping_1023assembly.txt","r")
changedrecord=open("/ifs/scratch/biochem/ll_lab/yf2342/changedrecord.fasta","w")
capped_assembly=open("/ifs/scratch/biochem/ll_lab/yf2342/capped_dedupe_98_Euplotes_woodruffi_trinity_spades_cap3_1023.fasta","w")
CAmax={}
GTmax={}
for line in cappingCA:
	contigname=line.split("\t")[0]
	TAS=int(line.split("\t")[1])
	CAmax[contigname]=TAS

for line in cappingGT:
	contigname=line.split("\t")[0]
	TAS=int(line.split("\t")[1])
	GTmax[contigname]=TAS



for record in SeqIO.parse(assembly,"fasta"):
	name=record.id
	if name in CAmax and record.id in GTmax:
		record.seq="CCCCAAAACCCCAAAACCCCAAAACCCC"+record.seq[CAmax[name]:GTmax[name]]+"GGGGTTTTGGGGTTTTGGGGTTTTGGGG"
		SeqIO.write(record,changedrecord,"fasta")
		SeqIO.write(record,capped_assembly,"fasta")
	elif name in CAmax:
		record.seq="CCCCAAAACCCCAAAACCCCAAAACCCC"+record.seq[CAmax[name]:]
		SeqIO.write(record,changedrecord,"fasta")
		SeqIO.write(record,capped_assembly,"fasta")
	elif name in GTmax:
		record.seq=record.seq[:GTmax[name]]+"GGGGTTTTGGGGTTTTGGGGTTTTGGGG"
		SeqIO.write(record,changedrecord,"fasta")
		SeqIO.write(record,capped_assembly,"fasta")
	else:
		SeqIO.write(record,capped_assembly,"fasta")


assembly.close()
changedrecord.close()
capped_assembly.close()
cappingGT.close()
cappingCA.close()