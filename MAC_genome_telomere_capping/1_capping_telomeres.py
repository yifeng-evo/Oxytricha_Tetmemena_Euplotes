from Bio import SeqIO
CAtelopoint=open("/ifs/scratch/biochem/ll_lab/yf2342/candidate_CA_telomere_addition_sites_for_1023assembly","r")
GTtelopoint=open("/ifs/scratch/biochem/ll_lab/yf2342/candidate_GT_telomere_addition_sites_for_1023assembly","r")
TASinthemiddle=open("/ifs/scratch/biochem/ll_lab/yf2342/TASinthemiddle_1023assembly.txt","w")
CAtelo={}
GTtelo={}
CAfrequency={}
GTfrequency={}
CAmax={}
GTmax={}

assembly=open("/ifs/scratch/biochem/ll_lab/yf2342/dedupe_98_Euplotes_woodruffi_trinity_spades_cap3_1023.fasta","r")
cappingCA=open("/ifs/scratch/biochem/ll_lab/yf2342/CAcapping_1023assembly.txt","w")
cappingGT=open("/ifs/scratch/biochem/ll_lab/yf2342/GTcapping_1023assembly.txt","w")

contiglengh={}
for record in SeqIO.parse(assembly,"fasta"):
	contiglengh[record.id]=len(record.seq)
	CAtelo[record.id]={}
	GTtelo[record.id]={}


for line in CAtelopoint:
	contigname=line.split("\t")[0]
	TAS=int(line.split("\t")[1])
	frequency=int(line.split("\t")[2])
	if TAS <100: # only check TAS at ends
		CAtelo[contigname][TAS]=frequency
		if contigname in CAfrequency:
			CAfrequency[contigname].append(frequency)
		else:
			CAfrequency[contigname]=[frequency]
	else:#record telo in the middle:
		TASinthemiddle.write(line[:-1]+"\tCA\n")


for contigname in CAtelo:
	for TAS in CAtelo[contigname]:
		if max(CAfrequency[contigname])==CAtelo[contigname][TAS]:
			if contigname in CAmax: #if two TAS tied 
				if CAmax[contigname]>=TAS: #leave the smaller one as the authentic telomere addition site
					CAmax[contigname]=TAS
			else:
				CAmax[contigname]=TAS


for line in GTtelopoint:
	contigname=line.split("\t")[0]
	TAS=int(line.split("\t")[1])
	frequency=int(line.split("\t")[2])
	if TAS > contiglengh[contigname]-100: # only check TAS at ends
		GTtelo[contigname][TAS]=frequency
		if contigname in GTfrequency:
			GTfrequency[contigname].append(frequency)
		else:
			GTfrequency[contigname]=[frequency]
	else:#record telo in the middle:
		TASinthemiddle.write(line[:-1]+"\tGT\n")


for contigname in GTtelo:
	for TAS in GTtelo[contigname]:
		if max(GTfrequency[contigname])==GTtelo[contigname][TAS]:
			if contigname in GTmax: #if two TAS tied 
				if GTmax[contigname]<=TAS: #leave the smaller one as the authentic telomere addition site
					GTmax[contigname]=TAS
			else:
				GTmax[contigname]=TAS

for key in CAmax:
	cappingCA.write(key+"\t"+str(CAmax[key])+"\n")
for key in GTmax:
	cappingGT.write(key+"\t"+str(GTmax[key])+"\n")

	
CAtelopoint.close()
GTtelopoint.close()
assembly.close()
cappingCA.close()
cappingGT.close()








	




