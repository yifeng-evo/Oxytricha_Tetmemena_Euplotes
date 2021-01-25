from Bio import SeqIO
f=open("/ifs/scratch/biochem/ll_lab/yf2342/Ewoo_MAC_telo_reads_to_assembly_1023.psl","r")
f1=open("/ifs/scratch/biochem/ll_lab/yf2342/Euplotes_woo_raw_reads/corrected/telo_reads_renamed_trimmed_cutadapt.corrected.fasta","r")
CAtelopoint=open("/ifs/scratch/biochem/ll_lab/yf2342/candidate_CA_telomere_addition_sites_for_1023assembly","w")
GTtelopoint=open("/ifs/scratch/biochem/ll_lab/yf2342/candidate_GT_telomere_addition_sites_for_1023assembly","w")

GTreads=[]
CAreads=[]
CAtelo={}
GTtelo={}
contiglength={}


for record in SeqIO.parse(f1,"fasta"):
	if "GGGGTTTTGGGG" in record.seq:
		GTreads.append(record.id)
	if "CCCCAAAACCCC" in record.seq:
		CAreads.append(record.id)


f.readline()
f.readline()
f.readline()
f.readline()
f.readline()
while True:
	line=f.readline()
	if len(line)<3:
		break
	mismatches=float(line.split("\t")[1])
	qsize=float(line.split("\t")[10])
	qname=line.split("\t")[9]
	qstart=int(line.split("\t")[11])
	qend=int(line.split("\t")[12])
	tname=line.split("\t")[13]
	tstart=int(line.split("\t")[15])
	tend=int(line.split("\t")[16])
	if mismatches < 0.05*qsize: # only look at matches over 95% identical
		if line.split("\t")[8]=="+": #plus strand
			if qname in CAreads: # contain ccccaaaa
				if qstart == 0: 
					if tname not in CAtelo:
						CAtelo[tname]={}
					if tstart in CAtelo[tname]:
						CAtelo[tname][tstart]=CAtelo[tname][tstart]+1
					else:
						CAtelo[tname]={tstart:1}

			if qname in GTreads:
				if qend == qsize:
					if tname not in GTtelo:
						GTtelo[tname]={}
					if tend in GTtelo[tname]:
						GTtelo[tname][tend]=GTtelo[tname][tend]+1
					else:
						GTtelo[tname]={tend:1}
						
		elif line.split("\t")[8]=="-":
			if qname in CAreads: # contain ccccaaaa
				if qend == qsize:
					if tname not in GTtelo:
						GTtelo[tname]={}
					if tend in GTtelo[tname]:
						GTtelo[tname][tend]=GTtelo[tname][tend]+1
					else:
						GTtelo[tname]={tend:1}
									
				
			if qname in GTreads:
				if qstart == 0: 
					if tname not in CAtelo:
						CAtelo[tname]={}
					if tstart in CAtelo[tname]:
						CAtelo[tname][tstart]=CAtelo[tname][tstart]+1
					else:
						CAtelo[tname]={tstart:1}
						

print("start printing")

for contig in CAtelo:
	for i in CAtelo[contig]:
		if CAtelo[contig][i]>=5:
			CAtelopoint.write(contig+"\t"+str(i)+"\t"+str(CAtelo[contig][i])+"\n")

for contig in GTtelo:
	for i in GTtelo[contig]:
		if GTtelo[contig][i]>=5:
			GTtelopoint.write(contig+"\t"+str(i)+"\t"+str(GTtelo[contig][i])+"\n")

f.close()
f1.close()
CAtelopoint.close()
GTtelopoint.close()










