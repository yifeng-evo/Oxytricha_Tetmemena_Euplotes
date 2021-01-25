
# coding: utf-8

# In[ ]:


from Bio import AlignIO
import itertools
import os.path
from os import path
OG_file=open("./Orthogroups/Orthogroups_SingleCopyOrthologues.txt","r")
#first generate a gap gff


#combine list into ranges


#generate gap gff
def gap_gff(seq):
    i=0
    gap_ls=[]
    def ranges(i):
        for a, b in itertools.groupby(enumerate(i), lambda pair: pair[1] - pair[0]):
            b = list(b)
            yield b[0][1], b[-1][1]+1
    
    while i<len(seq):
        if seq[i]!="-":
            gap_ls.append(i)
	i=i+1
    return list(ranges(gap_ls))
    
            
for line in OG_file:
    OG=line[:-1]
    if os.stat("./Single_Copy_Orthologue_Sequences/"+OG+".pal2nal.out").st_size !=0:
        print ("./Single_Copy_Orthologue_Sequences/"+OG+".pal2nal.out")
        alignment=AlignIO.read("./Single_Copy_Orthologue_Sequences/"+OG+".pal2nal.out","clustal")
        for record in alignment:
            if "evm.model" in record.id:
                ewoo_id=record.id
                ewoo_seq=record.seq
            elif "TET" in record.id:
                tet_id=record.id
                tet_seq=record.seq
            else:
                oxy_id=record.id
                oxy_seq=record.seq
        ewoo_gap_range=gap_gff(ewoo_seq)
        oxy_gap_range=gap_gff(oxy_seq)
        tet_gap_range=gap_gff(tet_seq)
        f=open("./Single_Copy_Orthologue_gap_gff/"+OG+".gap.gff","w")
	f.write("#ewoo\t"+str(len(ewoo_gap_range))+"\n")
	f.write("#tet\t"+str(len(tet_gap_range))+"\n")
	f.write("#oxy\t"+str(len(oxy_gap_range))+"\n")
        for region in ewoo_gap_range:
            f.write(OG+"\tewoo\t"+ewoo_id+"\t"+str(len(ewoo_seq))+"\t"+str(region[0])+"\t"+str(region[1])+"\n")
        for region in oxy_gap_range:
            f.write(OG+"\toxy\t"+oxy_id+"\t"+str(len(oxy_seq))+"\t"+str(region[0])+"\t"+str(region[1])+"\n")
        for region in tet_gap_range:
            f.write(OG+"\ttet\t"+tet_id+"\t"+str(len(tet_seq))+"\t"+str(region[0])+"\t"+str(region[1])+"\n")
            
        f.close()

OG_file.close()

