
# coding: utf-8

# In[52]:


#bout should be sepearated bout for 3 orfs
bout=open("57kd_oxyMIC.bout","r")
output=open("TBE_57kd_oxyMIC.txt","w") #1based position
MIC_dict={}

import itertools

def convert(i):
    for a, b in itertools.groupby(enumerate(i), lambda pair: pair[1] - pair[0]):
        b = list(b)
        yield [b[0][1], b[-1][1]]



for line in bout:
    MIC=line.split("\t")[1]
    sstart=int(line.split("\t")[8])
    send=int(line.split("\t")[9])
        
    real_start=min(send,sstart)-1 # invert when 
    real_end=max(sstart,send)

    if MIC not in MIC_dict:
        MIC_dict[MIC]=list(range(real_start,real_end))
    else:
        MIC_dict[MIC].extend(list(range(real_start,real_end)))



for MIC in MIC_dict:
    tot_range=list(convert(sorted(set(MIC_dict[MIC]))))
    if len(tot_range)==1:
        output.write(MIC+"\t%d\t%d\n" % (tot_range[0][0]+1,tot_range[0][1]+1))
    else:
        block=[]
        start=tot_range[0][0]
        end=tot_range[0][1]+1 # initiation	itertools right close!!!!!
        block=[[start,end]]
        i=1

        while i <len(tot_range):	
            s=tot_range[i][0]
            e=tot_range[i][1]+1
            for region in block:
                v=0
                start=region[0]
                end=region[1]
                if abs(s-end)<=200:
                    region[1]=e
                    v=1
                    break
                else:
                    v=0
            if v==0:
                block.append([s,e])	
            i=i+1

        for k in block:
            output.write(MIC+"\t%d\t%d\n" % (k[0]+1,k[1]))


bout.close()
output.close()


# In[55]:


f=open("TBE_42kd_oxyMIC.txt","r")
bout=open("42kd_oxyMIC.bout","r")
output=open("TBE_42kd_oxyMIC_with_orientation.txt","w")
abnormal=open("TBE_42kd_oxyMIC_opposite_orientation.abnormal.txt","w")


dic={}
for line in f:
    rname=line.split("\t")[0]
    start=int(line.split("\t")[1])
    end=int(line.split("\t")[2])
    if rname not in dic:
        dic[rname]=[[start,end]]
    else:
        dic[rname].append([start,end])
        
for line in bout:
    MIC=line.split("\t")[1]
    #qstart=int(line.split("\t")[6])
    #qend=int(line.split("\t")[7])
    sstart=int(line.split("\t")[8])
    send=int(line.split("\t")[9])
    real_start=min(send,sstart) # invert when 
    real_end=max(sstart,send)
    frame=int(line.split("\t")[15])
    for region in dic[MIC]:
        if real_start >=region[0] and real_end<=region[1]:
            region=region.append(frame)
            break
         

           
f.close()
bout.close()
for MIC in dic:
    for region in dic[MIC]:
        f_count=len(region)-2
        if region[2]>0:
            initial=1
            ori="+"
        else:
            initial=-1
            ori="-"
        i=3
        while i<len(region):
            if region[i]>0:
                k=1
            else:
                k=-1
            if initial * k==-1:
                abnormal.write(MIC+"\t"+str(region)+"\n")
            i=i+1
        
        output.write(MIC+"\t"+str(region[0])+"\t"+str(region[1])+"\t"+str(f_count)+"\t"+ori+"\n")
                     

output.close()
abnormal.close()
    
    


# In[57]:


f42kd=open("TBE_42kd_oxyMIC_with_orientation.txt","r")
f22kd=open("TBE_22kd_oxyMIC_with_orientation.txt","r")
f57kd=open("TBE_57kd_oxyMIC_with_orientation.txt","r")
output=open("TBE_oxy_MIC_complete.txt","w")
#output2=open("TBE_oxy_MIC_complete_only42_22.txt","w")
linker1=open("42kd_22kd_linker_oxy.txt","w")
linker2=open("22kd_57kd_linker_oxy.txt","w")
#output3=open("OXY_other22kd_inregion.txt","w")
#output4=open("OXY_other57kd_inregion.txt","w")
windowsize=2000
lower_window=-100

dict_42={}
for line in f42kd:
    MIC=line.split("\t")[0]
    start=int(line.split("\t")[1])
    end=int(line.split("\t")[2])
    ori=line.split("\t")[4][:-1]
    if MIC not in dict_42:
        dict_42[MIC]=[["42kd",start,end,ori]]
    else:
        dict_42[MIC].append(["42kd",start,end,ori])

dict_22={}
for line in f22kd:
    MIC=line.split("\t")[0]
    start_22=int(line.split("\t")[1])
    end_22=int(line.split("\t")[2])
    ori_22=line.split("\t")[4][:-1]
    if MIC in dict_42:
        for region in dict_42[MIC]:
            start_42=region[1]
            end_42=region[2]
            ori_42=region[3]
            if lower_window<=start_22-end_42<=windowsize and ori_42=="+" and ori_22=="-":
#                 if len(region)==8:
#                     output3.write(line)
#                 else:
                linker1.write(str(start_22-end_42)+"\n")
                region.extend(["22kd",start_22,end_22,ori_22])
#                     if MIC not in dict_22:
#                         dict_22[MIC]=[region]
#                     else:
#                         dict_22[MIC].append(region)
            elif lower_window<=start_42-end_22<=windowsize and ori_42=="-" and ori_22=="+":
#                 if len(region)==8:
#                     output3.write(line)
#                 else:
                linker1.write(str(start_42-end_22)+"\n")
                region.extend(["22kd",start_22,end_22,ori_22])
#                     if MIC not in dict_22:
#                         dict_22[MIC]=[region]
#                     else:
#                         dict_22[MIC].append(region)

#don'breank

# for MIC in dict_22:
#     for region in dict_22[MIC]:
#         output2.write(MIC+"\t"+str(region)+"\n")

for line in f57kd:
    MIC=line.split("\t")[0]
    start_57=int(line.split("\t")[1])
    end_57=int(line.split("\t")[2])
    ori_57=line.split("\t")[4][:-1]
    if MIC in dict_42:
        for region in dict_42[MIC]:
            if "42kd" in region and "22kd" in region:
                start_22=region[5]
                end_22=region[6]
                ori_22=region[7]
                if lower_window<=start_57-end_22 <=windowsize and ori_57=="+" and ori_22=="-":
#                     if len(region)==12:
#                         output4.write(line)
#                     else:
                    linker2.write(str(start_57-end_22)+"\n")
                    region.extend(["57kd",start_57,end_57,ori_57])

                elif lower_window<=start_22-end_57 <=windowsize and ori_57=="-" and ori_22=="+":
#                     if len(region)==9:
#                         output4.write(line)
#                     else:
                    linker2.write(str(start_22-end_57)+"\n")
                    region.extend(["57kd",start_57,end_57,ori_57])
                
for MIC in dict_42: 
    for region in dict_42[MIC]:
        if "57kd" in region:
            i=0
            while i <len(region):
                if region[i]=="42kd":
                    output.write(MIC+"\t42kd_%d_%d_%s" % (region[i+1],region[i+2],region[i+3]))
                elif region [i]=="22kd":
                    output.write("\t22kd_%d_%d_%s" % (region[i+1],region[i+2],region[i+3]))
                else:
                    output.write("\t57kd_%d_%d_%s" % (region[i+1],region[i+2],region[i+3]))
                i=i+4 
            output.write("\n")
                
                    
output.close()
f42kd.close()
f22kd.close()
f57kd.close()
#output2.close()
linker1.close()
linker2.close()
# output3.close()
# output4.close()


# In[22]:


f=open("sorted_TBE_oxy_MIC_complete.txt","r") #####check tail!!!!!
output=open("final_TBE_oxy_MIC_complete.txt","w")
output2=open("complete_TBE_orf_oxy_MIC.txt","w")

#MIC_pre=""
overlap=[]
list_complete=[]

for line in f:
    line=line[:-1]
    MIC=line.split("\t")[0]
    for region in line.split("\t")[1:]:
        if MIC+"_"+region in list_complete:
            overlap.append(MIC+"_"+region)
        else:
            list_complete.append(MIC+"_"+region)
        
f.close()
f=open("sorted_TBE_oxy_MIC_complete.txt","r")

overlap_line=[]
for line in f:
    line=line[:-1]
    MIC=line.split("\t")[0]
    x=0
    for region in line.split("\t")[1:]:
        v=0
        if MIC+"_"+region in overlap: #overlap region
            x=1
            for group in overlap_line:
                if line in group:
                    v=1
                    break
                else:
                    for line1 in group:
                        if region in line1:
                            group = group.append(line)
                            v=1
                            break
                    if v==1:
                        break
                       
            if v==0:
                overlap_line.append([line])
            break
            
    if x==0:
        output.write(line+"\n")
        
             
                
f.close()
for region in list_complete:
    output2.write(region+"\n")

output2.close()

for group in overlap_line:
    orf=[]
    for line in group:
        orf=orf+line.split("\t")[1:]
    orf=list(set(orf))
    orf.sort()
    output.write(line.split("\t")[0]+"\t"+("\t").join(orf)+"\n")
    
output.close()   


# In[27]:


str22="TBE_22kd_oxyMIC_with_orientation.txt"
f22=open(str22,"r")
f42=open("TBE_42kd_oxyMIC_with_orientation.txt","r")
f57=open("TBE_57kd_oxyMIC_with_orientation.txt","r")

f_42_22=open("TBE_oxy_MIC_42_22.txt","w")
f_22_57=open("TBE_oxy_MIC_22_57.txt","w")
windowsize=2000
lower_window=-100

dict_42={}
for line in f42:
    MIC=line.split("\t")[0]
    start=int(line.split("\t")[1])
    end=int(line.split("\t")[2])
    ori=line.split("\t")[4][:-1]
    if MIC not in dict_42:
        dict_42[MIC]=[["42kd",start,end,ori]]
    else:
        dict_42[MIC].append(["42kd",start,end,ori])

f42.close()
for line in f22:
    MIC=line.split("\t")[0]
    start_22=int(line.split("\t")[1])
    end_22=int(line.split("\t")[2])
    ori_22=line.split("\t")[4][:-1]
    if MIC in dict_42:
        for region in dict_42[MIC]:
            start_42=region[1]
            end_42=region[2]
            ori_42=region[3]
            if lower_window<=start_22-end_42<=windowsize and ori_42=="+" and ori_22=="-":
                region.extend(["22kd",start_22,end_22,ori_22])
            elif lower_window<=start_42-end_22<=windowsize and ori_42=="-" and ori_22=="+":
                region.extend(["22kd",start_22,end_22,ori_22])
f22.close()


for MIC in dict_42: 
    for region in dict_42[MIC]:
        if "22kd" in region:
            i=0
            while i <len(region):
                if region[i]=="42kd":
                    f_42_22.write(MIC+"\t42kd_%d_%d_%s" % (region[i+1],region[i+2],region[i+3]))
                elif region [i]=="22kd":
                    f_42_22.write("\t22kd_%d_%d_%s" % (region[i+1],region[i+2],region[i+3]))
                i=i+4 
            f_42_22.write("\n")
                
                
f22=open(str22,"r")
dict_22={}
for line in f22:
    MIC=line.split("\t")[0]
    start=int(line.split("\t")[1])
    end=int(line.split("\t")[2])
    ori=line.split("\t")[4][:-1]
    if MIC not in dict_22:
        dict_22[MIC]=[["22kd",start,end,ori]]
    else:
        dict_22[MIC].append(["22kd",start,end,ori])  
        
for line in f57:
    MIC=line.split("\t")[0]
    start_57=int(line.split("\t")[1])
    end_57=int(line.split("\t")[2])
    ori_57=line.split("\t")[4][:-1]
    if MIC in dict_22:
        for region in dict_22[MIC]:
            start_22=region[1]
            end_22=region[2]
            ori_22=region[3]
            if lower_window<=start_57-end_22<=windowsize and ori_57=="+" and ori_22=="-":
                region.extend(["57kd",start_57,end_57,ori_57])
            elif lower_window<=start_22-end_57<=windowsize and ori_57=="-" and ori_22=="+":
                region.extend(["57kd",start_57,end_57,ori_57])
                

for MIC in dict_22: 
    for region in dict_22[MIC]:
        if "57kd" in region:
            i=0
            while i <len(region):
                if region[i]=="22kd":
                    f_22_57.write(MIC+"\t22kd_%d_%d_%s" % (region[i+1],region[i+2],region[i+3]))
                elif region [i]=="57kd":
                    f_22_57.write("\t57kd_%d_%d_%s" % (region[i+1],region[i+2],region[i+3]))
                i=i+4 
            f_22_57.write("\n")


f42.close()
f22.close()
f57.close()
f_42_22.close()
f_22_57.close()
    


# In[24]:


f=open("sorted_TBE_oxy_MIC_42_22.txt","r") #####check tail!!!!!
output=open("final_TBE_oxy_MIC_42_22.txt","w")
comp_orf=open("complete_TBE_orf_oxy_MIC.txt","r")
output2=open("only42_22_TBE_orf_oxy_MIC.txt","w")

comp_list=[]
for line in comp_orf:
    comp_list.append(line[:-1])

list_42_22=[]
overlap=[]
for line in f:
    line=line[:-1]
    MIC=line.split("\t")[0]
    for region in line.split("\t")[1:]:
        if MIC+"_"+region in comp_list:
            pass
        else:
            if MIC+"_"+region in list_42_22:
                overlap.append(MIC+"_"+region)
            else:
                list_42_22.append(MIC+"_"+region)

f.close()

f=open("sorted_TBE_oxy_MIC_42_22.txt","r")             
overlap_line=[]
for line in f:
    line=line[:-1]###very important!!!
    MIC=line.split("\t")[0]
    x=0
    for region in line.split("\t")[1:]:
        v=0
        if str(MIC+"_"+region) in overlap: #overlap region
            x=1
            for group in overlap_line:
                if line in group:
                    v=1
                    break
                else:
                    for line1 in group:
                        if region in line1:
                            group = group.append(line)
                            v=1     
                            break
                        if v==1:
                            break
            if v==0:
                overlap_line.append([line])
            break
            
    if x==0:
        output.write(line+"\n")
        
             
                
f.close()
for region in list_42_22:
    output2.write(region+"\n")

output2.close()

for group in overlap_line:
    orf=[]
    for line in group:
        orf=orf+line.split("\t")[1:]
    orf=list(set(orf))
    orf.sort()
    output.write(line.split("\t")[0]+"\t"+("\t").join(orf)+"\n")
    
output.close()   
comp_orf.close()
    


# In[29]:


f=open("sorted_TBE_oxy_MIC_22_57.txt","r") #####check tail!!!!!
output=open("final_TBE_oxy_MIC_22_57.txt","w")
comp_orf=open("complete_TBE_orf_oxy_MIC.txt","r")
output2=open("only22_57_TBE_orf_oxy_MIC.txt","w")

comp_list=[]
for line in comp_orf:
    comp_list.append(line[:-1])

list_42_22=[]
overlap=[]
for line in f:
    line=line[:-1]
    MIC=line.split("\t")[0]
    for region in line.split("\t")[1:]:
        if MIC+"_"+region in comp_list:
            pass
        else:
            if MIC+"_"+region in list_42_22:
                overlap.append(MIC+"_"+region)
            else:
                list_42_22.append(MIC+"_"+region)

f.close()

f=open("sorted_TBE_oxy_MIC_22_57.txt","r")             
overlap_line=[]
for line in f:
    line=line[:-1]###very important!!!
    MIC=line.split("\t")[0]
    x=0
    for region in line.split("\t")[1:]:
        v=0
        if str(MIC+"_"+region) in overlap: #overlap region
            x=1
            for group in overlap_line:
                if line in group:
                    v=1
                    break
                else:
                    for line1 in group:
                        if region in line1:
                            group = group.append(line)
                            v=1     
                            break
                        if v==1:
                            break
            if v==0:
                overlap_line.append([line])
            break
            
    if x==0:
        output.write(line+"\n")
        
             
                
f.close()
for region in list_42_22:
    output2.write(region+"\n")

output2.close()

for group in overlap_line:
    orf=[]
    for line in group:
        orf=orf+line.split("\t")[1:]
    orf=list(set(orf))
    orf.sort()
    output.write(line.split("\t")[0]+"\t"+("\t").join(orf)+"\n")
    
output.close()   
comp_orf.close()
    
    


# In[11]:


from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO
bout=open("22kd_tetMIC.bout","r")
output=open("TetMIC_22kDa.faa","w")
for line in bout:
    MIC=line.split("\t")[1]
    qstart=int(line.split("\t")[6])
    qend=int(line.split("\t")[7])
    if (qend-qstart+1)>=0.9*193:
        sstart=int(line.split("\t")[8])
        send=int(line.split("\t")[9])
        sseq=line.split("\t")[13]
        record=SeqRecord(Seq(sseq,IUPAC.protein),id=MIC+"_%d_%d" % (sstart,send))
        SeqIO.write(record,output,"fasta")
        
bout.close()
output.close()
        
        
    


# In[6]:


from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO
bout=open("57kd_tetMIC.bout","r")
output=open("TetMIC_57kDa.faa","w")
for line in bout:
    MIC=line.split("\t")[1]
    qstart=int(line.split("\t")[6])
    qend=int(line.split("\t")[7])
    if (qend-qstart+1)>=0.7*482:
        sstart=int(line.split("\t")[8])
        send=int(line.split("\t")[9])
        sseq=line.split("\t")[13]
        record=SeqRecord(Seq(sseq,IUPAC.protein),id=MIC+"_%d_%d" % (sstart,send))
        SeqIO.write(record,output,"fasta")
        
bout.close()
output.close()


# In[10]:


from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO
bout=open("42kd_tetMIC.bout","r")
output=open("TetMIC_42kDa.faa","w")
for line in bout:
    MIC=line.split("\t")[1]
    qstart=int(line.split("\t")[6])
    qend=int(line.split("\t")[7])
    if (qend-qstart+1)>=0.9*354:
        sstart=int(line.split("\t")[8])
        send=int(line.split("\t")[9])
        sseq=line.split("\t")[13]
        record=SeqRecord(Seq(sseq,IUPAC.protein),id=MIC+"_%d_%d" % (sstart,send))
        SeqIO.write(record,output,"fasta")
        
bout.close()
output.close()


# In[1]:


from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO
bout=open("57kd_to_oxytri_MIC.bout","r")
output=open("OxyMIC_57kDa.faa","w")
for line in bout:
    MIC=line.split("\t")[1]
    qstart=int(line.split("\t")[6])
    qend=int(line.split("\t")[7])
    if (qend-qstart+1)>=0.9*482:
        sstart=int(line.split("\t")[8])
        send=int(line.split("\t")[9])
        sseq=line.split("\t")[13]
        record=SeqRecord(Seq(sseq,IUPAC.protein),id=MIC+"_%d_%d" % (sstart,send))
        SeqIO.write(record,output,"fasta")
        
bout.close()
output.close()



# In[34]:


import re
f42kd=open("TBE_42kd_oxyMIC_with_orientation.txt","r")
f22kd=open("TBE_22kd_oxyMIC_with_orientation.txt","r")
f57kd=open("TBE_57kd_oxyMIC_with_orientation.txt","r")
comp=open("final_TBE_oxy_MIC_complete.txt","r")
f22_57=open("final_TBE_oxy_MIC_22_57.txt","r")
f42_22=open("final_TBE_oxy_MIC_42_22.txt","r")
orflist1=open("complete_TBE_orf_oxy_MIC.txt","r")
orflist2=open("only42_22_TBE_orf_oxy_MIC.txt","r")
orflist3=open("only22_57_TBE_orf_oxy_MIC.txt","r")
output=open("TBE_oxy_region_merged.txt","w")


def get_region(line):
    ls=[]
    for orf in line.split("\t")[1:]:
        s=int(orf.split("_")[1])
        e=int(orf.split("_")[2])
        ls.append(s)
        ls.append(e)
    return min(ls),max(ls)

def get_orflist(f):
    ls=[]
    for line in f:
        ls.append(line[:-1]) # remove\n
    return ls
        
comp_orf=get_orflist(orflist1)
orf42_22=get_orflist(orflist2)
orf22_57=get_orflist(orflist3)

for line in comp:
    MIC=line.split("\t")[0]
    start,end=get_region(line)
    output.write(MIC+"\t%d\t%d\tcomplete\n" % (start,end))
               
for line in f42_22:
    line=line[:-1]
    MIC=line.split("\t")[0]
    for region in line.split("\t")[1:]:
        if MIC+"_"+region in orf42_22:
            start,end=get_region(line)
            output.write(MIC+"\t%d\t%d\t42_22\n" % (start,end))
            break

                
for line in f22_57:
    line=line[:-1]
    MIC=line.split("\t")[0]
    for region in line.split("\t")[1:]:
        if MIC+"_"+region in orf22_57:
            start,end=get_region(line)
            output.write(MIC+"\t%d\t%d\t22_57\n" % (start,end))
            break
    

for line in f42kd:
    MIC=line.split("\t")[0] 
    start=line.split("\t")[1]
    end=line.split("\t")[2]
    ori=line.split("\t")[4][:-1]
    orf=MIC+"_42kd_"+start+"_"+end+"_"+ori
    if orf not in comp_orf and orf not in orf42_22:
        output.write(MIC+"\t%s\t%s\t42kd\n" % (start,end))


for line in f22kd:
    MIC=line.split("\t")[0] 
    start=line.split("\t")[1]
    end=line.split("\t")[2]
    ori=line.split("\t")[4][:-1]
    orf=MIC+"_22kd_"+start+"_"+end+"_"+ori
    if orf not in comp_orf and orf not in orf42_22 and orf not in orf22_57:
        output.write(MIC+"\t%s\t%s\t22kd\n" % (start,end))

for line in f57kd:
    MIC=line.split("\t")[0] 
    start=line.split("\t")[1]
    end=line.split("\t")[2]
    ori=line.split("\t")[4][:-1]
    orf=MIC+"_57kd_"+start+"_"+end+"_"+ori
    if orf not in comp_orf and orf not in orf22_57:
        output.write(MIC+"\t%s\t%s\t57kd\n" % (start,end))

        
f42kd.close()
f22kd.close()
f57kd.close()
comp.close()
f22_57.close()
f42_22.close()
output.close()
orflist1.close()
orflist2.close()
orflist3.close()




