


f=open("./sdrap_ewoo_11032020_pid95_add90/sdrap_ewoo_11032020_pid95_add90_properties.tsv","r")
MDS=open("./sdrap_ewoo_11032020_pid95_add90/sdrap_ewoo_11032020_pid95_add90_prod_intervals.bed","r")
pri_MDS=open("./sdrap_ewoo_11032020_pid95_add90/Primary_sdrap_ewoo_11032020_pid95_add90_prod_intervals.bed","w")

pair={}
for line in f:
    if line[0:4]!="prec":
        MIC=line.split("\t")[0]
        MAC=line.split("\t")[1]
        cov=float(line.split("\t")[2])
        if MAC not in pair:
            pair[MAC]=[MIC,cov]
        else:
            if cov>pair[MAC][1]:
                pair[MAC]=[MIC,cov]
            
for line in MDS:
    if line[0]!="#":
        MAC=line.split("\t")[0]
        MIC_info=line.split("\t")[3]
        MIC_name_ls=MIC_info.split("_")[2:-10]
        MIC="_".join(MIC_name_ls)
        if MAC in pair and pair[MAC][0]==MIC:
            pri_MDS.write(line)
        
f.close()
MDS.close()
pri_MDS.close()
    


# In[19]:


import matplotlib.pyplot as plt
import statistics
tet_mds=open("./sdrap_tet_10272020_pid95_add90/Primary_sdrap_tet_10272020_pid95_add90_prod_intervals.bed","r")
oxy_mds=open("./sdrap_oxy_1025_pid95_add90/Primary_sdrap_oxy_1025_pid95_add90_prod_intervals.bed","r")
ewoo_mds=open("./sdrap_ewoo_11032020_pid95_add90/Primary_sdrap_ewoo_11032020_pid95_add90_prod_intervals.bed","r")

def MDS_len_ls (file):
    mds_ls=[]
    for line in file:
        mds_s=int(line.split("\t")[1])
        mds_e=int(line.split("\t")[2])
        mds_len=mds_e-mds_s+1
        mds_ls.append(mds_len)
    return mds_ls

def print_info (name,ls):
    print (name+"_MDS:%d_median:%d_average:%d_min:%d_max:%d" % (len(ls),statistics.median(ls),sum(ls)/len(ls),min(ls),max(ls)))
    

tet_ls=MDS_len_ls(tet_mds)
oxy_ls=MDS_len_ls(oxy_mds)
ewoo_ls=MDS_len_ls(ewoo_mds)

plt.style.use('default')    
plt.hist([tet_ls,oxy_ls,ewoo_ls],bins=50, density=True, range=(0,1501),color=["orange","gray","#00b9f1"],label=["Tetmemena","Oxytricha","Euplotes"])
plt.legend(loc="best")
#plt.title("MDS length distribution")
plt.ylabel("Frequency")
plt.xlabel("MDS length (bp)")
plt.savefig('./Figure1/3species_MDS_distribution.png', format='png', bbox_inches='tight',dpi=300)
print_info("tet",tet_ls)
print_info("oxy",oxy_ls)
print_info("ewoo",ewoo_ls)

plt.show()
tet_mds.close()
oxy_mds.close()
ewoo_mds.close()


# In[17]:


f=open("junction_scramble_sdrap_ewoo_11032020_pid95_add90_prec_prod_mdsindex_primaryonly.txt","r")
pointer=open("./sdrap_ewoo_11032020_pid95_add90/sdrap_ewoo_11032020_pid95_add90_prod_pointers.bed","r")
output=open("Scramble_sdrap_ewoo_11032020_pid95_add90_prod_pointers.bed","w")

junction_dic={}
for line in f:
    if "single_mds" not in line:
        MIC=line.split("\t")[0]
        MAC=line.split("\t")[1]
        junction=line.split("\t")[4][:-2]#remove |\n in the end
        if MAC in junction_dic:
            junction_dic[MAC][MIC]={}
        else:
            junction_dic[MAC]={}
            junction_dic[MAC][MIC]={}
        for info in junction.split("|"):
            junc=info.split(":")[0]
            feature=info.split(":")[1]
            junction_dic[MAC][MIC][junc]=feature
            
for line in pointer:
    if line[0]!="#":
        MAC=line.split("\t")[0]
        info=line.split("\t")[3]
        MIC_ls=info.split("_")[1:-11]
        MIC=("_").join(MIC_ls)
        mds_1=int(info.split("_")[-11])
        mds_2=int(info.split("_")[-8])
        junc=str(min(mds_1,mds_2))+"_"+str(max(mds_1,mds_2))
        output.write(line[:-1]+"\t"+junction_dic[MAC][MIC][junc]+"\n")
        
f.close()
pointer.close()
output.close()
            
        
    


# In[18]:


f=open("Scramble_sdrap_ewoo_11032020_pid95_add90_prod_pointers.bed","r")
output=open("master_Scramble_sdrap_ewoo_11032020_pid95_add90_prod_pointers.bed","w")

pointer_dic={}
for line in f:
    MAC=line.split("\t")[0]
    start=line.split("\t")[1]
    end=line.split("\t")[2]
    feature=line.split("\t")[-1][:-1] #remove \n
    if MAC+"\t"+start+"\t"+end not in pointer_dic:
        pointer_dic[MAC+"\t"+start+"\t"+end]=feature
    else:
        if feature!=pointer_dic[MAC+"\t"+start+"\t"+end]:
            pointer_dic[MAC+"\t"+start+"\t"+end]="either"
            
for pointer in pointer_dic:
    output.write(pointer+"\t"+pointer_dic[pointer]+"\n")
    
output.close()
f.close()
    
                
            
            
    
    


# In[7]:


f=open("master_Scramble_sdrap_tet_10272020_pid95_add90_prod_pointers.bed","r")
output=open("master_Scramble_sdrap_tet_10272020_pid95_add90_MAC_contigs.txt","w")

contig_dic={}
for line in f:
    MAC=line.split("\t")[0]
    feature=line.split("\t")[-1][:-1]
    if MAC not in contig_dic:
        contig_dic[MAC]=[feature]
    else:
        contig_dic[MAC].append(feature)
a=0
b=0
c=0
d=0
e=0
for contig in contig_dic:
    contig_set=sorted(set(contig_dic[contig]))
    if contig_set ==["nonscrambled"]:
        a=a+1
        output.write(contig+"\t"+"All_nonscrambled\n")
    elif contig_set ==["scrambled"]:
        b=b+1
        output.write(contig+"\t"+"All_scrambled\n")
    elif contig_set==["either","nonscrambled"]:
        c=c+1
        output.write(contig+"\t"+"low_confidence_nonscrambled\n")
    elif contig_set==["either"]:
        d=d+1
        output.write(contig+"\t"+"low_confidence_nonscrambled\n")
    else:
        e=e+1
        output.write(contig+"\t"+"Scrambled\n")
        
print (a,b,c,d,e)
f.close()
output.close()
    


# In[18]:


import matplotlib.pyplot as plt
import statistics
f=open("master_Scramble_sdrap_tet_10272020_pid95_add90_prod_pointers.bed","r")
nonscrambled_ls=[]
scrambled_ls=[]
either_ls=[]
for line in f:
    MAC=line.split("\t")[0]
    start=int(line.split("\t")[1])
    end=int(line.split("\t")[2])
    length=end-start+1
    feature=line.split("\t")[3][:-1]
    if length <=30:
        if feature=="nonscrambled":
            nonscrambled_ls.append(length)
        elif feature=="either":
            either_ls.append(length)
        else:
            scrambled_ls.append(length)
           
f.close()

plt.style.use('default')    
plt.hist([nonscrambled_ls,scrambled_ls],bins=30, density=True, range=(0,30.1),color=["blue","red"],label=["nonscrambled","scrambled"])
plt.legend(loc="best")
#plt.title("Euplotes pointer length distribution")
plt.ylabel("Frequency")
plt.xlabel("pointer length (bp)")
plt.savefig('./Figure1/Tetmemena_pointer_distribution_below30bp.pdf', format='pdf', bbox_inches='tight')
print ("nonscrambled:%d_median:%d_average:%d" % (len(nonscrambled_ls),statistics.median(nonscrambled_ls),sum(nonscrambled_ls)/len(nonscrambled_ls)))
print ("scrambled:%d_median:%d_average:%d" % (len(scrambled_ls),statistics.median(scrambled_ls),sum(scrambled_ls)/len(scrambled_ls)))
plt.show()    
    



    

