f=open("../SDRAP/sdrap_oxy_1025_pid95_add90/sdrap_oxy_1025_pid95_add90_prec_intervals.bed","r")
for line in f:
    if line[0]!="#":
        MIC=line.split("\t")[0]
        MAC_info=line.split("\t")[3]
        MAC=MAC_info.split("_")[2]
        split=open("./oxy_split_prec_interval/%s__%s" % (MIC,MAC),"a")
        split.write(line)
        split.close()
    
import os
directory = os.fsencode("./oxy_split_prec_interval")
output=open("sdrap_oxy_1025_pid95_add90_prec_prod_mdsindex.txt","w")
for file in os.listdir(directory):
    filename = os.fsdecode(file)
    MIC=filename.split("__")[0]
    MAC=filename.split("__")[1]
    index_ls=[]
    f=open("./oxy_split_prec_interval/"+filename,"r")
    for line in f:
        info=line.split("\t")[3]
        index=int(info.split("_")[-8])
        index_ls.append(index)
    output.write("%s\t%s\t%d\t%s\n" % (MIC,MAC,len(index_ls),str(index_ls)))
    f.close()
    

output.close()
