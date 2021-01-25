import os
f=open("../SDRAP/sdrap_tet_10272020_pid95_add90/sdrap_tet_10272020_pid95_add90_prec_intervals.bed","r")
for line in f:
    if line[0]!="#":
        MIC=line.split("\t")[0]
        MAC_info=line.split("\t")[3]
        MAC=MAC_info.split("_")[2]
        split=open("./tet_split_prec_interval/%s__%s" % (MIC,MAC),"a")
        split.write(line)
        split.close()


directory = os.fsencode("./tet_split_prec_interval")
output=open("sdrap_tet_10272020_pid95_add90_prec_prod_mdsindex_primaryonly.txt","w")
for file in os.listdir(directory):
    filename = os.fsdecode(file)
    MIC=filename.split("__")[0]
    MAC=filename.split("__")[1]
    index_ls=[]
    f=open("./tet_split_prec_interval/"+filename,"r")
    pre_ls=[]
    for line in f:
        info=line.split("\t")[3]
        if info[0:3]=="pre":
            start=int(info.split("_")[3])
            end=int(info.split("_")[4])
            pre_ls.append([start,end])
    f.close()
    f=open("./tet_split_prec_interval/"+filename,"r")
    for line in f:
        info=line.split("\t")[3]
        start=int(info.split("_")[3])
        end=int(info.split("_")[4])
        if [start,end] in pre_ls:
            index=int(info.split("_")[-8])#+"p"
            index_ls.append(index)
        #else:
        #    index=int(info.split("_")[-8])
        #    index_ls.append(index)
    output.write("%s\t%s\t%d\t%s\n" % (MIC,MAC,len(index_ls),str(index_ls)))
    f.close()

output.close()
