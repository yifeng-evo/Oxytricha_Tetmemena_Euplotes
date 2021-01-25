f=open("sdrap_ewoo_11032020_pid95_add90_prec_prod_mdsindex_primaryonly.txt","r")
output=open("junction_scramble_sdrap_ewoo_11032020_pid95_add90_prec_prod_mdsindex_primaryonly.txt","w")

def get_index(ls,element):
    i=0
    index_ls=[]
    while i< len(ls):
        if ls[i]==element:
            index_ls.append(i)
        i=i+1
    return index_ls

def combine_adjacent_dup(ls):
    b=[ls[0]]
    i=1
    while i <len(ls):
        if ls[i]==ls[i-1]:
            pass
        else:
            b.append(ls[i])
        i=i+1

    return b


for line in f:
    index_ls=line.split("\t")[3][1:-2]
    digit_index=[int(i) for i in index_ls.split(",")]
    dedup_ls=combine_adjacent_dup(digit_index)
    abs_digit_index=[abs(i) for i in dedup_ls]
    max_index=max(abs_digit_index) # without direction
    if max_index>1: #multimds
        junction_dic={}
        mds_1=1
        mds_2=2
        while mds_2<=max_index:
            junction=str(mds_1)+"_"+str(mds_2)
            pos_1=get_index(abs_digit_index,mds_1)
            pos_2=get_index(abs_digit_index,mds_2)
            ordered=0
            scrambled=0
            for m in pos_1:
                for n in pos_2:
                    if dedup_ls[n]>0 and dedup_ls[m]>0 and n-m==1:#[1,2]
                        ordered=1
                    elif dedup_ls[n]<0 and dedup_ls[m]<0 and m-n==1: #[-2,-1]
                        ordered=1
                    else:
                        scrambled=1
            if ordered==1 and scrambled==1:
                feature="either"
            elif ordered==1 and scrambled==0:
                feature="nonscrambled"
            else:
                feature="scrambled"
            junction_dic[junction]=feature
            mds_1=mds_1+1
            mds_2=mds_2+1
        #if only one mds, don't write
        output.write (line[:-1]+"\t")
        for junction in junction_dic:
            output.write("%s:%s|"% (junction,junction_dic[junction]))
        output.write("\n")
    else:
        output.write(line[:-1]+"\tsingle_mds\n")

f.close()
output.close()
