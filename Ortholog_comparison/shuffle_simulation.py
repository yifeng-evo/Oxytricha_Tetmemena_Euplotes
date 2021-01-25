import random
import os.path
from os import path
import sys
#seed=int(sys.argv[1])
OG_file=open("./Orthogroups/Orthogroups_SingleCopyOrthologues_3species.txt","r")
#random.seed(seed)

def convert2aln(pointer,spe_rl):
    if pointer!=-1:
        head=0
        tail=0
        for gap in spe_rl:
            gap_len=abs(gap[1]-gap[0])
            if head+gap_len<pointer:
                head=head+gap_len
            else:
                tail=pointer-head
                break
        pointer_aln=gap[0]+tail
        return pointer_aln
    else:
        return -1


def circle_shuffle(random_num,ori_pos,cds_len):
    new_pos=ori_pos+random_num
    if new_pos>cds_len:
        new_pos=new_pos-cds_len
    return new_pos

def check_junction (old_ls,shuffle_ls,random_int,file):
    i=0
    v=0
    while i<len(old_ls):
        if old_ls[i][1]-old_ls[i][0]!=shuffle_ls[i][1]-shuffle_ls[i][0]:# they span the circle junction
            v=1
#            print ("regenerating",random_int, file)
            break
        i=i+1
    return v


########################################prepare of shuffle    
def shuffle_pointer(file_pos_on_cds):

    pos_on_cds=open(file_pos_on_cds,"r")
    cds_len=int(pos_on_cds.readline().split("\t")[3])
    pos_on_cds.close()
    right_shift=random.randint(1,cds_len)

    pos_on_cds=open(file_pos_on_cds,"r")
    ori_ls=[]
    new_ls=[]

    #pointer_used_in_shuffle=open(file_pos_on_cds+"_used","w")
    for line in pos_on_cds:
        pointer1=int(line.split("\t")[4])
        pointer2=int(line.split("\t")[5])
        feature=line[:-1].split("\t")[-1]
        info1=line.split("\t")[1]
        info2=line.split("\t")[2]
        info3=line.split("\t")[6]
        info4=line.split("\t")[7]
        info=("\t").join([info1,info2,info3,info4])
        if pointer1!=-1 and pointer2!=-1:# 
            pointer_s=min(pointer1,pointer2) #1based
            pointer_e=max(pointer1,pointer2)#1based
            ori_ls.append([pointer_s,pointer_e,feature,info])

            new_pt_s=circle_shuffle(right_shift,pointer_s,cds_len)
            new_pt_e=circle_shuffle(right_shift,pointer_e,cds_len)
            new_ls.append([new_pt_s,new_pt_e,feature,info])

    pos_on_cds.close()


    v=check_junction(ori_ls,new_ls,right_shift,file_pos_on_cds)
    
    k=0          
    while v==1 and k<=1000:
        right_shift=random.randint(1,cds_len)
        new_ls=[]
        for pointer in ori_ls:
            new_pt_s=circle_shuffle(right_shift,pointer[0],cds_len)
            new_pt_e=circle_shuffle(right_shift,pointer[1],cds_len)
            new_ls.append([new_pt_s,new_pt_e,pointer[2],pointer[3]])
        v=check_junction(ori_ls,new_ls,right_shift,file_pos_on_cds)
        k=k+1

    if v==0:
        return new_ls,cds_len
    else:
        print (file_pos_on_cds)

        # output=open(file_pos_on_cds+"test_simulation","w")
        # i=0
        # while i <len(ori_ls):
        #     output.write(file_pos_on_cds+"\t%d\t%d\t%d\t%d\n" % (ori_ls[i][0],ori_ls[i][1],new_ls[i][0],new_ls[i][1]))
        #     i=i+1
        # output.close()

    #return None

    
    
for line in OG_file:
    OG=line.split("\t")[0]
    #print (OG)
    oxy_shuffle_pointer_ls=[]
    ewoo_shuffle_pointer_ls=[]
    tet_shuffle_pointer_ls=[]
    if path.exists("./Oxy_pos_on_cds/"+OG+".pointer_pos_on_cds.txt"):
        oxy_shuffle_pointer_ls,oxy_cds_len=shuffle_pointer("./Oxy_pos_on_cds/"+OG+".pointer_pos_on_cds.txt")#pointer_pos_ls=[] if pointer1 or pointer2=-1
    if path.exists("./Ewoo_pos_on_cds/"+OG+".pointer_pos_on_cds.txt"):
        ewoo_shuffle_pointer_ls,ewoo_cds_len=shuffle_pointer("./Ewoo_pos_on_cds/"+OG+".pointer_pos_on_cds.txt")
    if path.exists("./Tet_pos_on_cds/"+OG+".pointer_pos_on_cds.txt"):
        tet_shuffle_pointer_ls,tet_cds_len=shuffle_pointer("./Tet_pos_on_cds/"+OG+".pointer_pos_on_cds.txt")
        

    #######################convert new position to gapped position
    if path.exists("./Single_Copy_Orthologue_gap_gff/"+OG+".gap.gff"):
        gap_gff=open("./Single_Copy_Orthologue_gap_gff/"+OG+".gap.gff","r")
        oxy_rl=[]
        ewoo_rl=[]
        tet_rl=[]
        for line in gap_gff:
            if line[0]!="#":
                gap_1=int(line.split("\t")[4])
                gap_2=int(line.split("\t")[5])
                if line.split("\t")[1]=="oxy":
                    oxy_pro=line.split("\t")[2]
                    oxy_len=int(line.split("\t")[3])
                    oxy_rl.append([gap_1,gap_2])
                elif line.split("\t")[1]=="tet":
                    tet_pro=line.split("\t")[2]
                    tet_rl.append([gap_1,gap_2])
                    tet_len=int(line.split("\t")[3])
                elif line.split("\t")[1]=="ewoo":
                    ewoo_pro=line.split("\t")[2]
                    ewoo_rl.append([gap_1,gap_2])
                    ewoo_len=int(line.split("\t")[3])
                    
                    
   #############output#################                 

        f=open("./shuffle_simulation/tmp/"+OG+"_oxy.txt","w")
        for pointer in oxy_shuffle_pointer_ls:
            pt_s=pointer[0]-1
            pt_e=pointer[1]
            feature=pointer[2]
            info=pointer[3]
            pointer_aln_1=convert2aln(pt_s,oxy_rl)
            pointer_aln_2=convert2aln(pt_e,oxy_rl)
            #mid=(pointer_aln_1+pointer_aln_2)/2
            f.write("%s\t%d\t%d\t%s\n" % (info,pointer_aln_1,pointer_aln_2,feature))
        f.close()
        
        f=open("./shuffle_simulation/tmp/"+OG+"_ewoo.txt","w")
        for pointer in ewoo_shuffle_pointer_ls:
            pt_s=pointer[0]-1
            pt_e=pointer[1]
            feature=pointer[2]
            info=pointer[3]
            pointer_aln_1=convert2aln(pt_s,ewoo_rl)
            pointer_aln_2=convert2aln(pt_e,ewoo_rl)
            #mid=(pointer_aln_1+pointer_aln_2)/2
            f.write("%s\t%d\t%d\t%s\n" % (info,pointer_aln_1,pointer_aln_2,feature))
        f.close()
        
        f=open("./shuffle_simulation/tmp/"+OG+"_tet.txt","w")
        for pointer in tet_shuffle_pointer_ls:
            pt_s=pointer[0]-1
            pt_e=pointer[1]
            feature=pointer[2]
            info=pointer[3]
            pointer_aln_1=convert2aln(pt_s,tet_rl)
            pointer_aln_2=convert2aln(pt_e,tet_rl)
            #mid=(pointer_aln_1+pointer_aln_2)/2
            f.write("%s\t%d\t%d\t%s\n" % (info,pointer_aln_1,pointer_aln_2,feature))
        f.close()

        gap_gff.close()

OG_file.close()
