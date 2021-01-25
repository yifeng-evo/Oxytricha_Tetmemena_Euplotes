import matplotlib .pyplot as plt
import os.path
from os import path
from matplotlib.backends.backend_pdf import PdfPages
OG_file=open("./Orthogroups/Orthogroups_SingleCopyOrthologues_3species.txt","r")

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

def draw_pointer_species(folder,species,ls):
    if species =="oxy":
        v=1
    elif species=="tet":
        v=2
    else:
        v=3
    if path.exists(folder+OG+".pointer_pos_on_cds.txt"):
        pointer_on_cds=open(folder+OG+".pointer_pos_on_cds.txt","r")
        for line in pointer_on_cds:
            pointer_1=int(line.split("\t")[4])
            pointer_2=int(line.split("\t")[5])
            feature=line.split("\t")[8][:-1]
            if feature=="scrambled":
                fc="red"#orange
            elif feature=="nonscrambled":
                fc="blue"
            else:
                fc="purple" #purple for either

            if pointer_1!=-1 and pointer_2!=-1:
                pointer_s=min(pointer_1,pointer_2)-1#convert to 0 based
                pointer_e=max(pointer_1,pointer_2)
                pointer_aln_1=convert2aln(pointer_s,ls)
                pointer_aln_2=convert2aln(pointer_e,ls)
                f=open("./shuffle_simulation/real/"+OG+"_"+species+".txt","a")
                f.write(line.split("\t")[1]+"\t"+line.split("\t")[2]+"\t"+line.split("\t")[6]+"\t"+line.split("\t")[7]+"\t"+str(pointer_aln_1)+"\t"+str(pointer_aln_2)+"\t"+feature+"\n")#$1\t$2\toriginal pointer start\t original pointer end\t gff pointer 0
                f.close()
                x=[pointer_aln_1,pointer_aln_2]
                y=[v,v]
                plt.plot(x,y,alpha=0.5,color=fc,markersize=15,marker="|",linewidth=12)
            elif pointer_1==-1 and pointer_2==-1:
                pass
            else:
                pointer_s=max(pointer_1,pointer_2)-1
                pointer_e=max(pointer_1,pointer_2)
                pointer_aln_1=convert2aln(pointer_s,ls)
                pointer_aln_2=convert2aln(pointer_e,ls)
                x=[pointer_aln_1,pointer_aln_2]
                y=[v,v]
                plt.plot(x,y,alpha=0.5,color=fc,markersize=15,marker="d",linewidth=12)


def draw_intron_species(folder,species,ls):
    if species =="oxy":
        v=1
#        sc="gray"
    elif species=="tet":
        v=2
#        sc="gray"
    else:
        v=3
#        sc="gray"
    sc="green"
    if path.exists(folder+OG+".exon_pos_on_cds.txt"):
        exon_on_cds=open(folder+OG+".exon_pos_on_cds.txt","r")
        for line in exon_on_cds:
            exon_1=int(line.split("\t")[4])
            exon_2=int(line.split("\t")[5])  
            exon_s=min(exon_1,exon_2)-1
            exon_e=exon_s+1          
            exon_aln_1=convert2aln(exon_s,ls)  
            exon_aln_2=convert2aln(exon_e,ls)  
            x=[exon_aln_1,exon_aln_2]
            y=[v,v]           
            plt.plot(x,y,alpha=0.5,color=sc,markersize=15,marker="o",linewidth=12)



for line in OG_file:
    OG=line.split("\t")[0]
    ###################gap info 
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
        fig= plt.figure(figsize=(20,2))
        plt.axes(frameon=False)
        plt.yticks([])#hide y axis
        plt.xlim(-0.5,oxy_len+0.5)#set x axis range
        plt.ylim(0.5,3.5)
         #hiding axes
        plt.text(0.001,0.15,s="Oxy:"+oxy_pro,transform=plt.gcf().transFigure)
        plt.text(0.001,0.5,s="Tet:"+tet_pro,transform=plt.gcf().transFigure)
        plt.text(0.001,0.8,s="Ewoo:"+ewoo_pro,transform=plt.gcf().transFigure)
        plt.hlines(y=1,xmin=0,xmax=oxy_len,alpha=0.5,color="gray")#xmax =alignment length
        plt.hlines(y=2,xmin=0,xmax=tet_len,alpha=0.5,color="gray")
        plt.hlines(y=3,xmin=0,xmax=ewoo_len,alpha=0.5,color="gray")
        ####what if oxy_rl =[]?
        for gap in oxy_rl:
            plt.hlines(y=1,xmin=gap[0],xmax=gap[1],alpha=0.5,color="gray",linewidth=10)
        for gap in tet_rl:
            plt.hlines(y=2,xmin=gap[0],xmax=gap[1],alpha=0.5,color="gray",linewidth=10)
        for gap in ewoo_rl:
            plt.hlines(y=3,xmin=gap[0],xmax=gap[1],alpha=0.5,color="gray",linewidth=10)
        if oxy_rl==[]:
            plt.hlines(y=1,xmin=0,xmax=oxy_len,alpha=0.5,color="gray",linewidth=10)
        if tet_rl==[]:
            plt.hlines(y=2,xmin=0,xmax=tet_len,alpha=0.5,color="gray",linewidth=10)
        if ewoo_rl==[]:
            plt.hlines(y=3,xmin=0,xmax=ewoo_len,alpha=0.5,color="gray",linewidth=10)

            
    ####################transfer cds pos info to alignment pos ############################
        draw_pointer_species("./Ewoo_pos_on_cds/","ewoo",ewoo_rl)
        draw_intron_species("./Ewoo_pos_on_cds/","ewoo",ewoo_rl)
                
        draw_pointer_species("./Oxy_pos_on_cds/","oxy",oxy_rl)
        draw_intron_species("./Oxy_pos_on_cds/","oxy",oxy_rl)
        
        draw_pointer_species("./Tet_pos_on_cds/","tet",tet_rl)
        draw_intron_species("./Tet_pos_on_cds/","tet",tet_rl)     

            
        # with PdfPages("./Ortholog_comparison_visualization/all_OG.pdf") as pdf:
        #     plt.title(OG)
        #     pdf.sacefig()  
        plt.savefig("./Ortholog_comparison_visualization/"+OG+".pdf",format="pdf",bbox_inches="tight")  
        plt.close()
    #else:
     #   print ("no gap file for "+OG) # if 
        ####

OG_file.close()
