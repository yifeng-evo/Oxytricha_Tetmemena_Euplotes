import sys
f=open(sys.argv[1],"r")
output=open("concise_below30"+sys.argv[1],"w")

info_dic={}
for line in f:
	MAC=line.split("\t")[0]
	start=int(line.split("\t")[1])
	end=int(line.split("\t")[2])
	length=end-start+1
	feature=line.split("\t")[3][:-1] #remove\n
	if length<=30:
		if MAC not in info_dic:
			info_dic[MAC]=[[start,end,feature]]
		else:
			v=0
			for pointer in info_dic[MAC]:
				s=pointer[0]
				e=pointer[1]
				l=e-s+1
				if s<=start<=e or s<=end<=e: #overlap
					if length>l: #keep the longest pointer, throw away the shorter one
						info_dic[MAC].remove(pointer)
						info_dic[MAC].append([start,end,feature])
						print (MAC,pointer)
						v=1
						break
					else:
						v=1
						break # keep the original 
				elif start<=s and e<=end:
					info_dic[MAC].remove(pointer)
					info_dic[MAC].append([start,end,feature])
					print (MAC,pointer)
					v=1
					break
			if v==0:
				info_dic[MAC].append([start,end,feature])

for MAC in info_dic:
	for pointer in info_dic[MAC]:
		output.write("%s\t%d\t%d\t%s\n" % (MAC,pointer[0],pointer[1],pointer[2]))

output.close()
f.close()
