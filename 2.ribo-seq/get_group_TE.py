#!/usr/bin/env python

def main(group_list,exp,outpre):
	groups,group_dic = [],{}
	f1 = open(group_list)
	for line in f1:
		line = line.strip().split("\t")
		if line[1] not in groups:
			groups.append(line[1])
			group_dic[line[1]] = [line[0]]
		else :
			group_dic[line[1]].append(line[0])
	f1.close()

	tmp = open(exp)
	group_exp = outpre+".xls"
	tmp_exp = outpre+".tmp"
	out = open(group_exp,'w')
	out2 = open(tmp_exp,'w')
	i = 0
	group_index = {}
	for line in tmp:
		line = line.strip().split("\t")
		i += 1
		if i == 1:
			out.write(line[0]+"\t"+"\t".join(groups)+"\n")
			out2.write(line[0]+"\t"+"\t".join(groups)+"\n")
			for group in groups:
				samples = group_dic[group]
				indexs = [line.index(sample) for sample in samples]
				group_index[group] = indexs
			continue
		group_values,tmp_values = [],[]
		for group in groups:
			inds = group_index[group]
			values = []
			for m in inds:
				if line[m] != "NA":
					values.append(float(line[m]))
			if values:
				mean_value = "%.6f" % (sum(values)/len(values))
				tmp_values.append(mean_value)
			else :
				mean_value = "NA"
				tmp_values.append("0")
			group_values.append(mean_value)
		out.write(line[0]+"\t"+"\t".join(group_values)+"\n")
		out2.write(line[0]+"\t"+"\t".join(tmp_values)+"\n")
	tmp.close();out.close();out.close()

if __name__=="__main__":
	import sys,re,os

	if len(sys.argv) != 4:
		print "USage : python get_group_TE.py <group_list> <TE_sample> <TE_group>\n"
		sys.exit()
	
	main(*sys.argv[1:])
