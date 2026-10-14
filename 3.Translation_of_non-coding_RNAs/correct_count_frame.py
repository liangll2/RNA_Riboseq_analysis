#!/usr/bin/python
#-*-coding:utf-8-*-
#Author: zhang wei
#date : 2019-03-12

def choice_len(cds_bed,bam):
	bed = {}
	with open(cds_bed) as f:
		for line in f:
			line = line.strip().split('\t')
			bed[line[0]]=[int(line[1]),int(line[2])]

	dic = {}
	for i in range(25,36):
		dic[i] = [0,0,0]
	if not os.path.exists(bam):
		print "error : %s is not found" % bam
		sys.exit()
	else :
		pass
	
	f1 = os.popen("/Bio/bin/samtools view -F 4 %s" % bam)
	for line in f1:
		line = line.strip()
		each = line.split('\t')
		if each[2] not in bed:
			continue
		start = bed[each[2]][0]+1
		end = bed[each[2]][1]
		length = len(each[9])
		if int(each[3]) >=start-16 and int(each[3]) <=end+3 and length in dic:
			n = (int(each[3])-start)%3
			if n==0:
				dic[length][0] +=1
			if n==1:
				dic[length][1] +=1
			if n==2:
				dic[length][2] +=1
	f1.close()
	
	dic_frame = {}
	for reads_len in dic:
		max_number = max(dic[reads_len])
		sum_number = sum(dic[reads_len])
		if max_number >= 0.5*sum_number:
			max_index = dic[reads_len].index(max_number)
			if max_index == 0:
				dic_frame[reads_len] = 1
			elif max_index == 1:
				dic_frame[reads_len] = 2
			else :
				dic_frame[reads_len] = 3

	return dic_frame

def codon_count(cds_bed,orf_bed,bam,outfile):
	dic_bed = defaultdict(dict)
	dic_codon = {}
	bed = open(orf_bed)
	for line in bed:
		line = line.strip().split("\t")
		dic_bed[line[0]][line[3]] = [int(line[1]),int(line[2])]
		dic_codon[line[3]] = [0,0,0]
	bed.close()

	dic_frame = choice_len(cds_bed,bam)
	print dic_frame

	tmp = os.popen("/Bio/bin/samtools view -F 4 %s" % bam)
	for line in tmp:
		line = line.strip().split("\t")
		if len(line[9]) not in dic_frame:
			continue
		if line[2] not in dic_bed:
			continue
		pos = int(line[3])
		for orf in dic_bed[line[2]]:
			start = dic_bed[line[2]][orf][0]+1
			end = dic_bed[line[2]][orf][1]
			if pos < start-16 or pos > end:
				continue
			codon = (pos-start)%3
			if codon == 0 and dic_frame[len(line[9])] == 1:
				dic_codon[orf][0] += 1
			if codon == 0 and dic_frame[len(line[9])] == 2:
				dic_codon[orf][2] += 1
			if codon == 0 and dic_frame[len(line[9])] == 3:
				dic_codon[orf][1] += 1
			if codon == 1 and dic_frame[len(line[9])] == 1:
				dic_codon[orf][1] += 1
			if codon == 1 and dic_frame[len(line[9])] == 2:
				dic_codon[orf][0] += 1
			if codon == 1 and dic_frame[len(line[9])] == 3:
				dic_codon[orf][2] += 1
			if codon == 2 and dic_frame[len(line[9])] == 1:
				dic_codon[orf][2] += 1
			if codon == 2 and dic_frame[len(line[9])] == 2:
				dic_codon[orf][1] += 1
			if codon == 2 and dic_frame[len(line[9])] == 3:
				dic_codon[orf][0] += 1
	tmp.close()
	
	orf_info = open(orf_bed)
	out = open(outfile,'w')
	out.write("id\torf_type\tsource_trans\tstart\tend\tlength\tframe1_count\tframe2_count\tframe3_count\n")
	for each in orf_info:
		each = each.strip().split("\t")
		start = str(int(each[1])+1)
		orf_len = str(int(each[2]) - int(each[1]))
		if each[3] not in dic_codon:
			frame = "0\t0\t0"
		else :
			frame = "\t".join([str(x) for x in dic_codon[each[3]]])
		out.write(each[3]+"\t"+each[4]+"\t"+each[0]+"\t"+start+"\t"+each[2]+"\t"+orf_len+"\t"+frame+"\n")
	orf_info.close()
	out.close()


if __name__ == "__main__":
	import sys,os,re
	from collections import Counter
	from collections import defaultdict
	if len(sys.argv) != 5:
		print "Usage: python correct_count_frame.py <cds_bed> <orf_bed> <bam> <outfile>\n"
		sys.exit(0)

	codon_count(*sys.argv[1:])

