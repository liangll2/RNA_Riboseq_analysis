#!/usr/bin/python
#-*-coding:utf-8-*-
#Author : zhang wei
#date : 2019-03-26

def ORFscore(frame_count):
	if frame_count == [0,0,0]:
		orf_score = "NA"
	else :
		mean = np.mean(frame_count)
		value = sum([float((x-mean)**2)/mean for x in frame_count]) +1
		if frame_count[0] < frame_count[1] or frame_count[0] < frame_count[2]:
			orf_score = math.log(value,2)*-1
		else :
			orf_score = math.log(value,2)*1

	return orf_score

def combine_score(sample_lis,indir,outfile):
	dic_count = {}
	keys = []
	out = open(outfile,'w')
	out2 = open(outfile+".withframe", 'w')
	i = 0
	for name in sample_lis:
		line = indir+"/%s.orf_frame.xls" % name
		if not os.path.exists(line):
			print "error there is no such file : %s" % line
			sys.exit()
		else :
			pass
		i += 1
		with open(line) as tmp:
			for each in tmp:
				each = each.strip().split("\t")
				key = "\t".join(each[:6])
				if i == 1:
					keys.append(key)
					if key.startswith("id"):
						dic_count[key] = "\t".join(each[6:])
					else :
						dic_count[key] = [int(each[6]),int(each[7]),int(each[8])]
				else :
					if key.startswith("id"):
						continue
					dic_count[key][0] += int(each[6])
					dic_count[key][1] += int(each[7])
					dic_count[key][2] += int(each[8])

	for x in keys:
		if x.startswith("id"):
			out.write(x+"\tORFscore\n")
			out2.write(x+"\t"+ dic_count[x]+"\tORFscore\n")
			continue
		orf_score = ORFscore(dic_count[x])
		out.write(x+"\t"+str(orf_score)+"\n")
		out2.write(x+"\t"+"\t".join([str(y) for y in dic_count[x]])+"\t"+str(orf_score)+"\n")
	out.close()
	out2.close()

def plot(prefix):
	script = '''#!/Bio/bin/Rscript
data<-read.table("%s.xls",header = T)
data_uORF<-data[which(data$orf_type=="uORF"),]
data_dORF<-data[which(data$orf_type=="dORF"),]
data_mORF<-data[which(data$orf_type=="mORF"),]
xmax<-max(data[which(data$ORFscore!="NA"),]$ORFscore)
xmin<-min(data[which(data$ORFscore!="NA"),]$ORFscore)
uORF_score<-density(data_uORF[which(data_uORF$ORFscore!="NA"),]$ORFscore)
dORF_score<-density(data_dORF[which(data_dORF$ORFscore!="NA"),]$ORFscore)
mORF_score<-density(data_mORF[which(data_mORF$ORFscore!="NA"),]$ORFscore)
ymax<-max(uORF_score$y,dORF_score$y,mORF_score$y)
ymin<-min(uORF_score$y,dORF_score$y,mORF_score$y)
pdf(file='%s.pdf',width = 10,height=8)
par(mar=c(4,5,1,1))
plot(uORF_score,xlim=c(xmin,xmax),ylim=c(ymin,ymax),type="n", xlab="ORFscore", ylab="density",main="",cex.lab=2,cex.axis=1.8)
polygon(uORF_score,col=rgb(241,64,59,100,maxColorValue = 255))
polygon(dORF_score,col=rgb(116,198,239,150,maxColorValue = 255))
polygon(mORF_score,col=rgb(141,88,165,150,maxColorValue = 255))
legend("topleft",c("uORF","dORF","mORF"),pch=20,col=c(rgb(241,64,59,150,maxColorValue = 255),rgb(116,198,239,150,maxColorValue = 255),rgb(141,88,165,150,maxColorValue = 255)),cex=2)
dev.off()
png(file='%s.png',width = 800,height=640)
par(mar=c(4,5,1,1))
plot(uORF_score,xlim=c(xmin,xmax),ylim=c(ymin,ymax),type="n", xlab="ORFscore", ylab="density",main="",cex.lab=2,cex.axis=1.8)
polygon(uORF_score,col=rgb(241,64,59,100,maxColorValue = 255))
polygon(dORF_score,col=rgb(116,198,239,150,maxColorValue = 255))
polygon(mORF_score,col=rgb(141,88,165,150,maxColorValue = 255))
legend("topleft",c("uORF","dORF","mORF"),pch=20,col=c(rgb(241,64,59,150,maxColorValue = 255),rgb(116,198,239,150,maxColorValue = 255),rgb(141,88,165,150,maxColorValue = 255)),cex=2)
dev.off()
''' % (prefix,prefix,prefix)
	with open(prefix+".draw.r",'w') as f:
		f.write(script)
	assert not os.system("/Bio/bin/Rscript %s.draw.r" % prefix)

def plot2(prefix):
	script = '''#!/Bio/bin/Rscript

file <- "%s.xls"
outp <- "%s"

data <- read.table(file, header = T)
data_lncORF <- data[which(data$orf_type=="lncORF"),]
data_uORF   <- data[which(data$orf_type=="uORF"),]
data_dORF   <- data[which(data$orf_type=="dORF"),]
data_mORF   <- data[which(data$orf_type=="mORF"),]

xmax <- max(data[which(data$ORFscore!="NA"),]$ORFscore)
xmin <- min(data[which(data$ORFscore!="NA"),]$ORFscore)
uORF_score <- density(data_uORF[which(data_uORF$ORFscore!="NA"),]$ORFscore)
dORF_score <- density(data_dORF[which(data_dORF$ORFscore!="NA"),]$ORFscore)
mORF_score <- density(data_mORF[which(data_mORF$ORFscore!="NA"),]$ORFscore)
if(nrow(data_lncORF) >= 2){
    lncORF_score <- density(data_lncORF[which(data_lncORF$ORFscore!="NA"),]$ORFscore)
    ymax <- max(uORF_score$y,dORF_score$y,mORF_score$y,lncORF_score$y)
    ymin <- min(uORF_score$y,dORF_score$y,mORF_score$y,lncORF_score$y)
} else {
    lncORF_score <- c()
    ymax <- max(uORF_score$y,dORF_score$y,mORF_score$y)
    ymin <- min(uORF_score$y,dORF_score$y,mORF_score$y)
}

#colors <- c(unlist(strsplit('#3A973699,#D3242199,#2771A799,#C6AFD199', ',')))
colors <- c(unlist(strsplit('#00468B99,#935E9F99,#42C1BB99,#EDE54699', ',')))

pdf(file = paste0(outp, ".pdf"), width = 10, height = 8)
par(mar = c(4,5,1,1))
if(nrow(data_lncORF) >= 2) {
    plot(mORF_score, xlim=c(xmin,xmax), ylim=c(ymin,ymax), type="n", xlab="ORFscore", 
    ylab="density", main="", cex.lab=2, cex.axis=1.8)
    polygon(mORF_score,col=colors[4],border=NA)
    polygon(dORF_score,col=colors[3],border=NA)
    polygon(uORF_score,col=colors[2],border=NA)
    polygon(lncORF_score,col=colors[1],border=NA)
    legend("topleft", c("lncORF","uORF","dORF","mORF"), pch=20, col=colors, cex=2)
} else {
    plot(mORF_score, xlim=c(xmin,xmax), ylim=c(ymin,ymax), type="n", xlab="ORFscore", 
    ylab="density", main="", cex.lab=2, cex.axis=1.8)
    polygon(mORF_score,col=colors[4],border=NA)
    polygon(dORF_score,col=colors[3],border=NA)
    polygon(uORF_score,col=colors[2],border=NA)
    legend("topleft", c("uORF","dORF","mORF"), pch=20, col=colors[2:4], cex=2)
}
dev.off()
system(paste0("/Bio/usr/bin/convert -density 300 ", outp, ".pdf", " ", outp, ".png"))

''' % (prefix,prefix)
	with open(prefix+".draw.r",'w') as f:
		f.write(script)
	assert not os.system("/Bio/bin/Rscript %s.draw.r" % prefix)

def main(groupfile,indir,outdir):
	group_dic = {}
	group,sample = [],[]
	with open(groupfile) as tmp:
		for line in tmp:
			line = line.strip().split("\t")
			group.append(line[0])
			sample.append(line[1])
			if line[0] not in group_dic:
				group_dic[line[0]] = [line[1]]
			else :
				group_dic[line[0]].append(line[1])
	
	all_file = outdir+"/allsample.ORFscore.xls"
	combine_score(sample,indir,all_file)
	plot2(outdir+"/allsample.ORFscore")
	for item in group:
		sample_lis = group_dic[item]
		outfile = outdir+"/%s.ORFscore.xls" % item
		combine_score(sample_lis,indir,outfile)
if __name__ == "__main__":
	import sys,os,re
	import math
	import numpy as np
	if len(sys.argv) != 4:
		print "Usage : python combine_frame_ORFscore.py <Groups.txt> <indir> <outdir>\n"
		sys.exit()

	main(*sys.argv[1:])

