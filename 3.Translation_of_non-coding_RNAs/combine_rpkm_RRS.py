#!/usr/bin/python
#-*-coding:utf-8-*-
#Author : zhang wei
#date : 2019-03-26

def get_map_reads(sample_lis,indir):
	map_reads = 0
	for name in sample_lis:
		stat = indir+"/bams/%s.stat.xls" % (name)
		tmp = open(stat)
		i = 0
		for each in tmp:
			each = each.strip("\n")
			i += 1
			if i == 2:
				map_reads += int(each.split("\t")[1])
		tmp.close()
	
	return map_reads

def get_rpkm(count,total,length):
	rpkm = round(float(10**6*count)/(total*(float(length)/1000)),2)
	return rpkm

def combine_rss(sample_lis,indir,outfile):
	total = get_map_reads(sample_lis,indir)
	
	dic_count,dic_len = {},{}
	dic_info = {}
	i = 0
	for name in sample_lis:
		i += 1
		sample_exp = indir+"/exp/%s.ORF.exp.xls" % (name)
		tmp = open(sample_exp)
		for line in tmp:
			line = line.strip()
			each = line.split("\t")
			if line.startswith("id"):
				dic_count["id"] = line
				orf = each.index("orf_count")
				utr = each.index("3utr_count")
				orflen_index = each.index("orf_length")
				utrlen_index = each.index("3utr_length")
				continue
			key = each[0]
			if i == 1:
				dic_len[key] = [int(each[orflen_index]),int(each[utrlen_index])]
				dic_count[key] = [int(each[orf]),int(each[utr])]
				dic_info[key] = "\t".join(each[:utrlen_index+1])
			else :
				dic_count[key][0] += int(each[orf])
				dic_count[key][1] += int(each[utr])
		tmp.close()

	expfile = indir+"/RRS_ORFscore/allsample.ORFscore.xls"
	if not os.path.exists(expfile):
		print "error : there is no such file: %s" % expfile
		sys.exit()
	else :
		pass
	f = open(expfile)
	out = open(outfile,'w')
	for line in f:
		line = line.strip()
		each = line.split("\t")
		if line.startswith("id"):
			out.write(dic_count["id"]+"\n")
			continue
		orf_count = dic_count[each[0]][0]
		orf_length = dic_len[each[0]][0]
		utr_count = dic_count[each[0]][1]
		utr_length = dic_len[each[0]][1]
		orf_rpkm = get_rpkm(orf_count,total,orf_length)
		if utr_count != 0 and utr_length != 0:
			utr_rpkm = get_rpkm(utr_count,total,utr_length)
			if utr_rpkm == 0:
				rrs = orf_rpkm
			else :
				rrs = orf_rpkm/utr_rpkm
		elif utr_count == 0 and utr_length != 0 :
			tmp_utr_count = 1
			utr_rpkm = get_rpkm(tmp_utr_count,total,utr_length)
			if utr_rpkm == 0 :
				rrs = orf_rpkm
			else :
				rrs = orf_rpkm/utr_rpkm
			utr_rpkm = 0.00
		else :
			utr_rpkm = 0.00
			rrs = orf_rpkm
		out.write(dic_info[each[0]]+"\t"+str(orf_count)+"\t"+str(utr_count)+"\t"+str(orf_rpkm)+"\t"+str(utr_rpkm)+"\t"+str(rrs)+"\n")
	out.close()

def plot(prefix):
	script = '''#!/Bio/bin/Rscript
data<-read.table("%s.xls",header = T)
data<-data[which(data$RRS>0),]
data_uORF<-data[which(data$orf_type=="uORF"),]
data_dORF<-data[which(data$orf_type=="dORF"),]
data_mORF<-data[which(data$orf_type=="mORF"),]
ymax<-max(log10(data[which(data$RRS!="NA"),]$RRS))
ymin<-min(log10(data[which(data$RRS!="NA"),]$RRS))
uORF_rrs<-density(log10(data_uORF[which(data_uORF$RRS!="NA"),]$RRS))
dORF_rrs<-density(log10(data_dORF[which(data_dORF$RRS!="NA"),]$RRS))
mORF_rrs<-density(log10(data_mORF[which(data_mORF$RRS!="NA"),]$RRS))
pdf(file='%s.pdf',width = 10,height=8)
par(mar=c(4,5,1,1))
ya<-max(uORF_rrs$y,dORF_rrs$y,mORF_rrs$y)
yn<-min(uORF_rrs$y,dORF_rrs$y,mORF_rrs$y)
plot(uORF_rrs,ylim=c(yn,ya),xlim=c(ymin,ymax),type="n", xlab="RRS", ylab="density",main="",cex.lab=2,cex.axis=1.8)
polygon(uORF_rrs,col=rgb(241,64,59,100,maxColorValue = 255))
polygon(dORF_rrs,col=rgb(116,198,239,150,maxColorValue = 255))
polygon(mORF_rrs,col=rgb(141,88,165,150,maxColorValue = 255))
legend("topleft",c("uORF","dORF","mORF"),pch=20,col=c(rgb(241,64,59,150,maxColorValue = 255),rgb(116,198,239,150,maxColorValue = 255),rgb(141,88,165,150,maxColorValue = 255)),cex=2)
dev.off()
png(file='%s.png',width = 800,height=640)
par(mar=c(4,5,1,1))
ya<-max(uORF_rrs$y,dORF_rrs$y,mORF_rrs$y)
yn<-min(uORF_rrs$y,dORF_rrs$y,mORF_rrs$y)
plot(uORF_rrs,ylim=c(yn,ya),xlim=c(ymin,ymax),type="n", xlab="RRS", ylab="density",main="",cex.lab=2,cex.axis=1.8)
polygon(uORF_rrs,col=rgb(241,64,59,100,maxColorValue = 255))
polygon(dORF_rrs,col=rgb(116,198,239,150,maxColorValue = 255))
polygon(mORF_rrs,col=rgb(141,88,165,150,maxColorValue = 255))
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

xmax <- max(log10(data[which(data$RRS!="NA"),]$RRS))
xmin <- min(log10(data[which(data$RRS!="NA" & data$RRS!=0),]$RRS))
uORF_rrs <- density(log10(data_uORF[which(data_uORF$RRS!="NA"),]$RRS))
dORF_rrs <- density(log10(data_dORF[which(data_dORF$RRS!="NA"),]$RRS))
mORF_rrs <- density(log10(data_mORF[which(data_mORF$RRS!="NA"),]$RRS))
if(nrow(data_lncORF) >= 2){
    lncORF_rrs <- density(log10(data_lncORF[which(data_lncORF$RRS!="NA"),]$RRS))
    ymax <- max(uORF_rrs$y,dORF_rrs$y,mORF_rrs$y,lncORF_rrs$y)
    ymin <- min(uORF_rrs$y,dORF_rrs$y,mORF_rrs$y,lncORF_rrs$y)
} else {
    lncORF_rrs <- c()
    ymax <- max(uORF_rrs$y,dORF_rrs$y,mORF_rrs$y)
    ymin <- min(uORF_rrs$y,dORF_rrs$y,mORF_rrs$y)
}

#colors <- c(unlist(strsplit('#3A973699,#D3242199,#2771A799,#C6AFD199', ',')))
colors <- c(unlist(strsplit('#00468B99,#935E9F99,#42C1BB99,#EDE54699', ',')))

pdf(file = paste0(outp, ".pdf"), width = 10, height = 8)
par(mar = c(4,5,1,1))
if(nrow(data_lncORF) >= 2) {
    plot(mORF_rrs, xlim=c(xmin,xmax), ylim=c(ymin,ymax), type="n", xlab="log10(RRS)", 
    ylab="density", main="", cex.lab=2, cex.axis=1.8)
    polygon(mORF_rrs,col=colors[4],border=NA)
    polygon(dORF_rrs,col=colors[3],border=NA)
    polygon(uORF_rrs,col=colors[2],border=NA)
    polygon(lncORF_rrs,col=colors[1],border=NA)
    legend("topleft", c("lncORF","uORF","dORF","mORF"), pch=20, col=colors, cex=2)
} else {
    plot(mORF_rrs, xlim=c(xmin,xmax), ylim=c(ymin,ymax), type="n", xlab="log10(RRS)", 
    ylab="density", main="", cex.lab=2, cex.axis=1.8)
    polygon(mORF_rrs,col=colors[4],border=NA)
    polygon(dORF_rrs,col=colors[3],border=NA)
    polygon(uORF_rrs,col=colors[2],border=NA)
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

	all_file = outdir+"/allsample.RRS.xls"
	combine_rss(sample,indir,all_file)
	plot2(outdir+"/allsample.RRS")
	for item in group:
		sample_lis = group_dic[item]
		outfile = outdir+"/%s.RRS.xls" % item
		combine_rss(sample_lis,indir,outfile)
if __name__ == "__main__":
	import sys,os
	if len(sys.argv) != 4:
		print "Usage : pythpn combine_rpkm_RRS.py <Groups.txt> <indir> <outdir>\n"
		sys.exit()

	main(*sys.argv[1:])

