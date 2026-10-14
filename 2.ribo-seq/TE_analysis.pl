#! /usr/bin/perl

use utf8;
use strict;
use warnings;
#use Getopt::Long;
#use Bio::SeqIO;
#use Bio::Seq;
#use List::Util qw/sum min max/;
#use List::MoreUtils qw/uniq/;
#use File::Basename qw/basename dirname/;
#use File::Spec::Functions qw/rel2abs/;
#use FindBin qw/$Bin $Script/;
#use lib $Bin;
#Author: jhguo@genedenovo.com
die "perl $0 <RNA.exp.xls> <RIB.exp.xls> <compairs> <output_prefix>\n" unless(@ARGV eq 4);
open RNA,$ARGV[0] or die "failed to open RNA exps file";
my %cal;
my $head = <RNA>;
chomp($head);
my @heads = split("\t",$head);
my %rna_exp;
while(<RNA>){
	chomp;
	my @as = split("\t",$_);
	for (my $i = 1;$i<@as;$i++){
		next if ($heads[$i] =~ /_count/);
		$rna_exp{$as[0]}{$heads[$i]} = $as[$i];
	}
	$cal{$as[0]} ++;
}
close RNA;
my %rib_exp;
open RIB,$ARGV[1] or die "failed to open RIB exps file";
$head = <RIB>;
chomp($head);
@heads = split("\t",$head);
while(<RIB>){
	chomp;
	my @as = split("\t",$_);
	for (my $i = 1;$i<@as;$i++){
		next if ($heads[$i] =~ /_count/);
		$rib_exp{$as[0]}{$heads[$i]} = $as[$i];
	}
	$cal{$as[0]} ++;
}
close RIB;
open CP,$ARGV[2] or die $!;
my %cps;
my %ofhs;
while(<CP>){
	chomp;
	my @as = split("\t",$_);
	$cps{$as[0]} = $as[1];
	open my $ofh,">$ARGV[3].$as[1].fordraw.txt" or die $!;
	print $ofh "gene\tlog2RNA\tlog2TE\n";
	$ofhs{$as[1]}=$ofh;
}
close CP;
open OUT,">$ARGV[3].tmp" or die "failed to write output";
my $headflag = 0;
foreach my $key(sort{$a cmp $b} keys %cal){
	if($cal{$key} == 2){
		#print head
		if($headflag == 0){
			print OUT "genes";
			foreach my $tekey (sort{$a cmp $b} keys %cps){
				my ($ribo,$rna) = split("-vs-",$tekey);
				print OUT "\t$ribo\_RIB_fpkm\t$rna\_RNA_fpkm\t$cps{$tekey}"
			}
			print OUT "\n";
			$headflag = 1;
		}
		print OUT $key;
		foreach my $tekey(sort {$a cmp $b} keys %cps){
			my ($ribo,$rna) = split("-vs-",$tekey);
			my $te ;
#			$rib_exp{$key}{"$ribo\_fpkm"} = 0.01 unless ($rib_exp{$key}{"$ribo\_fpkm"} > 0);
#			$rna_exp{$key}{"$rna\_fpkm"} = 0.01 unless ($rna_exp{$key}{"$rna\_fpkm"} > 0);
			if($rib_exp{$key}{"$ribo\_fpkm"} > 0 and $rna_exp{$key}{"$rna\_fpkm"} >0){
				$te = $rib_exp{$key}{"$ribo\_fpkm"}/$rna_exp{$key}{"$rna\_fpkm"};
				print OUT "\t".$rib_exp{$key}{"$ribo\_fpkm"}."\t".$rna_exp{$key}{"$rna\_fpkm"}."\t".sprintf("%.6f",$te);
#				print OUT "\t".sprintf("%.6f",$te);
				my $fs = $ofhs{$cps{$tekey}};
				print $fs $key."\t".sprintf("%.6f",log($rna_exp{$key}{"$rna\_fpkm"})/log(2))."\t".sprintf("%.6f",log($te)/log(2))."\n";
			}else {
				$te = "NA";
#				print OUT "\t$te";
				print OUT "\t".$rib_exp{$key}{"$ribo\_fpkm"}."\t".$rna_exp{$key}{"$rna\_fpkm"}."\t$te";
			}
		}
		print OUT "\n";
	}
}
foreach my $cp(keys %ofhs){
	open TMP,">$ARGV[3].r" or die $!;
	print TMP "
		data<-read.table(\"$ARGV[3].$cp.fordraw.txt\",header = T)
		cor <- round(cor(data\$log2RNA, data\$log2TE, method = \"pearson\"), 4)
		pdf(\"$ARGV[3].$cp.pdf\")
		par(mar=c(5,5,4,2))
		plot(data\$log2RNA,data\$log2TE,pch=16,cex=0.2,col=\"red\",xlab=\"log2(TE)\",ylab=\"log2 mRNA abundance (fpkm)\",cex.axis=2,cex.lab=2)
		legend(\"topleft\", legend = paste(\"Pearson correlation: \",cor,sep=\"\"), col =rgb(255,0,0,max=255,alpha=255), pch=20,cex=1.2)
		png(\"$ARGV[3].$cp.png\")
		par(mar=c(5,5,4,2))
		plot(data\$log2RNA,data\$log2TE,pch=16,cex=0.2,col=\"red\",xlab=\"log2 (TE)\",ylab=\"log2 mRNA abundance (fpkm)\",cex.axis=2,cex.lab=2)
		legend(\"topleft\", legend = paste(\"Pearson correlation: \",cor,sep=\"\"), col =rgb(255,0,0,max=255,alpha=255), pch=20,cex=1.2)
	";
	close TMP;
	`Rscript $ARGV[3].r`;
}
