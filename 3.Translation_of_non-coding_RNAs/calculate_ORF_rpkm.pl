#! /usr/bin/perl

use utf8;
use strict;
use warnings;
use lib "/home/guojianhui/lib";
use Fzuguo;
use File::Basename;
my $pwd = dirname(__FILE__);
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
die "perl $0 <all_ORF.bed> <bam> <stat.xls> <outdir/prefix>\n" unless(@ARGV eq 4);
open BED,$ARGV[0] or die $!;
my %bed_uorf;
my %bed_uorf_utr;
my %bed_dorf;
my %bed_dorf_utr;
my %bed_morf;
my %bed_morf_utr;
my %exist_uorf;
my %exist_uorf_utr;
my %exist_dorf;
my %exist_dorf_utr;
my %exist_morf;
my %exist_morf_utr;
## bed format:
## bed{transid}[0] -> [start,end] ,......
while(<BED>){
	chomp;
	my @as = split("\t",$_);
	if($as[3] eq "uORF"){
		unless($exist_uorf{$as[0]."\t".$as[1]."\t".$as[2]}){
			push(@{$bed_uorf{$as[0]}},[$as[1],$as[2]]);
			$exist_uorf{$as[0]."\t".$as[1]."\t".$as[2]} = 1;
		}
		unless($exist_uorf_utr{$as[0]."\t".$as[4]."\t".$as[5]}){
			push(@{$bed_uorf_utr{$as[0]}},[$as[4],$as[5]]);
			$exist_uorf_utr{$as[0]."\t".$as[4]."\t".$as[5]} = 1;
		}
	}elsif($as[3] eq "dORF"){
		unless($exist_dorf{$as[0]."\t".$as[1]."\t".$as[2]}){
			push(@{$bed_dorf{$as[0]}},[$as[1],$as[2]]);
			$exist_dorf{$as[0]."\t".$as[1]."\t".$as[2]} = 1;
		}
		unless($exist_dorf_utr{$as[0]."\t".$as[4]."\t".$as[5]}){
			push(@{$bed_dorf_utr{$as[0]}},[$as[4],$as[5]]);
			$exist_dorf_utr{$as[0]."\t".$as[4]."\t".$as[5]} = 1;
		}
	}elsif($as[3] eq "mORF"){
		unless($exist_morf{$as[0]."\t".$as[1]."\t".$as[2]}){
			push(@{$bed_morf{$as[0]}},[$as[1],$as[2]]);
			$exist_morf{$as[0]."\t".$as[1]."\t".$as[2]} = 1;
		}
		unless($exist_morf_utr{$as[0]."\t".$as[4]."\t".$as[5]}){
			push(@{$bed_morf_utr{$as[0]}},[$as[4],$as[5]]);
			$exist_morf_utr{$as[0]."\t".$as[4]."\t".$as[5]} = 1;
		}
	}
}
close BED;


open FILE,$ARGV[2] or die "failed to open file1";
<FILE>;
my $map_gene_reads;
while(<FILE>){
	chomp;
	my @as = split("\t",$_);
	$map_gene_reads = $as[1];
}
close FILE;
if($ARGV[1] =~ /bam$/){
	open BAM,"samtools view $ARGV[1]|" or die $!;
}else{
	open BAM,$ARGV[1] or die $!;
}
my %ct_uorf;
my %ct_uorf_utr;
my %ct_dorf;
my %ct_dorf_utr;
my %ct_morf;
my %ct_morf_utr;
#不考虑平均分配
my $flag = "";
while(<BAM>){
	chomp;
	my @as = split("\t",$_);
	if(defined $bed_uorf{$as[2]}){
		my $start = $as[3];
		my $maplength = 0;
		while($as[5] =~ /(\d+)[IM]/g){
			$maplength += $1;
		}
		my $end = $start + $maplength - 1;
		foreach my $id(@{$bed_uorf{$as[2]}}){
			if(($start > $id->[0] and $start <= $id->[1] ) or ($end > $id->[0] and $end <= $id->[1])){
				my $ids = $as[2]."\t".$id->[0]."\t".$id->[1];
				$ct_uorf{$ids}++;
			}
		}
		foreach my $id(@{$bed_uorf_utr{$as[2]}}){
			if( $start > $id->[0] and $end <= $id->[1]){
				my $ids = $as[2]."\t".$id->[0]."\t".$id->[1];
				$ct_uorf_utr{$ids}++;
			}
		}
	}
	if(defined $bed_dorf{$as[2]}){
		my $start = $as[3];
		my $maplength = 0;
		while($as[5] =~ /(\d+)[IM]/g){
			$maplength += $1;
		}
		my $end = $start + $maplength - 1;
		foreach my $id(@{$bed_dorf{$as[2]}}){
			if(($start > $id->[0] and $start <= $id->[1] ) or ($end > $id->[0] and $end <= $id->[1])){
				my $ids = $as[2]."\t".$id->[0]."\t".$id->[1];
				$ct_dorf{$ids}++;
			}
		}
		foreach my $id(@{$bed_dorf_utr{$as[2]}}){
			if( $start > $id->[0] and $end <= $id->[1]){
				my $ids = $as[2]."\t".$id->[0]."\t".$id->[1];
				$ct_dorf_utr{$ids}++;
			}
		}
	}
	if(defined $bed_morf{$as[2]}){
		my $start = $as[3];
		my $maplength = 0;
		while($as[5] =~ /(\d+)[IM]/g){
			$maplength += $1;
		}
		my $end = $start + $maplength - 1;
		foreach my $id(@{$bed_morf{$as[2]}}){
			if(($start > $id->[0] and $start <= $id->[1] ) or ($end > $id->[0] and $end <= $id->[1])){
				my $ids = $as[2]."\t".$id->[0]."\t".$id->[1];
				$ct_morf{$ids}++;
			}
		}
		foreach my $id(@{$bed_morf_utr{$as[2]}}){
			if( $start > $id->[0] and $end <= $id->[1]){
				my $ids = $as[2]."\t".$id->[0]."\t".$id->[1];
				$ct_morf_utr{$ids}++;
			}
		}
	}

}
open OUT,">$ARGV[3].ORF.exp.xls" or die $!;
print OUT "id\torf_type\tsource_trans\torf_start\torf_end\torf_length\t3utr_start\t3utr_end\t3utr_length\torf_count\t3utr_count\torf_fpkm\t3utr_fpkm\tRRS\n";
open BED,$ARGV[0] or die $!;
my $id = 1;
while(<BED>){
	chomp;
	my @as = split("\t",$_);
	if($as[3] eq "uORF"){
		my $count_orf = defined $ct_uorf{join("\t",@as[0..2])} ? $ct_uorf{join("\t",@as[0..2])} : 0;
		my $count_utr = defined $ct_uorf_utr{$as[0]."\t".$as[4]."\t".$as[5]} ?  $ct_uorf_utr{$as[0]."\t".$as[4]."\t".$as[5]} :  0;
		my $fpkm_orf = sprintf("%.2f",$count_orf*10**6/($map_gene_reads*(($as[2]-$as[1])/1000)));
		my $fpkm_utr;
		my $rrs;
		if($count_utr != 0 and ($as[5]-$as[4])!=0){
			$fpkm_utr = sprintf("%.2f",$count_utr*10**6/($map_gene_reads*(($as[5]-$as[4])/1000)));
			if($fpkm_utr == 0){
				$rrs = $fpkm_orf;
			}else{
				$rrs = $fpkm_orf/$fpkm_utr;
			}
		}else{
			$fpkm_utr = 0.01;
			$rrs = $fpkm_orf;
		}
		$rrs = 0.01 unless (defined $rrs);
		$rrs = sprintf("%.2f",$rrs);
		$count_orf = 0 unless (defined $count_orf);
		$count_utr = 0 unless (defined $count_utr);
		$fpkm_orf = 0.00 unless (defined $fpkm_orf);
		$fpkm_utr = 0.00 unless (defined  $fpkm_utr);
		print OUT "ORF$id\t$as[3]\t$as[0]\t".($as[1]+1)."\t$as[2]\t".($as[2]-$as[1])."\t".($as[4]+1)."\t".$as[5]."\t".($as[5]-$as[4])."\t$count_orf\t$count_utr\t$fpkm_orf\t$fpkm_utr\t$rrs\n";
	}elsif($as[3] eq "dORF"){
		my $count_orf = defined $ct_dorf{join("\t",@as[0..2])} ? $ct_dorf{join("\t",@as[0..2])} : 0;
		my $count_utr = defined $ct_dorf_utr{$as[0]."\t".$as[4]."\t".$as[5]} ? $ct_dorf_utr{$as[0]."\t".$as[4]."\t".$as[5]} : 0; 
		my $fpkm_orf = sprintf("%.2f",$count_orf*10**6/($map_gene_reads*(($as[2]-$as[1])/1000)));
		my $fpkm_utr;
		my $rrs;
		if($count_utr != 0 and ($as[5]-$as[4])!=0){
			$fpkm_utr = sprintf("%.2f",$count_utr*10**6/($map_gene_reads*(($as[5]-$as[4])/1000)));
			if($fpkm_utr == 0){
				$rrs = $fpkm_orf;
			}else{
				$rrs = $fpkm_orf/$fpkm_utr;
			}
		}else{
			$fpkm_utr = 0.01;
			$rrs = $fpkm_orf;
		}
		$rrs = 0.01 unless (defined $rrs);
		$rrs = sprintf("%.2f",$rrs);
		$count_orf = 0 unless (defined $count_orf);
		$count_utr = 0 unless (defined $count_utr);
		$fpkm_orf = 0.00 unless (defined $fpkm_orf);
		$fpkm_utr = 0.00 unless (defined  $fpkm_utr);
		print OUT "ORF$id\t$as[3]\t$as[0]\t".($as[1]+1)."\t$as[2]\t".($as[2]-$as[1])."\t".($as[4]+1)."\t".$as[5]."\t".($as[5]-$as[4])."\t$count_orf\t$count_utr\t$fpkm_orf\t$fpkm_utr\t$rrs\n";
	}else{
		my $count_orf = defined $ct_morf{join("\t",@as[0..2])} ? $ct_morf{join("\t",@as[0..2])} : 0;
		my $count_utr = defined $ct_morf_utr{$as[0]."\t".$as[4]."\t".$as[5]} ? $ct_morf_utr{$as[0]."\t".$as[4]."\t".$as[5]} : 0;
		my $fpkm_orf = sprintf("%.2f",$count_orf*10**6/($map_gene_reads*(($as[2]-$as[1])/1000)));
		my $fpkm_utr;
		my $rrs;
		if($count_utr != 0 and ($as[5]-$as[4])!=0){
			$fpkm_utr = sprintf("%.2f",$count_utr*10**6/($map_gene_reads*(($as[5]-$as[4])/1000)));
			if($fpkm_utr == 0){
				$rrs = $fpkm_orf;
			}else{
				$rrs = $fpkm_orf/$fpkm_utr;
			}
		}else{
			$fpkm_utr = 0.01;
			$rrs = $fpkm_orf;
		}
		$rrs = 0.01 unless (defined $rrs);
		$rrs = sprintf("%.2f",$rrs);
		$count_orf = 0 unless (defined $count_orf);
		$count_utr = 0 unless (defined $count_utr);
		$fpkm_orf = 0.00 unless (defined $fpkm_orf);
		$fpkm_utr = 0.00 unless (defined  $fpkm_utr);
		print OUT "ORF$id\t$as[3]\t$as[0]\t".($as[1]+1)."\t$as[2]\t".($as[2]-$as[1])."\t".($as[4]+1)."\t".$as[5]."\t".($as[5]-$as[4])."\t$count_orf\t$count_utr\t$fpkm_orf\t$fpkm_utr\t$rrs\n";
	}
	$id++;
}
