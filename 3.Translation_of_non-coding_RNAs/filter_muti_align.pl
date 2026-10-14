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
die "perl $0 <gene_trans.map> <sort.byreas.bam> <output_prefix>\n" unless(@ARGV eq 3);
open FILE,$ARGV[0] or die "failed to open file1";
my %tran2gene;
while(<FILE>){
	chomp;
	my @as = split("\t",$_);
	$tran2gene{$as[1]} = $as[0];
}
close FILE;
open BAM,"samtools view -h $ARGV[1]|" or die "faile to read $ARGV[1]";
open OUT,"|samtools view -Sb  - > $ARGV[2].filter.bam" or die "failed to write $ARGV[2].filter.bam";
my $flag="";
my $raw=0;
my $left = 0;
my %genehash=();
while(<BAM>){
	if(/^@/){
		print OUT $_;
	}else{
		my @as = split("\t",$_);
		if($flag eq ""){
			push(@{$genehash{$tran2gene{$as[2]}}},$_);
			$flag = $as[0];
		}elsif($flag eq $as[0]){
			if(scalar(keys %genehash)>1){
				next;
			}
			push(@{$genehash{$tran2gene{$as[2]}}},$_);
		}else{
			if(scalar(keys %genehash) == 1){
				my ($id) = keys %genehash;
				foreach my $order(@{$genehash{$id}}){
					print OUT $order;
				}
				$left++;
			}
			$raw++;
			%genehash=();
			push(@{$genehash{$tran2gene{$as[2]}}},$_);
			$flag = $as[0];
		}
	}
}
if(scalar(keys %genehash) == 1){
	my ($id) = keys %genehash;
	foreach my $order(@{$genehash{$id}}){
		print OUT $order;
	}
	$left++;
}

$raw++;
close BAM;
open STAT,">$ARGV[2].stat.xls" or die "failed to write $ARGV[2].stat.xls";
print STAT "ALL\tfiltered\n$raw\t$left";
close STAT;
