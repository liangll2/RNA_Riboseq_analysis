#! /usr/bin/perl

use utf8;
use strict;
use warnings;
use Getopt::Long;
#use Bio::SeqIO;
use Bio::Seq;
#use List::Util qw/sum min max/;
#use List::MoreUtils qw/uniq/;
#use File::Basename qw/basename dirname/;
#use File::Spec::Functions qw/rel2abs/;
#use FindBin qw/$Bin $Script/;
#use lib $Bin;

my ($help,$chr,$bed,$gff,$target_fa,$pre,$nuc_length,$skip,$nuc_length_min);
GetOptions(
	"help"=>\$help,
	"chr:s" =>\$chr,
	"bed:s" =>\$bed,
	"gff:s" =>\$gff,
	"target_fa:s" =>\$target_fa,
	"pre:s" =>\$pre,
	"maxlength:s" =>\$nuc_length,
	"minlength:s" =>\$nuc_length_min,
	"skip:s" => \$skip
);
my $usage = "
	[usage]
	perl $0 -chr chr.fa -bed chr.bed -pre output.prefix -maxlength 300[default:300] -minlength 0[default:9]
	perl $0 -target_fa xxx.fa -pre output.prefix -maxlength 300
	perl $0 -help
	perl $0 -chr chr.fa -gff chr.gff -pre output.prefix -maxlength 300
	perl $0 -chr chr.fa -bed chr.bed -gff chr.gff -pre output.prefix -maxlength 300

	[parameter help]
	-chr         : fa format file which contains target sequence, must used with option -bed or -gff
	-bed         : record target sequece position messages,format: chr|start_loc|end_loc|name|score|strand
	-gff         : input gff/gtf file #MUST contain gene locs position 
	-target_fa   : directly input target sequence fa file; use this option;the -chr -bed -gtf option will be ignored;
	-pre         : output prefix
	-skip        : if your intput file is chr & gff,the program will default run intergenic,utr,intron sORF predict process
	               if you want to skip some process,please use this parameter : -skip intergenic,utr or -skip intron or etc.
	-help        : show this help message
";
die $usage if($help || !defined($pre));
my $mark_target = 0;
my $mark_gff = 0;
if(!defined($nuc_length)){
	$nuc_length = 300-3;
}else{
	$nuc_length = $nuc_length-3;
}
if(!defined($nuc_length_min)){
	$nuc_length_min = 9;
}
#print undef($skip);
if(!defined($skip)){
	$skip = "";
}
if(defined($target_fa)){
	$mark_target = 1;
	predict_orf($target_fa,$pre,$nuc_length_min);
	exit();
}elsif(defined($bed) && defined($chr)){
	`bedtools getfasta  -fi $chr -bed $bed -name -s > $pre.fa`;
	predict_orf("$pre.fa",$pre,$nuc_length_min);
}elsif(defined($gff) && defined($chr)){
#	open GENOME,$chr or die"cannot open:$!";
	open GFF,$gff or die "cant open gfffile;";
	$mark_gff = 1;
}else{
	die $usage;
}


if($mark_gff == 1){
	if($skip !~ /utr/){
		#############################################################
		########### utr predict ##########################
		#############################################################
		`gffread $gff -g $chr -w $pre.exon.fa -W`;
		`perl /home/miaoxin/Pipeline/Small_RNA/SmallRNA_Programs/ref/extract_utr_from_exon.pl $pre.exon.fa $pre`;
		predict_orf("$pre\_5UTR.fa","$pre\_5UTR",$nuc_length_min);
		predict_orf("$pre\_3UTR.fa","$pre\_3UTR",$nuc_length_min);
	}
if($skip !~ /intergenic/){
#############################################################
########### intergenic ORF predict ##########################
#############################################################

	my %fa;my $key;my %falength;
	if(-s "$chr.fai"){
		open FAI,"$chr.fai" or die "failed to open fai file";
		while(<FAI>){
			my @ar = split("\t",$_);
			$falength{$ar[0]} = $ar[1];
		}
	}else{
		`samtools faidx $chr`;
		open FAI,"$chr.fai" or die "failed to open fai file";
		while(<FAI>){
			my @ar = split("\t",$_);
			$falength{$ar[0]} = $ar[1];
		}
	}
	
	#begin to search gene locs;
	open GENLOC,">./geneloc.bed" or die "failed to write geneloc.bed files";
	while(<GFF>){
		next if /^#/;
		my @ar = split(/\t/,$_);
		if($ar[2] eq "gene"){
			print GENLOC $ar[0]."\t".($ar[3]-1)."\t".$ar[4]."\t\t\t$ar[6]\n"; 
		}
	}
	close GFF;
	`cat geneloc.bed|sort -k1,1 -k2,2n -k3,3n|bedtools merge -S + -i - > ./$pre.geneloc.merge.plus.bed`;
	`cat geneloc.bed|sort -k1,1 -k2,2n -k3,3n|bedtools merge -S - -i - > ./$pre.geneloc.merge.minus.bed`;
	#begin to get intergenic locs
	open BED,"$pre.geneloc.merge.plus.bed" or die "failed to open geneloc.merge.plus.bed";
	open INTERBED,">$pre.intergenic.plus.bed" or die "failed to write intergeic bed file";
	my $chr_flag="";
	my $cal=0;
	while(<BED>){
		chomp;
		my @a=split("\t",$_);
		if($chr_flag eq ""){
			#line1
			if($a[1] != 0){
				$cal++;
				print INTERBED $a[0]."\t0\t$a[1]\tPLUS_INTER$cal\t\t+\n";
			}
			$chr_flag = $a[0];
			print INTERBED $a[0]."\t$a[2]\t";
		}else{
			if($chr_flag eq $a[0]){
				$cal++;
				print INTERBED $a[1]."\tPLUS_INTER$cal\t\t+\n";
				print INTERBED $a[0]."\t".$a[2]."\t";
				$chr_flag = $a[0];
			}else{
				$cal++;
				print INTERBED $falength{$chr_flag}."\tPLUS_INTER$cal\t\t+\n";
				if($a[1] != 0){
					$cal++;
					print INTERBED $a[0]."\t0\t$a[1]\tPLUS_INTTER$cal\t\t+\n";
				}
				$chr_flag = $a[0];
				print INTERBED $a[0]."\t$a[2]\t";
			}
		}
	}
	$cal++;
	print INTERBED $falength{$chr_flag}."\tPLUS_INTER$cal\t\t+";
	close BED;
	close INTERBED;
	open BED,"$pre.geneloc.merge.minus.bed" or die "failed to open geneloc.merge.minus.bed";
	open INTERBED,">$pre.intergenic.minus.bed" or die "failed to write intergeic bed file";
	$chr_flag="";
	$cal = 0;
	while(<BED>){
		chomp;
		my @a=split("\t",$_);
		if($chr_flag eq ""){
			#line1
			if($a[1] != 0){
				$cal++;
				print INTERBED $a[0]."\t0\t$a[1]\tMINUS_INTER$cal\t\t-\n";
			}
			$chr_flag = $a[0];
			print INTERBED $a[0]."\t$a[2]\t";
		}else{
			if($chr_flag eq $a[0]){
				$cal++;
				print INTERBED $a[1]."\tMINUS_INTER$cal\t\t-\n";
				print INTERBED $a[0]."\t".$a[2]."\t";
				$chr_flag = $a[0];
			}else{
				$cal++;
				print INTERBED $falength{$chr_flag}."\tMINUS_INTER$cal\t\t-\n";
				if($a[1] != 0){
					$cal++;
					print INTERBED $a[0]."\t0\t$a[1]\tMINUS_INTER$cal\t\t-\n";
				}
				$chr_flag = $a[0];
				print INTERBED $a[0]."\t$a[2]\t";
			}
		}
	}
	$cal++;
	print INTERBED $falength{$chr_flag}."\tMINUS_INTER$cal\t\t-";
	close BED;
	close INTERBED;
	#begin to predict Intergenic ORF;
	#plus
	bed_to_fa("$pre.intergenic.plus.bed",$chr,"$pre.intergenic.plus.fa");
	predict_orf("$pre.intergenic.plus.fa","$pre.intergenic.plus",$nuc_length_min);
	#minus
	bed_to_fa("$pre.intergenic.minus.bed",$chr,"$pre.intergenic.minus.fa");
	predict_orf("$pre.intergenic.minus.fa","$pre.intergenic.minus",$nuc_length_min);
}
if($skip !~ /intron/){
#################################################
########### intron predict ######################
#################################################

	open GFF,$gff or die "failed to open gff file";
	my %exons;
	my %strand;
	if($gff =~ /\.gff$/){
		while(<GFF>){
			next if /^#/;
			chomp;
			my @ar = split("\t",$_);
			next if($ar[2] ne "exon");
			my $line = $_.";";
			my ($transcript) = $line =~ /transcript_id=(.*?);/;
			push (@{$exons{$transcript}},$ar[0]."\t".$ar[3]."\t".$ar[4]."\t".$ar[6]);
			$strand{$transcript} .= $ar[6];
		}
	}elsif($gff =~ /\.gtf$/){
		while(<GFF>){
			next if /^#/;
			chomp;
			my @ar = split("\t",$_);
			next if($ar[2] ne "exon");
			my ($transcript) = $_ =~ /transcript_id \"(.*?)\"/;
			push (@{$exons{$transcript}},$ar[0]."\t".$ar[3]."\t".$ar[4]."\t".$ar[6]);
			$strand{$transcript} .= $ar[6];
		}
	}
	open INTRON,">$pre.fakeintron.gtf" or die "failed to write fake intron gff file";
	foreach my $key(keys %exons){
		if(@{$exons{$key}}>1){
			print INTRON get_intron_gtf(\@{$exons{$key}},$key);
		}
	}
	#merge introns;
	`cut -f 1,4,5,7 $pre.fakeintron.gtf |sort -k1,1 -k2,2n -k3,3n|awk '{print \$1"\\t"\$2"\\t"\$3"\\t\\t\\t"\$4}' > ./$pre.intron.sort.bed`;
	`bedtools merge -i $pre.intron.sort.bed -s > $pre.intron.merge.temp`;
	`awk '{print \$1"\\t"\$2-1"\\t"\$3"\\tINTRON"NR"\\t\\t"\$4}' ./$pre.intron.merge.temp > ./$pre.intron.merge.bed`;
	bed_to_fa("$pre.intron.merge.bed",$chr,"$pre.intron.merge.fa");
	predict_orf("$pre.intron.merge.fa","$pre.intron",$nuc_length_min);
}
sub get_intron_gtf{
	my $index = $_[0];
	my $trans = $_[1];
	my @intron;
	my $s="";
	my $flag = 0;
	foreach my $v(sort{(split("\t",$a))[1] <=> (split("\t",$b))[1]} @{$index}){
		my @a = split("\t",$v);
		if($flag == 0){
			$s = $s."$a[0]\tFakeIntron\texon\t".($a[2]+1)."\t";
			$flag++;
		}else{
			$flag++;
			$s = $s.($a[1]-1)."\t.\t$a[3]\t.\ttranscript_id \"$trans\_intron\"\n";
			if($flag < @{$index}){
				$s = $s."$a[0]\tFakeIntron\texon\t".($a[2]+1)."\t";
			}
		}
	}
	return($s);
}


}



sub nc2pep{
	my $seq = $_[0];
	my $obj = Bio::Seq -> new(-seq => $seq);
	return $obj->translate->seq;
}

sub bed_to_fa{
	my $bedfile = $_[0];
	my $fafile = $_[1];
	my $outfile = $_[2];
	`bedtools getfasta -fi $fafile -bed $bedfile -s -name -fo $outfile`;
}

sub predict_orf{
	my $fafile = $_[0];
	my $prefix = $_[1];
	my $length = $_[2];
	my %fa;my $key;
	open INFA,$fafile or die "failed to open target fa file";
	while(<INFA>){
		chomp;
		if(/>/){
			s/>//;     my @a=split(/\s+/,$_);    $key=$a[0];
		}else{
			$fa{$key}.=uc($_);
		}
	}
	close INFA;
	my %cuts;
	#only find terminal sites
	open DETAIL,">$prefix.detail.xls" or die "can't open $prefix.detail.xls";
	print DETAIL "id\tstart\tend\tnuc\tpep\tpep_length\n";
	foreach my $chr(sort keys %fa){
		next if (length($fa{$chr}) <= 9);
		while($fa{$chr}=~m/(TAG|TAA|TGA)/ig){
			#absolute position
			my $match_end_pos = pos($fa{$chr});
			#getnear start locs
			$cuts{$chr}{$match_end_pos} = "$chr\t$match_end_pos\t$1";
			my $start = $match_end_pos-$nuc_length-3;
			if($start < 0){
				$start = 0;
			}
			my $gap = $match_end_pos-$start;
			my $seqs = substr($fa{$chr},$start,$match_end_pos-$start);
			for(my $i=length($seqs)-3;$i>=3;$i=$i-3){
				if(substr($seqs,$i-3,3) =~ /(TAG|TAA|TGA)/){
					last;
				}
				if(substr($seqs,$i-3,3) eq "ATG"){
					my $seqnuc = substr($fa{$chr},$i+($match_end_pos-$gap)-3,$match_end_pos-($i+($match_end_pos-$gap)-3));
					if(length($seqnuc) >= $length ){
						my $pep = nc2pep($seqnuc);
						print DETAIL $chr."\t".($i+($match_end_pos-$gap)-3)."\t".$match_end_pos."\t$seqnuc\t".$pep."\t".(length($pep)-1)."\n";
					}
				}
			}
		}
	}
}

