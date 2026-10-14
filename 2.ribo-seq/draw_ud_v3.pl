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
die "perl $0 <dir/*.bam> <gtf>\n" unless(@ARGV eq 2);
my $basename = basename($ARGV[0]);
my $path = $ARGV[0];
$path =~ s/$basename$//;
$basename =~ s/\.bam//g;
`mkdir -p $path/$basename`;
`ln -s $ARGV[0] $path/$basename/$basename.bam`;
print ("source /state/partition1/WORK/lib-devel/working_library/anaconda2/bin/activate /state/partition1/WORK/lib-devel/working_library/anaconda2/envs/R363;/state/partition1/WORK/lib-devel/working_library/anaconda2/envs/R363/bin//Rscript /Bio/Bin/pipeline/Ribo_seq/v6.0/Script/do.r $basename $path $ARGV[1]");
`source /state/partition1/WORK/lib-devel/working_library/anaconda2/bin/activate /state/partition1/WORK/lib-devel/working_library/anaconda2/envs/R363;/state/partition1/WORK/lib-devel/working_library/anaconda2/envs/R363/bin/Rscript /Bio/Bin/pipeline/Ribo_seq/v6.0/Script/do.r $basename $path $ARGV[1]`;

