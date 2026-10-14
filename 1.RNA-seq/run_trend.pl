#!/usr/bin/perl
use strict;
use warnings;
use FindBin qw($Bin $RealBin);
use Getopt::Long;
use Carp qw(carp confess croak);
use Data::Dumper;
use lib "$Bin/../";
use Support_Program;

usage() if(@ARGV < 2);

my ($diff_dir,$diff_glist,$conf_file,$dir_conf,$group_diff_dir);
my $type = "union"; # union or intersection
my $run = 1;
my ($pro_num,$ratio,$pvalue) = (20,2,0.05);
my ($id_list,$group,$exp_file,$ref_pfx,$outdir,$diff_list);
my $run_enrich = "yes";
my $pipe_type = "RNAseq";
my ($run_ko,$run_go,$run_do,$run_reactome) = (1,1,0,0);
my $mirna_target;
GetOptions( "diff_dir=s"=>\$diff_dir,
        "group_diff_dir=s" => \$group_diff_dir,
        "info_conf=s" => \$conf_file,
        "dirinfo=s" => \$dir_conf,
        "diff_list=s" => \$diff_list,
        "diff_glist=s"=>\$diff_glist,
        "id_list=s" => \$id_list,
        "group_file=s" => \$group,
        "exp=s" => \$exp_file,
        "ref_pfx=s" => \$ref_pfx,
        "outdir=s" => \$outdir,
        "enrich=s" => \$run_enrich,
        "run_ko=s" => \$run_ko,
        "run_go=s" => \$run_go,
        "run_do=s" => \$run_do,
        "run_reactome=s" => \$run_reactome,
        "type=s" => \$type,
        "pro=s"  => \$pro_num,
        "ratio=s" => \$ratio,
        "pavlue=s" => \$pvalue,
        "mirna_target|source_gene=s" => \$mirna_target,
        "pipe_type=s" => \$pipe_type,
        "run=s"=>\$run);

my $src = "$RealBin/pre_src";
my $get_glist = "perl $src/get_glist.pl";
my $id_get_glist = "perl $src/id_get_glist.pl";
my $get_exp = "perl $src/get_exp.pl";
my $trend = "perl $RealBin/Trend.pl";

mkdir $outdir if(!-e $outdir);
InitOpts(run=>$run);
$ratio = sprintf("%.4f",log($ratio)/log(2));

if(defined $conf_file){
    my %dir_info;
    my (%exp_file,%lnc_exp_file);

    if(defined $dir_conf){
        %dir_info = ReadDir($dir_conf);
        $group = "$dir_info{conf_file}/group.list" if(!defined $group);
        if(!defined $ref_pfx){
            $ref_pfx = exists $dir_info{result_annot} ? $dir_info{result_annot} : "$dir_info{annotation}/all";
        }
        if($pipe_type eq 'lncRNA'){
            %exp_file = GetExpFilebyConf(\%dir_info,'RNAseq');
            %lnc_exp_file = GetExpFilebyConf(\%dir_info,"lncRNA");
        }else{
            %exp_file = GetExpFilebyConf(\%dir_info,$pipe_type);
        }
        $mirna_target = $dir_info{result_mirna_target} if(exists $dir_info{result_mirna_target});
    }else{
        %exp_file = (sample => $exp_file, group=>$exp_file);
        $group_diff_dir //= $diff_dir;
        %dir_info = (group_diff => $group_diff_dir, sample_diff => $diff_dir);
    }

#    print Dumper(%dir_info);

    my %cc = qw/sample sample_diff   group group_diff/;
    if($pipe_type eq 'DIA' || $pipe_type eq 'itraq'){
        %cc = qw/sample diff group diff/;
    }
    open CONF,$conf_file or die $!;
    my $i = 0;
    my $no_exp = 0 if(!defined $exp_file);
    while(<CONF>){
        chomp;
        next if(/^#/);
        my ($sample_group,@samples) = split/\t/,$_;
        Error2("No [$sample_group] exp file") if(!exists $exp_file{$sample_group} || $exp_file{$sample_group} eq 'none');
        $exp_file = $exp_file{$sample_group};
        Error2("No dir info of [$cc{$sample_group}] in conf file [$dir_conf]") if(!exists $dir_info{$cc{$sample_group}});
        my $diff_dir = $dir_info{$cc{$sample_group}};
        $diff_dir = "$diff_dir/all_mirna" if($pipe_type eq 'smallRNA');
        TrendByConf($exp_file,$diff_dir,$outdir,$type,"trend",$pipe_type,$mirna_target,@samples);
        
        if(exists $dir_info{"$cc{$sample_group}_lnc"}){
            $diff_dir = $dir_info{"$cc{$sample_group}_lnc"};
            $exp_file = $lnc_exp_file{$sample_group};
            TrendByConf($exp_file,$diff_dir,$outdir,$type,"trend_lnc",$pipe_type,$mirna_target,@samples);
        }
    }
    close CONF;
}elsif(defined $diff_dir){ # get glist from diff dir
    if(defined $diff_list){
        Excu2("$get_glist $diff_list $diff_dir $outdir/gene.glist $type");
        print STDERR "No...\n";
    }
    else{
        Excu2("$id_get_glist $id_list $diff_dir $outdir/gene.glist $type");
#        Trend("$outdir/gene.glist");
        print STDERR "No...\n";
    }
}else{
    print STDERR "No...\n";
}

sub TrendByConf{
    my $exp_file = shift;
    my $diff_dir = shift;
    my $outdir = shift;
    my $type = shift;
    my $name = shift;  ## trend,  trend_lnc
    my $pipe_type = shift;
    my $mirna_target = shift;
    my @samples = @_;

    my $name2 = ucfirst($name);
    my $result_name = join("_","$name2",@samples);
    my $id_list = join(",",@samples);
    my $pipe = ($pipe_type eq 'lncRNA' && $name ne 'trend_lnc') ? "RNAseq" : $pipe_type;
#    my $pipe = $name eq 'trend_lnc' ? 'lncRNA' : "mRNA";

    my $name_tmp = join("_","$name",@samples);
    my $outdir_tmp = "$outdir/$name_tmp";#.join("_","trend",@samples);
    mkdir $outdir_tmp if(!-e $outdir_tmp);
#        CheckDir($outdir_tmp,"$outdir/bak","Trend$i");
    Excu2("$id_get_glist $id_list $diff_dir $outdir_tmp/gene.glist $type");   
    Trend("$outdir_tmp/gene.glist",$id_list,$exp_file,$outdir_tmp,pipe=>$pipe,target=>$mirna_target);
    Excu2("ln -sf --no-dereference $name_tmp/result $outdir/$result_name");

    return;
}

sub Trend{
    my $diff_glist = shift;
    my $id_list = shift;
    my $exp_file = shift;
    my $outdir = shift;
    my %opts = @_;

    $opts{pipe} //= "RNAseq";

    my $ref_pfx_tmp = $ref_pfx;
    my $desc = "$ref_pfx.Annot.txt";
    if($opts{pipe} eq 'lncRNA'){
        $ref_pfx_tmp = "none";
        $desc = "$ref_pfx.ts.symbol";
    }
    my %key_word_cc = qw/lncRNA lncRNA RNAseq Genes DIA Protein itraq protein smallRNA miRNAs 
                        Metabolome Metabolites circRNA circRNAs isoseq_RNAseq Genes
                        Denovo Genes/;
    Error2("No key_word  info of pipe [$opts{pipe}]") if(!exists $key_word_cc{$opts{pipe}});
    my $key_word = $key_word_cc{$opts{pipe}};
    my $target = defined $opts{target} ? "-t $opts{target}" : "";
    Error2("Group info is not defined, please add -group_file") if(!defined $group);

    Excu2("$get_exp $diff_glist $id_list $group $exp_file >$outdir/exp.xls");
    Excu2("$trend -i $outdir/exp.xls -o $outdir/  -key $key_word -ref $ref_pfx_tmp -desc $desc -pro $pro_num -p $pvalue -ratio $ratio $target -run_ko $run_ko -run_go $run_go -run_do $run_do -run_reactome $run_reactome -run $run");
    
    return;
}

sub CheckDir{
    my $dir = shift;
    my $bak_dir = shift;
    my $dir_name = shift;

    if(-e $dir){
        my $time = time();
        mkdir $bak_dir if(!-e $bak_dir);
        Excu2("mv $dir $bak_dir/$dir_name.$time");
        Excu2("mkdir $dir");
    }else{
        mkdir $dir;
    }

    return;
}

sub ReadDir{
    my $list = shift;

    open DIR,$list or die $!;
    my %hash;
    while(<DIR>){
        chomp;
        my @aa = split/\t/,$_;
        $hash{$aa[0]} = $aa[1];
    }
    close DIR;

    return %hash;
}

sub GetExpFilebyConf{
    my $conf_ref = shift;
    my $type = shift;

    $type //= "RNAseq";

    my %hash;
    if($type eq 'DIA'){
        $hash{sample} = $conf_ref->{result_all_protein} || "none";
        $hash{group} = $conf_ref->{result_all_group_exp} || "none";
    }elsif($type eq 'lncRNA'){
        $hash{sample} = $conf_ref->{result_lnc_exp} || "none";
        $hash{group} = $conf_ref->{result_lnc_group_exp} || "none";
    }elsif(exists $conf_ref->{result_all_exp}){
        $hash{sample} = $conf_ref->{result_all_exp};
        $hash{group} = $conf_ref->{result_all_group_exp} || $hash{sample};
    }else{
        my $exp_dir;
        if(exists $conf_ref->{RSEM}){
            $exp_dir = $conf_ref->{RSEM};
        }else{
            $exp_dir = $conf_ref->{expression};
        }
        $hash{sample} = "$exp_dir/all.genes.expression.annot.xls";
        $hash{group} = "$exp_dir/all.genes.group.expression.annot.xls";
    }

    return %hash;
}

sub GetFile{
    my $pfx = shift;
    my @adds = @_;

    my $out;
    for my $add(@adds){
        next if(!-e "$pfx$add");
        $out = "$pfx$add";
        last;
    }
    die "No file of $pfx [@adds]\n" if(!defined $out);

    return $out;
}

sub usage{
    
    my $usage =<<USAGE; 
    
    Usage: perl $0 -id_list <id_list> -group_file <group_file> -exp <exp_file> -ref_pfx <ref_pfx> -outdir <outdir>
          -conf_file -dir_conf -ref_pfx -outdir 
                
                -id_list     id list for trend
                -group_file  group_file for get mean
                -exp         exp file 
                -ref_pfx     ref pfx for enrich 
                -outdir      output dir

                -diff_dir    diff result dir  (*.filter.xls)
                -diff_list   diff list 
                
                -diff_glist  glist file

                -enrich      run enrich or not, [yes]/no

         trend  parameter 
                
                -pro         profile num [20]
                -ratio       ratio [2] 
                -pvalue      pavlue [0.05]
USAGE

    print "$usage\n";
    exit(1);
}
