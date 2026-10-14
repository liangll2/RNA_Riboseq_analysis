#!/usr/bin/perl

use strict;
use warnings;
#use Getopt::Std;
use Getopt::Long;
use Cwd 'abs_path';
use FindBin qw($RealBin);
use List::Util qw/sum/;
use Tie::File;
use Data::Dumper;
use lib"$RealBin/lib";
use DEBUG2;
use Carp qw(carp confess croak);

my @ARGV_bak = @ARGV;

my %compare;
my %opts = qw/d no     g nogroup    p q    
              t 0.05   l 2          n no   m no
              soft  edgeR           key_name Genes          run 1
              run_go 1  run_ko 1  run_do 0  run_reactome 0  run_GSEA 0
              run_MF 1  run_BP 1  run_CC 1
              enrich_summary 0  
              nomap_die  1
              e no    r no    0 yes    x 0/;  ## e => enrich  r => report

#getopts('a:c:d:e:f:g:l:m:n:o:p:r:t:0:hxyz',\%opts);
GetOptions(\%opts,"a=s","c=s","d=s","e=s","f=s","g=s","l|fc=s","m=s","n=s","o=s","p|pq=s","r:s",
            "t|pq_cut=s","0=s","h","x","y","z","soft=s","Rscript=s","key_name=s",
            "enrich_summary=s","pair=s",
            "peptide=s",
            "run_MF=s","run_BP=s","run_CC=s",
            "nomap_die=s","source_gene=s",
            "run_go=s","run_ko=s","run_do=s","run_reactome=s","run_GSEA=s","run=s");

&usage if($opts{h});
&usage unless ($opts{f} || $opts{c});
$opts{0} = "no"  if($opts{run} == 0);

die "ERROR: your matrix file [$opts{f}] is not exists!\n" unless -e $opts{f};

# softwares
my $Rscript = `which Rscript`;chomp($Rscript);
$Rscript = $opts{Rscript} if(exists $opts{Rscript});


my $src = "$RealBin/src";
my $soft_adddesc = "perl $src/add_desc.pl";
my $g_annot_flag = 0;

my $pheatmap = "Rscript $src/pheatmap.r";
#my $enrich = "perl $Bin/../enrich/enrich.pl";
#my $enrich = "sh $Bin/src/enrich.sh";
my $enrich = "perl $src/enrich.pl";
my $edgeR = "$Rscript $src/edgeR.r";
my $DESeq2 = "$Rscript $src/DEseq2.r";
my $edgeR_stat = "perl $src/edgeR.stat.pl";
my $run_GSEA = "perl $src/run_GSEA.pl";
my $R_plot = "perl $src/R_plot.pl";
my $Ttest = "$Rscript $src/Ttest.r";
my $DIA = "perl $src/DIA_diff.pl";
my $snp = "perl $src/cal_snp_diff.pl";
my $get_source_gene = "perl $src/source_gene.pl";

## get error msg
my $g_exp_type = "fpkm";
my $msg;
my $DOWARN =1;
BEGIN { $SIG{__WARN__} = sub { $msg = $_[0];warn $_[0] if $DOWARN } };

## parameter deal
my $log = $opts{y}?0:1;##default yes
my $dispersion = $opts{d};
$opts{p} = lc($opts{p});
my $outdir = $opts{o} ? $opts{o} : 'edgeR' . time();
mkdir $outdir unless(-e $outdir);
mkdir "$outdir/enrich" if($opts{e} eq 'yes');

$opts{z} = 1;
## bak the cmd
RunCmd("diff.task","perl $0 @ARGV_bak",$opts{z});
openLOG("$outdir/out.log",\%opts,'diff') if($log);
LOG("CheckParameter");
&CheckParameter(\%opts) if ($log);

# read pair info 
my %pair_info = ReadPairInfo($opts{pair}) if(exists $opts{pair});

## filter duplicate row
LOG("Read Input files");
$opts{f} = FilterDuplicateRow($opts{f},$opts{o});

# read the column ids
my $colids = read_matrixF_header($opts{f},$opts{g});
#print Dumper(%$colids)."\n";

# get the col part data and the compare data 
_fetch_col_data($opts{g},$opts{c},\%compare);

## run the cmd
#my $Rscript = "/Bio/bin/Rscript";
my $err_flag = 0;
my $draw_flag = 0;
my %col = (p=>1,q=>2);
my @file_list;
LOG("Compare diffs");
#for my $name(keys %{$compare{'file'}}) {
for my $name(@{$compare{order}}){
    &LOG("Compare: $name");
    ## check matrix file, false return 1
    next if($opts{0} eq 'yes' && &CheckFileNumber($compare{'file'}{$name},1,1));
    my $ret = 0;
    my $pair_test = exists $opts{pair} ? "yes" : "no";
    if($opts{soft} eq 'edgeR'){
        if(($compare{'num'}{$name} eq '1,1' || $compare{'num'}{$name} =~ /^1,/ || $compare{'num'}{$name} =~ /,1$/) && $dispersion eq 'no'){
            $dispersion = 0.01;
        }
        $ret = Excu("$edgeR $compare{'file'}{$name} $compare{'num'}{$name} $compare{'name'}{$name}  $outdir/$name $dispersion $col{$opts{p}} $opts{t} $opts{l} $opts{n} $compare{'fpkm'}{$name} $pair_test");
    }elsif($opts{soft} eq 'DESeq2'){
        $ret = Excu("$DESeq2 $compare{'file'}{$name} $compare{'num'}{$name} $compare{'name'}{$name} $outdir/$name $col{$opts{p}} $opts{t} $opts{l} $pair_test");
    }elsif($opts{soft} eq 't-test'){
        $compare{'fpkm'}{$name} //= 'none';
        $ret = Excu("$Ttest $compare{'file'}{$name} $compare{'num'}{$name} $compare{'name'}{$name} $outdir/$name $col{$opts{p}} $opts{t} $opts{l} $pair_test $compare{'fpkm'}{$name}");
    }elsif($opts{soft} eq 'DIA'){
        $ret = Excu("$DIA -protein $opts{f} -peptide $opts{peptide} -group $opts{g} -diff_list $name -fc $opts{l} -pq_cut $opts{t} -pq $opts{p} -outdir $outdir");
    }elsif($opts{soft} eq 'snp'){
        $ret = Excu("$snp $compare{'fpkm'}{$name} $compare{'num'}{$name} $compare{'name'}{$name} $outdir/$name.all.xls");
        next;
    }elsif($opts{soft} eq 'OPLS-DA'){
        
    }else{
        LOG("Unknown soft type [$opts{soft}]\n");
        die "Unknown soft type [$opts{soft}]\n";
    }
    Excu("perl $src/combine_edgeR.pl $outdir/$name.all.xls $compare{'fpkm'}{$name} $compare{'num'}{$name} $opts{p} $opts{t} $opts{l} $opts{x}") if($opts{m} eq 'yes');

    if($ret != 0){
        $err_flag = 1;
        &Warn("There are some error in your input file! If you use Chinese characters or Greek letter, please instead of it. ");
        ## if failed, next
        next;
    }
    $draw_flag = 1;
    push @file_list,"$outdir/$name.filter.xls";

    ## run volcano
    if($opts{0} eq 'yes' || -e "$outdir/$name.all.xls"){
        open FILE,"$outdir/$name.all.xls" or die $!;
        my $head = <FILE>;
        close FILE;

#	my @cols = split/,/,$compare{'num'}{$name};
#	my $xcol = 1+$cols[0]*2+$cols[1]*2+2+1;
#	my $ycol = 1+$cols[0]*2+$cols[1]*2+2+1+$col{$opts{p}};
    
        my %pq_name = ('p'=>'PValue','q'=>'FDR');
        my $xcol = &GetColsNum($head,'log2(fc)');
        my $ycol = &GetColsNum($head,$pq_name{$opts{p}});

        Excu("$R_plot volcano -infile $outdir/$name.all.xls -xcol $xcol -ycol $ycol -xcut $opts{l} -ycut $opts{t} -title \"$name volcano plot\" -outpfx $outdir/$name.DE.volcano");
        if($opts{e} eq 'yes'){ ## enrich glist 
            if(exists $opts{source_gene}){
                Excu("$get_source_gene $outdir/$name.filter.xls $opts{source_gene} $outdir/enrich/$name.glist");
            }else{
                Excu("cut -f1,$xcol $outdir/$name.filter.xls | sed 1d >$outdir/enrich/$name.glist");
            }
        }
    }
    if($opts{r} eq 'yes' && $opts{m} eq 'yes'){
        Excu("$R_plot pheatmap -infile $outdir/$name.filter.xls -outdir $outdir -outpfx $name.heatmap");
    }
    if($opts{r} eq 'yes' && $opts{soft} eq 'DIA'){
        my $col = "1-".($compare{'allnum'}{$name}+1);
        Excu("$R_plot pheatmap -infile $outdir/$name.filter.xls -col $col -outdir $outdir -outpfx $name.heatmap ");
    }

    if($opts{run_GSEA} == 1){
        die "No ref info for GSEA\n" if(!exists $opts{a});;
        my $gen2sym = "$opts{a}.gen2sym";
        die "No ref [$gen2sym] for GSEA\n" if(!-e $gen2sym);
        Excu("$run_GSEA -diff $outdir/$name.all.xls -idmap $gen2sym -outdir $outdir/GSEA -outpfx $name");
    }
}

## add annot 
if($g_annot_flag == 1){
    Excu("$soft_adddesc $outdir/annot.xls 1 1 $outdir/*all.xls $outdir/*filter.xls");
}

## draw hist stat figure
my $draw = "$edgeR_stat ".join(",",@file_list)." $outdir/diff.stat $opts{key_name}";
Excu($draw);

## error cmd
ErrorCmd("error.task","perl $0 @ARGV_bak",$opts{z}) if($err_flag != 0 && $log);

## rm 
Excu("rm $opts{o}/tmp.xls") if(!$opts{x});
Excu("rm $opts{o}/*.fpkm.tmp") if(!$opts{x} && $opts{m} eq 'yes');
Excu("rm $opts{o}/*count.tmp") if(!$opts{x});
Excu("rm $opts{o}/annot.xls") if(!$opts{x} && $g_annot_flag);

if($opts{e} eq 'yes' && exists $opts{a}){
    my $ref_pfx = "$opts{a}";
    my $ko = "$ref_pfx.kopath";
    my $bgl = "$ref_pfx.bgl";
    my $go = "$ref_pfx";
    my $goannot = "$ref_pfx.annot";
    my $annot_txt = "$ref_pfx.Annot.txt";

    my $type = exists $opts{source_gene} ? "nodiff" : "diff";

    Excu("$enrich -i $outdir/enrich -ref $ref_pfx -type $type -run_go $opts{run_go} -run_ko $opts{run_ko} -run_do $opts{run_do} -run_reactome $opts{run_reactome} -run_MF $opts{run_MF} -run_BP $opts{run_BP} -run_CC $opts{run_CC} -summary $opts{enrich_summary} -nomap_die $opts{nomap_die}");
    
    if(-s $annot_txt && !exists $opts{source_gene}){
        Excu("perl $Bin/src/add_desc.pl $annot_txt 1 1 $outdir/*-vs-*all.xls $outdir/*-vs-*filter.xls");
    }
}



complete() if($log);
#print "GeneDenovo\n";

sub ReadPairInfo{
    my $file = shift;

    local($_);
    my %hash;
    open my $in_fh,$file or die $!;
    while(<$in_fh>){
        chomp;
        my @aa = split/\t/,$_;
        for my $i(0..$#aa){
            for my $y(0..$#aa){
                $hash{$aa[$i]}{$aa[$y]} = 1;
            }
        }
    }
    close $in_fh;
#    print Dumper(%hash);

    return %hash;
}


sub read_matrixF_header
{
    my $file = shift;
    my $sample_list = shift;
    my $cols;

    open IN,$file or die $!;
    
    my $header = <IN>; chomp $header;
    my @colids = split /\t/ , $header;
    shift @colids;
    %$cols = map { $colids[$_] => $_ + 2 } 0 .. $#colids;
    if($opts{m} eq 'yes'){
        my $annot_col = "1,";
        my $ii = 1;
        for my $i(@colids){
            $ii++;
            if($i =~ /_rpkm/){
                $g_exp_type = "rpkm";
            }
            if($i =~ /_rpm/){
                $g_exp_type = "rpm";
            }
            next if($i =~ /_count/ || $i=~ /_fpkm/ || $i =~ /_rpkm/ || $i =~ /_rpm/);
            $annot_col .= "$ii,";
        }
        if($annot_col ne "1,"){
            AwkPrint($opts{f},$annot_col,"$outdir/annot.xls");
            $g_annot_flag = 1;
        }
    }

    if($opts{soft} eq 'DIA'){
        open SAMPLE,$sample_list or die $!;
        my @samples;
        while(<SAMPLE>){
            chomp;
            my @aa = split/\t/,$_;
            push @samples,$aa[0];
        }
        close SAMPLE;
        my %sample = map{$_=>1} @samples;
        my $annot_col = "1,";
        my $nn = 1;
        for my $i(@colids){
            $nn++;
            next if(exists $sample{$i});
            $annot_col .= "$nn,";
        }
        if($annot_col ne "1,"){
            AwkPrint($opts{f},$annot_col,"$outdir/annot.xls");
            $g_annot_flag = 1;
        }
    }

    return $cols;
}

sub _fetch_col_data
{
    my $col_list = shift;
    my $compare = shift;
    my $hash = shift;

    my @temp;
    my @samples;
    my %samples;
    my %sample_type;
    my %style;
    my $group = 0;

    ## get the group column
    if(-s $col_list){
        open IN,$col_list or die $!;
        my %flag;
        while(<IN>)
        {
            chomp;
            chop $_ if ($_ =~ /\r$/);
            next if ($_ eq "" || $_ =~ /^#/);
            s/\s+$//;
            next if(/^\s*$/);
            my ($id,$group_id) = split/\t/,$_;##id && group
            next if(exists $flag{$group_id});  ## group err, next it 
            my $key = $id;
            if($opts{m} eq 'yes'){ ## os, id_count, id_fpkm
                $key = "$id\_count";
            }
            if(! $colids->{$key}){ ## check the sample id
                delete $$hash{'list'}{$group_id} if(exists $$hash{'list'}{$group_id});
                delete $$hash{'num'}{$group_id} if(exists $$hash{'num'}{$group_id});
                if($opts{m} eq 'yes'){
                    Warn("Can't find column name [$key] in the source file! Your file may be not from GeneDenovo! So filter group [$group_id]") if($log);
                }else{
                    Warn("There is not a column name [$id] in the source file! So filter group [$group_id]") if($log);
                }
                $flag{$group_id} = 1;
                next;
            }
            $hash->{count_col}{$id} = $colids->{$key};
            $$hash{'list'}{$group_id} .= "$colids->{$key},";  ## group member cols
            $$hash{name_list}{$group_id} .= "$id,";
            $$hash{'num'}{$group_id}++;  ## group member num
            if(!exists $hash->{list}{$id}){
                $hash->{list}{$id} = "$colids->{$key},";
                $$hash{name_list}{$id} .= "$id,";
                $$hash{'num'}{$id} = 1;
            }
            if($opts{m} eq 'yes'){
                my $kk = "$id\_$g_exp_type";
                if(! $colids->{$kk}){
                    stop("Your data is not genedenovo, please check your") if($log);
                }
                $$hash{'list2'}{$group_id} .= "$colids->{$kk},";
                $hash->{exp_col}{$id} = $colids->{$kk};
                if(!exists $hash->{list2}{$id}){
                    $$hash{'list2'}{$id} .= "$colids->{$kk},";
                }
            }
        }
        close IN;
    }else{
        ## no group info 
        #  group is always exist
        for my $id (sort keys %{$colids}){
            my $name = $id;
            if($id =~ /_count/){
                $name =~ s/_count//;
                $hash->{list}{$name} = "$colids->{$id},";
                $hash->{num}{$name} = 1;
            }
            if($id =~ /_fpkm/ || $id =~ /_rpkm/ || $id =~ /_rpm/){
                $name =~ s/_fpkm$//;
                $name =~ s/_rpkm$//;
                $name =~ s/_rpm$//;
                $hash->{list2}{$name} = "$colids->{$id},";
            }
            $hash->{list}{$id} = "$colids->{$id},";
            $hash->{num}{$id} = 1;
        }
    }
#    print Dumper($hash);

    ## get the compare file and compare command
    open IN2, $compare or die $!;
    while(<IN2>){
        chomp;
        chop $_ if ($_ =~ /\r$/);
        next if ($_ eq "" || $_ =~ /^\s*#/);
        s/\s+$//;
        my @arr = split/\t/,$_;##control && treat
        my $name = join("-vs-",@arr);
        my $ok = 1;
        for my $n(@arr){
            if(! $$hash{'list'}{$n}){
                Warn("There is not a group named [$n]! So filter the compare [$name]") if($log);
                $ok = 0;
            }
        }
        next if($ok == 0);
        
        my $col_list = join("",map{$hash->{list}{$_}}@arr);
        my $num_list = join(",",map{$hash->{num}{$_}}@arr);
        ($col_list,$num_list) = GetPairColList($hash,\%pair_info,@arr,'count_col') if(exists $opts{pair});
        AwkPrint($opts{f},"1,$col_list","$outdir/$name.count.tmp");
        if($opts{m} eq 'yes'){
            for my $n(@arr){
                Error("No info of [$n]") if(!exists $$hash{'list2'}{$n});
            }
            my $col_list = join("",map{$hash->{list2}{$_}}@arr);
            ($col_list) =  GetPairColList($hash,\%pair_info,@arr,'exp_col') if(exists $opts{pair});
            AwkPrint($opts{f},"1,$col_list","$outdir/$name.fpkm.tmp");
        }
        
        $$hash{'file'}{$name} = "$outdir/$name.count.tmp";
        $$hash{'fpkm'}{$name} = 'none';
        $$hash{'fpkm'}{$name} = "$outdir/$name.fpkm.tmp" if($opts{m} eq 'yes');
        $$hash{'num'}{$name} = "$num_list";
        $$hash{'allnum'}{$name} = $$hash{'num'}{$arr[0]}+$$hash{'num'}{$arr[1]};
        $$hash{'name'}{$name} = "$arr[0],$arr[1]";
        push @{$$hash{order}},$name;
    }
    close IN2;
}

sub GetPairColList{
    my $hash_ref = shift;
    my $pair_info_ref = shift;
    my $control = shift;
    my $treat = shift;
    my $key = shift;

    $key //= 'count_col';

    my @control_name = split/,/,$hash_ref->{name_list}{$control};
    my @treat_name = split/,/,$hash_ref->{name_list}{$treat};
    my %treat_list = map{$_=>1}@treat_name;

    my @control_col;
    my @treat_col;
#    print Dumper($pair_info_ref);
#    print Dumper($hash_ref);
    for my $i(0..$#control_name){
        for my $y(0..$#treat_name){
#            print "$control_name[$i]\t$treat_name[$y]\n";
            ## find pair  && get pair col
            if(exists $pair_info_ref->{$control_name[$i]}{$treat_name[$y]}){
                die "No col info of [$control_name[$i]] [$key]" if(!exists $hash_ref->{$key}{$control_name[$i]});
                push @control_col,$hash_ref->{$key}{$control_name[$i]};
                push @treat_col,$hash_ref->{$key}{$treat_name[$y]};
                last;
            }
        }
    }

    my $list = join(",",@control_col,@treat_col);
    my $num_list = join(",",scalar(@control_col),scalar(@treat_col));

    print "@control_name\t@control_col\t=\t@treat_name\t@treat_col\t$list\t$num_list\n";

    return ("$list,",$num_list);
}

sub array_sum{
        my $array = shift;
        my $line = shift;
        my $sum;

        my $num = @{$array} -1;
#		print "$num:@{$array}\n";
        for (0..$num){
            my $i = $_+1;
            next if ($$array[$_] eq 'NA');
            if($$array[$_] eq ""){
                $$array[$_] = 0;
                Warn("In line $line column $i: there is value not defined, use 0 instead!");
            }
            if(isNumber($$array[$_])){
                &stop("In line $line column $i: [$$array[$_]] is not numeric, please check your file");
            }
            $sum += $$array[$_];
        }
        return $sum;
}

sub isNumber{ # yes => 0,  no => 1
    my $out = 0;
        
    my $i = 1;
    for(@_){
        next if($_ eq 'Filtered' || $_ eq 'NA');
        my $a = $_ + 0;
        if($_ =~ /inf/i || (defined $msg && $msg =~ /isn't numeric/)){
            $out = $i;
            last;
        }
        $i++;
    }
    return $out;
}

sub CheckFileNumber{## file, colname, head
    my $file = shift;
    my $col = shift;
    my $head = shift;

    my $num;
#    print "$file\n";
    open FILE,$file or die $!;
    while(<FILE>){
        chomp;
        my @aa = split/\t/,$_;
        $num = @aa if($. == 1);
        my $num2 = @aa;
        if($num2 != $num){
            &Warn("There is $num columns in the first line, but only [$num2] in line $.");
            return 1;
        }
        next if($. == 1 && $head == 1);
        shift @aa if($col ==1);
        my $err = &isNumber(@aa);
        if($err){
            &Warn("In line $. column $err, value [$aa[$err-1]] is not numeric.");
            return 1;
        }
    }
    close FILE;
    return 0;
}

sub CheckParameter{
    my $hash = shift;
    
    my $msg = "please check your file!\nIf you save the file from excel, please confirm what you save is not =>unicode txt<=";
    if(-s $$hash{f} && -B $$hash{f}){
        &stop("Input count file is not txt file, $msg");
    }
    if(-s $$hash{g} && -B $$hash{g}){
        &stop("Group file is not txt file, $msg");
    }
    if(-s $$hash{c} && -B $$hash{c}){
        &stop("Compare file is not txt file, $msg");
    }

    if($$hash{d} ne 'no' && &isNumber($$hash{d})){
        &stop("Parameter dispersion must be numeric or no, [$$hash{d}]");
    }
    if($$hash{p} ne 'p' && $$hash{p} ne 'q'){
        &stop("Parameter test value must be p or q, [$$hash{p}]");
    }
    if(&isNumber($$hash{t})){
        &stop("Parameter  pqcutoff must be numeric, [$$hash{t}]");
    }
    if(&isNumber($$hash{l}) || $$hash{l}<0){
        &stop("Parameter log2fc cutoff must be numeric and must bigger than 0, [$$hash{l}]");
    }
    
    my $log = '';
    $log .= "Data type(genedenovo): $opts{m}\n";
    $log .= "Group Info:\n".ReadFile($hash->{g})."\n" if($hash->{g} ne 'nogroup');
    $log .= "Compare Info:\n".ReadFile($hash->{c})."\n";
    $log .= "log2fc cutoff: $hash->{l}\n";
    $log .= "pqcut: $hash->{p} $hash->{t}\n";
    $log .= "edgeR dispersion : $hash->{d}\n" if($hash->{soft} eq 'edgeR');
    $log .= "Use edgeR NormFactors: $hash->{n}\n" if($hash->{soft} eq 'edgeR');
    LOG($log);
    
    return;
}

sub ReadFile{
    my $file = shift;
    open FILE,$file or die $!;
    my @aa = <FILE>;
    close FILE;

    return join("",@aa);
}

sub GetColsNum{
    my $line = shift;
    my $key = shift;

    chomp($line);
    my @cols = split/\t/,$line;
    my $out = 0;
    for my $i(0...$#cols){
#        print "$cols[$i]\n";
        if($cols[$i] =~ /\Q$key\E/i){
            return $i+1;
        }
    }
    print STDERR"Can't found key[$key] in the file head, please check your file!\n";
    exit(0);
}

sub FilterDuplicateRow{
    my $file = shift;
    my $outdir = shift;

    my %hash;
    open FILE,$file or die $!;
    my $out  = '';
    my $flag = 0;
    while(<FILE>){
        chomp;
        my @aa = split/\t/,$_;
        if(exists $hash{$aa[0]}){
            $flag = 1;
            next;
        }
        $hash{$aa[0]} = 1;
        $out .= "$_\n";
    }
    close FILE;
    &Warn("There are some same rownames in your file and will filter it, please check your file!") if($flag);
    open OUT,">","$outdir/tmp.xls" or die $!;
    print OUT"$out";
    close OUT;

    return "$outdir/tmp.xls";
}

sub Excu{
    my $cmd = shift;

    print "$cmd\n";
    my $ret = 0;
    if($opts{0} eq 'yes'){
        $ret = system($cmd);
        if($ret != 0){
            Error("ret$ret\t$cmd");
        }
    }

    return $ret;
}

sub Error{

    confess("Error:@_\n");
    return;
}

sub AwkPrint{
    my $infile = shift;
    my $cols = shift;
    my $outfile = shift;

    chop($cols);
    my @cols = split/,/,$cols;
    $cols =~ s/,/,\$/g;
    my $cmd = "awk -F\$'\\t' -vOFS='\\t' '{print \$$cols}' $infile >$outfile\n";
#    print "$cmd\n";
    Excu($cmd) ;

    return;
}

# help info.
sub usage
{
    print <<USAGE;
Version: 0.3
Date:    2014-04-15
Contact: pai\@genedenovo.com

Usage:   perl $0 [options]

Options: -f INT    the count matrix table file, separated by <TAB>
         -g INT    group info, format: sampleid\\tgroupid
         -c INT    compare list, format: group 1\\tgroup 2
         -d INT    dispersion, experiment:0.4, model organism: 0.1, technical repeat: 0.01
         -o INT    outpur dir
         -p INT    pvalue or fdr p/q [q]
         -t INT    pvalue or fdr cutoff [0.05]
         -l INT    FC cutoff [2],  will do log2(fc) in r

         -n STR    use NormFactors or not [no]/yes
         -m STR    data from genedenovo or not [no]/yes
         -r STR    report or not [no]/yes  ## run pheatmap
         -x        test, no rm

         -e STR    run enrich or not [no]/yes
         -a STR    enrich ref [go=>.annot, ko=>.kopath, GSEA=>.gen2sym, reactome=>.reactome, do=>.do]
         -run_go, -run_ko, -run_GSEA, -run_reactome -run_do 

         -pair     pair info, one pair a line,  A1<tab>B1<tab>C1
                                                A2<tab>B2<tab>C2


         -soft     [edgeR],DEseq2 
         -Rscript  Rscript path, DEseq2 use Rscript3

Example: perl $0 -f all.count.xls -c column.lis -r row.lis  -o <outdir> -p p -t 0.05 -l 1

     ...
USAGE
    exit;
}
