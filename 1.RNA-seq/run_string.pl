#!/usr/bin/env perl 

use strict;
use warnings;
use FindBin qw/$RealBin $Script/;
#use File::Basename qw/basename dirname/;
use Getopt::Long;
#use lib $Bin;
use Data::Dumper;
use lib "$RealBin/..";
use Cwd 'abs_path';
use Carp qw(carp confess croak);
use utf8;
use Support_Program;

my $version = "v1.0";
#my $db = "$RealBin/db";
my $src = "$RealBin/src";

my $usage =<<USAGE;

    Desc: use for xxx

    version: $version

    perl $0 <> 
        db 
        string 
        combineidmap


    e.g: perl $0 <> <>
         perl $0 db -infile hsa.gen2sym -idmap ../hsa.id_mapping_file -outpfx hsa 
         perl $0 db -idmap hsa.ensemble.idmap -outpfx hsa
         perl $0 db -pep t.fa -class animal  -outpfx hsa
         perl $0 db -pep t.fa -class animal  -outpfx hsa -ref_sp hsa.string.sp
USAGE

die "$usage" if(@ARGV<1); 
my $type = shift;
my %opts = ("run"=>1,"type"=>"gene","cytoscape"=>1);
GetOptions(\%opts,"infile=s","idmap=s","outpfx=s","class=s","db=s","pep=s","newgene=s","run=s","symbol=s","desc=s","ref_sp=s","ref_idmap=s","new_idmap=s","outfile=s","dna=s","type=s","cytoscape=s");

$opts{idmap} = abs_path($opts{idmap}) if(exists $opts{idmap});
$opts{infile} = abs_path($opts{infile}) if(exists $opts{infile});

#my $diamond = "/Bio/Linux-src-files/NGS/Align/diamond/diamond-0.9.22/diamond";
my %apps; 
#Support_Program("$RealBin/../",\%apps,1);
Support_Program_conf("$RealBin/../",\%apps,"$RealBin/../Small_pipe.soft.conf");
my $diamond = $apps{diamond};
my %conf = Conf($apps{string_db}); ## db 

## db 
if($type eq 'db'){
    FilePrepare(%opts);
}

if($type eq 'string'){
    String(%opts);
}

if($type eq 'combineidmap'){
    CombineIdmap(%opts);
}


sub FilePrepare{
    my %opts = @_;
    
    ## check db 
    my $sp;
    if(exists $opts{newgene}){ ## 新基因，直接通过pep查找
        $sp = MapDb(%opts,'class'=>'all');  
        exit(0);
    }elsif(defined $opts{idmap} && -s $opts{idmap}){ ##通过idmap 文件查找
        my $col = GetCol($opts{idmap},0,"ProteinName"); ## col //= 2
        $col //= 2;
        $sp = CheckDB($opts{idmap},$col);
        if($sp){
            my $id_key = $opts{type} eq 'gene' ? "GeneID" : "TranscriptID";
            my $id_col = GetCol($opts{idmap},1,"$id_key");
            Excu("awk '{print \$$id_col\"\\t\"\$$col}' $opts{idmap} >$opts{outpfx}.string.idmap");
#            Excu("ln -sf $opts{idmap} $opts{outpfx}.string.idmap");
        }else{
            print "idmap can't find sp\n";
        }
    }elsif(exists $opts{infile}){ ## 通过infile id 去找
        $sp = CheckDB($opts{infile},1);
        if($sp){
            Excu("ln -sf $opts{infile} $opts{outpfx}.string.idmap");
        }
    }else{ ## 

    }

    if($sp){
        print "found [$sp]\n";
    }else{ ## class=> sp , pep序列去找 对应关系
        $sp = MapDb(%opts);
    }
    Excu("ln -sf $conf{sp_dir}/$sp.* $opts{outpfx}.string.db");
    if(defined $opts{ref_sp}){
        Excu("cp $opts{ref_sp} $opts{outpfx}.string.sp");# if($opts{ref_sp} ne "$opts{outpfx}.string.sp");
    }else{
        open OUT,">$opts{outpfx}.string.sp" or die $!;
        print OUT "$sp\n";
        close OUT;
    }

    return;
}

sub MapDb{
    my %opts = @_;

    die "No db info of [$opts{class}]\n" if(!exists $conf{$opts{class}});
    if(exists $opts{newgene}){
        $opts{dna} = $opts{pep};
    }
    if(exists $opts{dna}){ ## dna 
        Excu("$diamond blastx -k 0 -f 6 -p 4 -e 1e-5 --salltitles -q $opts{dna} -d $conf{$opts{class}} -o $opts{outpfx}.blast.out --tmpdir /tmp");       
    }else{
        Excu("$diamond blastp -k 0 -f 6 -p 4 -e 1e-5 --salltitles -q $opts{pep} -d $conf{$opts{class}} -o $opts{outpfx}.blast.out --tmpdir /tmp");       
    }
    my $ref_sp = exists $opts{ref_sp} ? "$opts{ref_sp}" : "";
    Excu("perl $src/blast2idmap.pl $opts{outpfx}.blast.out $opts{outpfx}.string.idmap $ref_sp");

    return $conf{name}{$opts{class}};
}

sub String{
    my %opts = @_;

#    my %db = ReadDB($opts{db});
    my %idmap = ReadIdmap($opts{idmap},0);
    my %symbol = ReadSymbol($opts{symbol}) if(exists $opts{symbol});
    my %annot = ReadAnnot($opts{desc}) if(exists $opts{desc});

    open IN,$opts{infile} or die "$!:$opts{infile}";
    `rm -f $opts{outpfx}.node.txt` if(-l "$opts{outpfx}.node.txt");
    open NODE,">$opts{outpfx}.node.txt" or die $!;
    print NODE "id\ttype".AddAnnot('head',\%annot)."\n";
    if(exists $opts{symbol}){
        open NODES,">$opts{outpfx}.node.symbol.txt" or die $!;
        print NODES "id\ttype".AddAnnot('head',\%annot)."\n";
    }
    my %idlist;
    my %e_node;
    while(<IN>){
        chomp;
        my @aa = split;
        next if($aa[0] eq 'id');
        $idlist{$aa[0]} = 1;
        my $id = $aa[0];
        if(!exists $e_node{$id}){
            print NODE "$id\tgene".AddAnnot($id,\%annot)."\n";
            if(SymbolOK($id,\%symbol)){
                print NODES "$symbol{$id}\tgene".AddAnnot($id,\%annot)."\n";
            }
        }else{
            $e_node{$id} = 1;
        }
    }
    close IN;

    open OUT,">$opts{outpfx}.string.txt" or die $!;
    print OUT "protein1\tprotein2\tcombined_score\n";
    if(exists $opts{symbol}){
        open OUTS,">$opts{outpfx}.string.symbol.txt" or die $!;
        print OUTS "protein1\tprotein2\tcombined_score\n";
    }
    my %e;
    open DB,$opts{db} or die "$!:$opts{db}";
    <DB>;
    while(<DB>){
        chomp;
        my ($source,$target,$score) = split/\t/,$_;
        $source =~ s/^\d+\.//;
        $target =~ s/^\d+\.//;
        next if(exists $e{$source}{$target} || exists $e{$target}{$source});
        $e{$source}{$target} = 1;

        my @source = exists $idmap{$source} ? @{$idmap{$source}} : ($source);
        my @target = exists $idmap{$target} ? @{$idmap{$target}} : ($target);

        for my $source0(@source){
            next if(!exists $idlist{$source0});
            for my $target0 (@target){
                next if(!exists $idlist{$target0});
                print OUT "$source0\t$target0\t$score\n";
                if(SymbolOK($source0,\%symbol) && SymbolOK($target0,\%symbol)){
                    print OUTS "$symbol{$source0}\t$symbol{$target0}\t$score\n";
                }
            }
        }
    }
    close OUT;
    close NODE;
    if(exists $opts{symbol}){
        close OUTS;
        close NODES;
    }

    my $outdir;
    my $pfx;
    if($opts{outpfx} =~ /(.*)\/(.*)/){
        $outdir = $1;
        $pfx = $2;
    }else{
        $outdir = "./";
        $pfx = $opts{outpfx};
    }

    if($opts{cytoscape} == 1){
        Excu("perl $src/make_cytoscape.pl $opts{outpfx}.string.txt $opts{outpfx}.node.txt $outdir $pfx euler");
        Excu("perl $src/make_cytoscape.pl $opts{outpfx}.string.symbol.txt $opts{outpfx}.node.symbol.txt $outdir $pfx.symbol euler") if(exists $opts{symbol});
    }

    return;
}

sub CombineIdmap{
    my %opts = @_;

    my %idmap = ReadIdmap($opts{ref_idmap},0);

    open OUT,">$opts{outfile}" or die $!;
    print OUT "GeneID\tProteinName\n";
    for my $id(sort keys %idmap){
        for my $source_id(@{$idmap{$id}}){
            next if($source_id eq 'GeneID');
            print OUT "$source_id\t$id\n";
        }
    }
    open NEW,$opts{new_idmap} or die $!;
    while(<NEW>){
        next if(/GeneID/);
        print OUT "$_";
    }
    close NEW;
    close OUT;

    return;
}

sub SymbolOK{
    my $id = shift;
    my $symbol_ref = shift;
    
#    print Dumper($symbol_ref);die;
    if(exists $symbol_ref->{$id} && $symbol_ref->{$id} ne '--' && $symbol_ref->{$id} ne '-'){
        return 1;
    }
    
    return 0;
}

sub ReadDB{
    my $file = shift;
    
    local($_);
    my %hash;
    open FILE,$file or die "$!:$file\n";
    while(<FILE>){
        chomp;
        my @aa = split;
        $aa[0] =~ s/^\d+\.//;
        $aa[1] =~ s/^\d+\.//;
        $hash{$aa[0]}{$aa[1]} = $aa[2];
    }
    close FILE;

    return  %hash;
}

sub ReadSymbol{
    my $file = shift;
    local($_);
    my %hash;
    open FILE,$file or die "$!:$file\n";
    while(<FILE>){
        chomp;
        my @aa = split/\t/,$_;
        next if(@aa < 2);
        $aa[1] =~ s/\'//g;
        $hash{$aa[0]} = $aa[1];
    }
    close FILE;

    return %hash;
}

sub ReadIdmap{
    my $file = shift;
    my $strict = shift;
    $strict //= 1;

    local($_);
    my %hash;
    my $col1 = GetCol($file,$strict,'GeneID');
    my $col2 = GetCol($file,$strict,'ProteinName');
    $col1 //= 1;
    $col2 //= 2;
    $col1 -= 1;
    $col2 -= 1;

    open FILE,$file or die "$!:$file\n";
    while(<FILE>){
        chomp;
        my @aa = split;
        if(@aa > 1){
            $aa[$col2] =~ s/^\d+\.//;
            push @{$hash{$aa[$col2]}},$aa[$col1];
        }else{
            @{$hash{$aa[$col1]}} = ($aa[$col1]);
        }
    }
    close FILE;
#    print Dumper(%hash);

    return %hash;
}


sub CheckDB{
    my $file = shift;
    my $col = shift;

    $col //= 1;
    $col -= 1;

    open NAME,$conf{name_list} or die "$!: ,$conf{name_list}\n";
    my %name_list;
    while(<NAME>){
        chomp;
        my @aa = split;
        my ($sp,$gene) = split/\./,$aa[0],2;
        $name_list{$gene} = $sp;
    }
    close NAME;
    
    open FILE,$file or die $!;
    while(<FILE>){
        chomp;
        my @aa = split;
        $aa[$col] =~ s/\..*$//;## protein version
        if(exists $name_list{$aa[$col]}){
            return $name_list{$aa[$col]};
        }
    }
    
    return;
}

sub GetCol{
    my $file = shift;
    my $strict = shift;
    my $name = shift;

    $strict //= 1;

    open  FILE,$file or die "$!:$file\n";
    my $head = <FILE>;
    chomp($head);
    my @head = split/\t/,$head;
    for my $i(0..$#head){
        if($head[$i] eq $name){
            return $i+1;
        }
    }
    close FILE;
    Error("Can't find name [$name] in [$head] of [$file") if($strict == 1);

    return;
}

sub ReadAnnot{
    my $annot = shift;

    my %hash;
    local($_);
    open ANNOT,$annot or die "$!:$annot";
    while(<ANNOT>){
        chomp;
        my @aa = split/\t/,$_;
        my $id = shift @aa;
        if($. == 1){
            $id = "head";
            $hash{col_num} = scalar(@aa);
        }
        $hash{$id} = join("\t",@aa);
    }
    close ANNOT;

    return %hash;
}

sub AddAnnot{
    my $col_id = shift;
    my $annot_ref = shift;;

    my $out;
    if(exists $annot_ref->{col_num}){
        $annot_ref->{$col_id} //= join("\t",("-")x$annot_ref->{col_num});
        $out = "\t$annot_ref->{$col_id}";    
    }else{
        $out = "";
    }
    return $out;
}

sub Conf{
    my $dir = shift;

    my %hash;
    $hash{name_list} = "$dir/all.10.names.xls";
    $hash{sp_dir} = "$dir/String-db/protein.links.sorted.by.species/";
    $hash{animal} = "$dir/9606.dmnd";
    $hash{plant} = "$dir/3702.dmnd";
    $hash{fungi} = "$dir/4932.dmnd";
    $hash{mirco} = "$dir/511145.dmnd";
    $hash{all} = "$dir/all.dmnd";
    $hash{name}{animal} = "9606";
    $hash{name}{plant} = "3702";
    $hash{name}{fungi} = "4932";
    $hash{name}{micro} = "511145";
    $hash{name}{all} = "all";

    $hash{vertebrates} = $hash{animal};
    $hash{invertebrates} = $hash{animal};
    $hash{bacteria} = $hash{mirco};
    $hash{virus} = $hash{mirco};
    $hash{name}{vertebrates} = $hash{name}{animal};
    $hash{name}{invertebrates} = $hash{name}{animal};
    $hash{name}{bacteria} = $hash{name}{mirco};
    $hash{name}{virus} = $hash{name}{mirco};



    return %hash;
}

sub Excu{
    my $cmd = shift;

    print "$cmd\n";
    if($opts{run} == 1){
        my $ret = system($cmd);
        Error("Error[$ret]:$cmd") if($ret != 0);
    }

    return;
}

sub Error{

    confess("Error:@_\n");
    return;
}
