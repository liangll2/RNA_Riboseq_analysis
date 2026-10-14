# Step1: Filter raw data with fastp v0.18.0
fastp -a AGATCGGAAGAGC -q 20 -u 50 -n 15 -l 50 -i MG63_1.fq.gz -I MG63_2.fq.gz -o clean_fq_dir/MG63_1.fq.gz -O clean_fq_dir/MG63_2.fq.gz -j clean_fq_dir/MG63.json -h clean_fq_dir/MG63.html
fastp -a AGATCGGAAGAGC -q 20 -u 50 -n 15 -l 50 -i Soas2_1.fq.gz -I Soas2_2.fq.gz -o clean_fq_dir/Soas2_1.fq.gz -O clean_fq_dir/Soas2_2.fq.gz -j clean_fq_dir/Soas2.json -h clean_fq_dir/Soas2.html
fastp -a AGATCGGAAGAGC -q 20 -u 50 -n 15 -l 50 -i h143B-1_1.fq.gz -I h143B-1_2.fq.gz -o clean_fq_dir/h143B-1_1.fq.gz -O clean_fq_dir/h143B-1_2.fq.gz -j clean_fq_dir/h143B-1.json -h clean_fq_dir/h143B-1.html
fastp -a AGATCGGAAGAGC -q 20 -u 50 -n 15 -l 50 -i h143B-2_1.fq.gz -I h143B-2_2.fq.gz -o clean_fq_dir/h143B-2_1.fq.gz -O clean_fq_dir/h143B-2_2.fq.gz -j clean_fq_dir/h143B-2.json -h clean_fq_dir/h143B-2.html
fastp -a AGATCGGAAGAGC -q 20 -u 50 -n 15 -l 50 -i h143B-3_1.fq.gz -I h143B-3_2.fq.gz -o clean_fq_dir/h143B-3_1.fq.gz -O clean_fq_dir/h143B-3_2.fq.gz -j clean_fq_dir/h143B-3.json -h clean_fq_dir/h143B-3.html
fastp -a AGATCGGAAGAGC -q 20 -u 50 -n 15 -l 50 -i hFOB-1_1.fq.gz -I hFOB-1_2.fq.gz -o clean_fq_dir/hFOB-1_1.fq.gz -O clean_fq_dir/hFOB-1_2.fq.gz -j clean_fq_dir/hFOB-1.json -h clean_fq_dir/hFOB-1.html
fastp -a AGATCGGAAGAGC -q 20 -u 50 -n 15 -l 50 -i hFOB-2_1.fq.gz -I hFOB-2_2.fq.gz -o clean_fq_dir/hFOB-2_1.fq.gz -O clean_fq_dir/hFOB-2_2.fq.gz -j clean_fq_dir/hFOB-2.json -h clean_fq_dir/hFOB-2.html
fastp -a AGATCGGAAGAGC -q 20 -u 50 -n 15 -l 50 -i hFOB-3_1.fq.gz -I hFOB-3_2.fq.gz -o clean_fq_dir/hFOB-3_1.fq.gz -O clean_fq_dir/hFOB-3_2.fq.gz -j clean_fq_dir/hFOB-3.json -h clean_fq_dir/hFOB-3.html

# Step2: Align clean reads to rRNA with bowtie2 v2.4.4
bowtie2 --local -x $rRNA_index -1 clean_fq_dir/MG63_1.fq.gz -2 clean_fq_dir/MG63_2.fq.gz --un-conc-gz rRNAremove/MG63_%.fq.gz
bowtie2 --local -x $rRNA_index -1 clean_fq_dir/Soas2_1.fq.gz -2 clean_fq_dir/Soas2_2.fq.gz --un-conc-gz rRNAremove/Soas2_%.fq.gz
bowtie2 --local -x $rRNA_index -1 clean_fq_dir/h143B-1_1.fq.gz -2 clean_fq_dir/h143B-1_2.fq.gz --un-conc-gz rRNAremove/h143B-1_%.fq.gz
bowtie2 --local -x $rRNA_index -1 clean_fq_dir/h143B-2_1.fq.gz -2 clean_fq_dir/h143B-2_2.fq.gz --un-conc-gz rRNAremove/h143B-2_%.fq.gz
bowtie2 --local -x $rRNA_index -1 clean_fq_dir/h143B-3_1.fq.gz -2 clean_fq_dir/h143B-3_2.fq.gz --un-conc-gz rRNAremove/h143B-3_%.fq.gz
bowtie2 --local -x $rRNA_index -1 clean_fq_dir/hFOB-1_1.fq.gz -2 clean_fq_dir/hFOB-1_2.fq.gz --un-conc-gz rRNAremove/hFOB-1_%.fq.gz
bowtie2 --local -x $rRNA_index -1 clean_fq_dir/hFOB-2_1.fq.gz -2 clean_fq_dir/hFOB-2_2.fq.gz --un-conc-gz rRNAremove/hFOB-2_%.fq.gz
bowtie2 --local -x $rRNA_index -1 clean_fq_dir/hFOB-3_1.fq.gz -2 clean_fq_dir/hFOB-3_2.fq.gz --un-conc-gz rRNAremove/hFOB-3_%.fq.gz

# Step3: Align clean reads to genome using HISAT2 v2.1.0
hisat2 --rna-strandness RF -x $genome_index -q -1 rRNAremove/MG63_1.fq.gz -2 rRNAremove/MG63_2.fq.gz --known-splicesite-infile $genome.splices_sites --dta --summary-file align_genome/MG63.stat --new-summary >align_genome/MG63.sam
hisat2 --rna-strandness RF -x $genome_index -q -1 rRNAremove/Soas2_1.fq.gz -2 rRNAremove/Soas2_2.fq.gz --known-splicesite-infile $genome.splices_sites --dta --summary-file align_genome/Soas2.stat --new-summary >align_genome/Soas2.sam
hisat2 --rna-strandness RF -x $genome_index -q -1 rRNAremove/h143B-1_1.fq.gz -2 rRNAremove/h143B-1_2.fq.gz --known-splicesite-infile $genome.splices_sites --dta --summary-file align_genome/h143B-1.stat --new-summary >align_genome/h143B-1.sam
hisat2 --rna-strandness RF -x $genome_index -q -1 rRNAremove/h143B-2_1.fq.gz -2 rRNAremove/h143B-2_2.fq.gz --known-splicesite-infile $genome.splices_sites --dta --summary-file align_genome/h143B-2.stat --new-summary >align_genome/h143B-2.sam
hisat2 --rna-strandness RF -x $genome_index -q -1 rRNAremove/h143B-3_1.fq.gz -2 rRNAremove/h143B-3_2.fq.gz --known-splicesite-infile $genome.splices_sites --dta --summary-file align_genome/h143B-3.stat --new-summary >align_genome/h143B-3.sam
hisat2 --rna-strandness RF -x $genome_index -q -1 rRNAremove/hFOB-1_1.fq.gz -2 rRNAremove/hFOB-1_2.fq.gz --known-splicesite-infile $genome.splices_sites --dta --summary-file align_genome/hFOB-1.stat --new-summary >align_genome/hFOB-1.sam
hisat2 --rna-strandness RF -x $genome_index -q -1 rRNAremove/hFOB-2_1.fq.gz -2 rRNAremove/hFOB-2_2.fq.gz --known-splicesite-infile $genome.splices_sites --dta --summary-file align_genome/hFOB-2.stat --new-summary >align_genome/hFOB-2.sam
hisat2 --rna-strandness RF -x $genome_index -q -1 rRNAremove/hFOB-3_1.fq.gz -2 rRNAremove/hFOB-3_2.fq.gz --known-splicesite-infile $genome.splices_sites --dta --summary-file align_genome/hFOB-3.stat --new-summary >align_genome/hFOB-3.sam

# Step4: Transcriptome reconstruction and annotation using stringtie v1.3.4 
stringtie --rf -f 0.3 align_genome/MG63.sort.bam -G $genome.gtf -o assembly/MG63.gtf
stringtie --rf -f 0.3 align_genome/Soas2.sort.bam -G $genome.gtf -o assembly/Soas2.gtf
stringtie --rf -f 0.3 align_genome/h143B-1.sort.bam -G $genome.gtf -o assembly/h143B-1.gtf
stringtie --rf -f 0.3 align_genome/h143B-2.sort.bam -G $genome.gtf -o assembly/h143B-2.gtf
stringtie --rf -f 0.3 align_genome/h143B-3.sort.bam -G $genome.gtf -o assembly/h143B-3.gtf
stringtie --rf -f 0.3 align_genome/hFOB-1.sort.bam -G $genome.gtf -o assembly/hFOB-1.gtf
stringtie --rf -f 0.3 align_genome/hFOB-2.sort.bam -G $genome.gtf -o assembly/hFOB-2.gtf
stringtie --rf -f 0.3 align_genome/hFOB-3.sort.bam -G $genome.gtf -o assembly/hFOB-3.gtf
ls $PWD/assembly/*.gtf  > assembly/assembly_GTF.list
stringtie --merge -f 0.3 -m 200 -G $genome.gtf -o assembly/merge.gtf assembly/assembly_GTF.list
cuffcompare -r $genome.gtf -o assembly/merged assembly/merged.gtf
# Based on the cuffcompare results, the novel transcripts are those with class_code u, i, j, x, e and having more than one exon.
cat $genome.gtf assembly/newgene.tmp.gtf > assembly/exp_ref.gtf

# Step5: Expression analysis
# The align_and_estimate_abundance.pl script in Trinity v2.0.6 utilizes RSEM v1.2.19 software to calculate expression.
gffread assembly/exp_ref.gtf -g $genome.fa -w RSEM/all.exon.fa
grep '^>' RSEM/all.exon.fa| awk '{sub(/>/,"",$1);sub(/gene=/,"",$2;print $2"\t"$1)}' > RSEM/gene_trans_map
perl align_and_estimate_abundance.pl --transcripts RSEM/all.exon.fa --gene_trans_map RSEM/gene_trans_map --est_method RSEM --aln_method bowtie2 --seqType fq --left Rremove/MG63_1.fq.gz --right Rremove/MG63_2.fq.gz --output_dir RSEM/MG63 --output_prefix MG63 --quality --phred33-quals --SS_lib_type RF
perl align_and_estimate_abundance.pl --transcripts RSEM/all.exon.fa --gene_trans_map RSEM/gene_trans_map --est_method RSEM --aln_method bowtie2 --seqType fq --left Rremove/Soas2_1.fq.gz --right Rremove/Soas2_2.fq.gz --output_dir RSEM/Soas2 --output_prefix Soas2 --quality --phred33-quals --SS_lib_type RF
perl align_and_estimate_abundance.pl --transcripts RSEM/all.exon.fa --gene_trans_map RSEM/gene_trans_map --est_method RSEM --aln_method bowtie2 --seqType fq --left Rremove/h143B-1_1.fq.gz --right Rremove/h143B-1_2.fq.gz --output_dir RSEM/h143B-1 --output_prefix h143B-1 --quality --phred33-quals --SS_lib_type RF
perl align_and_estimate_abundance.pl --transcripts RSEM/all.exon.fa --gene_trans_map RSEM/gene_trans_map --est_method RSEM --aln_method bowtie2 --seqType fq --left Rremove/h143B-2_1.fq.gz --right Rremove/h143B-2_2.fq.gz --output_dir RSEM/h143B-2 --output_prefix h143B-2 --quality --phred33-quals --SS_lib_type RF
perl align_and_estimate_abundance.pl --transcripts RSEM/all.exon.fa --gene_trans_map RSEM/gene_trans_map --est_method RSEM --aln_method bowtie2 --seqType fq --left Rremove/h143B-3_1.fq.gz --right Rremove/h143B-3_2.fq.gz --output_dir RSEM/h143B-3 --output_prefix h143B-3 --quality --phred33-quals --SS_lib_type RF
perl align_and_estimate_abundance.pl --transcripts RSEM/all.exon.fa --gene_trans_map RSEM/gene_trans_map --est_method RSEM --aln_method bowtie2 --seqType fq --left Rremove/hFOB-1_1.fq.gz --right Rremove/hFOB-1_2.fq.gz --output_dir RSEM/hFOB-1 --output_prefix hFOB-1 --quality --phred33-quals --SS_lib_type RF
perl align_and_estimate_abundance.pl --transcripts RSEM/all.exon.fa --gene_trans_map RSEM/gene_trans_map --est_method RSEM --aln_method bowtie2 --seqType fq --left Rremove/hFOB-2_1.fq.gz --right Rremove/hFOB-2_2.fq.gz --output_dir RSEM/hFOB-2 --output_prefix hFOB-2 --quality --phred33-quals --SS_lib_type RF
perl align_and_estimate_abundance.pl --transcripts RSEM/all.exon.fa --gene_trans_map RSEM/gene_trans_map --est_method RSEM --aln_method bowtie2 --seqType fq --left Rremove/hFOB-3_1.fq.gz --right Rremove/hFOB-3_2.fq.gz --output_dir RSEM/hFOB-3 --output_prefix hFOB-3 --quality --phred33-quals --SS_lib_type RF

# Step6: Differential expression analysis use DESeq2 v1.24; enrich analysis use kegg Release 94.0 and GO 3.8.2
perl run_RNAseq_diff.pl -f genes.expression.xls -g group.list -c group_diff.list -d 0.01 -p q -t 0.05  -l 2 -m yes -r yes -a 8.Annotation/1.annotation/all -e yes -run_go 1 -run_ko 1 -run_do 1 -run_reactome 1 -soft DESeq2 -Rscript Rscript.3.6.0 -o group_diff

# Step7: string analysis use STRING v10 and Cytoscape v3.7.1
perl run_string.pl string -infile hFOB-vs-h143B.filter.xls -db all.string.db -idmap all.string.idmap -symbol all.gen2sym -outpfx hFOB-vs-h143B
perl run_string.pl string -infile hFOB-vs-MG63.filter.xls -db all.string.db -idmap all.string.idmap -symbol all.gen2sym -outpfx hFOB-vs-MG63
perl run_string.pl string -infile hFOB-vs-Soas2.filter.xls -db all.string.db -idmap all.string.idmap -symbol all.gen2sym -outpfx hFOB-vs-Soas2

# Step8: trend analysis use STEM v1.0
perl run_trend.pl -info_conf Trend.list -dirinfo pipe_dir.list -pipe_type lncRNA -ratio 2 -run_do 1 -run_reactome 1 -outdir trend
