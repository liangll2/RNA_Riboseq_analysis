# Step1: Filter raw data with fastp v0.18.0
fastp -a TGGAATTCT -q 20 -u 50 -w 4 -l 10 -i h143B-1r.fq.gz -o clean_fq_dir/h143B-1r.fq.gz -j clean_fq_dir/h143B-1r.json -h clean_fq_dir/h143B-1r.html
fastp -a TGGAATTCT -q 20 -u 50 -w 4 -l 10 -i h143B-2r.fq.gz -o clean_fq_dir/h143B-2r.fq.gz -j clean_fq_dir/h143B-2r.json -h clean_fq_dir/h143B-2r.html
fastp -a TGGAATTCT -q 20 -u 50 -w 4 -l 10 -i h143B-3r.fq.gz -o clean_fq_dir/h143B-3r.fq.gz -j clean_fq_dir/h143B-3r.json -h clean_fq_dir/h143B-3r.html
fastp -a TGGAATTCT -q 20 -u 50 -w 4 -l 10 -i hFOB-1r.fq.gz -o clean_fq_dir/hFOB-1r.fq.gz -j clean_fq_dir/hFOB-1r.json -h clean_fq_dir/hFOB-1r.html
fastp -a TGGAATTCT -q 20 -u 50 -w 4 -l 10 -i hFOB-2r.fq.gz -o clean_fq_dir/hFOB-2r.fq.gz -j clean_fq_dir/hFOB-2r.json -h clean_fq_dir/hFOB-2r.html
fastp -a TGGAATTCT -q 20 -u 50 -w 4 -l 10 -i hFOB-3r.fq.gz -o clean_fq_dir/hFOB-3r.fq.gz -j clean_fq_dir/hFOB-3r.json -h clean_fq_dir/hFOB-3r.html
fastp -a TGGAATTCT -q 20 -u 50 -w 4 -l 10 -i MG63r.fq.gz -o clean_fq_dir/MG63r.fq.gz -j clean_fq_dir/MG63r.json -h clean_fq_dir/MG63r.html
fastp -a TGGAATTCT -q 20 -u 50 -w 4 -l 10 -i Soas2r.fq.gz -o clean_fq_dir/Soas2r.fq.gz -j clean_fq_dir/Soas2r.json -h clean_fq_dir/Soas2r.html

# Step2: rRNA removal with bowtie2 v2.2.8
bowtie2 --local -x $rRNA_index -U clean_fq_dir/h143B-1r.fq.gz --un-gz rRNAremove/h143B-1r.fq.gz
bowtie2 --local -x $rRNA_index -U clean_fq_dir/h143B-2r.fq.gz --un-gz rRNAremove/h143B-2r.fq.gz
bowtie2 --local -x $rRNA_index -U clean_fq_dir/h143B-3r.fq.gz --un-gz rRNAremove/h143B-3r.fq.gz
bowtie2 --local -x $rRNA_index -U clean_fq_dir/hFOB-1r.fq.gz --un-gz rRNAremove/hFOB-1r.fq.gz
bowtie2 --local -x $rRNA_index -U clean_fq_dir/hFOB-2r.fq.gz --un-gz rRNAremove/hFOB-2r.fq.gz
bowtie2 --local -x $rRNA_index -U clean_fq_dir/hFOB-3r.fq.gz --un-gz rRNAremove/hFOB-3r.fq.gz
bowtie2 --local -x $rRNA_index -U clean_fq_dir/MG63r.fq.gz --un-gz rRNAremove/MG63r.fq.gz
bowtie2 --local -x $rRNA_index -U clean_fq_dir/Soas2r.fq.gz --un-gz rRNAremove/Soas2r.fq.gz

# Step3: tRNA sRNA removal with bowtie2 v2.2.8
bowtie2 --local -x $index -U rRNAremove/h143B-1r.fq.gz --un-gz sRNAremove/h143B-1r.fq.gz
bowtie2 --local -x $index -U rRNAremove/h143B-2r.fq.gz --un-gz sRNAremove/h143B-2r.fq.gz
bowtie2 --local -x $index -U rRNAremove/h143B-3r.fq.gz --un-gz sRNAremove/h143B-3r.fq.gz
bowtie2 --local -x $index -U rRNAremove/hFOB-1r.fq.gz --un-gz sRNAremove/hFOB-1r.fq.gz
bowtie2 --local -x $index -U rRNAremove/hFOB-2r.fq.gz --un-gz sRNAremove/hFOB-2r.fq.gz
bowtie2 --local -x $index -U rRNAremove/hFOB-3r.fq.gz --un-gz sRNAremove/hFOB-3r.fq.gz
bowtie2 --local -x $index -U rRNAremove/MG63r.fq.gz --un-gz sRNAremove/MG63r.fq.gz
bowtie2 --local -x $index -U rRNAremove/Soas2r.fq.gz --un-gz sRNAremove/Soas2r.fq.gz

# Step4: Extraction of reads with lengths between 20 and 40 bp
zcat sRNAremove/h143B-1r.fq.gz | awk 'NR%4==1{id=$0;flag=0}; (NR%4==2 && length($0) >= 20 && length($0) <= 40){flag=1;print id"\n"$0"\n+"}; (NR%4==0 && flag==1){print $0}'| gzip > sRNAremove/h143B-1r.fq.gz
zcat sRNAremove/h143B-2r.fq.gz | awk 'NR%4==1{id=$0;flag=0}; (NR%4==2 && length($0) >= 20 && length($0) <= 40){flag=1;print id"\n"$0"\n+"}; (NR%4==0 && flag==1){print $0}'| gzip > sRNAremove/h143B-2r.fq.gz
zcat sRNAremove/h143B-3r.fq.gz | awk 'NR%4==1{id=$0;flag=0}; (NR%4==2 && length($0) >= 20 && length($0) <= 40){flag=1;print id"\n"$0"\n+"}; (NR%4==0 && flag==1){print $0}'| gzip > sRNAremove/h143B-3r.fq.gz
zcat sRNAremove/hFOB-1r.fq.gz | awk 'NR%4==1{id=$0;flag=0}; (NR%4==2 && length($0) >= 20 && length($0) <= 40){flag=1;print id"\n"$0"\n+"}; (NR%4==0 && flag==1){print $0}'| gzip > sRNAremove/hFOB-1r.fq.gz
zcat sRNAremove/hFOB-2r.fq.gz | awk 'NR%4==1{id=$0;flag=0}; (NR%4==2 && length($0) >= 20 && length($0) <= 40){flag=1;print id"\n"$0"\n+"}; (NR%4==0 && flag==1){print $0}'| gzip > sRNAremove/hFOB-2r.fq.gz
zcat sRNAremove/hFOB-3r.fq.gz | awk 'NR%4==1{id=$0;flag=0}; (NR%4==2 && length($0) >= 20 && length($0) <= 40){flag=1;print id"\n"$0"\n+"}; (NR%4==0 && flag==1){print $0}'| gzip > sRNAremove/hFOB-3r.fq.gz
zcat sRNAremove/MG63r.fq.gz | awk 'NR%4==1{id=$0;flag=0}; (NR%4==2 && length($0) >= 20 && length($0) <= 40){flag=1;print id"\n"$0"\n+"}; (NR%4==0 && flag==1){print $0}'| gzip > sRNAremove/MG63r.fq.gz
zcat sRNAremove/Soas2r.fq.gz | awk 'NR%4==1{id=$0;flag=0}; (NR%4==2 && length($0) >= 20 && length($0) <= 40){flag=1;print id"\n"$0"\n+"}; (NR%4==0 && flag==1){print $0}'| gzip > sRNAremove/Soas2r.fq.gz

# Step5: Align filtered reads to genome using STAR v020201
STAR --runMode alignReads --readFilesCommand zcat --twopassMode Basic --outSAMtype BAM Unsorted --outSAMunmapped None --genomeDir $genomeDir --readFilesIn sRNAremove/h143B-1r.fq.gz --outFileNamePrefix align/h143B-1r
STAR --runMode alignReads --readFilesCommand zcat --twopassMode Basic --outSAMtype BAM Unsorted --outSAMunmapped None --genomeDir $genomeDir --readFilesIn sRNAremove/h143B-2r.fq.gz --outFileNamePrefix align/h143B-2r
STAR --runMode alignReads --readFilesCommand zcat --twopassMode Basic --outSAMtype BAM Unsorted --outSAMunmapped None --genomeDir $genomeDir --readFilesIn sRNAremove/h143B-3r.fq.gz --outFileNamePrefix align/h143B-3r
STAR --runMode alignReads --readFilesCommand zcat --twopassMode Basic --outSAMtype BAM Unsorted --outSAMunmapped None --genomeDir $genomeDir --readFilesIn sRNAremove/hFOB-1r.fq.gz --outFileNamePrefix align/hFOB-1r
STAR --runMode alignReads --readFilesCommand zcat --twopassMode Basic --outSAMtype BAM Unsorted --outSAMunmapped None --genomeDir $genomeDir --readFilesIn sRNAremove/hFOB-2r.fq.gz --outFileNamePrefix align/hFOB-2r
STAR --runMode alignReads --readFilesCommand zcat --twopassMode Basic --outSAMtype BAM Unsorted --outSAMunmapped None --genomeDir $genomeDir --readFilesIn sRNAremove/hFOB-3r.fq.gz --outFileNamePrefix align/hFOB-3r
STAR --runMode alignReads --readFilesCommand zcat --twopassMode Basic --outSAMtype BAM Unsorted --outSAMunmapped None --genomeDir $genomeDir --readFilesIn sRNAremove/MG63r.fq.gz --outFileNamePrefix align/MG63r
STAR --runMode alignReads --readFilesCommand zcat --twopassMode Basic --outSAMtype BAM Unsorted --outSAMunmapped None --genomeDir $genomeDir --readFilesIn sRNAremove/Soas2r.fq.gz --outFileNamePrefix align/Soas2r

# Step6:  Trinucleotide periodicity analysis using riboWaltz v1.0
bowtie2 -k 10 --no-unal -x $genome_exon_index -U sRNAremove/h143B-1r.fq.gz | samtools view -Sb -| samtools sort - align/exon/h143B-1r
bowtie2 -k 10 --no-unal -x $genome_exon_index -U sRNAremove/h143B-2r.fq.gz | samtools view -Sb -| samtools sort - align/exon/h143B-2r
bowtie2 -k 10 --no-unal -x $genome_exon_index -U sRNAremove/h143B-3r.fq.gz | samtools view -Sb -| samtools sort - align/exon/h143B-3r
bowtie2 -k 10 --no-unal -x $genome_exon_index -U sRNAremove/hFOB-1r.fq.gz | samtools view -Sb -| samtools sort - align/exon/hFOB-1r
bowtie2 -k 10 --no-unal -x $genome_exon_index -U sRNAremove/hFOB-2r.fq.gz | samtools view -Sb -| samtools sort - align/exon/hFOB-2r
bowtie2 -k 10 --no-unal -x $genome_exon_index -U sRNAremove/hFOB-3r.fq.gz | samtools view -Sb -| samtools sort - align/exon/hFOB-3r
bowtie2 -k 10 --no-unal -x $genome_exon_index -U sRNAremove/MG63r.fq.gz | samtools view -Sb -| samtools sort - align/exon/MG63r
bowtie2 -k 10 --no-unal -x $genome_exon_index -U sRNAremove/Soas2r.fq.gz | samtools view -Sb -| samtools sort - align/exon/Soas2r
perl draw_ud_v3.pl align/exon/h143B-1r.bam $genome.gtf
perl draw_ud_v3.pl align/exon/h143B-2r.bam $genome.gtf
perl draw_ud_v3.pl align/exon/h143B-3r.bam $genome.gtf
perl draw_ud_v3.pl align/exon/hFOB-1r.bam $genome.gtf
perl draw_ud_v3.pl align/exon/hFOB-2r.bam $genome.gtf
perl draw_ud_v3.pl align/exon/hFOB-3r.bam $genome.gtf
perl draw_ud_v3.pl align/exon/MG63r.bam $genome.gtf
perl draw_ud_v3.pl align/exon/Soas2r.bam $genome.gtf

# Step7: Gene expression at the translational level
# The align_and_estimate_abundance.pl script in Trinity utilizes RSEM v1.2.19 software to calculate expression.
perl align_and_estimate_abundance.pl --transcripts $genome_cds_imdex --gene_trans_map $gene_trans_map --est_method RSEM --aln_method bowtie2 --seqType fq --single sRNAremove/h143B-1r.fq.gz --output_dir EXPS/h143B-1r --output_prefix h143B-1r --quality --phred33-quals --rsem_add_opts --seed-length 10
perl align_and_estimate_abundance.pl --transcripts $genome_cds_imdex --gene_trans_map $gene_trans_map --est_method RSEM --aln_method bowtie2 --seqType fq --single sRNAremove/h143B-2r.fq.gz --output_dir EXPS/h143B-2r --output_prefix h143B-2r --quality --phred33-quals --rsem_add_opts --seed-length 10
perl align_and_estimate_abundance.pl --transcripts $genome_cds_imdex --gene_trans_map $gene_trans_map --est_method RSEM --aln_method bowtie2 --seqType fq --single sRNAremove/h143B-3r.fq.gz --output_dir EXPS/h143B-3r --output_prefix h143B-3r --quality --phred33-quals --rsem_add_opts --seed-length 10
perl align_and_estimate_abundance.pl --transcripts $genome_cds_imdex --gene_trans_map $gene_trans_map --est_method RSEM --aln_method bowtie2 --seqType fq --single sRNAremove/hFOB-1r.fq.gz --output_dir EXPS/hFOB-1r --output_prefix hFOB-1r --quality --phred33-quals --rsem_add_opts --seed-length 10
perl align_and_estimate_abundance.pl --transcripts $genome_cds_imdex --gene_trans_map $gene_trans_map --est_method RSEM --aln_method bowtie2 --seqType fq --single sRNAremove/hFOB-2r.fq.gz --output_dir EXPS/hFOB-2r --output_prefix hFOB-2r --quality --phred33-quals --rsem_add_opts --seed-length 10
perl align_and_estimate_abundance.pl --transcripts $genome_cds_imdex --gene_trans_map $gene_trans_map --est_method RSEM --aln_method bowtie2 --seqType fq --single sRNAremove/hFOB-3r.fq.gz --output_dir EXPS/hFOB-3r --output_prefix hFOB-3r --quality --phred33-quals --rsem_add_opts --seed-length 10
perl align_and_estimate_abundance.pl --transcripts $genome_cds_imdex --gene_trans_map $gene_trans_map --est_method RSEM --aln_method bowtie2 --seqType fq --single sRNAremove/MG63r.fq.gz --output_dir EXPS/MG63r --output_prefix MG63r --quality --phred33-quals --rsem_add_opts --seed-length 10
perl align_and_estimate_abundance.pl --transcripts $genome_cds_imdex --gene_trans_map $gene_trans_map --est_method RSEM --aln_method bowtie2 --seqType fq --single sRNAremove/Soas2r.fq.gz --output_dir EXPS/Soas2r --output_prefix Soas2r --quality --phred33-quals --rsem_add_opts --seed-length 10

# Step8: Differential expression analysis use edgeR v3.12.1
perl run_DE_analysis.pl --matrix all.genes.counts.matrix --method edgeR --output genes_edgeR_result --samples_file Groups.txt --contrasts Contrast.txt

# Step9: Translational Efficiency Calculation
perl TE_analysis.pl RNAseq_all.genes.expression.xls Riboseq_all.genes.expression.xls TE_sample.list TE_calculate/TE_sample
python get_group_TE.py sample_group.list TE_calculate/TE_sample.xls TE_calculate/TE_group
