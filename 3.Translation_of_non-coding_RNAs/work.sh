# Step1: sORF prediction
perl orf_predict.pl -target_fa $longest_5utr.fa -pre sORF/prepare/uORF -minlength 63 -maxlength 453
perl orf_predict.pl -target_fa $longest_3utr.fa -pre sORF/prepare/dORF -minlength 63 -maxlength 453
cat prepare/uORF.bed prepare/dORF.bed prepare/longest_mORF.bed > sORF/prepare/all_ORF.bed

# Step2: Quantification of sORF translation abundance
perl filter_muti_align.pl EXPS/gene_trans.map align/exon/h143B-1r.bam sORF/bams/h143B-1r
perl filter_muti_align.pl EXPS/gene_trans.map align/exon/h143B-2r.bam sORF/bams/h143B-2r
perl filter_muti_align.pl EXPS/gene_trans.map align/exon/h143B-3r.bam sORF/bams/h143B-3r
perl filter_muti_align.pl EXPS/gene_trans.map align/exon/hFOB-1r.bam sORF/bams/hFOB-1r
perl filter_muti_align.pl EXPS/gene_trans.map align/exon/hFOB-2r.bam sORF/bams/hFOB-2r
perl filter_muti_align.pl EXPS/gene_trans.map align/exon/hFOB-3r.bam sORF/bams/hFOB-3r
perl filter_muti_align.pl EXPS/gene_trans.map align/exon/MG63r.bam sORF/bams/MG63r
perl filter_muti_align.pl EXPS/gene_trans.map align/exon/Soas2r.bam sORF/bams/Soas2r
perl calculate_ORF_rpkm.pl sORF/prepare/all_ORF.bed  sORF/bams/h143B-1r.filter.bam  sORF/bams/h143B-1r.stat.xls sORF/bams/h143B-1r
perl calculate_ORF_rpkm.pl sORF/prepare/all_ORF.bed  sORF/bams/h143B-1r.filter.bam  sORF/bams/h143B-1r.stat.xls sORF/bams/h143B-1r
perl calculate_ORF_rpkm.pl sORF/prepare/all_ORF.bed  sORF/bams/h143B-1r.filter.bam  sORF/bams/h143B-1r.stat.xls sORF/bams/h143B-1r
perl calculate_ORF_rpkm.pl sORF/prepare/all_ORF.bed  sORF/bams/hFOB-1r.filter.bam  sORF/bams/hFOB-1r.stat.xls sORF/bams/hFOB-1r
perl calculate_ORF_rpkm.pl sORF/prepare/all_ORF.bed  sORF/bams/hFOB-1r.filter.bam  sORF/bams/hFOB-1r.stat.xls sORF/bams/hFOB-1r
perl calculate_ORF_rpkm.pl sORF/prepare/all_ORF.bed  sORF/bams/hFOB-1r.filter.bam  sORF/bams/hFOB-1r.stat.xls sORF/bams/hFOB-1r
perl calculate_ORF_rpkm.pl sORF/prepare/all_ORF.bed  sORF/bams/MG63r.filter.bam  sORF/bams/MG63r.stat.xls sORF/bams/MG63r
perl calculate_ORF_rpkm.pl sORF/prepare/all_ORF.bed  sORF/bams/Soas2r.filter.bam  sORF/bams/Soas2r.stat.xls sORF/bams/Soas2r

# Step3: Differential analysis use edgeR v3.12.1
perl run_DE_analysis.pl --matrix all.ORFs.counts.filter.matrix --method edgeR --output orfs_edgeR_result --samples_file Groups.txt --contrasts Contrast.txt

# Step4: Prediction of coding potential
# Calculation of ORFscore
python correct_count_frame.py prepare/longest_mORF.bed $exp_filter.mORF_sORF.bed sORF/bams/h143B-1r.filter.bam sORF/RRS_ORFscore/h143B-1r.orf_frame.xls
python correct_count_frame.py prepare/longest_mORF.bed $exp_filter.mORF_sORF.bed sORF/bams/h143B-2r.filter.bam sORF/RRS_ORFscore/h143B-2r.orf_frame.xls
python correct_count_frame.py prepare/longest_mORF.bed $exp_filter.mORF_sORF.bed sORF/bams/h143B-3r.filter.bam sORF/RRS_ORFscore/h143B-3r.orf_frame.xls
python correct_count_frame.py prepare/longest_mORF.bed $exp_filter.mORF_sORF.bed sORF/bams/hFOB-1r.filter.bam sORF/RRS_ORFscore/hFOB-1r.orf_frame.xls
python correct_count_frame.py prepare/longest_mORF.bed $exp_filter.mORF_sORF.bed sORF/bams/hFOB-2r.filter.bam sORF/RRS_ORFscore/hFOB-2r.orf_frame.xls
python correct_count_frame.py prepare/longest_mORF.bed $exp_filter.mORF_sORF.bed sORF/bams/hFOB-3r.filter.bam sORF/RRS_ORFscore/hFOB-3r.orf_frame.xls
python correct_count_frame.py prepare/longest_mORF.bed $exp_filter.mORF_sORF.bed sORF/bams/MG63r.filter.bam sORF/RRS_ORFscore/MG63r.orf_frame.xls
python correct_count_frame.py prepare/longest_mORF.bed $exp_filter.mORF_sORF.bed sORF/bams/Soas2r.filter.bam sORF/RRS_ORFscore/Soas2r.orf_frame.xls
python combine_frame_ORFscore.py Groups.txt sORF/RRS_ORFscore sORF/RRS_ORFscore
# Calculation of RSS
python combine_rpkm_RRS.py Groups.txt sORF/ sORF/RRS_ORFscore
# Calculating Fickett and Hexamer Scores using CPAT 1.2.2
python make_hexamer_tab.py -c $genome_cds.fa -n $genome_lnc.fa > sORF/Fickett_Hexamer//ORF.Hexamer.tsv
python make_logitModel.py -x sORF/Fickett_Hexamer//ORF.Hexamer.tsv -c $_refGeneMrna.fa  -n $genome_lnc.fa -o sORF/Fickett_Hexamer//ORF
python cpat.py -g sORF/all_ORF.DNA.fasta -d sORF/Fickett_Hexamer//ORF.logit.RData -x sORF/Fickett_Hexamer//ORF.Hexamer.tsv -o ORF.result.xls
