# Usage: bash rna-seq-map-Bd.sh BdHS01-1TP5_S1
# 
# Purpose: map Brachy RNA-seq, pair-ended.
#
# Note: For example, this script reads BdHS01-1TP5_S1_R1_raw.fastq and BdHS01-1TP5_S1_R2_raw.fastq,
#       and map these reads to Brachy transcriptome.
#
# 7 OCT 2016 hui

filename=$1

# edit the following upper-case variables if necessary
GFF3_FILE=/media/pw_synology3/PW_HiSeq_data/Reference_genomes/Bdistachyon_v3.1_Phytozome11_downloaded03062016/v3.1/annotation/Bdistachyon_314_v3.1.gene_exons.gff3
GTF_FILE=/home/hui/rna-seq-map/Bdistachyon_314_v3.1.gene_exons.gtf  # produced by gffread with GFF3_FILE as input, prepared for htseq-count
TRANSCRIPTOME_INDEX=/home/hui/rna-seq-map/transcriptome/gene_index # will be produced by tophat in the first run
BOWTIE2INDEX_PATH=/media/pw_synology3/PW_HiSeq_data/Reference_genomes/Bdistachyon_v3.1_Phytozome11_downloaded03062016/v3.1/Bowtie2Index # must be prepared before running this script
GENOME_LENGTH=271923306 # number obtained from http://plants.ensembl.org/Brachypodium_distachyon/Info/Annotation/
CHROMOSOME_INFO_FILE=/media/pw_synology3/PW_HiSeq_data/Reference_genomes/Bdistachyon_v3.1_Phytozome11_downloaded03062016/v3.1/ChromInfo.txt # manually created
LIBRARY_TYPE=fr-firststrand  # options for tophat and cufflinks


# remove adaptors
echo "...removing adaptors" > run_"$filename".log
java -jar /home/Program_NGS_sl-pw-srv01/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads 6 \
-trimlog trimmolog1.txt \
"$filename"_R1_raw.fastq \
"$filename"_R2_raw.fastq \
"$filename"_R1_raw_trimmo_paired_truseq3-2_2_10_5_1.fastq \
"$filename"_R1_raw_trimmo_unpaired_truseq3-2_2_10_5_1.fastq \
"$filename"_R2_raw_trimmo_paired_truseq3-2_2_10_5_1.fastq \
"$filename"_R2_raw_trimmo_unpaired_truseq3-2_2_10_5_1.fastq \
ILLUMINACLIP:/home/Program_NGS_sl-pw-srv01/Trimmomatic-0.32/adapters/TruSeq3-PE-2.fa:2:10:5:1 2> trimmolog2.txt


# map with tophat
mkdir -p tophat_results
echo "...mapping with tophat" >> run_"$filename".log
tophat \
--max-multihits 1 \
--num-threads 6 \
--library-type "$LIBRARY_TYPE" \
--GTF "$GFF3_FILE" \
--transcriptome-index "$TRANSCRIPTOME_INDEX" \
--output-dir tophat_results/"$filename"_trimmo_nomixed/ \
--no-mixed \
"$BOWTIE2INDEX_PATH"/bd \
"$filename"_R1_raw_trimmo_paired_truseq3-2_2_10_5_1.fastq "$filename"_R2_raw_trimmo_paired_truseq3-2_2_10_5_1.fastq 2> tophatlog1.txt

# rename and sort
mv tophat_results/"$filename"_trimmo_nomixed/accepted_hits.bam \
   tophat_results/"$filename"_trimmo_nomixed/"$filename"_trimmo_paired_2_10_5_1_tophat_Bd_nomixed.bam

echo "...sorting mapped reads" >> run_"$filename".log
samtools sort tophat_results/"$filename"_trimmo_nomixed/"$filename"_trimmo_paired_2_10_5_1_tophat_Bd_nomixed.bam \
              tophat_results/"$filename"_trimmo_nomixed/"$filename"_trimmo_paired_2_10_5_1_tophat_Bd_nomixed_sorted 2>> run_"$filename".log


# mark and remove duplicates 
echo "...removing duplicates" >> run_"$filename".log
java -Xmx4g -jar /home/Program_NGS_sl-pw-srv01/picard-tools-1.103/MarkDuplicates.jar \
INPUT=tophat_results/"$filename"_trimmo_nomixed/"$filename"_trimmo_paired_2_10_5_1_tophat_Bd_nomixed_sorted.bam \
OUTPUT=tophat_results/"$filename"_trimmo_nomixed/"$filename"_trimmo_paired_2_10_5_1_tophat_Bd_nomixed_sorted_rmdup_picard.bam \
METRICS_FILE=dup.txt VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true 2> markdup_stderr.txt


# indexing
samtools index tophat_results/"$filename"_trimmo_nomixed/"$filename"_trimmo_paired_2_10_5_1_tophat_Bd_nomixed_sorted_rmdup_picard.bam

# count number of reads in raw data
raw_line=$(wc -l "$filename"_R1_raw.fastq "$filename"_R2_raw.fastq | grep 'total' | grep -o '[0-9]*') 
raw_count=$(echo "$raw_line/4" | bc -l)
echo "number of raw reads: $raw_count" >> run_"$filename".log

# get flagstat
echo "number of clean reads mapped:" >> run_"$filename".log
samtools flagstat tophat_results/"$filename"_trimmo_nomixed/"$filename"_trimmo_paired_2_10_5_1_tophat_Bd_nomixed_sorted_rmdup_picard.bam \
>> run_"$filename".log

# estimate genome average and normalise
echo "...normalising reads" >> run_"$filename".log
genomeCoverageBed -split -bg -ibam tophat_results/"$filename"_trimmo_nomixed/"$filename"_trimmo_paired_2_10_5_1_tophat_Bd_nomixed_sorted_rmdup_picard.bam \
-g $CHROMOSOME_INFO_FILE \
> tophat_results/"$filename"_trimmo_nomixed/"$filename"_trimmo_paired_2_10_5_1_tophat_Bd_nomixed_sorted_rmdup_picard.bedgraph


sum=$(samtools depth tophat_results/"$filename"_trimmo_nomixed/"$filename"_trimmo_paired_2_10_5_1_tophat_Bd_nomixed_sorted_rmdup_picard.bam | awk '{sum+=$3;cnt++}END{printf "%.0f", sum}')
sum_norm=$(echo "$sum/$GENOME_LENGTH" | bc -l)
echo "genome normalised coverage: $sum_norm" >> run_"$filename".log

export MYVAR=$sum_norm
perl -e 'print $ENV{MYVAR}."\n"'

# normalise read counts by genome-wide coverage
perl -ne 'chomp($_); @a=split(/\t/,$_);print $a[0]."\t".$a[1]."\t".$a[2]."\t".$a[3]/$ENV{MYVAR}."\n";' \
tophat_results/"$filename"_trimmo_nomixed/"$filename"_trimmo_paired_2_10_5_1_tophat_Bd_nomixed_sorted_rmdup_picard.bedgraph \
> tophat_results/"$filename"_trimmo_nomixed/"$filename"_trimmo_paired_2_10_5_1_tophat_Bd_nomixed_sorted_rmdup_picard_genomenorm.bedgraph

# convert bedgraph to bigwig
bedGraphToBigWig tophat_results/"$filename"_trimmo_nomixed/"$filename"_trimmo_paired_2_10_5_1_tophat_Bd_nomixed_sorted_rmdup_picard_genomenorm.bedgraph \
$CHROMOSOME_INFO_FILE \
tophat_results/"$filename"_trimmo_nomixed/"$filename"_trimmo_paired_2_10_5_1_tophat_Bd_nomixed_sorted_rmdup_picard_genomenorm.bw

# cufflinks
echo "...calling FPKMs" >> run_"$filename".log
cufflinks -q --output-dir cufflinks_results/cufflinks_"$filename"_trimmo_nomixed/ \
--GTF "$GFF3_FILE" \
--num-threads 10 \
--frag-bias-correct "$BOWTIE2INDEX_PATH"/bd.fa \
--multi-read-correct \
--library-type "$LIBRARY_TYPE" \
tophat_results/"$filename"_trimmo_nomixed/"$filename"_trimmo_paired_2_10_5_1_tophat_Bd_nomixed_sorted_rmdup_picard.bam

# obtain raw reads using HTseq-count
echo "...calling raw reads" >> run_"$filename".log
samtools sort -n tophat_results/"$filename"_trimmo_nomixed/"$filename"_trimmo_paired_2_10_5_1_tophat_Bd_nomixed_sorted_rmdup_picard.bam \
tophat_results/"$filename"_trimmo_nomixed/"$filename"_trimmo_paired_2_10_5_1_tophat_Bd_nomixed_sorted_rmdup_picard_name_sorted 2>> run_"$filename".log

htseq-count -q -r name -s no -f bam -t exon -i gene_id \
tophat_results/"$filename"_trimmo_nomixed/"$filename"_trimmo_paired_2_10_5_1_tophat_Bd_nomixed_sorted_rmdup_picard_name_sorted.bam \
"$GTF_FILE" \
> tophat_results/"$filename"_trimmo_nomixed/"$filename"_trimmo_paired_2_10_5_1_tophat_Bd_nomixed_sorted_rmdup_picard_name_sorted_htseq_count.ct

# combine FPKM/TPM/RAW
echo "...combining FPKM/TPM/Raw" >> run_"$filename".log
cp cufflinks_results/cufflinks_"$filename"_trimmo_nomixed/genes.fpkm_tracking .
cut -f1,4,7,10 genes.fpkm_tracking > genes.fpkm_tracking_fpkm
cp tophat_results/"$filename"_trimmo_nomixed/"$filename"_trimmo_paired_2_10_5_1_tophat_Bd_nomixed_sorted_rmdup_picard_name_sorted_htseq_count.ct genes.raw

Rscript -e 'genes.fpkm <- read.table("genes.fpkm_tracking_fpkm", sep="\t", quote="", header=T, colClasses = c(rep("character", 3), "numeric"))' \
-e 'genes.fpkm$TPM <-genes.fpkm$FPKM/sum(genes.fpkm$FPKM)*10^6' \
-e 'genes.raw <- read.table("genes.raw", sep="\t", quote="", header=F, colClasses = c("character", "numeric"))' \
-e 'names(genes.raw) <- c("tracking_id", "Raw")' \
-e 'genes.fpkm.raw <- merge(genes.fpkm, genes.raw, by = "tracking_id")' \
-e 'write.table(genes.fpkm.raw, "genes_fpkm_tpm_raw.txt", sep="\t", quote=F, row.names=F)' 

mv genes_fpkm_tpm_raw.txt "$filename"_trimmo_paired_2_10_5_1_tophat_Bd_nomixed_sorted_rmdup_picard_combined_read.txt # final gene expression table containing FPKM, TPM and raw counts

# remove files
rm -f trimmolog1.txt
rm -f *.bedgraph
rm -f genes.fpkm_tracking
rm -f genes.fpkm_tracking_fpkm
rm -f genes.raw
