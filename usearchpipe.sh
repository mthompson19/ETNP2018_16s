#!/bin/bash
#SBATCH --job-name=usearch16s
#SBATCH -n 2
#SBATCH -N 1
#SBATCH --output usearch16s.out
#SBATCH --error usearch16s.err
#SBATCH -p penglab-48core

#Active Qiime
source activate qiime2-2021.8

# This is a script to run 16S analysis using USEARCH

# Navigate to working directory
cd /work/mat19/usearchtest 

# Set directory to TMPDIR 

export TMPDIR='/usearchtest' 
echo $TMPDIR

#Naviage to raw_reads data

#Quality 

echo "Checking the quality of raw reads..."

mkdir ./fastq_eestats

for fq in *.fastq; do 

usearch -fastq_eestats2 $fq -output ./fastq_eestats/$fq -length_cutoffs 50,300,10
done

#Trim forward and reverse reads to 200bp

echo "Trimming reads to 200bp..."
for fq in *.fastq; do

NAME=$(basename $fq .fastq)
usearch -fastx_truncate ${fq} -trunclen 200 -fastqout ${NAME}_200bp.fastq
done

#Merge reads - 
echo "Merge reads"

usearch -fastq_mergepairs *R1.fastq -relabel @ -fastq_maxdiffs 10 -fastq_pctid 80 -fastqout ./200bp_merged.fq -report ./merge_200bp_report.txt

#Filter reads
echo "filter reads"
usearch -fastq_filter 200bp_merged.fq -fastq_maxee 1.0 -fastq_minlen 150 -fastaout filtered.fa

#Dereplicate reads 
echo "duplicate reads"
usearch -fastx_uniques filtered.fa -fastaout uniques.fa -sizeout -relabel Uniq

#Denoise 
echo "denoise" 
usearch -unoise3 uniques.fa -zotus zotus.fa
# or run this
usearch -cluster_otus uniques.fa -otus otus.fa 

#ZOTU table 
echo "ZOTU table" 
usearch -otutab 200bp_merged.fq -zotus otus.fa -otutabout otutab.txt -biomout otutab.json -mapout zmap.txt -notmatched unmapped.fa -dbmatched otus_with_sizes.fa -sizeout -threads 10
## go into zotus.fa and remove the methanobacterium and bacillus (you have to merge the ZOTU and Taxonomy in phyloseq to figure out which ZOTU is each 
### they are contaminants and only found in abundance in Stn3Cast9370m02um 

#Import ZOTU to Qiime
echo "Import ZOTU to Qiime"
qiime tools import --input-path otus.fa --output-path otus.qza --type 'FeatureData[Sequence]'


#Qiime Taxonomy 
echo "Qiime Taxonomy"
qiime feature-classifier classify-sklearn \
 --i-classifier silva-138-99-515-806-nb-classifier.qza \
 --i-reads otus.qza \
 --o-classification otu_16S_taxonomy.qza

#Export to tsv
echo "export to tsv"
qiime tools export --input-path otu_18S_taxonomy.qza --output-path 18S_taxonomy_tsv

#done
