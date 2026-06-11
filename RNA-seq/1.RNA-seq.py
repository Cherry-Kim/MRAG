#!/usr/bin/env python3
import glob
import os

REF_DIR = './Reference/hg38/'
GENOME = 'GRCh38'
num_threads = 96

def step1_preprocessing(sample, num_threads, out_dir = 'fastqc'):
    os.system('/BIO1/MRAG/HYK/ATAC/TrimGalore-0.6.11/trim_galore --paired  --gzip -o trim_galore/ '+sample+'_1.fastq.gz '+sample+'_2.fastq.gz')

    os.makedirs(out_dir, exist_ok=True)
    os.system('fastqc -t '+str(num_threads)+' --nogroup -o '+out_dir+' trim_galore/'+sample+'_1_val_1.fq.gz')
    os.system('fastqc -t '+str(num_threads)+' --nogroup -o '+out_dir+' trim_galore/'+sample+'_2_val_2.fq.gz')

def step2_STAR(sample, num_threads, REF_DIR):
    os.system('/BIO1/MRAG/HYK/RNA/STAR-2.7.11b/bin/Linux_x86_64/STAR --runThreadN '+str(num_threads)+' --genomeDir '+REF_DIR+' --readFilesIn trim_galore/'+sample+'_1_val_1.fq.gz trim_galore/'+sample+'_2_val_2.fq.gz --readFilesCommand zcat --outTmpDir '+'TEMP_'+sample+' --outFileNamePrefix '+sample+' --quantMode TranscriptomeSAM')
    os.remove(sample+'Aligned.out.sam')

def step3_RSEM(sample, REF_DIR, GENOME, num_threads):
    os.system('/BIO1/MRAG/HYK/RNA/RSEM/rsem-calculate-expression --paired-end --bam -p '+str(num_threads)+' '+sample+'Aligned.toTranscriptome.out.bam '+REF_DIR+GENOME+' '+sample)
 

def Main():
    file_list = os.listdir('./')
    f_list=[file for file in file_list if file.endswith("_1.fastq.gz")]
    for fname in f_list:
        sample = fname.replace('_1.fastq.gz','')
        print(sample)
        #step1_preprocessing(sample, num_threads)
        #step2_STAR(sample, num_threads, REF_DIR)
        step3_RSEM(sample, REF_DIR, GENOME, num_threads)
Main()
