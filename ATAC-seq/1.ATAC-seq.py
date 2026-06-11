#!/usr/bin/env python3
import glob
import os

REF_DIR = './Reference/hg38/'
ALIGN_DIR = './alignments/'
GENOME = 'hg38'
GTF_FILE = os.path.join(REF_DIR, 'hg38.refGene.gtf')
num_threads = 96

def step0_preprocessing(sample, num_threads, out_dir = 'fastqc'):
    os.system('trim_galore --paired  --gzip -o trim_galore/ '+sample+'_R1.fastq.gz '+sample+'_R2.fastq.gz')
    os.makedirs(out_dir, exist_ok=True)
    os.system('fastqc -t '+str(num_threads)+' --nogroup -o '+out_dir+' trim_galore/'+sample+'_R1_val_1.fq.gz')
    os.system('fastqc -t '+str(num_threads)+' --nogroup -o '+out_dir+' trim_galore/'+sample+'_R2_val_2.fq.gz')

def step1_align(sample, REF_DIR, GENOME, num_threads, ALIGN_DIR):
    os.makedirs(ALIGN_DIR, exist_ok=True)
    os.system('bowtie2 -p '+str(num_threads)+' -x '+REF_DIR+GENOME+' -1 trim_galore/'+sample+'_R1_val_1.fq.gz -2 trim_galore/'+sample+'_R2_val_2.fq.gz -S '+ALIGN_DIR+sample+'.sam 2>'+ALIGN_DIR+sample+'.bowtie2Log')
    os.system('samtools view -bS '+ALIGN_DIR+sample+'.sam > '+ALIGN_DIR+sample+'.bam') 
    os.system('samtools sort -@ '+str(num_threads)+' '+ALIGN_DIR+sample+'.bam > '+ALIGN_DIR+sample+'.sorted.bam')
    os.system('samtools index '+ALIGN_DIR+sample+'.sorted.bam')
    os.remove(ALIGN_DIR+sample+'.sam')

def step2_filter(sample, REF_DIR, GENOME, ALIGN_DIR, picard_dir = 'TEMP_PICARD/'):
    os.system('bamtools filter -in '+ALIGN_DIR+sample+'.sorted.bam -out '+ALIGN_DIR+sample+'.filterd.bam -mapQuality ">=30" -isProperPair true')
    os.system('samtools view '+ALIGN_DIR+sample+'.filterd.bam | egrep -v chrM | samtools view -bT '+REF_DIR+GENOME+'.fa -o '+ALIGN_DIR+sample+'.filterd.-chrM.bam')

    os.makedirs(picard_dir, exist_ok=True)
    os.system('java -Xmx8g -jar /BIO1/picard.jar MarkDuplicates TMP_DIR='+picard_dir+' VALIDATION_STRINGENCY=LENIENT I='+ALIGN_DIR+sample+'.filterd.-chrM.bam O='+ALIGN_DIR+sample+'.dedup.sam M='+sample+'.duplicate_metrics REMOVE_DUPLICATES=true AS=true')
    os.system('java -Xmx8g -jar /BIO1/picard.jar SortSam TMP_DIR='+picard_dir+' VALIDATION_STRINGENCY=LENIENT SO=coordinate I='+ALIGN_DIR+sample+'.dedup.sam O='+ALIGN_DIR+sample+'.dedup.bam CREATE_INDEX=true')
    os.system('mv *_metrics '+ALIGN_DIR)
    os.remove(ALIGN_DIR+sample+'.dedup.sam')
   
    os.system('bamCoverage -b '+ALIGN_DIR+sample+'.dedup.bam -o '+sample+'.bw')

def step3_peakcall(sample, ALIGN_DIR, REF_DIR, SP = 'hs'):
    os.system('bedtools bamtobed -i '+ALIGN_DIR+sample+'.dedup.bam  > '+ALIGN_DIR+sample+'.bed')
    os.system('macs2 callpeak -t '+ALIGN_DIR+sample+'.bed -n '+sample+' -g '+SP+' -B -q 0.05 --outdir macs2_'+sample+' --nomodel --shift -100 --extsize 200')

    os.system('bedtools intersect -v -a macs2_'+sample+'/'+sample+'_peaks.narrowPeak -b '+REF_DIR+'GRCh38_unified_blacklist.bed > '+sample+'_peaks.narrowPeak.filt')

def step4_peak_annotation(sample, REF_DIR, GTF_FILE):
    os.system('awk -F "\t" \'{print $4 "\t" $1 "\t" $2 "\t" $3 "\t" "+"}\' '+sample+'_peaks.narrowPeak.filt > '+sample+'.input.bed')
    os.system('perl /BIO1/MRAG/HYK/homer/bin/annotatePeaks.pl '+sample+'.input.bed '+REF_DIR+' -gtf '+GTF_FILE+' > '+sample+'_annotations.txt')

def step5_DASTK(num_threads, REF_DIR, GENOME):
    treatment_dir = "DASTK/treatment"
    control_dir = "DASTK/control"
    os.makedirs(treatment_dir, exist_ok=True)
    os.makedirs(control_dir, exist_ok=True)

    treatment = ["1_pLenti_ONECUT2_OE_1a_merge_peaks.narrowPeak.filt",
                 "2_pLenti_ONECUT2_OE_1b_merge_peaks.narrowPeak.filt"
                        ]
    control = [  "3_pLenti_cMyc_DDK_2a_merge_peaks.narrowPeak.filt",
                 "4_pLenti_cMyc_DDK_2b_merge_peaks.narrowPeak.filt"
                        ]

    os.system(f"cat {' '.join(treatment)} > {treatment_dir}/treatment.narrowPeak.filt")
    os.system(f"cat {' '.join(control)} > {control_dir}/control.narrowPeak.filt")
        
    os.system('python ./DAStk/DAStk/process_atac.py -t '+str(num_threads)+' -e '+treatment_dir+'/treatment.narrowPeak.filt  -m '+REF_DIR+'grch38_motifs/ -g '+GENOME+' -o '+treatment_dir)
    os.system('python ./DAStk/DAStk/process_atac.py -t '+str(num_threads)+' -e '+control_dir+'/control.narrowPeak.filt  -m '+REF_DIR+'grch38_motifs/ -g '+GENOME+' -o '+control_dir)

    os.system('python ./DAStk/DAStk/differential_md_score.py -1 '+treatment_dir+'/treatment.narrowPeak_md_scores.txt -2 '+control_dir+'/control.narrowPeak_md_scores.txt   -m "treatment" -n "control" -b -o DASTK')


def Main():
    fastp_files = glob.glob('*_R1.fastq.gz')
    for fname in fastp_files:
        sample = fname.replace('_R1.fastq.gz', '')
        print(sample)
        step0_preprocessing(sample, num_threads)
        step1_align(sample, REF_DIR, GENOME, num_threads, ALIGN_DIR)
        step2_filter(sample, REF_DIR, GENOME, ALIGN_DIR)
        step3_peakcall(sample, ALIGN_DIR, REF_DIR)
        step4_peak_annotation(sample, REF_DIR, GTF_FILE)
    step5_DASTK(num_threads, REF_DIR, GENOME)
Main()
