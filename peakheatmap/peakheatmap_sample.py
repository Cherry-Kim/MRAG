import os
## Install deepTools ###
#git clone https://github.com/deeptools/deepTools.git
#cd deepTools
#pip install .

def Step1_bw(bam_DIR, sample):
    bw = sample + '.CPM.bw'

    cmd = f'bamCoverage --bam {bam_DIR}{sample}.dedup.bam \
    --normalizeUsing CPM \
    -o {bw} \
    --effectiveGenomeSize 2308125299 \
    --exactScaling \
    --binSize 10 \
    --extendReads 100'

    os.system(cmd)


def Step2_matrix():
    bw_files = [
    '1_pLenti_ONECUT2_OE_1a_merge.CPM.bw',
    '2_pLenti_ONECUT2_OE_1b_merge.CPM.bw',
    '3_pLenti_cMyc_DDK_2a_merge.CPM.bw',
    '4_pLenti_cMyc_DDK_2b_merge.CPM.bw'
    ]
    bw_string = ' '.join(bw_files)

    cmd = f'computeMatrix reference-point --referencePoint center \
            -R /BIO1/MRAG/HYK/Plot/CistromeDB/up.bed /BIO1/MRAG/HYK/Plot/CistromeDB/down.bed \
            -S {bw_string} \
            -b 3000 -a 3000 \
            --skipZeros \
            -p 96 \
            -o OC2vsCyc_matrix.gz'
    
    os.system(cmd)

def Step3_plot():
    cmd = f'plotHeatmap -m OC2vsCyc_matrix.gz \
            -out OC2vsCyc_heatmap.png \
            --regionsLabel UP DOWN \
            --xAxisLabel ""  \
            --samplesLabel \
            ONECUT2_1a \
            ONECUT2_1b \
            cMyc_2a \
            cMyc_2b\
            --colorList \
            "white,red" \
            "white,red" \
            "white,darkblue" \
            "white,darkblue"'
    
    os.system(cmd)


def Main():
    bam_DIR = '/BIO1/MRAG/HYK/ATAC/alignments/'
    files = ['1_pLenti_ONECUT2_OE_1a_merge', '2_pLenti_ONECUT2_OE_1b_merge' ]
    #files = ['1_pLenti_ONECUT2_OE_1a_merge', '2_pLenti_ONECUT2_OE_1b_merge', '3_pLenti_cMyc_DDK_2a_merge', '4_pLenti_cMyc_DDK_2b_merge']
    for sample in files:
        print(sample)
        #Step1_bw(bam_DIR, sample)
    #Step2_matrix()
    Step3_plot()
Main()
