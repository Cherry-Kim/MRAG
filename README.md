# MRAG: An integrative data analysis pipeline for modeling transcriptional regulatory networks in accessible genome using ATAC-seq and RNA-seq data


Modeling transcriptional regulatory networks (TRNs) in response to a perturbation provides insights into cellular behavior by describing interactions between transcription factors (TFs) and their transcriptional targets. Here, we present **MRAG (Master Regulator analysis of Accessible Genome)**, which reconstructs transcriptional-regulatory networks the integration of the assay for transposase-accessible chromatin with high-throughput sequencing (ATAC-seq) and RNA sequencing (RNA-seq). MRAG first identifies enriched TF motifs within a list of accessible chromatin regions from multiple replicates and experimental conditions, as well as differential expression by a perturbation. TF activity score is computed by combining differential motif score (ATAC-seq) and TF-target enrichment score (RNA-seq) to selected master TFs. MRAG also prioritizes potential gene targets of the master TFs with differential ATAC-seq peaks on the regulatory region with master TFs’ motif coupled with significant differential expression by a perturbation. As an integrative analysis tool of chromatin accessibility and gene expression profiles, application of MRAG will improve TRN inference in humans and animal models.

![image](https://github.com/Cherry-Kim/MRAG/assets/64776690/ad8b9b28-5068-42ca-a23e-6e4c60157fff)

 
## Installation
CellNeighborEX requires Python version >=3.8, <3.11. We recommend using conda environment to avoid dependency conflicts. The dependencies are listed in requirements.txt.
```
conda create -n myenv python=3.10
conda activate myenv
```
```
function test() {
  console.log("notice the blank line before this function?");
}
```

Python API documentation and tutorials
Please see this Read the Docs.



