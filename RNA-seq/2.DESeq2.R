#!/usr/bin/env Rscript
Step1_tximport <- function(){	
	library(tximport)
	files <- c("SRR26175961.genes.results", "SRR26175962.genes.results", "SRR26175963.genes.results", "SRR26175964.genes.results")
	names(files) <- c("SRR26175961", "SRR26175962", "SRR26175963", "SRR26175964")
#                SRR26175961                 SRR26175962
#"SRR26175961.genes.results" "SRR26175962.genes.results"

	txi <- tximport(files,
                type = "rsem",  txIn = FALSE, txOut = FALSE)
        write.table(txi$counts, "txi.rsem.counts.txt", col.names=NA, row.names=T, quote=F,sep='\t')
}


Step2_DESeq2 <- function(){
	library(DESeq2)
	library(tibble)
	rawCounts <- txi$counts
	rawCounts <- as.data.frame(rawCounts)
	rawCounts <- rownames_to_column(rawCounts, var = "GENE")
	#> head(rawCounts)
        #        GENE SRR26175961 SRR26175962 SRR26175963 SRR26175964
	#1 ENSG00000000003.15       52.00        43.0      129.00      142.00

	A <- rawCounts$GENE
	B <- rawCounts[, 2:ncol(rawCounts)]
	C <- cbind(A,B)
	counts <- as.matrix(C[,-1])
	rownames(counts) <- C[,1]

	Sample <- factor(c("SRR26175961", "SRR26175962",
                   "SRR26175963", "SRR26175964"))
	Group <- factor(c("T","T","N","N"))
	Group <- relevel(Group, ref = "N")
	sampleData <- data.frame(Sample=Sample, Group)

	dds <- DESeqDataSetFromMatrix(
	  countData = round(counts),
	  colData = sampleData,
	  design = ~ Group
	)
	dds <- DESeq(dds)

	dds_results <- results(dds)
	dds_results <- dds_results[order(dds_results$padj),]
	write.csv(dds_results, file="DESeq2.csv", quote=F)
}

Step3_gene <- function(){
	library('biomaRt')
	dds_results <- read.csv('DESeq2.csv')
        genes <-  dds_results$X
	genes_clean <- sub("\\..*", "", genes)
	
	library(org.Hs.eg.db)
	library(AnnotationDbi)
	gene_symbol <- mapIds(
	  org.Hs.eg.db,
	  keys = genes_clean,
	  column = "SYMBOL",
	  keytype = "ENSEMBL",
	  multiVals = "first"
	)
        dds_results$GENE <- gene_symbol[genes_clean]
	dds_results <- dds_results[, c("X", "GENE", setdiff(colnames(dds_results), c("X", "GENE")))]
	colnames(dds_results)[colnames(dds_results) == "X"] <- "ID"
        write.csv(dds_results, "DESeq2.gene.csv", quote = FALSE, row.names = FALSE)
}
