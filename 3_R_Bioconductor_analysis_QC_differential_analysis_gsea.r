# load packages
library(MAGeCKFlute)
library(clusterProfiler)
library(ggplot2)
library(ggrepel)

set.seed(171336)

# increase size of plots from default
options(repr.plot.width = 14, 
        repr.plot.height = 14) # from 7, 7

## QC
# count summary from mageck count
file1 <- "gw01.countsummary.txt"
countsummary <- read.delim(file1, check.names = FALSE)
countsummary

# Gini index - smaller value indicates more eveness of the count distribution. 
# (Recommended: around 0.1 for plasmid or initial state samples, and around 0.2-0.3 for negative selection samples)
BarView(countsummary, x = "Label", y = "GiniIndex",
        ylab = "Gini index", main = "Evenness of sgRNA reads")

# Missed sgRNAs: sgRNAs that have 0 counts (Recommended: no more than 1%)
countsummary$Missed = 100*(countsummary$Zerocounts/countsummary$TotalsgRNAs)
BarView(countsummary, x = "Label", y = "Missed", fill = "#394E80",
        ylab = "% missed gRNAs", main = "% Missed sgRNAs")

# Read mapping
# Reads: (Recommended: 100~300 times the number of sgRNAs)
# mapping percentage: calculated as Mapped/Reads (Recommended: at least 60%)
MapRatesView(countsummary)

# read data from mageck test
## path to the gene summary file
file1 <- "gw01_test.gene_summary.txt"
## path to the sgRNA summary file 
file2 <- "gw01_test.sgrna_summary.txt"
gene_data <- ReadRRA(file1)
head(gene_data)
sgrna_data <- ReadsgRRA(file2)
head(sgrna_data)

# volcano plot
gene_data$LogFDR <- -log10(gene_data$FDR)
p1 <- ScatterView(gene_data, x = "Score", y = "LogFDR", label = "id", 
                 model = "volcano", top = 10, 
                size = 3)
print(p1)

# Rank all the genes based on their scores and label genes in the rank plot
gene_data$Rank <- rank(gene_data$Score)
p1 <- ScatterView(gene_data, x = "Rank", y = "Score", label = "id", 
                 top = 10, auto_cut_y = TRUE, ylab = "Log2FC", 
                 groups = c("top", "bottom"), size = 3)
print(p1)

# positively selected genes
gene_data$RandomIndex <- sample(1:nrow(gene_data), nrow(gene_data))
gene_data <- gene_data[order(-gene_data$Score), ]
gg <- gene_data[gene_data$Score > 0, ]
p1 <- ScatterView(gg, x = "RandomIndex", y = "Score", label = "id",
                 y_cut = CutoffCalling(gene_data$Score,2), 
                 groups = "top", top = 10, ylab = "Log2FC", size = 4)
p1

# negatively selected genes
gg <- gene_data[gene_data$Score < 0, ]
p2 <- ScatterView(gg, x = "RandomIndex", y = "Score", label = "id",
                 y_cut = CutoffCalling(gene_data$Score, 2), 
                 groups = "bottom", top = 10, ylab = "Log2FC", size = 4)
p2

# visualize the rank of sgRNAs targeting top selected genes.
p2 <- sgRankView(sgrna_data, top = 10, bottom = 10)
print(p2)

# create weighted gene list with scores
geneList = gene_data$Score
names(geneList) = gene_data$id
geneList

# GSEA using KEGG database; limiting to gene sets size 15-500
gseRes <- enrich.GSE(geneList, type = "KEGG", limit = c(15, 500), eps = 0)
gseRes@result

# look at top 25 KEGG hits
gseRes@result[1:25, c('Description', 'NES', 'p.adjust', 'Count')]

## visualize enrichment results
# barplot
df <- gseRes@result
df$logFDR <- -log10(df$p.adjust)
p <- BarView(df[1:10, ], "Description", 'logFDR')
p <- p + labs(x = NULL) + coord_flip()
p

# dot plot
EnrichedView(gseRes, bottom = 10, mode = 1)

#gseaplot - T cell receptor signaling pathway
gseaplot(gseRes, geneSetID = 7, title = gseRes$Description[7])

#gseaplot - PD-L1 expression and PD-1 checkpoint pathway in cancer
gseaplot(gseRes, geneSetID = 6, title = gseRes$Description[6])

#gseaplot - Th17 cell differentiation
gseaplot(gseRes, geneSetID = 13, title = gseRes$Description[13])

#gseaplot - Th1 and Th2 cell differentiation
gseaplot(gseRes, geneSetID = 14, title = gseRes$Description[14])

sessionInfo()
