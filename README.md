# Analysis of genome-wide CRISPR screen in primary human T cells

Here I analyze [data](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA487352) from a [genome-wide CRISPR screen in primary human T cells](https://www.sciencedirect.com/science/article/pii/S0092867418313333) downloaded from the NCBI Sequence Read Archive (SRA) database. Primary human CD8+ T cells were isolated from four donors, stimulated, transduced with the Brunello single-guide RNA (sgRNA) library, electroporated with Cas9, and then restimulated. Cells were then sorted into 'non-proliferating' cells and highly proliferating cells followed by sequencing of the two populations.

The analysis happens in two notebooks: a Jupyter Python notebook for going from the raw FASTQ reads to counts, followed by an analysis in the R/Bioconductor notebook of the counts.

The steps followed in the Jupyter Python notebook are:
- downloading of the raw FASTQ data from the SRA database ([prefetch and fasterq-dump](https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump))
- quantifying the abundance of the 77,441 sgRNAs in the 8 samples (4 donors, 2 cell populations each)
- statistical tests to assess the significance of sgRNA or gene knockout effects between 'non-proliferating' cells and highly proliferating cells

[MAGeCK](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0554-4) is used for the last two steps.

Downstream analysis of the output of MAGeCK is subsequently performed in the R/Bioconductor notebook (and corresponding script). The script performs quality assessment, conducts differential analysis, and performs Gene Set Enrichment Analysis (GSEA) on the MAGeCK output.

### Required Input Files:
- per-sample count summary information from MAGeCK count.
- gene summary file from MAGeCK test.
- sgRNA summary file from MAGeCK test.

### Required Packages:
- `MAGeCKFlute`: For analysis of MAGeCK output.
- `clusterProfiler`: For conducting GSEA.
- `ggplot2`, `ggrepel`: For data visualization.

### Outputs:
- Various visualizations such as bar plots, scatter plots (volcano, rank, etc.), and GSEA plots.
- Summary tables of enrichment analysis results.
- Diagnostic plots for assessing the quality of gene sets and differential gene expression.

### Steps Included:
1. Quality checks on count data (Gini index, missed sgRNAs, read mapping).
2. Exploratory data analysis using `ScatterView` for differential gene expression.
3. Visualization of positively and negatively selected genes.
4. Examination of sgRNA rankings targeting top-selected genes.
5. GSEA analysis using the KEGG database.
6. Visualization of GSEA results through various plots.
