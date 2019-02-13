### Acknowledgements ###

Code developed when providing analysis support from the City of Hope Integrative Genomics Core with requests / suggestions from *Yuan Yuan MD/PhD* (re-analysis with different scripts after other IGC staff), *Susan Neuhausen* (early re-analysis with different scripts, after/before other IGC staff), and *[External Collaborator B]* (currently, did not want to be explicitly acknowledged here).

I would also like to thank the following users that gave us permission to acknowledge these templates (or scripts similar to these templates) were used for analysis of their data: *Luke Gruenert* ([Onconetics](http://onconetics.com/), collaboration with David Salzman at [sRNAlytics](https://www.srnalytics.com/) helped me realize original, untrimmed reads had barcode to alleviate ligation bias, as described in [Sun et al. 2011](https://www.ncbi.nlm.nih.gov/pubmed/22016383)), *Shiuan Chen* (collaboration with Susan Neuhausen), *Rahul Jandial* (with lab member Khairul Ansari), and *Sharad Singhal*.  There was also 1 other external collaborator that requested to not be acknowledged here (for this particular data type, at this time).

There was a previous IGC function that I did use for quantification, and that is why I say “Among templates provided here, I would recommend using run_featureCounts.py.”  However, I did end up using featureCounts for quantification in some of the later projects.

In general, some projects have been analyzed or re-analyzed with different staff members (sometimes in different labs/institutes, sometimes within the IGC).  So, the use of data with modified versions of these templates may not be what is eventually used in the associated publication.

Xiwei Wu (the Director of the Integrative Genomics Core) was kind enough to review this acknowledgement, and he approves of the content.

### \~Temporary Note\~ ###
**I apologize for the confusion, but I would like to emphasize that these are called “templates” because I almost always have to modify the code for each project (beyond the parameter files), meaning they will be more difficult for other people to use in the same way.**  This was unfortunately not immediately clear to me when I created the templates.

I also believe that the process of writing the scripts for analysis (such as the templates) is very important for the learning process, and it is very important that you understand all the steps for analysis before presenting them in a paper.

So, I will post an update when more specific guidance / suggestions can be provided.  **In the meantime, I apologize, but I cannot provide user support for the templates.**

# miRNAseq_template
Template for analysis of Illumina miRNA-Seq data

### Order to Run Scripts ###

#### 0) Prepare miRbase reference GTF

If you wish to directly compare alignments for miRNA expression for an earlier build (such as hg19),
you'll need to liftOver the annotation file (which you can do using `liftOver_GFF.py`).

The script `extract_mature_miRNA.py` will filter precursor miRNA and convert the GFF to a GTF (with precurser ID 
in gene_id, mature gene id for transcript_id, and the gene name in gene_name).  If you perform quantification using HT-Seq,
you need to filter the primary transcripts to avoid flagging the mature miRNAs as ambiguous (some miRNA have identical sequences
in multiple locations, but the script will give you an idea of how many sequence you'd be excluding if you require unique mappings).
For some programs, you may also need annotations in GTF format rather than GFF format.  Also, by default, Bowtie wil only report 1
alignment per read (but "-m 1" would need to be added for that alignment to be unique).

mirBase also provides FASTA files for "high-confidence" miRNA.  You can use `extract_mature_miRNA_high_confidence.py`
to create a .gtf file that only contains the "high-confidence" miRNA for your mapping.

Reference preparation scripts don't use parameter file, so you'll need to edit the code directly.

#### 1) Trim first 3 nt and remove adapters using `cutadapt_filter.py` or `cutadapt_filter_cluster.py`.

#### 2) Calculate cutadapt statistics using `cutadapt_filter_statistics.py`

You'll need to parse FastQC output to calculate the cutadapt statistics.  If those files are not already available, they can be produced using `cluster_FastQC.py`.

#### 3) Align reads using Bowtie `align_Bowtie_cluster.py`

#### 4) Calculate alignment stats using `collect_flagstat_stats.py`

#### 5) Quantify read counts.  Among templates provided here, I would recommend using `run_featureCounts.py`.

These scripts are slightly different than those used in [TopHat_Workflow](https://github.com/cwarden45/RNAseq_templates/tree/master/TopHat_Workflow).  More specifically, only counts are tabulated (without RPKM values) and two types of counts are provided: mature miRNA and "other" genes.  For miRNAs, HT-Seq is switched to "intersection-strict" and featureCounts is run with "-M --fracOverlap 1".  The script `combine_gtf.py` can create the "other" gene table as a mix of UCSC RepeatMasker rRNA (with download instructions in comments), UCSC tRNA, and the gene annotations from [Genome_Ref_Code](https://github.com/cwarden45/RNAseq_templates/tree/master/Genome_Ref_Code) (includes primary miRNA and snoRNAs, for example).  All quantifications assume gene on same strand as read. 

#### 6) Calculate Count-Per-Million (CPM) abundance.  Among templates provided here, I would recommend using `featureCounts_reformat.R`.  Downstream analysis can work with any CPM table formatted with the miRNA in the first column and the short sample ID (desired sample label, "userID" in sample description file) in the remaining columns.

#### 7) QC plots using `qc.R`

#### 8) Differential expression using `DEG_miRNA.R`

### Dependencies (some optional) ###

Most Python scripts can be run using this [Docker image](https://hub.docker.com/r/cwarden45/rnaseq-dependencies/)

*Adapter Removal*

cutadapt: http://cutadapt.readthedocs.io/en/stable/guide.html

*Alignment*

Bowtie: http://bowtie-bio.sourceforge.net/index.shtml

*Read Counting*

HTseq: http://www-huber.embl.de/HTSeq/doc/install.html#install

featureCounts: http://bioinf.wehi.edu.au/featureCounts/

GenomicAlignments: https://bioconductor.org/packages/release/bioc/html/GenomicAlignments.html

*Differential Expression*

edgeR: https://bioconductor.org/packages/release/bioc/html/edgeR.html

limma-voom: https://bioconductor.org/packages/release/bioc/html/limma.html

DESeq2: https://bioconductor.org/packages/release/bioc/html/DESeq2.html

qvalue: https://bioconductor.org/packages/release/bioc/html/qvalue.html

*Visualization*

gplots: https://cran.r-project.org/web/packages/gplots/index.html

RColorBrewer: https://cran.r-project.org/web/packages/RColorBrewer/index.html

heatmap.3: https://github.com/obigriffith/biostar-tutorials/blob/master/Heatmaps/heatmap.3.R

heatmap.3 example: https://www.biostars.org/p/18211/

IGV / igvtools: http://software.broadinstitute.org/software/igv/

*Reference Build Conversion*

liftOver command line: https://genome-store.ucsc.edu/

genome build chain files (for genome of interest): http://hgdownload.cse.ucsc.edu/downloads.html

*miRNA Targets*

TargetScan: http://www.targetscan.org/vert_71/

tarBase: http://diana.imis.athena-innovation.gr/DianaTools/index.php?r=tarbase/index

IPA microRNA Target Filter: https://www.qiagenbioinformatics.com/products/features/microrna-target-filter/

*Miscellaneous*

FastQC: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

Biopython: http://biopython.org/wiki/Biopython

### Parameter Values ###
| Parameter | Value|
|---|---|
|comp_name	| Name of differential expression comparison (used to name output file)
|plot_groups | Names of columns in *sample_description_file* to be plotted in QC and differential expression plots.  Use commas to plot multiple groups|
|deg_groups|Names of columns in *sample_description_file* to be plotted in QC and differential expression plots.  Use commas to include multiple variables (for multivariate model or gene list filtering)|
|treatment_group|Treatment group for primary variable; enter *continuous* for a continuous variable and a correlation will be provided instead of a fold-change value.|
|Raw_Code_PC|Path to output folder for most results|
|Result_Folder|Path to output folder for selected, final results|
|Alignment_Folder|Path to TopHat Alignments|
|Reads_Folder|Path to Reads for TopHat Alignment|
|Cluster_Email|If running alignment on a cluster, e-mail for notifications|
|pvalue_method|Method to Calculate P-value.  Can be *edgeR*, *limma-voom*, *DESeq2*, *lm* (linear regression), or *aov* (ANOVA)|
|fdr_method|Method to Calculate FDR.  Can be *BH* (Benjamini and Hochberg),*q-value*, or *q-lfdr*|
|genome|Name of genome build|
|Bowtie_Ref| Path to Bowtie ref|
|Threads|Number of Threads for TopHat Alignment|
|miRNA_GTF|miRbase mature miRNA GTF|
|other_GTF|GTF with quantifications for other genes in separate quantification (no overlap between mature and primary miRNA annotations)|
|sample_description_file|Name of Sample Description File|
|total_counts_file|Name of File to Contain Total Read Counts and Cutadapt-Filtered Read Counts|
|aligned_stats_file|Name of File to Contain Aligned Read Counts|
|cluster_distance| Distance metric for dendrogram.  Can be *Euclidean* or *Pearson_Dissimilarity*|
|cpm_file|Name of File to Contain Rounded CPM Expression Values|
|CPM_norm|How to count number of aligned reads: *aligned*, *miRNA*, or *quantified-other*|
|counts_file|Name of File to Contain Read Counts Per Gene|
|cpm_expression_cutoff|Minimum Robust Expression Level (used for filtering genes for differential expression and calculating robust miRNA percent)|
|minimum_fraction_expressed|Minimum fraction of samples with expression above *cpm_expression_cutoff*. Filter for differential expression analysis.|
|fold_change_cutoff|Minimum fold-change difference to consider a gene differentially expressed.  For miRNA-Seq, fold-change values are calculated from log2(CPM + 1) values, and *cpm_expression_cutoff* is only used for filtering miRNAs.|
|cor_cutoff|If using a continuous variable, minimum absolute correlation to consider a gene differentially expressed|
|pvalue_cutoff|Maximum p-value to consider a gene differenitally expressed|
|fdr_cutoff|Maximum FDR to consider a gene differentially expressed|
|interaction| Method for comparing an interaction of two variables.  Can be *model* or *no*|
|secondary_plot|Set to "continuous" if you wish to plot the second variable as a continuous variable in the differentially expressed miRNA heatmap|
