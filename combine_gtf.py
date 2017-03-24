import os

rRNA_gtf = "path/to/hg19.repeatMasker.rRNA.gtf"
tRNA_gtf = "path/to/hg19_UCSC_tRNA.gtf"
known_gene_gtf = "../TxDb_hg19_gene.gtf"
combined_gtf = "hg19_combined_gene_YYMMDD.gtf"

#to specifically download rRNA Repeat Masker annotations:
#go to Table Browser (https://genome.ucsc.edu/cgi-bin/hgTables) for your genome of interest (I've tested hg19)
#go to "Repeats" and select "Repeat Masker"
#table should be set to "rmsk"
#select "filter", and set "repClass does match rRNA" (you have to type in "rRNA")
#solution courtesy of http://seqanswers.com/forums/archive/index.php/t-41868.html

#for tRNA, using Table Browser, go to "Genes and Gene Predictions" and select "tRNA genes" as the track

#for known gene .gtf creation, please see https://github.com/cwarden45/RNAseq_templates/tree/master/Genome_Ref_Code

command = "cat " + rRNA_gtf + " " + tRNA_gtf + " " + known_gene_gtf + " > " + combined_gtf
os.system(command)