input.file = "sample_description.txt"
output.file = "sample_description_v2.txt"
alignment.folder = "../Bowtie_Alignment"

input.table = read.table(input.file, head=T, sep="\t")

miRNA.htseq.files = paste(alignment.folder,"/",input.table$sampleID,"/HTseq_",input.table$sampleID,"_mature_miRNA_counts.txt",sep="")
other.htseq.files = paste(alignment.folder,"/",input.table$sampleID,"/HTseq_",input.table$sampleID,"_other_gene_counts.txt",sep="")

miRNA.featureCounts.files = paste(alignment.folder,"/",input.table$sampleID,"/featureCounts_",input.table$sampleID,"_mature_miRNA_counts.txt",sep="")
other.featureCounts.files = paste(alignment.folder,"/",input.table$sampleID,"/featureCounts_",input.table$sampleID,"_other_gene_counts.txt",sep="")

output.table = data.frame(input.table, miRNA.htseq.files, other.htseq.files, miRNA.featureCounts.files, other.featureCounts.files)
write.table(output.table,output.file,sep="\t", quote=F, row.names=F)