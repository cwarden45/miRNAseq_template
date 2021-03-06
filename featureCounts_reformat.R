parameter.file = "parameters.txt"

##### Ideally, you don't have to modify code below this point #####

normalizeTotalExpression = function (geneExpr, totalReads) {
	return(geneExpr / totalReads)
}#end def normalizeTotalExpression

count.defined.values = function(arr, expr.cutoff)
{
	sig.figures = 1
	if (expr.cutoff > 0)
		sig.figures = 0
	expr.cutoff = round(expr.cutoff, digits=sig.figures)
	arr = round(arr, digits=sig.figures)
	return(length(arr[arr > expr.cutoff]))
}#end def count.values

trimmed.counts = function(counts, min.percent, max.percent)
{
	total.counts = sum(counts)
	counts = counts[order(counts)]
	min.index = min.percent * length(counts)
	max.index = max.percent * length(counts)
	counts = counts[min.index:max.index]
	trimmed.counts = sum(counts)
	trimmed.percent = round(100 * trimmed.counts/total.counts, digits=1)
	return(trimmed.percent)
}#end def trimmed.counts

param.table = read.table(parameter.file, header=T, sep="\t")
comp.name=as.character(param.table$Value[param.table$Parameter == "comp_name"])
genome=as.character(param.table$Value[param.table$Parameter == "genome"])
min.expression = as.numeric(as.character(param.table$Value[param.table$Parameter == "cpm_expression_cutoff"]))
user.folder = as.character(param.table$Value[param.table$Parameter == "Result_Folder"])
sample.file = as.character(param.table$Value[param.table$Parameter == "sample_description_file"])
total.reads.file = as.character(param.table$Value[param.table$Parameter == "total_counts_file"])
aligned.stat.file = as.character(param.table$Value[param.table$Parameter == "aligned_stats_file"])
counts.file = as.character(param.table$Value[param.table$Parameter == "counts_file"])
cpm.file = as.character(param.table$Value[param.table$Parameter == "cpm_file"])
aligned.type = as.character(param.table$Value[param.table$Parameter == "CPM_norm"])

sample.table = read.table(sample.file, header=T, sep="\t")
sampleID = as.character(sample.table$sampleID)
sample.label = as.character(sample.table$userID)
dash.flag = grep("-",sample.label)
if(length(dash.flag) > 0){
	print(paste(paste(sample.label[dash.flag],collapse=",")," samples labels have dashes in their labels",sep=""))
}
num.flag = grep("^[0-9]",sample.label)
if(length(num.flag) > 0){
	print(paste(paste(sample.label[num.flag],collapse=",")," samples labels start with numbers",sep=""))
}

if((length(dash.flag) > 0)|(length(num.flag) > 0)){
	stop("Please make sure sample labels do not have dashes and do not start with numbers")
}

total.reads.table = read.table(total.reads.file, header=T, sep="\t")
totalID = as.character(total.reads.table$SampleID)
total.reads.table = total.reads.table[match(sampleID, totalID),]
print(total.reads.table$TotalReads)

aligned.stat.table = read.table(aligned.stat.file, header=T, sep="\t")
alignedID = as.character(aligned.stat.table$SampleID)
aligned.stat.table = aligned.stat.table[match(sampleID, alignedID),]
total.reads = as.numeric(total.reads.table$TotalReads)
percent.cutadapt.pass = as.character(total.reads.table$Percent.Cutadapt.Filtered)
percent.aligned.reads = as.character(aligned.stat.table$Percent.Aligned)

#other counts
other.counts.file = gsub(".txt$","_other.txt",counts.file)

count.files = as.character(sample.table$other.featureCounts.files)

temp.file = count.files[[1]]
temp.table = read.table(temp.file, sep="\t", header=T)
genes = as.character(temp.table[[1]])

other.count.mat = matrix(nrow=nrow(temp.table), ncol=length(sampleID))
colnames(other.count.mat) = sample.label

matched.genes = c()

for (i in 1:length(sampleID)){
	temp.file = count.files[[i]]
	temp.table = read.table(temp.file, sep="\t", header=T)
	temp.genes = as.character(temp.table[[1]])
	
	if(length(matched.genes) != 0){
		matched.genes = temp.genes[match(matched.genes, temp.genes, nomatch=0)]
	} else if(!identical(temp.genes, genes)){
		print("Genes not in same order - most likely, different quantification method was used for different samples")
		userAns = readline(prompt="Do you wish to proceed with a subset of matched gene symbols? (y/n): ")
		userAns = tolower(substr(userAns, 1, 1))
		if (userAns != "y"){
			stop("Please re-run featureCounts with all of your samples")
		} else {
			matched.genes = temp.genes[match(genes, temp.genes, nomatch=0)]
		}#end else
	}else {
		other.count.mat[,i] = temp.table[[7]]
	}#end else
}#end for (i in 1:length(sampleID))

if (length(matched.genes) != 0){
	other.count.mat = matrix(nrow=length(matched.genes), ncol=length(sampleID))
	colnames(other.count.mat) = sample.label
	
	for (i in 1:length(sampleID)){
	temp.file = count.files[[i]]
	temp.table = read.table(temp.file, sep="\t", header=F)
	temp.genes = as.character(temp.table[[1]])
	temp.counts = temp.table[[7]]
	
	other.count.mat[,i] = temp.counts[match(matched.genes, temp.genes, nomatch=0)]
	}#end for (i in 1:length(sampleID))
	
	genes = matched.genes
}#end if (length(matched.genes) != 0)

annotated.counts = data.frame(genes, other.count.mat)
write.table(annotated.counts, file = other.counts.file, sep="\t", row.names=F, quote=T)

result.file = paste(user.folder,other.counts.file,sep="/")
write.table(annotated.counts, file=result.file, row.names=F, quote=F, sep="\t")

#miRNA counts
count.files = as.character(sample.table$miRNA.featureCounts.files)

temp.file = count.files[[1]]
temp.table = read.table(temp.file, sep="\t", header=T)
miRNA = as.character(temp.table[[1]])

miRNA.count.mat = matrix(nrow=nrow(temp.table), ncol=length(sampleID))
colnames(miRNA.count.mat) = sample.label

matched.genes = c()

for (i in 1:length(sampleID)){
	temp.file = count.files[[i]]
	temp.table = read.table(temp.file, sep="\t", header=T)
	temp.genes = as.character(temp.table[[1]])
	
	if(length(matched.genes) != 0){
		matched.genes = temp.genes[match(matched.genes, temp.genes, nomatch=0)]
	} else if(!identical(temp.genes, miRNA)){
		print("Genes not in same order - most likely, different quantification method was used for different samples")
		userAns = readline(prompt="Do you wish to proceed with a subset of matched gene symbols? (y/n): ")
		userAns = tolower(substr(userAns, 1, 1))
		if (userAns != "y"){
			stop("Please re-run featureCounts with all of your samples")
		} else {
			matched.genes = temp.genes[match(miRNA, temp.genes, nomatch=0)]
		}#end else
	} else {
		miRNA.count.mat[,i] = temp.table[[7]]
	}#end else
}#end for (i in 1:length(sampleID))

if (length(matched.genes) != 0){
	count.mat = matrix(nrow=length(matched.genes), ncol=length(sampleID))
	colnames(count.mat) = sample.label
	
	for (i in 1:length(sampleID)){
	temp.file = count.files[[i]]
	temp.table = read.table(temp.file, sep="\t", header=F)
	temp.genes = as.character(temp.table[[1]])
	temp.counts = temp.table[[7]]
	
	miRNA.count.mat[,i] = temp.counts[match(matched.genes, temp.genes, nomatch=0)]
	}#end for (i in 1:length(sampleID))
	
	miRNA = matched.genes
}#end if (length(matched.genes) != 0)

annotated.counts = data.frame(miRNA,miRNA.count.mat)
write.table(annotated.counts, file = counts.file, sep="\t", row.names=F, quote=T)

result.file = paste(user.folder,counts.file,sep="/")
write.table(annotated.counts, file=result.file, row.names=F, quote=F, sep="\t")

#aligned reads calculation
if(aligned.type == "aligned"){
	aligned.reads = as.numeric(aligned.stat.table$Bowtie.Aligned.Reads)
} else if(aligned.type =="miRNA"){
	aligned.reads=apply(miRNA.count.mat, 2, sum)
} else if(aligned.type =="quantified-other"){
	aligned.reads=apply(other.count.mat, 2, sum)
}else {
	stop("Print CPM_norm must be either 'aligned', 'miRNA', or 'quantified-other'")
}#end else

total.million.aligned.reads = aligned.reads / 1000000

miRNA.reads = apply(miRNA.count.mat, 2, sum)
percent.miRNA.reads = round(100 * miRNA.reads / aligned.reads, digits=1)

other.reads = apply(other.count.mat, 2, sum)
percent.other.reads = round(100 * other.reads / aligned.reads, digits=1)

#other CPM
other.cpm.file = gsub(".txt$","_other.txt",cpm.file)

other.CPM = round(t(apply(other.count.mat, 1, normalizeTotalExpression, totalReads = total.million.aligned.reads)))
colnames(other.CPM) = sample.label

annotated.cpm = data.frame(genes, other.CPM)
write.table(annotated.cpm, file = other.cpm.file, sep="\t", row.names=F, quote=T)

result.file = paste(user.folder, other.cpm.file, sep="/")
write.table(annotated.cpm, file=result.file, row.names=F, quote=F, sep="\t")

#miRNA CPM
miRNA.CPM = round(t(apply(miRNA.count.mat, 1, normalizeTotalExpression, totalReads = total.million.aligned.reads)))
colnames(miRNA.CPM) = sample.label

trimmed.percent = apply(miRNA.CPM, 2, trimmed.counts, min.percent=0.3, max.percent=0.95)

expressed.gene.counts = apply(miRNA.CPM, 2, count.defined.values, expr.cutoff = min.expression)
percent.expressed.genes = round( 100 * expressed.gene.counts / nrow(miRNA.CPM), digits=1)
coverage.table = data.frame(Sample = sample.label, total.reads = total.reads, percent.cutadapt.pass=percent.cutadapt.pass,
							aligned.reads=aligned.reads, percent.aligned.reads =percent.aligned.reads,
							percent.miRNA.reads = paste(percent.miRNA.reads,"%",sep=""),
							percent.exon.reads = paste(percent.other.reads,"%",sep=""),
							robust.miRNA.count = expressed.gene.counts, robust.miRNA.percent = paste(percent.expressed.genes,"%",sep=""),
							trimmed.percent=paste(trimmed.percent,"%",sep=""))
write.table(coverage.table, file="miRNA_coverage_stats.txt", quote=F, row.names=F, sep="\t")

#tables have different file formats for downstream R analysis versus other applications that involve parsing text files
#	--> don't set the Result folder to the working directory, or you may skip genes during DEG analysis
annotated.cpm = data.frame(miRNA, miRNA.CPM)
write.table(annotated.cpm, file = cpm.file, sep="\t", row.names=F, quote=T)

result.file = paste(user.folder, cpm.file, sep="/")
write.table(annotated.cpm, file=result.file, row.names=F, quote=F, sep="\t")
