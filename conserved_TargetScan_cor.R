deg.file = "[comp.id]_[miRNA-Seq criteria].txt"
result.file = "[comp.id]_TargetScan_target_cor.txt"


meta.file = "sample_description.txt"
miRNA.file = "CPM.txt"
#code should work with either RNA-Seq or microarray data, but you may need to modify code for specific labels
mRNA.file = "log2_fpkm.txt"
target.file = "Predicted_Targets_Info.default_predictions.txt"

meta.table = read.table(meta.file, head=T, sep="\t")
miRNA.samples = as.character(meta.table$userID)

miRNA.table = read.table(miRNA.file, head=T, sep="\t")
miRNA = as.character(miRNA.table$miRNA)
miRNA = sub("hsa-","", miRNA)
miRNA.expr = log2(miRNA.table[,match(miRNA.samples,names(miRNA.table))]+1)
miRNA.samples = names(miRNA.expr)

mRNA.table = read.table(mRNA.file, head=T, sep="\t")
genes = as.character(mRNA.table$symbol)
mRNA.expr = mRNA.table[,match(miRNA.samples, names(mRNA.table))]

target.table = read.table(target.file, head=T, sep="\t")
paired.miRNA = target.table$miR.Family
paired.mRNA = target.table$Gene.Symbol

pairIDs = paste(paired.miRNA, paired.mRNA, sep="-")
print(dim(target.table))
target.table = target.table[match(unique(pairIDs),pairIDs),]
print(dim(target.table))
target.miRNAs = as.character(target.table$miR.Family)
target.miRNAs = sub("hsa-","", target.miRNAs)
target.mRNAs = as.character(target.table$Gene.Symbol)

paired.miRNA = c()
paired.mRNA = c()
miRNA.trend = c()
target.pvalues = c()
target.cor = c()

deg.table = read.table(deg.file, sep="\t", header=T)
deg.table = deg.table[deg.table$status != "No Change",]
print(dim(deg.table))

for (i in 1:nrow(deg.table)){
	print(deg.table[i,])
	temp.miRNA = as.character(deg.table$miRNA[i])
	temp.miRNA = sub("hsa-","", temp.miRNA)
	temp.trend = as.character(deg.table$status[i])
	if (temp.miRNA %in% target.miRNAs){
		temp.paired.miRNA = as.character(target.miRNAs[target.miRNAs == temp.miRNA])
		temp.paired.mRNA = as.character(target.mRNAs[target.miRNAs == temp.miRNA])
		print(temp.paired.miRNA)
		print(temp.paired.mRNA)
		paired.miRNA = c(paired.miRNA, temp.paired.miRNA)
		paired.mRNA = c(paired.mRNA, temp.paired.mRNA)
		miRNA.trend = c(miRNA.trend, rep(temp.trend, times=length(temp.paired.miRNA )))
		
		temp.cor = rep(NA, times=length(temp.paired.miRNA ))
		temp.pvalue = rep(NA, times=length(temp.paired.miRNA))
		
		for (j in 1:length(temp.paired.miRNA)){
			if((temp.paired.miRNA[j] %in% miRNA) & (temp.paired.mRNA[j] %in% genes)){
				temp.miRNA.expr = as.numeric(miRNA.expr[miRNA == temp.paired.miRNA[j],])
				temp.mRNA.expr = mRNA.expr[genes == temp.paired.mRNA[j],]
				if(nrow(temp.mRNA.expr) == 1){
					temp.mRNA.expr = as.numeric(temp.mRNA.expr)
				} else {
					stop("Shouldn't have to average genes with this summarization")
					#temp.mRNA.expr = apply(temp.mRNA.expr, 2, mean)
				}

				result = cor.test(as.numeric(temp.miRNA.expr), as.numeric(temp.mRNA.expr))
				temp.cor[j] = result$estimate
				temp.pvalue[j] = result$p.value				
			} else if(temp.paired.miRNA[j] %in% miRNA){
				multimapped.probes = genes[grep(" // ",genes)]
				potential.matches = multimapped.probes[grep(temp.paired.mRNA[j], multimapped.probes)]
				if(length(potential.matches > 1)){
					match.flag = rep(0, times=length(potential.matches))
					for (k in 1:length(potential.matches)){
						temp.gene.names = unlist(strsplit(potential.matches[k], split= " // "))
						if(temp.paired.mRNA[j] %in% temp.gene.names){
							match.flag[k] = 1
						}
					}#end for (k in 1:length(potential.matches))
					if(sum(match.flag) >= 1){
						print(paste("Using multimapped probe for ",temp.paired.mRNA[j],sep=""))
						use.probes = potential.matches[match.flag == 1]
						
						temp.miRNA.expr = as.numeric(miRNA.expr[miRNA == temp.paired.miRNA[j],])
						temp.mRNA.expr = mRNA.expr[match(use.probes, genes),]
						if(nrow(temp.mRNA.expr) == 1){
							temp.mRNA.expr = as.numeric(temp.mRNA.expr)
						} else {
							temp.mRNA.expr = apply(temp.mRNA.expr, 2, mean)
						}

						result = cor.test(as.numeric(temp.miRNA.expr), as.numeric(temp.mRNA.expr))
						temp.cor[j] = result$estimate
						temp.pvalue[j] = result$p.value							
					}#end if(sum(match.flag) >= 1)

				}#end if(length(potential.matches > 1))
			}#end else if(as.character(temp.pairs$miR.Family[j]) %in% miRNA)
		}#end for (j in 1:nrow(temp.pairs))
		
		target.pvalues = c(target.pvalues, temp.pvalue)
		target.cor = c(target.cor, temp.cor)
	}#end if (temp.miRNA %in% target.table$miR.Family)
	
}#end for (i in 1:nrow(deg.table))

target.fdr = p.adjust(target.pvalues, "fdr")

output.table = data.frame(miRNA=paired.miRNA, miRNA.trend=miRNA.trend, target=paired.mRNA,
						cor.coef = target.cor, cor.lm.pvalue = target.pvalues, cor.lm.fdr = target.fdr)
output.table = output.table[!is.na(output.table$cor.lm.pvalue),]
write.table(output.table, file="TargetScan_Conserved_cor.txt", sep="\t", row.names=F)

summary.miRNA = as.character(levels(as.factor(as.character(output.table$miRNA))))
num.neg.cor = rep(0, times=length(summary.miRNA))
num.targets = rep(0, times=length(summary.miRNA))
num.pos.cor = rep(0, times=length(summary.miRNA))
pos.cor.rate = rep(NA, times=length(summary.miRNA))

for (i in 1:length(summary.miRNA)){
	print(summary.miRNA[i])
	temp.table = output.table[output.table$miRNA == summary.miRNA[i],]
	print(dim(temp.table))
	num.targets[i] = nrow(temp.table)
	num.neg.cor[i] = length(temp.table$miRNA[(temp.table$cor.coef < 0) & (temp.table$cor.lm.pvalue < 0.05)])
	num.pos.cor[i] = length(temp.table$miRNA[(temp.table$cor.coef > 0) & (temp.table$cor.lm.pvalue < 0.05)])
	if (num.pos.cor[i] >= num.neg.cor[i]){
		pos.cor.rate[i] = 1
	} else{
		pos.cor.rate[i] = num.pos.cor[i] / num.neg.cor[i]
	}
}#end for (i in 1:length(summary.miRNA))

summary.table = data.frame(miRNA= summary.miRNA, neg.cor = num.neg.cor, total.cor = num.targets,
							pos.cor = num.pos.cor, miRNA.estimated.fdr = pos.cor.rate)
write.table(summary.table, file=result.file, sep="\t", row.names=F)
