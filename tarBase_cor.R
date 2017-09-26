
tarbase.target.mapping.file = "selected_tarBase_pairs.txt"
cor.file = "../../../../Result/Target_Integration/tarBase/combined_tarbase_target_cor.txt"
result.file = "../../../../Result/Target_Integration/tarBase/combined_tarbase_target_summary.txt"

meta.file = "sample_description.txt"
miRNA.file = "CPM.txt"
#code should work with RNA-Seq or microarray data
mRNA.file = "log2_fpkm.txt"


meta.table = read.table(meta.file, head=T, sep="\t")
sample.labels = as.character(meta.table$userID)

miRNA.table = read.table(miRNA.file, head=T, sep="\t")
miRNA = as.character(miRNA.table$miRNA)
miRNA.expr = log2(miRNA.table[,match(sample.labels,names(miRNA.table))]+1)
miRNA.samples = names(miRNA.expr)

mRNA.table = read.table(mRNA.file, head=T, sep="\t")
genes = as.character(mRNA.table$symbol)
mRNA.expr = mRNA.table[,match(sample.labels,names(mRNA.table))]

target.table = read.table(tarbase.target.mapping.file, head=T, sep="\t")
interactions = paste(target.table[,1],target.table[,2])
print(dim(target.table))
target.table = target.table[match(unique(interactions),interactions),]
print(dim(target.table))

target.pvalues = rep(NA, times=nrow(target.table))
target.cor = rep(NA, times=nrow(target.table))

for (i in 1:nrow(target.table)){
	print(target.table[i,])
	temp.miRNA = as.character(target.table[i,1])
	temp.mRNA = as.character(target.table[i,2])
	
	if((temp.miRNA %in% miRNA) & (temp.mRNA %in% genes)){
		temp.miRNA.expr = as.numeric(miRNA.expr[miRNA == temp.miRNA,])
		temp.mRNA.expr = mRNA.expr[genes == temp.mRNA,]
		if(nrow(temp.mRNA.expr) == 1){
			temp.mRNA.expr = as.numeric(temp.mRNA.expr)
		} else {
			stop("Shouldn't have to average genes with this summarization")
			#temp.mRNA.expr = apply(temp.mRNA.expr, 2, mean)
		}

		result = cor.test(as.numeric(temp.miRNA.expr), as.numeric(temp.mRNA.expr))
		target.cor[i] = result$estimate
		target.pvalues[i] = result$p.value				
	} else if(temp.miRNA %in% miRNA){
		multimapped.probes = genes[grep(" // ",genes)]
		potential.matches = multimapped.probes[grep(temp.mRNA, multimapped.probes)]
		if(length(potential.matches > 1)){
			match.flag = rep(0, times=length(potential.matches))
			for (k in 1:length(potential.matches)){
				temp.gene.names = unlist(strsplit(potential.matches[k], split= " // "))
				if(temp.mRNA %in% temp.gene.names){
					match.flag[k] = 1
				}
			}#end for (k in 1:length(potential.matches))
			if(sum(match.flag) >= 1){
				print(paste("Using multimapped probe for ",temp.mRNA,sep=""))
				use.probes = potential.matches[match.flag == 1]
						
				temp.miRNA.expr = as.numeric(miRNA.expr[miRNA == temp.miRNA,])
				temp.mRNA.expr = mRNA.expr[match(use.probes, genes),]
				if(nrow(temp.mRNA.expr) == 1){
					temp.mRNA.expr = as.numeric(temp.mRNA.expr)
				} else {
					stop("Shouldn't have to average genes with this summarization")
					#temp.mRNA.expr = apply(temp.mRNA.expr, 2, mean)
				}

				result = cor.test(as.numeric(temp.miRNA.expr), as.numeric(temp.mRNA.expr))
				target.cor[i] = result$estimate
				target.pvalues[i] = result$p.value						
			}#end if(sum(match.flag) >= 1)

		}#end if(length(potential.matches > 1))
	}#end else if(as.character(temp.pairs$miRNA[j]) %in% miRNA)
}#end for (i in 1:nrow(target.table))

target.fdr = p.adjust(target.pvalues, "fdr")


output.table = data.frame(target.table, cor.coef = target.cor, 	cor.pvalue = target.pvalues, cor.fdr = target.fdr)
output.table = output.table[!is.na(output.table$cor.pvalue),]
write.table(output.table, file=cor.file, sep="\t", row.names=F)

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
	num.neg.cor[i] = length(temp.table$miRNA[(temp.table$cor.coef < 0) & (temp.table$cor.pvalue < 0.05)])
	num.pos.cor[i] = length(temp.table$miRNA[(temp.table$cor.coef > 0) & (temp.table$cor.pvalue < 0.05)])
	if (num.pos.cor[i] >= num.neg.cor[i]){
		pos.cor.rate[i] = 1
	} else{
		pos.cor.rate[i] = num.pos.cor[i] / num.neg.cor[i]
	}
}#end for (i in 1:length(summary.miRNA))

summary.table = data.frame(miRNA= summary.miRNA, neg.cor = num.neg.cor, total.cor = num.targets,
							pos.cor = num.pos.cor, miRNA.estimated.fdr = pos.cor.rate)
write.table(summary.table, file=result.file, sep="\t", row.names=F)
