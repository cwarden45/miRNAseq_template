normalizeTotalExpression = function (geneExpr, totalReads) {
	return(geneExpr / totalReads)
}#end def normalizeTotalExpression

avgGroupExpression = function (geneExpr, groups) {
	avg.expr = tapply(geneExpr, groups, mean)
	return(avg.expr)
}#end def avgGroupExpression

standardize.arr = function(arr)
{
	center.arr = as.numeric(arr) - mean(as.numeric(arr), na.rm=T)
	norm.arr = center.arr / sd(center.arr, na.rm=T)
	return(norm.arr)
}#end def standardize.arr

count.defined.values = function(arr, expr.cutoff)
{
	return(length(arr[arr > expr.cutoff]))
}#end def count.values

count.na.values = function(arr)
{
	return(length(arr[is.na(arr)]))
}#end def count.values

ratio2fc = function(value)
{
	if(value >= 0){
		return(2^value)
	} else {
		return (-2^(-value))
	}
}#end def ratio2fc

gene.lm = function(arr, var1, var2=c(), var3=c())
{	
	if (length(var2) == 0){
		fit = lm(as.numeric(arr) ~ var1)
		result = summary(fit)
		pvalue = result$coefficients[2,4]
	} else if (length(var3) == 0){
		fit = lm(as.numeric(arr) ~ var1 + var2)
		result = summary(fit)
		pvalue = result$coefficients[2,4]
	} else {
		fit = lm(as.numeric(arr) ~ var2*var3 + var2 + var3)
		result = summary(fit)
		pvalue = result$coefficients[4,4]
	}
	return(pvalue)
}#end def gene.lm

gene.aov = function(arr, var1, var2=c(), var3=c())
{	
	if (length(var2) == 0){
		fit = aov(as.numeric(arr) ~ var1)
		result = summary(fit)
		aov.pvalue = result[[1]][['Pr(>F)']][1]
	} else if (length(var3) == 0){
		fit = aov(as.numeric(arr) ~ var1 + var2)
		result = summary(fit)
		aov.pvalue = result[[1]][['Pr(>F)']][1]
	} else {
		fit = aov(as.numeric(arr) ~ var2*var3 + var2 + var3)
		result = summary(fit)
		aov.pvalue = result[[1]][['Pr(>F)']][3]
	}
	return(aov.pvalue)
}#end def gene.aov

calc.gene.cor = function(arr, indep.var)
{	
	na.count = length(arr[!is.na(arr)])
	if((na.count >= 3) & (sd(arr) != 0)){
		gene.cor.coef = cor(arr,indep.var)
	} else {
		gene.cor.coef = NA
	}
	return(gene.cor.coef)
}#end def calc.gene.cor

library(gplots)
fixed.color.palatte = c("green","orange","purple","cyan","pink","maroon","yellow","grey","red","blue","black","darkgreen","thistle1","tan","orchid1",colors())

param.table = read.table("parameters.txt", header=T, sep="\t")
comp.name=as.character(param.table$Value[param.table$Parameter == "comp_name"])
min.expression = as.numeric(as.character(param.table$Value[param.table$Parameter == "cpm_expression_cutoff"]))
min.fraction.expressed = as.numeric(as.character(param.table$Value[param.table$Parameter == "minimum_fraction_expressed"]))
fc.cutoff = as.numeric(as.character(param.table$Value[param.table$Parameter == "fold_change_cutoff"]))
cor.cutoff = as.numeric(as.character(param.table$Value[param.table$Parameter == "cor_cutoff"]))
pvalue.cutoff = as.numeric(as.character(param.table$Value[param.table$Parameter == "pvalue_cutoff"]))
fdr.cutoff = as.numeric(as.character(param.table$Value[param.table$Parameter == "fdr_cutoff"]))
deg.groups = unlist(strsplit(as.character(param.table$Value[param.table$Parameter == "deg_groups"]), split=","))
plot.groups = unlist(strsplit(as.character(param.table$Value[param.table$Parameter == "plot_groups"]), split=","))
trt.group = as.character(param.table$Value[param.table$Parameter == "treatment_group"])
interaction.flag = as.character(param.table$Value[param.table$Parameter == "interaction"])
interaction.flag[interaction.flag == "none"]="no"
pvalue.method = as.character(param.table$Value[param.table$Parameter == "pvalue_method"])
fdr.method = as.character(param.table$Value[param.table$Parameter == "fdr_method"])
output.folder = as.character(param.table$Value[param.table$Parameter == "Raw_Code_PC"])
user.folder = as.character(param.table$Value[param.table$Parameter == "Result_Folder"])
sample.description.file = as.character(param.table$Value[param.table$Parameter == "sample_description_file"])
cpm.file = as.character(param.table$Value[param.table$Parameter == "cpm_file"])
trt.group2 = as.character(param.table$Value[param.table$Parameter == "secondary_plot"])

setwd(output.folder)

sample.description.table = read.table(sample.description.file, sep="\t", header=T)
longID = sample.description.table$sampleID
sample.label = sample.description.table$userID

deg.group.table = sample.description.table[,deg.groups]
if (length(deg.groups) == 1){
	deg.meta = sample.description.table[!is.na(deg.group.table),]
} else {
	deg.grp.na.counts = apply(deg.group.table, 1, count.na.values)
	deg.meta = sample.description.table[deg.grp.na.counts == 0,]
}

cpm.table = read.table(cpm.file, head=T, sep="\t")
CPM= cpm.table[,match(sample.label,names(cpm.table))]
CPM = matrix(as.numeric(unlist(CPM)), ncol=length(sample.label))

CPM =CPM[rowSums(!is.na(CPM)) == length(sample.label),]
cpm.table = cpm.table[rowSums(!is.na(CPM)) == length(sample.label),]
print(dim(cpm.table))
print(dim(CPM))

print(dim(CPM))
expressed.sample.count = apply(CPM, 1, count.defined.values, expr.cutoff = min.expression)
cpm.table = cpm.table[expressed.sample.count >= round(min.fraction.expressed * ncol(CPM)),]
CPM = CPM[expressed.sample.count >= round(min.fraction.expressed * ncol(CPM)),]
print(dim(CPM))

miRNA = cpm.table$miRNA
colnames(CPM) = as.character(sample.label)
rownames(CPM) = as.character(miRNA)

log2.CPM = log2(CPM + 1)

if(length(plot.groups) == 1){
	print("Averaging Expression for One Variable (for plot.groups)")
	grp = sample.description.table[,plot.groups]
} else if ((length(plot.groups) == 2)&(interaction.flag == "no")){
	print("Averaging Expression for First Variable (for plot.groups)")
	grp = sample.description.table[,plot.groups[1]]
} else if (length(plot.groups) == 2){
	print("Averaging Expression for Interaction Variable (for plot.groups)")
	grp = paste(sample.description.table[,plot.groups[1]],sample.description.table[,plot.groups[2]],sep=":")
} else {
	stop("Code only compatible with 2 variables (with or without a 3rd interaction variable")
}

groupIDs = as.character(levels(as.factor(grp)))
average.cpm = data.frame(t(apply(log2.CPM, 1, avgGroupExpression, groups = grp)))
if(length(groupIDs) == 1){
	average.cpm = t(average.cpm)
} else {
	average.cpm = average.cpm
}
colnames(average.cpm) = paste("avg.log2.cpm", sub("-",".",groupIDs), sep=".")


#remove undefined group IDs (so, you can visualize samples that aren't in your comparison)
if(length(deg.groups) == 1){
	var1 = sample.description.table[,deg.groups]
	deg.log2.CPM = log2.CPM[,!is.na(var1)]
	deg.CPM = CPM[,!is.na(var1)]
	var1 = var1[!is.na(var1)]
	if (trt.group != "continuous"){
		var1 = as.factor(as.character(var1[!is.na(var1)]))
	}
} else if (length(deg.groups) == 2){
		var1 = sample.description.table[,deg.groups[1]]
		var2 = sample.description.table[,deg.groups[2]]
		deg.samples = !is.na(var1)&!is.na(var2)
		deg.log2.CPM = log2.CPM[,deg.samples]
		deg.CPM = CPM[,deg.samples]
		var1 = var1[deg.samples]
		if (trt.group != "continuous"){
			var1 = as.factor(as.character(var1[!is.na(var1)]))
		}
		var2 = var2[deg.samples]
		if (trt.group2 != "continuous"){
			var2 = as.factor(as.character(var2[!is.na(var2)]))
		}
} else {
	stop("Code currently doesn't support more than 2 group model for DEG (with or without interaction)")
}

if(length(deg.groups) == 1){
	print("Averaging Expression for One Variable (for deg.groups)")
	contrast.grp = var1
} else if ((length(deg.groups) == 2)&(interaction.flag == "no")){
	print("Averaging Expression for First Variable (for deg.groups)")
	contrast.grp = var1
} else if (length(deg.groups) == 2){
	print("Averaging Expression for Interaction Variable (for deg.groups)")
	contrast.grp = paste(var1,var2,sep=":")
} else {
	stop("Code only compatible with 2 variables (with or without a 3rd interaction variable")
}

if (trt.group == "continuous"){
	contrast.grp = as.numeric(contrast.grp)
	
	gene.cor = apply(deg.log2.CPM, 1, calc.gene.cor, indep.var=contrast.grp)
	
	fc.table = data.frame(cor=gene.cor)
} else {
	groupIDs = as.character(levels(as.factor(contrast.grp)))
	contrast.cpm = data.frame(t(apply(deg.log2.CPM, 1, avgGroupExpression, groups = contrast.grp)))
	colnames(contrast.cpm) = paste("avg.log2.cpm", sub("-",".",groupIDs), sep=".")
}#end else

if((interaction.flag == "no") & (trt.group != "continuous")){
	print("Calculating fold-change for primary variable")
	trt.expr = contrast.cpm[,paste("avg.log2.cpm", sub("-",".",trt.group), sep=".")]
	cntl.expr = contrast.cpm[,paste("avg.log2.cpm", sub("-",".",groupIDs[groupIDs != trt.group]), sep=".")]

	log2ratio = round(trt.expr - cntl.expr, digits = 2)
	fc = round(sapply(log2ratio, ratio2fc), digits = 2)
	fc.table = data.frame(log2ratio=log2ratio, fold.change=fc)
} else if (interaction.flag == "model"){
	if ((trt.group == "continuous")&(trt.group2 == "continuous")){
		print("Calculating  correlation for secondary variable")
		sec.contrast.grp = as.numeric(var2)
		
		gene.cor2 = apply(deg.log2.CPM, 1, calc.gene.cor, indep.var=sec.contrast.grp)
		
		fc.table = data.frame(prim.cor=gene.cor, sec.cor = gene.cor2)
	} else if (trt.group == "continuous"){
		print("Fold-change / correlation cutoff not used for mixed variable analysis")
		print("NOTE: 'Up-Regulated' R output refers to genes that vary with FDR and p-value cutoffs")
		print("However, fold-change / correlation values for each separate variable are still provided")

		sec.groupIDs = var2
		sec.groups = as.character(levels(as.factor(sec.groupIDs)))
		sec.contrast.cpm = data.frame(t(apply(deg.log2.CPM, 1, avgGroupExpression, groups = sec.groupIDs)))
		colnames(sec.contrast.cpm) = paste("avg.log2.cpm", sub("-",".",groupIDs), sep=".")
		sec.trt.expr = sec.contrast.cpm[,paste("avg.log2.cpm", sub("-",".",trt.group2), sep=".")]
		sec.cntl.expr = sec.contrast.cpm[,paste("avg.log2.cpm", sub("-",".",sec.groups[sec.groups != trt.group2]), sep=".")]

		sec.log2ratio = round(sec.trt.expr - sec.cntl.expr, digits = 2)
		sec.fc = round(sapply(sec.log2ratio, ratio2fc), digits = 2)
		
		fc.table = data.frame(prim.cor=gene.cor, sec.fc = sec.fc)
	} else if (trt.group2 == "continuous"){	
		print("Fold-change / correlation cutoff not used for mixed variable analysis")
		print("NOTE: 'Up-Regulated' R output refers to genes that vary with FDR and p-value cutoffs")
		print("However, fold-change / correlation values for each separate variable are still provided")

		prim.groupIDs = var1
		prim.groups = as.character(levels(as.factor(prim.groupIDs)))
		prim.contrast.cpm = data.frame(t(apply(deg.log2.CPM, 1, avgGroupExpression, groups = prim.groupIDs)))
		colnames(prim.contrast.cpm) = paste("avg.log2.cpm", sub("-",".",prim.groups), sep=".")
		prim.trt = trt.group
		prim.cntl = prim.groups[prim.groups != trt.group]
		prim.trt.expr = prim.contrast.cpm[,paste("avg.log2.cpm", sub("-",".",prim.trt), sep=".")]
		prim.cntl.expr = prim.contrast.cpm[,paste("avg.log2.cpm", sub("-",".",prim.cntl), sep=".")]

		prim.log2ratio = round(prim.trt.expr - prim.cntl.expr, digits = 2)
		prim.fc = round(sapply(prim.log2ratio, ratio2fc), digits = 2)
		
		sec.contrast.grp = as.numeric(var2)
		gene.cor2 = apply(deg.log2.CPM, 1, calc.gene.cor, indep.var=sec.contrast.grp)
		
		fc.table = data.frame(prim.fc=prim.fc, sec.cor = gene.cor2)
	} else {
		print("Calculating fold-change table for primary variables (within subsets of secondary variable)")
		prim.groups = paste(var1,var2,sep=":")
		prim.trt = paste(trt.group,trt.group2,sep=":")
		prim.cntl = paste(prim.groups[prim.groups != trt.group],trt.group2,sep=":")
		prim.trt.expr = contrast.cpm[,paste("avg.log2.cpm", sub("-",".",prim.trt), sep=".")]
		prim.cntl.expr = contrast.cpm[,paste("avg.log2.cpm", sub("-",".",prim.cntl), sep=".")]

		prim.log2ratio = round(prim.trt.expr - prim.cntl.expr, digits = 2)
		prim.fc = round(sapply(prim.log2ratio, ratio2fc), digits = 2)

		sec.groups = as.character(levels(as.factor(sample.description.table[,deg.groups[2]])))
		sec.trt = paste(trt.group, sec.groups[sec.groups != trt.group2], sep=":")
		sec.cntl = paste(prim.groups[prim.groups != trt.group], sec.groups[sec.groups != trt.group2], sep=":")
		sec.trt.expr = contrast.cpm[,paste("avg.log2.cpm", sub("-",".",sec.trt), sep=".")]
		sec.cntl.expr = contrast.cpm[,paste("avg.log2.cpm", sub("-",".",sec.cntl), sep=".")]

		sec.log2ratio = round(sec.trt.expr - sec.cntl.expr, digits = 2)
		sec.fc = round(sapply(sec.log2ratio, ratio2fc), digits = 2)

		overall.log2ratio = prim.log2ratio - sec.log2ratio
		overall.fc = round(sapply(overall.log2ratio, ratio2fc), digits = 2)

		fc.table = data.frame(fc1 = prim.fc, fc2=sec.fc, fc3=overall.fc)
		colnames(fc.table) = c(paste("fold.change",trt.group,":",trt.group2,sep="."),
								paste("fold.change",trt.group,":",sec.groups[sec.groups != trt.group2], sep="."),
								"overall.fold.change")
	}#end else
}else if(trt.group == "continuous"){
	print("Skipping fold-change calculation for continuous variable")
}else{
		stop("interaction must be \"no\", \"model\", or \"filter-overlap\"")
}#end else

rep.check = 1
for (i in 1:length(deg.groups)){
	deg.group = deg.groups[i]
	
	if((i == 1) & (trt.group != "continuous")){
		deg.group.values = as.factor(as.character(deg.meta[,deg.group]))
		min.reps = min(table(deg.group.values))
		if (min.reps < 2){
			rep.check=0
			print("There are not at least 2 samples per-group in order to calculate p-value.")
			print("In the future, please make sure you at least have duplicate samples.")
		}#end if (min.reps < 2)
	} else if ((i == 2) & (trt.group2 != "continuous")){
		deg.group.values = as.factor(as.character(deg.meta[,deg.group]))
		min.reps = min(table(deg.group.values))
		if (min.reps < 2){
			rep.check=0
			print("There are not at least 2 samples per-group in order to calculate p-value.")
			print("In the future, please make sure you at least have duplicate samples.")
		}#end if (min.reps < 2)
	} else if (i > 2){
		stop("Workflow currently doesn't support use of more than 2 variables")
	}
}#end for (deg.group in deg.groups)

if(rep.check == 1){
	#start p-value calculation
	if (pvalue.method == "edgeR"){
		library(edgeR)
		
		#can test removing ', lib.size = rep(1000000,ncol(deg.CPM))'
		y = DGEList(counts=deg.CPM, genes=miRNA, lib.size = rep(1000000,ncol(deg.CPM)))
		y = estimateCommonDisp(y)
		if (length(deg.groups) == 1){
			print("edgeR with 1 variable")
			if (trt.group == "continuous"){
				var1 = as.numeric(var1)
			}
			design = model.matrix(~var1)
			fit = glmFit(y, design)
			lrt = glmLRT(fit, coef=2)
			test.pvalue = lrt$table$PValue
		} else if ((length(deg.groups) == 2)&(interaction.flag == "no")){
			print("edgeR with 2 variables")

			if (trt.group == "continuous"){
				var1 = as.numeric(var1)
			} else{
				var1 = as.factor(var1)
			}

			if (trt.group2 == "continuous"){
				var2 = as.numeric(var2)
			} else{
				var2 = as.factor(var2)
			}
			design = model.matrix(~var1 + var2)
			fit = glmFit(y, design)
			lrt = glmLRT(fit, coef=2)
			test.pvalue = lrt$table$PValue
		} else if ((length(deg.groups) == 2)&(interaction.flag == "model")){
			print("edgeR with 2 variables plus interaction")

			if (trt.group == "continuous"){
				var1 = as.numeric(var1)
			}

			if (trt.group2 == "continuous"){
				var2 = as.numeric(var2)
			}
			design = model.matrix(~var1*var2 + var1 + var2)
			fit = glmFit(y, design)
			lrt = glmLRT(fit, coef=4)
			test.pvalue = lrt$table$PValue
		}
	} else if (pvalue.method == "DESeq2"){
		library(DESeq2)
		
		if (length(deg.groups) == 1){
			print("DESeq2 with 1 variable")
			if (trt.group == "continuous"){
				var1 = as.numeric(var1)
			}

			colData = data.frame(var1=var1)
			rownames(colData) = colnames(deg.CPM)
			dds= DESeqDataSetFromMatrix(countData = deg.CPM,
							colData = colData,
							design = ~ var1)
			dds = DESeq(dds)
			res = results(dds)
			test.pvalue = res$pvalue
		} else if ((length(deg.groups) == 2)&(interaction.flag == "no")){
			print("DESeq2 (LRT) with 2 variables")

			if (trt.group == "continuous"){
				var1 = as.numeric(var1)
			} else{
				var1 = as.factor(var1)
			}

			if (trt.group2 == "continuous"){
				var2 = as.numeric(var2)
			} else{
				var2 = as.factor(var2)
			}

			colData = data.frame(var1=var1, var2=var2)
			rownames(colData) = colnames(deg.CPM)
			dds = DESeqDataSetFromMatrix(countData = deg.CPM,
							colData = colData,
							design = ~ var1 + var2)
			Wald.flag = TRUE
			if (Wald.flag){
				dds = DESeq(dds)
				if (trt.group == "continuous"){
					res = results(dds, name = "var1")
				} else{
					other.groups = as.character(levels(as.factor(as.character(var1[var1 != trt.group]))))
					if(length(other.groups) > 1){
						print("DESeq2 Wald-test will look at differences between two groups.")
						print("You can manually switch code to use LRT instead of Wald test (set Wald.flag = FALSE).")
						#The paired sample design in the DESeq2 manual uses the Wald test
						stop("Or, please consider using an interaction model with LRT if your primary variable has more than two groups.")
					}
					res = results(dds, contrast = c("var1", trt.group, other.groups))
				}
			} else {
				dds = DESeq(dds, test="LRT", reduced = ~ var2)
				res = results(dds)
				#print(head(res))
			}
			test.pvalue = res$pvalue
		} else if ((length(deg.groups) == 2)&(interaction.flag == "model")){
			print("DESeq2 (LRT) with 2 variables plus interaction")

			if (trt.group == "continuous"){
				var1 = as.numeric(var1)
			}

			if (trt.group2 == "continuous"){
				var2 = as.numeric(var2)
			}

			colData = data.frame(var1=var1, var2=var2)
			rownames(colData) = colnames(deg.CPM)
			dds = DESeqDataSetFromMatrix(countData = deg.CPM,
							colData = colData,
							design = ~ var1*var2 + var1 + var2)
			dds = DESeq(dds, test="LRT", reduced = ~ var1 + var2)
			res = results(dds)
			test.pvalue = res$pvalue
		}
	} else if (pvalue.method == "limma-voom"){
		library(edgeR)
		library(limma)
		
		#can add ', lib.size = rep(1000000,ncol(deg.CPM))' to account for use of CPM values
		y = DGEList(counts=deg.CPM, genes=miRNA)
		if (length(deg.groups) == 1){
			print("limma-voom with 1 variable")

			if (trt.group == "continuous"){
				var1 = as.numeric(var1)
			}
			design = model.matrix(~var1)
			png(paste(comp.name,"voom_plot.png",sep="_"))
			v = voom(y,design,plot=TRUE)
			dev.off()
			fit = lmFit(v,design)
			fit = eBayes(fit)
			pvalue.mat = data.frame(fit$p.value)
			test.pvalue = pvalue.mat[,2]
		} else if ((length(deg.groups) == 2)&(interaction.flag == "no")){
			print("limma-voom with 2 variables")

			if (trt.group == "continuous"){
				var1 = as.numeric(var1)
			} else{
				var1 = as.factor(var1)
			}

			if (trt.group2 == "continuous"){
				var2 = as.numeric(var2)
			} else{
				var2 = as.factor(var2)
			}
			design = model.matrix(~var1 + var2)
			png(paste(comp.name,"voom_plot.png",sep="_"))
			v = voom(y,design,plot=TRUE)
			dev.off()
			fit = lmFit(v,design)
			fit = eBayes(fit)
			pvalue.mat = data.frame(fit$p.value)
			test.pvalue = pvalue.mat[,2]
		} else if ((length(deg.groups) == 2)&(interaction.flag == "model")){
			print("limma-voom with 2 variables plus interaction")

			if (trt.group == "continuous"){
				var1 = as.numeric(var1)
			}

			if (trt.group2 == "continuous"){
				var2 = as.numeric(var2)
			}
			design = model.matrix(~var1*var2 + var1 + var2)
			png(paste(comp.name,"voom_plot.png",sep="_"))
			v = voom(y,design,plot=TRUE)
			dev.off()
			fit = lmFit(v,design)
			fit = eBayes(fit)
			pvalue.mat = data.frame(fit$p.value)
			test.pvalue = pvalue.mat[,4]
		}
	} else if (pvalue.method == "lm"){
		if (length(deg.groups) == 1){
			print("log2(CPM+1) linear regression with 1 variable")

				if (trt.group == "continuous"){
					var1 = as.numeric(var1)
				}
				test.pvalue = apply(deg.log2.CPM, 1, gene.lm, var1=var1)
		} else if ((length(deg.groups) == 2)&(interaction.flag == "no")){
			print("log2(CPM+1) linear regression with 2 variables")

			if (trt.group == "continuous"){
				var1 = as.numeric(var1)
			}
				
			if (trt.group2 == "continuous"){
				var2 = as.numeric(var2)
			}
			test.pvalue = apply(deg.log2.CPM, 1, gene.lm, var1=var1, var2=var2)
		} else if ((length(deg.groups) == 2)&(interaction.flag == "model")){
			print("log2(CPM+1) linear regression with 2 variables plus interaction")
			var3 = as.factor(paste(var1,var2,sep=":"))

			if (trt.group == "continuous"){
				var1 = as.numeric(var1)
			}

			if (trt.group == "continuous"){
				var2 = as.numeric(var2)
			}
			test.pvalue = apply(deg.log2.CPM, 1, gene.lm, var1=var3, var2=var1, var3=var2)
		}
	} else if (pvalue.method == "ANOVA"){
		if (length(deg.groups) == 1){
			print("log2(CPM+1) ANOVA with 1 variable")
				
			if (trt.group == "continuous"){
				var1 = as.numeric(var1)
			}
			test.pvalue = apply(deg.log2.CPM, 1, gene.aov, var1=var1)
		} else if ((length(deg.groups) == 2)&(interaction.flag == "no")){
			print("log2(CPM+1) ANOVA with 2 variables")

			if (trt.group == "continuous"){
				var1 = as.numeric(var1)
			}

			if (trt.group2 == "continuous"){
				var2 = as.numeric(var2)
			}
			test.pvalue = apply(deg.log2.CPM, 1, gene.aov, var1=var1, var2=var2)
		} else if ((length(deg.groups) == 2)&(interaction.flag == "model")){
			print("log2(CPM+1) ANOVA with 2 variables plus interaction")
			var3 = as.factor(paste(var1,var2,sep=":"))

			if (trt.group == "continuous"){
				var1 = as.numeric(var1)
			}

			if (trt.group2 == "continuous"){
				var2 = as.numeric(var2)
			}
			test.pvalue = apply(deg.log2.CPM, 1, gene.aov, var1=var3, var2=var1, var3=var2)
		}
	} else{
		stop("pvalue_method must be \"edgeR\", \"limma-voom\", \"DESeq2\", \"lm\", or \"ANOVA\"")
	}
} else{
	test.pvalue = rep(1,times=length(miRNA))
	prim.pvalue = rep(1,times=length(miRNA))
	sec.pvalue = rep(1,times=length(miRNA))
}#end else


if (trt.group == "continuous"){
	upID = "Increased Expression"
	downID = "Decreased Expression"
} else {
	upID = paste(trt.group," Up",sep="")
	downID = paste(trt.group," Down",sep="")	
}


if (interaction.flag == "no"){
	if (fdr.method == "BH"){
		fdr = p.adjust(test.pvalue, "fdr")
	} else if (fdr.method == "q-value"){
		library(qvalue)
		qobj <- qvalue(p = test.pvalue)
		fdr = qobj$qvalue
		png(paste(comp.name,"_",pvalue.method,"_qvalue_plot.png",sep=""))
		qHist = hist(qobj)
		print(qHist)
		dev.off()
	} else if (fdr.method == "q-lfdr"){
		library(qvalue)
		qobj <- qvalue(p = test.pvalue)
		fdr = qobj$lfdr
		png(paste(comp.name,"_",pvalue.method,"_qvalue_plot.png",sep=""))
		qHist = hist(qobj)
		print(qHist)
		dev.off()
	} else {
		stop("fdr_method must be \"BH\", \"q-value\", or \"q-lfdr\"")
	}
	status = rep("No Change", times=length(fdr))
	if (trt.group == "continuous"){
		status[(gene.cor >= cor.cutoff) & (test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = upID
		status[(gene.cor <= -cor.cutoff) & (test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = downID
		pvalue.table = data.frame(p.value = test.pvalue, FDR = fdr)
	} else{
		status[(fc >= fc.cutoff) & (test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = upID
		status[(fc <= -fc.cutoff) & (test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = downID
		pvalue.table = data.frame(p.value = test.pvalue, FDR = fdr)
	}#end else
} else{
	trt.group = prim.trt
	if(interaction.flag == "model"){
		if (fdr.method == "BH"){
			fdr = p.adjust(test.pvalue, "fdr")
		} else if (fdr.method == "q-value"){
			library(qvalue)
			qobj <- qvalue(p = test.pvalue)
			fdr = qobj$qvalue
			png(paste(comp.name,"_",pvalue.method,"_qvalue_plot.png",sep=""))
			qHist = hist(qobj)
			print(qHist)
			dev.off()
		} else if (fdr.method == "q-lfdr"){
			library(qvalue)
			qobj <- qvalue(p = test.pvalue)
			fdr = qobj$lfdr
			png(paste(comp.name,"_",pvalue.method,"_qvalue_plot.png",sep=""))
			qHist = hist(qobj)
			print(qHist)
			dev.off()
		} else {
			stop("fdr_method must be \"BH\", \"q-value\", or \"q-lfdr\"")
		}
		status = rep("No Change", times=length(fdr))
		if ((trt.group == "continuous")&(trt.group2 == "continuous")){
			status[(gene.cor.int >= cor.cutoff) & (test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = upID
			status[(gene.cor.int <= -cor.cutoff) & (test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = downID
			pvalue.table = data.frame(p.value = test.pvalue, FDR = fdr)
		} else if ((trt.group != "continuous")&(trt.group2 != "continuous")){
			status[(overall.fc >= fc.cutoff) & (test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = upID
			status[(overall.fc <= -fc.cutoff) & (test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = downID
			pvalue.table = data.frame(p.value = test.pvalue, FDR = fdr)
		} else {
			upID = "Variable Expression"
			status[(test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = upID
			pvalue.table = data.frame(p.value = test.pvalue, FDR = fdr)
		}#end else
	} else{
		stop("interaction must be \"no\" or \"model\"")
	}#end else
}#end else

print(paste("Up-Regulated: ",length(status[status == upID]),sep=""))
print(paste("Down-Regulated: ",length(status[status == downID]),sep=""))

if(rep.check == 1){
	deg.table = data.frame(miRNA = miRNA,
							average.cpm, fc.table,
							pvalue.table, status = status)
} else {
	deg.table = data.frame(miRNA = miRNA,
							average.cpm, fc.table, status = status)	
}#end else

deg.file = paste(comp.name,"_",pvalue.method,"_DEG_fc_",fc.cutoff,"_fdr_",fdr.cutoff,"_pval_",pvalue.cutoff,".txt",sep="")
deg.file = gsub(":",".",deg.file)
write.table(deg.table, file=deg.file, row.names=F, quote=F, sep="\t")

final.deg.file = paste(user.folder,"/DEG/",comp.name,"_DEG_stats.txt",sep="")
write.table(deg.table, file=final.deg.file, row.names=F, quote=F, sep="\t")

temp.cpm = log2.CPM
temp.cpm = temp.cpm[status != "No Change", ]
deg.miRNA = miRNA[status != "No Change"]

if(length(deg.miRNA) > 1){
	if(length(plot.groups) > 1){
		source("heatmap.3.R")
		if((trt.group != "continuous")&(trt.group2 != "continuous")){
			grp1 = as.character(sample.description.table[,plot.groups[1]])
			grp2 = as.character(sample.description.table[,plot.groups[2]])

			group.levels = c(levels(as.factor(grp1)),levels(as.factor(grp2)))

			color.palette <- fixed.color.palatte[1:length(group.levels)]
			labelColors1 = rep("black",times=length(sample.label))
			for (i in 1:length(group.levels)){
				labelColors1[grp1 == as.character(group.levels[i])] = color.palette[i]
			}#end for (i in 1:length(group.levels))
			labelColors2 = rep("black",times=length(sample.label))
			for (i in 1:length(group.levels)){
				labelColors2[grp2 == as.character(group.levels[i])] = color.palette[i]
			}#end for (i in 1:length(group.levels))
		}else if((trt.group == "continuous")&(trt.group2 == "continuous")){
			stop("Add code for two continuous variables")
		}else if(trt.group == "continuous"){
			grp1 = as.numeric(sample.description.table[,plot.groups[1]])
			grp2 = as.character(sample.description.table[,plot.groups[2]])
		
			labelColors1 = rep("black",times=length(sample.label))
			library("RColorBrewer")
			continuous.color.breaks = 10
				
			plot.var = as.numeric(grp1)
			plot.var.min = min(plot.var, na.rm=T)
			plot.var.max = max(plot.var, na.rm=T)
				
			plot.var.range = plot.var.max - plot.var.min
			plot.var.interval = plot.var.range / continuous.color.breaks
				
			color.range = colorRampPalette(c("green","black","orange"))(n = continuous.color.breaks)
			plot.var.breaks = plot.var.min + plot.var.interval*(0:continuous.color.breaks)
			for (j in 1:continuous.color.breaks){
				#print(paste(plot.var.breaks[j],"to",plot.var.breaks[j+1]))
				labelColors1[(plot.var >= plot.var.breaks[j]) &(plot.var <= plot.var.breaks[j+1])] = color.range[j]
			}#end for (j in 1:continuous.color.breaks)
			
			group.levels = c(levels(as.factor(grp2)))
			color.palette <- fixed.color.palatte[3:(2+length(group.levels))]
			labelColors2 = rep("black",times=length(sample.label))
			for (i in 1:length(group.levels)){
				labelColors2[grp2 == as.character(group.levels[i])] = color.palette[i]
			}#end for (i in 1:length(group.levels))
		}else{
			grp1 = as.character(sample.description.table[,plot.groups[1]])
			grp2 = as.numeric(sample.description.table[,plot.groups[2]])

			group.levels = c(levels(as.factor(grp1)))
			color.palette <- fixed.color.palatte[1:(length(group.levels))]
			labelColors1 = rep("black",times=length(sample.label))
			for (i in 1:length(group.levels)){
				labelColors1[grp1 == as.character(group.levels[i])] = color.palette[i]
			}#end for (i in 1:length(group.levels))
			
			labelColors2 = rep("black",times=length(sample.label))
			library("RColorBrewer")
			continuous.color.breaks = 10
				
			plot.var = as.numeric(grp2)
			plot.var.min = min(plot.var, na.rm=T)
			plot.var.max = max(plot.var, na.rm=T)
				
			plot.var.range = plot.var.max - plot.var.min
			plot.var.interval = plot.var.range / continuous.color.breaks
				
			color.range = colorRampPalette(c("purple","black","cyan"))(n = continuous.color.breaks)
			plot.var.breaks = plot.var.min + plot.var.interval*(0:continuous.color.breaks)
			for (j in 1:continuous.color.breaks){
				#print(paste(plot.var.breaks[j],"to",plot.var.breaks[j+1]))
				labelColors2[(plot.var >= plot.var.breaks[j]) &(plot.var <= plot.var.breaks[j+1])] = color.range[j]
			}#end for (j in 1:continuous.color.breaks)
		}
		
		std.expr = apply(temp.cpm, 1, standardize.arr)
		if(length(deg.miRNA) < 25){
			colnames(std.expr) = deg.miRNA
		} else {
			colnames(std.expr) = rep("", length(deg.miRNA))
		}
		rownames(std.expr) = sample.label

		column_annotation <- as.matrix(deg.miRNA)
		colnames(column_annotation) <- c("")

		row_annotation <- data.frame(label1 = labelColors1, label2 = labelColors2)
		row_annotation = as.matrix(t(row_annotation))
		rownames(row_annotation) <- c(plot.groups)

		heatmap.file <- paste(comp.name,"_",pvalue.method,"_DEG_fc_",fc.cutoff,"_fdr_",fdr.cutoff,"_pval_",pvalue.cutoff,".png",sep="")
		heatmap.file = gsub(":",".",heatmap.file)
		png(file = heatmap.file)
		heatmap.3(std.expr, col=colorpanel(33, low="blue", mid="black", high="red"), density.info="none", key=TRUE,
					RowSideColors=row_annotation, trace="none", margins = c(10,15),RowSideColorsSize=4, dendrogram="both")
		if((trt.group != "continuous")&(trt.group2 != "continuous")){
					legend("topright", legend=group.levels,
							col=color.palette,
							pch=15, cex=0.7)
		}else if((trt.group == "continuous")&(trt.group2 == "continuous")){
			stop("Add code for two continuous variables")
		}else if(trt.group == "continuous"){
			legend("right",legend=c(round(plot.var.max,digits=1),rep("",length(color.range)-2),round(plot.var.min,digits=1)),
								col=rev(color.range),  pch=15, y.intersp = 0.4, cex=0.8, pt.cex=1.5)
			legend("topright", legend=group.levels, col=color.palette, pch=15, cex=0.7)
		}else{
			legend("right",legend=c(round(plot.var.max,digits=1),rep("",length(color.range)-2),round(plot.var.min,digits=1)),
								col=rev(color.range),  pch=15, y.intersp = 0.4, cex=0.8, pt.cex=1.5)
			legend("topright", legend=group.levels, col=color.palette, pch=15, cex=0.7)
		}
		dev.off()
			
		if(interaction.flag != "no"){
			temp.fc.table = as.matrix(fc.table)
			if (((trt.group == "continuous") & (trt.group2 == "continuous")) | ((trt.group != "continuous") & (trt.group2 != "continuous"))){
				temp.fc.table = temp.fc.table[,-ncol(temp.fc.table)]
			}
			temp.fc.table = temp.fc.table[status != "No Change", ]
			if(length(deg.miRNA) < 25){
				rownames(temp.fc.table) = deg.miRNA
			} else {
				rownames(temp.fc.table) = rep("",times=length(deg.miRNA))
			}
			colnames(temp.fc.table) = gsub(".:.",":",gsub("fold.change.","",colnames(temp.fc.table)))
		
			temp.fc.table[temp.fc.table < -10] = -10
			temp.fc.table[temp.fc.table > 10] = 10
		
			heatmap.file <- paste("fold_change_",comp.name,"_",pvalue.method,"_DEG_fc_",fc.cutoff,"_fdr_",fdr.cutoff,"_pval_",pvalue.cutoff,".png",sep="")
			heatmap.file = gsub(":",".",heatmap.file)
			png(file = heatmap.file)
			heatmap.2(temp.fc.table, col=colorpanel(33, low="blue", mid="black", high="red"), density.info="none", key=TRUE,
						trace="none", margins = c(20,5), cexCol=1.5)
			dev.off()
		}#end if(interaction.flag != "no")
		
	} else {
		labelColors = rep("black",times=length(sample.label))
		if(trt.group == "continuous"){
			library("RColorBrewer")
			continuous.color.breaks = 10
			
			plot.var = as.numeric(sample.description.table[,plot.groups])
			plot.var.min = min(plot.var, na.rm=T)
			plot.var.max = max(plot.var, na.rm=T)
			
			plot.var.range = plot.var.max - plot.var.min
			plot.var.interval = plot.var.range / continuous.color.breaks
			
			color.range = colorRampPalette(c("green","black","orange"))(n = continuous.color.breaks)
			plot.var.breaks = plot.var.min + plot.var.interval*(0:continuous.color.breaks)
			for (j in 1:continuous.color.breaks){
				#print(paste(plot.var.breaks[j],"to",plot.var.breaks[j+1]))
				labelColors[(plot.var >= plot.var.breaks[j]) &(plot.var <= plot.var.breaks[j+1])] = color.range[j]
			}#end for (j in 1:continuous.color.breaks)
		}else{
			group.levels = levels(as.factor(sample.description.table[,plot.groups]))
			color.palette = fixed.color.palatte[1:length(group.levels)]
			for (i in 1:length(group.levels)){
				labelColors[grp == as.character(group.levels[i])] = color.palette[i]
			}#end for (i in 1:length(group.levels))
		}

		std.expr = apply(temp.cpm, 1, standardize.arr)
		if(length(deg.miRNA) < 25){
			colnames(std.expr) = deg.miRNA
		} else {
			colnames(std.expr) = rep("", length(deg.miRNA))
		}
		rownames(std.expr) = sample.label
		
		heatmap.file <- paste(comp.name,"_DEG_fc_",fc.cutoff,"_fdr_",fdr.cutoff,"_pval_",pvalue.cutoff,".png",sep="")
		heatmap.file = gsub(":",".",heatmap.file)
		png(file = heatmap.file)
		heatmap.2(std.expr, col=colorpanel(33, low="blue", mid="black", high="red"), density.info="none", key=TRUE,
					 RowSideColors=labelColors, trace="none", margins = c(10,15))

		if(trt.group == "continuous"){
			legend("right",legend=c(round(plot.var.max,digits=1),rep("",length(color.range)-2),round(plot.var.min,digits=1)),
								col=rev(color.range),  pch=15, y.intersp = 0.4, cex=0.8, pt.cex=1.5)
		}else{
			legend("topright", group.levels, col=color.palette, pch=15)
		}
		dev.off()
	}#end else
}#end if(length(deg.miRNA) > 1)
