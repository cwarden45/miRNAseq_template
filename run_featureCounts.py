import sys
import re
import os

parameterFile = "parameters.txt"
threads = 4
finishedSamples = ()

alignmentFolder = ""
miRNA_gtf = ""
other_gtf = ""

inHandle = open(parameterFile)
lines = inHandle.readlines()
			
for line in lines:
	line = re.sub("\n","",line)
	line = re.sub("\r","",line)
	
	lineInfo = line.split("\t")
	param = lineInfo[0]
	value = lineInfo[1]
	
	if param == "Alignment_Folder":
		alignmentFolder = value
		
	if param == "miRNA_GTF":
		miRNA_gtf = value

	if param == "other_GTF":
		other_gtf = value		

		
fileResults = os.listdir(alignmentFolder)

for file in fileResults:
	result = re.search(".bam$",file)
	fullPath = os.path.join(alignmentFolder, file)
	
	if result:
		sample = re.sub(".bam$","",file)
		sortResult = re.search(".name.sort.bam",file)
		
		if (sample not in finishedSamples) and (not sortResult):
			print sample
			
			countsFile = alignmentFolder + "/" + sample + "/featureCounts_" + sample+ "_mature_miRNA_counts.txt"
			command = "/opt/subread-1.5.2-source/bin/featureCounts -M --fracOverlap 1 -s 1 -g gene_name -T "+str(threads)+" -a " + miRNA_gtf + " -o " + countsFile + " " + fullPath
			os.system(command)
			
			countsFile = alignmentFolder + "/" + sample + "/featureCounts_" + sample + "_other_gene_counts.txt"
			command = "/opt/subread-1.5.2-source/bin/featureCounts -s 1 -T "+str(threads)+" -a " + other_gtf + " -o " + countsFile+ " " + fullPath
			os.system(command)