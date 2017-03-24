import sys
import re
import os
import subprocess

parameterFile = "parameters.txt"

statFile = ""
readsFolder = ""

inHandle = open(parameterFile)
lines = inHandle.readlines()
			
for line in lines:
	line = re.sub("\n","",line)
	line = re.sub("\r","",line)
	
	lineInfo = line.split("\t")
	param = lineInfo[0]
	value = lineInfo[1]

	if param == "Reads_Folder":
		readsFolder = value

	if param == "total_counts_file":
		statFile = value
		
if (readsFolder == "") or (readsFolder == "[required]"):
	print "Need to enter a value for 'Reads_Folder'!"
	sys.exit()
	
if (statFile == "") or (statFile == "[required]"):
	print "Need to enter a value for 'total_counts_file'!"
	sys.exit()
	
cutadaptFolder = readsFolder + "/Trimmed_Reads"

statHandle = open(statFile,"w")
text = "SampleID\tTotalReads\tCutadapt.Reads\tPercent.Cutadapt.Filtered\n"
statHandle.write(text)
	
fastqcFolder = readsFolder + "/QC"
fileResults = os.listdir(readsFolder)

for file in fileResults:
	result = re.search("(.*)_L\d{3}_R1_001.fastq.gz$",file)
	
	if result:
		sample = result.group(1)
		print sample
		
		#get total reads from FastQC
		fastqcPrefix = re.sub(".fastq.gz","",file)
		fastQCtext = fastqcFolder + "/" + fastqcPrefix + "_fastqc/fastqc_data.txt"
		
		inHandle = open(fastQCtext)
		line = inHandle.readline()
		
		lineCount = 0
		
		while line:
			line = re.sub("\n","",line)
			line = re.sub("\r","",line)
			
			lineCount += 1
			
			if lineCount == 7:
				
				totalResult = re.search("Total Sequences\t(\d+)",line)
				if totalResult:
					totalReads = int(totalResult.group(1))
				else:
					print "Problem parsing FastQC file!\n"
					sys.exit()
			
			line = inHandle.readline()
		
		inHandle.close()
		
		#get number of filtered forward reads
		cutadaptFQ = cutadaptFolder + "/" + sample + "_R1_cutadapt.fastq"
		
		command = "wc -l " + cutadaptFQ
		wcText = subprocess.check_output(command, shell=True)
		countResult = re.search("^\s+(\d+)",wcText)
		if countResult:
			fqLines = int(countResult.group(1))
		else:
			print "Problem parsing " + wcText
			sys.exit()
		cutadaptReads = fqLines / 4
		percentKept = 100 * float(cutadaptReads)/float(totalReads)
		
		text = sample + "\t" + str(totalReads) + "\t" + str(cutadaptReads) + "\t" + '{0:.2f}'.format(percentKept) + "%\n"
		statHandle.write(text)