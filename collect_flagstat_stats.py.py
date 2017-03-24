import sys
import re
import os
import subprocess

parameterFile = "parameters.txt"

statFile = ""
alignmentFolder = ""

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

	if param == "aligned_stats_file":
		statFile = value
		
if (alignmentFolder == "") or (alignmentFolder == "[required]"):
	print "Need to enter a value for 'Alignment_Folder'!"
	sys.exit()

if (statFile == "") or (statFile == "[required]"):
	print "Need to enter a value for 'aligned_stats_file'!"
	sys.exit()
	
statHandle = open(statFile,"w")
text = "SampleID\tBowtie.Cutadapt.Reads\tBowtie.Aligned.Reads\tPercent.Aligned\n"
statHandle.write(text)
	
fileResults = os.listdir(alignmentFolder)

for file in fileResults:
	result = re.search("(.*).bam$",file)
	
	if result:
		sample = result.group(1)
		print sample
		
		flagstat_file = alignmentFolder + "/" + sample + "/flagstat.txt"
		inHandle = open(flagstat_file)
		line = inHandle.readline()
		
		lineCount = 0
		
		totalReads = 0
		alignedReads = 0
		
		while line:
			line = re.sub("\n","",line)
			line = re.sub("\r","",line)
			
			lineCount += 1
			
			if lineCount == 1:
				countResult = re.search("^(\d+)",line)
				
				if countResult:
					totalReads = int(countResult.group(1))
				else:
					print "Problem parsing total reads from line:\n"
					print line
					sys.exit()
			elif lineCount == 3:
				countResult = re.search("^(\d+)",line)
				
				if countResult:
					alignedReads = int(countResult.group(1))
				else:
					print "Problem parsing aligned reads from line:\n"
					print line
					sys.exit()
					
			line = inHandle.readline()
		
		inHandle.close()
		
		if alignedReads != 0:
			percentAligned = 100 * float(alignedReads)/float(totalReads)
		else:
			percentAligned=0
		text = sample + "\t" + str(totalReads) + "\t" + str(alignedReads) + "\t" + '{0:.2f}'.format(percentAligned) + "%\n"
		statHandle.write(text)