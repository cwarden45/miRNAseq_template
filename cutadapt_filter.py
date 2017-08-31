import os
import sys
import re
from Bio.Seq import Seq

parameterFile = "parameters.txt"
finishedSamples = ()

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
		
if (readsFolder == "") or (readsFolder == "[required]"):
	print "Need to enter a value for 'Reads_Folder'!"
	sys.exit()

cutadaptFolder = readsFolder + "/Trimmed_Reads"
command = "mkdir " + cutadaptFolder
os.system(command)

fileResults = os.listdir(readsFolder)

jobCount = 0
for file in fileResults:
	result = re.search("(.*)_S\d+_L\d{3}_R1_001.fastq$",file)
	
	if result:
		jobCount += 1
		sample = result.group(1)
		if sample not in finishedSamples:
			print sample
			
			read1 = readsFolder + "/" + file

			trim1 = cutadaptFolder + "/"+ sample + "_R1_removeFirst3.fastq"
			
			command= "cutadapt --cut=3 " + read1 + " > "+trim1
			os.system(command)
			
			trim2 = cutadaptFolder + "/" + sample + "_R1_cutadapt.fastq"
	
			#can add '--max-n 0' to remove degenerate nucleotides
			command = "cutadapt -a TCTGGAATTCTCGGGTGCCAAGGAACTCC -m 16 " +  trim1 + " > " + trim2
			os.system(command)

			command = "rm " +  trim1
			os.system(command)

			command = "gzip " + read1
			os.system(command)		
