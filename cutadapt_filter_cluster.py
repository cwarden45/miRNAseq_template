import os
import sys
import re
from Bio.Seq import Seq

#you might have to modify the code if using cutadapt on a cluster where you don't have admin rights
user = ""

parameterFile = "parameters.txt"
finishedSamples = ()

readsFolder = ""
email = ""

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

	if param == "Cluster_Email":
		email = value
		
if (readsFolder == "") or (readsFolder == "[required]"):
	print "Need to enter a value for 'Reads_Folder'!"
	sys.exit()
	
if (email == "") or (email == "[required]"):
	print "Need to enter a value for 'Cluster_Email'!"
	sys.exit()

cutadaptFolder = readsFolder + "/Trimmed_Reads"
command = "mkdir " + cutadaptFolder
os.system(command)

submitAll = "serial_cutadapt_qsub.sh"
outHandle = open(submitAll,"w")
text = "#!/bin/bash\n"
text = text + "#$ -M "+email+"\n"
text = text + "#$ -m bea\n"
text = text + "#$ -N miRNAcut\n"
text = text + "#$ -q short.q\n"
text = text + "#$ -l vf=1G\n"
text = text + "#$ -j yes\n"
text = text + "#$ -o miRNAcut.log\n"
text = text + "#$ -cwd\n"
text = text + "#$ -V\n"
outHandle.write(text)

#uncomment out if you need to install the latest version of cutadapt
text = "pip install --user --upgrade cutadapt\n"
#outHandle.write(text)

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
			
			text= "/home/"+user+"/.local/bin/cutadapt  --cut=3 " + read1 + " > "+trim1+"\n"
			outHandle.write(text)
			
			trim2 = cutadaptFolder + "/" + sample + "_R1_cutadapt.fastq"
	
			#can add '--max-n 0' to remove degenerate nucleotides
			text = "/home/"+user+"/.local/bin/cutadapt -a TCTGGAATTCTCGGGTGCCAAGGAACTCC -m 16 " +  trim1 + " > " + trim2  + "\n"
			outHandle.write(text)

			text = "rm " +  trim1 + "\n"
			outHandle.write(text)
			
			text = "gzip " + read1 + "\n"
			outHandle.write(text)
			
