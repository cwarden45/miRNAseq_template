import sys
import re
import os
import subprocess

parameterFile = "parameters.txt"
finishedSamples = ()

ref = ""
alignmentFolder = ""
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
	
	if param == "Alignment_Folder":
		alignmentFolder = value
		
	if param == "Bowtie_Ref":
		ref = value
		
	if param == "Cluster_Email":
		email = value

if (email == "") or (email == "[required]"):
	print "Need to enter a value for 'Cluster_Email'!"
	sys.exit()
	
if (ref == "") or (ref == "[required]"):
	print "Need to enter a value for 'Bowtie_Ref'!"
	sys.exit()

if (readsFolder == "") or (readsFolder == "[required]"):
	print "Need to enter a value for 'Reads_Folder'!"
	sys.exit()
	
if (alignmentFolder == "") or (alignmentFolder == "[required]"):
	print "Need to enter a value for 'Alignment_Folder'!"
	sys.exit()

cutadaptFolder = readsFolder + "/Trimmed_Reads"
fileResults = os.listdir(cutadaptFolder)

submitAll = "master_Bowtie_queue.sh"
masterHandle = open(submitAll,"w")
text = "#!/bin/bash\n"
masterHandle.write(text)

jobCount = 0
for file in fileResults:
	result = re.search("(.*)_R1_cutadapt.fastq$",file)
	
	if result:
		sample = result.group(1)
		jobCount += 1
		
		if (sample not in finishedSamples):
			print sample
			shellScript = sample + ".sh"
			text = "sbatch " + shellScript + "\n"
			masterHandle.write(text)

			outHandle = open(shellScript, "w")
			text = "#!/bin/bash\n"
			text = text + "#SBATCH -J miB"+str(jobCount)+"\n"
			text = text + "#SBATCH --mail-type=ALL\n"
			text = text + "#SBATCH --mail-user="+email+"\n"
			text = text + "#SBATCH -n 1\n"#one thread
			text = text + "#SBATCH -N 1\n"
			text = text + "#SBATCH --mem=4g\n"
			text = text + "#SBATCH --time=24:00:00\n"
			text = text + "#SBATCH --output=miB"+str(jobCount)+".log\n\n"
			
			text = text + "module load Bowtie/1.2.2\n"
			text = text + "module load samtools/1.6\n\n"
			outHandle.write(text)
				
			outputSubfolder = alignmentFolder +"/" + sample
			text = "mkdir " + outputSubfolder + "\n"
			outHandle.write(text)
										
			read1 =cutadaptFolder + "/" + file

			#could also trim reads at this step
			alnSam = outputSubfolder + "/aligned.sam"	
			text = "bowtie --best -S " + ref + "  " + read1 + "  " + alnSam+ "\n" 
			outHandle.write(text)
									
			alnBam = outputSubfolder + "/aligned.bam"
			#can add "-F 4" to exclude unaligned reads
			text = "/opt/SAMtools/1.6/bin/samtools view -bS " + alnSam + " > " + alnBam + "\n"
			outHandle.write(text)
			
			userBam = alignmentFolder + "/" + sample + ".bam"
			text = "/opt/SAMtools/1.6/bin/samtools sort -o " + userBam + " " + alnBam + "\n"
			outHandle.write(text)
			
			text = "rm " + alnSam + "\n"
			outHandle.write(text)

			text = "rm " + alnBam + "\n"
			outHandle.write(text)
			
			text = "/opt/SAMtools/1.6/bin/samtools index " + userBam + "\n"
			outHandle.write(text)
			
			statFile = outputSubfolder + "/flagstat.txt"	
			text = "/opt/SAMtools/1.6/bin/samtools flagstat " + userBam + " > "+statFile+"\n"
			outHandle.write(text)
			
			text = "gzip " + read1
			outHandle.write(text)
