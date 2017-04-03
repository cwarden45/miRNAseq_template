import re

miRNA = ""
input_file = "_RAW.txt"
output_file = "_REFORMAT.txt"

outHandle = open(output_file,"w")

inHandle = open(input_file)
line = inHandle.readline()

while line:
	line = re.sub("\n","",line)
	line = re.sub("\r","",line)
	
	geneResult = re.search("(\S+)",line)
	if geneResult:
		print "|" + geneResult.group(1) + "|"
		text = miRNA + "\t" + geneResult.group(1) + "\n"
		outHandle.write(text)
	line = inHandle.readline()