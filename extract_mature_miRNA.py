import sys
import re
import os

input = "hg38_to_hg19_mirBase_v21.gff"
description = "mirbase_hg38_liftover"
output = "hg38_to_hg19_mirBase_v21.gtf"

outHandle = open(output, 'w')

inHandle = open(input)
lines = inHandle.readlines()

miRNAhash = {}

for line in lines:
	line = re.sub("\n","",line)
	line = re.sub("\r","",line)
	
	chrResult = re.search("^chr",line)
	if chrResult:
		lineInfo = line.split("\t")
		chr = lineInfo[0]
		type = lineInfo[2]
		start = lineInfo[3]
		stop = lineInfo[4]
		strand = lineInfo[6]
		annText = lineInfo[8]
		
		if type == "miRNA":
			nameResult = re.search(";Name=(.*);Derives_from", annText)
			geneIDresult = re.search("ID=(.*);Alias", annText)
			preIDresult = re.search("Derives_from=(.*)$", annText)

			matureID = nameResult.group(1)
			
			if matureID in miRNAhash:
				miRNAhash[matureID] = miRNAhash[matureID] + 1
			else:
				miRNAhash[matureID] = 1
			
			preID = preIDresult.group(1)
			geneID = geneIDresult.group(1)
			#geneID is really primary miRNA
			#transcriptID is geneID for mature transcript
			#use gene_name for quantification
			text = chr + "\t"+description+"\texon\t" + start + "\t" + stop + "\t0.0\t" + strand + "\t.\tgene_id \"" +preID+ "\"; transcript_id \"" +geneID+ "\"; gene_name \"" +matureID+ "\"\n"
			outHandle.write(text)
			
multiMapMirna = 0

for miRNA in miRNAhash:
	if miRNAhash[miRNA] > 1:
		print miRNA + " has " + str(miRNAhash[miRNA]) + " genomic locations"
		multiMapMirna +=1
		
		
multiCopyPercent = 100 * float(multiMapMirna)/float(len(miRNAhash.keys()))		
print "Total miRNA with multiple locations: " + str(multiMapMirna) + " (" +'{0:.2f}'.format(multiCopyPercent)+ "%)"
print str(len(miRNAhash.keys())) + " mature miRNA"