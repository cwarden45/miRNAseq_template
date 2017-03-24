import sys
import re
import os

input_file = "hsa_hg38_mirBase_v21.gff3"
output_file = "hg38_to_hg19_mirBase_v21.gff"
chain_file = "hg38ToHg19.over.chain"
unmapped_file = "hg38_v21_unmapped.gff"

command = "liftOver -gff " + input_file + " " + chain_file + " " + output_file + " " + unmapped_file
os.system(command)