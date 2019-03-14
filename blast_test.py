from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
from collections import defaultdict
import csv
import sys 


result_handle = NCBIWWW.qblast("blastp","nr",open(sys.argv[1]).read(),megablast=False,expect=1000,word_size=6,hitlist_size=10,threshold=21)

for record in NCBIXML.parse(result_handle):
    if record.alignments:
        #filtering first hit base on bit_score
        record.alignments.sort(key = lambda align: -max(hsp.score for hsp in align.hsps))




