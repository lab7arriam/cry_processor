#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 12 12:39:48 2018

@author: anton
"""

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
from collections import defaultdict
import csv


def blast_to_id (record):    
    result_handle = NCBIWWW.qblast("blastp", "nt", record.seq)
    blast_records = NCBIXML.parse(result_handle)
    print(blast_records)
    try:
        for blast_record in blast_records:
            print(blast_record.alignments[0])
            if blast_record.alignments:
                ID = " ".join(list(list(str(blast_record.alignments[0].title).split('|'))[4].split(' ')[i] for i in [1, 2]))
                break
    except:
        ID = "No match"
    return(ID)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='De_Bruijn_Graph visualization')
    parser.add_argument('-i', help='Please enter input sequence file', type=str)
    parser.add_argument('-o', help='Please enter output sequence file', type=str)
    parser.set_defaults(feature=True)
    args = parser.parse_args()
    i, o = args.i,args.o
    
    record_dict = defaultdict(list)
    info_dict = defaultdict(list)
    
    with open(i, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            record_dict[blast_to_id(record)].append(record) 
            info_dict[blast_to_id(record)].append(record.id)
            print(record_dict)
    New_records = []
    print(New_records)

    for  ID, record in record_dict.items():
        if len (record) == 1:
            New_records.append(SeqRecord((record_dict[ID][0].seq), id=ID))
        else:
            for contig in record:
                New_records.append(SeqRecord((contig.seq), id=ID))
    SeqIO.write(New_records, o, "fasta")

    info_file = open('/home/anton/cry_processor/info_log.tsv', 'w',newline='')
    my_writer = csv.writer(info_file, delimiter='\t') 
    
    for ID, record in info_dict.items():
        row = [None, None]
        try:
            row[0], row[1] = ID, '; '.join(record)
        except:
            row[0], row[1] = ID, record[0]
        my_writer.writerow(row)

    
    
    
