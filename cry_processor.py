import subprocess
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_protein
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIWWW, NCBIXML
from collections import defaultdict
import sys
import csv


class CryProcessor:
    def __init__(self, cry_quiery, hmmer_dir, processing_flag):
        cmd_sto = subprocess.call('if [ ! -d cry_extraction ]; then mkdir cry_extraction; fi; if [ ! -d cry_extraction/logs ]; then mkdir cry_extraction/logs; fi', shell=True)
        self.cry_quiery = cry_quiery
        self.hmmer_dir = hmmer_dir
        self.processing_flag = processing_flag
        self.init_count = len(list(SeqIO.parse(open(self.cry_quiery),"fasta")))
        self.one_dom_count = 0
        self.two_dom_count = 0
        self.three_dom_count = 0

    def run_hmmer(self,queiry,out_index,model_type,dir_flag,log):
        cmd_sto = subprocess.call('cd cry_extraction; if [ ! -d {0} ]; then mkdir {0}; fi'.format(dir_flag), shell=True)
        cmd_sto = subprocess.call('{0}/binaries/hmmsearch -A {1} ./data/models/{2} {3} >> cry_extraction/logs/{4}.log'.format(self.hmmer_dir,'cry_extraction/'+dir_flag + '/'+queiry.split('.')[0]+out_index,model_type,queiry,log), shell=True) 
        cmd_fasta = subprocess.call('{0}/binaries/esl-reformat fasta {1} > {2}'.format(self.hmmer_dir,'cry_extraction/'+dir_flag+ '/'+queiry.split('.')[0]+out_index, 'cry_extraction/'+dir_flag+ '/'+queiry.split('.')[0]+out_index.replace('sto','fasta')), shell=True) 

    def find_cry(self):
        print('Searching for unprocessed cry toxins')
        self.run_hmmer(str(self.cry_quiery),'_full_extracted.sto','Complete.hmm','full_toxins', 'full_extraction')
       
    def find_domains(self):
        for i in range(1,4):
            print('Searching for domain {} of cry toxins'.format(i))
            self.run_hmmer(str(self.cry_quiery),'_D{}_extracted.sto'.format(i),'D{}.hmm'.format(i),'domains','domains_extraction')

    def cry_3D_ids_extractor(self):
        print('Exctracting cry toxins with 3 domains')
        self.dom_dict=defaultdict(list)
        for i in range(1,4):
            for record in SeqIO.parse(open('cry_extraction/{0}/{1}'.format('domains',self.cry_quiery.split('.')[0]+'_D{}_extracted.fasta'.format(i))),"fasta"):
                if '|' in record.id:
                    name = record.id.split('|')[1].split('/')[0] + '|' +'|'.join('|'.join(record.description.split('] ')[1].split(' ')).split('_'))
                else:
                    name = record.id.split('/')[0] + '|' +'|'.join('|'.join(record.description.split('] ')[1].split(' ')).split('_'))
                self.dom_dict[name].append('D{}'.format(i))
        self.one_dom_count = len([key for key in self.dom_dict if len(set(self.dom_dict[key]))==1])
        self.two_dom_count = len([key for key in self.dom_dict if len(set(self.dom_dict[key]))==2])
        return([key for key in self.dom_dict if len(set(self.dom_dict[key]))==3])


    def cry_digestor(self):
        self.id_list = self.cry_3D_ids_extractor()
        self.three_dom_count = len(self.id_list)
        print("Performing cry toxins processing")
        self.coordinate_dict = defaultdict(list)
        dom_start=int(self.processing_flag)
        for i in range(dom_start,4):
            for record in SeqIO.parse(open('cry_extraction/{0}/{1}'.format('domains',self.cry_quiery.split('.')[0]+'_D{}_extracted.fasta'.format(i))),"fasta"):
                if '|' in record.id:
                    name = record.id.split('|')[1].split('/')[0] + '|' +'|'.join('|'.join(record.description.split('] ')[1].split(' ')).split('_'))
                else:
                    name = record.id.split('/')[0] + '|' +'|'.join('|'.join(record.description.split('] ')[1].split(' ')).split('_'))
                if name in self.id_list:
                    self.coordinate_dict[name].extend(record.id.split('/')[1].split('-'))
        new_rec_list=list()        
        for record in SeqIO.parse(open(self.cry_quiery),"fasta"):
            if record.description != record.id:
                name = record.id + '|' +'|'.join('|'.join(record.description.split(' ')[1:]).split('_'))
            else:
                name= record.id + '|' +'|'.join('|'.join(record.description.split(' ')[0:]).split('_'))
            if name in self.id_list:
                 start = min([int(x) for x in self.coordinate_dict[name]])-1
                 stop = max([int(x) for x in self.coordinate_dict[name]])-1
                 new_rec_list.append(SeqRecord(Seq(str(record.seq[start:stop]),generic_protein),id=name.split('|')[0], description=" ".join(name.split('|')[1:])))
        SeqIO.write(new_rec_list,'cry_extraction/raw_processed_{}.fasta'.format(str(self.cry_quiery).split('.')[0]), "fasta")
        print('{} sequences recieved'.format(self.init_count))
        print('{} potential cry-toxins found'.format(self.one_dom_count+self.two_dom_count+self.three_dom_count))
        print('{} toxins with one domain'.format(self.one_dom_count))
        print('{} toxins with two domains'.format(self.two_dom_count))
        print('{} toxins with three domains'.format(self.three_dom_count))

    def blast_raw_output(self):
        id_match_table = defaultdict(list)
        new_ids=list()
        anotated_rec_list=list()
        print("Blasting raw output")   
        result_handle = NCBIWWW.qblast("blastp","nr",open('cry_extraction/raw_processed_{}.fasta'.format(str(self.cry_quiery).split('.')[0])).read(),megablast=False,expect=1000,word_size=6,hitlist_size=10,threshold=21)
        blast_stat = open('cry_extraction/logs/blast_alignments_{}.log'.format(str(self.cry_quiery).split('.')[0]),'w')
        for record in NCBIXML.parse(result_handle):
            blast_stat.write('********************************************************\n')
            if record.alignments:
                for alignment in record.alignments:
                    print >> blast_stat,alignment
                    blast_stat.write('\n')
                    for hsp in alignment.hsps:
                        print >> blast_stat,hsp
                        blast_stat.write('\n')
            #filtering first hit base on bit_score
                record.alignments.sort(key = lambda align: -max(hsp.score for hsp in align.hsps))
                new_ids.append([str(record.alignments[0].hit_id.split('|')[3]).strip(), str(record.alignments[0].hit_id).strip()]) 
            else:
                new_ids.append(["NO BLAST MATCH", '-']) 
        blast_stat.close()
        for init_rec, new_rec in zip(SeqIO.parse(open('cry_extraction/raw_processed_{}.fasta'.format(str(self.cry_quiery).split('.')[0])),"fasta"),new_ids):
            id_match_table[init_rec.id].append(init_rec.description)
            id_match_table[init_rec.id].extend(new_rec)
            anotated_rec_list.append(SeqRecord(Seq(str(init_rec.seq),generic_protein),id=new_rec[0]))
        SeqIO.write(anotated_rec_list,'cry_extraction/raw_blasted_{}.fasta'.format(str(self.cry_quiery).split('.')[0]), "fasta")
        csvwfile = open('cry_extraction/logs/id_match_table_{}.tsv'.format(str(self.cry_quiery).split('.')[0]), 'w')
        my_writer = csv.writer(csvwfile, delimiter='\t')
        my_writer.writerow(['Initial ID', 'Description', 'Blast hit ID', 'Description'])
        print(len(new_ids))
        print(len(id_match_table))
        for init_id in id_match_table:
            my_writer.writerow([init_id, id_match_table[init_id][0], id_match_table[init_id][1],id_match_table[init_id][2]])
        csvwfile.close()
          

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='cry_processor')
    parser.add_argument('-fi', help='Please enter full path to fasta file', metavar='File',
                        type=str, required=True)
    parser.add_argument('-hm', help='Please specify path to hmmer', metavar='Hmmer_directory',
                        type=str, required=True)
    parser.add_argument('-pr', help='Please choose processig type: 1 for extracting all domains, 2 for extratins 2-3 domains', metavar='Int',type=str, default=2)
    parser.set_defaults(feature=True)
    args = parser.parse_args()
    fi,hm,pr = args.fi, args.hm, args.pr
    pr = CryProcessor(fi, hm,pr)
    pr.find_cry()
    pr.find_domains()
    pr.cry_digestor()
    pr.blast_raw_output()
