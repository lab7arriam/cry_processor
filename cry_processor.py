import subprocess
import argparse
from Bio import SeqIO, Entrez
from Bio.Seq import Seq
from Bio.Alphabet import generic_protein
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import csv
import time



class CryProcessor:
    def __init__(self, cry_quiery, hmmer_dir, processing_flag, hm_threads):
        cmd_sto = subprocess.call('if [ ! -d cry_extraction ]; then mkdir cry_extraction; fi; if [ ! -d cry_extraction/logs ]; then mkdir cry_extraction/logs; fi', shell=True)
        self.cry_quiery = cry_quiery
        self.hmmer_dir = hmmer_dir
        self.processing_flag = processing_flag
        self.hm_threads = hm_threads
        self.init_count = len(list(SeqIO.parse(open(self.cry_quiery),"fasta")))
        self.one_dom_count = 0
        self.two_dom_count = 0
        self.three_dom_count = 0

    def run_hmmer(self,queiry,out_index,model_type,dir_flag,log,hm_threads):
        cmd_sto = subprocess.call('cd cry_extraction; if [ ! -d {0} ]; then mkdir {0}; fi'.format(dir_flag), shell=True)
        cmd_sto = subprocess.call('{0}/binaries/hmmsearch --cpu {5} -A {1} ./data/models/{2} {3} >> cry_extraction/logs/{4}.log'.format(self.hmmer_dir,'cry_extraction/'+dir_flag + '/'+queiry.split('.')[0]+out_index,model_type,queiry,log, hm_threads), shell=True) 
        cmd_fasta = subprocess.call('{0}/binaries/esl-reformat fasta {1} > {2}'.format(self.hmmer_dir,'cry_extraction/'+dir_flag+ '/'+queiry.split('.')[0]+out_index, 'cry_extraction/'+dir_flag+ '/'+queiry.split('.')[0]+out_index.replace('sto','fasta')), shell=True) 

    def find_cry(self):
        print('Searching for unprocessed cry toxins')
        self.run_hmmer(str(self.cry_quiery),'_full_extracted.sto','Complete.hmm','full_toxins', 'full_extraction',self.hm_threads)
       
    def find_domains(self):
        for i in range(1,4):
            print('Searching for domain {} of cry toxins'.format(i))
            self.run_hmmer(str(self.cry_quiery),'_D{}_extracted.sto'.format(i),'D{}.hmm'.format(i),'domains','domains_extraction',self.hm_threads)

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
        print('{} potential cry toxins found'.format(self.one_dom_count+self.two_dom_count+self.three_dom_count))
        print('{} toxins with one domain'.format(self.one_dom_count))
        print('{} toxins with two domains'.format(self.two_dom_count))
        print('{} toxins with three domains'.format(self.three_dom_count))

    def annotate_raw_output(self):
        new_records = list()
        new_ids=dict()
        un_count=0
        total_count=0
        print("Annotating raw output with diamond")  
        cmd_pre_dia = subprocess.call('cd cry_extraction; cp ../data/diamond_data/cry_nomenclature.dmnd .; touch diamond.log', shell=True) 
        cmd_dia = subprocess.call('cd cry_extraction;../diamond blastp -d cry_nomenclature -q raw_processed_{0}.fasta -o diamond_matches_{0}.txt --al aligned_{0}.fa --un unaligned_{0}.fa --max-target-seqs 1 --log --verbose 2>> diamond.log; rm cry_nomenclature.dmnd; mv *.log logs/; mv *gned* logs/'.format(str(self.cry_quiery).split('.')[0]), shell=True)
        with open("cry_extraction/diamond_matches_{}.txt".format(str(self.cry_quiery).split('.')[0]), 'r') as csv_file:
            my_reader = csv.reader(csv_file, delimiter='\t') 
            for row in my_reader:
                total_count+=1
                if float(row[2])<100.0:
                    new_ids[row[0]]=row[1]+'|'+ str(row[2])
                    un_count+=1
        for init_rec in SeqIO.parse(open('cry_extraction/raw_processed_{}.fasta'.format(str(self.cry_quiery).split('.')[0])),"fasta"):
           if init_rec.id in new_ids.keys():
               new_records.append(SeqRecord(Seq(str(init_rec.seq),generic_protein),id=new_ids[init_rec.id],description=init_rec.description))
        print('{} sequences matched with database'.format(total_count))
        print('{} toxins different from database found'.format(un_count))
        SeqIO.write(new_records,'cry_extraction/unique_{}.fasta'.format(str(self.cry_quiery).split('.')[0]), "fasta")
       
    def make_summary_table(self):
        print('Searching for the metadata')
        summary_dict=defaultdict(dict)
        with open("cry_extraction/diamond_matches_{}.txt".format(str(self.cry_quiery).split('.')[0]), 'r') as csv_file:
            my_reader = csv.reader(csv_file, delimiter='\t') 
            for row in my_reader:
                summary_dict[row[0]]=defaultdict(list)
                summary_dict[row[0]]['init']=row[1:3]
        for init_rec in SeqIO.parse(open('cry_extraction/raw_processed_{}.fasta'.format(str(self.cry_quiery).split('.')[0])),"fasta"):
           if init_rec.id in summary_dict.keys():
               summary_dict[init_rec.id]['init'].append(init_rec.description)
        for key in summary_dict:
           handle = Entrez.efetch(db="protein",rettype='ipg',retmode='text', id =key)
           handle_list=[el.split('\t') for el in handle.read().split('\n')]
           hit_counter=0
           for i in range(len(handle_list)-1):
              if len(handle_list[i+1])>2:
                  summary_dict[key]['hit'+str(hit_counter)]=handle_list[i+1]
                  hit_counter+=1
           time.sleep(3)
        with open("cry_extraction/annotation_table_{}.tsv".format(str(self.cry_quiery).split('.')[0]), 'w') as csv_file:
            my_writer = csv.writer(csv_file, delimiter='\t') 
            init_row = ['protein_id', 'initial_description', 'top_cry_hit', 'cry_identity', 'source', 'nucl_accession', 'start','stop', 'strand','ipg_prot_id','ipg_prot_name', 'organism', 'strain','assembly']
            my_writer.writerow(init_row)
            for key in summary_dict:
                iter_num=len(summary_dict[key].keys())
                for i in range(iter_num-1):
                    if i==0:
                        row=[key]+[summary_dict[key]['init'][2]]+[summary_dict[key]['init'][0]]+[summary_dict[key]['init'][1]]+summary_dict[key]['hit'+str(i)][1:]
                    else:
                        row=['--']*4+summary_dict[key]['hit'+str(i)][1:]
                    my_writer.writerow(row)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='cry_processor')
    parser.add_argument('-fi', help='Please enter full path to fasta file', metavar='File',
                        type=str, required=True)
    parser.add_argument('-hm', help='Please specify path to hmmer', metavar='Hmmer_directory',
                        type=str, required=True)
    parser.add_argument('-pr', help='Please choose processig type: 1 for extracting all domains, 2 for extratins 2-3 domains', metavar='Int',type=str, default=1)
    parser.add_argument('-th', help='Please specify number of threads for hmmer', metavar='Int',type=str, default=1)
    parser.set_defaults(feature=True)
    args = parser.parse_args()
    fi,hm,pr,th = args.fi, args.hm, args.pr,args.th
    pr = CryProcessor(fi, hm,pr, th)
    #pr.find_cry()
    #pr.find_domains()
    #pr.cry_digestor()
    #pr.annotate_raw_output()
    pr.make_summary_table()
