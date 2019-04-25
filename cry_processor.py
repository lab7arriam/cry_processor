import subprocess
import argparse
from Bio import SeqIO, Entrez
from Bio.Seq import Seq
from Bio.Alphabet import generic_protein
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import csv
import time
import os
import sys
import re

class CryProcessor:
    def __init__(self,quiery_dir, cry_quiery, hmmer_dir, processing_flag, hm_threads, email, regime):
        self.home_dir = ('/').join(os.path.realpath(__file__).split('/')[0:len(os.path.realpath(__file__).split('/'))-1])
        self.cry_quiery = cry_quiery
        self.hmmer_dir = hmmer_dir
        self.processing_flag = processing_flag
        self.hm_threads = hm_threads
        self.quiery_dir = quiery_dir
        self.init_count = len(list(SeqIO.parse(open(self.cry_quiery),"fasta")))
        self.one_dom_count = 0
        self.two_dom_count = 0
        self.three_dom_count = 0
        self.email = email
        self.regime = regime
        cmd_init = subprocess.call('if [ ! -d $PWD/{0} ]; then mkdir $PWD/{0}; fi; if [ ! -d $PWD/{0}/cry_extraction ]; then mkdir $PWD/{0}/cry_extraction; fi; if [ ! -d $PWD/{0}/cry_extraction/logs ]; then mkdir $PWD/{0}/cry_extraction/logs; fi'.format(self.quiery_dir), shell=True)

    def run_hmmer(self,queiry,out_index,model_type,dir_flag,log,hm_threads,quiery_dir):
        if self.hmmer_dir:
            cmd_make_dirs = subprocess.call('cd $PWD/{1}/cry_extraction; if [ ! -d {0} ]; then mkdir {0}; fi'.format(dir_flag, quiery_dir), shell=True)

            cmd_search = subprocess.call('{0}/binaries/hmmsearch --cpu {5} -A {1} {6}/data/models/{2} {3} >> $PWD/{7}/cry_extraction/logs/{4}.log'.format(self.hmmer_dir,'$PWD/{0}/cry_extraction/'.format(self.quiery_dir)+dir_flag + '/'+queiry.split('/')[len(queiry.split('/'))-1].split('.')[0]+out_index,model_type,queiry,log, hm_threads, self.home_dir, self.quiery_dir), shell=True) 

            cmd_fasta = subprocess.call('{0}/binaries/esl-reformat fasta {1} > {2}'.format(self.hmmer_dir,'$PWD/{0}/cry_extraction/'.format(self.quiery_dir)+dir_flag+ '/'+queiry.split('/')[len(queiry.split('/'))-1].split('.')[0]+out_index, '$PWD/{0}/cry_extraction/'.format(self.quiery_dir)+dir_flag+ '/'+queiry.split('/')[len(queiry.split('/'))-1].split('.')[0]+out_index.replace('sto','fasta')), shell=True) 

        else:
            cmd_make_dirs = subprocess.call('cd $PWD/{1}/cry_extraction; if [ ! -d {0} ]; then mkdir {0}; fi'.format(dir_flag, quiery_dir), shell=True)

            cmd_search = subprocess.call('hmmsearch --cpu {4} -A {0} {5}/data/models/{1} {2} >> $PWD/{6}/cry_extraction/logs/{3}.log'.format('$PWD/{0}/cry_extraction/'.format(self.quiery_dir)+dir_flag + '/'+queiry.split('/')[len(queiry.split('/'))-1].split('.')[0]+out_index,model_type,queiry,log, hm_threads, self.home_dir, self.quiery_dir), shell=True) 

            cmd_fasta = subprocess.call('esl-reformat fasta {0} > {1}'.format('$PWD/{0}/cry_extraction/'.format(self.quiery_dir)+dir_flag+ '/'+queiry.split('/')[len(queiry.split('/'))-1].split('.')[0]+out_index, '$PWD/{0}/cry_extraction/'.format(self.quiery_dir)+dir_flag+ '/'+queiry.split('/')[len(queiry.split('/'))-1].split('.')[0]+out_index.replace('sto','fasta')), shell=True) 

    def find_cry(self, qiuery):
        print('Searching for unprocessed cry toxins')
        self.run_hmmer(str(qiuery),'_full_extracted.sto','Complete.hmm','full_toxins', 'full_extraction',self.hm_threads, self.quiery_dir)
       
    def find_domains(self, qiuery):
        for i in range(1,4):
            print('Searching for domain {} of cry toxins'.format(i))
            self.run_hmmer(str(qiuery),'_D{}_extracted.sto'.format(i),'D{}.hmm'.format(i),'domains','domains_extraction', self.hm_threads, self.quiery_dir)

    def cry_3D_ids_extractor(self):
        print('Exctracting cry toxins with 3 domains')
        self.dom_dict=defaultdict(list)
        for i in range(1,4):
            for record in SeqIO.parse(open('{3}/{2}/cry_extraction/{0}/{1}'.format('domains',self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0]+'_D{}_extracted.fasta'.format(i),self.quiery_dir,subprocess.Popen(['pwd'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0].strip())),"fasta"):
                if '|' in record.id:
                    name = record.id.split('|')[1].split('/')[0] + '|' +'|'.join('|'.join(record.description.split('] ')[1].split(' ')).split('_'))
                else:
                    name = record.id.split('/')[0] + '|' +'|'.join('|'.join(record.description.split('] ')[1].split(' ')).split('_'))
                self.dom_dict[name].append('D{}'.format(i))
        self.one_dom_count = len([key for key in self.dom_dict if len(set(self.dom_dict[key]))==1])
        self.two_dom_count = len([key for key in self.dom_dict if len(set(self.dom_dict[key]))==2])
        return([key for key in self.dom_dict if len(set(self.dom_dict[key]))==3])

    def cry_digestor(self):
        if self.regime == 'do':
            self.find_domains(self.cry_quiery)
            self.id_list = self.cry_3D_ids_extractor()
            self.three_dom_count = len(self.id_list)
            print("Performing cry toxins processing")
            self.coordinate_dict = defaultdict(list)
            dom_start=int(self.processing_flag)
            for i in range(dom_start,4):
                for record in SeqIO.parse(open('{3}/{2}/cry_extraction/{0}/{1}'.format('domains',self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0]+'_D{}_extracted.fasta'.format(i),self.quiery_dir,subprocess.Popen(['pwd'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0].strip())),"fasta"):
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
            
            SeqIO.write(new_rec_list,'/{0}/{1}/cry_extraction/raw_processed_{2}.fasta'.format(subprocess.Popen(['pwd'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0].strip(),self.quiery_dir,self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0]), "fasta")
            print('{} sequences recieved'.format(self.init_count))
            print('{} potential cry toxins found'.format(self.one_dom_count+self.two_dom_count+self.three_dom_count))
            print('{} toxins with one domain'.format(self.one_dom_count))
            print('{} toxins with two domains'.format(self.two_dom_count))
            print('{} toxins with three domains'.format(self.three_dom_count))

        if self.regime == 'fd':
            self.find_cry(self.cry_quiery)
            cmd_take_file = subprocess.call("cp {0}/{1}/cry_extraction/full_toxins/{2}_full_extracted.fasta {0}/{1}/cry_extraction/; mv {0}/{1}/cry_extraction/{2}_full_extracted.fasta {0}/{1}/cry_extraction/{2}.fasta; sed -i 's/[0-9]*|//g' {0}/{1}/cry_extraction/{2}.fasta; sed -i 's/\/[0-9]*-[0-9]*//g' {0}/{1}/cry_extraction/{2}.fasta; sed -i 's/ \[subseq from]\ / /g' {0}/{1}/cry_extraction/{2}.fasta".format(subprocess.Popen(['pwd'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0].strip(),self.quiery_dir,self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0]),shell=True)
            self.find_domains('{0}/{1}/cry_extraction/{2}.fasta'.format(subprocess.Popen(['pwd'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0].strip(),self.quiery_dir,self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0]))
            self.id_list = self.cry_3D_ids_extractor()
            self.three_dom_count = len(self.id_list)
            print("Performing cry toxins processing")
            self.coordinate_dict = defaultdict(list)
            self.first_filter_count = len(list(SeqIO.parse(open("{0}/{1}/cry_extraction/{2}.fasta".format(subprocess.Popen(['pwd'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0].strip(),self.quiery_dir,self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0])),"fasta")))
            dom_start=int(self.processing_flag)
            for i in range(dom_start,4):
                for record in SeqIO.parse(open('{3}/{2}/cry_extraction/{0}/{1}'.format('domains',self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0]+'_D{}_extracted.fasta'.format(i),self.quiery_dir,subprocess.Popen(['pwd'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0].strip())),"fasta"):
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
            SeqIO.write(new_rec_list,'/{0}/{1}/cry_extraction/raw_processed_{2}.fasta'.format(subprocess.Popen(['pwd'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0].strip(),self.quiery_dir,self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0]), "fasta")
            print('{} sequences recieved'.format(self.init_count))
            print('{} sequences after first search'.format(self.first_filter_count))
            print('{} potential cry toxins found'.format(self.one_dom_count+self.two_dom_count+self.three_dom_count))
            print('{} toxins with one domain'.format(self.one_dom_count))
            print('{} toxins with two domains'.format(self.two_dom_count))
            print('{} toxins with three domains'.format(self.three_dom_count))

    def annotate_raw_output(self):
        new_records = list()
        self.new_ids=dict()
        un_count=0
        total_count=0
        print("Annotating raw output with diamond")  
        cmd_pre_dia = subprocess.call('cd /{0}/{1}/cry_extraction; cp {2}/data/diamond_data/cry_nomenclature.dmnd .; touch diamond.log'.format(subprocess.Popen(['pwd'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0].strip(),self.quiery_dir,self.home_dir), shell=True) 
        cmd_dia = subprocess.call('cd /{1}/{2}/cry_extraction;{3}/diamond blastp -d cry_nomenclature -q raw_processed_{0}.fasta -o diamond_matches_{0}.txt --al aligned_{0}.fa --un unaligned_{0}.fa --max-target-seqs 1 --log --verbose 2>> diamond.log; rm cry_nomenclature.dmnd; mv *.log logs/; mv *gned* logs/'.format(self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0],subprocess.Popen(['pwd'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0].strip(),self.quiery_dir,self.home_dir), shell=True)
        with open("/{0}/{1}/cry_extraction/diamond_matches_{2}.txt".format(subprocess.Popen(['pwd'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0].strip(),self.quiery_dir,self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0]), 'r') as csv_file:
            my_reader = csv.reader(csv_file, delimiter='\t') 
            for row in my_reader:
                total_count+=1
                if float(row[2])<100.0:
                    self.new_ids[row[0]]=row[1]+'|'+ str(row[2])
                    un_count+=1
        for init_rec in SeqIO.parse(open('/{0}/{1}/cry_extraction/raw_processed_{2}.fasta'.format(subprocess.Popen(['pwd'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0].strip(),self.quiery_dir,self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0])),"fasta"):
           if init_rec.id in self.new_ids.keys():
               new_records.append(SeqRecord(Seq(str(init_rec.seq),generic_protein),id=self.new_ids[init_rec.id],description=init_rec.description))
        print('{} sequences matched with database'.format(total_count))
        print('{} toxins different from database found'.format(un_count))
        SeqIO.write(new_records,'/{0}/{1}/cry_extraction/unique_{2}.fasta'.format(subprocess.Popen(['pwd'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0].strip(),self.quiery_dir,self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0]), "fasta")
        if self.regime == 'fd':
            cmd_clean_up = subprocess.call("rm {0}/{1}/cry_extraction/{2}.fasta".format(subprocess.Popen(['pwd'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0].strip(),self.quiery_dir,self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0]),shell=True)
       
    def make_summary_table(self):
        print('Searching for the metadata')
        Entrez.email = "{}".format(self.email)
        summary_dict=defaultdict(dict)
        with open("/{0}/{1}/cry_extraction/diamond_matches_{2}.txt".format(subprocess.Popen(['pwd'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0].strip(),self.quiery_dir,self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0]), 'r') as csv_file:
            my_reader = csv.reader(csv_file, delimiter='\t') 
            for row in my_reader:
                summary_dict[row[0]]=defaultdict(list)
                summary_dict[row[0]]['init']=row[1:3]
        making_smb = subprocess.call("touch /{0}/{1}/cry_extraction/annotation_table_{2}.tsv".format(subprocess.Popen(['pwd'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0].strip(),self.quiery_dir,self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0]), shell = True)
        for init_rec in SeqIO.parse(open('/{0}/{1}/cry_extraction/raw_processed_{2}.fasta'.format(subprocess.Popen(['pwd'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0].strip(),self.quiery_dir,self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0])),"fasta"):
           if init_rec.id in summary_dict.keys():
               summary_dict[init_rec.id]['init'].append(init_rec.description)
        init_row = ['protein_id', 'initial_description', 'top_cry_hit', 'cry_identity', 'source', 'nucl_accession', 'start','stop', 'strand','ipg_prot_id','ipg_prot_name', 'organism', 'strain','assembly\n']
        f = open("/{0}/{1}/cry_extraction/annotation_table_{2}.tsv".format(subprocess.Popen(['pwd'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0].strip(),self.quiery_dir,self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0]),"w")
        f.write('\t'.join(init_row))
        f.close() 
        for key in summary_dict:
            handle = Entrez.efetch(db="protein",rettype='ipg',retmode='text', id =key)
            handle_list=[el.split('\t') for el in handle.read().split('\n')]
            hit_counter=0
            for i in range(len(handle_list)-1):
               if len(handle_list[i+1])>2:
                   summary_dict[key]['hit'+str(hit_counter)]=handle_list[i+1]
                   hit_counter+=1
            iter_num=len(summary_dict[key].keys())
            for i in range(iter_num-1):
                if i==0:
                    row=[key]+[summary_dict[key]['init'][2]]+[summary_dict[key]['init'][0]]+[summary_dict[key]['init'][1]]+summary_dict[key]['hit'+str(i)][1:]
                else:
                    row=['--']*4+summary_dict[key]['hit'+str(i)][1:]
                f = open("/{0}/{1}/cry_extraction/annotation_table_{2}.tsv".format(subprocess.Popen(['pwd'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0].strip(),self.quiery_dir,self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0]),"a+")
                f.write('\t'.join(row) + '\n') 
                f.close()
            time.sleep(3)

    def upload_nucl(self):
        print('uploading nucleotide sequences')
        Entrez.email = "{}".format(self.email)
        keys_for_nucl = defaultdict(list)
        with open("/{0}/{1}/cry_extraction/annotation_table_{2}.tsv".format(subprocess.Popen(['pwd'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0].strip(),self.quiery_dir,self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0]), 'r') as csvfile:
             my_reader = csv.reader(csvfile, delimiter='\t') 
             for row in my_reader:
                 if row[0] in self.new_ids.keys():
                     keys_for_nucl[row[0]].append('|'.join('|'.join((row[1]).split(' ')).split('_')))
                     keys_for_nucl[row[0]].extend([row[5], row[6], row[7], row[8]])
        print(keys_for_nucl)
        for key in keys_for_nucl:
            handle = Entrez.efetch(db="nucleotide",rettype='fasta',retmode='text', id = keys_for_nucl[key][1])
            fasta_rec = SeqIO.read(handle, "fasta")
            handle.close()
            if keys_for_nucl[key][4] == '+':
                print(fasta_rec.seq[int(keys_for_nucl[key][2])-1:int(keys_for_nucl[key][3])])
                print(len(fasta_rec.seq[int(keys_for_nucl[key][2])-1:int(keys_for_nucl[key][3])]),len(fasta_rec.seq[int(keys_for_nucl[key][2])-1:int(keys_for_nucl[key][3])])/3)
            elif keys_for_nucl[key][4] == '-':
                print(fasta_rec.seq[int(keys_for_nucl[key][2])-1:int(keys_for_nucl[key][3])].reverse_complement())
                print(len(fasta_rec.seq[int(keys_for_nucl[key][2])-1:int(keys_for_nucl[key][3])].reverse_complement()),len(fasta_rec.seq[int(keys_for_nucl[key][2])-1:int(keys_for_nucl[key][3])].reverse_complement())/3)
            print(fasta_rec.id)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='cry_processor')
    parser.add_argument('-fi', help='Please enter full path to fasta file', metavar='File',
                        type=str, required=True)
    parser.add_argument('-hm', help='Please specify path to hmmer if your hmmer is installed locally', metavar='Hmmer_directory',
                        type=str, default='')
    parser.add_argument('-pr', help='Please choose processig type: 1 for extracting all domains, 2 for extrating 2-3 domains', metavar='Int',type=str, default=1)
    parser.add_argument('-th', help='Please specify number of threads for hmmer', metavar='Int',type=str, default=1)
    parser.add_argument('-ma', help='Please enter e-mail address for NCBI annotation', metavar='Str',type=str, default='')
    parser.add_argument('-od', help='Please specify output directory', metavar='Str',type=str, required=True)
    parser.add_argument('-r', help='Please choose pipeline type: do - domain only search with subsequent unioning; fd - searching for potential cry-toxins with subsequent processing', metavar='Str',type=str, default='do')
    parser.add_argument('--annotate', '-a',action='store_true',help='make final NCBI annotation')
    parser.add_argument('-nu', help='Please specify the way to upload nucleotide records: fn - uploading full sequences, pn - uploading processed subsequences', metavar='Nucl_uploading_type', type=str, default='')
    parser.set_defaults(feature=True)
    args = parser.parse_args()
    od,fi,hm,pr,th, ma, r, a, nu = args.od, args.fi, args.hm, args.pr,args.th, args.ma, args.r, args.annotate, args.nu
    pr = CryProcessor(od, fi, hm,pr, th, ma, r)
    pr.cry_digestor()
    pr.annotate_raw_output()
    #if a: 
    #    pr.make_summary_table()
    if nu:
        pr.upload_nucl()
