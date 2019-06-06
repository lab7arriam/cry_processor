import subprocess
import argparse
from Bio import SeqIO, Entrez
from Bio.Seq import Seq
from Bio.Alphabet import generic_protein, generic_dna
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import csv
import time
import os
import sys
import re

class CryProcessor:
    def __init__(self,quiery_dir, cry_quiery, hmmer_dir, processing_flag, hm_threads, email, regime, nucl_type, annot_flag,kmer_size,meta_flag,forw,rev):
        self.home_dir = ('/').join(os.path.realpath(__file__).split('/')[0:len(os.path.realpath(__file__).split('/'))-1])
        self.cry_quiery = cry_quiery
        self.hmmer_dir = hmmer_dir
        self.processing_flag = processing_flag
        self.hm_threads = hm_threads
        self.quiery_dir = quiery_dir
        if self.cry_quiery:
            self.init_count = re.sub("b",'',re.sub("\'",'',str(subprocess.check_output("grep '>' {} | wc -l".format(self.cry_quiery), shell =True).strip())))
            self.one_dom_count = 0
            self.two_dom_count = 0
            self.three_dom_count = 0
        self.email = email
        self.regime = regime
        self.nucl_type = nucl_type
        self.annot_flag = annot_flag
        self.kmer_size = kmer_size
        self.meta_flag = meta_flag
        self.forw = forw
        self.rev = rev
        cmd_init = subprocess.call('if [ ! -d $PWD/{0} ]; then mkdir $PWD/{0}; fi; if [ ! -d $PWD/{0}/cry_extraction ]; then mkdir $PWD/{0}/cry_extraction; fi; if [ ! -d $PWD/{0}/cry_extraction/logs ]; then mkdir $PWD/{0}/cry_extraction/logs; fi'.format(self.quiery_dir), shell=True)
        self.main_dir = re.sub("b",'',re.sub("\'",'',str(subprocess.Popen(['pwd'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0].strip())))


    def run_hmmer(self,queiry,out_index,model_type,dir_flag,log,hm_threads,quiery_dir):
        if self.hmmer_dir:
            cmd_make_dirs = subprocess.call('cd $PWD/{1}/cry_extraction; if [ ! -d {0} ]; then mkdir {0}; fi'.format(dir_flag, quiery_dir), shell=True)

            cmd_search = subprocess.call('{0}/binaries/hmmsearch --cpu {5} -A {1} {6}/data/models/{2} {3} >> $PWD/{7}/cry_extraction/logs/{4}.log'.format(self.hmmer_dir,'$PWD/{0}/cry_extraction/'.format(self.quiery_dir)+dir_flag + '/'+queiry.split('/')[len(queiry.split('/'))-1].split('.')[0]+out_index,model_type,queiry,log, hm_threads, self.home_dir, self.quiery_dir), shell=True) 

            cmd_fasta = subprocess.call('{0}/binaries/esl-reformat fasta {1} > {2}'.format(self.hmmer_dir,'$PWD/{0}/cry_extraction/'.format(self.quiery_dir)+dir_flag+ '/'+queiry.split('/')[len(queiry.split('/'))-1].split('.')[0]+out_index, '$PWD/{0}/cry_extraction/'.format(self.quiery_dir)+dir_flag+ '/'+queiry.split('/')[len(queiry.split('/'))-1].split('.')[0]+out_index.replace('sto','fasta')), shell=True) 

        else:
            cmd_make_dirs = subprocess.call('cd $PWD/{1}/cry_extraction; if [ ! -d {0} ]; then mkdir {0}; fi'.format(dir_flag, quiery_dir), shell=True)

            cmd_search = subprocess.call('{7}hmmsearch --cpu {4} -A {0} {5}/data/models/{1} {2} >> $PWD/{6}/cry_extraction/logs/{3}.log'.format('$PWD/{0}/cry_extraction/'.format(self.quiery_dir)+dir_flag + '/'+queiry.split('/')[len(queiry.split('/'))-1].split('.')[0]+out_index,model_type,queiry,log, hm_threads, self.home_dir, self.quiery_dir,self.home_dir+'/include/'), shell=True) 

            cmd_fasta = subprocess.call('{2}esl-reformat fasta {0} > {1}'.format('$PWD/{0}/cry_extraction/'.format(self.quiery_dir)+dir_flag+ '/'+queiry.split('/')[len(queiry.split('/'))-1].split('.')[0]+out_index, '$PWD/{0}/cry_extraction/'.format(self.quiery_dir)+dir_flag+ '/'+queiry.split('/')[len(queiry.split('/'))-1].split('.')[0]+out_index.replace('sto','fasta'),self.home_dir+'/include/'), shell=True) 

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
            for record in SeqIO.parse(open('{3}/{2}/cry_extraction/{0}/{1}'.format('domains',self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0]+'_D{}_extracted.fasta'.format(i),self.quiery_dir,self.main_dir)),"fasta"):
                if '|' in record.id:
                    name = record.id.split('|')[1].split('/')[0] + '|' +'|'.join('|'.join(record.description.split('] ')[1].split(' ')).split('_'))
                else:
                    name = record.id.split('/')[0] + '|' +'|'.join('|'.join(record.description.split('] ')[1].split(' ')).split('_'))
                if 'D{}'.format(i) not in self.dom_dict[name]:
                    self.dom_dict[name].append('D{}'.format(i))
        self.one_dom_count = len([key for key in self.dom_dict if len(set(self.dom_dict[key]))==1])
        self.two_dom_count = len([key for key in self.dom_dict if len(set(self.dom_dict[key]))==2])
        return([key for key in self.dom_dict if len(set(self.dom_dict[key]))>=3])

    def cry_digestor(self):
        if self.regime == 'do':
            self.find_domains(self.cry_quiery)
            self.id_list = self.cry_3D_ids_extractor()
            self.three_dom_count = len(self.id_list)
            print("Performing cry toxins processing")
            self.coordinate_dict = defaultdict(list)
            dom_start=int(self.processing_flag)
            pre_dict=defaultdict(dict)
            ids_check_set=set()
            check_flag=0
            for i in range(dom_start,4):
                for record in SeqIO.parse(open('{3}/{2}/cry_extraction/{0}/{1}'.format('domains',self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0]+'_D{}_extracted.fasta'.format(i),self.quiery_dir,self.main_dir)),"fasta"):
                    dom_list=list()
                    if '|' in record.id:
                        name = record.id.split('|')[1].split('/')[0] + '|' +'|'.join('|'.join(record.description.split('] ')[1].split(' ')).split('_'))
                    else:
                        name = record.id.split('/')[0] + '|' +'|'.join('|'.join(record.description.split('] ')[1].split(' ')).split('_'))
                    if name in self.id_list:
                        if 'D'+str(i) in pre_dict[name].keys():
                            pre_dict[name]['D'+str(i)].extend(record.id.split('/')[1].split('-'))
                        else:
                            pre_dict[name]['D'+str(i)]=record.id.split('/')[1].split('-')
            for name_key in pre_dict:
                for i in range(dom_start,4):
                    if len(pre_dict[name_key]['D'+str(i)])==2:
                        self.coordinate_dict[name_key].extend(pre_dict[name_key]['D'+str(i)])
                    else:
                        d=0
                        ind=0
                        for j in range(0,len(pre_dict[name_key]['D'+str(i)]),2):
                            if int(pre_dict[name_key]['D'+str(i)][j+1]) - int(pre_dict[name_key]['D'+str(i)][j]) > d:
                                d=int(pre_dict[name_key]['D'+str(i)][j+1]) - int(pre_dict[name_key]['D'+str(i)][j])
                                ind=j
                        self.coordinate_dict[name_key].extend(pre_dict[name_key]['D'+str(i)][ind:ind+2])
            new_rec_list=list()   
            full_rec_list=list()     
            for record in SeqIO.parse(open(self.cry_quiery),"fasta"):
                if record.description != record.id:
                    name = record.id + '|' +'|'.join('|'.join(record.description.split(' ')[1:]).split('_'))
                else:
                    name= record.id + '|' +'|'.join('|'.join(record.description.split(' ')[0:]).split('_'))
                if name in self.id_list:
                     if record.id in ids_check_set:
                         check_flag=1
                     ids_check_set.add(record.id)
                     start = min([int(x) for x in self.coordinate_dict[name]])-1
                     stop = max([int(x) for x in self.coordinate_dict[name]])
                     full_rec_list.append(SeqRecord(Seq(str(record.seq),generic_protein),id=name.split('|')[0], description=" ".join(name.split('|')[1:])))
                     new_rec_list.append(SeqRecord(Seq(str(record.seq[start:stop]),generic_protein),id=name.split('|')[0], description=" ".join(name.split('|')[1:])))
            
            SeqIO.write(new_rec_list,'/{0}/{1}/cry_extraction/raw_processed_{2}.fasta'.format(self.main_dir,self.quiery_dir,self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0]), "fasta")
            SeqIO.write(full_rec_list,'/{0}/{1}/cry_extraction/raw_full_{2}.fasta'.format(self.main_dir,self.quiery_dir,self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0]), "fasta")
            if check_flag==1:
                print('Warning! Identical ids in queiry, uncertain mapping can occur')
            print('{} sequences recieved'.format(self.init_count))
            print('{} potential cry toxins found'.format(self.one_dom_count+self.two_dom_count+self.three_dom_count))
            print('{} toxins with one domain'.format(self.one_dom_count))
            print('{} toxins with two domains'.format(self.two_dom_count))
            print('{} toxins with three domains'.format(self.three_dom_count))

        if self.regime == 'fd':
            self.find_cry(self.cry_quiery)
            cmd_take_file = subprocess.call("cp {0}/{1}/cry_extraction/full_toxins/{2}_full_extracted.fasta {0}/{1}/cry_extraction/; mv {0}/{1}/cry_extraction/{2}_full_extracted.fasta {0}/{1}/cry_extraction/{2}.fasta; sed -i 's/[0-9]*|//g' {0}/{1}/cry_extraction/{2}.fasta; sed -i 's/\/[0-9]*-[0-9]*//g' {0}/{1}/cry_extraction/{2}.fasta; sed -i 's/ \[subseq from]\ / /g' {0}/{1}/cry_extraction/{2}.fasta".format(self.main_dir,self.quiery_dir,self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0]),shell=True)
            self.find_domains('{0}/{1}/cry_extraction/{2}.fasta'.format(self.main_dir,self.quiery_dir,self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0]))
            self.id_list = self.cry_3D_ids_extractor()
            self.three_dom_count = len(self.id_list)
            print("Performing cry toxins processing")
            self.coordinate_dict = defaultdict(list)
            pre_dict=defaultdict(dict)
            ids_check_set=set()
            check_flag=0
            self.first_filter_count = len(list(SeqIO.parse(open("{0}/{1}/cry_extraction/{2}.fasta".format(self.main_dir,self.quiery_dir,self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0])),"fasta")))
            dom_start=int(self.processing_flag)
            for i in range(dom_start,4):
                for record in SeqIO.parse(open('{3}/{2}/cry_extraction/{0}/{1}'.format('domains',self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0]+'_D{}_extracted.fasta'.format(i),self.quiery_dir,self.main_dir)),"fasta"):
                    if '|' in record.id:
                        name = record.id.split('|')[1].split('/')[0] + '|' +'|'.join('|'.join(record.description.split('] ')[1].split(' ')).split('_'))
                    else:
                        name = record.id.split('/')[0] + '|' +'|'.join('|'.join(record.description.split('] ')[1].split(' ')).split('_'))
                    if name in self.id_list:
                        if 'D'+str(i) in pre_dict[name].keys():
                            pre_dict[name]['D'+str(i)].extend(record.id.split('/')[1].split('-'))
                        else:
                            pre_dict[name]['D'+str(i)]=record.id.split('/')[1].split('-')
            for name_key in pre_dict:
                for i in range(dom_start,4):
                    if len(pre_dict[name_key]['D'+str(i)])==2:
                        self.coordinate_dict[name_key].extend(pre_dict[name_key]['D'+str(i)])
                    else:
                        d=0
                        ind=0
                        for j in range(0,len(pre_dict[name_key]['D'+str(i)]),2):
                            if int(pre_dict[name_key]['D'+str(i)][j+1]) - int(pre_dict[name_key]['D'+str(i)][j]) > d:
                                d=int(pre_dict[name_key]['D'+str(i)][j+1]) - int(pre_dict[name_key]['D'+str(i)][j])
                                ind=j
                        self.coordinate_dict[name_key].extend(pre_dict[name_key]['D'+str(i)][ind:ind+2])
            new_rec_list=list()   
            full_rec_list=list()    
            for record in SeqIO.parse(open("{0}/{1}/cry_extraction/{2}.fasta".format(self.main_dir,self.quiery_dir,self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0])),"fasta"):
                if record.description != record.id:
                    name = record.id + '|' +'|'.join('|'.join(record.description.split(' ')[1:]).split('_'))
                else:
                    name= record.id + '|' +'|'.join('|'.join(record.description.split(' ')[0:]).split('_'))
                if name in self.id_list:
                     if record.id in ids_check_set:
                         check_flag=1
                     ids_check_set.add(record.id)
                     start = min([int(x) for x in self.coordinate_dict[name]])-1
                     stop = max([int(x) for x in self.coordinate_dict[name]])
                     new_rec_list.append(SeqRecord(Seq(str(record.seq[start:stop].upper()),generic_protein),id=name.split('|')[0], description=" ".join(name.split('|')[1:])))
                     full_rec_list.append(SeqRecord(Seq(str(record.seq),generic_protein),id=name.split('|')[0], description=" ".join(name.split('|')[1:])))
            SeqIO.write(new_rec_list,'/{0}/{1}/cry_extraction/raw_processed_{2}.fasta'.format(self.main_dir,self.quiery_dir,self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0]), "fasta")
            SeqIO.write(full_rec_list,'/{0}/{1}/cry_extraction/raw_full_{2}.fasta'.format(self.main_dir,self.quiery_dir,self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0]), "fasta")
            if check_flag==1:
                print('Warning! Identical ids in queiry, uncertain mapping can occur')            
            print('{} sequences recieved'.format(self.init_count))
            print('{} sequences after first search'.format(self.first_filter_count))
            print('{} potential cry toxins found'.format(self.one_dom_count+self.two_dom_count+self.three_dom_count))
            print('{} toxins with one domain'.format(self.one_dom_count))
            print('{} toxins with two domains'.format(self.two_dom_count))
            print('{} toxins with three domains'.format(self.three_dom_count))


        with open("/{0}/{1}/cry_extraction/proteins_domain_mapping_{2}.bed".format(self.main_dir,self.quiery_dir,self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0]), 'w') as csv_file:
            my_writer = csv.writer(csv_file, delimiter='\t') 
            init_row = ['id','domain' ,'start', 'stop', 'description']
            my_writer.writerow(init_row)
            for key in self.coordinate_dict:
                for i in range(0,5,2):
                    if i==0:
                        row=[key.split('|')[0],'domain {}'.format(1),int(self.coordinate_dict[key][i])-int(self.coordinate_dict[key][0]),int(self.coordinate_dict[key][i+1])-int(self.coordinate_dict[key][0])-1, " ".join(key.split('|')[1:])] 
                    else:
                        if i==2:
                            dn=2
                        else:
                            dn=3
                        row=[key.split('|')[0],'domain {}'.format(dn),int(self.coordinate_dict[key][i])-int(self.coordinate_dict[key][0])-1,int(self.coordinate_dict[key][i+1])-int(self.coordinate_dict[key][0])-1, " ".join(key.split('|')[1:])] 
                    my_writer.writerow(row)

        with open("/{0}/{1}/cry_extraction/coordinate_matches_{2}.txt".format(self.main_dir,self.quiery_dir,self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0]), 'w') as csv_file:
            my_writer = csv.writer(csv_file, delimiter='\t') 
            for key in self.coordinate_dict:
                my_writer.writerow([key]+self.coordinate_dict[key])

    def annotate_raw_output(self):
        new_records = list()
        self.new_ids=dict()
        un_count=0
        total_count=0
        print("Annotating raw output with diamond")  
        cmd_pre_dia = subprocess.call('cd /{0}/{1}/cry_extraction; cp {2}/data/diamond_data/cry_nomenclature.dmnd .; touch diamond.log'.format(self.main_dir,self.quiery_dir,self.home_dir), shell=True) 
        cmd_dia = subprocess.call('cd /{1}/{2}/cry_extraction;{3}/diamond blastp -d cry_nomenclature -q raw_full_{0}.fasta -o diamond_matches_{0}.txt --al aligned_{0}.fa --un unaligned_{0}.fa --max-target-seqs 1 --log --verbose 2>> diamond.log; rm cry_nomenclature.dmnd; mv *.log logs/; mv *gned* logs/'.format(self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0],self.main_dir,self.quiery_dir,self.home_dir+'/include'), shell=True)
        with open("/{0}/{1}/cry_extraction/diamond_matches_{2}.txt".format(self.main_dir,self.quiery_dir,self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0]), 'r') as csv_file:
            my_reader = csv.reader(csv_file, delimiter='\t') 
            for row in my_reader:
                total_count+=1
                if float(row[2])<100.0:
                    self.new_ids[row[0]]=row[1]+'|'+ str(row[2])
                    un_count+=1
        for init_rec in SeqIO.parse(open('/{0}/{1}/cry_extraction/raw_processed_{2}.fasta'.format(self.main_dir,self.quiery_dir,self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0])),"fasta"):
           if init_rec.id in self.new_ids.keys():
               new_records.append(SeqRecord(Seq(str(init_rec.seq),generic_protein),id=self.new_ids[init_rec.id].split('|')[0]+'(' + self.new_ids[init_rec.id].split('|')[1]+')'+'_'+init_rec.description.split()[0],description=init_rec.description))
        print('{} sequences matched with database'.format(total_count))
        print('{} toxins different from database found'.format(un_count))
        SeqIO.write(new_records,'/{0}/{1}/cry_extraction/unique_{2}.fasta'.format(self.main_dir,self.quiery_dir,self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0]), "fasta")
        if self.regime == 'fd':
            cmd_clean_up = subprocess.call("mv {0}/{1}/cry_extraction/{2}.fasta {0}/{1}/cry_extraction/first_search_{2}.fasta".format(self.main_dir,self.quiery_dir,self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0]),shell=True)
        if not self.nucl_type and not self.annot_flag:
            cmd_clean_up = subprocess.call('cd $PWD/{0}/cry_extraction; mv *coordinate_matches* logs/;mv *diamond_matches* logs/'.format(self.quiery_dir), shell=True)
       
    def make_summary_table(self):
        print('Searching for the metadata')
        Entrez.email = "{}".format(self.email)
        summary_dict=defaultdict(dict)
        with open("/{0}/{1}/cry_extraction/diamond_matches_{2}.txt".format(self.main_dir,self.quiery_dir,self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0]), 'r') as csv_file:
            my_reader = csv.reader(csv_file, delimiter='\t') 
            for row in my_reader:
                summary_dict[row[0]]=defaultdict(list)
                summary_dict[row[0]]['init']=row[1:3]
        making_smb = subprocess.call("touch /{0}/{1}/cry_extraction/annotation_table_{2}.tsv".format(self.main_dir,self.quiery_dir,self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0]), shell = True)
        for init_rec in SeqIO.parse(open('/{0}/{1}/cry_extraction/raw_processed_{2}.fasta'.format(self.main_dir,self.quiery_dir,self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0])),"fasta"):
           if init_rec.id in summary_dict.keys():
               summary_dict[init_rec.id]['init'].append(init_rec.description)
        init_row = ['protein_id', 'initial_description', 'top_cry_hit', 'cry_identity', 'source', 'nucl_accession', 'start','stop', 'strand','ipg_prot_id','ipg_prot_name', 'organism', 'strain','assembly\n']
        f = open("/{0}/{1}/cry_extraction/annotation_table_{2}.tsv".format(self.main_dir,self.quiery_dir,self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0]),"w")
        f.write('\t'.join(init_row))
        f.close() 
        for key in summary_dict:
            try:
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
                    f = open("/{0}/{1}/cry_extraction/annotation_table_{2}.tsv".format(self.main_dir,self.quiery_dir,self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0]),"a+")
                    f.write('\t'.join(row) + '\n') 
                    f.close()
                time.sleep(3)
            except:
                print('Warning! Bad link for!', key)
        if not self.nucl_type:
            cmd_clean_up = subprocess.call('cd $PWD/{0}/cry_extraction; mv *coordinate_matches* logs/;mv *diamond_matches* logs/'.format(self.quiery_dir), shell=True)

    def upload_nucl(self):
        print('Uploading nucleotide sequences')
        Entrez.email = "{}".format(self.email)
        keys_for_nucl = defaultdict(list)
        domain_coord_dict = defaultdict(list)
        with open("/{0}/{1}/cry_extraction/coordinate_matches_{2}.txt".format(self.main_dir,self.quiery_dir,self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0]), 'r') as csv_file:
            my_reader = csv.reader(csv_file, delimiter='\t')
            for row in my_reader:
                 domain_coord_dict[row[0]]=row[1:]
        init_count=0
        with open("/{0}/{1}/cry_extraction/annotation_table_{2}.tsv".format(self.main_dir,self.quiery_dir,self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0]), 'r') as csvfile:
             my_reader = csv.reader(csvfile, delimiter='\t') 
             for row in my_reader:
                 if row[0]!= '--' and row[0]!='protein_id':
                     init_count+=1
                     keys_for_nucl[row[0]].append('|'.join((row[1]).split(' ')))
                     keys_for_nucl[row[0]].extend([row[5], row[6], row[7], row[8], row[2], row[3]])
        if keys_for_nucl[list(keys_for_nucl.keys())[0]][0] not in domain_coord_dict:
            for key in keys_for_nucl:
                keys_for_nucl[key][0]=keys_for_nucl[key][0]+'|'+keys_for_nucl[key][0]
        f_nuc_recs = []
        p_nuc_recs = []
        if self.nucl_type == 'fn':
            for key in keys_for_nucl:
                try:
                    handle = Entrez.efetch(db="nucleotide",rettype='fasta',retmode='text', id = keys_for_nucl[key][1])
                    fasta_rec = SeqIO.read(handle, "fasta")
                    handle.close()
                    if keys_for_nucl[key][4] == '+':
                        f_nuc_recs.append(SeqRecord(Seq(str(fasta_rec.seq[int(keys_for_nucl[key][2])-1:int(keys_for_nucl[key][3])]),generic_dna),id=keys_for_nucl[key][5]+'('+keys_for_nucl[key][6]+')'+'_'+keys_for_nucl[key][0].split('|')[0], description=" ".join(keys_for_nucl[key][0].split('|'))))
                    elif keys_for_nucl[key][4] == '-':
                        f_nuc_recs.append(SeqRecord(Seq(str(fasta_rec.seq[int(keys_for_nucl[key][2])-1:int(keys_for_nucl[key][3])].reverse_complement()),generic_dna),id=keys_for_nucl[key][5]+'('+keys_for_nucl[key][6]+')'+'_'+keys_for_nucl[key][0].split('|')[0], description=" ".join(keys_for_nucl[key][0].split('|'))))
                except:
                    print('Warning! Download error for ', keys_for_nucl[key][2], '(' + keys_for_nucl[key][1]+ ')')

            SeqIO.write(f_nuc_recs,"/{0}/{1}/cry_extraction/{2}_full_nucl.fna".format(self.main_dir,self.quiery_dir,self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0]), "fasta")

        elif self.nucl_type == 'pn':
            for key in keys_for_nucl:
                try:
                    handle = Entrez.efetch(db="nucleotide",rettype='fasta',retmode='text', id = keys_for_nucl[key][1])
                    fasta_rec = SeqIO.read(handle, "fasta")
                    handle.close()
                    if keys_for_nucl[key][4] == '+':
                        dum_seq = str(fasta_rec.seq[int(keys_for_nucl[key][2])-1:int(keys_for_nucl[key][3])])
                        p_nuc_recs.append(SeqRecord(Seq(dum_seq[(min([int(x) for x in domain_coord_dict[keys_for_nucl[key][0]]])-1)*3:max([int(x) for x in domain_coord_dict[keys_for_nucl[key][0]]])*3],generic_dna),id=keys_for_nucl[key][5]+'('+keys_for_nucl[key][6]+')'+'_'+keys_for_nucl[key][0].split('|')[0], description=" ".join(keys_for_nucl[key][0].split('|'))))
                    elif keys_for_nucl[key][4] == '-':
                        dum_seq = str(Seq(str(fasta_rec.seq[int(keys_for_nucl[key][2])-1:int(keys_for_nucl[key][3])].reverse_complement()),generic_dna))
                        p_nuc_recs.append(SeqRecord(Seq(dum_seq[(min([int(x) for x in domain_coord_dict[keys_for_nucl[key][0]]])-1)*3:max([int(x) for x in domain_coord_dict[keys_for_nucl[key][0]]])*3],generic_dna),id=keys_for_nucl[key][5]+'('+keys_for_nucl[key][6]+')'+'_'+keys_for_nucl[key][0].split('|')[0], description=" ".join(keys_for_nucl[key][0].split('|'))))
                except:
                    print('Warning! Download error for ', keys_for_nucl[key][2], '(' + keys_for_nucl[key][1]+ ')')

        elif self.nucl_type == 'an':
            for key in keys_for_nucl:
                try:
                    handle = Entrez.efetch(db="nucleotide",rettype='fasta',retmode='text', id = keys_for_nucl[key][1])
                    fasta_rec = SeqIO.read(handle, "fasta")
                    handle.close()
                    if keys_for_nucl[key][4] == '+':
                        dum_seq = str(fasta_rec.seq[int(keys_for_nucl[key][2])-1:int(keys_for_nucl[key][3])])
                        p_nuc_recs.append(SeqRecord(Seq(dum_seq[(min([int(x) for x in domain_coord_dict[keys_for_nucl[key][0]]])-1)*3:max([int(x) for x in domain_coord_dict[keys_for_nucl[key][0]]])*3],generic_dna),id=keys_for_nucl[key][5]+'('+keys_for_nucl[key][6]+')'+'_'+keys_for_nucl[key][0].split('|')[0], description=" ".join(keys_for_nucl[key][0].split('|'))))
                        f_nuc_recs.append(SeqRecord(Seq(str(fasta_rec.seq[int(keys_for_nucl[key][2])-1:int(keys_for_nucl[key][3])]),generic_dna),id=keys_for_nucl[key][5]+'('+keys_for_nucl[key][6]+')'+'_'+keys_for_nucl[key][0].split('|')[0], description=" ".join(keys_for_nucl[key][0].split('|'))))
                    elif keys_for_nucl[key][4] == '-':
                        dum_seq = str(Seq(str(fasta_rec.seq[int(keys_for_nucl[key][2])-1:int(keys_for_nucl[key][3])].reverse_complement()),generic_dna))
                        p_nuc_recs.append(SeqRecord(Seq(dum_seq[(min([int(x) for x in domain_coord_dict[keys_for_nucl[key][0]]])-1)*3:max([int(x) for x in domain_coord_dict[keys_for_nucl[key][0]]])*3],generic_dna),id=keys_for_nucl[key][5]+'('+keys_for_nucl[key][6]+')'+'_'+keys_for_nucl[key][0].split('|')[0], description=" ".join(keys_for_nucl[key][0].split('|'))))
                        f_nuc_recs.append(SeqRecord(Seq(str(fasta_rec.seq[int(keys_for_nucl[key][2])-1:int(keys_for_nucl[key][3])].reverse_complement()),generic_dna),id=keys_for_nucl[key][5]+'('+keys_for_nucl[key][6]+')'+'_'+keys_for_nucl[key][0].split('|')[0], description=" ".join(keys_for_nucl[key][0].split('|'))))
                except:
                    print('Warning! Download error for ', keys_for_nucl[key][2], '(' + keys_for_nucl[key][1]+ ')')

            SeqIO.write(p_nuc_recs,"/{0}/{1}/cry_extraction/{2}_processed_nucl.fna".format(self.main_dir,self.quiery_dir,self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0]), "fasta")
        print('{} nucleotide sequences downloaded'.format(max([len(p_nuc_recs), len(f_nuc_recs)])))


    def map_nucl(self):
        keys_for_nucl = defaultdict(list)
        domain_coord_dict = defaultdict(list)       
        with open("/{0}/{1}/cry_extraction/coordinate_matches_{2}.txt".format(self.main_dir,self.quiery_dir,self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0]), 'r') as csv_file:
            my_reader = csv.reader(csv_file, delimiter='\t')
            for row in my_reader:
                 domain_coord_dict[row[0]]=row[1:]
        init_count=0
        with open("/{0}/{1}/cry_extraction/annotation_table_{2}.tsv".format(self.main_dir,self.quiery_dir,self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0]), 'r') as csvfile:
            my_reader = csv.reader(csvfile, delimiter='\t')
            for row in my_reader:
                if row[0]!= '--' and row[0]!='protein_id':                          
                    init_count+=1
                    keys_for_nucl[row[0]].append('|'.join((row[1]).split(' ')))
                    keys_for_nucl[row[0]].extend([row[5], row[6], row[7], row[8], row[2], row[3]])                              
        if keys_for_nucl[list(keys_for_nucl.keys())[0]][0] not in domain_coord_dict:
            for key in keys_for_nucl:
                keys_for_nucl[key][0]=keys_for_nucl[key][0]+'|'+keys_for_nucl[key][0]

        with open("/{0}/{1}/cry_extraction/nucl_domain_mapping_{2}.bed".format(self.main_dir,self.quiery_dir,self.cry_quiery.split('/')[len(self.cry_quiery.split('/'))-1].split('.')[0]), 'w') as csvfile:
             my_writer = csv.writer(csvfile, delimiter='\t') 
             init_row = ['id','domain', 'start', 'stop', 'description']
             my_writer.writerow(init_row)
             for key in keys_for_nucl:
                 try:
                     for i in range(0,5,2):
                         if i==0:
                            row=[keys_for_nucl[key][5]+'('+keys_for_nucl[key][6]+')'+'_'+keys_for_nucl[key][0].split('|')[0],'domain {}'.format(1),(int(domain_coord_dict[keys_for_nucl[key][0]][i])-1)*3-(int(domain_coord_dict[keys_for_nucl[key][0]][0])-1)*3,int(domain_coord_dict[keys_for_nucl[key][0]][i+1])*3-(int(domain_coord_dict[keys_for_nucl[key][0]][0])-1)*3-1, " ".join(keys_for_nucl[key][0].split('|'))]  
                         else:
                            if i==2:
                                dn=2
                            else:
                                dn=3
                            row=[keys_for_nucl[key][5]+'('+keys_for_nucl[key][6]+')'+'_'+keys_for_nucl[key][0].split('|')[0],'domain {}'.format(dn),(int(domain_coord_dict[keys_for_nucl[key][0]][i])-1)*3-(int(domain_coord_dict[keys_for_nucl[key][0]][0])-1)*3-1,int(domain_coord_dict[keys_for_nucl[key][0]][i+1])*3-(int(domain_coord_dict[keys_for_nucl[key][0]][0])-1)*3-1, " ".join(keys_for_nucl[key][0].split('|'))]      
                         my_writer.writerow(row)
                 except:
                     pass
        cmd_clean_up = subprocess.call('cd $PWD/{0}/cry_extraction; mv *coordinate_matches* logs/;mv *diamond_matches* logs/'.format(self.quiery_dir), shell=True)

    def launch_racer(self,kmer):
      print('Searching sequences from gfa file')
      cmd_race = subprocess.call('{0}/pathracer {6}/data/models/D1.hmm {1} {5} --output $PWD/{3}/cry_extraction/pathracer_output -t {2}  > /dev/null'.format(self.home_dir+'/include',self.cry_quiery,self.hm_threads,self.quiery_dir,'pathracer',kmer,self.home_dir), shell=True) 
      cmd_merge = subprocess.call('cd $PWD/{3}/cry_extraction/pathracer_output; cat *seqs.fa >> mearged.fasta; cp pathracer.log ../logs/'.format(self.home_dir+'/include',self.cry_quiery,self.hm_threads,self.quiery_dir,'pathracer',kmer,self.home_dir), shell=True) 


    def use_spades(self):
        print('Building assembly graph')
        if self.meta_flag:
            cmd_spades = subprocess.call('{0}/SPAdes-3.13.1-Linux/bin/spades.py --meta -1 {1} -2 {2} -o $PWD/{3}/cry_extraction/assembly -t {4} > /dev/null'.format(self.home_dir+'/include', self.forw, self.rev, self.quiery_dir, self.hm_threads), shell=True) 
            cmd_merge = subprocess.call('cd $PWD/{0}/cry_extraction/assembly; cp *.log ../logs/'.format(self.quiery_dir), shell=True)
        else:
            cmd_spades = subprocess.call('{0}/SPAdes-3.13.1-Linux/bin/spades.py -1 {1} -2 {2} -o $PWD/{3}/cry_extraction/assembly -t {4} > /dev/null'.format(self.home_dir+'/include', self.forw, self.rev,self.quiery_dir, self.hm_threads), shell=True) 
            cmd_merge = subprocess.call('cd $PWD/{0}/cry_extraction/assembly; cp *.log ../logs/'.format(self.quiery_dir), shell=True)       

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='cry_processor')
    parser.add_argument('-fi', help='Enter input file: fasta file or gfa file', metavar='Fasta file/GFA file',
                        type=str, default=None)
    parser.add_argument('-hm', help='Path to hmmer if hmmer is installed locally', metavar='Hmmer_directory',
                        type=str, default='')
    parser.add_argument('-pr', help='Processig type: 1 for extracting all domains, 2 for extrating 2-3 domains only', metavar='Int',type=str, default=1)
    parser.add_argument('-th', help='Number of threads for hmmer and pathracer', metavar='Int',type=str, default=8)
    parser.add_argument('-ma', help='e-mail address for NCBI annotation', metavar='Str',type=str, default='')
    parser.add_argument('-od', help='Output directory', metavar='Str',type=str, required=True)
    parser.add_argument('-r', help='Pipeline type: do - domain only search with subsequent unioning; fd - searching for potential cry-toxins with subsequent processing', metavar='Str',type=str, default='do')
    parser.add_argument('--annotate', '-a',action='store_true',help='make final output annotation with ipg')
    parser.add_argument('-nu', help='Uploading nucleotide records: fn - uploading full sequences, pn - uploading processed subsequences, an - both processed and unprocessed', metavar='Nucl_uploading_type', type=str, default='')
    parser.add_argument('--path_racer', '-pa',action='store_true',help='Searching for cry toxins from gfa file with pathracer')
    parser.add_argument('-fo', help='Forward illumina reads', metavar='Fasta file',type=str, default=None)
    parser.add_argument('-re', help='Reverse illumina reads', metavar='Fasta file',type=str, default=None)
    parser.add_argument('--meta','-m', action='store_true',help='Metagenomic regime for spades')
    parser.add_argument('-k', help='k-mer size for pathracer', metavar='Int',type=str, default=21)
    parser.set_defaults(feature=True)
    args = parser.parse_args()
    od,fi,hm,pr,th, ma, r, a, nu, mra,k,fr,rr,meta = args.od, args.fi, args.hm, args.pr,args.th, args.ma, args.r, args.annotate, args.nu,args.path_racer,args.k,args.fo,args.re,args.meta
    cr = CryProcessor(od, fi, hm,pr, th, ma, r, nu,a,k,meta,fr,rr)
    if not mra and not fr:
        cr.cry_digestor()
        cr.annotate_raw_output()
        if a: 
            cr.make_summary_table()
        if nu:
            cr.upload_nucl()
            cr.map_nucl()
    elif mra and not fr:
        cr.launch_racer(k)
        fi = od + '/cry_extraction/pathracer_output/mearged.fasta'
        cr = CryProcessor(od, fi, hm,pr, th, ma, r, nu,a,k,meta,fr,rr)
        cr.cry_digestor()
        cr.annotate_raw_output()
        if a: 
            cr.make_summary_table()
        if nu:
            cr.upload_nucl()
            cr.map_nucl()
    elif fr:
        cr.use_spades()
        fi = od + '/cry_extraction/assembly/assembly_graph_with_scaffolds.gfa'
        cr = CryProcessor(od, fi, hm,pr, th, ma, r, nu,a,k,meta,fr,rr)
        cr.launch_racer(k)
        fi = od + '/cry_extraction/pathracer_output/mearged.fasta'
        cr = CryProcessor(od, fi, hm,pr, th, ma, r, nu,a,k,meta,fr,rr)
        cr.cry_digestor()
        cr.annotate_raw_output()
        if a: 
            cr.make_summary_table()
        if nu:
            cr.upload_nucl()
            cr.map_nucl()
        
   

