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
import os.path
import sys
import re
import logging
import copy
import datetime

 
class CryProcessor:
    def __init__(self, 
                 query_dir, 
                 cry_query, 
                 hmmer_dir, 
                 processing_flag, 
                 hm_threads, 
                 email, 
                 mode, 
                 nucl_type, 
                 annot_flag,
                 kmer_size,
                 meta_flag,
                 forw,
                 rev,
                 silent_mode):
        #initialize the logger
        self.logger = logging.getLogger("cry_processor")
        self.logger.setLevel(logging.INFO)
        fh = logging.FileHandler("cry_processor.log")
        #set the formatter for the logger
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        fh.setFormatter(formatter)
        self.logger.addHandler(fh)
        self.logger.info("Initializing...")
        self.cry_query = cry_query
        self.hmmer_dir = hmmer_dir
        self.processing_flag = processing_flag
        self.hm_threads = hm_threads
        self.query_dir = query_dir
        if self.cry_query:
            #calculate the number of sequences (only if the fasta file is present)
            if 'gfa' not in self.cry_query:
                self.init_count = re.sub("b",'',
                                  re.sub("\'",'', 
                                  str(subprocess.check_output("grep '>' {} | wc -l".format(self.cry_query), 
                                  shell =True).strip())))
                self.one_dom_count = 0
                self.two_dom_count = 0
                self.three_dom_count = 0
        self.email = email
        self.mode = mode
        self.nucl_type = nucl_type
        self.annot_flag = annot_flag
        self.kmer_size = kmer_size
        self.meta_flag = meta_flag
        self.forw = forw
        self.rev = rev
        self.racer_flag = 0
        self.hmmer_result_flag = 0 
        self.silent_mode = silent_mode
        if os.path.exists(os.path.realpath(self.query_dir)):
            self.query_dir = self.query_dir + "_"+ str(datetime.datetime.now()).split('.')[0].replace(' ', '_').replace(':', '_')
        print(self.query_dir)
        #creating output directories
        cmd_init = subprocess.call('if [ ! -d {0} ]; \
                                     then mkdir {0}; \
                                                      fi; \
                if [ ! -d {0}/logs ]; \
                 then mkdir {0}/logs; \
                                                      fi'.format(os.path.realpath(
                                                                 self.query_dir)), 
                                                                 shell=True)


    def run_hmmer(self,query,out_index,model_type,dir_flag,log,hm_threads,query_dir):
        """
        Runs hmmsearch on query
        Args
        =====
        query: fasta file for performing search
        out_index: postfix for .sto output file
        model_type: which hmm-model to use
        dir_flag: output directory for search
        log: prefix for log-file
        hm_threads: number of threads to use
        query_dir: directory where query file is located
        """
        #if self.hmmer_dir is not specified use hmmsearch from the include directory
        if self.hmmer_dir:
            #check the output directory structure
            cmd_make_dirs = subprocess.call('cd {0}; \
                                            if [ ! -d {1} ]; \
                                            then mkdir {1}; \
                                            fi'.format(os.path.join(
                                            os.path.realpath(
                                            query_dir)),
                                            dir_flag), 
                                            shell=True)
            #perform hmmsearch 

            cmd_search = subprocess.call('{0} \
                                          -E 50 --nonull2 --max --incE 250 --domZ 0.1 --cpu {1} \
                                          -A {2} \
                                          {3} \
                                          {4} >> {5}'.format(os.path.join(
                                                             os.path.realpath(
                                                             self.hmmer_dir),
                                                             "binaries",
                                                             "hmmsearch"),
                                                              hm_threads,
                                                             os.path.join(
                                                             os.path.realpath(
                                                             self.query_dir),
                                                              dir_flag,
                                                             self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0]+ 
                                                             out_index),
                                                             os.path.join(
                                                             os.path.dirname
                                                             (__file__),
                                                              'data',
                                                              'models',
                                                              model_type),
                                                              os.path.realpath(
                                                              query), 
                                                              os.path.join(
                                                              os.path.realpath(
                                                              self.query_dir),
                                                              "logs", 
                                                              log+'.log')), 
                                                              shell=True) 
            #converting to fasta files

            cmd_fasta = subprocess.call('if [ -s {0} ];\
                                         then {2} fasta {0} > {1};\
                                         fi'.format(os.path.join(
                                              os.path.realpath(
                                              self.query_dir),
                                              dir_flag,
                                              self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0]+ 
                                              out_index), 
                                              os.path.join(
                                              os.path.realpath(
                                              self.query_dir),
                                              dir_flag,
                                              self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0]+ 
                                              out_index.replace('sto','fasta')),
                                              os.path.join(os.path.realpath(
                                              self.hmmer_dir),
                                              "binaries","esl-reformat")), 
                                              shell=True)
# -E 50 --nonull2 --max --incE 250 --domZ 0.1 

        else:
            cmd_make_dirs = subprocess.call('cd {0}; \
                                            if [ ! -d {1} ]; \
                                            then mkdir {1}; \
                                            fi'.format(os.path.join(
                                            os.path.realpath(
                                            query_dir)),
                                            dir_flag), 
                                            shell=True)

            cmd_search = subprocess.call('{0} \
                                          -E 50 --nonull2 --max --incE 250 --domZ 0.1 --cpu {1} \
                                          -A {2} \
                                          {3} \
                                          {4} >> {5}'.format(
                                          os.path.join(
                                          os.path.dirname(
                                          __file__),
                                          "include",
                                          "hmmsearch"),
                                          hm_threads,
                                          os.path.join(
                                          os.path.realpath(
                                          self.query_dir),
                                          dir_flag,
                                          self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0]+ 
                                          out_index),
                                          os.path.join(
                                          os.path.dirname(
                                          __file__),
                                          'data',
                                          'models',
                                          model_type),
                                          os.path.realpath(
                                          query),
                                          os.path.join(
                                          os.path.realpath(
                                          self.query_dir),
                                          "logs", 
                                          log+'.log')), 
                                          shell=True) 

            cmd_fasta = subprocess.call('if [ -s {0} ];\
                                         then {2} fasta {0} > {1};\
                                         fi'.format(os.path.join(
                                         os.path.realpath(
                                         self.query_dir),
                                         dir_flag,
                                         self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0]+ 
                                         out_index), 
                                         os.path.join(
                                         os.path.realpath(
                                         self.query_dir),
                                         dir_flag,
                                         self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0]+ 
                                         out_index.replace('sto','fasta')),
                                         os.path.join(
                                         os.path.dirname(
                                         __file__),
                                         "include",
                                         "esl-reformat")), 
                                         shell=True) 

    def find_cry(self, qiuery):
        """
        Launches run_hmmer function on full unprocessed toxins (Complete.hmm model in data/models)
        Args
        =====
        query: fasta file for performing search
        """
        if not self.silent_mode:
            print('Searching for the unprocessed Cry toxins')
        self.logger.info('Searching for the unprocessed Cry toxins')
        self.run_hmmer(str(qiuery),
                       '_full_extracted.sto',
                       'Complete.hmm',
                       'full_toxins', 
                       'full_extraction',
                       self.hm_threads, 
                       self.query_dir)
    
    def find_domains(self, qiuery):
        """
        Launches run_hmmer function on full domains of cry-toxins (D1.hmm, D2.hmm and D3.hmm models in data/models)
        Args
        =====
        query: fasta file for performing search
        """
        # loop over domain indicies
        for i in range(1,4):
            if not self.silent_mode:
                print('Searching for the {} domain  of the Cry toxins'.format(i))
            self.logger.info('Searching for the {} domain of the Cry toxins'.format(i))
            self.run_hmmer(str(qiuery),
                           '_D{}_extracted.sto'.format(i),
                           'D{}.hmm'.format(i),
                           'domains',
                           'domains_extraction', 
                            self.hm_threads, 
                            self.query_dir)
        self.run_hmmer(str(qiuery),
                           '_cry31.sto',
                           'cry31.hmm',
                           'domains',
                           'domains_extraction', 
                            self.hm_threads, 
                            self.query_dir)
        self.run_hmmer(str(qiuery),
                           '_cry_58.sto',
                           'cry_58.hmm',
                           'domains',
                           'domains_extraction', 
                            self.hm_threads, 
                            self.query_dir)
        self.run_hmmer(str(qiuery),
                           '_extra.sto',
                           'extra.hmm',
                           'domains',
                           'domains_extraction', 
                            self.hm_threads, 
                            self.query_dir)
        self.run_hmmer(str(qiuery),
                           '_Endotoxin_M.sto',
                           'Endotoxin_M.hmm',
                           'domains',
                           'domains_extraction', 
                            self.hm_threads, 
                            self.query_dir)
        self.run_hmmer(str(qiuery),
                           '_Endotoxin_mid.sto',
                           'Endotoxin_mid.hmm',
                           'domains',
                           'domains_extraction', 
                            self.hm_threads, 
                            self.query_dir)

        self.run_hmmer(str(qiuery),
                           '_long_32_11.sto',
                           'long_32_11.hmm',
                           'domains',
                           'domains_extraction', 
                            self.hm_threads, 
                            self.query_dir)
        self.run_hmmer(str(qiuery),
                           '_extra_3.sto',
                           'extra_3.hmm',
                           'domains',
                           'domains_extraction', 
                            self.hm_threads, 
                            self.query_dir)

        mearging_cmd = subprocess.call("cd {0}; ls | grep '*fasta';\
                                         cat $(ls | grep -E '*(_D2_extracted|_cry31|_long_32_11|_Endotoxin_mid|_Endotoxin_M|_extra|_cry_58).fasta') > tmp_2.fasta; \
                                         mv tmp_2.fasta {1};   \
                                         fi".format(os.path.join(
                                         os.path.realpath(
                                         self.query_dir),
                                         'domains'),
                                         self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0]+ 
                                         '_D2_extracted.fasta'), 
                                         shell=True)

        mearging_cmd = subprocess.call("cd {0};\
                                         cat $(ls | grep -E '*(_D3_extracted|_extra_3).fasta') > tmp_3.fasta; \
                                         mv tmp_3.fasta {1};   \
                                         fi".format(os.path.join(
                                         os.path.realpath(
                                         self.query_dir),
                                         'domains'),
                                         self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0]+ 
                                         '_D3_extracted.fasta'), 
                                         shell=True)
        exist_check_flag = re.sub("b",'',
                                  re.sub("\'",'', 
                                  str(subprocess.check_output("cd {0}; \
                                  if [ ! -s {1} ] || [ ! -s {2} ] || [ ! -s {3} ];\
                                  then echo 'no'; \
                                  fi".format(os.path.join(
                                  os.path.realpath(
                                  self.query_dir),
                                  'domains'),
                                  self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0]+ 
                                         '_D1_extracted.fasta',
                                  self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0]+ 
                                         '_D2_extracted.fasta',
                                  self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0]+ 
                                         '_D3_extracted.fasta'),
                                  shell =True).strip())))
        if exist_check_flag == 'no':
            self.hmmer_result_flag = 1

        

    def cry_3D_ids_extractor(self):
        """
        Extracts ids from query for sequences possesing 3 domains of cry toxins
        """
        if not self.silent_mode:
            print('Exctracting the Cry toxins with 3 domains')
        self.logger.info('Exctracting the Cry toxins with 3 domains')
        #dictionary of ids is used for the further extraction and processing
        self.dom_dict=defaultdict(list)
        for i in range(1,4):
            #open all the fasta-files after hmmsearch and save them to the dictionary
            for record in SeqIO.parse(os.path.join(
                                      os.path.realpath(
                                      self.query_dir),
                                      'domains',
                                      self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0]+
                                      '_D{}_extracted.fasta'.format(i)),
                                      "fasta"):
                #transform the sequence name by combining its id and description
                if '|' in record.id:
                    name = record.id.split('|')[1].split('/')[0] + '|' + \
                           '|'.join('|'.join(record.description.split('] ')[1].split(' ')).split('_'))
                else:
                    name = record.id.split('/')[0] + '|' + \
                           '|'.join('|'.join(record.description.split('] ')[1].split(' ')).split('_'))

                if 'D{}'.format(i) not in self.dom_dict[name]:
                    self.dom_dict[name].append('D{}'.format(i))
        self.one_dom_count = len([key for key in self.dom_dict if len(set(self.dom_dict[key]))==1])
        self.two_dom_count = len([key for key in self.dom_dict if len(set(self.dom_dict[key]))==2])
        #returning ids for the sequences with 3 domains (possesing D1,D2,D3 as dictionary items)
        return([key for key in self.dom_dict if len(set(self.dom_dict[key]))>=3])

    def cry_digestor(self):
        """
        Performs processing toxins with 3 domains: cuts left and right flanking sequences 
        with saving only domains and linkers between domains
        """
        if self.mode == 'do':
            #domain_only mode: perform hmmsearch only on the domain hmm-models
            self.find_domains(self.cry_query)
            if self.hmmer_result_flag ==1:
                if not self.silent_mode:
                    print('No toxins found')
                self.logger.info('No toxins found')
            #obtaining an id list with the cry_3D_ids_extractor function
            else:
                self.id_list = self.cry_3D_ids_extractor()
                final_id_list=list()
                if not self.silent_mode:
                    print("Performing the Cry toxins processing")
                self.logger.info('Performing the Cry toxins processing')
                self.coordinate_dict = defaultdict(list)
                #a starting domain (depends on the -pr flag: you can extract all 3 domains by default or extract 2 and 3 domains only)
                dom_start=int(self.processing_flag)
                #dictioinary for checking: is needed because sometimes sequences have more than only one match after hmmersearch for one domain (maybe duplication regions?)
                pre_dict=defaultdict(dict)
                ids_check_set=set()
                #check_flag for identical ids in the query
                check_flag=0
                for i in range(dom_start,4):
                    #iterate over fasta fles with domains
                    for record in SeqIO.parse(os.path.join(
                                              os.path.realpath(
                                              self.query_dir),
                                              'domains',
                                              self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0]+
                                              '_D{}_extracted.fasta'.format(i)),
                                              "fasta"):
                        dom_list=list()
                        if '|' in record.id:
                            name = record.id.split('|')[1].split('/')[0] + '|' + \
                               '|'.join('|'.join(record.description.split('] ')[1].split(' ')).split('_'))
                        else:
                            name = record.id.split('/')[0] + '|' + \
                               '|'.join('|'.join(record.description.split('] ')[1].split(' ')).split('_'))
                        if name in self.id_list:
                            #add all the hmm-matches to pre_dict for the particular sequence
                            if 'D'+str(i) in pre_dict[name].keys():
                            #add start and stop coordinates, based on hmmsearch for the each domain
                                pre_dict[name]['D'+str(i)].extend(record.id.split('/')[1].split('-'))
                            else:
                                pre_dict[name]['D'+str(i)]=record.id.split('/')[1].split('-')
                for name_key in pre_dict:
                #processing pre_dict to get the final domains mapping
                    for i in range(dom_start,4):
                        if len(pre_dict[name_key]['D'+str(i)])==2:
                            #processing pre_dict to get the final domains mappings
                            #if it contains only two items (start and stop), we will use it
                            self.coordinate_dict[name_key].extend(pre_dict[name_key]['D'+str(i)])
                        else:
                            #choose the longest domain sequence if multiple mathces occur
                            d=0
                            ind=0
                            #get an index of the longest match after looping over all the mathes
                            for j in range(0,len(pre_dict[name_key]['D'+str(i)]),2):
                                if int(pre_dict[name_key]['D'+str(i)][j+1]) - int(pre_dict[name_key]['D'+str(i)][j]) > d:
                                    d=int(pre_dict[name_key]['D'+str(i)][j+1]) - int(pre_dict[name_key]['D'+str(i)][j])
                                    ind=j
                            self.coordinate_dict[name_key].extend(pre_dict[name_key]['D'+str(i)][ind:ind+2])
                for key in self.coordinate_dict:
                    if int(self.coordinate_dict[key][1]) - int(self.coordinate_dict[key][0])>75 and int(self.coordinate_dict[key][3]) - int(self.coordinate_dict[key][2])>75 and int(self.coordinate_dict[key][5]) - int(self.coordinate_dict[key][4])>75:
                        if int(self.coordinate_dict[key][4]) <= int(self.coordinate_dict[key][3]) and int(self.coordinate_dict[key][3]) - int(self.coordinate_dict[key][4]) <=55:
                            self.coordinate_dict[key][3]=int(self.coordinate_dict[key][3])-(int(self.coordinate_dict[key][3]) - int(self.coordinate_dict[key][4]))-2
                        if int(self.coordinate_dict[key][2]) <= int(self.coordinate_dict[key][1]) and int(self.coordinate_dict[key][1]) - int(self.coordinate_dict[key][2]) <=55:
                            self.coordinate_dict[key][2]=int(self.coordinate_dict[key][2])+(int(self.coordinate_dict[key][1]) - int(self.coordinate_dict[key][2]))+2
                for key in self.coordinate_dict:
                    #check if the sequece with all 3 domains is growing monotonically
                    if all(int(x)<int(y) for x, y in zip(self.coordinate_dict[key], self.coordinate_dict[key][1:])) and int(self.coordinate_dict[key][1]) - int(self.coordinate_dict[key][0])>75 and int(self.coordinate_dict[key][3]) - int(self.coordinate_dict[key][2])>70 and int(self.coordinate_dict[key][5]) - int(self.coordinate_dict[key][4])>70:
                        final_id_list.append(key)
                #create records for the full and processed sequences with 3 domains
                new_rec_list=list()   
                full_rec_list=list()
                dummy_dict = copy.deepcopy(self.coordinate_dict)
                self.coordinate_dict = defaultdict(list)
                for key in final_id_list:
                    self.coordinate_dict[key]=dummy_dict[key]  
                self.three_dom_count = len(final_id_list)
                for record in SeqIO.parse(self.cry_query,"fasta"):
                    if record.description != record.id:
                        name = record.id + '|' +'|'.join('|'.join(record.description.split(' ')[1:]).split('_'))
                    else:
                        name= record.id + '|' +'|'.join('|'.join(record.description.split(' ')[0:]).split('_'))
                    if name in final_id_list:
                         #check if multiple identical ids are in the query to get warning
                         if record.id in ids_check_set:
                             check_flag=1
                         ids_check_set.add(record.id)
                         #substract 1 from the starting index because the hmmer output is 1-based
                         start = min([int(x) for x in self.coordinate_dict[name]])-1
                         stop = max([int(x) for x in self.coordinate_dict[name]])
                         full_rec_list.append(SeqRecord(Seq(str(record.seq),
                                           generic_protein),
                                           id=name.split('|')[0], 
                                           description=" ".join(name.split('|')[1:])))
                         new_rec_list.append(SeqRecord(Seq(str(record.seq[start:stop]),
                                           generic_protein),
                                           id=name.split('|')[0], 
                                            description=" ".join(name.split('|')[1:])))
                #save full and processed sequences
                SeqIO.write(new_rec_list,
                        os.path.join(
                        os.path.realpath(
                        self.query_dir),
                        'raw_processed_{}.fasta'.format(
                        self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0])), 
                        "fasta")
                SeqIO.write(full_rec_list,
                        os.path.join(
                        os.path.realpath(
                        self.query_dir),
                        'raw_full_{}.fasta'.format(
                        self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0])), 
                        "fasta")
                if check_flag==1:
                    if not self.silent_mode:
                        print('Warning! Identical ids in the query, the uncertain mapping can occur')
                    self.logger.warning('Warning! Identical ids in the query, the uncertain mapping can occur')
                if not self.silent_mode:
                    print('{} sequences recieved'.format(self.init_count))
                    print('{} potential Cry toxins found'.format(self.one_dom_count+self.two_dom_count+self.three_dom_count))
                    print('{} toxins with one domain'.format(self.one_dom_count))
                    print('{} toxins with two domains'.format(self.two_dom_count))
                    print('{} toxins with three domains'.format(self.three_dom_count))
                self.logger.info('{} sequences recieved'.format(self.init_count))
                self.logger.info('{} potential Cry toxins found'.format(self.one_dom_count+self.two_dom_count+self.three_dom_count))
                self.logger.info('{} toxins with one domain'.format(self.one_dom_count))
                self.logger.info('{} toxins with two domains'.format(self.two_dom_count))
                self.logger.info('{} toxins with three domains'.format(self.three_dom_count))

        if self.mode == 'fd':
            self.find_cry(self.cry_query)
            if self.hmmer_result_flag ==1:
                if not self.silent_mode:
                    print('No toxins found')
                self.logger.info('No toxins found')
            else:
                #find the domains mode: perform hmmsearch on the full model, then repeat it on domain models
                #launch the search for the full potential cry-toxins
                #extract the fasta-file after hmmsearch and use it as a query for the domain search 
                #process the hmmsearch output with sed
                cmd_take_file = subprocess.call("cp {0} \
                                            {1}; \
                                            mv {2} \
                                            {3}; \
                                            sed -i 's/[0-9]*|//g' {3}; \
                                            sed -i 's/\/[0-9]*-[0-9]*//g' {3}; \
                                            sed -i 's/ \[subseq from]\ / /g' \
                                            {3}".format(os.path.join(
                                            os.path.realpath(
                                            self.query_dir),
                                            'full_toxins',
                                           '{}_full_extracted.fasta'.format(
                                            self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0])),
                                            os.path.join(
                                            os.path.realpath(
                                            self.query_dir)), 
                                            os.path.join(
                                            os.path.realpath(
                                            self.query_dir),
                                            '{}_full_extracted.fasta'.format(
                                            self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0])), 
                                            os.path.join(
                                            os.path.realpath(
                                            self.query_dir),
                                            '{}.fasta'.format(
                                            self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0]))),
                                            shell=True)
                #lauch the domain search
                self.find_domains('{}'.format(
                                       os.path.join(
                                       os.path.realpath(
                                       self.query_dir),
                                       '{}.fasta'.format(
                                       self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0]))))
                #further steps are identical to do_mode and are described above
                self.id_list = self.cry_3D_ids_extractor()
                final_id_list=list()
                if not self.silent_mode:
                    print("Performing the Cry toxins processing")
                self.logger.info("Performing the Cry toxins processing")
                self.coordinate_dict = defaultdict(list)
                pre_dict=defaultdict(dict)
                ids_check_set=set()
                check_flag=0
                self.first_filter_count = len(list(SeqIO.parse('{}'.format(
                                              os.path.join(
                                              os.path.realpath(
                                              self.query_dir),
                                              '{}.fasta'.format(
                                              self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0]))),
                                              "fasta")))
                dom_start=int(self.processing_flag)
                for i in range(dom_start,4):
                    for record in SeqIO.parse(os.path.join(
                                              os.path.realpath(
                                              self.query_dir),
                                              'domains',
                                              self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0]+
                                              '_D{}_extracted.fasta'.format(i)),
                                              "fasta"):
                        if '|' in record.id:
                            name = record.id.split('|')[1].split('/')[0] + '|' +\
                               '|'.join('|'.join(record.description.split('] ')[1].split(' ')).split('_'))
                        else:
                            name = record.id.split('/')[0] + '|' +\
                               '|'.join('|'.join(record.description.split('] ')[1].split(' ')).split('_'))
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
                for key in self.coordinate_dict:
                    if all(int(x)<int(y) for x, y in zip(self.coordinate_dict[key], self.coordinate_dict[key][1:])) and int(self.coordinate_dict[key][1]) - int(self.coordinate_dict[key][0])>75:
                        final_id_list.append(key)
                new_rec_list=list()   
                full_rec_list=list()
                dummy_dict = copy.deepcopy(self.coordinate_dict)
                self.coordinate_dict = defaultdict(list)
                for key in final_id_list:
                    self.coordinate_dict[key]=dummy_dict[key] 
                self.three_dom_count = len(final_id_list)
                for record in SeqIO.parse('{}'.format(os.path.join(
                                                      os.path.realpath(
                                                      self.query_dir),
                                                      '{}.fasta'.format(
                                                      self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0]))),
                                                      "fasta"):
                    if record.description != record.id:
                        name = record.id + '|' +'|'.join('|'.join(record.description.split(' ')[1:]).split('_'))
                    else:
                        name= record.id + '|' +'|'.join('|'.join(record.description.split(' ')[0:]).split('_'))
                    if name in final_id_list:
                         if record.id in ids_check_set:
                             check_flag=1
                         ids_check_set.add(record.id)
                         start = min([int(x) for x in self.coordinate_dict[name]])-1
                         stop = max([int(x) for x in self.coordinate_dict[name]])
                         new_rec_list.append(SeqRecord(Seq(str(record.seq[start:stop].upper()),
                                          generic_protein),
                                          id=name.split('|')[0], 
                                          description=" ".join(name.split('|')[1:])))
                         full_rec_list.append(SeqRecord(Seq(str(record.seq),
                                          generic_protein),
                                          id=name.split('|')[0], 
                                          description=" ".join(name.split('|')[1:])))
                SeqIO.write(new_rec_list,
                        os.path.join(
                        os.path.realpath(
                        self.query_dir),
                        'raw_processed_{}.fasta'.format(
                        self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0])), 
                        "fasta")
                SeqIO.write(full_rec_list,
                        os.path.join(
                        os.path.realpath(
                        self.query_dir),
                        'raw_full_{}.fasta'.format(
                        self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0])), 
                        "fasta")
                if check_flag==1:
                    if not self.silent_mode:
                        print('Warning! Identical ids in the query, the uncertain mapping can occur')
                    self.logger.warning('Warning! Identical ids in the query, the uncertain mapping can occur')
                if not self.silent_mode:
                    print('{} sequences recieved'.format(self.init_count))
                    print('{} potential Cry toxins found'.format(self.one_dom_count+self.two_dom_count+self.three_dom_count))
                    print('{} toxins with one domain'.format(self.one_dom_count))
                    print('{} toxins with two domains'.format(self.two_dom_count))
                    print('{} toxins with three domains'.format(self.three_dom_count))
                self.logger.info('{} sequences recieved'.format(self.init_count))
                self.logger.info('{} potential Cry toxins found'.format(self.one_dom_count+self.two_dom_count+self.three_dom_count))
                self.logger.info('{} toxins with one domain'.format(self.one_dom_count))
                self.logger.info('{} toxins with two domains'.format(self.two_dom_count))
                self.logger.info('{} toxins with three domains'.format(self.three_dom_count))


        #create a file with the domain mappings for the each cry-toxin, posessing all 3 domains
        if self.hmmer_result_flag !=1:
            #create a mapping file only if the hmmsearch output is not empty
            with open(os.path.join(
                      os.path.realpath(
                      self.query_dir),
                      'proteins_domain_mapping_processed_{}.bed'.format(
                      self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0])),'w') as csv_file:
                my_writer = csv.writer(csv_file, delimiter='\t') 
                #create a bed-file for mapping to specify start and stop coordinates, id and description
                #init_row = ['id','domain' ,'start', 'stop', 'description']
                #my_writer.writerow(init_row)
                for key in self.coordinate_dict:
                    #iterate over starting indicies in the dictionary of coordinates
                    if self.processing_flag!="2":
                        for i in range(0,5,2):
                            if i==0:
                                row=[key.split('|')[0],
                                  int(self.coordinate_dict[key][i])-int(self.coordinate_dict[key][0]),
                                  int(self.coordinate_dict[key][i+1])-int(self.coordinate_dict[key][0])-1,
                                 'domain_{}'.format(1)+'_'+key.split('|')[0]] 
                            else:
                                if i==2:
                                    dn=2
                                else:
                                    dn=3
                                row=[key.split('|')[0],
                                 int(self.coordinate_dict[key][i])-int(self.coordinate_dict[key][0])-1,
                                 int(self.coordinate_dict[key][i+1])-int(self.coordinate_dict[key][0])-1,
                                 'domain_{}'.format(dn)+'_'+key.split('|')[0]] 
                            my_writer.writerow(row)
                    else:
                        for i in range(0,3,2):
                            if i==0:
                                row=[key.split('|')[0],
                                  int(self.coordinate_dict[key][i])-int(self.coordinate_dict[key][0]),
                                  int(self.coordinate_dict[key][i+1])-int(self.coordinate_dict[key][0])-1,
                                 'domain_{}'.format(2)+'_'+key.split('|')[0]]  
                            else:
                                if i==2:
                                    dn=3
                                row=[key.split('|')[0],
                                 int(self.coordinate_dict[key][i])-int(self.coordinate_dict[key][0])-1,
                                 int(self.coordinate_dict[key][i+1])-int(self.coordinate_dict[key][0])-1,
                                 'domain_{}'.format(dn)+'_'+key.split('|')[0]]  
                            my_writer.writerow(row)



            #create a bed-file with the full protein mappings
            with open(os.path.join(
                      os.path.realpath(
                      self.query_dir),
                      'proteins_domain_mapping_full_{}.bed'.format(
                        self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0])),'w') as csv_file:
                my_writer = csv.writer(csv_file, delimiter='\t') 
                #init_row = ['id','domain' ,'start', 'stop', 'description']
                #my_writer.writerow(init_row)
                for key in self.coordinate_dict:
                    if self.processing_flag!="2":
                        for i in range(0,5,2):
                            if i==0:
                                row=[key.split('|')[0],
                                  int(self.coordinate_dict[key][i]),
                                  int(self.coordinate_dict[key][i+1])-1,
                                 'domain_{}'.format(1)+'_'+key.split('|')[0]] 
                            else:
                                if i==2:
                                    dn=2
                                else:
                                    dn=3
                                row=[key.split('|')[0],
                                 int(self.coordinate_dict[key][i])-1,
                                 int(self.coordinate_dict[key][i+1])-1,
                                 'domain_{}'.format(dn)+'_'+key.split('|')[0]] 

                        
                            my_writer.writerow(row)
                    else:
                        for i in range(0,3,2):
                            if i==0:
                                row=[key.split('|')[0],
                                  int(self.coordinate_dict[key][i]),
                                  int(self.coordinate_dict[key][i+1])-1,
                                 'domain_{}'.format(2)+'_'+key.split('|')[0]] 
                            else:
                                if i==2:
                                    dn=3
                                row=[key.split('|')[0],
                                 int(self.coordinate_dict[key][i])-1,
                                 int(self.coordinate_dict[key][i+1])-1,
                                 'domain_{}'.format(dn)+'_'+key.split('|')[0]] 
                            my_writer.writerow(row)

            #save the original dictionary of coordinates for checking
            with open(os.path.join(
                      os.path.realpath(
                      self.query_dir),
                      'coordinate_matches_{}.txt'.format(
                      self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0])),'w') as csv_file:
                my_writer = csv.writer(csv_file, delimiter='\t') 
                for key in self.coordinate_dict:
                    my_writer.writerow([key]+self.coordinate_dict[key])

    def annotate_raw_output(self):
        """
        Performs annotating raw sequences with diamond ocer database from btnomenclature
        Uses full records for annotation to get potentially new toxins
        """
        #ids for toxins, differ from btnomenclature
        new_records = list()
        self.new_ids=dict()
        un_count=0
        total_count=0
        if not self.silent_mode:
            print("Annotating the raw output with diamond") 
        self.logger.info("Annotating the raw output with diamond") 
        
        #diamond should have a database file in the same directory with the query. 
        #copying database for pefrorming blastp
        cmd_pre_dia = subprocess.call('cd {0}; \
                                       cp {1} .; \
                                       touch diamond.log'.format(
                                       os.path.join(
                                       os.path.realpath(
                                       self.query_dir)),
                                       os.path.realpath(
                                       os.path.join(
                                       os.path.dirname(__file__),
                                       'data',
                                       "diamond_data",
                                       "cry_nomenclature.dmnd"))), 
                                       shell=True) 
        #performing diamond blastp, save only the first hit


        cmd_dia = subprocess.call('cd {0};\
                                   {1} blastp \
                                    -d cry_nomenclature \
                                    -q raw_full_{2}.fasta \
                                    -o diamond_matches_{2}.txt \
                                    --al aligned_{2}.fa \
                                    --un unaligned_{2}.fa \
                                    --max-target-seqs 1 \
                                    --log --verbose 2>> diamond.log; \
                                    rm cry_nomenclature.dmnd; \
                                    mv *.log logs/; \
                                    mv *gned* logs/'.format(
                                    os.path.join(
                                    os.path.realpath(
                                    self.query_dir)),
                                    os.path.realpath(
                                    os.path.join(
                                    os.path.dirname(__file__),
                                    "include",
                                    "diamond")),
                                    self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0]), 
                                    shell=True)
        #analyse diamond matches to get new toxins
        with open(os.path.join(
                  os.path.realpath(
                  self.query_dir),
                  'diamond_matches_{}.txt'.format(
                  self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0])),'r') as csv_file:
            my_reader = csv.reader(csv_file, delimiter='\t') 
            for row in my_reader:
                total_count+=1
                if float(row[2])<100.0:
                    self.new_ids[row[0]]=row[1]+'|'+ str(row[2])
                    un_count+=1
        #extract unique sequences from the raw sequences
        for init_rec in SeqIO.parse(os.path.join(
                                    os.path.realpath(
                                    self.query_dir),
                                    'raw_processed_{}.fasta'.format(
                                    self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0])), 
                                    "fasta"):
           if init_rec.id in self.new_ids.keys():
               new_records.append(SeqRecord(Seq(str(init_rec.seq),
                                   generic_protein),
                                   id=self.new_ids[init_rec.id].split('|')[0]+'(' + 
                                   self.new_ids[init_rec.id].split('|')[1]+')'+
                                   '_'+init_rec.description.split()[0],
                                   description=init_rec.description))
        if not self.silent_mode:
            print('{} sequences matched with the database'.format(total_count))
            print('{} toxins different from the database found'.format(un_count))
        self.logger.info('{} sequences matched with the database'.format(total_count)) 
        self.logger.info('{} toxins different from the database found'.format(un_count)) 
        SeqIO.write(new_records,os.path.join(
                                os.path.realpath(
                                self.query_dir),
                                'unique_{}.fasta'.format(
                                self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0])), 
                                "fasta")

        new_coord_dict = defaultdict(list)
        #print(self.coordinate_dict)
        #print(self.new_ids)
        for key in self.coordinate_dict:
            if key.split("|")[0] in self.new_ids.keys():
                #print(self.new_ids[key.split("|")[0]].split('|')[0] +"("+self.new_ids[key.split("|")[0]].split('|')[1]+")_"+key.split("|")[0])
                new_coord_dict[self.new_ids[key.split("|")[0]].split('|')[0] +"("+self.new_ids[key.split("|")[0]].split('|')[1]+")_"+key.split("|")[0]]=self.coordinate_dict[key]
        with open(os.path.join(
                      os.path.realpath(
                      self.query_dir),
                      'unique_proteins_domain_mapping_{}.bed'.format(
                        self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0])),'w') as csv_file:
                my_writer = csv.writer(csv_file, delimiter='\t') 
                for key in new_coord_dict:
                    if self.processing_flag!="2":
                        for i in range(0,5,2):
                            if i==0:
                                row=[key.split('|')[0],
                                  int(new_coord_dict[key][i])-int(new_coord_dict[key][0]),
                                  int(new_coord_dict[key][i+1])-1-int(new_coord_dict[key][0]),
                                 'domain_{}'.format(1)+'_'+key.split('|')[0]] 
                            else:
                                if i==2:
                                    dn=2
                                else:
                                    dn=3
                                row=[key.split('|')[0],
                                 int(new_coord_dict[key][i])-1-int(new_coord_dict[key][0]),
                                 int(new_coord_dict[key][i+1])-1-int(new_coord_dict[key][0]),
                                 'domain_{}'.format(dn)+'_'+key.split('|')[0]] 

                        
                            my_writer.writerow(row)
                    else:
                        for i in range(0,3,2):
                            if i==0:
                                row=[key.split('|')[0],
                                  int(new_coord_dict[key][i])-int(new_coord_dict[key][0]),
                                  int(new_coord_dict[key][i+1])-1-int(new_coord_dict[key][0]),
                                 'domain_{}'.format(2)+'_'+key.split('|')[0]] 
                            else:
                                if i==2:
                                    dn=3
                                row=[key.split('|')[0],
                                 int(new_coord_dict[key][i])-1-int(new_coord_dict[key][0]),
                                 int(new_coord_dict[key][i+1])-1-int(new_coord_dict[key][0]),
                                 'domain_{}'.format(dn)+'_'+key.split('|')[0]] 
                            my_writer.writerow(row)

        #move intermediate files to the log directory if the annotation flag is not specified
        if self.mode == 'fd':
            cmd_clean_up = subprocess.call("mv {0} \
                                               {1};".format(os.path.join(
                                                      os.path.realpath(
                                                      self.query_dir),
                                                      '{}.fasta'.format(
                                                      self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0])),
                                                      os.path.join(
                                                      os.path.realpath(
                                                      self.query_dir),
                                                      'first_search_{}.fasta'.format(
                                                      self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0]))),
                                               shell=True)
        if not self.nucl_type and not self.annot_flag:
            cmd_clean_up = subprocess.call('cd {0}; \
                                            mv *coordinate_matches* logs/;\
                                            mv *diamond_matches* logs/'.format(os.path.join(
                                                      os.path.realpath(
                                                      self.query_dir))), 
                                            shell=True)
       
    def make_summary_table(self):
        """
        Annotates raw sequences by searching metadata from NCBI ipg database
        """
        if not self.silent_mode:
            print('Searching for the metadata')
        self.logger.info('Searching for the metadata') 
        #it is better to specify your e-mail address for the correct ncbi searching
        Entrez.email = "{}".format(self.email)
        summary_dict=defaultdict(dict)
        #load information from the diamond blastp search
        with open(os.path.join(
                  os.path.realpath(
                  self.query_dir),
                  'diamond_matches_{}.txt'.format(
                   self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0])),'r') as csv_file:
            my_reader = csv.reader(csv_file, delimiter='\t') 
            for row in my_reader:
                summary_dict[row[0]]=defaultdict(list)
                summary_dict[row[0]]['init']=row[1:3]
        #create a tsv-file with metadata
        making_smb = subprocess.call("touch {}".format(os.path.join(
                                                       os.path.realpath(
                                                      self.query_dir),
                                                      'annotation_table_{}.tsv'.format(
                                                      self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0]))), 
                                                      shell = True)
        for init_rec in SeqIO.parse(os.path.join(
                                    os.path.realpath(
                                    self.query_dir),
                                    'raw_processed_{}.fasta'.format(
                                    self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0])), 
                                    "fasta"):
          #save initial information about the sequence
           if init_rec.id in summary_dict.keys():
               summary_dict[init_rec.id]['init'].append(init_rec.description)
        #save the header of the tsv-file
        init_row = ['protein_id', 
                    'initial_description', 
                    'top_cry_hit', 
                    'cry_identity', 
                    'source', 
                    'nucl_accession', 
                    'start',
                    'stop', 
                    'strand',
                    'ipg_prot_id',
                    'ipg_prot_name', 
                    'organism', 
                    'strain',
                    'assembly\n']
        #write rows for the each sequence 
        f = open(os.path.join(
                 os.path.realpath(
                 self.query_dir),
                 'annotation_table_{}.tsv'.format(
                 self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0])),
                 "w")
        f.write('\t'.join(init_row))
        f.close() 
        for key in summary_dict: 
            #checking if id is in ipg
            try:
                #getting an ipg-table for the query
                handle = Entrez.efetch(db="protein",
                                       rettype='ipg',
                                       retmode='text',
                                       id =key)
                #parse the ipg output
                handle_list=[el.split('\t') for el in handle.read().split('\n')]
                hit_counter=0
                #iterate over the hits from ipg
                for i in range(len(handle_list)-1):
                   if len(handle_list[i+1])>2:
                       summary_dict[key]['hit'+str(hit_counter)]=handle_list[i+1]
                       hit_counter+=1
                iter_num=len(summary_dict[key].keys())
                for i in range(iter_num-1):
                    #save all the data for the initial sequence from the query
                    if i==0:
                        row=[key]+[summary_dict[key]['init'][2]]+\
                            [summary_dict[key]['init'][0]]+\
                            [summary_dict[key]['init'][1]]+\
                            summary_dict[key]['hit'+str(i)][1:]
                    else:
                        #mark other hits as additional
                        row=['--']*4+summary_dict[key]['hit'+str(i)][1:]
                    f = open(os.path.join(
                             os.path.realpath(
                             self.query_dir),
                             'annotation_table_{}.tsv'.format(
                             self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0])),
                             "a+")
                    f.write('\t'.join(row) + '\n') 
                    f.close()
                time.sleep(3)
            except:
                if not self.silent_mode:
                    print('Warning! Bad link for!', key)
                self.logger.warning('Warning! Bad link for!', key) 
        #clean from intermediate files if no nucleotide uploadings are specified
        if not self.nucl_type:
            cmd_clean_up = subprocess.call('cd {0}; \
                                           mv *coordinate_matches* logs/;\
                                           mv *diamond_matches* logs/'.format(os.path.join(
                                           os.path.realpath(
                                           self.query_dir))), 
                                           shell=True)

           

    def upload_nucl(self):
        """
        Uploads nucleotide sequences based on previous annotation step
        """
        if not self.silent_mode:
            print('Uploading the nucleotide sequences')
        self.logger.info('Uploading the nucleotide sequences') 
        Entrez.email = "{}".format(self.email)
        #a dictionary for start and stop positions of the nucleotide sequences from the annotation tsv-file
        keys_for_nucl = defaultdict(list)
        #a dictionary for domain positions from the protein domain mappings
        domain_coord_dict = defaultdict(list)
        with open(os.path.join(
                  os.path.realpath(
                  self.query_dir),
                  'coordinate_matches_{}.txt'.format(
                   self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0])), 'r') as csv_file:
            my_reader = csv.reader(csv_file, delimiter='\t')
            for row in my_reader:
                 domain_coord_dict[row[0]]=row[1:]
        init_count=0
        with open(os.path.join(
                  os.path.realpath(
                  self.query_dir),
                  'annotation_table_{}.tsv'.format(
                  self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0])), 'r') as csvfile:
             my_reader = csv.reader(csvfile, delimiter='\t') 
             for row in my_reader:
                 #search nucleotide sequences only for the initial sequences from the tsv-file
                 if row[0]!= '--' and row[0]!='protein_id':
                     init_count+=1
                     keys_for_nucl[row[0]].append('|'.join((row[1]).split(' ')))
                     #add info about the accession, start, stop an strand orientation
                     keys_for_nucl[row[0]].extend([row[5], row[6], row[7], row[8], row[2], row[3]])
        if keys_for_nucl[list(keys_for_nucl.keys())[0]][0] not in domain_coord_dict:
            for key in keys_for_nucl:
                keys_for_nucl[key][0]=keys_for_nucl[key][0]+'|'+keys_for_nucl[key][0]
        #upload full records if the fn flag is specified
        f_nuc_recs = []
        p_nuc_recs = []
        if self.nucl_type == 'fn':
            for key in keys_for_nucl:
                #checking the link from the annotation table
                try:
                    #search for the fasta file 
                    handle = Entrez.efetch(db="nucleotide",
                                            rettype='fasta',
                                            retmode='text', 
                                            id = keys_for_nucl[key][1])
                    fasta_rec = SeqIO.read(handle, "fasta")
                    handle.close()
                    #perform reverse complement if strand is -
                    if keys_for_nucl[key][4] == '+':
                        f_nuc_recs.append(SeqRecord(Seq(str(fasta_rec.seq[int(keys_for_nucl[key][2])-
                                          1:int(keys_for_nucl[key][3])]),
                                          generic_dna),
                                          id=keys_for_nucl[key][5]+
                                          '('+keys_for_nucl[key][6]+
                                          ')'+'_'+keys_for_nucl[key][0].split('|')[0], 
                                          description=" ".join(keys_for_nucl[key][0].split('|'))))
                    elif keys_for_nucl[key][4] == '-':
                        f_nuc_recs.append(SeqRecord(Seq(str(fasta_rec.seq[int(keys_for_nucl[key][2])-
                                         1:int(keys_for_nucl[key][3])].reverse_complement()),
                                         generic_dna),
                                         id=keys_for_nucl[key][5]+
                                         '('+keys_for_nucl[key][6]+
                                         ')'+'_'+keys_for_nucl[key][0].split('|')[0], 
                                         description=" ".join(keys_for_nucl[key][0].split('|'))))
                    time.sleep(3)
                except:
                    if not self.silent_mode:
                        print('Warning! Download error for ', 
                               keys_for_nucl[key][2], 
                               '(' + keys_for_nucl[key][1]+ ')')
                    self.logger.warning('Warning! Download error for ', 
                                         keys_for_nucl[key][2], 
                                         '(' + keys_for_nucl[key][1]+ ')') 

            SeqIO.write(f_nuc_recs,os.path.join(
                                   os.path.realpath(
                                   self.query_dir),
                                   '{}_full_nucl.fna'.format(
                                   self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0])),
                                   "fasta")

        #upload full processed records if the pn flag is specified
        elif self.nucl_type == 'pn':
            for key in keys_for_nucl:
                try:
                    handle = Entrez.efetch(db="nucleotide",
                                            rettype='fasta',
                                            retmode='text', 
                                            id = keys_for_nucl[key][1])
                    fasta_rec = SeqIO.read(handle, "fasta")
                    handle.close()
                    if keys_for_nucl[key][4] == '+':
                        #dum_seq is a full nucleotide sequence
                        dum_seq = str(fasta_rec.seq[int(keys_for_nucl[key][2])-1:int(keys_for_nucl[key][3])])
                        #cut left and rigth flanking sequenses from the duf_seq according to the protein domain mappings
                        #transform to nucleotide coordinates by substracting 1 and multiplying by 3
                        p_nuc_recs.append(SeqRecord(
                                        Seq(dum_seq[(
                                        min([int(x) for x in domain_coord_dict[keys_for_nucl[key][0]]])-1)*3:
                                        max([int(x) for x in domain_coord_dict[keys_for_nucl[key][0]]])*3],
                                        generic_dna),
                                        id=keys_for_nucl[key][5]+
                                        '('+keys_for_nucl[key][6]+
                                        ')'+'_'+keys_for_nucl[key][0].split('|')[0], 
                                        description=" ".join(keys_for_nucl[key][0].split('|'))))
                    elif keys_for_nucl[key][4] == '-':
                        dum_seq = str(Seq(
                                      str(fasta_rec.seq[int(keys_for_nucl[key][2])-1:
                                      int(keys_for_nucl[key][3])].reverse_complement()),
                                      generic_dna))
                        p_nuc_recs.append(SeqRecord(
                                      Seq(dum_seq[(
                                      min([int(x) for x in domain_coord_dict[keys_for_nucl[key][0]]])-1)*3:
                                      max([int(x) for x in domain_coord_dict[keys_for_nucl[key][0]]])*3],
                                      generic_dna),
                                      id=keys_for_nucl[key][5]+
                                     '('+keys_for_nucl[key][6]+
                                     ')'+'_'+keys_for_nucl[key][0].split('|')[0], 
                                     description=" ".join(keys_for_nucl[key][0].split('|'))))
                    time.sleep(3)
                except:
                    if not self.silent_mode:
                        print('Warning! Download error for ', 
                               keys_for_nucl[key][2], 
                               '(' + keys_for_nucl[key][1]+ ')')
                    self.logger.warning('Warning! Download error for ', 
                                         keys_for_nucl[key][2], 
                                         '(' + keys_for_nucl[key][1]+ ')') 
            SeqIO.write(p_nuc_recs,os.path.join(
                                   os.path.realpath(
                                   self.query_dir),
                                   '{}_processed_nucl.fna'.format(
                                   self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0])),
                                   "fasta")
        

        #upload both full and processed records if the an flag is specified
        elif self.nucl_type == 'an':
            for key in keys_for_nucl:
                try:
                    handle = Entrez.efetch(db="nucleotide",
                                           rettype='fasta',
                                           retmode='text',
                                           id = keys_for_nucl[key][1])
                    fasta_rec = SeqIO.read(handle, "fasta")
                    handle.close()
                    if keys_for_nucl[key][4] == '+':
                        dum_seq = str(fasta_rec.seq[int(keys_for_nucl[key][2])-1:
                                  int(keys_for_nucl[key][3])])
                        p_nuc_recs.append(SeqRecord(
                                       Seq(dum_seq[(
                                       min([int(x) for x in domain_coord_dict[keys_for_nucl[key][0]]])-1)*3:
                                       max([int(x) for x in domain_coord_dict[keys_for_nucl[key][0]]])*3],
                                       generic_dna),
                                       id=keys_for_nucl[key][5]+
                                       '('+keys_for_nucl[key][6]+
                                       ')'+'_'+keys_for_nucl[key][0].split('|')[0], 
                                       description=" ".join(keys_for_nucl[key][0].split('|'))))
                        f_nuc_recs.append(SeqRecord(
                                       Seq(str(fasta_rec.seq[int(keys_for_nucl[key][2])-1:
                                       int(keys_for_nucl[key][3])]),
                                       generic_dna),
                                       id=keys_for_nucl[key][5]+
                                       '('+keys_for_nucl[key][6]+
                                       ')'+'_'+keys_for_nucl[key][0].split('|')[0], 
                                       description=" ".join(keys_for_nucl[key][0].split('|'))))
                    elif keys_for_nucl[key][4] == '-':
                        dum_seq = str(Seq(str(fasta_rec.seq[int(keys_for_nucl[key][2])-1:
                                      int(keys_for_nucl[key][3])].reverse_complement()),
                                      generic_dna))
                        p_nuc_recs.append(SeqRecord(
                                      Seq(dum_seq[(
                                      min([int(x) for x in domain_coord_dict[keys_for_nucl[key][0]]])-1)*3:
                                      max([int(x) for x in domain_coord_dict[keys_for_nucl[key][0]]])*3],
                                      generic_dna),
                                      id=keys_for_nucl[key][5]+
                                      '('+keys_for_nucl[key][6]+
                                      ')'+'_'+keys_for_nucl[key][0].split('|')[0], 
                                      description=" ".join(keys_for_nucl[key][0].split('|'))))
                        f_nuc_recs.append(SeqRecord(
                                      Seq(str(fasta_rec.seq[int(
                                      keys_for_nucl[key][2])-1:
                                      int(keys_for_nucl[key][3])].reverse_complement()),
                                      generic_dna),
                                      id=keys_for_nucl[key][5]+
                                      '('+keys_for_nucl[key][6]+
                                      ')'+'_'+keys_for_nucl[key][0].split('|')[0], 
                                      description=" ".join(keys_for_nucl[key][0].split('|')))) 
                    time.sleep(3)
                except:
                    if not self.silent_mode:
                        print('Warning! Download error for ', 
                               keys_for_nucl[key][2], 
                               '(' + keys_for_nucl[key][1]+ ')')
                    self.logger.warning('Warning! Download error for ', 
                                         keys_for_nucl[key][2], 
                                         '(' + keys_for_nucl[key][1]+ ')') 
                     

            SeqIO.write(p_nuc_recs,os.path.join(
                                   os.path.realpath(
                                   self.query_dir),
                                   '{}_processed_nucl.fna'.format(
                                   self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0])),
                                   "fasta")

            SeqIO.write(f_nuc_recs,os.path.join(
                                   os.path.realpath(
                                   self.query_dir),
                                   '{}_full_nucl.fna'.format(
                                   self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0])),
                                   "fasta")
        if not self.silent_mode:
            print('{} nucleotide sequences downloaded'.format(max([len(p_nuc_recs), 
                                                       len(f_nuc_recs)])))
        self.logger.info('{} nucleotide sequences downloaded'.format(max([len(p_nuc_recs), 
                                                       len(f_nuc_recs)])))
        self.nuc_count=max([len(p_nuc_recs),len(f_nuc_recs)])

    def map_nucl(self):
        """
        Create mapping files for nucleotide sequences, requires downloaded sequences
        """
        keys_for_nucl = defaultdict(list)
        domain_coord_dict = defaultdict(list)   
        #read information about othe riginal protein coordinates for unprocessed sequences     
        with open(os.path.join(
                  os.path.realpath(
                  self.query_dir),
                  'coordinate_matches_{}.txt'.format(
                  self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0])), 
                  'r') as csv_file:
            my_reader = csv.reader(csv_file, delimiter='\t')
            for row in my_reader:
                 domain_coord_dict[row[0]]=row[1:]
        init_count=0
        #add info about stop and start of nucleotide sequences
        with open(os.path.join(
                  os.path.realpath(
                  self.query_dir),
                  'annotation_table_{}.tsv'.format(
                  self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0])),
                  'r') as csvfile:
            my_reader = csv.reader(csvfile, delimiter='\t')
            for row in my_reader:
                if row[0]!= '--' and row[0]!='protein_id':                          
                    init_count+=1
                    keys_for_nucl[row[0]].append('|'.join((row[1]).split(' ')))
                    keys_for_nucl[row[0]].extend([row[5], row[6], row[7], row[8], row[2], row[3]])                              
        if keys_for_nucl[list(keys_for_nucl.keys())[0]][0] not in domain_coord_dict:
            for key in keys_for_nucl:
                keys_for_nucl[key][0]=keys_for_nucl[key][0]+'|'+keys_for_nucl[key][0]
        #mappings for processed sequences
        with open(os.path.join(
                  os.path.realpath(
                  self.query_dir),
                  'nucl_domain_mapping_processed_{}.bed'.format(
                        self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0])), 
                   'w') as csvfile:
             my_writer = csv.writer(csvfile, delimiter='\t') 
             #init_row = ['id','domain', 'start', 'stop', 'description']
             #my_writer.writerow(init_row)
             for key in keys_for_nucl:
                 try:
                     if self.processing_flag!="2":
                         for i in range(0,5,2):
                             if i==0:
                            #substract coordinates of the first domain start to get the mappings for processed sequences
                                row=[keys_for_nucl[key][5]+
                                 '('+keys_for_nucl[key][6]+')'+'_'+
                                 keys_for_nucl[key][0].split('|')[0],
                                 'domain {}'.format(1)+'_'+key.split('|')[0],
                                 (int(domain_coord_dict[keys_for_nucl[key][0]][i])-1)*3-
                                 (int(domain_coord_dict[keys_for_nucl[key][0]][0])-1)*3,
                                 int(domain_coord_dict[keys_for_nucl[key][0]][i+1])*3-
                                 (int(domain_coord_dict[keys_for_nucl[key][0]][0])-1)*3-1,
                                  " ".join(keys_for_nucl[key][0].split('|'))]  
                             else:
                                if i==2:
                                    dn=2
                                else:
                                    dn=3
                                row=[keys_for_nucl[key][5]+
                                '('+keys_for_nucl[key][6]+')'+'_'+
                                keys_for_nucl[key][0].split('|')[0],
                                'domain {}'.format(dn)+'_'+key.split('|')[0],
                                (int(domain_coord_dict[keys_for_nucl[key][0]][i])-1)*3-
                                (int(domain_coord_dict[keys_for_nucl[key][0]][0])-1)*3-1,
                                int(domain_coord_dict[keys_for_nucl[key][0]][i+1])*3-
                                (int(domain_coord_dict[keys_for_nucl[key][0]][0])-1)*3-1, 
                                " ".join(keys_for_nucl[key][0].split('|'))]      
                             my_writer.writerow(row)
                     else:
                         for i in range(0,3,2):
                             if i==0:
                                row=[keys_for_nucl[key][5]+
                                 '('+keys_for_nucl[key][6]+')'+'_'+
                                 keys_for_nucl[key][0].split('|')[0],
                                 'domain {}'.format(2)+'_'+key.split('|')[0],
                                 (int(domain_coord_dict[keys_for_nucl[key][0]][i])-1)*3-
                                 (int(domain_coord_dict[keys_for_nucl[key][0]][0])-1)*3,
                                 int(domain_coord_dict[keys_for_nucl[key][0]][i+1])*3-
                                 (int(domain_coord_dict[keys_for_nucl[key][0]][0])-1)*3-1,
                                  " ".join(keys_for_nucl[key][0].split('|'))]  
                             else:
                                if i==2:
                                    dn=3
                                row=[keys_for_nucl[key][5]+
                                '('+keys_for_nucl[key][6]+')'+'_'+
                                keys_for_nucl[key][0].split('|')[0],
                                'domain {}'.format(dn)+'_'+key.split('|')[0],
                                (int(domain_coord_dict[keys_for_nucl[key][0]][i])-1)*3-
                                (int(domain_coord_dict[keys_for_nucl[key][0]][0])-1)*3-1,
                                int(domain_coord_dict[keys_for_nucl[key][0]][i+1])*3-
                                (int(domain_coord_dict[keys_for_nucl[key][0]][0])-1)*3-1, 
                                " ".join(keys_for_nucl[key][0].split('|'))]      
                             my_writer.writerow(row)
                 except:
                     pass
        #save mappings for full nucleotide sequences
        with open(os.path.join(
                  os.path.realpath(
                  self.query_dir),
                  'nucl_domain_mapping_full_{}.bed'.format(
                   self.cry_query.split('/')[len(self.cry_query.split('/'))-1].split('.')[0])), 
                   'w') as csvfile:
             my_writer = csv.writer(csvfile, delimiter='\t') 
             #init_row = ['id','domain', 'start', 'stop', 'description']
             #my_writer.writerow(init_row)
             for key in keys_for_nucl:
                 try:
                     if self.processing_flag!="2":
                         for i in range(0,5,2):
                             if i==0:
                                row=[keys_for_nucl[key][5]+
                                 '('+keys_for_nucl[key][6]+')'+'_'+
                                 keys_for_nucl[key][0].split('|')[0],
                                 'domain_{}'.format(1)+'_'+key.split('|')[0],
                                 (int(domain_coord_dict[keys_for_nucl[key][0]][i])-1)*3,
                                 int(domain_coord_dict[keys_for_nucl[key][0]][i+1])*3-1,
                                  " ".join(keys_for_nucl[key][0].split('|'))]  
                             else:
                                if i==2:
                                    dn=2
                                else:
                                    dn=3
                                row=[keys_for_nucl[key][5]+
                                '('+keys_for_nucl[key][6]+')'+'_'+
                                keys_for_nucl[key][0].split('|')[0],
                                'domain_{}'.format(dn)+'_'+key.split('|')[0],
                                (int(domain_coord_dict[keys_for_nucl[key][0]][i])-1)*3-1,
                                int(domain_coord_dict[keys_for_nucl[key][0]][i+1])*3-1, 
                                " ".join(keys_for_nucl[key][0].split('|'))]      
                             my_writer.writerow(row)
                     else:
                         for i in range(0,5,2):
                             if i==0:
                                row=[keys_for_nucl[key][5]+
                                 '('+keys_for_nucl[key][6]+')'+'_'+
                                 keys_for_nucl[key][0].split('|')[0],
                                 'domain_{}'.format(2)+'_'+key.split('|')[0],
                                 (int(domain_coord_dict[keys_for_nucl[key][0]][i])-1)*3,
                                 int(domain_coord_dict[keys_for_nucl[key][0]][i+1])*3-1,
                                  " ".join(keys_for_nucl[key][0].split('|'))]  
                             else:
                                if i==2:
                                    dn=3
                                row=[keys_for_nucl[key][5]+
                                '('+keys_for_nucl[key][6]+')'+'_'+
                                keys_for_nucl[key][0].split('|')[0],
                                'domain_{}'.format(dn)+'_'+key.split('|')[0],
                                (int(domain_coord_dict[keys_for_nucl[key][0]][i])-1)*3-1,
                                int(domain_coord_dict[keys_for_nucl[key][0]][i+1])*3-1, 
                                " ".join(keys_for_nucl[key][0].split('|'))]      
                             my_writer.writerow(row)
                 except:
                     pass

        cmd_clean_up = subprocess.call('cd {}; \
                                        mv *coordinate_matches* logs/;\
                                        mv *diamond_matches* logs/'.format(os.path.join(
                                                      os.path.realpath(
                                                      self.query_dir))), 
                                                      shell=True)

    def launch_racer(self,kmer):
        """
        launches pathracer on gfa-file
        Args
        =====
        kmer: kmer size for pathracer
        """
        if not self.silent_mode:
            print('Searching for the sequences from the gfa file')
        self.logger.info('Searching for the sequences from the gfa file')
        #launch pathracer
        cmd_race = subprocess.call('{0} \
                                  {1} \
                                  {2} {3} \
                                  -l 0.3 \
                                  --output {4} \
                                  -t {5}  > /dev/null'.format(
                                          os.path.join(
                                          os.path.dirname(
                                          __file__),
                                          "include",
                                          "pathracer"),
                                          os.path.join(
                                          os.path.dirname(
                                          __file__),
                                          'data',
                                          'models',
                                          'Complete.hmm'),
                                           os.path.realpath(
                                           self.cry_query),
                                           kmer,
                                           os.path.join(
                                           os.path.realpath(
                                           self.query_dir),
                                           'pathracer_output'),
                                           self.hm_threads), 
                                           shell=True)


        #save all the sequences and the pathracer log
        cmd_merge = subprocess.call('cd {0}; \
                                   if [  -f *seqs.fa ]; \
                                   then cat *seqs.fa >> mearged.fasta;\
                                   fi;\
                                   cp pathracer.log ../logs/;\
                                   '.format(os.path.join(
                                           os.path.realpath(
                                           self.query_dir),
                                           'pathracer_output')), 
                                           shell=True) 
        exist_check_flag = re.sub("b",'',
                                  re.sub("\'",'', 
                                  str(subprocess.check_output("cd {}; \
                                  if [ ! -f mearged.fasta ];\
                                  then echo 'no'; \
                                  fi".format(os.path.join(
                                           os.path.realpath(
                                           self.query_dir),
                                           'pathracer_output')),
                                            shell =True).strip())))
        if exist_check_flag == 'no':
            self.racer_flag = 1
        if self.racer_flag == 0:
            #rename ids from the pathracer output
            #double {} is used in awk to exclude it from string formating
            cmd_rename = subprocess.call("cd {}; \
                                   if [  -s mearged.fasta ]; \
                                   then awk \'/^>/{{print \">seq\" ++i; next}}{{print}}\' mearged.fasta > tmp ;\
                                   rm mearged.fasta;\
                                   mv tmp mearged.fasta;\
                                   echo 'done';\
                                   fi".format(os.path.join(
                                           os.path.realpath(
                                           self.query_dir),
                                           'pathracer_output')), 
                                            shell=True)
      

    def use_spades(self):
        """
        launches spades on illumina reads, usese mrtaspades if --meta flag is specified
        """
        if not self.silent_mode:
            print('Building the assembly graph')
        self.logger.info('Building the assembly graph')
        if str(self.meta_flag)=='True': 
            #use mataspades if the --meta flag is specified
            cmd_spades = subprocess.call('{0} \
                                         --meta \
                                         -1 {1} -2 {2} \
                                         -o {3} \
                                         -t {4} > /dev/null'.format(
                                          os.path.join(
                                          os.path.dirname(
                                          __file__),
                                          "include",
                                          "SPAdes-3.13.1-Linux",
                                          "bin",
                                          "spades.py"),
                                          os.path.realpath(
                                          self.forw),
                                          os.path.realpath(
                                          self.rev),
                                          os.path.join(
                                          os.path.realpath(
                                          self.query_dir),
                                          'assembly'),
                                          self.hm_threads),
                                          shell=True) 

            cmd_merge = subprocess.call('cd {}; \
                                         cp *.log ../logs/'.format(
                                          os.path.join(
                                          os.path.realpath(
                                          self.query_dir),
                                          'assembly')),
                                          shell=True)
        else:
            cmd_spades = subprocess.call('{0} \
                                         -1 {1} -2 {2} \
                                         -o {3} \
                                         -t {4} > /dev/null'.format(
                                          os.path.join(
                                          os.path.dirname(
                                          __file__),
                                          "include",
                                          "SPAdes-3.13.1-Linux",
                                          "bin",
                                          "spades.py"),
                                          os.path.realpath(
                                          self.forw),
                                          os.path.realpath(
                                          self.rev),
                                          os.path.join(
                                          os.path.realpath(
                                          self.query_dir),
                                          'assembly'),
                                          self.hm_threads), 
                                          shell=True) 
            cmd_merge = subprocess.call('cd {}; \
                                         cp *.log ../logs/'.format(os.path.join(
                                          os.path.realpath(
                                          self.query_dir),
                                          'assembly')), 
                                          shell=True)       


class Crylauncher:
    def __init__(self):
        pass
    def LaunchProcessor(od,fi,hm,pr,th, ma, r, a, nu, mra,k,fr,rr,meta,s):
    #check if the input file exists
        if fi:
            if not os.path.isfile(fi):
                raise Exception('The input file does not exist')  
            if not mra:
                fasta_flag = re.sub("b",'',
                                re.sub("\'",'', 
                                str(subprocess.check_output("grep '>' {} | wc -l".format(os.path.realpath(fi)), 
                                shell =True).strip())))
                if int(fasta_flag)==0:
                    raise Exception('No records are present in the fasta file')  
                
        #check if modes are not mixed 
        if (mra and fr) or (fr and fi) or (fi and re and fr) or (fi and re and fr and mra):
            raise Exception('You should not mix the --pathracer mode with the assembly mode and the fasta-serching mode; choose only the one option')
        #check if only the one file with reads is specified
        if fr and not rr or rr and not fr:
            raise Exception('Specify both forward and reverse reads')
        if fr:
            #check if files with reads are present
            if not os.path.isfile(rr) and not os.path.isfile(fr):
                raise Exception('The Files with reverse and forward reads do not exist')
            elif not os.path.isfile(rr):
                raise Exception('The File with reverse reads does not exists')
            elif not os.path.isfile(fr):
                raise Exception('The File with forward reads does not exists')
        #initialize CryProcessor class
        cr = CryProcessor(od, fi, hm,pr, th, ma, r, nu,a,k,meta,fr,rr,s)
        if not mra and not fr:
            #pipeline for the protein fasta file
            cr.cry_digestor()
            if cr.hmmer_result_flag != 1:
            #go to the next steps only if the hmmsearch result is not empty
                cr.annotate_raw_output()
                if a: 
                    cr.make_summary_table()
                if nu:
                    cr.upload_nucl()
                    if cr.nuc_count!=0:
                        cr.map_nucl()
        elif mra and not fr:
        #pipeline for the gfa file
            cr.launch_racer(k)
            if cr.racer_flag != 1:
            #go to the next step only if the pathracer output is not empty
                fi = os.path.join(os.path.realpath(
                              od),
                              'pathracer_output',
                              'mearged.fasta')
                cr = CryProcessor(od, fi, hm,pr, th, ma, r, nu,a,k,meta,fr,rr,s)
                cr.cry_digestor()
                if cr.hmmer_result_flag != 1:
                    cr.annotate_raw_output()
                    if a: 
                        cr.make_summary_table()
                    if nu:
                        cr.upload_nucl()
                        if cr.nuc_count!=0:
                            cr.map_nucl()
            else:
                if not s:
                    print('No toxins found in the assembly graph')
                cr.logger.info('No toxins found in the assembly graph')

        elif fr:
        #full pipeline with spades assembly
            cr.use_spades()
            fi = os.path.join(
                          os.path.realpath(
                          od),
                          'assembly',
                          'assembly_graph_with_scaffolds.gfa') 
            cr = CryProcessor(od, fi, hm,pr, th, ma, r, nu,a,k,meta,fr,rr,s)
            cr.launch_racer(k)
            if cr.racer_flag != 1:
                fi = os.path.join(
                              os.path.realpath(
                              od),
                              'pathracer_output',
                              'mearged.fasta')
                cr = CryProcessor(od, fi, hm,pr, th, ma, r, nu,a,k,meta,fr,rr,s)
                cr.cry_digestor()
                if cr.hmmer_result_flag != 1:
                    cr.annotate_raw_output()
                    if a: 
                        cr.make_summary_table()
                    if nu:
                        cr.upload_nucl()
                        if cr.nuc_count!=0:
                            cr.map_nucl()
            else:
                if not s:
                    print('No toxins found in the assembly graph')
                cr.logger.info('No toxins found in the assembly graph')
        if not s:
            print('CryProcessor has finished, thanks for crying (with) us!')
        cr.logger.info('CryProcessor has finished, thanks for crying (with) us!')
        move_log = subprocess.call('mv cry_processor.log {}'.format(os.path.join(
                              os.path.realpath(
                              od),
                              'logs')), 
                               shell=True)
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Example usage: pyhton3 cry_processor.py -fi <input.fasta> -od <the output directory>')
    parser.add_argument('-fi', 
                        help='Enter the input file: the fasta file or the gfa file', 
                        metavar='Fasta file/GFA file',
                        type=str, 
                        default=None)
    parser.add_argument('-hm', 
                        help='Path to Hmmer if Hmmer is installed locally', 
                        metavar='Hmmer_directory',
                        type=str, 
                        default='')
    parser.add_argument('-pr', 
                        help='The processig type: 1 for extracting all the domains, 2 for extrating 2-3 domains only', 
                        metavar='Int',
                        type=str, 
                        default=1)
    parser.add_argument('-th', 
                        help='Number of threads for Hmmer/PathRacer/SPAdes', 
                        metavar='Int',
                        type=str, 
                        default=8)
    parser.add_argument('-ma', 
                        help='e-mail address for the NCBI annotation', 
                        metavar='Str',
                        type=str, 
                        default='')
    parser.add_argument('-od', 
                         help='The output directory', 
                         metavar='Str',
                         type=str, 
                         required=True)
    parser.add_argument('-r', 
                         help='The pipeline type: do - domain only search with the subsequent combining the searching results;\
                         fd - searching for the potential Cry toxins with the subsequent processing',
                          metavar='Str',
                         type=str, 
                         default='do')
    parser.add_argument('--annotate', 
                        '-a',
                         action='store_true',
                         help='make the final output annotation with the IPG database')
    parser.add_argument('-nu', 
                         help='Uploading the nucleotide records: fn - uploading the full sequences, \
                         pn - uploading the processed subsequences, \
                         an - the both processed and unprocessed sequences',
                         metavar='Nucl_uploading_type', 
                         type=str, default='')
    parser.add_argument('--pathracer',
                        '-pa',
                         action='store_true',
                         help='Searching for the Cry toxins from the gfa file via PathRacer')
    parser.add_argument('-fo', 
                        help='The forward illumina reads', 
                        metavar='Fastq file',
                        type=str, 
                        default=None)
    parser.add_argument('-re', 
                        help='The reverse illumina reads', 
                        metavar='Fastq file',
                        type=str, 
                        default=None)
    parser.add_argument('--meta',
                        '-m', 
                        action='store_true',
                        help='The metagenomic mode for SPAdes')
    parser.add_argument('-k', 
                        help='The k-mer size for PathRacer', 
                        metavar='Int',
                        type=str, 
                        default=55)
    parser.add_argument('--silent',
                        '-s', 
                        action='store_true',
                        help='The silent mode')

    parser.set_defaults(feature=True)
    args = parser.parse_args()
    od,fi,hm,pr,th, ma, r, a, nu, mra,k,fr,rr,meta,s = args.od, args.fi, args.hm, args.pr,args.th, args.ma, args.r, args.annotate, args.nu,args.pathracer,args.k,args.fo,args.re,args.meta,args.silent
    Crylauncher.LaunchProcessor(od,fi,hm,pr,th, ma, r, a, nu, mra,k,fr,rr,meta,s)

