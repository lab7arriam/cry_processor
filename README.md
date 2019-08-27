[![Build Status](https://travis-ci.com/anton-shikov/cry_processor.svg?branch=master)](https://travis-ci.com/anton-shikov/cry_processor) ![python3.5](https://img.shields.io/badge/language-python%203.5-blue.svg)
![image](https://img.shields.io/badge/Platforms-linux--64-orange.svg)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# CryProcessor
CryProcessor is a high-troughtput tool for the Cry toxins mining from the fasta-files or directly from the illumina reads.

## About CryProcessor

CryProcessor is a python-written tool for searching and extracting Cry toxins from illumina sequence data or from the protein fasta files. It includes several parts: an hmm-based scanning for potential Cry toxins, obtaining information about the domains, extracting Cry toxins  with 3 domains only and comparing found toxins with Bt nomenclature. <br> <br> The mode for performing the toxins search directly from the illumina reads implies building an assembly graph (using <i>SPAdes</i>) and the subsequent mining toxins directly from the obtained assebmly graph. 

## CryProcessor Pipeline

The following text stands for the full pipeline description (for the illumina reads). To start, <i>SPAdes</i> (http://cab.spbu.ru/software/spades/) or <i>metaSPAdes</i> (http://cab.spbu.ru/software/meta-spades/) are implemented to get the assembly graph from the fastq-files. After that, the potential Cry toxins (with at least 30% identity to the hmm-consensus) are extracted from the assembly paths via <i>PathRacer</i> (http://cab.spbu.ru/software/pathracer/). Then <i>hmmsearch</i> (http://hmmer.org/) is used to find Cry toxin domains in the obtained sequences. In the next step, the results of <i>hmmsearch</i> are combined to get the toxins that posses all three domains. <br> <br>The coordinates of the domains are used to cut flanking sequences and save the domains with the corresponding linkers. The full sequences (without processing procedure) are used to compare the obtained toxins with the Bt nomenclature database via <i>diamond blastp</i> (https://github.com/bbuchfink/diamond). The non-identical sequences are extracted and marked as the potentially new toxins. <br><br>For all the found sequences (both identical to presented in Bt nomenclature and the novel sequences) an online ipg-annotation (Identical Protein Group) is performed (to see the annotation output read the annotation output section below). Finally, nucleotide sequences, corresponding to the protein sequences of the found toxins, are downloaded. Metadata will be uploaded only if the accession numbers are present in the query.

## Installation and Usage
### Prerequisites
<ul>
  <li>python (3.7 or higher); </li>
  <li>Biopython (1.66 or higher).
  </li>
</ul>

To install <i>Biopython</i> use the following command:

```
~$ pip install biopython
```

### Installing

To install CryProcessor clone git repository to your PC:

```
~$ git clone https://github.com/lab7arriam/cry_processor
```
After downloading, CryProcessor is ready to use.

You can also add CryProcessor to the PATH typing:

```
~$ PATH=$PATH:/path/to/install
```

If you want to use prebuild Docker container with CryProcessor, pull it using the following command:

```
~$ docker pull lab7arriam/cry_processor
```

## Running and Using Tool

### Quick Usage
To extract Cry toxins from the protein fasta file simply execute the following command:
```
~$ /path/to/install/cry_processor.py -fi input.faa  -od output_dir
```
This command will automatically search for the Cry toxins in the fasta file with amino acid sequences.

If you use Docker container, you should use following command:

```
~$ docker run --rm -v /path/to/data:/data lab7arriam/cry_processor cry_processor -fi /data/<input file> -od /data/<output_dir> [args]
```

### Supported Input Formats

<ul>
  <li>fasta files with protein sequences; </li>
  <li>gfa files (representing genome assebly graph);</li>
  <li>forward and reverse illumina reads.</li>
</ul>

### Tool Options: 
The full list of tool options:
```
-fi <input.fasta> or <input.gfa> an input file in the fasta format of in the gfa format
-hm <path to hmmer directory> path to the hmmer directory if you want to use local hmmer
-pr <1 or 2> the processig type: 1 for extracting all the domains, 2 for extrating 2-3 domains only (default 1)
-th the nubmer of threads for hmmer/SPAdes/PathRacer (default 8)
-ma <e-mail> an e-mail address for the connecting to NCBI
-od <output dir> the output directory
-r <do or fd> the running mode: do - the domain only search only; fd - the full toxins mining with the subsequent domain search
-a (--annotate) perform the data anotation with ipg (only if the accession numbers are present)
-nu <fn or pn or an> upload the nucleotide sequences: fn - the full sequences, pn - the processed sequences, an - the both variants
-pa (--pathracer) - launching PathRacer on the gfa file
-fo <input_1.fastq> - forward illumina reads
-re <input_2.fastq> - reverse illumina reads
-n (--meta) - the flag for specifying the metagenomic mode forSPAdes
-k <number> - the K-mer size for PathRacer (default 21)
-s (--silent) - disable the console output
```

### Using Tool for the Fasta Files
To use the tool for the files in the fasta format execute the command, presented in quick usage, you can also specify the annotation (writing an e-mail address is strongly recommended):
```
~$ /path/to/install/cry_processor.py -fi input.faa  -od output_dir -ma <e-mail address> --annotate
```
Use the <i>-nu</i> flag to download nucleotide sequences:

```
~$ /path/to/install/cry_processor.py -fi input.faa  -od output_dir -ma <e-mail address> --annotate -nu pn
```
The pipeline of searching could be performed in two modes:
<ul>
  <li><i>do</i> - the domain only mode. Searches the Cry toxin domains from the full queiry, then combines this data to extract the toxins with 3 domains; </li>
  <li><i>fd</i> - the find domains mode. At the begining, <i>hmmsearch</i> with the full Cry toxin model is performed. Next, the domains are extracted and combined. This mode is quicker, as far as the domain search is performed on previously extracted sequences. </li>
</ul>

### Using Tool for the .gfa Files
You can apply the Cry toxins search directly from the assembly graph in the gfa format with the following commad:

```
~$ /path/to/install/cry_processor.py -fi input.gfa  -od output_dir --path_racer
```
### Using Tool for the Illumina Reads
This mode includes the reads assembly with SPAdes and the subsequent hmm-based toxins mining. To implement this use the following command:

```
~$ /path/to/install/cry_processor.py -fo forward_reads.fastq -re reverse_reads.fastq -od output_dir 
```
If you want to search for the Cry toxins from metagenomic reads specify <i>--meta</i> flag:

```
~$ /path/to/install/cry_processor.py -fo forward_reads.fastq -re reverse_reads.fastq -od output_dir --meta
```
### Tips for the Correct Running

Note that you cannot mix modes:
<ul>
  <li>Do not use the <i>--pathracer</i> flag with the illumina fastq files query flags (<i>-fo</i> and <i>-re</i>); </li>
  <li>Do not mix the <i>-fi</i> agrument with <i>-fo</i> and <i>-re</i> arguments; </li>
  <li>Do not mix the <i>--pathracer</i> flag and the <i>--meta</i> flag; </li>
  <li>Do not specify the <i>--meta</i> mode with the <i>-fi</i> agrument; </li>
  <li>You should use both the <i>-fo</i> and <i>-re</i> argumens together; </li>
</ul>
Using the <i>-nu</i> flag is possible only if the <i>--annotate</i> flag is specified. <br> Note that performing the anotation is not recommended for the gfa and assembly modes, because the online annotation is impossible without the accession numbers in the query.


### The Annotation Output
Using the <i>--annotate</i> flag will perform the NCBI-search in the ipg database for the submitted accession numbers within the query and return gathered information in tsv-format with the following structure:
<table>
  <tr>
    <td>protein_id</td>
    <td>initial_description</td>
    <td>top_cry_hit</td>
    <td>cry_identity</td>
    <td>source</td>
    <td>nucl_accession</td>
    <td>start</td>
    <td>stop</td>
    <td>strand</td>
    <td>ipg_prot_id</td>
    <td>ipg_prot_name</td>
    <td>organism</td>
    <td>strain</td>
    <td>assembly</td>
  </tr>
  <tr>
    <td>AFU17015.1</td>
    <td>AFU17015.1 pesticidal crystal proteinBacillus thuringiensis MC28</td>
    <td>Cry4Cc1</td>
    <td>99.4</td>
    <td>INSDC</td>
    <td>CP003690.1</td>
    <td>58993</td>
    <td>62628</td>
    <td>+</td>
    <td>AFU17015.1</td>
    <td>pesticidal crystal protein</td>
    <td>Bacillus thuringiensis MC28</td>
    <td>MC28</td>
    <td>GCA_000300475.1</td>
  </tr>
</table>

Columns description:
<ul>
  <li>protein_id - initial sequence id in the query </li>
  <li>initial_description - protein description in the query </li>
  <li>top_cry_hit - Cry protein from Bt nomenclature with the best identity score according to <i>diamond</i></li>
  <li>cry_identity - identity score with the closest Cry protein from Bt nomenclature </li>
  <li>source - database source of the sequence </li>
  <li>nucl_accession - accession number for the nucleotide sequence </li>
  <li>start - starting position in the full nucleotide sequence </li>
  <li>stop - ending position in the full nucleotide sequence </li> 
  <li>strand - DNA strand orientation for the nucleotide sequence </li>
  <li>ipg_prot_id - protein id according to the ipg database </li>
  <li>ipg_prot_name - protein description to the ipg database </li>
  <li>organism - taxon description </li>
  <li>strain - strain of the organism </li>
  <li>assembly - assembly ids where the nucleotide sequence is present</li>
</ul>
Rows possesing the first four coloums refer to the initial sequences (which are found in the query), results for all the identical proteins are marked with --. <br>Here's an expamle of such rows:

<table>
  <tr>
    <td> -- </td>
    <td> -- </td>
    <td> -- </td>
    <td> -- </td>
    <td>RefSeq</td>
    <td>NZ_CP010111.1</td>
    <td>142760</td>
    <td>146194</td>
    <td>-</td>
    <td>WP_080989235.1</td>
    <td>pesticidal protein</td>
    <td>Bacillus thuringiensis serovar indiana</td>
    <td>HD521</td>
    <td>GCF_001183785.1</td>
  </tr>
</table>
Note that all the rows, marked with -- in the first four columns, are identical to the initial sequenses in the queiry with the full columns, located above.

### The Output Files Structure
In the output directory, specified with the <i>-od</i> flag, the <i>cry_extraction</i> directory is created. It can contain the following subdirectories with the following files:

<ul>
  <li>cry_extraction/assembly - this directory is created if the assembly mode is enabled. It contains the files and directories that refer to the <i>SPAdes</i> output, for instanse, the assembly graph with the contigs in the gfa format. To see the full list of the <i>SPAdes</i> output look at the <i>SPAdes</i> manual (http://cab.spbu.ru/software/spades/); </li>
  <li>cry_extraction/pathracer - this directory includes the <i>PathRacer</i> output. It will be created if the <i>pathracer</i> mode or the assembly mode are enabled. To see the full list of the <i>PathRacer</i> output read the manual (http://cab.spbu.ru/software/pathracer/); </li>
  <li>cry_extraction/full_toxins -this directory is created if the <i>fd</i> mode is enabled. It contains the following files: 
   	<ul>
           <li> &lt;input&gt;_full_extracted.sto - the result of <i>hmmsearch</i> with the full-toxin model in the sto format; </li>
           <li> &lt;input&gt;_full_extracted.fasta - the result of <i>hmmsearch</i> with the full-toxin model in the fasta format;  </li>
         </ul>
  <li>cry_extraction/domains. This directory contains the results of <i>hmmsearch</i> for all the domain models in the sto and fasta formats: 
   	<ul>
           <li> &lt;input&gt;_D1_extracted.sto; </li>
           <li> &lt;input&gt;_D1_extracted.fasta;  </li>
           <li> &lt;input&gt;_D2_extracted.sto; </li>
           <li> &lt;input&gt;_D2_extracted.fasta;  </li>
           <li> &lt;input&gt;_D3_extracted.sto; </li>
           <li> &lt;input&gt;_D3_extracted.fasta;  </li>
         </ul>
  </li>
  <li>cry_extraction/logs. This directory contains the log files for the different stages of the pipeline: 
   	<ul>
           <li> spades.log - the log file for the assembly process; </li>
           <li> pathracer.log - the log file for the hmm-based search for the Cry toxins from the gfa file;  </li>
           <li> full_extraction.log - the log file for performing <i>hmmsearch</i> on the full-toxin model;  </li>
           <li> domains_extraction.log - the log file for performing <i>hmmsearch</i> on the domain models;  </li>
           <li> coordinate_matches_&lt;input&gt;.txt - the dictionary of the domain mappings for the query; </li>
           <li> diamond.log - the log file for performing <i>diamond blastp</i> agaist the Bt nomenclature database;  </li>
           <li> aligned_&lt;input&gt;.fa - the aligned sequences from <i>diamond blastp</i>;  </li>
           <li> unaligned_&lt;input&gt;.fa - the unligned sequences from <i>diamond blastp</i>;  </li>
           <li> diamond_matches_&lt;input&gt;.txt - the results of <i>diamond blastp</i> in the table form;  </li>
           <li> cry_processor.log - the main CryProcessor logfile;  </li>
         </ul>
  </li>
  <li>cry_extraction/first_search_&lt;input&gt;.fasta - the result of <i>hmmsearch</i> with the full model on the query in the fasta format. This file is later used as a queiry for the domains search; </li>
  <li>cry_extraction/raw_full_&lt;input&gt;.fasta. This file contains the full protein sequences of the Cry toxins possesing all three domains; </li>
  <li>cry_extraction/raw_processed_&lt;input&gt;.fasta -this file contains the processed protein sequences (without left and right lanking parts) of the Cry toxins possesing all three domains; </li>
  <li>cry_extraction/proteins_domain_mapping_full_&lt;input&gt;.bed - the file in the bed format, containing the domain coordinate mappings for the uprocessed protein sequences; </li>
  <li>cry_extraction/proteins_domain_mapping_processed_&lt;input&gt;.bed - the file in the bed format, containing the domain coordinate mappings for the processed protein sequences; </li>
  <li>cry_extraction/unique_&lt;input&gt;.fasta - this file contains the potentially new Cry toxins, differing from Bt nnomenclature afrer <i>diamond blastp</i>; </li>
  <li>cry_extraction/annotation_table_&lt;input&gt;.tsv - the result of the ipg annotation of the found toxins (the structure of the file is described above); </li>
  <li>cry_extraction/&lt;input&gt;_full_nucl.fna - the full nucleotide sequences for the found toxins; </li>
  <li>cry_extraction/&lt;input&gt;_processed_nucl.fna - the processed nucleotide sequences for the found toxins; </li>
  <li>cry_extraction/nucl_domain_mapping_full_&lt;input&gt;.bed - the file in the bed format, containing the domain coordinate mappings for the uprocessed nucleotide sequences; </li>
  <li>cry_extraction/nucl_domain_mapping_processed_&lt;input&gt;.bed - the file in the bed format, containing the domain coordinate mappings for the processed nucleotide sequences; </li>
</ul>

Note, that in the files, marked with the <i>raw</i> prefix the initial accession numbers from the query are used as ids, while in the files, obtained after <i>diamond blastp</i> (&lt;input&gt;_full_nucl.fna, &lt;input&gt;_processed_nucl.fna and cry_extraction/unique_&lt;input&gt;.fasta) the id structure is modified in the following way: <br>
<ul>
<li>Cry53Aa1(40.7)_WP_103591149.1 </li>
<li>&lt;top Cry protein hit from Bt nomenclature&gt;(&lt;the identity score with this hit&gt;)_&lt;the initial accession number&gt;</li>
</ul>

## References

1. Bankevich A., Nurk S., Antipov D., Gurevich A., Dvorkin M., Kulikov A. S., Lesin V., Nikolenko S., Pham S., Prjibelski A., Pyshkin A., Sirotkin A., Vyahhi N., Tesler G., Alekseyev M. A., Pevzner P. A. SPAdes: A New Genome Assembly Algorithm and Its Applications to Single-Cell Sequencing. Journal of Computational Biology, 2012 
2. Buchfink B.,  Xie C. &  Huson D.H. Fast and Sensitive Protein Alignment using DIAMOND. Nature Methods, 2015
3. Finn R.D., Clements J.,  Eddy S.R. HMMER Web Server: Interactive Sequence Similarity Searching. Nucleic Acids Research, 2011.
4. Nurk S., Meleshko D., Korobeynikov A., Pevzner P. A. metaSPAdes: a new versatile de novo metagenomics assembler. Genome Research, 2017

## Author

* **Anton Shikov** - *Initial work* - [anton-shikov](https://github.com/anton-shikov)


