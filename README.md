# Cry Processor
Cry processor is high-troughtput tool for mining cry-toxins from fasta-files or directly from illumina reads.

## About Cry processor

Cry processor is a python-written tool for searching cry toxins from sequence data or from protein fasta files. It Includes hmm-based searching of potential toxins,  obtaining information about domains, extracting cry toxins with 3 domains only and comparing fonuded toxins with btnomenclature. A regime for performing toxins search directly from Illumina reads iplies asembly graph building using spades and subsequent toxins mining directly from assebmly graph.

## Cry processor pipeline

Here full pipeline is described (from illumina reads). By first spades or metaspades is performed to get assembly from fastq-files. After that potential cry toxins (with at least 30% identity to hmm-consensus) are extracted from assembly paths via pathracer. Then hmmsearch is implemented to find cry toxins domains in obtained sequences. Searching result are than combined to get toxins possesing all three domains. Coordinates of domains are used to cut flanking sequences and save domains and linkers bertween domains. Full sequences are used to compare obtained toxins with btnomenclature database using diamond blastp. Non-identical sequences then are extracted as pitenial new toxins. For all found sequences (both identical to presented in bt_nomenclature and novel) online ipg-annotation is performed (to see anotation output watch sections below). Finally, nucleotide sequences corresponding to protein sequences of found toxins are downloaded. Metadata is upload only if accession numbers are present in the quiery.

### Prerequisites
<ul>
  <li>python (3.7 or higher) </li>
  <li>Biopython (1.66 or higher)
  </li>
</ul>

To install Biopython use the folloinw comand

```
~$ pip install biopython

```


### Installing

To install this tool clone this repository to your PC.

```
~$ git clone https://github.com/cry_Processor
```
After downloading a tool is ready to use.

## Running and using tool

## Quick usage
To extract cry toxins from the protein fasta file simply run the following command:
```
~$ python3 cry_processor.py -fi input.faa  -od output_dir
```
This command will automatically search for cry-toxins in fasta-file with amino acid sequences

## 


### Annotation output
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

## Tool options: 
```
-fi <input.fasta> or <input.gfa> Input file in fasta-format of in gfa-format
-hm <path to hmmer directory> Path to hmmemer directory if you want to use local hmmer
-pr <1 or 2> Processig type: 1 for extracting all domains, 2 for extrating 2-3 domains only
-th Nubmer of threads for hmmer/spades/pathracer
-ma <e-mail> e-mail address for connecting to NCBI
-od <output dir> output directory
-r <do or fd> Working regime: do - domain only search only; fd full toxins mining with subsequent domain search
-a (--annotate) Perform data anotation with ipg (only if accession numbers are present)
-nu <fn or pn or an> Upload nucleotide sequences: fn - full sequences, pn -processed sequences, an - both variants
-pa (--path_racer) - Launching pathracer on gfa file
-fo <input_1.fastq> - forward illumina reads
-re <input_2.fastq> - reverse illumina reads
-n (--meta) - flag for  specifying metagenomic spades regime
-k <number> - K-mer size for pathracer
```


## Author

* **Anton Shikov** - *Initial work* - [anton-shikov](https://github.com/anton-shikov)


