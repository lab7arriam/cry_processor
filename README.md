# Cry Processor
Cry processor is high-troughtput tool for mining cry-toxins from fasta-files or directly from illumina reads.

## About Cry processor

Cry processor is a python-written tool for searching cry toxins from sequence data or from protein fasta files. It Includes hmm-based searching of potential toxins,  obtaining information about domains, extracting cry toxins with 3 domains only and comparing fonuded toxins with btnomenclature. A regime for performing toxins search directly from Illumina reads iplies asembly graph building using spades and subsequent toxins mining directly from assebmly graph.

## Cry processor pipeline

Here full pipeline is described (from illumina reads). By first spades or metaspades is performed to get assembly from fastq-files. After that potential cry toxins (with at least 30% identity to hmm-consensus) are extracted from assembly paths via pathracer. Then hmmsearch is implemented to find cry toxins domains in obtained sequences. Searching result are than combined to get toxins possesing all three domains. Coordinates of domains are used to cut flanking sequences and save domains and linkers bertween domains. Full sequences are used to compare obtained toxins with btnomenclature database using diamond blastp. Non-identical sequences then are extracted as pitenial new toxins. For all found sequences (both identical to presented in bt_nomenclature and novel) online ipg-annotation is performed (to see anotation output watch sections below). Finally, nucleotide sequences corresponding to protein sequences of found toxins are downloaded. Metadata is upload only if accession numbers are present in the quiery.

##Installation and usage
### Prerequisites
<ul>
  <li>python (3.7 or higher) </li>
  <li>Biopython (1.66 or higher)
  </li>
</ul>

To install Biopython use the following comand

```
~$ pip install biopython

```

### Installing

To install this tool clone tool repository to your PC.

```
~$ git clone https://github.com/cry_Processor
```
After downloading cry processor is ready to use.

## Running and using tool

### Quick usage
To extract cry toxins from the protein fasta file simply run the following command:
```
~$ python3 cry_processor.py -fi input.faa  -od output_dir
```
This command will automatically search for cry-toxins in fasta-file with amino acid sequences 

### Supported input formats

<ul>
  <li>fasta files with protein sequences </li>
  <li>gfa file (genome assebly graph)</li>
  <li>forward and reverse illumina reads</li>
</ul>

### Tool options: 
Full list of tool options:
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

## Using tool for fasta files
## Using tool for gfa files
## Using tool for illumina reads

### Annotation output
Using --annotating flag would perform NCBI-search in ipg(identical protein group) database for submitted accession number and would return information in tsv-format with the following structure:
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

Coloums description:
<ul>
  <li>protein_id - initial id in quiery </li>
  <li>initial_description - protein description in quiery </li>
  <li>top_cry_hit - Cry-protein from btnomenclature with the biggest identity score </li>
  <li>cry_identity - identity score with closest Cry-protein from btnomenclature </li>
  <li>source - Database source of sequence </li>
  <li>nucl_accession - accession number for nucleotide sequence </li>
  <li>start - starting position in full nucleotide sequence </li>
  <li>stop - ending position in full nucleotide sequence </li>
  <li>strand - DNA strand orientation for nucleotide sequence </li>
  <li>ipg_prot_id - protein according to ipg database </li>
  <li>ipg_prot_name - protein description to ipg database </li>
  <li>organism - taxon description </li>
  <li>strain - strain of the organism </li>
  <li>assembly - assembly id where nucleotide sequence is present</li>
</ul>
Rows possesing first 4 coloums fefer to initial sequences (those found in quiery), results for all identical proteins are marked with --. Hee is an expamle of such rows:

<table>
  <tr>
    <td>--</td>
    <td>--</td>
    <td>--</td>
    <td>--</td>
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
Note that all rows marked with -- in first four rows are identical to initial sequense in queiry with full rows

### Output files structure



## Author

* **Anton Shikov** - *Initial work* - [anton-shikov](https://github.com/anton-shikov)


