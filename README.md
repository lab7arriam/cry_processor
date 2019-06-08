# CryProcessor
CryProcessor is a high-troughtput tool for cry-toxins mining from fasta-files or directly from illumina reads.

## About CryProcessor

CryProcessor is a python-written tool for searching and extracting cry toxins from illumina sequence data or from protein fasta files. It includes several parts: hmm-based searching of potential cry-toxins, obtaining information about domains, extracting cry-toxins only with 3 domains and comparing found toxins with btnomenclature. A regime for performing toxins search directly from illumina reads implies building assembly graph, using  <i> SPAdes <i>, and subsequent mining toxin directly from obtained assebmly graph. 

## CryProcessor pipeline

The fololowing text stands for full pipeline (for illumina reads). To start, <i>SPAdes<i> or <i>metaSPAdes<i> are implemented to get assembly graph from the fastq-files. After that, potential cry toxins (with at least 30% identity to the hmm-consensus) are extracted from assembly paths via <i>PathRacer<i>. Then <i>hmmsearch<i> is used to find cry toxins domains in the obtained sequences.In the next step results of <i>hmmsearch<i> are combined to get toxins that posses all three domains. <br> Coordinates of the domains are used to cut flanking sequences and save domains with corresponding linkers. Full sequences (without processing procedure) are used to compare obtained toxins with Bt nomenclature database via <i>diamond blastp<i>. Non-identical sequences are extracted and marked as the potentially new toxins. For all the found sequences (both identical to presented in Bt nomenclature and novel) an online ipg(Identical Protein Group)-annotation is performed (to see the annotation output reed the annotation output section below). Finally, nucleotide sequences, corresponding to protein sequences of the found toxins, are downloaded. Metadata could be uploaded only if accession numbers are present in the quiery.

## Installation and usage
### Prerequisites
<ul>
  <li>python (3.7 or higher) </li>
  <li>Biopython (1.66 or higher)
  </li>
</ul>

To install Biopython use the following command:

```
~$ pip install biopython
```

### Installing

To install CryProcessor clone git repository to your PC.

```
~$ git clone https://github.com/cry_Processor
```
After downloading CryProcessor is ready to use.

## Running and using tool

### Quick usage
To extract cry toxins from the protein fasta file simply execute the following command:
```
~$ python3 cry_processor.py -fi input.faa  -od output_dir
```
This command will automatically search for cry-toxins in fasta-file with amino acid sequences.

### Supported input formats

<ul>
  <li>fasta files with protein sequences </li>
  <li>gfa file (genome assebly graph)</li>
  <li>forward and reverse illumina reads</li>
</ul>

### Tool options: 
Full list of tool options:
```
-fi <input.fasta> or <input.gfa> input file in fasta-format of in gfa-format
-hm <path to hmmer directory> path to hmmemer directory if you want to use local hmmer
-pr <1 or 2> processig type: 1 for extracting all domains, 2 for extrating 2-3 domains only (default 1)
-th nubmer of threads for hmmer/spades/pathracer (default 8)
-ma <e-mail> e-mail address for connecting to NCBI
-od <output dir> output directory
-r <do or fd> working regime: do - domain only search only; fd - full toxins mining with subsequent domain search
-a (--annotate) perform data anotation with ipg (only if accession numbers are present)
-nu <fn or pn or an> upload nucleotide sequences: fn - full sequences, pn -processed sequences, an - both variants
-pa (--pathracer) - launching pathracer on gfa file
-fo <input_1.fastq> - forward illumina reads
-re <input_2.fastq> - reverse illumina reads
-n (--meta) - flag for  specifying metagenomic spades regime
-k <number> - K-mer size for pathracer (default 21)
-s (--silent) - disable console output
```

### Using tool for fasta files
To use tool for fasta file execute the command, presented in quick usage, you can also specify annotation (writing e-nail ):
```
~$ python3 cry_processor.py -fi input.faa  -od output_dir -ma <e-mail address> --annotate
```
Use -nu flag to download nucleotide sequences:

```
~$ python3 cry_processor.py -fi input.faa  -od output_dir -ma <e-mail address> --annotate -nu pn
```
Pipeline of searching could be performed in two regimes:
<ul>
  <li>do - domain only regime. Searches cry-toxin domains from full queiry, then combines this data to extract toxins with 3 domains </li>
  <li>fd - find domains regime. At the begining hmmsearch using full cry toxin model is performed. Domains are extracted next, this regime is quicker as far as domain search is performed on already extracted sequences </li>
</ul>

### Using tool for gfa files
You can apply cry toxins search directly from assembly graph in gfa format with the following commad:

```
~$ python3 cry_processor.py -fi input.gfa  -od output_dir --path_racer
```
### Using tool for illumina reads
This regime includes read assembly with SPades and subsequent hmm-based toxins mining, to impliment this use thw following command:

```
~$ python3 cry_processor.py -fo forward_reads.fastq -re reverse_reads.fastq -od output_dir 
```
If you want to search for cry-toxins from metagenomic reads specify --meta flag
Note that you cannot mix regimes. Do not use --pathracer flag with illumina fastq-files quiery, do not mix -fi agrument with -fo and -re arguments for the correct working.

### Annotation output
Using --annotate flag would perform NCBI-search in ipg(identical protein group) database for submitted accession numbers within quiery and would return information in tsv-format with the following structure:
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
  <li>protein_id - initial sequence id in the quiery </li>
  <li>initial_description - protein description in quiery </li>
  <li>top_cry_hit - cry-protein from btnomenclature with the best identity score according to diamond</li>
  <li>cry_identity - identity score with the closest cry-protein from btnomenclature </li>
  <li>source - database source of the sequence </li>
  <li>nucl_accession - accession number for nucleotide sequence </li>
  <li>start - starting position in full nucleotide sequence </li>
  <li>stop - ending position in full nucleotide sequence </li>
  <li>strand - DNA strand orientation for nucleotide sequence </li>
  <li>ipg_prot_id - protein id according to ipg database </li>
  <li>ipg_prot_name - protein description to ipg database </li>
  <li>organism - taxon description </li>
  <li>strain - strain of the organism </li>
  <li>assembly - assembly id where nucleotide sequence is present</li>
</ul>
Rows possesing first 4 coloums fefer to initial sequences (those found in quiery), results for all identical proteins are marked with --. Here is an expamle of such rows:

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
Note that all rows marked with -- in first four coloumns are identical to initial sequense in queiry with full coloumns located above.

### Output files structure



## Author

* **Anton Shikov** - *Initial work* - [anton-shikov](https://github.com/anton-shikov)


