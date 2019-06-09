# CryProcessor
CryProcessor is a high-troughtput tool for the cry-toxins mining from fasta-files or directly from illumina reads.

## About CryProcessor

CryProcessor is a python-written tool for searching and extracting cry toxins from illumina sequence data or from protein fasta files. It includes several parts: hmm-based searching of potential cry-toxins, obtaining information about domains, extracting cry-toxins only with 3 domains and comparing found toxins with BT nomenclature. <br> <br> A regime for the performing toxins search directly from illumina reads implies building an assembly graph (using  <i> SPAdes </i>) and the subsequent mining toxins directly from the obtained assebmly graph. 

## CryProcessor Pipeline

The following text stands for the full pipeline description (for illumina reads). To start, <i>SPAdes</i> or <i>metaSPAdes</i> are implemented to get assembly graph from the fastq-files. After that, potential cry toxins (with at least 30% identity to the hmm-consensus) are extracted from assembly paths via <i>PathRacer</i>. Then <i>hmmsearch</i> is used to find cry toxins domains in the obtained sequences. In the next step, results of <i>hmmsearch</i> are combined to get toxins that posses all three domains. <br> <br>Coordinates of the domains are used to cut flanking sequences and save domains with corresponding linkers. Full sequences (without processing procedure) are used to compare obtained toxins with Bt nomenclature database via <i>diamond blastp</i>. Non-identical sequences are extracted and marked as the potentially new toxins. <br><br>For all the found sequences (both identical to presented in Bt nomenclature and novel sequences) an online ipg-annotation (Identical Protein Group) is performed (to see the annotation output reed the annotation output section below). Finally, nucleotide sequences, corresponding to the protein sequences of the found toxins, are downloaded. Metadata will be uploaded only if accession numbers are present in the quiery.

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
~$ git clone https://github.com/cry_Processor
```
After downloading, CryProcessor is ready to use.

## Running and Using Tool

### Quick Usage
To extract cry toxins from the protein fasta file simply execute the following command:
```
~$ python3 cry_processor.py -fi input.faa  -od output_dir
```
This command will automatically search for cry-toxins in the fasta file with amino acid sequences.

### Supported Input Formats

<ul>
  <li>fasta files with protein sequences; </li>
  <li>gfa file (genome assebly graph);</li>
  <li>forward and reverse illumina reads.</li>
</ul>

### Tool Options: 
The full list of tool options:
```
-fi <input.fasta> or <input.gfa> an input file in fasta format of in gfa format
-hm <path to hmmer directory> path to the hmmer directory if you want to use local hmmer
-pr <1 or 2> processig type: 1 for extracting all the domains, 2 for extrating 2-3 domains only (default 1)
-th nubmer of threads for hmmer/spades/pathracer (default 8)
-ma <e-mail> e-mail address for connecting to NCBI
-od <output dir> the output directory
-r <do or fd> working regime: do - domain only search only; fd - full toxins mining with subsequent domain search
-a (--annotate) perform data anotation with ipg (only if accession numbers are present)
-nu <fn or pn or an> upload nucleotide sequences: fn - full sequences, pn - processed sequences, an - both variants
-pa (--pathracer) - launching PathRacer on gfa file
-fo <input_1.fastq> - forward illumina reads
-re <input_2.fastq> - reverse illumina reads
-n (--meta) - flag for specifying metagenomic SPAdes regime
-k <number> - K-mer size for PathRacer (default 21)
-s (--silent) - disable console output
```

### Using Tool for Fasta Files
To use tool for the files in fasta format execute the command, presented in quick usage, you can also specify annotation (writing e-mail is strongly recommended):
```
~$ python3 cry_processor.py -fi input.faa  -od output_dir -ma <e-mail address> --annotate
```
Use <i>-nu</i> flag to download nucleotide sequences:

```
~$ python3 cry_processor.py -fi input.faa  -od output_dir -ma <e-mail address> --annotate -nu pn
```
The pipeline of searching could be performed in two regimes:
<ul>
  <li><i>do</i> - domain only regime. Searches cry-toxin domains from the full queiry, then combines this data to extract toxins with 3 domains; </li>
  <li><i>fd</i> - find domains regime. At the begining, <i>hmmsearch</i> with the full cry-toxin model is performed. Next, domains are extracted and combined. this regime is quicker, as far as the domain search is performed on previously extracted sequences. </li>
</ul>

### Using Tool for .gfa Files
You can apply the cry toxins search directly from the assembly graph in gfa format with the following commad:

```
~$ python3 cry_processor.py -fi input.gfa  -od output_dir --path_racer
```
### Using Tool for Illumina Reads
This regime includes the reads assembly with SPAdes and the subsequent hmm-based toxins mining. To impliment this use the following command:

```
~$ python3 cry_processor.py -fo forward_reads.fastq -re reverse_reads.fastq -od output_dir 
```
If you want to search for cry-toxins from metagenomic reads specify <i>--meta</i> flag.
Note that you cannot mix regimes. <br>Do not use the <i>--pathracer</i> flag with the illumina fastq files quiery (<i>-fo</i> and <i>-re</i>), do not mix the <i>-fi</i> agrument with <i>-fo</i> and <i>-re</i> arguments for the correct running.

### Annotation Output
Using the <i>--annotate</i> flag will perform the NCBI-search in the ipg database for submitted accession numbers within the quiery and return gathered information in tsv-format with the following structure:
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
  <li>initial_description - protein description in the quiery </li>
  <li>top_cry_hit - cry-protein from BT nomenclature with the best identity score according to <i>diamond</i></li>
  <li>cry_identity - identity score with the closest cry-protein from btnomenclature </li>
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
Rows possesing the first four coloums refer to the initial sequences (which are found in the quiery), results for all the identical proteins are marked with --. <br>Here's an expamle of such rows:

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
Note that all the rows, marked with -- in the first four coloumns, are identical to the initial sequense in the queiry with the full coloumns, located above.

### Output Files Structure



## Author

* **Anton Shikov** - *Initial work* - [anton-shikov](https://github.com/anton-shikov)


