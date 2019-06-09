# CryProcessor
CryProcessor is a high-troughtput tool for the Cry toxins mining from the fasta-files or directly from the illumina reads.

## About CryProcessor

CryProcessor is a python-written tool for searching and extracting Cry toxins from illumina sequence data or from the protein fasta files. It includes several parts: an hmm-based searching of potential Cry toxins, obtaining information about domains, extracting Cry toxins only with 3 domains and comparing found toxins with BT nomenclature. <br> <br> A regime for  performing the toxins search directly from the illumina reads implies building an assembly graph (using  <i> SPAdes </i>) and the subsequent mining toxins directly from the obtained assebmly graph. 

## CryProcessor Pipeline

The following text stands for the full pipeline description (for the illumina reads). To start, <i>SPAdes</i> (http://cab.spbu.ru/software/spades/) or <i>metaSPAdes</i> (http://cab.spbu.ru/software/meta-spades/) are implemented to get the assembly graph from the fastq-files. After that, the potential Cry toxins (with at least 30% identity to the hmm-consensus) are extracted from the assembly paths via <i>PathRacer</i> (http://cab.spbu.ru/software/pathracer/). Then <i>hmmsearch</i> (http://hmmer.org/) is used to find Cry toxin domains in the obtained sequences. In the next step, the results of <i>hmmsearch</i> are combined to get the toxins that posses all three domains. <br> <br>The coordinates of the domains are used to cut flanking sequences and save the domains with the corresponding linkers. The full sequences (without processing procedure) are used to compare the obtained toxins with Bt nomenclature database via <i>diamond blastp</i> (https://github.com/bbuchfink/diamond). Non-identical sequences are extracted and marked as the potentially new toxins. <br><br>For all the found sequences (both identical to presented in Bt nomenclature and novel sequences) an online ipg-annotation (Identical Protein Group) is performed (to see the annotation output reed the annotation output section below). Finally, nucleotide sequences, corresponding to the protein sequences of the found toxins, are downloaded. Metadata will be uploaded only if the accession numbers are present in the quiery.

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
To extract Cry toxins from the protein fasta file simply execute the following command:
```
~$ python3 cry_processor.py -fi input.faa  -od output_dir
```
This command will automatically search for the Cry toxins in the fasta file with amino acid sequences.

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
-pr <1 or 2> processig type: 1 for extracting all the domains, 2 for extrating 2-3 domains only (default 1)
-th nubmer of threads for hmmer/SPAdes/PathRacer (default 8)
-ma <e-mail> e-mail address for connecting to NCBI
-od <output dir> the output directory
-r <do or fd> the working regime: do - domain only search only; fd - the full toxins mining with subsequent domain search
-a (--annotate) perform data anotation with ipg (only if accession numbers are present)
-nu <fn or pn or an> upload nucleotide sequences: fn - full sequences, pn - processed sequences, an - both variants
-pa (--pathracer) - launching PathRacer on the gfa file
-fo <input_1.fastq> - forward illumina reads
-re <input_2.fastq> - reverse illumina reads
-n (--meta) - the flag for specifying the metagenomic regime forSPAdes
-k <number> - the K-mer size for PathRacer (default 21)
-s (--silent) - disable console output
```

### Using Tool for Fasta Files
To use the tool for the files in the fasta format execute the command, presented in quick usage, you can also specify the annotation (writing an e-mail address is strongly recommended):
```
~$ python3 cry_processor.py -fi input.faa  -od output_dir -ma <e-mail address> --annotate
```
Use the <i>-nu</i> flag to download nucleotide sequences:

```
~$ python3 cry_processor.py -fi input.faa  -od output_dir -ma <e-mail address> --annotate -nu pn
```
The pipeline of searching could be performed in two regimes:
<ul>
  <li><i>do</i> - domain only regime. Searches the Cry toxin domains from the full queiry, then combines this data to extract the toxins with 3 domains; </li>
  <li><i>fd</i> - find domains regime. At the begining, <i>hmmsearch</i> with the full Cry toxin model is performed. Next, the domains are extracted and combined. This regime is quicker, as far as the domain search is performed on previously extracted sequences. </li>
</ul>

### Using Tool for .gfa Files
You can apply the Cry toxins search directly from the assembly graph in the gfa format with the following commad:

```
~$ python3 cry_processor.py -fi input.gfa  -od output_dir --path_racer
```
### Using Tool for Illumina Reads
This regime includes the reads assembly with SPAdes and the subsequent hmm-based toxins mining. To implement this use the following command:

```
~$ python3 cry_processor.py -fo forward_reads.fastq -re reverse_reads.fastq -od output_dir 
```
If you want to search for the Cry toxins from metagenomic reads specify <i>--meta</i> flag.
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
  <li>top_cry_hit - Cry protein from BT nomenclature with the best identity score according to <i>diamond</i></li>
  <li>cry_identity - identity score with the closest Cry protein from btnomenclature </li>
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


## References

1.  Bankevich A., Nurk S., Antipov D., Gurevich A., Dvorkin M., Kulikov A. S., Lesin V., Nikolenko S., Pham S., Prjibelski A., Pyshkin A., Sirotkin A., Vyahhi N., Tesler G., Alekseyev M. A., Pevzner P. A. SPAdes: A New Genome Assembly Algorithm and Its Applications to Single-Cell Sequencing. Journal of Computational Biology, 2012 
2.  Buchfink B.,  Xie C. &  Huson D.H. Fast and Sensitive Protein Alignment using DIAMOND. Nature Methods, 2015
3. Nurk S., Meleshko D., Korobeynikov A., Pevzner P. A. metaSPAdes: a new versatile de novo metagenomics assembler. Genome Research, 2017 
4. Finn R.D., Clements J.,  Eddy S.R. HMMER Web Server: Interactive Sequence Similarity Searching. Nucleic Acids Research, 2011.

## Author

* **Anton Shikov** - *Initial work* - [anton-shikov](https://github.com/anton-shikov)


