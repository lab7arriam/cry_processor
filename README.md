# Cry Processor
Cry processor is high-troughtput tool for mining cry-toxins from fasta-files or directly from illumina reads.

## Getting Started

Tool description

### Prerequisites
<ul>
  <li>python (3.7 or higher) </li>
  <li>Biopython (1.66 or higher)
  </li>
</ul>

### Installing

To install this tool clone this repository to your PC.

```
~$ git clone https://github.com/cry_Processor
```

## Running and using tool

## Example usage
```
~$ python3 cry_processor.py -fi input.faa  -od output_dir
```
This command will automatically search for cry-toxins in fasta-file with amino acid sequences


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
-fi <input.fasta>/<input.gfa> Input file in fasta-format of in gfa-format
-hm <path to hmmer directory> Path to hmmemer directory if you want to use local hmmer
-pr 1/2 Processig type: 1 for extracting all domains, 2 for extrating 2-3 domains only
-th Nubmer of threads for hmmer/spades/pathracer
-ma <e-mail> e-mail address for connecting to NCBI
-od <output dir> output directory
-r 

```


## Author

* **Anton Shikov** - *Initial work* - [anton-shikov](https://github.com/anton-shikov)


