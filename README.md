tbAttSeq is Tome's internal pipeline to process Att-Seq (in-vitro recombination with known, short oligos) data. 

**Install Dependencies**

Make sure you have conda/mamba/etc installed on your system

Clone the repo

`cd tbAttSeq`

`conda env create -f environment.yml`

when its finished creating the environment:

`conda activate tb-attseq`

**Samplesheet**

The samplesheet should have the following format, comma separated

| sample_name | fastq_dir              | oligos            | group  |
|-------------|------------------------|-------------------|--------|
| NC          | input/reads/NC*        | input/oligos.csv  | NC     |
| WT          | input/reads/WT*        | input/oligos.csv  | WT     |
| PL1750      | input/reads/PL1750*    | input/oligos.csv  | PL1750 |
| PL1768      | input/reads/PL1768*    | input/oligos.csv  | PL1768 |
| PL2412      | input/reads/PL2412*    | input/oligos.csv  | PL2412 |
| PL2421      | input/reads/PL2421*    | input/oligos.csv  | PL2421 |
| PL2648      | input/reads/PL2648*    | input/oligos.csv  | PL2648 |

**Oligos**

In the samplesheet, one of the columns links to a file containing information about the known oligos (usually attB sequences) and should look like this:

| name    | sequence                                           |
|---------|----------------------------------------------------|
| CAS421  | GTTTGGTTTGTTTGCAACGGCAGTGACGGAGGTTGGGAGCCAGGCT   |
| CAS008  | TCCAGGCTTGTCCACAGTTGTAGTCTCCATGGGAGAAGCAGCTTCT   |
| CAS023  | CCGTGCAGCTTCAGCCTCTGCCGTCTCAGCTGAGCAAGCAGGAGGC   |
| CAS007  | CTTCTGCTTGCTGACAGCTGAGGTTCCTATTGCCAAGATCTCTGCA   |
| CAS004  | AGAAACAGCTCCAACTTCAGAGGTTGCTGCAGGATAAACTCAGAAA   |
| CAS1056 | CAACAGTTTGTGTGTAACAGGCGTTACAGTGGGGAGAAGCCAGGGTC |

How to run:

nextflow run pipeline.nf -c input/nextflow.config
