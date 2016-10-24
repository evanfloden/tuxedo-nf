# Tuxedo-NF

A Nextflow implementation of the Tuxedo Suite of Tools Workflow is based on the 2016 *Nature Protocols* publication: ["Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown"](http://www.nature.com/nprot/journal/v11/n9/full/nprot.2016.095.html)

[![nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.22.0-brightgreen.svg)](http://nextflow.io)

## Quick start 

Make sure you have all the required dependencies listed in the last section.

Install the Nextflow runtime by running the following command:

    $ curl -fsSL get.nextflow.io | bash


When done, you can launch the pipeline execution by entering the command shown below:

    $ nextflow run skptic/tuxedo-nf
    

By default the pipeline is executed against the provided example dataset. 
Check the *Pipeline parameters* section below to see how enter your data on the program 
command line.

All parameters can be specified at the command line or alternatively specified in a parameters config file. 

Default parameters can be found in the params_default.config file.     
    

## Pipeline parameters

#### `--reads` 
   
* Specifies the local location of the reads *fastq* file(s).
* Multiple files can be specified using the usual wildcards (*, ?), in this case make sure to surround the parameter string
  value by single quote characters (see the example below)
* It must end in `.fastq` or fastq.gz.
* See `--sra_ids` and `--use_sra` to pull reads directly from the NCBI SRA.
* Involved in the task mapping.
* By default it is set to the Tuxedo-NF's location: `./example-data/reads/SRR*_*_{1,2}.fastq.gz`


Example: 

    $ nextflow run skptic/tuxedo-nf --reads '/home/dataset/*.fastq'

This will handle each fastq file as a seperate sample.

Read pairs of samples can be specified using the glob file pattern. Consider a more complex situation where there are three samples (A, B and C), with A and B being paired reads and C being single ended. The read files could be:
    
    sample_A_1.fastq
    sample_A_2.fastq
    sample_B_1.fastq
    sample_B_2.fastq 
    sample_C_1.fastq

The reads may be specified as below:

    $ nextflow run skptic/tuxedo-nf --reads '/home/dataset/sample_*_{1,2}.fastq'    


  
#### `--genome`

* The location of the genome multi-fasta file.
* It should end in `.fa`
* See `--genome_address` and `--download_genome` to pull a genome directly.
* Involved in the task: genome_index.
* By default it is set to the Tuxedo-NF's localization: `./example-data/genome/genome.fa`

Example:

    $ nextflow run skptic/tuxedo-nf --genome /home/user/my_genome/example.fa


#### `--index`

* The location of a HISAT2 Index.
* It should point to the location of the index (without the `.X.ht2` suffix)
* See `--run_index` to have the index generated from a genome.
* Involved in the task: genome_index.
* By default it is set to the Tuxedo-NF's localization: `./example-data/index/genome`

Example:

    $ nextflow run skptic/tuxedo-nf --genome /home/user/my_genome_index/example



#### `--pheno`

* The location of the phenotype/group description file.
* It should end in `.csv` or `.txt`.
* This file defines the groups that each sample belong to.
* ids should be the sample name, for example the SRA id or the name of the fastq (without `_1.fastq` the suffix)
* Usually a tab or comma delimited file, for example:
	
	ids	groups
	SRR	control
	SRR	control
	SRR	disease
	SRR	disease

* Involved in the task: ballgown.
* By default it is set to the Tuxedo-NF's localization: `./example-data/pheno/exp-info.txt`

Example:

    $ nextflow run skptic/tuxedo-nf --transcriptome /home/user/my_exp/exp-info.txt


#### `--download_genome`

* Boolean value [ true || false ]
* If `true`, the file specified by `--genome_address` will be downloaded and used for subsequent tasks
* Involved in the tasks: download_genome
* By default `--download_genome=false`

Example:

    $ nextflow run skptic/tuxedo-nf --download-genome=true 
or equivalently just
    $ nextflow run skptic/tuxedo-nf --download-genome


#### `--genome_address`

* http or ftp address specifying a genome in fasta format.
* It should end in `.gz` or `.tar.gz`
* Sources for genomes include for [ensembl](ftp://ftp.ensembl.org/pub/release-86/fasta/) and [UCSC](ftp://hgdownload.cse.ucsc.edu/goldenPath/currentGenomes/).
* If `--download_genome` is `true` then the genome will be downloaded from the above address and used for subsequent tasks
* Involved in the tasks: download_genome
* By default `--genome_address=ftp://ftp.ensembl.org/pub/release-86/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz`

Example:

    $ nextflow run skptic/tuxedo-nf --genome_address=ftp://hgdownload.cse.ucsc.edu/goldenPath/currentGenomes/Homo_sapiens/bigZips/chromFa.tar.gz


#### `--download_annotation`

* Boolean value [ true || false ]
* If `true`, the file specified by `--annotation_address` will be downloaded and used for subsequent tasks
* Involved in the tasks: download_annotation 
* By default `--download_annotation=false`

Example:

    $ nextflow run skptic/tuxedo-nf --download_annotation=true
or equivalently just
    $ nextflow run skptic/tuxedo-nf --download_annotation


#### `--annotation_address`

* http or ftp address specifying an annotation gtf format
* It should end in `.gz` or `tar.gz`.
* Sources for genomes include for [ensembl](ftp://ftp.ensembl.org/pub/release-86/gtf/) and [UCSC](http://hgdownload.cse.ucsc.edu/downloads.html).
* If `--download-annotation` is `true` then the annotation will be downloaded from the above address and used for subsequent tasks
* Involved in the tasks: download_annotation
* By default `--annotation_address=ftp://ftp.ensembl.org/pub/release-86/gtf/homo_sapiens/Homo_sapiens.GRCh38.86.gtf.gz`

Example:

    $ nextflow run skptic/tuxedo-nf --genome_address=ftp://hgdownload.cse.ucsc.edu/goldenPath/currentGenomes/Homo_sapiens/bigZips/chromFa.tar.gz


#### `--run_index`

* Boolean value [ true || false ]
* If `true`, and the file specified by `--genome` or `--genome_address` will be indexed and used for subsequent tasks
* If `--download_genome` is true, then `--run_index` must be true
* Involved in the tasks: genome_index
* By default `--run_index=true`

Example:

    $ nextflow run skptic/tuxedo-nf --run_index=false


#### `--use_sra`

* Boolean value [ true || false ]
* If `true`, the SRA IDs specified by `--sra_ids` will be prefectched and used for subsequent tasks
* Involved in the tasks: sra_prefetch
* By default `--use_sra=false`

Example:

    $ nextflow run skptic/tuxedo-nf --use_sra=true
or equivalently just
    $ nextflow run skptic/tuxedo-nf --use_sra


#### `--sra_ids`

* Comma seperated list of SRA reads using *Run* accession id (ussually SRR or HRR)
* See the [NCBI Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra/) to more information or to browse data.
* If `--use-sra` is `true` then the SRA ids listed will be prefetched from either the local cache (see `--cache`) or downloaded.
* Involved in the tasks: sra_prefetch
* By default `--sra_ids="ERR188044,ERR188104,ERR188234,ERR188245,ERR188257,ERR188273,ERR188337,ERR188383,ERR188401,ERR188428,ERR188454,ERR204916"`

Example:

    $ nextflow run skptic/tuxedo-nf --sra_ids=`SRR349706,SRR349707,SRR349708`


#### `--cache`

* Location of NCBI cache. 
* If `--use-sra` is `true` then the SRA ids listed will be prefetched and stored in the `--cache` location.
* It should contain a directory called `ncbi` which contains the sra-tools/vdb database. 
* Can be useful if several pipelines use the same input sequences, saving on storage and bandwidth.
* Involved in the tasks: sra_prefetch
* By default `--cache=./cache`

Example:

    $ nextflow run skptic/tuxedo-nf --cache=`/your/ncbi_cache_location`



#### `--output` 
   
* Specifies the folder where the results will be stored for the user.  
* It does not matter if the folder does not exist.
* By default is set to Tuxedo-NF's folder: `./results` 

Example: 

    $ nextflow run skptic/tuxedo-nf --output /home/user/my_results 
  


## Cluster support

Tuxedo-NF execution relies on [Nextflow](http://www.nextflow.io) framework which provides an 
abstraction between the pipeline functional logic and the underlying processing system.

Thus it is possible to execute it on your computer or any cluster resource
manager without modifying it.

Currently the following platforms are supported:

  + Oracle/Univa/Open Grid Engine (SGE)
  + Platform LSF
  + SLURM
  + PBS/Torque


By default the pipeline is parallelized by spanning multiple threads in the machine where the script is launched.

To submit the execution to a SGE cluster create a file named `nextflow.config`, in the directory
where the pipeline is going to be launched, with the following content:

    process {
      executor='sge'
      queue='<your queue name>'
    }

In doing that, tasks will be executed through the `qsub` SGE command, and so your pipeline will behave like any
other SGE job script, with the benefit that *Nextflow* will automatically and transparently manage the tasks
synchronisation, file(s) staging/un-staging, etc.

Alternatively the same declaration can be defined in the file `$HOME/.nextflow/config`.

To lean more about the avaible settings and the configuration file read the 
[Nextflow documentation](http://www.nextflow.io/docs/latest/config.html).
  
  
Dependencies 
------------

 * Java 7+ 
 * Nextflow (0.22.0 or higher)
 * Docker
