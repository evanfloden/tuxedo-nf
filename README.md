# Tuxedo-NF

A Nextflow implementation of the Tuxedo Suite of Tools

Based on Nature Protocols: ["Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown"](http://www.nature.com/nprot/journal/v11/n9/full/nprot.2016.095.html)

[![nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.22.0-brightgreen.svg)](http://nextflow.io)

## Quick start 

Download the example dataset using the following commands:

    $ wget -O- http://genome.crg.es/~efloden/tuxexdo-nf/example-data.tar.gz | tar xzf


Make sure you have all the required dependencies listed in the last section.

Install the Nextflow runtime by running the following command:

    $ curl -fsSL get.nextflow.io | bash


When done, you can launch the pipeline execution by entering the command shown below:

    $ nextflow run cbcrg/tuxedo-nf
    

By default the pipeline is executed against the provided example dataset. 
Check the *Pipeline parameters*  section below to see how enter your data on the program 
command line.     
    


## Pipeline parameters

#### `--reads` 
   
* Specifies the location of the reads *fastq* file(s).
* Multiple files can be specified using the usual wildcards (*, ?), in this case make sure to surround the parameter string
  value by single quote characters (see the example below)
* It must end in `.fastq` or fastq.gz.
* Involved in the task: mapping.
* By default it is set to the Tuxedo-NF's location: `./example-data/ERR*_*_{1,2}.fastq.gz`

Example: 

    $ nextflow run cbcrg/tuxedo-nf --reads '/home/dataset/*.fastq'

This will handle each fastq file as a seperate sample.

Read pairs of samples can be specified using the glob file pattern. Consider a more complex situation where there are three samples (A, B and C), with A and B being paired reads and C being single ended. The read files could be:
    
    sample_A_1.fastq
    sample_A_2.fastq
    sample_B_1.fastq
    sample_B_2.fastq 
    sample_C_1.fastq

The reads may be specified as below:

    $ nextflow run cbcrg/tuxedo-nf --reads '/home/dataset/sample_*_{1,2}.fastq'    

  
#### `--genome`

* The location of the genome multi-fasta file.
* It should end in `.fa`
* Involved in the task: index.
* By default it is set to the Tuxedo-NF's localization: `./example-data/genome/chrX.fa`

Example:

    $ nextflow run cbcrg/tuxedo-nf --transcriptome /home/user/my_genome/example.fa


#### `--output` 
   
* Specifies the folder where the results will be stored for the user.  
* It does not matter if the folder does not exist.
* By default is set to Tuxedo-NF's folder: `./results` 

Example: 

    $ nextflow run cbcrg/tuxedo-nf --output /home/user/my_results 
  


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
