/*
 * Copyright (c) 2016, Centre for Genomic Regulation (CRG) and the authors.
 *
 *   This file is part of 'Tuxedo-NF'.
 *
 *   Tuxedo-NF is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Tuxedo-NF is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Tuxedo-NF.  If not, see <http://www.gnu.org/licenses/>.
 */

/* 
 * Main Tuxedo-NF pipeline script
 *
 * @authors
 * Evan Floden <evanfloden@gmail.com> 
 */

log.info "T U X E D O - N F  ~  version 0.1"
log.info "====================================="
log.info "reads                  : ${params.reads}"
log.info "genome                 : ${params.genome}"
log.info "index                  : ${params.index}"
log.info "phenotype info         : ${params.pheno}"               
log.info "annotation             : ${params.annotation}"
log.info "download genome        : ${params.download_genome}"
log.info "download annotation    : ${params.download_annotation}"
log.info "UCSC annotation        : ${params.UCSC_annotation}"
log.info "run index              : ${params.run_index}"
log.info "use SRA                : ${params.use_sra}"
log.info "SRA ids                : ${params.sra_ids}"
log.info "NCBI reads cache       : ${params.cache}"
log.info "output                 : ${params.output}"
log.info "\n"

/*
 * Input parameters validation
 */

annotation_file               = file(params.annotation)
genome_file                   = file(params.genome)
UCSC_annotation               = params.UCSC_annotation.toString().toUpperCase()
pheno_file                    = file(params.pheno)
index_file                    = file(params.index)
index_file1                   = index_file + ".1.ht2"
index_name                    = index_file.getFileName()
index_dir                     = index_file.getParent()
sra_ids_list                  = params.sra_ids.tokenize(",") 
cache                         = file(params.cache)


/*
 * validate and create a channel for genome/index input files
 */

if( !params.download_genome && !genome_file.exists() && params.run_index && !index_file1.exists() ) 
	exit 1, "Missing genome file: ${genome_file}"

if( !params.download_genome && !index_file1.exists() && !params.run_index) 
        exit 1, "Missing index file: ${index_file}"

if( params.download_genome ) {
    process download_genome {

            input:
            val (params.genome_address)

            output:
            file "*.fa" into genomes

            script:
            //
            // Genome Download
            //
            """
            wget ${params.genome_address}
            gunzip -d *.gz
            """
    }
} else { 
    Channel 
        .fromPath { genome_file}
        .set { genomes }
}


/*
 * validate and create a channel for annotation input files
 */

if( params.download_annotation ) {
  process download_annotation {

            input:
            val (params.annotation_address)

            output:
            file "*.gtf" into annotations1, annotations2, annotations3, annotations4

            script:
            //
            // Annotation Download
            //
            """
            wget ${params.annotation_address}
            gunzip -d *.gz
            """
    }
} else {    
    Channel    
        .fromPath ( annotation_file )
        .into { annotations1; annotations2; annotations3; annotations4}
}


/*
 * Create a channel for read files 
 */

if ( !params.use_sra ) {
    Channel
        .fromFilePairs( params.reads, size: -1 , flat: true)
        .ifEmpty { error "Cannot find any reads matching: ${params.seqs}" and use_sra is false}
        .set { read_files } 
}

/*
 * Create a channel for SRA IDs
 */

Channel
    .from( sra_ids_list )
    .set { sra_read_ids }


// GENOME INDEXING
// ===============

if (params.run_index) {
    process genome_index {
        input:
        file genome_file from genomes

        output:
        set val ("genome_index"), file("index_dir") into genome_index

        script:
        //
        // HISAT2 genome index
        //
        """
        hisat2-build ${genome_file} genome_index
        mkdir index_dir
        mv genome_index* index_dir/.
        """
    }
}
else {
    process premade_index {
        input:
        file(index_dir)
        val(index_name)

        output:
        set val(index_name), file("index_dir") into genome_index

        script:
        //
        // Premade HISAT2 genome index
        //
        """
        mkdir index_dir
        cp ${index_dir}/${index_name}.*.ht2 index_dir/.
        """
    }
}    



if (params.use_sra) {
    process sra_prefetch {

        publishDir = [path: {params.cache}, mode: 'copy', overwrite: 'true' ]
        tag "sra_id: $sra_id"

        input:
        val(sra_id) from sra_read_ids

        output:
        file "ncbi/**" optional true into sra_cache_elements
        val (sra_id) into prefetched_sras

        script:
        //
        // SRA Cache Check and Download
        //
        """
        prefetch -a "/home/sra_user/.aspera/connect/bin/ascp|/home/sra_user/.aspera/connect/etc/asperaweb_id_dsa.openssh" -t fasp ${sra_id}        
        """
    }
}

if (params.use_sra) {  
    process sra_mapping {
        tag "reads: $sra_id"

        input:
        set val (index_name), file(index_dir) from genome_index.first()
        val(sra_id) from prefetched_sras

        output:
        set val(sra_id), file("${sra_id}.sam") into hisat2_sams

        script:
        //
        // HISAT2 mapper using SRAToolkit with ncbi-vdb support
        //
        """
        fastq-dump --split-files ${sra_id}
        hisat2 -x ${index_dir}/${index_name} -1 ${sra_id}_1.fastq -2 ${sra_id}_2.fastq -S ${sra_id}.sam
        """
    }
}
else {
    process mapping {
        tag "reads: $name"

        input:
        set val(index_name), file(index_dir) from genome_index.first()
        set val(name), file(reads) from read_files

        output:
        set val(name), file("${name}.sam") into hisat2_sams

        script:
        //
        // HISAT2 mapper
        //
        def single = reads instanceof Path
        if( single ) 
            """
            hisat2 -x ${index_dir}/${index_name} -U ${reads[0]} -S ${name}.sam
            """
        else 
            """
            hisat2 -x ${index_dir}/${index_name} -1 ${reads[0]} -2 ${reads[1]} -S ${name}.sam
            """
    }
}


process sam2bam {
    tag "sam2bam: $name"

    input:
    set val(name), file(sam) from hisat2_sams

    output:
    set val(name), file("${name}.bam") into hisat2_bams

    script:
    //
    // SAM to sorted BAM files
    //
    """   
    samtools view -S -b ${sam} | samtools sort -o ${name}.bam -    
    """
}


hisat2_bams.into { hisat2_bams1; hisat2_bams2 } 

process stringtie_assemble_transcripts {
    tag "stringtie assemble transcripts: $name"

    input:
    set val(name), file(bam) from hisat2_bams1
    file (annotation_f) from annotations1.first()

    output:
    file("${name}.gtf") into hisat2_transcripts

    script:
    //
    // Assemble Transcripts per sample
    //
    """   
    stringtie -p ${task.cpus} -G ${annotation_f} -o ${name}.gtf -l ${name} ${bam}
    """
}   


hisat2_transcripts.into { hisat2_transcripts1; hisat2_transcripts2 }

hisat2_transcripts2
    .toList()
    .set { grouped_transcripts }

hisat2_transcripts1
  .collectFile () { file ->  ['gtf_filenames.txt', file.name + '\n' ] }
  .set { GTF_filenames }

process merge_stringtie_transcripts {
    tag "merge stringtie transcripts"

    input:
    file (merge_list) from GTF_filenames
    file (gtfs) from grouped_transcripts
    file (annotation_f) from annotations2

    output:
    file("stringtie_merged.gtf") into merged_transcripts

    script:
    //
    // Merge all stringtie transcripts
    //
    """
    stringtie --merge  -p ${task.cpus}  -G ${annotation_f} -o stringtie_merged.gtf  ${merge_list}
    """
}

merged_transcripts.into { merged_transcripts1; merged_transcripts2 }

process transcript_abundance {
    tag "reads: $name"

    input:
    set val(name), file(bam) from hisat2_bams2
    file merged_transcript_file from merged_transcripts1.first()

    output:
    file("${name}") into ballgown_data

    script:
    //
    // Estimate abundances of merged transcripts in each sample
    //
    """
    stringtie -e -B -p ${task.cpus} -G ${merged_transcript_file} -o ${name}/${name}_abundance.gtf ${bam}
    """
}

process gffcompare {
    tag "gffcompare"
    
    input:
    file (annotation_f) from annotations3
    file (merged_transcripts) from merged_transcripts2

    output:
    file("merged_gffcompare") into gffcompare

    shell:
    //
    // Compare merged stringtie transcripts with annotation
    //
    """
    gffcompare -r ${annotation_f} -G -o compare_merged ${merged_transcripts}
    mkdir merged_gffcompare
    mv compare_merged* merged_gffcompare/.
    """
}	


process ballgown {
    tag "ballgown"

    input:
    file pheno_file
    file annotation_f from annotations4
    file ballgown_dir from ballgown_data.toList()

    output:
    file("genes_results.csv") into sig_genes
    file("transcripts_results.csv") into sig_transcripts


    shell:
    //
    // Merge all stringtie transcripts
    //
    '''
    #!/usr/bin/env Rscript
 
    pheno_data_file <- "!{pheno_file}"

    library(ballgown)
    library(RSkittleBrewer)
    library(genefilter)
    library(dplyr)
    library(devtools)

    pheno_data <- read.table(pheno_data_file, header=TRUE, colClasses = c("factor","factor","factor"))

    bg <- ballgown(dataDir = ".", samplePattern="ERR", pData=pheno_data)

    # Get Gene Symbols of Transcripts 
    gene_symbols <- getGenes("!{annotation_f}", structure(bg)$trans, UCSC=!{UCSC_annotation}, attribute="gene_name")
    vector1 <- sapply(gene_symbols, function(x){as.vector(x)})
    gene_symbols_vector <- sapply(vector1, function(x){toString(x)})
    expr(bg)$trans$gene_name = gene_symbols_vector

    results_transcripts <-  stattest(bg, feature='transcript', covariate='sex',
                            getFC=TRUE, meas='FPKM')

    results_genes <-  stattest(bg, feature='gene', covariate='sex',
                      getFC=TRUE, meas='FPKM')

    results_transcripts <- data.frame(geneName=ballgown::geneNames(bg),
                           geneID=ballgown::geneIDs(bg), results_transcripts)

    results_transcripts <- arrange(results_transcripts, pval)
    results_genes <-  arrange(results_genes, pval)

    write.csv(results_transcripts, "transcripts_results.csv", row.names=FALSE)
    write.csv(results_genes, "genes_results.csv", row.names=FALSE)

    '''
}


