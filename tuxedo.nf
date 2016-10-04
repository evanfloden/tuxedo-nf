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
log.info "reads                  : ${params.seqs}"
log.info "genome                 : ${params.genome}"
log.info "index                  : ${params.index}"
log.info "annotation             : ${params.annotation}"
log.info "run index              : ${params.run_index}"
log.info "output                 : ${params.output}"
log.info "\n"

/*
 * Input parameters validation
 */

genome_file                   = file(params.genome)
annotation_file               = file(params.annotation) 
index_file                    = file(params.index)
index_file1                   = index_file.fileName() + ".1.ht2"

/*
 * validate input files
 */
if( !genome_file.exists() ) exit 1, "Missing genome file: ${genome_file}"

if( !annotation_file.exists() ) exit 1, "Missing annotation file: ${annotation_file}"

/*
 * Create a channel for read files 
 */
 
Channel
    .fromFilePairs( params.seqs, size: -1 )
    .ifEmpty { error "Cannot find any reads matching: ${params.seqs}" }
    .set { read_files } 



/*
 * Prepare index files if required
 */
if( !params.run_index ) {
    if ( index.file1.exists() ) {
        premade_genome_index = index_file
    }
}
     



// GENOME INDEXING
// ===============

if (params.run_index) {
    process genome_index {
        input:
        file genome_file

        output:
        file "index_dir" into genome_index

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
    process prepare_index {

    }
}    

process mapping {
    tag "reads: $name"

    input:
    file index_dir from genome_index.first()
    set val(name), file(reads) from read_files

    output:
    file "${name}.sam" into hisat2_sams 

    script:
    //
    // HISAT2 mapper
    //
    def single = reads instanceof Path
    if( single ) {
        """
        mv ${index_dir}/* .
        hisat2 -x genome_index -U ${reads} -S ${name}.sam 2> ${name}.alnstats

        """
    }  
    else {
        """
        mv ${index_dir}/* .
        hisat2 -x genome_index -1 ${reads}[0] -2 ${reads}[1] -S ${name}.sam 2> ${name}.alnstats
        """
    }

}
