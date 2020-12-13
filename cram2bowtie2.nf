#!/usr/bin/env nextflow

/*parameters that are specified at the command line or via config file*/
params.newGenome                = false          /*genome fasta file, must specify complete path. Required parameters*/
params.originalGenome        = false          /*genome fasta file, must specify complete path. Required parameters*/
params.samplePath            = false          /*input folder, must specify complete path. Required parameters*/
params.outputDir             = false          /*output folder, must specify complete path. Required parameters*/
params.singleEnd             = false          /*options: true|false. true = the input type is single end reads; false = the input type is paired reads. Default is false*/

/*output folder paths*/
readPrepPath                 = "${params.outputDir}/read_prep"
alnPath                      = "${params.outputDir}/alns"

/*cluster parameters */
myExecutor                   = 'slurm'
params.myQueue               = 'normal'
defaultCPU                   = '1'
defaultMemory                = '20'
assemblerCPU                 = '12'
assemblerMemory              = '100'

/*software stack*/
params.bowtieMod             = 'Bowtie2/2.4.2-IGB-gcc-8.2.0'
params.samtoolsMod           = 'SAMtools/1.10-IGB-gcc-8.2.0'

/*Prepare input*/
orig_genome_file                  = file(params.originalGenome)
orig_genomeStore                  = orig_genome_file.getParent()
new_genome_file                  = file(params.newGenome)
new_genomeStore                  = new_genome_file.getParent()

// if( !genome_file.exists() ) exit 1, "Missing reference genome file: ${genome_file}"
CRAM_Ch1 = Channel.fromFilePairs("${params.samplePath}", size: 1)

// TEMP Workaround to test bowtie2 step, uncomment above and below processes when fixed
// fq_pe_ch = Channel.fromFilePairs("${params.samplePath}")

/*

  prepare_genome 
  This process is executed only once. This is the original reference genome from the CRAM file

*/


process prepare_genome{
    tag                    { genome }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   defaultCPU
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 params.samtoolsMod
    storeDir               orig_genomeStore
    validExitStatus        0
    
    input:
    file genome from orig_genome_file

    output:
    file "*.fai" into orig_genome_index_ch
    
    script:
    """
    samtools faidx ${genome}
    """
}

/*

  prepare_new_genome 
  This process is executed only once; this is the genome to be mapped to

*/


process prepare_new_genome{
    tag                    { genome }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   defaultCPU
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 params.samtoolsMod
    storeDir               new_genomeStore
    validExitStatus        0
    
    input:
    file genome from new_genome_file

    output:
    file "*.fai" into new_genome_index_ch
    
    script:
    """
    samtools faidx ${genome}
    """
}

/*

  prepare_bowtie2 
  This process is executed only once.  Generate new bowtie2 index if needed

*/


process prepare_bowtie2 {
    tag                    { genome }
    executor               myExecutor
    cpus                   defaultCPU
    memory                 '72 GB'
    cpus                   24
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 params.bowtieMod
    storeDir               new_genomeStore
    validExitStatus        0
    
    input:
    file genome from new_genome_file

    output:
    file "*.bt2" into bowtie2Index
    
    script:
    """
    bowtie2-build --threads ${task.cpus} ${genome} ${genome.getSimpleName()}
    """
}

/*

  qc_input 
  This process mainly checks the inputs for improperly formed CRAM/BAM input files

*/


process qc_input {
    tag                    { name }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   defaultCPU
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 params.samtoolsMod
    // publishDir             readPrepPath, mode: "copy"	    
    validExitStatus        0,1
    errorStrategy          'ignore'
    stageOutMode           'copy'
    
    input:
    set val(name), file(CRAM) from CRAM_Ch1	

    output:
    set val(name), file('*_ok.cram') optional true into extract_ch
    
    script:
    """
    samtools quickcheck ${CRAM}
    if [ \$? -eq 0 ]
    then
        cp ${CRAM} ${name}_ok.cram
    fi
    """

}


/*

   Read extraction.  

*/

process extract_fastq {
    tag                    { name }
    executor               myExecutor
    cpus                   4
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 params.samtoolsMod
    publishDir             readPrepPath, mode: "copy"
    scratch                "/scratch"    
    
    input:
    set val(id), file(cram) from extract_ch	
    file genome from orig_genome_file
    file index from orig_genome_index_ch

    output:
    set val(id), file("*.R{1,2}.fastq.gz") optional true into fq_pe_ch
    
    script:    
    """
    samtools collate -@ ${task.cpus} -f --reference ${genome} -o tmp.bam ${cram}
    
    samtools fastq -@ ${task.cpus} \\
        -1 ${id}.R1.fastq.gz -2 ${id}.R2.fastq.gz tmp.bam
        
    rm tmp.bam
    """

}

/*

   bowtie2 alignment

*/

process bowtie2_aln {
    tag                    { name }
    executor               myExecutor
    cpus                   24
    memory                 "48 GB"
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 params.samtoolsMod,params.bowtieMod
    publishDir             alnPath, mode: "copy"// 
    scratch                "/scratch"
    
    input:
    set val(id), file(alnReads) from fq_pe_ch
    file genome from new_genome_file
    file idx from bowtie2Index
    
    output:
    set val(id), file("*.bam") into bam_ch
    
    script:    
    """
    bowtie2 -p ${task.cpus - 2} \
        -x ${new_genome_file.getParent()}/${new_genome_file.getBaseName()} \
        -1 ${alnReads[0]} -2 ${alnReads[1]} 2> ${id}.runlog \
        | samtools sort -T . -@ 2 -o ${id}.bowtie2.bam 
    """

}

/*

   BAM to CRAM

*/


process samtools_cram {
    tag                    { name }
    executor               myExecutor
    cpus                   8
    memory                 "48 GB"
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 params.samtoolsMod
    publishDir             alnPath, mode: "copy"// 
    scratch                "/scratch"
    
    input:
    set val(id), file(bam) from bam_ch
    file genome from new_genome_file
    
    output:
    file("*.cram")
    
    script:    
    """           
    samtools view -C -h -@ ${task.cpus} -o ${id}.bowtie2.cram \
          --reference ${genome} ${bam}
    """

}