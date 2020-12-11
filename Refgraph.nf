#!/usr/bin/env nextflow

/*parameters that are specified at the command line or via config file*/
params.genome                = false          /*genome fasta file, must specify complete path. Required parameters*/
params.samplePath            = false          /*input folder, must specify complete path. Required parameters*/
params.outputDir             = false          /*output folder, must specify complete path. Required parameters*/
params.singleEnd             = false          /*options: true|false. true = the input type is single end reads; false = the input type is paired reads. Default is false*/
params.assembler             = 'megahit'      /*options: megahit|masurca. Default is megahit*/
params.skipKraken2           = true           /*options: true|false. Default is true which means that kraken2 will be skipped*/

/* parameters for readprep = qctrimming and adapter removal */
params.skipTrim              = false           /*qc-trimming of reads. options: true|false. Default is false*/     
params.min_read_length       = '20'            /*minimum length of read to be kept after trimming for downstream analysis. Default is 20*/
params.min_base_quality      = '20'            /*minimum base quality. Default is 20*/
params.guess_adapter         = true            /*auto-detect adapter from input file. options: true|false. Default is true*/
params.forward_adapter       = false           /*adapter sequence to be clipped off (forward). */
params.reverse_adapter       = false           /*adapter sequence to be clipped off (reverse). Used for paired reads only*.*/


/*output folder paths*/
readPrepPath                 = "${params.outputDir}/read_prep"
trimPath                     = "${params.outputDir}/trimmed"
megahitPath                  = "${params.outputDir}/megahit"
masurcaPath                  = "${params.outputDir}/masurca"
kraken2Path                  = "${params.outputDir}/kraken2"
multiqcPath                  = "${params.outputDir}/multiqc"
metricsPath                  = "${params.outputDir}/assembly_metrics"

/*cluster parameters */
myExecutor                   = 'slurm'
params.myQueue               = 'normal'
defaultCPU                   = '1'
defaultMemory                = '20'
assemblerCPU                 = '12'
assemblerMemory              = '100'
// params.clusterAcct           = " -A h3bionet "

/*software stack*/
params.perlMod               = 'Perl/5.24.1-IGB-gcc-4.9.4'
params.fastpMod              = 'fastp/0.20.0-IGB-gcc-4.9.4'
params.samtoolsMod           = 'SAMtools/1.10-IGB-gcc-8.2.0'
params.megahitMod            = 'MEGAHIT/1.2.9-IGB-gcc-8.2.0'
params.assemblathon          = "/home/groups/hpcbio/apps/FAlite/assemblathon_stats.pl"
params.multiqcMod            = "MultiQC/1.7-IGB-gcc-4.9.4-Python-3.6.1"

/*Prepare input*/
genome_file                  = file(params.genome)
genomeStore                  = genome_file.getParent()
if( !genome_file.exists() ) exit 1, "Missing reference genome file: ${genome_file}"
CRAM_Ch1 = Channel.fromFilePairs("${params.samplePath}", size: 1)

/*

  prepare_genome 
  This process is executed only once

*/


process prepare_genome{
    tag                    { genome }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   defaultCPU
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 params.samtoolsMod
    storeDir               genomeStore
    validExitStatus        0
    
    input:
    file genome from genome_file

    output:
    file "*.fai" into genome_index_ch
    
    script:
    """
    samtools faidx ${genome}
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

   Read extraction.  This step is tricky.
   
   For single end reads, we only need to extract unmapped as there isn't a mapped mate.  
   However this also means that the additional workflows downstream cannot be run as they 
   paired end read data (location placement requires having mapped mates)
   
   samtools view -f 4 > se_unmapped.bam
   For paired end reads, we need to extract three classes of reads, keeping the mates
   
   1. Both unmapped
   2. R1 mapped, R2 unmapped
   3. R2 mapped, R1 unmapped
   
   This can be accomplished with samtools though the logic is a little tricky with the bit
   flags.  Best way would be to extract pairs where any of the reads are unmapped, then 
   pull out subsets using flags. Otherwise we're running through the 
   
   # both unmapped; -f means unmapped, 
   # bit flag 12  = both reads unmapped (bit flag 4 & 8)
   # bit flag 2304 = not primary alignment (this removes these), not supplemental
   samtools view -hb -f 12 -F 2304 alignments.bam > both_unmapped.bam
    
   # R1 mapped, R2 not
   # bit flag 4   = R1 unmapped
   # bit flag 2312 = Mate unmapped and not primary alignment (removes these), not supplemental
   #                Note this is to make sure we're not keeping reads *also* in the first  
   #                set
   samtools view -hb -f 4 -F 2312 alignments.bam  > R1_unmapped.bam
   
   # R2 mapped, R1 not
   # bit flag 8    = R2 (mate) unmapped
   # bit flag 2308 = R1 (read) unmapped and not primary alignment (removes these), not supplemental
   #                Note this is to make sure we're not keeping reads *also* in the first  
   #                set
   samtools view -hb -f 8 -F 2308 alignments.bam  > R2_unmapped.bam

*/

/*

  extract_unmap 
  The input files are in CRAM format which have been previously aligned to genome prepared in previous step

*/

process extract_unmapped {
    tag                    { name }
    executor               myExecutor
    cpus                   4
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 params.samtoolsMod
    publishDir             readPrepPath, mode: "copy"
    validExitStatus        0,1
    errorStrategy          'finish'
    scratch                '/scratch'
    stageOutMode           'copy'
    
    input:
    set val(id), file(cram) from extract_ch	
    file genome from genome_file
    file index from genome_index_ch

    output:
    set val(id), file("${id}.both-unmapped.R{1,2}.fastq") optional true into fq_pe_ch
    set val(id), file("${id}.orphans.unmapped.fastq") optional true into fq_se_ch
    file "${id}.improper.bam" optional true // all improper pairs; we save this for now, might be useful later
    set val(id), file("${id}.unmapped.bam") optional true into unmapped_bam_ch  // all unaligned + mates
    
    script:
    if(params.singleEnd) {
    
    """
    # TODO: UNTESTED!!!!
    # we only need to extract reads that are unmapped, no worries about pairing
    samtools view -@ ${defaultCPU} -hbt ${index} -f 4 -o ${id}.unmapped.bam ${cram} 
    
    # convert to FASTQ
    samtools fastq -@ ${defaultCPU} ${id}.SE.unmapped.bam > ${id}.orphans.unmapped.fastq
    """
    
    } else {
    
    """
    # two stages; grab any non-properly paired reads (includes unmapped)
    samtools view -@ ${task.cpus} -hbt ${index} -G 2 -o ${id}.improper.bam ${cram}
    
    # now capture the three classes; this is much faster than parsing the full BAM each time
    
    # both unmapped
    samtools view -@ ${task.cpus} -hbt ${index} \\
        -f 12 -F 2304 -o ${id}.both-unmapped.bam ${id}.improper.bam

    samtools fastq -@ ${task.cpus} ${id}.both-unmapped.bam \\
        -1 ${id}.both-unmapped.R1.fastq -2 ${id}.both-unmapped.R2.fastq
        
    # R1 only unmapped
    samtools view -@ ${task.cpus} -hbt ${index} \\
        -f 4 -F 2312 -o ${id}.R1-unmapped.bam ${id}.improper.bam

    samtools fastq -@ ${task.cpus} ${id}.R1-unmapped.bam \\
        -1 ${id}.R1-unmapped.R1.fastq -2 ${id}.R1-unmapped.R2.fastq

    # R2 only unmapped
    samtools view -@ ${task.cpus} -hbt ${index} \\
        -f 8 -F 2308 -o ${id}.R2-unmapped.bam ${id}.improper.bam
        
    samtools fastq -@ ${task.cpus} ${id}.R1-unmapped.bam \\
        -1 ${id}.R2-unmapped.R1.fastq -2 ${id}.R2-unmapped.R2.fastq
    
    # combine unmapped BAM files for later analysis
    samtools merge -@ ${task.cpus} ${id}.unmapped.bam ${id}.both-unmapped.bam ${id}.R1-unmapped.bam ${id}.R2-unmapped.bam
    
    # combined SE unmapped FASTQ into one file for assembly
    cat ${id}.R1-unmapped.R1.fastq ${id}.R2-unmapped.R2.fastq >> ${id}.orphans.unmapped.fastq
    
    # we could combine the mapped reads here, but we'll use the original BAM instead
    ## cat ${id}.R1-unmapped.R2.fastq ${id}.R2-unmapped.R1.fastq >> ${id}.orphans.mapped.fastq
    
    """
    
    }

}


/*

  trimming

*/

process trimming {
    tag                    { name }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   2
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    publishDir             trimPath, mode: "copy"
    module                 params.fastpMod
//     validExitStatus        0,1
//     errorStrategy          'finish'
//     scratch                '/scratch'
//     stageOutMode           'copy'

    input:
    set val(name), file(reads) from fq_pe_ch
    
    output:
    set val(name), file('*.PE.R{1,2}.trimmed.fq'), file('*.unpR{1,2}.trimmed.fq') optional true into trim_pe_ch
    set val(name), file('*.json') optional true into multiqc_pe_ch
    file '*'
	    
    script:
    trimOptions      = params.skipTrim ? ' ' :  ' -l 20 -q 20  --cut_right --cut_right_window_size 3 --cut_right_mean_quality 20 '
    adapterOptionsSE = params.guess_adapter ? ' ' : " --adapter_sequence=${params.forward_adapter} "
    adapterOptionsPE = params.guess_adapter ? ' --detect_adapter_for_pe ' : " --adapter_sequence=${params.forward_adapter}  --adapter_sequence_r2=${params.reverse_adapter} "

    if(params.singleEnd){
    """
    fastp --in1 ${reads[0]} --out1 "${name}.SE.R1.trimmed.fq" \\
        ${adapterOptionsSE} ${trimOptions} --thread ${task.cpus} \\
        -w ${task.cpus} --html "${name}"_SE_fastp.html \\
        --json "${name}"_SE_fastp.json
    """
    } else {
    """
    fastp --in1 ${reads[0]} --in2 ${reads[1]} \\
        --out1 "${name}.PE.R1.trimmed.fq"  --out2 "${name}.PE.R2.trimmed.fq" \\
        --unpaired1 "${name}.unpR1.trimmed.fq" --unpaired2 "${name}.unpR2.trimmed.fq" \\
        ${adapterOptionsPE}  ${trimOptions} --thread ${task.cpus} \\
        -w ${task.cpus}  --html "${name}"_PE_fastp.html \\
        --json "${name}"_PE_fastp.json
    """
    }
}


process trimming_orphans {
    tag                    { name }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   2
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    publishDir             trimPath, mode: "copy"
    module                 params.fastpMod

    input:
    set val(name), file(reads) from fq_se_ch
    
    output:
    set val(name), file('*.orphans.trimmed.fq') optional true into trim_orphan_ch
    set val(name), file('*.json') optional true into multiqc_orphan_ch
    file '*'
	    
    script:
    trimOptions      = params.skipTrim ? ' ' :  ' -l 20 -q 20  --cut_right --cut_right_window_size 3 --cut_right_mean_quality 20 '
    adapterOptions   = params.guess_adapter ? ' ' : ' --adapter_fasta="adapters.fa" '
    """
    # we generate an adapter sequence file on the fly here (should move this into a stored process)
    
    cat << ADAPTERS > adapters.fa
    >forward_adapter
    ${params.forward_adapter}
    >reverse_adapter
    ${params.reverse_adapter}
    ADAPTERS
    
    fastp --in1 ${reads[0]} --out1 "${name}.orphans.trimmed.fq" \\
        ${adapterOptions} ${trimOptions} --thread ${task.cpus} \\
        -w ${task.cpus} \\
        --html "${name}"_orphans_fastp.html \\
        --json "${name}"_orphans_fastp.json
    """
}



/*
  *Megahit for different input data types:
  *megahit -1 pe_1.fq -2 pe_2.fq -o out  # 1 paired-end library
  *megahit --12 interleaved.fq -o out # one paired & interleaved paired-end library
  *megahit -1 a1.fq,b1.fq,c1.fq -2 a2.fq,b2.fq,c2.fq -r se1.fq,se2.fq -o out # 3 paired-end libraries + 2 SE libraries
*/


// we will want to split this by SE CRAM and PE CRAM at some point; 
// logic is different enough to be an issue

// TODO: at the moment SE CRAM assembly is **not** supported, we will need an if/else
// to handle this one

process megahit_assemble {
    tag                    { name }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   assemblerCPU
    queue                  params.myQueue
    memory                 "$assemblerMemory GB"
    module                 params.megahitMod 
    publishDir             megahitPath , mode:'copy'
    validExitStatus        0,1
    errorStrategy          'finish'
    scratch                '/scratch'
    stageOutMode           'copy'
    
    // TODO: We need to sanity check that the channels are matched by the initial value (name)
    // See the join command: https://www.nextflow.io/docs/latest/operator.html?highlight=map#join
    input:
    set val(name), file(pefastqs), file(sefastqs) from trim_pe_ch
    set val(name2), file(orphans) from trim_orphan_ch

    output:
    set val(name), file('*.stats') optional true into metrics_ch
    file '*'

    script:
    if(params.singleEnd){
    """
    megahit -1 ${fastqs[0]} -o ${name}.megahit_results
    
    perl $params.assemblathon ${name}.megahit_results/final.contigs.fa > ${name}.megahit_results/final.contigs.fa.stats
    """
    } else {
    """
    megahit -1 ${pefastqs[0]} -2 ${pefastqs[1]} \\
        -r ${sefastqs[0]},${sefastqs[1]},${orphans} \\
        -o ${name}.megahit_results

    # megahit -1 ${pefastqs[0]} -2 ${pefastqs[1]}  -o ${name}.megahit_results

    perl $params.assemblathon ${name}.megahit_results/final.contigs.fa > ${name}.final.contigs.fa.stats

    """
    }
}

/*

  masurca
  To be added later

*/

/*

  kraken2
  To be added later

*/

/*

  QC all samples

*/
process MultiQC_readPrep {
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   2
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 params.multiqcMod
    publishDir             multiqcPath, mode: 'copy', overwrite: true
 
    input:
    file('*') from multiqc_pe_ch.collect()
    file('*') from multiqc_orphan_ch.collect()

    output:
    file "*"

    """
    multiqc .
    """
} 

process Assembly_metrics {
    executor               myExecutor
    cpus                   2
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 params.perlMod
    publishDir             metricsPath, mode: 'copy', overwrite: true
 
    input:
    file('*') from metrics_ch.collect()

    output:
    file "*"

    """
    cat *.stats > assembly_metrics.txt
    """
} 

