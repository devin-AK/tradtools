#!/usr/bin/env nextflow

/*
 * see https://github.com/CRG-CNAG/CalliNGS-NF/ for inspiration
 * https://github.com/czbiohub/nf-bowtie/blob/master/main.nf
 */
// https://carpentries-incubator.github.io/workflows-nextflow/05-processes-part1/index.html
// To Do: set cores to a variable (cpu.tasks)
// To Do: configure cpus & memory for individual processes / executors 
// To Do: pre-flight checks (e.g. can tradtoolsR be loaded)
// To Do: option to keep or discard intermediate files
// To Do: make sure base name of bt2 index is "index"
// To Do: include fastqc step
// To Do: make sure chromosome names from the two genome fasta files are unique

params.R1         = "$projectDir/example_data/small_R1.fastq.gz"
params.R2         = "$projectDir/example_data/small_R2.fastq.gz"
params.genome1    = "$projectDir/example_data/hg19_chr17.fa.gz"
params.genome2    = "$projectDir/example_data/lambda.fa.gz"
params.blacklist1 = "$projectDir/example_data/hg19-blacklist.v2.bed.gz"
params.blacklist2 = null
params.adapters   = "$projectDir/example_data/adapters2.fa"
params.results    = "$launchDir/results"

log.info """\

TRADtools v 0.1
================================================================================

|--------|
| INPUTS |
|--------|

  R1 fastq       : $params.R1
  R2 fastq       : $params.R2
  genome1 fasta  : $params.genome1
  genome2 fasta  : $params.genome2
  blacklist1 bed : $params.blacklist1
  blacklist2 bed : $params.blacklist2
  adapters fasta : $params.adapters
  results dir    : $params.results

================================================================================
"""

//WORKFLOW

workflow {
  //INPUTS
    //genomes
    Channel.fromPath( [params.genome1, params.genome2], checkIfExists: true)
      .collect()
      .set {genomes}
    Channel.fromPath( params.genome1 )
      .set {genome1}
    Channel.fromPath( params.genome2 )
      .set {genome2}
    //fastq reads
    Channel.fromPath( [params.R1, params.R2], checkIfExists: true)
      .collect()
      .set {reads}
    //adapters
    Channel.fromPath( params.adapters, checkIfExists: true)
      .set {adapters}
    //blacklist files
    if ( !params.blacklist1 ) {
      Channel.fromPath( "NULL" )
        .set {blacklist1}
    } else {
      Channel.fromPath( params.blacklist1 )
        .set {blacklist1}
    }
    if ( !params.blacklist2 ) {
      Channel.fromPath( "NULL" )
        .set {blacklist2}
    } else {
      Channel.fromPath( params.blacklist2 )
        .set {blacklist2}
    }



  //PIPELINE

  PREPARE_GENOMES(genomes)
      PREPARE_GENOMES.out.genomeIndex
        .collect()
        .set {index}
      PREPARE_GENOMES.out.chrNameFiles
        .collect()
        .set {chrnames}

  DEMULTIPLEX(reads, adapters) // Demultiplex the raw Fastq files using adapters
      DEMULTIPLEX.out.fq1Files // The following operators transform the demux R1 fastq files into key-value pairs to use with .combine()
        .flatten()
        .map(x -> tuple(x.simpleName,x))
        .set {demux_R1}
      DEMULTIPLEX.out.fq2Files
        .flatten()
        .map(x -> tuple(x.simpleName,x))
        .set {demux_R2}
      demux_R1.combine(demux_R2, by: 0) // Use combine to generate a channel "demux" with correctly linked demux R1 and R2 fastq files
        .set {demux}

  ALIGN(index, demux)
      ALIGN.out.bamFiles
        .flatten()
        .map(x -> tuple(x.simpleName,x))
        .set {bams} // bams is channel with bam files

  SPLIT_BY_SPECIES(bams, chrnames)
      SPLIT_BY_SPECIES.out.spec1BamFiles
        .flatten()
        .map(x -> tuple(x.simpleName,x))
        .combine(blacklist1)
        .combine(genome1)
        .set {bams1} // bams1 is channel with bam files from species 1
      SPLIT_BY_SPECIES.out.spec2BamFiles
        .flatten()
        .map(x -> tuple(x.simpleName,x))
        .combine(blacklist2)
        .combine(genome2)
        .set {bams2} // bams2 is channel with bam files from species 2
      bams1.concat(bams2)
        .set {bamsfinal} // the final list of bam files with associated blacklist files

  SIGNAL_TRACK(bamsfinal)

}



// PROCESSES

process PREPARE_GENOMES {
  publishDir "${params.results}", mode: "copy"
  cpus 8 

  input:
    tuple path(fasta1), path(fasta2)
  output:
    path "genome/STARgenome", emit: genomeIndex
    path "genome/fasta/*.fa", emit: genomeFastaFiles
    path "genome/*_chr.bed", emit: chrNameFiles

  shell:
  """
    mkdir -p genome/bowtie2_index
    mkdir -p genome/fasta
    gunzip -c ${fasta1} > genome/fasta/spec1.fa
    gunzip -c ${fasta2} > genome/fasta/spec2.fa
    cat genome/fasta/spec1.fa genome/fasta/spec2.fa > genome/fasta/combined_genomes.fa
    samtools faidx genome/fasta/spec1.fa
    samtools faidx genome/fasta/spec2.fa
    awk '{printf("%s\\t0\\t%s\\n",\$1,\$2);}' genome/fasta/spec1.fa.fai > genome/spec1_chr.bed
    awk '{printf("%s\\t0\\t%s\\n",\$1,\$2);}' genome/fasta/spec2.fa.fai > genome/spec2_chr.bed 
    STAR  --runMode genomeGenerate \
          --runThreadN ${task.cpus} \
          --genomeDir genome/STARgenome \
          --genomeFastaFiles genome/fasta/combined_genomes.fa
  """
}


process DEMULTIPLEX {
  publishDir "${params.results}", mode: "copy"
  cpus 16 

  input:
    tuple path(fastq_R1), path(fastq_R2)
    path adapters
  output:
    path "fastq/*R1.fastq", emit: fq1Files
    path "fastq/*R2.fastq", emit: fq2Files
    path "log/demux/*"

  """
    mkdir -p fastq
    mkdir -p log/demux
    cutadapt --json=log/demux/fragment_construction.cutadapt.json \
      --cores=${task.cpus} -m 50:55 -a CTGTCTCTTATACACATCTAGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
      -A "CACTGCNNNNNNNNAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT;max_error_rate=0.15;min_overlap=11" \
      -o R1_trim.fastq -p R2_trim.fastq ${fastq_R1} ${fastq_R2}
    cutadapt --cores=${task.cpus} --no-indels --discard-untrimmed -g file:${adapters} \
      -G AGATGTGTATAAGAGACAG \
      -o fastq/{name}.R1.fastq -p fastq/{name}.R2.fastq \
      R1_trim.fastq R2_trim.fastq
  """
}


process ALIGN {
  publishDir "${params.results}", mode: "copy"
  cpus 8 

  input:
    path genome_DIR
    tuple val(library_ID), path(demux_R1), path(demux_R2)
  output:
    path "bam/initial/*_possort_markdup.bam", emit: bamFiles
    path "log/align/*"
    path "log/duplicates/*"

  """
    mkdir -p tmp_sam
    mkdir -p bam/initial
    mkdir -p log/align
    mkdir -p log/duplicates

    STAR --runThreadN ${task.cpus} \
          --genomeDir ${genome_DIR} \
          --readFilesIn ${demux_R1} ${demux_R2} \
          --readFilesCommand awk "'NR%4==1{print\\\$1\\\":\\\"substr(\\\$2,7)}NR%4!=1{print\\\$0}'" \
          --alignIntronMax 1 \
          --alignEndsType EndToEnd \
          --outFilterScoreMinOverLread 0.3 \
          --outFilterMatchNminOverLread 0.3 \
          --outFilterMatchNmin 36 \
          --outSAMtype SAM \
          --outFileNamePrefix tmp_sam/${library_ID}_
    cp tmp_sam/*.out log/align

    samtools fixmate --threads ${task.cpus} -m -u tmp_sam/${library_ID}_Aligned.out.sam - | samtools sort --threads ${task.cpus} -u - | \
    samtools markdup --threads ${task.cpus} -f tmp_sam/${library_ID}_markdup_stats.txt -s -c -d 2500 --barcode-rgx \"^[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:[0-9]+:[0-9]+:[0-9]+:([A|C|G|T|N]{6,8})\$\" - bam/initial/${library_ID}_possort_markdup.bam

    Rscript -e "tradtoolsR::parse_samtools_markdup_stats(input_file='tmp_sam/${library_ID}_markdup_stats.txt',output_file='log/duplicates/${library_ID}_trad_markdups.csv')"

  """
}


process SPLIT_BY_SPECIES {
  publishDir "${params.results}", mode: "copy"

  input:
    tuple val(library_ID), path(bam_file)
    tuple path(chrname1), path(chrname2)
  output:
    path "bam/species1/*.bam", emit: spec1BamFiles
    path "bam/species2/*.bam", emit: spec2BamFiles
  
  """
    mkdir -p bam/species1
    mkdir -p bam/species2
    samtools view -L ${chrname1} -o bam/species1/${library_ID}_spec1.bam ${bam_file}
    samtools view -L ${chrname2} -o bam/species2/${library_ID}_spec2.bam ${bam_file}
  """
}


process SIGNAL_TRACK {
  debug true
  publishDir "${params.results}", mode: "copy"

  input:
    tuple val(library_ID), path(input_bam), path(input_blacklist), path(reference_fasta)
  output:
    path "bed/*.bed.gz", emit: bedFiles
    path "bw/*.bw", emit: bwFiles
    path "log/signal_track/*"
  
  """
  #!/usr/bin/env Rscript
  system('mkdir -p bed')
  system('mkdir -p bw')
  system('mkdir -p log/signal_track')
  library_ID <- '${library_ID}'
  input_bam <- '${input_bam}'
  input_blacklist <- '${input_blacklist}'
  reference_fasta <- '${reference_fasta}'
  st <- tradtoolsR::build_signal_track(input_bam = input_bam,
    input_blacklist = input_blacklist,
    reference_fasta = reference_fasta,
    output_bed = paste0('bed/',library_ID,'.bed.gz'),
    output_bw = paste0('bw/',library_ID,'.bw'),
    export_stats = TRUE,
    ranges_operation = 'default',
    mapQ = 10,
    standard_chromosomes = FALSE,
    overwrite = TRUE,
    verbose = TRUE)
  ID <- tools:::file_path_sans_ext(basename(input_bam))
  dinuc_log <- paste0(ID,'_trad_dinucs.csv')
  count_log <- paste0(ID,'_trad_readcounts.csv')
  invisible(file.copy(dinuc_log,'log/signal_track/'))
  invisible(file.copy(count_log,'log/signal_track/'))
  """
}






