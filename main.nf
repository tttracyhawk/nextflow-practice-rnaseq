params.reads = "${projectDir}/data/ggal/*_{1,2}.fq"
params.transcriptome_file = "${projectDir}/data/ggal/transcriptome.fa"
params.multiqc = "${projectDir}/multiqc"
params.outdir = "results"

log.info """
    RNASEQ-NF PIPELINE
    ==================
    transcriptome:  ${params.transcriptome_file}
    reads:          ${params.reads}
    outdir:         ${params.outdir}
""".stripIndent(true)


//define the `INDEX` process that creates a binary index
//given the transcriptome file

process INDEX {
    label 'process_long'

    input:
    path transcriptome

    output:
    path 'salmon_index'

    script:
    """
    salmon index --threads $task.cpus -t $transcriptome -i salmon_index
    """
}

process QUANTIFICATION {
    label 'process_long'
    tag "$sample_id"
    publishDir "$params.outdir/salmon_quant", mode: 'copy'

    input:
    path salmon_index
    tuple val(sample_id), path(reads)

    output:
    path "$sample_id/quant.sf", emit: quant
    path "$sample_id/logs_$sample_id", emit: logs

    script:
    """
    salmon quant \
        --threads $task.cpus \
        --libType=U \
        -i $salmon_index \
        -1 ${reads[0]} \
        -2 ${reads[1]} \
        -o $sample_id
    
    mv $sample_id/logs $sample_id/logs_$sample_id
    """
}

process FASTQC {
    label 'process_short'
    tag "FASTQC on $sample_id"
    publishDir "$params.outdir/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*.html", emit: html
    path "*.zip", emit: zip

    script:
    """
    fastqc $reads
    """
}

process MULTIQC {
    label 'process_short'
    publishDir "$params.outdir/multiqc", mode: 'copy'

    input:
    path '*'

    output:
    path 'multiqc_report.html'

    script:
    """
    multiqc .
    """
}

workflow {
    // Set up input files channels
    read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)

    index_ch = INDEX(params.transcriptome_file)
    quant_ch = QUANTIFICATION(index_ch, read_pairs_ch)

    fastqc_ch = FASTQC(read_pairs_ch)
    multiqc_in_ch = fastqc_ch.zip.mix(quant_ch.logs).collect()
    MULTIQC(multiqc_in_ch)
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc/multiqc_report.html\n" : "Oops .. something went wrong" )
}