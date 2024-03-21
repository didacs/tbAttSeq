nextflow.enable.dsl=2

params.samplesheet = '' // Assuming you have this parameter to specify the path to the samplesheet
params.attp_oligo = ''
params.outdir = 'results' // Assuming an output directory parameter

process ADAPTER_AND_POLY_G_TRIM {
    cache 'lenient'

    input:
        tuple val(sample_name), path(R1), path(R2), path(oligos), val(group)
    output:
        tuple val(sample_name), path("${sample_name}_trimmed.fastq.gz"), emit: fastq
        path("${sample_name}_fastp.json"), emit: fastp_stats

    script:
        """
        fastp -m -c --include_unmerged --dont_eval_duplication --low_complexity_filter --overlap_len_require 10 -i ${R1} -I ${R2} --merged_out ${sample_name}_trimmed.fastq.gz -w 16 -g -j ${sample_name}_fastp.json
        """
}

process CREATE_AMPLICONS {
    cache 'lenient'

    input:
        tuple val(sample_name), path(R1), path(R2), path(oligos), val(group)
    output:
        tuple val(sample_name), path("${sample_name}_amplicons.txt")

    script:
        """
        create_amplicon_files.py --attb_list ${oligos} --attp ${params.attp_oligo} --output ${sample_name}_amplicons.txt
        """
}

process CREATE_OLIGOS_REF {
    input:
        path oligos_path
        val attb_flank_left
        val attb_flank_right
    output:
        path "oligos.fa", emit: index
        path "oligos.fa.{amb,ann,bwt,pac,sa}", emit: ref_files
    script:
        """
        create_oligo_fasta.py \
            --attb_oligos ${oligos_path} \
            --attb_flank_left ${attb_flank_left} \
            --attb_flank_right ${attb_flank_right} \
            --output_fasta oligos.fa && \
        bwa index oligos.fa
        """
}

process ALIGN_READS_TO_OLIGOS {
    input:
        tuple val(sample_name), path(trimmed_reads)
        path oligos_index
        path oligos_ref_files
    output:
        tuple val(sample_name), path("${sample_name}.align_reads_to_oligos.bam")
    script:
    """
    bwa mem -p ${oligos_index} ${trimmed_reads} |\
        samtools view -b |\
        samtools sort --write-index -o ${sample_name}.align_reads_to_oligos.bam
    """
}

process DETECT_DSB {
    publishDir "${params.outdir}/${sample_name}", mode: "copy"
    input:
        tuple val(sample_name), path(aligned_bam)
    output:
        path "${sample_name}.dsb_counts.csv"
        path "${sample_name}.read_boundaries.csv"
    script:
    """
    dsb_quantification.py \
        --sample_name ${sample_name} \
        --bam ${aligned_bam} \
        --dinucleotide_position 44
    """
}


process DIRECT_SEARCH {
    publishDir "${params.outdir}/${sample_name}", mode: 'copy'
    input:
        tuple val(sample_name), path(merged_reads), path(amplicons)
        val(attb_left_flank)
        val(attb_right_flank)
        val(attp_left_flank)
        val(attp_right_flank)
    output:
        val(sample_name), emit: sample_name
        path("${sample_name}_recombination_data.csv"), emit: recombination_data

    script:
        """
        direct_search.py --fastq_file ${merged_reads} --attb_flank_left ${attb_left_flank} --attb_flank_right ${attb_right_flank} --attp_flank_left ${attp_left_flank} --attp_flank_right ${attp_right_flank} --amplicons_file ${amplicons} --sample_name ${sample_name}
        """
}

process COLLATE_RESULTS {
    publishDir "${params.outdir}", mode: 'copy'
    input:
        val(sample_names)
        path(recombination_files)
    output:
        path "results.csv"
    script:
        """
        collate_results.py -f ${recombination_files} -n ${sample_names} -o results.csv
        """
}

process MULTIQC {
    cache 'lenient'
    publishDir "${params.outdir}", mode: 'copy'
    input:
        path fastp_jsons
    output:
        file "multiqc_report.html"
        file "multiqc_data"
    script:
    """
        multiqc .
    """
}

workflow {
    // Create a channel for the input sample sheet
    Channel.fromPath(params.samplesheet)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def r1 = "${launchDir}/${row.fastq_dir}/*R1*.fastq.gz"
            def r2 = "${launchDir}/${row.fastq_dir}/*R2*.fastq.gz"
            def oligos = "${launchDir}/${row.oligos}"

            tuple(row.sample_name, file(r1), file(r2), file(oligos), row.group)
        }
        .set { combined_ch }

    trimmed_reads = ADAPTER_AND_POLY_G_TRIM(combined_ch)

    amplicons = CREATE_AMPLICONS(combined_ch)

    // Create reference from oligos
    create_oligo_ref_ch = Channel.fromPath(params.oligos_path)
    oligos_ref = CREATE_OLIGOS_REF(create_oligo_ref_ch, params.oligo_attb_flank_left, params.oligo_attb_flank_right)

    // align reads to oligos (substrates)
    aligned_bam = ALIGN_READS_TO_OLIGOS(trimmed_reads.fastq, oligos_ref.index.first(), oligos_ref.ref_files.first())

    // DSB detection
    DETECT_DSB(aligned_bam)

    trimmed_reads.fastq
        .combine(amplicons, by: 0)
        .set{trimmed_reads_and_amplicons}

    recombination_files = DIRECT_SEARCH(trimmed_reads_and_amplicons, params.attb_flank_left, params.attb_flank_right, params.attp_flank_left, params.attp_flank_right)

    recombination_files.sample_name.collect().map { it.join(' ') }.view()

    COLLATE_RESULTS(recombination_files.sample_name.collect().map { it.join(' ') },recombination_files.recombination_data.collect())

    MULTIQC(trimmed_reads.fastp_stats
            .flatten()  // Flatten the list
            .filter { it.toString().endsWith('.json') }  // Filter out only the paths ending with .json
            .collect()
            .ifEmpty([]))
}

workflow.onComplete {
    println "Pipeline execution summary:"
    println "Completed at: ${workflow.complete}"
    println "Duration    : ${workflow.duration}"
    println "Success     : ${workflow.success}"
    println "workDir     : ${workflow.workDir}"
    println "exit status : ${workflow.exitStatus}"
}
