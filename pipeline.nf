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
