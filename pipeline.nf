nextflow.enable.dsl=2

params.samplesheet = '' // Assuming you have this parameter to specify the path to the samplesheet
params.attp_oligo = ''
params.method = 'cs2'
params.outdir = 'results' // Assuming an output directory parameter

process TRIM_AND_CREATE_AMPLICONS {
    cache 'lenient'

    input:
        tuple val(sample_name), path(R1), path(R2), path(oligos), val(group)
        val(method)
    output:
        tuple val(sample_name), path("${sample_name}_trimmed.fastq.gz"), path("${sample_name}_amplicons.txt")

    script:
        """
        fastp -m -c --include_unmerged --dont_eval_duplication --low_complexity_filter \
            --overlap_len_require 10 -i ${R1} -I ${R2} \
            --merged_out ${sample_name}_trimmed.fastq.gz -w 16 -g -j ${sample_name}_fastp.json
        create_amplicon_files.py --attb_list ${oligos} --attp ${params.attp_oligo} --output ${sample_name}_amplicons.txt --method ${method}
        """
}

// process CS2_POOLED {
//     input:
//         tuple val(sample_name), path(merged_reads)
//         path(amplicons)

//     output:
//         // Assuming there's an output you want to emit, add here

//     script:
//         """
//         CRISPRessoPooled -r1 ${merged_reads} -f ${amplicons} --min_reads_to_use_region 1  -p 8 -o ./ -n ${sample_name} --write_detailed_allele_table --bam_output --place_report_in_output_folder --suppress_report --suppress_plots --limit_open_files_for_demux
//         """
// }

process DIRECT_SEARCH {
    publishDir "${params.outdir}/${sample_name}", mode: 'copy'
    input:
        tuple val(sample_name), path(merged_reads), path(amplicons)

        val(attb_left_flank)
        val(attb_right_flank)
        val(attp_left_flank)
        val(attp_right_flank)
    output:
        path("${sample_name}_recombination_data.csv")

    script:
        """
        direct_search.py --fastq_file ${merged_reads} --attb_flank_left ${attb_left_flank} --attb_flank_right ${attb_right_flank} --attp_flank_left ${attp_left_flank} --attp_flank_right ${attp_right_flank} --amplicons_file ${amplicons} --sample_name ${sample_name}
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

    trimmed_reads_and_amplicons = TRIM_AND_CREATE_AMPLICONS(combined_ch, params.method)

    DIRECT_SEARCH(trimmed_reads_and_amplicons, params.attb_flank_left, params.attb_flank_right, params.attp_flank_left, params.attp_flank_right)

    // if (params.method == 'cs2') {
    //     CS2_POOLED(amplicons.amplicons, amplicons)
    // } else {
    //     DIRECT_SEARCH(trimmed_reads, params.attb_flank_left, params.attb_flank_right, params.attp_flank_left, params.attp_flank_right, amplicons)
    // }
}

workflow.onComplete {
    println "Pipeline execution summary:"
    println "Completed at: ${workflow.complete}"
    println "Duration    : ${workflow.duration}"
    println "Success     : ${workflow.success}"
    println "workDir     : ${workflow.workDir}"
    println "exit status : ${workflow.exitStatus}"
}
