nextflow.enable.dsl=2

params.R1 = ''
params.R2 = ''
params.attb_list = ''
params.attp_oligo = ''
params.sample_name = ''
params.method = 'cs2'

params.attb_flank_left = ''
params.attb_flank_right = ''
params.attp_flank_left = ''
params.attp_flank_right = ''


process ADAPTER_AND_POLY_G_TRIM {
    cache 'lenient'

    input:
        tuple path(R1), path(R2)
        val sample_name
    output:
        path "${sample_name}_trimmed.fastq.gz", emit: trimmed_fastq_reads
        val sample_name

    script:
        """
        fastp -m -c --include_unmerged --dont_eval_duplication --low_complexity_filter --overlap_len_require 10 -i ${R1} -I ${R2} --merged_out ${sample_name}_trimmed.fastq.gz -w 16 -g -j ${sample_name}_fastp.json
        """
}

process CREATE_AMPLICONS {
    cache 'lenient'    

    input:
        path attb_list
        path attp_oligo
        val method
        val sample_name
    output:
        path "${sample_name}_amplicons.txt", emit: amplicons

    script:

    """
    create_amplicon_files.py --attb_list ${attb_list} --attp ${attp_oligo} --output ${sample_name}_amplicons.txt --method ${method}
    """  
}

process CS2_POOLED {

    input:
        path merged_reads
        path amplicons
        val sample_name

    output:

    script:

    """
    CRISPRessoPooled -r1 ${merged_reads} -f ${amplicons} --min_reads_to_use_region 1  -p 8 -o ./ -n ${sample_name} --write_detailed_allele_table --bam_output --place_report_in_output_folder --suppress_report --suppress_plots --limit_open_files_for_demux
    """

}

process DIRECT_SEARCH {
    publishDir "${params.outdir}/${sample_name}/"
    input:
        path merged_reads
        val attb_left_flank
        val attb_right_flank
        val attp_left_flank
        val attp_right_flank    
        path amplicons
        val sample_name
    output:
        path("${sample_name}_recombination_data.csv")
        val sample_name
    script:
    """
    direct_search.py --fastq_file ${merged_reads} --attb_flank_left ${attb_left_flank} --attb_flank_right ${attb_right_flank} --attp_flank_left ${attp_left_flank} --attp_flank_right ${attp_right_flank} --amplicons_file ${amplicons} --sample_name ${sample_name}
    """
}


workflow {
    
    Channel.fromPath(params.samplesheet)
    .splitCsv(header: true, sep: ',')
    .map { row -> 
        row.sample_name  // Emit only the sample_name
    }
    .set { samplename_ch }

    Channel.fromPath(params.samplesheet)
    .splitCsv(header: true, sep: ',')
    .map { row -> 
        def r1 = "${launchDir}/${row.fastq_dir}/*R1*.fastq.gz"
        def r2 = "${launchDir}/${row.fastq_dir}/*R2*.fastq.gz"
        tuple(
            file(r1), 
            file(r2), 
        )
    }
    .set { reads_ch }

    Channel.fromPath(params.samplesheet)
    .splitCsv(header: true, sep: ',')
    .map { row -> 
        def oligos = "${launchDir}/${row.oligos}"
        file(oligos)
    }
    .set { oligos_ch }

    Channel.fromPath(params.samplesheet)
    .splitCsv(header: true, sep: ',')
    .map { row -> 
        row.group 
    }
    .set { group_ch }

    trimmed_reads = ADAPTER_AND_POLY_G_TRIM(reads_ch, samplename_ch)
    if (params.method == 'cs2'){
        amplicons = CREATE_AMPLICONS(oligos_ch,params.attp_oligo,params.method,samplename_ch)
        cs2_output = CS2_POOLED(trimmed_reads.trimmed_fastq_reads, amplicons.amplicons,samplename_ch)
    } else {
        amplicons = CREATE_AMPLICONS(oligos_ch,params.attp_oligo,params.method,samplename_ch)
        output = DIRECT_SEARCH(trimmed_reads.trimmed_fastq_reads, params.attb_flank_left, params.attb_flank_right, params.attp_flank_left, params.attp_flank_right, amplicons.amplicons,samplename_ch)
    }
}

workflow.onComplete {

    def msg = """\
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """
        .stripIndent()
}
