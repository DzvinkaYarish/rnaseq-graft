process bam_to_junc {
    container = 'quay.io/eqtlcatalogue/leafcutter:v22.03.p4'
    publishDir "${params.outdir}/leafcutter/juncs", mode: 'copy', enabled: params.saveIndividualQuants

    input:
    tuple file(bam), file(bam_index)

    output:
    path "${bam.baseName}.junc", emit: junc
    path "${bam.baseName}_hg38.junc", emit: hg38_junc
    path "${bam.baseName}_gg6.junc", emit: gg6_junc
    
    script:
    // If confused about strands check this: https://rnabio.org/module-09-appendix/0009/12/01/StrandSettings/
    def leafcutter_strand = 0
    if (params.forward_stranded && !params.unstranded) {
        leafcutter_strand = 1
    } else if (params.reverse_stranded && !params.unstranded){
        leafcutter_strand = 2
    }
    """
    regtools junctions extract -s $leafcutter_strand -a 8 -m ${params.leafcutter_min_intron_length} -M ${params.leafcutter_max_intron_length} $bam -o ${bam.baseName}.junc
    #sed -i 's/gg6_/chr/' ${bam.baseName}.junc

    grep '^hg38' ${bam.baseName}.junc > ${bam.baseName}_hg38.junc
    grep '^gg6' ${bam.baseName}.junc > ${bam.baseName}_gg6.junc

    sed -i 's/gg6_/chr/' ${bam.baseName}_gg6.junc
    sed -i 's/hg38_/chr/' ${bam.baseName}_hg38.junc

    """
}

process cluster_introns {
    container = 'quay.io/eqtlcatalogue/leafcutter:v22.03.p4'
    tag "${junc_files.baseName}"
    publishDir "${params.outdir}/leafcutter", mode: 'copy'

    input:
    path junc_files
    path junc_files_hg38
    path junc_files_gg6

    output:
    path "leafcutter*.gz", emit: perind_counts
    path "*_refined", emit: refined

    script:
    """
    echo $junc_files
    leafcutter_cluster_regtools.py -j $junc_files -m ${params.leafcutter_min_split_reads} -o leafcutter -l ${params.leafcutter_max_intron_length}
    zcat leafcutter_perind_numers.counts.gz | sed '1s/^/phenotype_id /' | sed 's/.sorted//g' | sed -e 's/ /\t/g' | gzip -c > leafcutter_perind_numers.counts.formatted.gz

    echo $junc_files_hg38
    leafcutter_cluster_regtools.py -j $junc_files_hg38 -m ${params.leafcutter_min_split_reads} -o leafcutter_hg38 -l ${params.leafcutter_max_intron_length}
    zcat leafcutter_hg38_perind_numers.counts.gz | sed '1s/^/phenotype_id /' | sed 's/.sorted//g' | sed -e 's/ /\t/g' | gzip -c > leafcutter_hg38_perind_numers.counts.formatted.gz

    echo $junc_files_gg6
    leafcutter_cluster_regtools.py -j $junc_files_gg6 -m ${params.leafcutter_min_split_reads} -o leafcutter_gg6 -l ${params.leafcutter_max_intron_length}
    zcat leafcutter_gg6_perind_numers.counts.gz | sed '1s/^/phenotype_id /' | sed 's/.sorted//g' | sed -e 's/ /\t/g' | gzip -c > leafcutter_gg6_perind_numers.counts.formatted.gz

    """
}