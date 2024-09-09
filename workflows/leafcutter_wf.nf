nextflow.enable.dsl=2

include { bam_to_junc; cluster_introns} from '../modules/leafcutter'

workflow quant_leafcutter {
    take:
        bam_sorted_indexed
    
    main:
        bam_to_junc(bam_sorted_indexed)

        cluster_introns(bam_to_junc.out.junc.map { it.toString() }.collectFile(name: 'junction_files.txt', newLine: true),
        bam_to_junc.out.hg38_junc.map { it.toString() }.collectFile(name: 'junction_files_hg38.txt', newLine: true),
        bam_to_junc.out.gg6_junc.map { it.toString() }.collectFile(name: 'junction_files_gg6.txt', newLine: true)
        )
//         cluster_introns(bam_to_junc.out.junc.map{it.toString()}.collectFile(name: 'junction_files.txt', newLine: true))
}

