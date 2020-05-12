version 1.0

import "../tasks/tasks_read_utils.wdl" as read_utils

workflow merge_bams {
    call read_utils.merge_and_reheader_bams
    output {
        File   out_bam            = merge_and_reheader_bams.out_bam
        String viral_core_version = merge_and_reheader_bams.viralngs_version
    }   
}
