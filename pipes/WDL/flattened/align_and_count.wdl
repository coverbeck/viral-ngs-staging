version 1.0



workflow align_and_count_report {
    call reports__align_and_count as align_and_count
    output {
        File report               = align_and_count.report
        File report_top_hits      = align_and_count.report_top_hits
        String viral_core_version = align_and_count.viralngs_version
    }
}



task reports__align_and_count {
  input {
    File    reads_bam
    File    ref_db
    Int?    minScoreToFilter
    Int?    topNHits = 3

    Int?    machine_mem_gb
    String  docker="quay.io/broadinstitute/viral-core:2.0.21"
  }

  String  reads_basename=basename(reads_bam, ".bam")
  String  ref_basename=basename(ref_db, ".fasta")

  command {
    set -ex -o pipefail

    read_utils.py --version | tee VERSION

    ln -s ${reads_bam} ${reads_basename}.bam
    read_utils.py bwamem_idxstats \
      ${reads_basename}.bam \
      ${ref_db} \
      --outStats ${reads_basename}.count.${ref_basename}.txt.unsorted \
      ${'--minScoreToFilter=' + minScoreToFilter} \
      --loglevel=DEBUG

      sort -b -r -n -k3 ${reads_basename}.count.${ref_basename}.txt.unsorted > ${reads_basename}.count.${ref_basename}.txt
      head -n ${topNHits} ${reads_basename}.count.${ref_basename}.txt > ${reads_basename}.count.${ref_basename}.top_${topNHits}_hits.txt
  }

  output {
    File   report           = "${reads_basename}.count.${ref_basename}.txt"
    File   report_top_hits  = "${reads_basename}.count.${ref_basename}.top_${topNHits}_hits.txt"
    String viralngs_version = read_string("VERSION")
  }

  runtime {
    memory: select_first([machine_mem_gb, 7]) + " GB"
    cpu: 4
    docker: "${docker}"
    disks: "local-disk 375 LOCAL"
    dx_instance_type: "mem1_ssd1_v2_x4"
  }
}


