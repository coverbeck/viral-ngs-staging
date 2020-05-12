version 1.0



workflow coverage_table {
    call reports__coverage_report as coverage_report
    output {
        File   coverage_report_txt = coverage_report.coverage_report
        String viral_core_version  = coverage_report.viralngs_version
    }
}



task reports__coverage_report {
  input {
    Array[File]+ mapped_bams
    Array[File]  mapped_bam_idx # optional.. speeds it up if you provide it, otherwise we auto-index
    String       out_report_name="coverage_report.txt"

    String       docker="quay.io/broadinstitute/viral-core:2.0.21"
  }

  command {
    reports.py --version | tee VERSION
    reports.py coverage_only \
      ${sep=' ' mapped_bams} \
      ${out_report_name} \
      --loglevel DEBUG
  }

  output {
    File   coverage_report  = "${out_report_name}"
    String viralngs_version = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: "2 GB"
    cpu: 2
    disks: "local-disk 375 LOCAL"
    dx_instance_type: "mem1_ssd2_v2_x4"
  }
}


