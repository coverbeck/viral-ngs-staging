version 1.0



workflow fastq_to_ubam {
    call tasks_read_utils__FastqToUBAM as FastqToUBAM
    output {
        File unmapped_bam = FastqToUBAM.unmapped_bam
    }
}



task tasks_read_utils__FastqToUBAM {
  meta {
    description: "Converts FASTQ (paired or single) to uBAM and adds read group information."
  }
  input {
    File    fastq_1
    File?   fastq_2
    String  sample_name
    String  library_name
    String? readgroup_name
    String? platform_unit
    String? run_date
    String? platform_name
    String? sequencing_center

    String  docker="quay.io/broadinstitute/viral-core:2.0.21"
  }
  parameter_meta {
    fastq_1: { description: "Unaligned read1 file in fastq format", patterns: ["*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"] }
    fastq_2: { description: "Unaligned read2 file in fastq format. This should be empty for single-end read conversion and required for paired-end reads. If provided, it must match fastq_1 in length and order.", patterns: ["*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"] }
    sample_name: { description: "Sample name. This is required and will populate the 'SM' read group value and will be used as the output filename (must be filename-friendly)." }
    library_name: { description: "Library name. This is required and will populate the 'LB' read group value. SM & LB combinations must be identical for any sequencing reads generated from the same sequencing library, and must be distinct for any reads generated from different libraries." }
  }
  command {
      set -ex -o pipefail

      # find 90% memory
      mem_in_mb=`/opt/viral-ngs/source/docker/calc_mem.py mb 90`

      read_utils.py --version | tee VERSION

      picard -Xmx"$mem_in_mb"m \
        FastqToSam \
        FASTQ="~{fastq_1}" \
        ${"FASTQ2=" + fastq_2} \
        SAMPLE_NAME="${sample_name}" \
        LIBRARY_NAME="${library_name}" \
        OUTPUT="${sample_name}".bam \
        ${"READ_GROUP_NAME=" + readgroup_name} \
        ${"PLATFORM_UNIT=" + platform_unit} \
        ${"RUN_DATE=" + run_date} \
        ${"PLATFORM=" + platform_name} \
        ${"SEQUENCING_CENTER=" + sequencing_center}
  }
  runtime {
    docker: docker
    cpu: 2
    memory: "3 GB"
    disks: "local-disk 375 LOCAL"
    dx_instance_type: "mem1_ssd1_v2_x2"
  }
  output {
    File unmapped_bam = "~{sample_name}.bam"
  }
}


