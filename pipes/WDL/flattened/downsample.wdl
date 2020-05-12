version 1.0



workflow downsample {
    call reads__downsample_bams as downsample_bams

    output {
      Array[File] downsampled_bam  = downsample_bams.downsampled_bam
      String      viralngs_version = downsample_bams.viralngs_version
    }
}



task reads__downsample_bams {
  meta {
    description: "Downsample reads in a BAM file randomly subsampling to a target read count. Read deduplication can occur either before or after random subsampling, or not at all (default: not at all)."
  }

  input {
    Array[File]  reads_bam
    Int?         readCount
    Boolean?     deduplicateBefore=false
    Boolean?     deduplicateAfter=false

    Int?         machine_mem_gb
    String       docker="quay.io/broadinstitute/viral-core:2.0.21"
  }

  command {
    set -ex -o pipefail

    # find 90% memory
    mem_in_mb=`/opt/viral-ngs/source/docker/calc_mem.py mb 90`

    if [[ "${deduplicateBefore}" == "true" ]]; then
      DEDUP_OPTION="--deduplicateBefore"
    elif [[ "${deduplicateAfter}" == "true" ]]; then
      DEDUP_OPTION="--deduplicateAfter"
    fi

    if [[ "${deduplicateBefore}" == "true" && "${deduplicateAfter}" == "true" ]]; then
      echo "deduplicateBefore and deduplicateAfter are mutually exclusive. Only one can be used."
      exit 1
    fi
    
    read_utils.py --version | tee VERSION

    read_utils.py downsample_bams \
        ${sep=' ' reads_bam} \
        --outPath ./output \
        ${'--readCount=' + readCount} \
        $DEDUP_OPTION \
        --JVMmemory "$mem_in_mb"m
  }

  output {
    Array[File] downsampled_bam  = glob("output/*.downsampled-*.bam")
    String      viralngs_version = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: select_first([machine_mem_gb, 3]) + " GB"
    cpu:    4
    disks:  "local-disk 750 LOCAL"
    dx_instance_type: "mem1_ssd1_v2_x4"
  }
}


