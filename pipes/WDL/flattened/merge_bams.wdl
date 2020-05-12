version 1.0



workflow merge_bams {
    call read_utils__merge_and_reheader_bams as merge_and_reheader_bams
    output {
        File   out_bam            = merge_and_reheader_bams.out_bam
        String viral_core_version = merge_and_reheader_bams.viralngs_version
    }   
}



task read_utils__merge_and_reheader_bams {
    meta {
      description: "Merge and/or reheader bam files using a mapping table. This task can modify read group tags in a BAM header field for single BAM files or as part of a BAM merge operation. The output is a single BAM file (given one or more input BAMs) and a three-column tab delimited text table that defines: the field, the old value, and the new value (e.g. LB, old_lib_name, new_lib_name or SM, old_sample_name, new_sample_name)"
    }

    input {
      Array[File]+    in_bams
      String?         sample_name
      File?           reheader_table
      String          out_basename

      String          docker="quay.io/broadinstitute/viral-core:2.0.21"
    }

    command {
        set -ex -o pipefail

        read_utils.py --version | tee VERSION
        mem_in_mb=$(/opt/viral-ngs/source/docker/calc_mem.py mb 90)

        if [ ${length(in_bams)} -gt 1 ]; then
            read_utils.py merge_bams ${sep=' ' in_bams} merged.bam --JVMmemory="$mem_in_mb"m --loglevel DEBUG
        else
            echo "Skipping merge, only one input file"
            cp ${sep=' ' in_bams} merged.bam
        fi    

        # remap all SM values to user specified value
        if [ -n "${sample_name}" ]; then
          # create sample name remapping table based on existing sample names
          samtools view -H merged.bam | perl -n -e'/SM:(\S+)/ && print "SM\t$1\t'"${sample_name}"'\n"' | sort | uniq >> reheader_table.txt
        fi

        # remap arbitrary headers using user specified table
        if [[ -f "${reheader_table}" ]]; then
          cat "${reheader_table}" >> reheader_table.txt
        fi

        # reheader bam file if requested
        if [ -s reheader_table.txt ]; then
          read_utils.py reheader_bam merged.bam reheader_table.txt "${out_basename}.bam" --loglevel DEBUG
        else
          mv merged.bam "${out_basename}.bam"
        fi
    }

    output {
        File   out_bam          = "${out_basename}.bam"
        String viralngs_version = read_string("VERSION")
    }

    runtime {
        docker: "${docker}"
        memory: "3 GB"
        cpu: 2
        disks: "local-disk 750 LOCAL"
        dx_instance_type: "mem1_ssd2_v2_x4"
        preemptible: 0
    }
}


