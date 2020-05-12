version 1.0



workflow mafft {
    call interhost__multi_align_mafft as multi_align_mafft
    output {
        File        sampleNamesFile     = multi_align_mafft.sampleNamesFile
        Array[File] alignments_by_chr   = multi_align_mafft.alignments_by_chr
        String      viral_phylo_version = multi_align_mafft.viralngs_version
    }
}



task interhost__multi_align_mafft {
  input {
    Array[File]+   assemblies_fasta # fasta files, one per sample, multiple chrs per file okay
    String         out_prefix = "aligned"
    Int?           mafft_maxIters
    Float?         mafft_ep
    Float?         mafft_gapOpeningPenalty

    Int?           machine_mem_gb
    String         docker="quay.io/broadinstitute/viral-phylo:2.0.21.0"
  }

  command {
    interhost.py --version | tee VERSION
    interhost.py multichr_mafft \
      ${sep=' ' assemblies_fasta} \
      . \
      ${'--ep=' + mafft_ep} \
      ${'--gapOpeningPenalty=' + mafft_gapOpeningPenalty} \
      ${'--maxiters=' + mafft_maxIters} \
      --outFilePrefix ${out_prefix} \
      --preservecase \
      --localpair \
      --sampleNameListFile ${out_prefix}-sample_names.txt \
      --loglevel DEBUG
  }

  output {
    File        sampleNamesFile   = "${out_prefix}-sample_names.txt"
    Array[File] alignments_by_chr = glob("${out_prefix}*.fasta")
    String      viralngs_version  = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: select_first([machine_mem_gb, 7]) + " GB"
    cpu: 8
    dx_instance_type: "mem1_ssd1_v2_x8"
  }
}


