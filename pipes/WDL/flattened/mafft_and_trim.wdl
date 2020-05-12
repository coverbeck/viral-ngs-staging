version 1.0



workflow mafft_and_trim {
    call interhost__multi_align_mafft as mafft

    scatter(alignment in mafft.alignments_by_chr) {
        call interhost__trimal_clean_msa as trimal {
            input:
                in_aligned_fasta = alignment
        }
    }

    output {
        Array[File] alignments_by_chr          = mafft.alignments_by_chr
        Array[File] trimmed_alignements_by_chr = trimal.trimal_cleaned_fasta
        File        sampleNamesFile            = mafft.sampleNamesFile
        String      viral_phylo_version        = mafft.viralngs_version
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




task interhost__trimal_clean_msa {
  input {
    File     in_aligned_fasta

    Int?     machine_mem_gb
    String   docker="quay.io/biocontainers/trimal:1.4.1--h6bb024c_3"

    String   input_basename = basename(basename(in_aligned_fasta, ".fasta"), ".fa")
  }

  command {
    trimal -fasta -automated1 -in "${in_aligned_fasta}" -out "${input_basename}_trimal_cleaned.fasta"
  }

  output {
    File   trimal_cleaned_fasta = "${input_basename}_trimal_cleaned.fasta"
  }
  runtime {
    docker: "${docker}"
    memory: select_first([machine_mem_gb, 7]) + " GB"
    dx_instance_type: "mem1_ssd1_v2_x8"
  }
}


