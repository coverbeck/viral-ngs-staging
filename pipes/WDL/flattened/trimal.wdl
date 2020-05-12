version 1.0



workflow trimal {
    call interhost__trimal_clean_msa as trimal_clean_msa

    output {
        File trimmed_alignment = trimal_clean_msa.trimal_cleaned_fasta
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


