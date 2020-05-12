version 1.0




workflow genbank {

    input {
        File          reference_fasta
        Array[File]+  assemblies_fasta     # one per genome
    }

    call interhost__multi_align_mafft_ref as mafft {
        input:
            reference_fasta = reference_fasta,
            assemblies_fasta = assemblies_fasta
    }

    call ncbi__annot_transfer as annot {
        input:
            multi_aln_fasta = mafft.alignments_by_chr,
            reference_fasta = reference_fasta
    }
 
    call ncbi__prepare_genbank as prep_genbank {
        input:
            assemblies_fasta = assemblies_fasta,
            annotations_tbl = annot.transferred_feature_tables
    }

    output {
        Array[File] alignments_by_chr          = mafft.alignments_by_chr

        Array[File] transferred_feature_tables = annot.transferred_feature_tables

        Array[File] sequin_files               = prep_genbank.sequin_files
        Array[File] structured_comment_files   = prep_genbank.structured_comment_files
        Array[File] genbank_preview_files      = prep_genbank.genbank_preview_files
        Array[File] source_table_files         = prep_genbank.source_table_files
        Array[File] fasta_per_chr_files        = prep_genbank.fasta_per_chr_files
        Array[File] validation_files           = prep_genbank.validation_files
        File        errorSummary               = prep_genbank.errorSummary

        String      viral_phylo_version        = mafft.viralngs_version
    }

}



task interhost__multi_align_mafft_ref {
  input {
    File           reference_fasta
    Array[File]+   assemblies_fasta # fasta files, one per sample, multiple chrs per file okay
    String         fasta_basename = basename(reference_fasta, '.fasta')
    Int?           mafft_maxIters
    Float?         mafft_ep
    Float?         mafft_gapOpeningPenalty

    Int?           machine_mem_gb
    String         docker="quay.io/broadinstitute/viral-phylo:2.0.21.0"
  }

  command {
    interhost.py --version | tee VERSION
    interhost.py multichr_mafft \
      ${reference_fasta} ${sep=' ' assemblies_fasta} \
      . \
      ${'--ep=' + mafft_ep} \
      ${'--gapOpeningPenalty=' + mafft_gapOpeningPenalty} \
      ${'--maxiters=' + mafft_maxIters} \
      --outFilePrefix align_mafft-${fasta_basename} \
      --preservecase \
      --localpair \
      --sampleNameListFile align_mafft-${fasta_basename}-sample_names.txt \
      --loglevel DEBUG
  }

  output {
    #File         sampleNamesFile    = "align_mafft-${fasta_basename}-sample_names.txt"
    Array[File]+  alignments_by_chr  = glob("align_mafft-${fasta_basename}*.fasta")
    String        viralngs_version   = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: select_first([machine_mem_gb, 28]) + " GB"
    cpu: 8
    dx_instance_type: "mem2_ssd1_v2_x8"
  }
}




task ncbi__annot_transfer {
  input {
    Array[File]+ multi_aln_fasta
    File         reference_fasta
    Array[File]+ reference_feature_table

    Int?         machine_mem_gb
    String       docker="quay.io/broadinstitute/viral-phylo:2.0.21.0"
  }

  Array[Int]   chr_nums=range(length(multi_aln_fasta))

  parameter_meta {
    multi_aln_fasta:         { description: "fasta; multiple alignments of sample sequences for each chromosome" }
    reference_fasta:         { description: "fasta; all chromosomes in one file" }
    reference_feature_table: { description: "tbl; feature table corresponding to each chromosome in the alignment" }
  }

  command {
    set -ex -o pipefail
    ncbi.py --version | tee VERSION
    echo ${sep=' ' multi_aln_fasta} > alignments.txt
    echo ${sep=' ' reference_feature_table} > tbls.txt
    for i in ${sep=' ' chr_nums}; do
      _alignment_fasta=`cat alignments.txt | cut -f $(($i+1)) -d ' '`
      _feature_tbl=`cat tbls.txt | cut -f $(($i+1)) -d ' '`
      ncbi.py tbl_transfer_prealigned \
          $_alignment_fasta \
          ${reference_fasta} \
          $_feature_tbl \
          . \
          --oob_clip \
          --loglevel DEBUG
    done
  }

  output {
    Array[File] transferred_feature_tables = glob("*.tbl")
    String      viralngs_version           = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: select_first([machine_mem_gb, 3]) + " GB"
    cpu: 2
    dx_instance_type: "mem1_ssd1_v2_x2"
  }
}




task ncbi__prepare_genbank {
  input {
    Array[File]+ assemblies_fasta
    Array[File]  annotations_tbl
    File         authors_sbt
    File         biosampleMap
    File         genbankSourceTable
    File?        coverage_table # summary.assembly.txt (from Snakemake) -- change this to accept a list of mapped bam files and we can create this table ourselves
    String       sequencingTech
    String       comment # TO DO: make this optional
    String       organism
    String       molType = "cRNA"

    Int?         machine_mem_gb
    String       docker="quay.io/broadinstitute/viral-phylo:2.0.21.0"
  }

  command {
    set -ex -o pipefail
    cp ${sep=' ' annotations_tbl} .
    ncbi.py --version | tee VERSION
    ncbi.py prep_genbank_files \
        ${authors_sbt} \
        ${sep=' ' assemblies_fasta} \
        . \
        --mol_type ${molType} \
        --organism "${organism}" \
        --biosample_map ${biosampleMap} \
        --master_source_table ${genbankSourceTable} \
        ${'--coverage_table ' + coverage_table} \
        --comment "${comment}" \
        --sequencing_tech "${sequencingTech}" \
        --loglevel DEBUG
    mv errorsummary.val errorsummary.val.txt # to keep it separate from the glob
  }

  output {
    Array[File] sequin_files             = glob("*.sqn")
    Array[File] structured_comment_files = glob("*.cmt")
    Array[File] genbank_preview_files    = glob("*.gbf")
    Array[File] source_table_files       = glob("*.src")
    Array[File] fasta_per_chr_files      = glob("*.fsa")
    Array[File] validation_files         = glob("*.val")
    File        errorSummary             = "errorsummary.val.txt"
    String      viralngs_version         = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: select_first([machine_mem_gb, 3]) + " GB"
    cpu: 2
    dx_instance_type: "mem1_ssd1_v2_x2"
  }
}


