version 1.0




workflow isnvs_merge_to_vcf {
    input {
        File          reference_fasta
        Array[File]+  assemblies_fasta     # one per genome
    }

    call interhost__multi_align_mafft_ref as mafft {
        input:
            reference_fasta = reference_fasta,
            assemblies_fasta = assemblies_fasta
    }

    call tasks_intrahost__isnvs_vcf as isnvs_vcf {
        input:
            perSegmentMultiAlignments = mafft.alignments_by_chr,
            reference_fasta = reference_fasta
    }

    output {
        Array[File] alignments_by_chr   = mafft.alignments_by_chr
        File        isnvs_plain_vcf     = isnvs_vcf.isnvs_vcf
        File        isnvs_plain_vcf_idx = isnvs_vcf.isnvs_vcf_idx
        File        isnvs_annot_vcf     = isnvs_vcf.isnvs_annot_vcf
        File        isnvs_annot_vcf_idx = isnvs_vcf.isnvs_annot_vcf_idx
        File        isnvs_annot_txt     = isnvs_vcf.isnvs_annot_txt
        String      viral_phylo_version = mafft.viralngs_version
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




task tasks_intrahost__isnvs_vcf {
  input {
    Array[File]    vphaser2Calls
    Array[File]    perSegmentMultiAlignments
    File           reference_fasta

    Array[String]? snpEffRef
    Array[String]? sampleNames
    String?        emailAddress
    Boolean        naiveFilter=false

    Int?           machine_mem_gb
    String         docker="quay.io/broadinstitute/viral-phylo:2.0.21.0"
  }

  parameter_meta {
    vphaser2Calls:             { description: "vphaser output; ex. vphaser2.<sample>.txt.gz" }
    perSegmentMultiAlignments: { description: "aligned_##.fasta, where ## is segment number" }
    snpEffRef:                 { description: "list of accessions to build/find snpEff database" }
    sampleNames:               { description: "list of sample names" }
    emailAddress:              { description: "email address passed to NCBI if we need to download reference sequences" }
  }

  command {
    set -ex -o pipefail

    intrahost.py --version | tee VERSION

    SAMPLES="${sep=' ' sampleNames}"
    if [ -n "$SAMPLES" ]; then SAMPLES="--samples $SAMPLES"; fi

    providedSnpRefAccessions="${sep=' ' snpEffRef}"
    if [ -n "$providedSnpRefAccessions" ]; then 
      snpRefAccessions="$providedSnpRefAccessions";
    else
      snpRefAccessions="$(python -c "from Bio import SeqIO; print(' '.join(list(s.id for s in SeqIO.parse('${reference_fasta}', 'fasta'))))")"
    fi

    echo "snpRefAccessions: $snpRefAccessions"

    intrahost.py merge_to_vcf \
        ${reference_fasta} \
        isnvs.vcf.gz \
        $SAMPLES \
        --isnvs ${sep=' ' vphaser2Calls} \
        --alignments ${sep=' ' perSegmentMultiAlignments} \
        --strip_chr_version \
        ${true="--naive_filter" false="" naiveFilter} \
        --parse_accession
        
    interhost.py snpEff \
        isnvs.vcf.gz \
        $snpRefAccessions \
        isnvs.annot.vcf.gz \
        ${'--emailAddress=' + emailAddress}

    intrahost.py iSNV_table \
        isnvs.annot.vcf.gz \
        isnvs.annot.txt.gz
  }

  output {
    File        isnvs_vcf           = "isnvs.vcf.gz"
    File        isnvs_vcf_idx       = "isnvs.vcf.gz.tbi"
    File        isnvs_annot_vcf     = "isnvs.annot.vcf.gz"
    File        isnvs_annot_vcf_idx = "isnvs.annot.vcf.gz.tbi"
    File        isnvs_annot_txt     = "isnvs.annot.txt.gz"
    String      viralngs_version    = read_string("VERSION")
  }
  runtime {
    docker: "${docker}"
    memory: select_first([machine_mem_gb, 4]) + " GB"
    dx_instance_type: "mem1_ssd1_v2_x4"
  }
}


