version 1.0






workflow contigs {
    input {
        File         reads_unmapped_bam

        Array[File]  deplete_bmtaggerDbs = []
        Array[File]  deplete_blastDbs = []
        Array[File]  deplete_bwaDbs = []
    }

    if(length(deplete_bmtaggerDbs) + length(deplete_blastDbs) + length(deplete_bwaDbs) > 0) {
        call taxon_filter__deplete_taxa as deplete {
          input:
            raw_reads_unmapped_bam = reads_unmapped_bam,
            bmtaggerDbs = deplete_bmtaggerDbs,
            blastDbs = deplete_blastDbs,
            bwaDbs = deplete_bwaDbs
        }
    }

    call read_utils__rmdup_ubam as rmdup_ubam {
        input:
            reads_unmapped_bam = select_first([deplete.cleaned_bam, reads_unmapped_bam])
    }

    call assembly__assemble as spades {
        input:
            assembler = "spades",
            reads_unmapped_bam = rmdup_ubam.dedup_bam
    }

    # TO DO: taxonomic classification of contigs

    output {
        File?  cleaned_reads_unaligned_bam     = deplete.cleaned_bam
        File   deduplicated_reads_unaligned    = rmdup_ubam.dedup_bam
        File   contigs_fastas                  = spades.contigs_fasta

        Int?   read_counts_raw                 = deplete.depletion_read_count_pre
        Int?   read_counts_depleted            = deplete.depletion_read_count_post
        Int    read_counts_dedup               = rmdup_ubam.dedup_read_count_post
        Int    read_counts_prespades_subsample = spades.subsample_read_count

        String? deplete_viral_classify_version  = deplete.viralngs_version
        String  spades_viral_assemble_version   = spades.viralngs_version
   }
}



task taxon_filter__deplete_taxa {
  meta { description: "Runs a full human read depletion pipeline and removes PCR duplicates. Input database files (bmtaggerDbs, blastDbs, bwaDbs) may be any combination of: .fasta, .fasta.gz, or tarred up indexed fastas (using the software's indexing method) as .tar.gz, .tar.bz2, .tar.lz4, or .tar.zst." }

  input {
    File         raw_reads_unmapped_bam
    Array[File]? bmtaggerDbs
    Array[File]? blastDbs
    Array[File]? bwaDbs
    Int?         query_chunk_size
    Boolean?     clear_tags = false
    String?      tags_to_clear_space_separated = "XT X0 X1 XA AM SM BQ CT XN OC OP"

    Int?         cpu=8
    Int?         machine_mem_gb
    String       docker="quay.io/broadinstitute/viral-classify:2.0.21.3"
  }

  parameter_meta {
    raw_reads_unmapped_bam: { description: "unaligned reads in BAM format", patterns: ["*.bam"] }
    bmtaggerDbs: {
       description: "Optional list of databases to use for bmtagger-based depletion. Sequences in fasta format will be indexed on the fly, pre-bmtagger-indexed databases may be provided as tarballs.",
       patterns: ["*.fasta", "*.fasta.gz", "*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
    }
    blastDbs: {
      description: "Optional list of databases to use for blastn-based depletion. Sequences in fasta format will be indexed on the fly, pre-blast-indexed databases may be provided as tarballs.",
      patterns: ["*.fasta", "*.fasta.gz", "*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
    }
    bwaDbs: {
      description: "Optional list of databases to use for bwa mem-based depletion. Sequences in fasta format will be indexed on the fly, pre-bwa-indexed databases may be provided as tarballs.",
      patterns: ["*.fasta", "*.fasta.gz", "*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
    }
  }

  String       bam_basename = basename(raw_reads_unmapped_bam, ".bam")

  command {
    set -ex -o pipefail
    taxon_filter.py --version | tee VERSION

    if [ -z "$TMPDIR" ]; then
      export TMPDIR=$(pwd)
    fi

    # find memory thresholds
    mem_in_mb_50=`/opt/viral-ngs/source/docker/calc_mem.py mb 50`
    mem_in_mb_75=`/opt/viral-ngs/source/docker/calc_mem.py mb 75`

    # bmtagger and blast db args
    DBS_BMTAGGER="${sep=' ' bmtaggerDbs}"
    DBS_BLAST="${sep=' ' blastDbs}"
    DBS_BWA="${sep=' ' bwaDbs}"
    if [ -n "$DBS_BMTAGGER" ]; then DBS_BMTAGGER="--bmtaggerDbs $DBS_BMTAGGER"; fi
    if [ -n "$DBS_BLAST" ]; then DBS_BLAST="--blastDbs $DBS_BLAST"; fi
    if [ -n "$DBS_BWA" ]; then DBS_BWA="--bwaDbs $DBS_BWA"; fi
    
    if [[ "${clear_tags}" == "true" ]]; then
      TAGS_TO_CLEAR="--clearTags"
      if [[ -n "${tags_to_clear_space_separated}" ]]; then
        TAGS_TO_CLEAR="$TAGS_TO_CLEAR ${'--tagsToClear=' + tags_to_clear_space_separated}"
      fi
    fi

    # run depletion
    taxon_filter.py deplete \
      ${raw_reads_unmapped_bam} \
      tmpfile.raw.bam \
      tmpfile.bwa.bam \
      tmpfile.bmtagger_depleted.bam \
      ${bam_basename}.cleaned.bam \
      $DBS_BMTAGGER $DBS_BLAST $DBS_BWA \
      ${'--chunkSize=' + query_chunk_size} \
      $TAGS_TO_CLEAR \
      --JVMmemory="$mem_in_mb_50"m \
      --srprismMemory=$mem_in_mb_75 \
      --loglevel=DEBUG

    samtools view -c ${raw_reads_unmapped_bam} | tee depletion_read_count_pre
    samtools view -c ${bam_basename}.cleaned.bam | tee depletion_read_count_post
    reports.py fastqc ${bam_basename}.cleaned.bam ${bam_basename}.cleaned_fastqc.html --out_zip ${bam_basename}.cleaned_fastqc.zip
  }

  output {
    File   cleaned_bam               = "${bam_basename}.cleaned.bam"
    File   cleaned_fastqc            = "${bam_basename}.cleaned_fastqc.html"
    File   cleaned_fastqc_zip        = "${bam_basename}.cleaned_fastqc.zip"
    Int    depletion_read_count_pre  = read_int("depletion_read_count_pre")
    Int    depletion_read_count_post = read_int("depletion_read_count_post")
    String viralngs_version          = read_string("VERSION")
  }
  runtime {
    docker: "${docker}"
    memory: select_first([machine_mem_gb, 15]) + " GB"
    cpu: "${cpu}"
    disks: "local-disk 375 LOCAL"
    dx_instance_type: "mem1_ssd1_v2_x8"
    preemptible: 1
  }
}




task read_utils__rmdup_ubam {
  meta {
    description: "Perform read deduplication on unaligned reads."
  }

  input {
    File     reads_unmapped_bam
    String   method="mvicuna"

    Int?     machine_mem_gb
    String?  docker="quay.io/broadinstitute/viral-core:2.0.21"
  }

  parameter_meta {
    reads_unmapped_bam: { description: "unaligned reads in BAM format", patterns: ["*.bam"] }
    method:             { description: "mvicuna or cdhit" }
  }

  String reads_basename = basename(reads_unmapped_bam, ".bam")
  
  command {
    set -ex -o pipefail
    mem_in_mb=$(/opt/viral-ngs/source/docker/calc_mem.py mb 90)
    read_utils.py --version | tee VERSION

    read_utils.py rmdup_"${method}"_bam \
      "${reads_unmapped_bam}" \
      "${reads_basename}".dedup.bam \
      --JVMmemory "$mem_in_mb"m \
      --loglevel=DEBUG

    samtools view -c ${reads_basename}.dedup.bam | tee dedup_read_count_post
    reports.py fastqc ${reads_basename}.dedup.bam ${reads_basename}.dedup_fastqc.html --out_zip ${reads_basename}.dedup_fastqc.zip
  }

  output {
    File   dedup_bam             = "${reads_basename}.dedup.bam"
    File   dedup_fastqc          = "${reads_basename}.dedup_fastqc.html"
    File   dedup_fastqc_zip      = "${reads_basename}.dedup_fastqc.zip"
    Int    dedup_read_count_post = read_int("dedup_read_count_post")
    String viralngs_version      = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: select_first([machine_mem_gb, 7]) + " GB"
    cpu:    2
    disks:  "local-disk 375 LOCAL"
    dx_instance_type: "mem1_ssd1_v2_x2"
  }
}




task assembly__assemble {
    input {
      File     reads_unmapped_bam
      File     trim_clip_db

      Int?     trinity_n_reads=250000
      Int?     spades_n_reads=10000000
      Int?     spades_min_contig_len=0

      String?  assembler="trinity"  # trinity, spades, or trinity-spades
      Boolean? always_succeed=false

      # do this in two steps in case the input doesn't actually have "taxfilt" in the name
      String   sample_name = basename(basename(reads_unmapped_bam, ".bam"), ".taxfilt")

      Int?     machine_mem_gb
      String   docker="quay.io/broadinstitute/viral-assemble:2.0.21.0"
    }

    command {
        set -ex -o pipefail

        # find 90% memory
        mem_in_mb=`/opt/viral-ngs/source/docker/calc_mem.py mb 90`
        mem_in_gb=`/opt/viral-ngs/source/docker/calc_mem.py gb 90`

        assembly.py --version | tee VERSION

        if [[ "${assembler}" == "trinity" ]]; then
          assembly.py assemble_trinity \
            ${reads_unmapped_bam} \
            ${trim_clip_db} \
            ${sample_name}.assembly1-${assembler}.fasta \
            ${'--n_reads=' + trinity_n_reads} \
            ${true='--alwaysSucceed' false="" always_succeed} \
            --JVMmemory "$mem_in_mb"m \
            --outReads=${sample_name}.subsamp.bam \
            --loglevel=DEBUG

        elif [[ "${assembler}" == "spades" ]]; then
          assembly.py assemble_spades \
            ${reads_unmapped_bam} \
            ${trim_clip_db} \
            ${sample_name}.assembly1-${assembler}.fasta \
            ${'--nReads=' + spades_n_reads} \
            ${true="--alwaysSucceed" false="" always_succeed} \
            ${'--minContigLen=' + spades_min_contig_len} \
            --memLimitGb $mem_in_gb \
            --outReads=${sample_name}.subsamp.bam \
            --loglevel=DEBUG

        elif [[ "${assembler}" == "trinity-spades" ]]; then
          assembly.py assemble_trinity \
            ${reads_unmapped_bam} \
            ${trim_clip_db} \
            ${sample_name}.assembly1-trinity.fasta \
            ${'--n_reads=' + trinity_n_reads} \
            --JVMmemory "$mem_in_mb"m \
            --outReads=${sample_name}.subsamp.bam \
            ${true='--always_succeed' false='' always_succeed} \
            --loglevel=DEBUG
          assembly.py assemble_spades \
            ${reads_unmapped_bam} \
            ${trim_clip_db} \
            ${sample_name}.assembly1-${assembler}.fasta \
            --contigsUntrusted=${sample_name}.assembly1-trinity.fasta \
            ${'--nReads=' + spades_n_reads} \
            ${true='--alwaysSucceed' false='' always_succeed} \
            ${'--minContigLen=' + spades_min_contig_len} \
            --memLimitGb $mem_in_gb \
            --loglevel=DEBUG

        else
          echo "unrecognized assembler ${assembler}" >&2
          exit 1
        fi

        samtools view -c ${sample_name}.subsamp.bam | tee subsample_read_count >&2
    }

    output {
        File   contigs_fasta        = "${sample_name}.assembly1-${assembler}.fasta"
        File   subsampBam           = "${sample_name}.subsamp.bam"
        Int    subsample_read_count = read_int("subsample_read_count")
        String viralngs_version     = read_string("VERSION")
    }

    runtime {
        docker: "${docker}"
        memory: select_first([machine_mem_gb, 15]) + " GB"
        cpu: 4
        disks: "local-disk 375 LOCAL"
        dx_instance_type: "mem1_ssd1_v2_x8"
    }

}


