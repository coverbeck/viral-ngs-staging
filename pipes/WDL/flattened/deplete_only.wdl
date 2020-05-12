version 1.0



workflow deplete_only {
    call taxon_filter__deplete_taxa as deplete_taxa

    output {
        File   cleaned_bam               = deplete_taxa.cleaned_bam
        File   cleaned_fastqc            = deplete_taxa.cleaned_fastqc
        File   cleaned_fastqc_zip        = deplete_taxa.cleaned_fastqc_zip
        Int    depletion_read_count_pre  = deplete_taxa.depletion_read_count_pre
        Int    depletion_read_count_post = deplete_taxa.depletion_read_count_post
        String viral_classify_version    = deplete_taxa.viralngs_version
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


