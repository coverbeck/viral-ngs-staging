version 1.0



workflow classify_kraken2 {
    call metagenomics__kraken2 as kraken2

    output {
        File    kraken2_reads_report   = kraken2.kraken2_reads_report
        File    kraken2_summary_report = kraken2.kraken2_summary_report
        File    krona_report_html      = kraken2.krona_report_html
        String  viral_classify_version = kraken2.viralngs_version
    }
}



task metagenomics__kraken2 {
  meta {
    description: "Runs Kraken2 classification"
  }

  input {
    File     reads_bam
    File     kraken2_db_tgz         # {database.kdb,taxonomy}
    File     krona_taxonomy_db_tgz  # taxonomy.tab
    Float?   confidence_threshold
    Int?     min_base_qual

    Int?     machine_mem_gb
    String   docker="quay.io/broadinstitute/viral-classify:2.0.21.3"
  }

  parameter_meta {
    reads_bam: {
      description: "Reads to classify. May be unmapped or mapped or both, paired-end or single-end.",
      patterns: ["*.bam"] }
    kraken2_db_tgz: {
      description: "Pre-built Kraken database tarball containing three files: hash.k2d, opts.k2d, and taxo.k2d.",
      patterns: ["*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
    }
    krona_taxonomy_db_tgz: {
      description: "Krona taxonomy database containing a single file: taxonomy.tab, or possibly just a compressed taxonomy.tab",
      patterns: ["*.tab.zst", "*.tab.gz", "*.tab", "*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
    }
    confidence_threshold: {
      description: "Kraken2 confidence score threshold (0.0-1.0). See https://ccb.jhu.edu/software/kraken2/index.shtml?t=manual#confidence-scoring"
    }
    min_base_qual: {
      description: "Minimum base quality used in classification"
    }
  }

  String out_basename=basename(reads_bam, '.bam')

  command {
    set -ex -o pipefail

    if [ -z "$TMPDIR" ]; then
      export TMPDIR=$(pwd)
    fi
    DB_DIR=$(mktemp -d --suffix _db)
    mkdir -p $DB_DIR/kraken2 $DB_DIR/krona

    # decompress DB to $DB_DIR
    read_utils.py extract_tarball \
      ${kraken2_db_tgz} $DB_DIR/kraken2 \
      --loglevel=DEBUG
    du -hs $DB_DIR/kraken2

    # unpack krona taxonomy.tab
    if [[ ${krona_taxonomy_db_tgz} == *.tar.* ]]; then
      read_utils.py extract_tarball \
        ${krona_taxonomy_db_tgz} $DB_DIR/krona \
        --loglevel=DEBUG &  # we don't need this until later
    else
      if [[ "${krona_taxonomy_db_tgz}" == *.zst ]]; then
        cat "${krona_taxonomy_db_tgz}" | zstd -d > $DB_DIR/krona/taxonomy.tab &
      elif [[ "${krona_taxonomy_db_tgz}" == *.gz ]]; then
        cat "${krona_taxonomy_db_tgz}" | pigz -dc > $DB_DIR/krona/taxonomy.tab &
      elif [[ "${krona_taxonomy_db_tgz}" == *.bz2 ]]; then
        cat "${krona_taxonomy_db_tgz}" | bzip -dc > $DB_DIR/krona/taxonomy.tab &
      else
        cp "${krona_taxonomy_db_tgz}" $DB_DIR/krona/taxonomy.tab &
      fi
    fi

    metagenomics.py --version | tee VERSION

    metagenomics.py kraken2 \
      $DB_DIR/kraken2 \
      ${reads_bam} \
      --outReads   "${out_basename}".kraken2.reads.txt \
      --outReports "${out_basename}".kraken2.report.txt \
      ${"--confidence " + confidence_threshold} \
      ${"--min_base_qual " + min_base_qual} \
      --loglevel=DEBUG

    wait # for krona_taxonomy_db_tgz to download and extract
    pigz "${out_basename}".kraken2.reads.txt &

    metagenomics.py krona \
      "${out_basename}".kraken2.report.txt \
      $DB_DIR/krona \
      "${out_basename}".kraken2.krona.html \
      --sample_name "${out_basename}" \
      --noRank --noHits --inputType kraken2 \
      --loglevel=DEBUG

    wait # pigz reads.txt
  }

  output {
    File    kraken2_reads_report   = "${out_basename}.kraken2.reads.txt.gz"
    File    kraken2_summary_report = "${out_basename}.kraken2.report.txt"
    File    krona_report_html      = "${out_basename}.kraken2.krona.html"
    String  viralngs_version       = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: select_first([machine_mem_gb, 52]) + " GB"
    cpu: 8
    disks: "local-disk 750 LOCAL"
    dx_instance_type: "mem3_ssd1_v2_x8"
    preemptible: 2
  }
}


