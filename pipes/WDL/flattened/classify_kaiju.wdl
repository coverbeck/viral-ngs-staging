version 1.0



workflow classify_kaiju {
    call metagenomics__kaiju as kaiju
    output {
        File    kaiju_report           = kaiju.kaiju_report
        File    kaiju_reads            = kaiju.kaiju_reads
        File    krona_report_html      = kaiju.krona_report_html
        String  viral_classify_version = kaiju.viralngs_version
    }
}



task metagenomics__kaiju {
  input {
    File     reads_unmapped_bam
    File     kaiju_db_lz4  # <something>.fmi
    File     ncbi_taxonomy_db_tgz # taxonomy/{nodes.dmp, names.dmp}
    File     krona_taxonomy_db_tgz  # taxonomy/taxonomy.tab

    Int?     machine_mem_gb
    String   docker="quay.io/broadinstitute/viral-classify:2.0.21.3"
  }

  String   input_basename = basename(reads_unmapped_bam, ".bam")

  command {
    set -ex -o pipefail

    if [ -z "$TMPDIR" ]; then
      export TMPDIR=$(pwd)
    fi
    DB_DIR=$(mktemp -d --suffix _db)
    mkdir -p $DB_DIR/kaiju $DB_DIR/krona $DB_DIR/taxonomy

    lz4 -d ${kaiju_db_lz4} $DB_DIR/kaiju_db/kaiju.fmi

    read_utils.py extract_tarball \
      ${ncbi_taxonomy_db_tgz} $DB_DIR/taxonomy \
      --loglevel=DEBUG

    read_utils.py extract_tarball \
      ${krona_taxonomy_db_tgz} $DB_DIR/krona \
      --loglevel=DEBUG

    metagenomics.py --version | tee VERSION

    # classify contigs
    metagenomics.py kaiju \
      ${reads_unmapped_bam} \
      $DB_DIR/kaiju/kaiju.fmi \
      $DB_DIR/taxonomy \
      ${input_basename}.kaiju.summary_report.txt \
      --outReads ${input_basename}.kaiju.reads.txt.gz \
      --loglevel=DEBUG

    # run krona
    metagenomics.py krona \
      ${input_basename}.kaiju.summary_report.txt \
      $DB_DIR/krona \
      ${input_basename}.kaiju-krona.html \
      --inputType kaiju \
      --noRank --noHits \
      --loglevel=DEBUG
  }

  output {
    File    kaiju_report       = "${input_basename}.kaiju-summary_report.txt"
    File    kaiju_reads        = "${input_basename}.kaiju-reads.txt.gz"
    File    krona_report_html  = "${input_basename}.kaiju-krona.html"
    String  viralngs_version   = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: select_first([machine_mem_gb, 100]) + " GB"
    cpu: 16
    disks: "local-disk 375 LOCAL"
    dx_instance_type: "mem3_ssd1_v2_x16"
  }
}


