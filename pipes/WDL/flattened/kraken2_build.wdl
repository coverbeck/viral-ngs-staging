version 1.0



workflow kraken2_build {
    call metagenomics__build_kraken2_db as build_kraken2_db

    output {
        File kraken2_db   = build_kraken2_db.kraken2_db
        File taxdump_tgz  = build_kraken2_db.taxdump_tgz
        File krona_db     = build_kraken2_db.krona_db
    }
}



task metagenomics__build_kraken2_db {
  meta {
    description: "Builds a custom kraken2 database. Outputs tar.zst tarballs of kraken2 database, associated krona taxonomy db, and an ncbi taxdump.tar.gz. Requires live internet access if any standard_libraries are specified or if taxonomy_db_tgz is absent."
  }

  input {
    String        db_basename
    File?         taxonomy_db_tgz
    Array[String] standard_libraries = [
                      "archaea", "bacteria", "plasmid",
                      "viral", "human", "fungi", "protozoa",
                      "UniVec_Core"]
    Array[File]   custom_libraries = []
    Boolean       protein=false

    Int?        kmerLen
    Int?        minimizerLen
    Int?        minimizerSpaces
    Int?        maxDbSize
    Int?        zstd_compression_level

    Int?        machine_mem_gb
    String      docker="quay.io/broadinstitute/viral-classify:2.0.21.3"
  }

  parameter_meta {
    db_basename: { description: "A descriptive string used in output filenames. Outputs will be called kraken2-<db_basename>.tar.zst, krona-<db_basename>.tar.zst, and taxdump-<db_basename>.tar.gz" }
    taxonomy_db_tgz: {
       description: "Optional tarball of kraken2 taxonomy database directory. Omitting this input will cause a fresh download from NCBI at the time of build.",
       patterns: ["*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
    }
    standard_libraries: {
      description: "A list of 'standard' kraken2 databases to include in this build. Including any values here will cause fresh downloads of data at the time of build. A list of acceptable names is available at https://ccb.jhu.edu/software/kraken2/index.shtml?t=manual#custom-databases"
    }
    custom_libraries: {
      description: "A list of 'custom' kraken2 databases to include in this build. Headers must be formatted as described in the kraken2 documentation. These are fastas or tarball collections of such fastas--multiple may be provided here.",
      patterns: ["*.fasta", "*.fasta.zst", "*.fasta.gz", "*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
    }
    protein: {
      description: "Build a protein (translated search) database. Default is nucleotide."
    }
    kmerLen: {
      description: "(k) K-mer length in bp/aa (Kraken2 defaults: 35 nt, 15 aa)"
    }
    minimizerLen: {
      description: "(l) Minimizer length in bp/aa (Kraken2 defaults: 31 nt, 12 aa)"
    }
    minimizerSpaces: {
      description: "(s) Number of bp/aa in minimizer that are ignored in comparisons (Kraken2 defaults: 7 nt, 0 aa)"
    }
  }

  command {
    set -ex -o pipefail

    if [ -z "$TMPDIR" ]; then
      export TMPDIR=$(pwd)
    fi
    TAXDB_DIR=$(mktemp -d)
    FASTAS_DIR=$(mktemp -d)
    KRONA_DIR=$(mktemp -d)
    DB_DIR=$(mktemp -d)

    metagenomics.py --version | tee VERSION

    # prep input taxonomy db, if specified
    if [ -n "${taxonomy_db_tgz}" ]; then
      read_utils.py extract_tarball \
        ${taxonomy_db_tgz} $TAXDB_DIR \
        --loglevel=DEBUG &
      TAX_INPUT_CMD="--tax_db=$TAXDB_DIR"
    else
      TAX_INPUT_CMD=""
    fi

    # prep input custom fastas, if specified
    CUSTOM_INPUT_CMD=""
    if [ -n "${sep=' ' custom_libraries}" ]; then
      CUSTOM_INPUT_CMD="--custom_libraries "
      for TGZ in ${sep=' ' custom_libraries}; do
        if [[ ($TGZ == *.tar.*) || ($TGZ == *.tgz) ]]; then
          read_utils.py extract_tarball \
            $TGZ $FASTAS_DIR \
            --loglevel=DEBUG &
        else
          if [[ $TGZ == *.zst ]]; then
            cat $TGZ | zstd -d > $FASTAS_DIR/$(basename $TGZ .zst) &
          elif [[ $TGZ == *.gz ]]; then
            cat $TGZ | pigz -dc > $FASTAS_DIR/$(basename $TGZ .gz) &
          elif [[ $TGZ == *.bz2 ]]; then
            cat $TGZ | bzip -dc > $FASTAS_DIR/$(basename $TGZ .bz2) &
          else
            CUSTOM_INPUT_CMD="$CUSTOM_INPUT_CMD $TGZ"
          fi
        fi
      done
      wait # wait for all decompressions to finish
      for FASTA in $FASTAS_DIR/*; do
        CUSTOM_INPUT_CMD="$CUSTOM_INPUT_CMD $FASTA"
      done
    fi

    # prep standard libraries, if specified
    STD_INPUT_CMD=""
    if [ -n "${sep=' ' standard_libraries}" ]; then
      STD_INPUT_CMD="--standard_libraries ${sep=' ' standard_libraries}"
    fi

    # build kraken2 database
    wait # wait for all decompressions to finish
    metagenomics.py kraken2_build \
      $DB_DIR \
      $TAX_INPUT_CMD \
      $STD_INPUT_CMD \
      $CUSTOM_INPUT_CMD \
      --taxdump_out "taxdump-${db_basename}.tar.gz" \
      ${true='--protein' false='' protein} \
      ${'--kmerLen=' + kmerLen} \
      ${'--minimizerLen=' + minimizerLen} \
      ${'--minimizerSpaces=' + minimizerSpaces} \
      ${'--maxDbSize=' + maxDbSize}
      --loglevel=DEBUG
    tar -c -C $DB_DIR . | zstd ${"-" + zstd_compression_level} > "kraken2-${db_basename}.tar.zst" &

    # build matching krona db
    metagenomics.py krona_build \
      $KRONA_DIR --taxdump_tar_gz "taxdump-${db_basename}.tar.gz"
    cat $KRONA_DIR/taxonomy.tab | zstd -19 > "krona-${db_basename}-taxonomy.tab.zst"

    wait # tar/zst of kraken2 db
  }

  output {
    File        kraken2_db       = "kraken2-${db_basename}.tar.zst"
    File        taxdump_tgz      = "taxdump-${db_basename}.tar.gz"
    File        krona_db         = "krona-${db_basename}-taxonomy.tab.zst"
    String      viralngs_version = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: select_first([machine_mem_gb, 100]) + " GB"
    disks: "local-disk 750 LOCAL"
    cpu: 16
    dx_instance_type: "mem3_ssd1_v2_x16"
    preemptible: 0
  }
}


