version 1.0







workflow classify_multi {
    meta {
         description: "Runs raw reads through taxonomic classification (Kraken2), human read depletion (based on Kraken2), de novo assembly (SPAdes), taxonomic classification of contigs (BLASTx), and FASTQC/multiQC of reads."
    }

    input {
        Array[File]+ reads_bams

        File  ncbi_taxdump_tgz

        File  spikein_db
        File  trim_clip_db

        File  kraken2_db_tgz
        File  krona_taxonomy_db_kraken2_tgz
        File? blast_db_tgz
        File? krona_taxonomy_db_blast_tgz
    }

    parameter_meta {
        reads_bams: {
          description: "Reads to classify. May be unmapped or mapped or both, paired-end or single-end.",
          patterns: ["*.bam"]
        }
        spikein_db: {
          description: "ERCC spike-in sequences",
          patterns: ["*.fasta", "*.fasta.gz", "*.fasta.zst"]
        }
        trim_clip_db: {
          description: "Adapter sequences to remove via trimmomatic prior to SPAdes assembly",
          patterns: ["*.fasta", "*.fasta.gz", "*.fasta.zst"]
        }
        kraken2_db_tgz: {
          description: "Pre-built Kraken database tarball containing three files: hash.k2d, opts.k2d, and taxo.k2d.",
          patterns: ["*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
        }
        krona_taxonomy_db_kraken2_tgz: {
          description: "Krona taxonomy database containing a single file: taxonomy.tab, or possibly just a compressed taxonomy.tab",
          patterns: ["*.tab.zst", "*.tab.gz", "*.tab", "*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
        }
        blast_db_tgz: {
          description: "Pre-built BLAST database tarball containing an indexed blast database named 'nr'",
          patterns: ["*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
        }
        krona_taxonomy_db_blast_tgz: {
          description: "Krona taxonomy database: a tarball containing a taxonomy.tab file as well as accession to taxid mapping (a kraken-based taxonomy database will not suffice).",
          patterns: ["*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
        }
        ncbi_taxdump_tgz: {
          description: "An NCBI taxdump.tar.gz file that contains, at the minimum, a nodes.dmp and names.dmp file.",
          patterns: ["*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
        }
    }

    scatter(raw_reads in reads_bams) {
        call reports__fastqc as fastqc_raw {
            input: reads_bam = raw_reads
        }
        call reports__align_and_count as spikein {
            input:
                reads_bam = raw_reads,
                ref_db = spikein_db
        }
    }

    scatter(raw_reads in reads_bams) {
        # separate scatter blocks speeds up the gathers in DNAnexus and provides independent failure blocks
        call metagenomics__kraken2 as kraken2 {
            input:
                reads_bam = raw_reads,
                kraken2_db_tgz = kraken2_db_tgz,
                krona_taxonomy_db_tgz = krona_taxonomy_db_kraken2_tgz
        }
        call metagenomics__filter_bam_to_taxa as deplete {
            input:
                classified_bam = raw_reads,
                classified_reads_txt_gz = kraken2.kraken2_reads_report,
                ncbi_taxonomy_db_tgz = ncbi_taxdump_tgz,
                exclude_taxa = true,
                taxonomic_names = ["Vertebrata"],
                out_filename_suffix = "hs_depleted"
        }
        call reports__fastqc as fastqc_cleaned {
            input: reads_bam = deplete.bam_filtered_to_taxa
        }
        call metagenomics__filter_bam_to_taxa as filter_acellular {
            input:
                classified_bam = raw_reads,
                classified_reads_txt_gz = kraken2.kraken2_reads_report,
                ncbi_taxonomy_db_tgz = ncbi_taxdump_tgz,
                exclude_taxa = true,
                taxonomic_names = ["Vertebrata", "other sequences", "Bacteria"],
                out_filename_suffix = "acellular"
        }
    }

    scatter(clean_reads in filter_acellular.bam_filtered_to_taxa) {
        call read_utils__rmdup_ubam as rmdup_ubam {
           input:
                reads_unmapped_bam = clean_reads
        }
        call assembly__assemble as spades {
            input:
                assembler = "spades",
                reads_unmapped_bam = rmdup_ubam.dedup_bam,
                trim_clip_db = trim_clip_db,
                always_succeed = true
        }
        if(defined(blast_db_tgz) && defined(krona_taxonomy_db_blast_tgz)) {
            call metagenomics__blastx as blastx {
                input:
                    contigs_fasta = spades.contigs_fasta,
                    blast_db_tgz = select_first([blast_db_tgz]),
                    krona_taxonomy_db_tgz = select_first([krona_taxonomy_db_blast_tgz])
            }
        }
    }

    call reports__MultiQC as multiqc_raw {
        input:
            input_files = fastqc_raw.fastqc_zip,
            file_name   = "multiqc-raw.html"
    }

    call reports__MultiQC as multiqc_cleaned {
        input:
            input_files = fastqc_cleaned.fastqc_zip,
            file_name   = "multiqc-cleaned.html"
    }

    call reports__MultiQC as multiqc_dedup {
        input:
            input_files = rmdup_ubam.dedup_fastqc_zip,
            file_name   = "multiqc-dedup.html"
    }

    call reports__align_and_count_summary as spike_summary {
        input:
            counts_txt = spikein.report
    }

    call reports__aggregate_metagenomics_reports as metag_summary_report {
        input:
            kraken_summary_reports = kraken2.kraken2_summary_report
    }

    call metagenomics__krona as krona_merge_kraken2 {
        input:
            reports_txt_gz = kraken2.kraken2_summary_report,
            krona_taxonomy_db_tgz = krona_taxonomy_db_kraken2_tgz,
            input_type = "kraken2",
            out_basename = "merged-kraken2.krona"
    }

    if(defined(blast_db_tgz) && defined(krona_taxonomy_db_blast_tgz)) {
        call metagenomics__krona_merge as krona_merge_blastx {
            input:
                krona_reports = select_all(blastx.krona_report_html),
                out_basename = "merged-spades-blastx.krona"
        }
    }

    output {
        Array[File] cleaned_reads_unaligned_bams = deplete.bam_filtered_to_taxa
        Array[File] deduplicated_reads_unaligned = rmdup_ubam.dedup_bam
        Array[File] contigs_fastas               = spades.contigs_fasta

        Array[Int]  read_counts_raw                 = deplete.classified_taxonomic_filter_read_count_pre
        Array[Int]  read_counts_depleted            = deplete.classified_taxonomic_filter_read_count_post
        Array[Int]  read_counts_dedup               = rmdup_ubam.dedup_read_count_post
        Array[Int]  read_counts_prespades_subsample = spades.subsample_read_count

        File        multiqc_report_raw     = multiqc_raw.multiqc_report
        File        multiqc_report_cleaned = multiqc_cleaned.multiqc_report
        File        multiqc_report_dedup   = multiqc_dedup.multiqc_report
        File        spikein_counts         = spike_summary.count_summary
        File        kraken2_merged_krona   = krona_merge_kraken2.krona_report_html
        File        kraken2_summary        = metag_summary_report.krakenuniq_aggregate_taxlevel_summary
        File?       blastx_merged_krona   = krona_merge_blastx.krona_report_html

        Array[File] kraken2_summary_reports = kraken2.kraken2_summary_report
        Array[File] kraken2_krona_by_sample = kraken2.krona_report_html
        Array[File] blastx_report_by_sample = select_all(blastx.blast_report)
        Array[File] blastx_krona_by_sample  = select_all(blastx.krona_report_html)

        String      kraken2_viral_classify_version = kraken2.viralngs_version[0]
        String      deplete_viral_classify_version    = deplete.viralngs_version[0]
        String      spades_viral_assemble_version     = spades.viralngs_version[0]
    }
}



task reports__fastqc {
  input {
    File     reads_bam

    String   docker="quay.io/broadinstitute/viral-core:2.0.21"
  }

  String   reads_basename=basename(reads_bam, ".bam")

  command {
    set -ex -o pipefail
    reports.py --version | tee VERSION
    reports.py fastqc ${reads_bam} ${reads_basename}_fastqc.html --out_zip ${reads_basename}_fastqc.zip
  }

  output {
    File   fastqc_html      = "${reads_basename}_fastqc.html"
    File   fastqc_zip      = "${reads_basename}_fastqc.zip"
    String viralngs_version = read_string("VERSION")
  }

  runtime {
    memory: "2 GB"
    cpu: 1
    docker: "${docker}"
    disks: "local-disk 375 LOCAL"
    dx_instance_type: "mem1_ssd1_v2_x2"
  }
}




task reports__align_and_count {
  input {
    File    reads_bam
    File    ref_db
    Int?    minScoreToFilter
    Int?    topNHits = 3

    Int?    machine_mem_gb
    String  docker="quay.io/broadinstitute/viral-core:2.0.21"
  }

  String  reads_basename=basename(reads_bam, ".bam")
  String  ref_basename=basename(ref_db, ".fasta")

  command {
    set -ex -o pipefail

    read_utils.py --version | tee VERSION

    ln -s ${reads_bam} ${reads_basename}.bam
    read_utils.py bwamem_idxstats \
      ${reads_basename}.bam \
      ${ref_db} \
      --outStats ${reads_basename}.count.${ref_basename}.txt.unsorted \
      ${'--minScoreToFilter=' + minScoreToFilter} \
      --loglevel=DEBUG

      sort -b -r -n -k3 ${reads_basename}.count.${ref_basename}.txt.unsorted > ${reads_basename}.count.${ref_basename}.txt
      head -n ${topNHits} ${reads_basename}.count.${ref_basename}.txt > ${reads_basename}.count.${ref_basename}.top_${topNHits}_hits.txt
  }

  output {
    File   report           = "${reads_basename}.count.${ref_basename}.txt"
    File   report_top_hits  = "${reads_basename}.count.${ref_basename}.top_${topNHits}_hits.txt"
    String viralngs_version = read_string("VERSION")
  }

  runtime {
    memory: select_first([machine_mem_gb, 7]) + " GB"
    cpu: 4
    docker: "${docker}"
    disks: "local-disk 375 LOCAL"
    dx_instance_type: "mem1_ssd1_v2_x4"
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




task metagenomics__filter_bam_to_taxa {
  input {
    File           classified_bam
    File           classified_reads_txt_gz
    File           ncbi_taxonomy_db_tgz # nodes.dmp names.dmp
    Array[String]? taxonomic_names
    Array[Int]?    taxonomic_ids
    Int?           minimum_hit_groups
    Boolean        withoutChildren=false
    Boolean        exclude_taxa=false
    String         out_filename_suffix = "filtered"

    String         docker="quay.io/broadinstitute/viral-classify:2.0.21.3"
  }

  String out_basename = basename(classified_bam, ".bam") + "." + out_filename_suffix

  command {
    set -ex -o pipefail
    if [ -z "$TMPDIR" ]; then
      export TMPDIR=$(pwd)
    fi

    # decompress taxonomy DB to CWD
    read_utils.py extract_tarball \
      ${ncbi_taxonomy_db_tgz} . \
      --loglevel=DEBUG
    if [ -d "taxonomy" ]; then mv taxonomy/* .; fi

    touch taxfilterargs
    TAXNAMELIST="${write_lines(select_first([taxonomic_names, []]))}"
    if [ -n "$(cat $TAXNAMELIST)" ]; then echo "--taxNames" >> taxfilterargs; fi
    cat $TAXNAMELIST >> taxfilterargs

    TAX_IDs="${sep=' ' taxonomic_ids}"
    if [ -n "$TAX_IDs" ]; then TAX_IDs="--taxIDs $TAX_IDs"; fi

    metagenomics.py --version | tee VERSION

    samtools view -c ${classified_bam} | tee classified_taxonomic_filter_read_count_pre &

    cat taxfilterargs | xargs -d '\n' metagenomics.py filter_bam_to_taxa \
      ${classified_bam} \
      ${classified_reads_txt_gz} \
      "${out_basename}.bam" \
      nodes.dmp \
      names.dmp \
      $TAX_IDs \
      ${true='--exclude' false='' exclude_taxa} \
      ${true='--without-children' false='' withoutChildren} \
      ${'--minimum_hit_groups=' + minimum_hit_groups} \
      --out_count COUNT \
      --loglevel=DEBUG

    samtools view -c "${out_basename}.bam" | tee classified_taxonomic_filter_read_count_post
    wait
  }

  output {
    File    bam_filtered_to_taxa                        = "${out_basename}.bam"
    Int     classified_taxonomic_filter_read_count_pre  = read_int("classified_taxonomic_filter_read_count_pre")
    Int     reads_matching_taxa                         = read_int("COUNT")
    Int     classified_taxonomic_filter_read_count_post = read_int("classified_taxonomic_filter_read_count_post")
    String  viralngs_version                            = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: "7 GB"
    disks: "local-disk 375 LOCAL"
    cpu: 2
    dx_instance_type: "mem1_ssd2_v2_x4"
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




task metagenomics__blastx {
  meta {
    description: "Runs BLASTx classification"
  }

  input {
    File     contigs_fasta
    File     blast_db_tgz
    File     krona_taxonomy_db_tgz

    Int?     machine_mem_gb
    String   docker="quay.io/broadinstitute/viral-classify:2.0.21.3"
  }

  parameter_meta {
    contigs_fasta: {
      description: "Sequences to classify. Use for a small number of longer query sequences (e.g. contigs)",
      patterns: ["*.fasta"] }
    blast_db_tgz: {
      description: "Pre-built BLAST database tarball containing an indexed blast database named 'nr'",
      patterns: ["*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
    }
    krona_taxonomy_db_tgz: {
      description: "Krona taxonomy database: a tarball containing a taxonomy.tab file as well as accession to taxid mapping (a kraken-based taxonomy database will not suffice).",
      patterns: ["*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
    }
  }

  String out_basename=basename(contigs_fasta, '.fasta')

  command {
    set -ex -o pipefail

    if [ -z "$TMPDIR" ]; then
      export TMPDIR=$(pwd)
    fi
    DB_DIR=$(mktemp -d --suffix _db)
    mkdir -p $DB_DIR/blast $DB_DIR/krona

    # decompress DB to $DB_DIR
    read_utils.py extract_tarball \
      ${blast_db_tgz} $DB_DIR/blast \
      --loglevel=DEBUG

    # unpack krona taxonomy database
    read_utils.py extract_tarball \
      ${krona_taxonomy_db_tgz} $DB_DIR/krona \
      --loglevel=DEBUG &  # we don't need this until later

    blastx -version | tee VERSION

    blastx \
      -query ${contigs_fasta} \
      -db $DB_DIR/blast/nr \
      -out "${out_basename}.blastx.contigs.txt" \
      -outfmt 7 \
      -num_threads `nproc`

    wait # for krona_taxonomy_db_tgz to download and extract

    ktImportBLAST \
      -i -k \
      -tax $DB_DIR/krona \
      -o "${out_basename}.blastx.krona.html" \
      "${out_basename}.blastx.contigs.txt","${out_basename}"

    pigz "${out_basename}".blastx.contigs.txt
  }

  output {
    File    blast_report       = "${out_basename}.blastx.contigs.txt.gz"
    File    krona_report_html  = "${out_basename}.blastx.krona.html"
    String  blastx_version     = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: select_first([machine_mem_gb, 8]) + " GB"
    cpu: 32
    disks: "local-disk 375 LOCAL"
    dx_instance_type: "mem1_ssd1_v2_x36"
    preemptible: 1
  }
}




task reports__MultiQC {
  input {
    Array[File]     input_files = []

    Boolean         force = false
    Boolean         full_names = false
    String?         title
    String?         comment
    String?         file_name
    String          out_dir = "./multiqc-output"
    String?         template
    String?         tag
    String?         ignore_analysis_files
    String?         ignore_sample_names
    File?           sample_names
    Array[String]+? exclude_modules
    Array[String]+? module_to_use
    Boolean         data_dir = false
    Boolean         no_data_dir = false
    String?         output_data_format
    Boolean         zip_data_dir = false
    Boolean         export = false
    Boolean         flat = false
    Boolean         interactive = true
    Boolean         lint = false
    Boolean         pdf = false
    Boolean         megaQC_upload = false # Upload generated report to MegaQC if MegaQC options are found
    File?           config  # directory
    String?         config_yaml

    String          docker = "quay.io/biocontainers/multiqc:1.8--py_2"
  }

  parameter_meta {
    output_data_format: { description: "[tsv|yaml|json] default:tsv" }
  }

  # get the basename in all wdl use the filename specified (sans ".html" extension, if specified)
  String report_filename = if (defined(file_name)) then basename(select_first([file_name]), ".html") else "multiqc"

  command {
      set -ex -o pipefail

      echo "${sep='\n' input_files}" > input-filenames.txt
      echo "" >> input-filenames.txt

      multiqc \
      --file-list input-filenames.txt \
      --dirs \
      --outdir "${out_dir}" \
      ${true="--force" false="" force} \
      ${true="--fullnames" false="" full_names} \
      ${"--title " + title} \
      ${"--comment " + comment} \
      ${"--filename " + file_name} \
      ${"--template " + template} \
      ${"--tag " + tag} \
      ${"--ignore " + ignore_analysis_files} \
      ${"--ignore-samples" + ignore_sample_names} \
      ${"--sample-names " + sample_names} \
      ${true="--exclude " false="" defined(exclude_modules)}${sep=" --exclude " exclude_modules} \
      ${true="--module " false="" defined(module_to_use)}${sep=" --module " module_to_use} \
      ${true="--data-dir" false="" data_dir} \
      ${true="--no-data-dir" false="" no_data_dir} \
      ${"--data-format " + output_data_format} \
      ${true="--zip-data-dir" false="" zip_data_dir} \
      ${true="--export" false="" export} \
      ${true="--flat" false="" flat} \
      ${true="--interactive" false="" interactive} \
      ${true="--lint" false="" lint} \
      ${true="--pdf" false="" pdf} \
      ${false="--no-megaqc-upload" true="" megaQC_upload} \
      ${"--config " + config} \
      ${"--cl-config " + config_yaml }

      if [ -z "${file_name}" ]; then
        mv "${out_dir}/${report_filename}_report.html" "${out_dir}/${report_filename}.html"
      fi

      tar -c "${out_dir}/${report_filename}_data" | gzip -c > "${report_filename}_data.tar.gz"
  }

  output {
      File multiqc_report            = "${out_dir}/${report_filename}.html"
      File multiqc_data_dir_tarball  = "${report_filename}_data.tar.gz"
  }

  runtime {
    memory: "2 GB"
    cpu: 1
    docker: "${docker}"
    disks: "local-disk 375 LOCAL"
    dx_instance_type: "mem1_ssd1_v2_x2"
  }
}




task reports__align_and_count_summary {
  input {
    Array[File]+  counts_txt

    String        docker="quay.io/broadinstitute/viral-core:2.0.21"
  }

  command {
    set -ex -o pipefail

    mkdir spike_summaries
    cp ${sep=' ' counts_txt} spike_summaries/

    reports.py --version | tee VERSION
    reports.py aggregate_spike_count spike_summaries/ count_summary.tsv \
      --loglevel=DEBUG
  }

  output {
    File   count_summary    = "count_summary.tsv"
    String viralngs_version = read_string("VERSION")
  }

  runtime {
    memory: "3 GB"
    cpu: 2
    docker: "${docker}"
    disks: "local-disk 50 HDD"
    dx_instance_type: "mem1_ssd1_v2_x2"
  }
}




task reports__aggregate_metagenomics_reports {
  input {
    Array[File]+ kraken_summary_reports 
    String       aggregate_taxon_heading_space_separated  = "Viruses"
    String       aggregate_taxlevel_focus                 = "species"
    Int?         aggregate_top_N_hits                     = 5

    String       docker="quay.io/broadinstitute/viral-classify:2.0.21.3"
  }

  parameter_meta {
    aggregate_taxon_heading_space_separated: { description: "The taxonomic heading to analyze. More than one can be specified." }
    aggregate_taxlevel_focus:                { description: "species,genus,family,order,class,phylum,kingdom,superkingdom" }
    aggregate_top_N_hits:                    { description: "only include the top N hits from a given sample in the aggregate report" }
  }

  String       aggregate_taxon_heading = sub(aggregate_taxon_heading_space_separated, " ", "_") # replace spaces with underscores for use in filename

  command {
    set -ex -o pipefail

    metagenomics.py --version | tee VERSION
    metagenomics.py taxlevel_summary \
      ${sep=' ' kraken_summary_reports} \
      --csvOut aggregate_taxa_summary_${aggregate_taxon_heading}_by_${aggregate_taxlevel_focus}_top_${aggregate_top_N_hits}_by_sample.csv \
      --noHist \
      --taxHeading ${aggregate_taxon_heading_space_separated} \
      --taxlevelFocus ${aggregate_taxlevel_focus} \
      --zeroFill --includeRoot --topN ${aggregate_top_N_hits} \
      --loglevel=DEBUG
  }

  output {
    File   krakenuniq_aggregate_taxlevel_summary = "aggregate_taxa_summary_${aggregate_taxon_heading}_by_${aggregate_taxlevel_focus}_top_${aggregate_top_N_hits}_by_sample.csv"
    String viralngs_version = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: "3 GB"
    cpu: 1
    disks: "local-disk 50 HDD"
    dx_instance_type: "mem1_ssd2_v2_x2"
    preemptible: 0
  }
}




task metagenomics__krona {
  input {
    Array[File]+  reports_txt_gz
    File          krona_taxonomy_db_tgz
    String        out_basename = basename(basename(reports_txt_gz[0], '.gz'), '.txt')

    String?  input_type
    Int?     query_column
    Int?     taxid_column
    Int?     score_column
    Int?     magnitude_column

    Int?     machine_mem_gb
    String   docker="quay.io/broadinstitute/viral-classify:2.0.21.3"
  }

  command {
    set -ex -o pipefail
    if [ -z "$TMPDIR" ]; then
      export TMPDIR=$(pwd)
    fi
    DB_DIR=$(mktemp -d --suffix _db)
    mkdir -p $DB_DIR/krona

    metagenomics.py --version | tee VERSION

    # unpack krona taxonomy.tab
    if [[ ${krona_taxonomy_db_tgz} == *.tar.* ]]; then
      read_utils.py extract_tarball \
        ${krona_taxonomy_db_tgz} $DB_DIR/krona \
        --loglevel=DEBUG
    else
      if [[ "${krona_taxonomy_db_tgz}" == *.zst ]]; then
        cat "${krona_taxonomy_db_tgz}" | zstd -d > $DB_DIR/krona/taxonomy.tab
      elif [[ "${krona_taxonomy_db_tgz}" == *.gz ]]; then
        cat "${krona_taxonomy_db_tgz}" | pigz -dc > $DB_DIR/krona/taxonomy.tab
      elif [[ "${krona_taxonomy_db_tgz}" == *.bz2 ]]; then
        cat "${krona_taxonomy_db_tgz}" | bzip -dc > $DB_DIR/krona/taxonomy.tab
      else
        cp "${krona_taxonomy_db_tgz}" $DB_DIR/krona/taxonomy.tab
      fi
    fi

    metagenomics.py krona \
      ${sep=' ' reports_txt_gz} \
      $DB_DIR/krona \
      ${out_basename}.html \
      ${'--inputType=' + input_type} \
      ${'--queryColumn=' + query_column} \
      ${'--taxidColumn=' + taxid_column} \
      ${'--scoreColumn=' + score_column} \
      ${'--magnitudeColumn=' + magnitude_column} \
      --noRank --noHits \
      --loglevel=DEBUG
  }

  output {
    File    krona_report_html  = "${out_basename}.html"
    String  viralngs_version   = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: "3 GB"
    cpu: 1
    disks: "local-disk 50 HDD"
    dx_instance_type: "mem1_ssd2_v2_x2"
  }
}




task metagenomics__krona_merge {
  input {
    Array[File]  krona_reports
    String       out_basename

    Int?         machine_mem_gb
    String       docker="biocontainers/krona:v2.7.1_cv1"
  }

  command {
    set -ex -o pipefail
    ktImportKrona | head -2 | tail -1 | cut -f 2-3 -d ' ' | tee VERSION
    ktImportKrona -o "${out_basename}.html" ${sep=' ' krona_reports}
  }

  output {
    File    krona_report_html = "${out_basename}.html"
    String  krona_version     = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: select_first([machine_mem_gb, 3]) + " GB"
    cpu: 1
    disks: "local-disk 50 HDD"
    dx_instance_type: "mem1_ssd2_v2_x2"
  }
}


