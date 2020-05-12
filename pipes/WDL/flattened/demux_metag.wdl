version 1.0








workflow demux_metag {
    input {
        File spikein_db
        File trim_clip_db
        Array[File]? bmtaggerDbs  # .tar.gz, .tgz, .tar.bz2, .tar.lz4, .fasta, or .fasta.gz
        Array[File]? blastDbs  # .tar.gz, .tgz, .tar.bz2, .tar.lz4, .fasta, or .fasta.gz
        Array[File]? bwaDbs

        File kraken2_db_tgz
        File krona_taxonomy_db_kraken2_tgz
        File blast_db_tgz
        File krona_taxonomy_db_blast_tgz
    }

    call demux__illumina_demux as illumina_demux

    scatter(raw_reads in illumina_demux.raw_reads_unaligned_bams) {
        call reports__align_and_count as spikein {
            input:
                reads_bam = raw_reads,
                ref_db = spikein_db
        }
        call taxon_filter__deplete_taxa as deplete {
            input:
                raw_reads_unmapped_bam = raw_reads,
                bmtaggerDbs = bmtaggerDbs,
                blastDbs = blastDbs,
                bwaDbs = bwaDbs
        }
        call read_utils__rmdup_ubam as rmdup_ubam {
           input:
                reads_unmapped_bam = deplete.cleaned_bam
        }
        call assembly__assemble as spades {
            input:
                assembler = "spades",
                reads_unmapped_bam = rmdup_ubam.dedup_bam,
                trim_clip_db = trim_clip_db,
                always_succeed = true
        }
        call metagenomics__kraken2 as kraken2 {
            input:
                reads_bam = raw_reads,
                kraken2_db_tgz = kraken2_db_tgz,
                krona_taxonomy_db_tgz = krona_taxonomy_db_kraken2_tgz
        }
        call metagenomics__blastx as blastx {
            input:
                contigs_fasta = spades.contigs_fasta,
                blast_db_tgz = blast_db_tgz,
                krona_taxonomy_db_tgz = krona_taxonomy_db_blast_tgz
        }
    }

    call reports__MultiQC as multiqc_raw {
        input:
            input_files = illumina_demux.raw_reads_fastqc_zip,
            file_name   = "multiqc-raw.html"
    }

    call reports__MultiQC as multiqc_cleaned {
        input:
            input_files = deplete.cleaned_fastqc_zip,
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

    call metagenomics__krona_merge as krona_merge_kraken2 {
        input:
            krona_reports = kraken2.krona_report_html,
            out_basename = "merged-kraken2.krona.html"
    }

    call metagenomics__krona_merge as krona_merge_blastx {
        input:
            krona_reports = blastx.krona_report_html,
            out_basename = "merged-spades-blastx.krona.html"
    }

    output {
        Array[File] raw_reads_unaligned_bams     = illumina_demux.raw_reads_unaligned_bams
        Array[File] cleaned_reads_unaligned_bams = deplete.cleaned_bam
        Array[File] deduplicated_reads_unaligned = rmdup_ubam.dedup_bam
        Array[File] contigs_fastas               = spades.contigs_fasta

        Array[Int]  read_counts_raw                 = deplete.depletion_read_count_pre
        Array[Int]  read_counts_depleted            = deplete.depletion_read_count_post
        Array[Int]  read_counts_dedup               = rmdup_ubam.dedup_read_count_post
        Array[Int]  read_counts_prespades_subsample = spades.subsample_read_count

        File        demux_metrics            = illumina_demux.metrics
        File        demux_commonBarcodes     = illumina_demux.commonBarcodes
        File        demux_outlierBarcodes    = illumina_demux.outlierBarcodes

        File        multiqc_report_raw     = multiqc_raw.multiqc_report
        File        multiqc_report_cleaned = multiqc_cleaned.multiqc_report
        File        multiqc_report_dedup   = multiqc_dedup.multiqc_report
        File        spikein_counts         = spike_summary.count_summary
        File        kraken2_merged_krona   = krona_merge_kraken2.krona_report_html
        File        kraken2_summary        = metag_summary_report.krakenuniq_aggregate_taxlevel_summary
        File        blastx_merged_krona   = krona_merge_blastx.krona_report_html

        Array[File] kraken2_summary_reports = kraken2.kraken2_summary_report
        Array[File] kraken2_krona_by_sample = kraken2.krona_report_html
        Array[File] blastx_report_by_sample = blastx.blast_report
        Array[File] blastx_krona_by_sample  = blastx.krona_report_html

        String      demux_viral_core_version          = illumina_demux.viralngs_version
        String      kraken2_viral_classify_version = kraken2.viralngs_version[0]
        String      deplete_viral_classify_version    = deplete.viralngs_version[0]
        String      spades_viral_assemble_version     = spades.viralngs_version[0]
    }
}



task demux__illumina_demux {
  input {
    File    flowcell_tgz
    Int?    lane=1
    File?   samplesheet
    File?   runinfo
    String? sequencingCenter

    String? flowcell
    Int?    minimumBaseQuality = 10
    Int?    maxMismatches = 0
    Int?    minMismatchDelta
    Int?    maxNoCalls
    String? readStructure
    Int?    minimumQuality
    Int?    threads
    String? runStartDate
    Int?    maxReadsInRamPerTile
    Int?    maxRecordsInRam
    Boolean? forceGC=true

    Int?    machine_mem_gb
    String  docker="quay.io/broadinstitute/viral-core:2.0.21"
  }

  command {
    set -ex -o pipefail

    # find N% memory
    mem_in_mb=`/opt/viral-ngs/source/docker/calc_mem.py mb 85`

    if [ -z "$TMPDIR" ]; then
      export TMPDIR=$(pwd)
    fi
    FLOWCELL_DIR=$(mktemp -d)

    read_utils.py --version | tee VERSION

    read_utils.py extract_tarball \
      ${flowcell_tgz} $FLOWCELL_DIR \
      --loglevel=DEBUG

    # if we are overriding the RunInfo file, use the path of the file provided. Otherwise find the file
    if [ -n "${runinfo}" ]; then
      RUNINFO_FILE="${runinfo}"
    else
      # full RunInfo.xml path
      RUNINFO_FILE="$(find $FLOWCELL_DIR -type f -maxdepth 3 -name RunInfo.xml | head -n 1)"
    fi
    
    # Parse the lane count & run ID from RunInfo.xml file
    lane_count=$(xmllint --xpath "string(//Run/FlowcellLayout/@LaneCount)" $RUNINFO_FILE)
    if [ -z "$lane_count" ]; then
        echo "Could not parse LaneCount from RunInfo.xml. Please check RunInfo.xml is properly formatted"
    fi

    surface_count=$(xmllint --xpath "string(//Run/FlowcellLayout/@SurfaceCount)" $RUNINFO_FILE)
    if [ -z "$surface_count" ]; then
        echo "Could not parse SurfaceCount from RunInfo.xml. Please check RunInfo.xml is properly formatted"
    fi

    swath_count=$(xmllint --xpath "string(//Run/FlowcellLayout/@SwathCount)" $RUNINFO_FILE)
    if [ -z "$swath_count" ]; then
        echo "Could not parse SwathCount from RunInfo.xml. Please check RunInfo.xml is properly formatted"
    fi

    tile_count=$(xmllint --xpath "string(//Run/FlowcellLayout/@TileCount)" $RUNINFO_FILE)
    if [ -z "$tile_count" ]; then
        echo "Could not parse TileCount from RunInfo.xml. Please check RunInfo.xml is properly formatted"
    fi

    # total data size more roughly tracks total tile count
    total_tile_count=$((lane_count*surface_count*swath_count*tile_count))

    demux_threads="$(nproc --all)"
    if [ "$total_tile_count" -le 50 ]; then
        echo "Detected $total_tile_count tiles, interpreting as MiSeq run."
    elif [ "$total_tile_count" -le 150 ]; then
        echo "Detected $total_tile_count tiles, interpreting as HiSeq2k run."
    elif [ "$total_tile_count" -le 288 ]; then
        # increase the number of reads in ram per-tile for NextSeq, since the tiles are larger
        # without this setting, reads will spill to disk and may read the limit
        # on the number of files that can be opened
        max_reads_in_ram_per_tile=1500000
        max_records_in_ram=2000000
        echo "Detected $total_tile_count tiles, interpreting as NextSeq (mid-output) run."
    elif [ "$total_tile_count" -le 624 ]; then
        demux_threads=32 # with NovaSeq-size output, OOM errors can sporadically occur with higher thread counts
        mem_in_mb=$(/opt/viral-ngs/source/docker/calc_mem.py mb 80)
        max_reads_in_ram_per_tile=600000 # reduce the number of reads per tile since the NovaSeq has so many
        max_records_in_ram=2000000
        echo "Detected $total_tile_count tiles, interpreting as NovaSeq SP run."
    elif [ "$total_tile_count" -le 864 ]; then
        # increase the number of reads in ram per-tile for NextSeq, since the tiles are larger
        # without this setting, reads will spill to disk and may read the limit
        # on the number of files that can be opened
        max_reads_in_ram_per_tile=1500000 # reduce the number of reads per tile since the NovaSeq has so many
        max_records_in_ram=2500000
        echo "Detected $total_tile_count tiles, interpreting as NextSeq (high-output) run."
    elif [ "$total_tile_count" -le 896 ]; then
        echo "Detected $total_tile_count tiles, interpreting as HiSeq4k run."
    elif [ "$total_tile_count" -le 1408 ]; then
        demux_threads=32 # with NovaSeq-size output, OOM errors can sporadically occur with higher thread counts
        mem_in_mb=$(/opt/viral-ngs/source/docker/calc_mem.py mb 80)
        max_reads_in_ram_per_tile=600000 # reduce the number of reads per tile since the NovaSeq has so many
        max_records_in_ram=2000000
        echo "Detected $total_tile_count tiles, interpreting as NovaSeq run."
        echo "  **Note: Q20 threshold used since NovaSeq with RTA3 writes only four Q-score values: 2, 12, 23, and 37.**"
        echo "    See: https://www.illumina.com/content/dam/illumina-marketing/documents/products/appnotes/novaseq-hiseq-q30-app-note-770-2017-010.pdf"
    elif [ "$total_tile_count" -gt 1408 ]; then
        demux_threads=30 # with NovaSeq-size output, OOM errors can sporadically occur with higher thread counts
        mem_in_mb=$(/opt/viral-ngs/source/docker/calc_mem.py mb 80)
        max_reads_in_ram_per_tile=600000 # reduce the number of reads per tile since the NovaSeq has so many
        max_records_in_ram=2000000
        echo "Tile count: $total_tile_count tiles (unknown instrument type)."
    fi

    # use the passed-in (or default) WDL value first, then fall back to the auto-scaled value
    # if the result of this is null (nothing is passed in, no autoscaled value, no param is passed to the command)
    if [ -n "${minimumBaseQuality}" ]; then demux_min_base_quality="${minimumBaseQuality}"; else demux_min_base_quality="$demux_min_base_quality"; fi
    if [ -n "$demux_min_base_quality" ]; then demux_min_base_quality="--minimum_base_quality=$demux_min_base_quality";fi
    
    if [ -n "${threads}" ]; then demux_threads="${threads}"; else demux_threads="$demux_threads"; fi
    if [ -n "$demux_threads" ]; then demux_threads="--threads=$demux_threads"; fi
    

    if [ -n "${maxReadsInRamPerTile}" ]; then max_reads_in_ram_per_tile="${maxReadsInRamPerTile}"; else max_reads_in_ram_per_tile="$max_reads_in_ram_per_tile"; fi
    if [ -n "$max_reads_in_ram_per_tile" ]; then max_reads_in_ram_per_tile="--max_reads_in_ram_per_tile=$max_reads_in_ram_per_tile"; fi
    
    if [ -n "${maxRecordsInRam}" ]; then max_records_in_ram="${maxRecordsInRam}"; else max_records_in_ram="$max_records_in_ram"; fi
    if [ -n "$max_records_in_ram" ]; then max_records_in_ram="--max_records_in_ram=$max_records_in_ram"; fi

    # note that we are intentionally setting --threads to about 2x the core
    # count. seems to still provide speed benefit (over 1x) when doing so.
    illumina.py illumina_demux \
      $FLOWCELL_DIR \
      ${lane} \
      . \
      ${'--sampleSheet=' + samplesheet} \
      ${'--runInfo=' + runinfo} \
      ${'--sequencing_center=' + sequencingCenter} \
      --outMetrics=metrics.txt \
      --commonBarcodes=barcodes.txt \
      ${'--flowcell=' + flowcell} \
      $demux_min_base_quality \
      ${'--max_mismatches=' + maxMismatches} \
      ${'--min_mismatch_delta=' + minMismatchDelta} \
      ${'--max_no_calls=' + maxNoCalls} \
      ${'--read_structure=' + readStructure} \
      ${'--minimum_quality=' + minimumQuality} \
      ${'--run_start_date=' + runStartDate} \
      $max_reads_in_ram_per_tile \
      $max_records_in_ram \
      --JVMmemory="$mem_in_mb"m \
      $demux_threads \
      ${true='--force_gc=true' false="--force_gc=false" forceGC} \
      --compression_level=5 \
      --loglevel=DEBUG

    illumina.py guess_barcodes --expected_assigned_fraction=0 barcodes.txt metrics.txt barcodes_outliers.txt

    mkdir -p unmatched
    mv Unmatched.bam unmatched/

    OUT_BASENAMES=bam_basenames.txt
    for bam in *.bam; do
      echo "$(basename $bam .bam)" >> $OUT_BASENAMES
    done

    FASTQC_HARDCODED_MEM_PER_THREAD=250 # the value fastqc sets for -Xmx per thread, not adjustable
    num_cpus=$(nproc)
    num_bam_files=$(cat $OUT_BASENAMES | wc -l)
    num_fastqc_jobs=1
    num_fastqc_threads=1
    total_ram_needed_mb=250

    # determine the number of fastq jobs
    while [[ $total_ram_needed_mb -lt $mem_in_mb ]] && [[ $num_fastqc_jobs -lt $num_cpus ]] && [[ $num_fastqc_jobs -lt $num_bam_files ]]; do
        num_fastqc_jobs=$(($num_fastqc_jobs+1))
        total_ram_needed_mb=$(($total_ram_needed_mb+$FASTQC_HARDCODED_MEM_PER_THREAD))
    done
    # determine the number of fastqc threads per job
    while [[ $(($total_ram_needed_mb)) -lt $mem_in_mb ]] && [[ $(($num_fastqc_jobs*$num_fastqc_threads)) -lt $num_cpus ]]; do
        if [[ $(( $num_fastqc_jobs * $(($num_fastqc_threads+1)) )) -le $num_cpus ]]; then
            num_fastqc_threads=$(($num_fastqc_threads+1))
            total_ram_needed_mb=$(($num_fastqc_jobs*($FASTQC_HARDCODED_MEM_PER_THREAD*$num_fastqc_threads)))
        else
            break
        fi
    done

    # GNU Parallel refresher:
    # ",," is the replacement string; values after ":::" are substituted where it appears
    parallel --jobs $num_fastqc_jobs -I ,, \
      "reports.py fastqc \
        ,,.bam \
        ,,_fastqc.html \
        --out_zip ,,_fastqc.zip \
        --threads $num_fastqc_threads" \
      ::: `cat $OUT_BASENAMES`
  }

  output {
    File        metrics                  = "metrics.txt"
    File        commonBarcodes           = "barcodes.txt"
    File        outlierBarcodes          = "barcodes_outliers.txt"
    Array[File] raw_reads_unaligned_bams = glob("*.bam")
    File        unmatched_reads_bam      = "unmatched/Unmatched.bam"
    Array[File] raw_reads_fastqc         = glob("*_fastqc.html")
    Array[File] raw_reads_fastqc_zip     = glob("*_fastqc.zip")
    String      viralngs_version         = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: select_first([machine_mem_gb, 200]) + " GB"
    cpu: 32
    disks: "local-disk 2625 LOCAL"
    dx_instance_type: "mem3_ssd2_v2_x32"
    dx_timeout: "20H"
    preemptible: 0  # this is the very first operation before scatter, so let's get it done quickly & reliably
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


