version 1.0




workflow demux_only {
    call tasks_demux__illumina_demux as illumina_demux

    call reports__MultiQC as MultiQC {
        input:
            input_files = illumina_demux.raw_reads_fastqc_zip
    }

    output {
        Array[File] raw_reads_unaligned_bams = illumina_demux.raw_reads_unaligned_bams
        File        demux_metrics            = illumina_demux.metrics
        File        demux_commonBarcodes     = illumina_demux.commonBarcodes
        File        demux_outlierBarcodes    = illumina_demux.outlierBarcodes
        File        multiqc_report_raw       = MultiQC.multiqc_report
        String      demux_viral_core_version = illumina_demux.viralngs_version
    }
}



task tasks_demux__illumina_demux {
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


