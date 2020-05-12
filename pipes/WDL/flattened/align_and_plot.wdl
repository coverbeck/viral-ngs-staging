version 1.0




workflow align_and_plot {
    input {
        String? aligner_options = "-r Random -l 30 -g 40 -x 20 -t 502"
    }

    call assembly__align_reads as align {
        input:
            aligner_options = aligner_options
    }
    call reports__plot_coverage as plot_coverage {
        input:
            aligned_reads_bam = align.aligned_only_reads_bam,
            sample_name = basename(basename(basename(align.aligned_only_reads_bam, ".bam"), ".mapped"), ".clean")
    }

    output {
        File   aligned_bam                   = align.aligned_bam
        File   aligned_bam_idx               = align.aligned_bam_idx
        File   aligned_bam_flagstat          = align.aligned_bam_flagstat
        File   aligned_only_reads_bam        = align.aligned_only_reads_bam
        File   aligned_only_reads_bam_idx    = align.aligned_only_reads_bam_idx
        File   aligned_only_reads_fastqc     = align.aligned_only_reads_fastqc
        File   aligned_only_reads_fastqc_zip = align.aligned_only_reads_fastqc_zip
        Int    reads_provided                = align.reads_provided
        Int    reads_aligned                 = align.reads_aligned
        Int    read_pairs_aligned            = align.read_pairs_aligned
        Int    bases_aligned                 = align.bases_aligned
        Float  mean_coverage                 = align.mean_coverage
        String align_viral_core_version      = align.viralngs_version
        File   coverage_plot                 = plot_coverage.coverage_plot
        File   coverage_tsv                  = plot_coverage.coverage_tsv
        Int    reference_length              = plot_coverage.assembly_length
        String plot_viral_core_version       = plot_coverage.viralngs_version
    }
}



task assembly__align_reads {
  meta {
    description: "Align unmapped reads to a reference genome, either using novoalign (default) or bwa. Produces an aligned bam file (including all unmapped reads), an aligned-only bam file, both sorted and indexed, along with samtools flagstat output, fastqc stats (on mapped only reads), and some basic figures of merit."
  }

  input {
    File     reference_fasta
    File     reads_unmapped_bam

    File?    novocraft_license

    String?  aligner="novoalign" # novoalign or bwa
    String?  aligner_options
    Boolean? skip_mark_dupes=false

    String   docker="quay.io/broadinstitute/viral-core:2.0.21"

    String   sample_name = basename(basename(basename(reads_unmapped_bam, ".bam"), ".taxfilt"), ".clean")
  }

  parameter_meta {
    aligner: { description: "Short read aligner to use: novoalign or bwa. (Default: novoalign)" }
  }
  
  command {
    set -ex -o pipefail

    read_utils.py --version | tee VERSION

    cp ${reference_fasta} assembly.fasta
    grep -v '^>' assembly.fasta | tr -d '\n' | wc -c | tee assembly_length

    if [ "$(cat assembly_length)" != "0" ]; then

      # only perform the following if the reference is non-empty

      if [ "${aligner}" == "novoalign" ]; then
        read_utils.py novoindex \
          assembly.fasta \
          ${"--NOVOALIGN_LICENSE_PATH=" + novocraft_license} \
          --loglevel=DEBUG
      fi
      read_utils.py index_fasta_picard assembly.fasta --loglevel=DEBUG
      read_utils.py index_fasta_samtools assembly.fasta --loglevel=DEBUG

      read_utils.py align_and_fix \
        ${reads_unmapped_bam} \
        assembly.fasta \
        --outBamAll "${sample_name}.all.bam" \
        --outBamFiltered "${sample_name}.mapped.bam" \
        --aligner ${aligner} \
        ${'--aligner_options "' + aligner_options + '"'} \
        ${true='--skipMarkDupes' false="" skip_mark_dupes} \
        --JVMmemory=3g \
        ${"--NOVOALIGN_LICENSE_PATH=" + novocraft_license} \
        --loglevel=DEBUG

    else
      touch "${sample_name}.all.bam" "${sample_name}.mapped.bam"

    fi
    samtools index ${sample_name}.mapped.bam

    # collect figures of merit
    grep -v '^>' assembly.fasta | tr -d '\nNn' | wc -c | tee assembly_length_unambiguous
    samtools view -c ${reads_unmapped_bam} | tee reads_provided
    samtools view -c ${sample_name}.mapped.bam | tee reads_aligned
    # report only primary alignments 260=exclude unaligned reads and secondary mappings
    samtools view -h -F 260 ${sample_name}.all.bam | samtools flagstat - | tee ${sample_name}.all.bam.flagstat.txt
    grep properly ${sample_name}.all.bam.flagstat.txt | cut -f 1 -d ' ' | tee read_pairs_aligned
    samtools view ${sample_name}.mapped.bam | cut -f10 | tr -d '\n' | wc -c | tee bases_aligned
    python -c "print (float("`cat bases_aligned`")/"`cat assembly_length_unambiguous`") if "`cat assembly_length_unambiguous`">0 else 0" > mean_coverage

    # fastqc mapped bam
    reports.py fastqc ${sample_name}.mapped.bam ${sample_name}.mapped_fastqc.html --out_zip ${sample_name}.mapped_fastqc.zip
  }

  output {
    File   aligned_bam                   = "${sample_name}.all.bam"
    File   aligned_bam_idx               = "${sample_name}.all.bai"
    File   aligned_bam_flagstat          = "${sample_name}.all.bam.flagstat.txt"
    File   aligned_only_reads_bam        = "${sample_name}.mapped.bam"
    File   aligned_only_reads_bam_idx    = "${sample_name}.mapped.bai"
    File   aligned_only_reads_fastqc     = "${sample_name}.mapped_fastqc.html"
    File   aligned_only_reads_fastqc_zip = "${sample_name}.mapped_fastqc.zip"
    Int    reads_provided                = read_int("reads_provided")
    Int    reads_aligned                 = read_int("reads_aligned")
    Int    read_pairs_aligned            = read_int("read_pairs_aligned")
    Int    bases_aligned                 = read_int("bases_aligned")
    Float  mean_coverage                 = read_float("mean_coverage")
    String viralngs_version              = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: "7 GB"
    cpu: 8
    disks: "local-disk 375 LOCAL"
    dx_instance_type: "mem1_ssd1_v2_x8"
    preemptible: 1
  }
}




task reports__plot_coverage {
  input {
    File     aligned_reads_bam
    String   sample_name

    Boolean? skip_mark_dupes=false
    Boolean? plot_only_non_duplicates=false
    Boolean? bin_large_plots=false
    String?  binning_summary_statistic="max" # max or min

    String   docker="quay.io/broadinstitute/viral-core:2.0.21"
  }
  
  command {
    set -ex -o pipefail

    read_utils.py --version | tee VERSION

    samtools view -c ${aligned_reads_bam} | tee reads_aligned
    if [ "$(cat reads_aligned)" != "0" ]; then
      samtools index -@ "$(nproc)" "${aligned_reads_bam}"

      PLOT_DUPE_OPTION=""
      if [[ "${skip_mark_dupes}" != "true" ]]; then
        PLOT_DUPE_OPTION="${true='--plotOnlyNonDuplicates' false="" plot_only_non_duplicates}"
      fi
      
      BINNING_OPTION="${true='--binLargePlots' false="" bin_large_plots}"

      # plot coverage
      reports.py plot_coverage \
        "${aligned_reads_bam}" \
        "${sample_name}.coverage_plot.pdf" \
        --outSummary "${sample_name}.coverage_plot.txt" \
        --plotFormat pdf \
        --plotWidth 1100 \
        --plotHeight 850 \
        --plotDPI 100 \
        $PLOT_DUPE_OPTION \
        $BINNING_OPTION \
        --binningSummaryStatistic ${binning_summary_statistic} \
        --plotTitle "${sample_name} coverage plot" \
        --loglevel=DEBUG

    else
      touch ${sample_name}.coverage_plot.pdf ${sample_name}.coverage_plot.txt
    fi

    # collect figures of merit
    samtools view -H ${aligned_reads_bam} | perl -n -e'/^@SQ.*LN:(\d+)/ && print "$1\n"' |  python -c "import sys; print(sum(int(x) for x in sys.stdin))" | tee assembly_length
    # report only primary alignments 260=exclude unaligned reads and secondary mappings
    samtools view -h -F 260 ${aligned_reads_bam} | samtools flagstat - | tee ${sample_name}.flagstat.txt
    grep properly ${sample_name}.flagstat.txt | cut -f 1 -d ' ' | tee read_pairs_aligned
    samtools view ${aligned_reads_bam} | cut -f10 | tr -d '\n' | wc -c | tee bases_aligned
    python -c "print (float("`cat bases_aligned`")/"`cat assembly_length`") if "`cat assembly_length`">0 else 0" > mean_coverage
  }

  output {
    File   coverage_plot                 = "${sample_name}.coverage_plot.pdf"
    File   coverage_tsv                  = "${sample_name}.coverage_plot.txt"
    Int    assembly_length               = read_int("assembly_length")
    Int    reads_aligned                 = read_int("reads_aligned")
    Int    read_pairs_aligned            = read_int("read_pairs_aligned")
    Int    bases_aligned                 = read_int("bases_aligned")
    Float  mean_coverage                 = read_float("mean_coverage")
    String viralngs_version              = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: "7 GB"
    cpu: 2
    disks: "local-disk 375 LOCAL"
    dx_instance_type: "mem1_ssd1_v2_x4"
    preemptible: 1
  }
}


