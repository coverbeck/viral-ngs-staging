version 1.0






workflow assemble_denovo {
  
  input {
    File         reads_unmapped_bam

    Array[File]+ reference_genome_fasta

    Array[File]  deplete_bmtaggerDbs = []
    Array[File]  deplete_blastDbs = []
    Array[File]  deplete_bwaDbs =[]

    File?        filter_to_taxon_db
    File         trim_clip_db

    File?        novocraft_license

    Boolean      call_isnvs=false

    String       assembler="spades"
    Float?       scaffold_min_length_fraction
    Float?       scaffold_min_unambig
    Int?         scaffold_replace_length=55
    Int?         nucmer_max_gap
    Int?         nucmer_min_match
    Int?         nucmer_min_cluster
    Float?       scaffold_min_pct_contig_aligned
  }

  parameter_meta {
    raw_reads_unmapped_bam: { description: "unaligned reads in BAM format", patterns: ["*.bam"] }
    deplete_bmtaggerDbs: {
       description: "Optional list of databases to use for bmtagger-based depletion. Sequences in fasta format will be indexed on the fly, pre-bmtagger-indexed databases may be provided as tarballs.",
       patterns: ["*.fasta", "*.fasta.gz", "*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
    }
    deplete_blastDbs: {
      description: "Optional list of databases to use for blastn-based depletion. Sequences in fasta format will be indexed on the fly, pre-blast-indexed databases may be provided as tarballs.",
      patterns: ["*.fasta", "*.fasta.gz", "*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
    }
    deplete_bwaDbs: {
      description: "Optional list of databases to use for bwa mem-based depletion. Sequences in fasta format will be indexed on the fly, pre-bwa-indexed databases may be provided as tarballs.",
      patterns: ["*.fasta", "*.fasta.gz", "*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
    }
    filter_to_taxon_db: {
      description: "Optional database to use to filter read set to those that match by LASTAL. Sequences in fasta format will be indexed on the fly.",
      patterns: ["*.fasta"]
    }
    reference_genome_fasta: {
      description: "After denovo assembly, large contigs are scaffolded against a reference genome to determine orientation and to join contigs together, before further polishing by reads. You must supply at least one reference genome (all segments/chromomes in a single fasta file). If more than one reference is provided, contigs will be scaffolded against all of them and the one with the most complete assembly will be chosen for downstream polishing.",
      patterns: ["*.fasta"]
    }
  }

  String sample_name = basename(basename(reads_unmapped_bam, ".bam"), ".cleaned")

  if(length(deplete_bmtaggerDbs) + length(deplete_blastDbs) + length(deplete_bwaDbs) > 0) {
    call taxon_filter__deplete_taxa as deplete_taxa {
      input:
        raw_reads_unmapped_bam = reads_unmapped_bam,
        bmtaggerDbs = deplete_bmtaggerDbs,
        blastDbs = deplete_blastDbs,
        bwaDbs = deplete_bwaDbs
    }
  }

  if(defined(filter_to_taxon_db)) {
    call taxon_filter__filter_to_taxon as filter_to_taxon {
      input:
        reads_unmapped_bam = select_first([deplete_taxa.cleaned_bam, reads_unmapped_bam]),
        lastal_db_fasta = select_first([filter_to_taxon_db])
    }
  }

  call read_utils__rmdup_ubam as rmdup_ubam {
    input:
      reads_unmapped_bam = select_first([filter_to_taxon.taxfilt_bam, deplete_taxa.cleaned_bam, reads_unmapped_bam])
  }

  call assembly__assemble as assemble {
    input:
      reads_unmapped_bam = rmdup_ubam.dedup_bam,
      trim_clip_db = trim_clip_db,
      always_succeed = true,
      assembler = assembler,
      sample_name = sample_name
  }

  call assembly__scaffold as scaffold {
    input:
      contigs_fasta = assemble.contigs_fasta,
      reads_bam = select_first([filter_to_taxon.taxfilt_bam, deplete_taxa.cleaned_bam, reads_unmapped_bam]),
      reference_genome_fasta = reference_genome_fasta,
      min_length_fraction = scaffold_min_length_fraction,
      min_unambig = scaffold_min_unambig,
      replace_length = scaffold_replace_length,
      nucmer_max_gap = nucmer_max_gap,
      nucmer_min_match = nucmer_min_match,
      nucmer_min_cluster = nucmer_min_cluster,
      scaffold_min_pct_contig_aligned = scaffold_min_pct_contig_aligned
  }

  call assembly__refine_2x_and_plot as refine_2x_and_plot {
    input:
      assembly_fasta = scaffold.scaffold_fasta,
      reads_unmapped_bam = select_first([deplete_taxa.cleaned_bam, reads_unmapped_bam]),
      novocraft_license = novocraft_license,
      sample_name = sample_name
  }

  if(call_isnvs) {
    call intrahost__isnvs_per_sample as isnvs_per_sample {
        input:
            assembly_fasta = refine_2x_and_plot.final_assembly_fasta,
            mapped_bam = refine_2x_and_plot.aligned_bam
    }
  }

  output {
    File  final_assembly_fasta        = refine_2x_and_plot.final_assembly_fasta
    File  aligned_only_reads_bam      = refine_2x_and_plot.aligned_only_reads_bam
    File  coverage_plot               = refine_2x_and_plot.coverage_plot
    Int   assembly_length             = refine_2x_and_plot.assembly_length
    Int   assembly_length_unambiguous = refine_2x_and_plot.assembly_length_unambiguous
    Int   reads_aligned               = refine_2x_and_plot.reads_aligned
    Float mean_coverage               = refine_2x_and_plot.mean_coverage

    File?  cleaned_bam               = select_first([deplete_taxa.cleaned_bam, reads_unmapped_bam])
    File?  cleaned_fastqc            = deplete_taxa.cleaned_fastqc
    Int?   depletion_read_count_pre  = deplete_taxa.depletion_read_count_pre
    Int?   depletion_read_count_post = deplete_taxa.depletion_read_count_post

    File?  taxfilt_bam               = filter_to_taxon.taxfilt_bam
    File?  taxfilt_fastqc            = filter_to_taxon.taxfilt_fastqc
    Int?   filter_read_count_post    = filter_to_taxon.filter_read_count_post

    File   dedup_bam                 = rmdup_ubam.dedup_bam
    File   dedup_fastqc              = rmdup_ubam.dedup_fastqc
    Int    dedup_read_count_post     = rmdup_ubam.dedup_read_count_post

    File   contigs_fasta             = assemble.contigs_fasta
    File   subsampBam                = assemble.subsampBam
    Int    subsample_read_count      = assemble.subsample_read_count

    File   scaffold_fasta                        = scaffold.scaffold_fasta
    File   intermediate_scaffold_fasta           = scaffold.intermediate_scaffold_fasta
    File   intermediate_gapfill_fasta            = scaffold.intermediate_gapfill_fasta
    Int    assembly_preimpute_length             = scaffold.assembly_preimpute_length
    Int    assembly_preimpute_length_unambiguous = scaffold.assembly_preimpute_length_unambiguous
    String scaffolding_chosen_ref_name           = scaffold.scaffolding_chosen_ref_name
    File   scaffolding_stats                     = scaffold.scaffolding_stats
    File   scaffolding_alt_contigs               = scaffold.scaffolding_alt_contigs

    File?  isnvsFile                   = isnvs_per_sample.isnvsFile

    File aligned_bam                   = refine_2x_and_plot.aligned_bam
    File aligned_only_reads_bam_idx    = refine_2x_and_plot.aligned_only_reads_bam_idx
    File aligned_only_reads_fastqc     = refine_2x_and_plot.aligned_only_reads_fastqc
    File coverage_tsv                  = refine_2x_and_plot.coverage_tsv
    Int  read_pairs_aligned            = refine_2x_and_plot.read_pairs_aligned
    Int  bases_aligned                 = refine_2x_and_plot.bases_aligned

    String? deplete_viral_classify_version  = deplete_taxa.viralngs_version
    String? taxfilt_viral_classify_version  = filter_to_taxon.viralngs_version
    String  assemble_viral_assemble_version = assemble.viralngs_version
    String  scaffold_viral_assemble_version = scaffold.viralngs_version
    String  refine_viral_assemble_version   = refine_2x_and_plot.viralngs_version
    String? isnvs_viral_phylo_version       = isnvs_per_sample.viralngs_version
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




task taxon_filter__filter_to_taxon {
  meta { description: "This step reduces the read set to a specific taxon (usually the genus level or greater for the virus of interest)" }

  input {
    File     reads_unmapped_bam
    File     lastal_db_fasta
    Boolean? error_on_reads_in_neg_control = false
    Int?     negative_control_reads_threshold = 0
    String?  neg_control_prefixes_space_separated = "neg water NTC"

    Int?     machine_mem_gb
    String   docker="quay.io/broadinstitute/viral-classify:2.0.21.3"
  }

  # do this in two steps in case the input doesn't actually have "cleaned" in the name
  String   bam_basename = basename(basename(reads_unmapped_bam, ".bam"), ".cleaned")

  command {
    set -ex -o pipefail
    taxon_filter.py --version | tee VERSION

    # find 90% memory
    mem_in_mb=`/opt/viral-ngs/source/docker/calc_mem.py mb 90`

    if [[ "${error_on_reads_in_neg_control}" == "true" ]]; then
      ERROR_ON_NEG_CONTROL_ARGS="--errorOnReadsInNegControl"
      if [[ -n "${negative_control_reads_threshold}" ]]; then
        ERROR_ON_NEG_CONTROL_ARGS="$ERROR_ON_NEG_CONTROL_ARGS ${'--negativeControlReadsThreshold=' + negative_control_reads_threshold}"
      fi
      if [[ -n "${neg_control_prefixes_space_separated}" ]]; then
        ERROR_ON_NEG_CONTROL_ARGS="$ERROR_ON_NEG_CONTROL_ARGS ${'--negControlPrefixes=' + neg_control_prefixes_space_separated}"
      fi      
    fi

    taxon_filter.py filter_lastal_bam \
      ${reads_unmapped_bam} \
      ${lastal_db_fasta} \
      ${bam_basename}.taxfilt.bam \
      $ERROR_ON_NEG_CONTROL_ARGS \
      --JVMmemory="$mem_in_mb"m \
      --loglevel=DEBUG

    samtools view -c ${bam_basename}.taxfilt.bam | tee filter_read_count_post
    reports.py fastqc ${bam_basename}.taxfilt.bam ${bam_basename}.taxfilt_fastqc.html --out_zip ${bam_basename}.taxfilt_fastqc.zip
  }

  output {
    File   taxfilt_bam            = "${bam_basename}.taxfilt.bam"
    File   taxfilt_fastqc         = "${bam_basename}.taxfilt_fastqc.html"
    File   taxfilt_fastqc_zip     = "${bam_basename}.taxfilt_fastqc.zip"
    Int    filter_read_count_post = read_int("filter_read_count_post")
    String viralngs_version       = read_string("VERSION")
  }
  runtime {
    docker: "${docker}"
    memory: select_first([machine_mem_gb, 15]) + " GB"
    cpu: 16
    disks: "local-disk 375 LOCAL"
    dx_instance_type: "mem1_ssd1_v2_x8"
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




task assembly__scaffold {
    input {
      File         contigs_fasta
      File         reads_bam
      Array[File]+ reference_genome_fasta

      String?      aligner
      Float?       min_length_fraction
      Float?       min_unambig
      Int?         replace_length=55

      Int?         nucmer_max_gap
      Int?         nucmer_min_match
      Int?         nucmer_min_cluster
      Float?       scaffold_min_pct_contig_aligned

      Int?         machine_mem_gb
      String       docker="quay.io/broadinstitute/viral-assemble:2.0.21.0"

      # do this in multiple steps in case the input doesn't actually have "assembly1-x" in the name
      String       sample_name = basename(basename(basename(contigs_fasta, ".fasta"), ".assembly1-trinity"), ".assembly1-spades")
    }

    command {
        set -ex -o pipefail

        # find 90% memory
        mem_in_gb=`/opt/viral-ngs/source/docker/calc_mem.py gb 90`

        assembly.py --version | tee VERSION

        assembly.py order_and_orient \
          ${contigs_fasta} \
          ${sep=' ' reference_genome_fasta} \
          ${sample_name}.intermediate_scaffold.fasta \
          ${'--maxgap=' + nucmer_max_gap} \
          ${'--minmatch=' + nucmer_min_match} \
          ${'--mincluster=' + nucmer_min_cluster} \
          ${'--min_pct_contig_aligned=' + scaffold_min_pct_contig_aligned} \
          --outReference ${sample_name}.scaffolding_chosen_ref.fasta \
          --outStats ${sample_name}.scaffolding_stats.txt \
          --outAlternateContigs ${sample_name}.scaffolding_alt_contigs.fasta \
          --loglevel=DEBUG

        grep '^>' ${sample_name}.scaffolding_chosen_ref.fasta | cut -c 2- | tr '\n' '\t' > ${sample_name}.scaffolding_chosen_ref.txt

        assembly.py gapfill_gap2seq \
          ${sample_name}.intermediate_scaffold.fasta \
          ${reads_bam} \
          ${sample_name}.intermediate_gapfill.fasta \
          --memLimitGb $mem_in_gb \
          --maskErrors \
          --loglevel=DEBUG

        grep -v '^>' ${sample_name}.intermediate_gapfill.fasta | tr -d '\n' | wc -c | tee assembly_preimpute_length
        grep -v '^>' ${sample_name}.intermediate_gapfill.fasta | tr -d '\nNn' | wc -c | tee assembly_preimpute_length_unambiguous

        assembly.py impute_from_reference \
          ${sample_name}.intermediate_gapfill.fasta \
          ${sample_name}.scaffolding_chosen_ref.fasta \
          ${sample_name}.scaffolded_imputed.fasta \
          --newName ${sample_name} \
          ${'--replaceLength=' + replace_length} \
          ${'--minLengthFraction=' + min_length_fraction} \
          ${'--minUnambig=' + min_unambig} \
          ${'--aligner=' + aligner} \
          --loglevel=DEBUG
    }

    output {
        File   scaffold_fasta                        = "${sample_name}.scaffolded_imputed.fasta"
        File   intermediate_scaffold_fasta           = "${sample_name}.intermediate_scaffold.fasta"
        File   intermediate_gapfill_fasta            = "${sample_name}.intermediate_gapfill.fasta"
        Int    assembly_preimpute_length             = read_int("assembly_preimpute_length")
        Int    assembly_preimpute_length_unambiguous = read_int("assembly_preimpute_length_unambiguous")
        String scaffolding_chosen_ref_name           = read_string("${sample_name}.scaffolding_chosen_ref.txt")
        File   scaffolding_chosen_ref                = "${sample_name}.scaffolding_chosen_ref.fasta"
        File   scaffolding_stats                     = "${sample_name}.scaffolding_stats.txt"
        File   scaffolding_alt_contigs               = "${sample_name}.scaffolding_alt_contigs.fasta"
        String viralngs_version                      = read_string("VERSION")
    }

    runtime {
        docker: "${docker}"
        memory: select_first([machine_mem_gb, 15]) + " GB"
        cpu: 4
        disks: "local-disk 375 LOCAL"
        dx_instance_type: "mem1_ssd1_v2_x8"
    }
}




task assembly__refine_2x_and_plot {
    meta {
      description: "This combined task exists just to streamline the two calls to assembly.refine and one call to reports.plot_coverage that many denovo assembly workflows use. It saves on instance spin up and docker pull times, file staging time, and all steps contained here have similar hardware requirements. The more atomic WDL tasks are still available for custom workflows (see refine, refine_assembly_with_aligned_reads, align_reads, etc)."
    }

    input {
      File    assembly_fasta
      File    reads_unmapped_bam

      File?   novocraft_license

      String? refine1_novoalign_options="-r Random -l 30 -g 40 -x 20 -t 502"
      Float?  refine1_major_cutoff=0.5
      Int?    refine1_min_coverage=2

      String? refine2_novoalign_options="-r Random -l 40 -g 40 -x 20 -t 100"
      Float?  refine2_major_cutoff=0.5
      Int?    refine2_min_coverage=3

      String? plot_coverage_novoalign_options="-r Random -l 40 -g 40 -x 20 -t 100 -k"

      Int?    machine_mem_gb
      String  docker="quay.io/broadinstitute/viral-assemble:2.0.21.0"

      # do this in two steps in case the input doesn't actually have "cleaned" in the name
      String  sample_name = basename(basename(reads_unmapped_bam, ".bam"), ".cleaned")
    }

    command {
        set -ex -o pipefail

        # find 90% memory
        mem_in_mb=`/opt/viral-ngs/source/docker/calc_mem.py mb 90`

        assembly.py --version | tee VERSION

        ln -s ${assembly_fasta} assembly.fasta
        read_utils.py novoindex \
        assembly.fasta \
        ${"--NOVOALIGN_LICENSE_PATH=" + novocraft_license} \
        --loglevel=DEBUG

        # refine 1
        assembly.py refine_assembly \
          assembly.fasta \
          ${reads_unmapped_bam} \
          ${sample_name}.refine1.fasta \
          --outVcf ${sample_name}.refine1.pre_fasta.vcf.gz \
          --min_coverage ${refine1_min_coverage} \
          --major_cutoff ${refine1_major_cutoff} \
          --novo_params="${refine1_novoalign_options}" \
          --JVMmemory "$mem_in_mb"m \
          ${"--NOVOALIGN_LICENSE_PATH=" + novocraft_license} \
          --loglevel=DEBUG

        # refine 2
        assembly.py refine_assembly \
          ${sample_name}.refine1.fasta \
          ${reads_unmapped_bam} \
          ${sample_name}.fasta \
          --outVcf ${sample_name}.refine2.pre_fasta.vcf.gz \
          --min_coverage ${refine2_min_coverage} \
          --major_cutoff ${refine2_major_cutoff} \
          --novo_params="${refine2_novoalign_options}" \
          ${"--NOVOALIGN_LICENSE_PATH=" + novocraft_license} \
          --JVMmemory "$mem_in_mb"m \
          --loglevel=DEBUG

        # final alignment
        read_utils.py align_and_fix \
          ${reads_unmapped_bam} \
          ${sample_name}.fasta \
          --outBamAll ${sample_name}.all.bam \
          --outBamFiltered ${sample_name}.mapped.bam \
          --aligner_options "${plot_coverage_novoalign_options}" \
          --JVMmemory "$mem_in_mb"m \
          ${"--NOVOALIGN_LICENSE_PATH=" + novocraft_license} \
          --loglevel=DEBUG

        # collect figures of merit
        grep -v '^>' ${sample_name}.fasta | tr -d '\n' | wc -c | tee assembly_length
        grep -v '^>' ${sample_name}.fasta | tr -d '\nNn' | wc -c | tee assembly_length_unambiguous
        samtools view -c ${sample_name}.mapped.bam | tee reads_aligned
        # report only primary alignments 260=exclude unaligned reads and secondary mappings
        samtools view -h -F 260 ${sample_name}.all.bam | samtools flagstat - | tee ${sample_name}.all.bam.flagstat.txt
        grep properly ${sample_name}.all.bam.flagstat.txt | cut -f 1 -d ' ' | tee read_pairs_aligned
        samtools view ${sample_name}.mapped.bam | cut -f10 | tr -d '\n' | wc -c | tee bases_aligned
        #echo $(( $(cat bases_aligned) / $(cat assembly_length) )) | tee mean_coverage
        python -c "print (float("`cat bases_aligned`")/"`cat assembly_length`") if "`cat assembly_length`">0 else 0" > mean_coverage

        # fastqc mapped bam
        reports.py fastqc ${sample_name}.mapped.bam ${sample_name}.mapped_fastqc.html --out_zip ${sample_name}.mapped_fastqc.zip

        # plot coverage
        if [ $(cat reads_aligned) != 0 ]; then
          reports.py plot_coverage \
            ${sample_name}.mapped.bam \
            ${sample_name}.coverage_plot.pdf \
            --outSummary "${sample_name}.coverage_plot.txt" \
            --plotFormat pdf \
            --plotWidth 1100 \
            --plotHeight 850 \
            --plotDPI 100 \
            --plotTitle "${sample_name} coverage plot" \
            --loglevel=DEBUG
        else
          touch ${sample_name}.coverage_plot.pdf ${sample_name}.coverage_plot.txt
        fi
    }

    output {
        File refine1_sites_vcf_gz          = "${sample_name}.refine1.pre_fasta.vcf.gz"
        File refine1_assembly_fasta        = "${sample_name}.refine1.fasta"
        File refine2_sites_vcf_gz          = "${sample_name}.refine2.pre_fasta.vcf.gz"
        File final_assembly_fasta          = "${sample_name}.fasta"
        File aligned_bam                   = "${sample_name}.all.bam"
        File aligned_bam_idx               = "${sample_name}.all.bai"
        File aligned_bam_flagstat          = "${sample_name}.all.bam.flagstat.txt"
        File aligned_only_reads_bam        = "${sample_name}.mapped.bam"
        File aligned_only_reads_bam_idx    = "${sample_name}.mapped.bai"
        File aligned_only_reads_fastqc     = "${sample_name}.mapped_fastqc.html"
        File aligned_only_reads_fastqc_zip = "${sample_name}.mapped_fastqc.zip"
        File coverage_plot                 = "${sample_name}.coverage_plot.pdf"
        File coverage_tsv                  = "${sample_name}.coverage_plot.txt"
        Int  assembly_length               = read_int("assembly_length")
        Int  assembly_length_unambiguous   = read_int("assembly_length_unambiguous")
        Int  reads_aligned                 = read_int("reads_aligned")
        Int  read_pairs_aligned            = read_int("read_pairs_aligned")
        Int  bases_aligned                 = read_int("bases_aligned")
        Float mean_coverage                = read_float("mean_coverage")
        String viralngs_version            = read_string("VERSION")
    }

    runtime {
        docker: "${docker}"
        memory: select_first([machine_mem_gb, 7]) + " GB"
        cpu: 8
        disks: "local-disk 375 LOCAL"
        dx_instance_type: "mem1_ssd1_v2_x8"
    }
}




task intrahost__isnvs_per_sample {
  input {
    File    mapped_bam
    File    assembly_fasta

    Int?    threads
    Int?    minReadsPerStrand
    Int?    maxBias

    Int?    machine_mem_gb
    String  docker="quay.io/broadinstitute/viral-phylo:2.0.21.0"

    String  sample_name = basename(basename(basename(mapped_bam, ".bam"), ".all"), ".mapped")
  }

  command {
    intrahost.py --version | tee VERSION
    intrahost.py vphaser_one_sample \
        ${mapped_bam} \
        ${assembly_fasta} \
        vphaser2.${sample_name}.txt.gz \
        ${'--vphaserNumThreads' + threads} \
        --removeDoublyMappedReads \
        ${'--minReadsEach' + minReadsPerStrand} \
        ${'--maxBias' + maxBias}
  }

  output {
    File   isnvsFile        = "vphaser2.${sample_name}.txt.gz"
    String viralngs_version = read_string("VERSION")
  }
  runtime {
    docker: "${docker}"
    memory: select_first([machine_mem_gb, 7]) + " GB"
    dx_instance_type: "mem1_ssd1_v2_x8"
  }
}


