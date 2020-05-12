version 1.0



workflow scaffold_and_refine {
    input {
        File reads_unmapped_bam
    }

    call assembly__scaffold as scaffold {
        input:
            reads_bam = reads_unmapped_bam
    }

    call assembly__refine_2x_and_plot as refine_2x_and_plot {
        input:
            assembly_fasta = scaffold.scaffold_fasta,
            reads_unmapped_bam = reads_unmapped_bam
    }

  output {
    File  final_assembly_fasta        = refine_2x_and_plot.final_assembly_fasta
    File  aligned_only_reads_bam      = refine_2x_and_plot.aligned_only_reads_bam
    File  coverage_plot               = refine_2x_and_plot.coverage_plot
    Int   assembly_length             = refine_2x_and_plot.assembly_length
    Int   assembly_length_unambiguous = refine_2x_and_plot.assembly_length_unambiguous
    Int   reads_aligned               = refine_2x_and_plot.reads_aligned
    Float mean_coverage               = refine_2x_and_plot.mean_coverage

    File   scaffold_fasta                        = scaffold.scaffold_fasta
    File   intermediate_scaffold_fasta           = scaffold.intermediate_scaffold_fasta
    File   intermediate_gapfill_fasta            = scaffold.intermediate_gapfill_fasta
    Int    assembly_preimpute_length             = scaffold.assembly_preimpute_length
    Int    assembly_preimpute_length_unambiguous = scaffold.assembly_preimpute_length_unambiguous
    String scaffolding_chosen_ref_name           = scaffold.scaffolding_chosen_ref_name
    File   scaffolding_stats                     = scaffold.scaffolding_stats
    File   scaffolding_alt_contigs               = scaffold.scaffolding_alt_contigs

    File aligned_bam                   = refine_2x_and_plot.aligned_bam
    File aligned_only_reads_bam_idx    = refine_2x_and_plot.aligned_only_reads_bam_idx
    File aligned_only_reads_fastqc     = refine_2x_and_plot.aligned_only_reads_fastqc
    File coverage_tsv                  = refine_2x_and_plot.coverage_tsv
    Int  read_pairs_aligned            = refine_2x_and_plot.read_pairs_aligned
    Int  bases_aligned                 = refine_2x_and_plot.bases_aligned

    String scaffold_viral_assemble_version = scaffold.viralngs_version
    String refine_viral_assemble_version   = refine_2x_and_plot.viralngs_version
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


