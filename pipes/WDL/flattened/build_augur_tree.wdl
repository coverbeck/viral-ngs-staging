version 1.0



workflow build_augur_tree {
    meta {
        description: "Align assemblies, build trees, and convert to json representation suitable for Nextstrain visualization. See https://nextstrain.org/docs/getting-started/ and https://nextstrain-augur.readthedocs.io/en/stable/"
    }

    input {
        Array[File]     assembly_fastas
        File            metadata
        String          virus
        File            ref_fasta
        File            genbank_gb
        Array[String]?  ancestral_traits_to_infer
    }

    parameter_meta {
        assembly_fastas: {
          description: "Set of assembled genomes to align and build trees. These must represent a single chromosome/segment of a genome only. Fastas may be one-sequence-per-individual or a concatenated multi-fasta (unaligned) or a mixture of the two. Fasta header records need to be pipe-delimited (|) for each metadata value.",
          patterns: ["*.fasta", "*.fa"]
        }
        metadata: {
          description: "Metadata in tab-separated text format. See https://nextstrain-augur.readthedocs.io/en/stable/faq/metadata.html for details.",
          patterns: ["*.txt", "*.tsv"]
        }
        virus: {
          description: "A filename-friendly string that is used as a base for output file names."
        }
        ref_fasta: {
          description: "A reference assembly (not included in assembly_fastas) to align assembly_fastas against. Typically from NCBI RefSeq or similar.",
          patterns: ["*.fasta", "*.fa"]
        }
        genbank_gb: {
          description: "A 'genbank' formatted gene annotation file that is used to calculate coding consequences of observed mutations. Must correspond to the same coordinate space as ref_fasta. Typically downloaded from the same NCBI accession number as ref_fasta.",
          patterns: ["*.gb", "*.gbf"]
        }
        ancestral_traits_to_infer: {
          description: "A list of metadata traits to use for ancestral node inference (see https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/traits.html). Multiple traits may be specified; must correspond exactly to column headers in metadata file. Omitting these values will skip ancestral trait inference, and ancestral nodes will not have estimated values for metadata."
        }
    }

    call nextstrain__concatenate as concatenate {
        input:
            infiles     = assembly_fastas,
            output_name = "all_samples_combined_assembly.fasta"
    }
    call nextstrain__augur_mafft_align as augur_mafft_align {
        input:
            sequences = concatenate.combined,
            ref_fasta = ref_fasta,
            basename  = virus
    }
    call nextstrain__draft_augur_tree as draft_augur_tree {
        input:
            aligned_fasta  = augur_mafft_align.aligned_sequences,
            basename       = virus
    }
    call nextstrain__refine_augur_tree as refine_augur_tree {
        input:
            raw_tree       = draft_augur_tree.aligned_tree,
            aligned_fasta  = augur_mafft_align.aligned_sequences,
            metadata       = metadata,
            basename       = virus
    }
    if(defined(ancestral_traits_to_infer) && length(select_first([ancestral_traits_to_infer,[]]))>0) {
        call nextstrain__ancestral_traits as ancestral_traits {
            input:
                tree           = refine_augur_tree.tree_refined,
                metadata       = metadata,
                columns        = select_first([ancestral_traits_to_infer,[]]),
                basename       = virus
        }
    }
    call nextstrain__ancestral_tree as ancestral_tree {
        input:
            refined_tree   = refine_augur_tree.tree_refined,
            aligned_fasta  = augur_mafft_align.aligned_sequences,
            basename       = virus
    }
    call nextstrain__translate_augur_tree as translate_augur_tree {
        input:
            basename       = virus,
            refined_tree   = refine_augur_tree.tree_refined,
            nt_muts        = ancestral_tree.nt_muts_json,
            genbank_gb     = genbank_gb
    }
    call nextstrain__export_auspice_json as export_auspice_json {
        input:
            refined_tree   = refine_augur_tree.tree_refined,
            metadata       = metadata,
            branch_lengths = refine_augur_tree.branch_lengths,
            traits         = ancestral_traits.node_data_json,
            nt_muts        = ancestral_tree.nt_muts_json,
            aa_muts        = translate_augur_tree.aa_muts_json,
            basename       = virus
    }

    output {
        File  combined_assembly_fasta    = concatenate.combined
        File  augur_aligned_fasta        = augur_mafft_align.aligned_sequences
        File  raw_tree                   = draft_augur_tree.aligned_tree
        File  refined_tree               = refine_augur_tree.tree_refined
        File  branch_lengths             = refine_augur_tree.branch_lengths
        File  json_nt_muts               = ancestral_tree.nt_muts_json
        File  ancestral_sequences_fasta  = ancestral_tree.sequences
        File  json_aa_muts               = translate_augur_tree.aa_muts_json
        File? json_ancestral_traits      = ancestral_traits.node_data_json
        File  auspice_input_json         = export_auspice_json.virus_json
    }
}



task nextstrain__concatenate {
    meta {
        description: "This is nothing more than unix cat."
    }
    input {
        Array[File] infiles
        String      output_name
    }
    command {
        cat ~{sep=" " infiles} > "${output_name}"
    }
    runtime {
        docker: "ubuntu"
        memory: "1 GB"
        cpu:    1
        disks: "local-disk 375 LOCAL"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
    output {
        File combined = "${output_name}"
    }
}




task nextstrain__augur_mafft_align {
    meta {
        description: "Align multiple sequences from FASTA. See https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/align.html"
    }
    input {
        File     sequences
        File     ref_fasta
        String   basename

        File?    existing_alignment
        Boolean  fill_gaps = true
        Boolean  remove_reference = true

        Int?     machine_mem_gb
        String   docker = "nextstrain/base:build-20200506T095107Z"
    }
    command {
        augur align --sequences ~{sequences} \
            --reference-sequence ~{ref_fasta} \
            --output ~{basename}_aligned.fasta \
            ~{true="--fill-gaps" false="" fill_gaps} \
            ~{"--existing-alignment " + existing_alignment} \
            ~{true="--remove-reference" false="" remove_reference} \
            --debug \
            --nthreads auto
        cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes > MAX_RAM
    }
    runtime {
        docker: docker
        memory: select_first([machine_mem_gb, 104]) + " GB"
        cpu :   16
        disks:  "local-disk 375 LOCAL"
        preemptible: 2
        dx_instance_type: "mem3_ssd2_v2_x16"
    }
    output {
        File aligned_sequences = "~{basename}_aligned.fasta"
        File align_troubleshoot = stdout()
        File max_ram_usage_in_bytes = "MAX_RAM"
    }
}




task nextstrain__draft_augur_tree {
    meta {
        description: "Build a tree using a variety of methods. See https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/tree.html"
    }
    input {
        File     aligned_fasta
        String   basename

        String   method = "iqtree"
        String   substitution_model = "GTR"
        File?    exclude_sites
        File?    vcf_reference
        String?  tree_builder_args

        Int?     machine_mem_gb
        String   docker = "nextstrain/base:build-20200506T095107Z"
    }
    command {
        augur tree --alignment ~{aligned_fasta} \
            --output ~{basename}_raw_tree.nwk \
            --method ~{default="iqtree" method} \
            --substitution-model ~{default="GTR" substitution_model} \
            ~{"--exclude-sites " + exclude_sites} \
            ~{"--vcf-reference " + vcf_reference} \
            ~{"--tree-builder-args " + tree_builder_args} \
            --nthreads auto
    }
    runtime {
        docker: docker
        memory: select_first([machine_mem_gb, 30]) + " GB"
        cpu :   16
        disks:  "local-disk 375 LOCAL"
        dx_instance_type: "mem1_ssd1_v2_x16"
        preemptible: 2
    }
    output {
        File aligned_tree = "~{basename}_raw_tree.nwk"
    }
}




task nextstrain__refine_augur_tree {
    meta {
        description: "Refine an initial tree using sequence metadata. See https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/refine.html"
    }
    input {
        File     raw_tree
        File     aligned_fasta
        File     metadata
        String   basename

        Int?     gen_per_year
        Float?   clock_rate
        Float?   clock_std_dev
        Boolean  keep_root = false
        String?  root
        Boolean? covariance
        Boolean  keep_polytomies = false
        Int?     precision
        Boolean  date_confidence = true
        String?  date_inference = "marginal"
        String?  branch_length_inference
        String?  coalescent
        Int?     clock_filter_iqd = 4
        String?  divergence_units
        File?    vcf_reference

        Int?     machine_mem_gb
        String   docker = "nextstrain/base:build-20200506T095107Z"
    }
    command {
        augur refine \
            --tree ~{raw_tree} \
            --alignment ~{aligned_fasta} \
            --metadata ~{metadata} \
            --output-tree ~{basename}_refined_tree.nwk \
            --output-node-data ~{basename}_branch_lengths.json \
            --timetree \
            ~{"--clock-rate " + clock_rate} \
            ~{"--clock-std-dev " + clock_std_dev} \
            ~{"--coalescent " + coalescent} \
            ~{"--clock-filter-iqd " + clock_filter_iqd} \
            ~{"--gen-per-year " + gen_per_year} \
            ~{"--root " + root} \
            ~{"--precision " + precision} \
            ~{"--date-inference " + date_inference} \
            ~{"--branch-length-inference " + branch_length_inference} \
            ~{"--divergence-units " + divergence_units} \
            ~{true="--covariance" false="--no-covariance" covariance} \
            ~{true="--keep-root" false="" keep_root} \
            ~{true="--keep-polytomies" false="" keep_polytomies} \
            ~{true="--date-confidence" false="" date_confidence} \
            ~{"--vcf-reference " + vcf_reference}
    }
    runtime {
        docker: docker
        memory: select_first([machine_mem_gb, 30]) + " GB"
        cpu :   16
        disks: "local-disk 375 LOCAL"
        dx_instance_type: "mem1_ssd1_v2_x16"
        preemptible: 2
    }
    output {
        File tree_refined  = "~{basename}_refined_tree.nwk"
        File branch_lengths = "~{basename}_branch_lengths.json"
    }
}




task nextstrain__ancestral_traits {
    meta {
        description: "Infer ancestral traits based on a tree. See https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/traits.html"
    }
    input {
        File           tree
        File           metadata
        Array[String]+ columns
        String         basename

        Boolean        confidence = true
        File?          weights
        Float?         sampling_bias_correction

        Int?     machine_mem_gb
        String   docker = "nextstrain/base:build-20200506T095107Z"
    }
    command {
        augur traits \
            --tree ~{tree} \
            --metadata ~{metadata} \
            --columns ~{sep=" " columns} \
            --output-node-data "~{basename}_nodes.json" \
            ~{"--weights " + weights} \
            ~{true="--confidence" false="" confidence}
    }
    runtime {
        docker: docker
        memory: select_first([machine_mem_gb, 3]) + " GB"
        cpu :   2
        disks: "local-disk 50 HDD"
        dx_instance_type: "mem1_ssd1_v2_x2"
        preemptible: 2
    }
    output {
        File node_data_json = "~{basename}_nodes.json"
    }
}




task nextstrain__ancestral_tree {
    meta {
        description: "Infer ancestral sequences based on a tree. See https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/ancestral.html"
    }
    input {
        File     refined_tree
        File     aligned_fasta
        String   basename

        String   inference = "joint"
        Boolean  keep_ambiguous = false
        Boolean  infer_ambiguous = false
        Boolean  keep_overhangs = false
        File?    vcf_reference
        File?    output_vcf

        Int?     machine_mem_gb
        String   docker = "nextstrain/base:build-20200506T095107Z"
    }
    command {
        augur ancestral \
            --tree ~{refined_tree} \
            --alignment ~{aligned_fasta} \
            --output-node-data ~{basename}_nt_muts.json \
            ~{"--vcf-reference " + vcf_reference} \
            ~{"--output-vcf " + output_vcf} \
            --output-sequences ~{basename}_ancestral_sequences.fasta \
            ~{true="--keep-overhangs" false="" keep_overhangs} \
            --inference ~{default="joint" inference} \
            ~{true="--keep-ambiguous" false="" keep_ambiguous} \
            ~{true="--infer-ambiguous" false="" infer_ambiguous}
    }
    runtime {
        docker: docker
        memory: select_first([machine_mem_gb, 7]) + " GB"
        cpu :   4
        disks: "local-disk 50 HDD"
        dx_instance_type: "mem1_ssd1_v2_x4"
        preemptible: 2
    }
    output {
        File nt_muts_json = "~{basename}_nt_muts.json"
        File sequences    = "~{basename}_ancestral_sequences.fasta"
    }
}




task nextstrain__translate_augur_tree {
    meta {
        description: "export augur files to json suitable for auspice visualization. See https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/translate.html"
    }
    input {
        String basename
        File   refined_tree
        File   nt_muts
        File   genbank_gb

        File?  genes
        File?  vcf_reference_output
        File?  vcf_reference

        Int?   machine_mem_gb
        String docker = "nextstrain/base:build-20200506T095107Z"
    }
    command {
        augur translate --tree ~{refined_tree} \
            --ancestral-sequences ~{nt_muts} \
            --reference-sequence ~{genbank_gb} \
            ~{"--vcf-reference-output " + vcf_reference_output} \
            ~{"--vcf-reference " + vcf_reference} \
            ~{"--genes " + genes} \
            --output-node-data ~{basename}_aa_muts.json
    }
    runtime {
        docker: docker
        memory: select_first([machine_mem_gb, 3]) + " GB"
        cpu :   2
        disks:  "local-disk 50 HDD"
        dx_instance_type: "mem1_ssd1_v2_x2"
        preemptible: 2
    }
    output {
        File aa_muts_json = "~{basename}_aa_muts.json"
    }
}




task nextstrain__export_auspice_json {
    meta {
        description: "export augur files to json suitable for auspice visualization. See https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/export.html"
    }
    input {
        File   auspice_config
        File   metadata
        File   refined_tree
        File?  branch_lengths
        File?  traits
        File?  nt_muts
        File?  aa_muts
        File?  lat_longs_tsv
        File?  colors_tsv
        String basename

        Int?   machine_mem_gb
        String docker = "nextstrain/base:build-20200506T095107Z"
    }
    command {

        NODE_DATA_FLAG=""
        if [[ -n "~{branch_lengths}" || -n "~{traits}" || -n "~{nt_muts}" || -n "~{aa_muts}" ]]; then
          NODE_DATA_FLAG="--node-data "
        fi

        augur export v2 --tree ~{refined_tree} \
            --metadata ~{metadata} \
            $NODE_DATA_FLAG ~{sep=' ' select_all([branch_lengths,traits,nt_muts,aa_muts])}\
            --auspice-config ~{auspice_config} \
            ~{"--lat-longs " + lat_longs_tsv} \
            ~{"--colors " + colors_tsv} \
            --output ~{basename}_auspice.json
    }
    runtime {
        docker: docker
        memory: select_first([machine_mem_gb, 3]) + " GB"
        cpu :   2
        disks:  "local-disk 100 HDD"
        dx_instance_type: "mem1_ssd1_v2_x2"
        preemptible: 2
    }
    output {
        File virus_json = "~{basename}_auspice.json"
    }
}


