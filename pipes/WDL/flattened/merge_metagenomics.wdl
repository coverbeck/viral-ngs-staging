version 1.0




workflow merge_metagenomics {

    input {
        Array[File] krakenuniq_summary_reports
        Array[File] krona_per_sample
    }

    call metagenomics__krona_merge as krona_merge {
        input:
            krona_reports = krona_per_sample,
            out_basename = "merged_krona"
    }

    call reports__aggregate_metagenomics_reports as aggregate_metagenomics_reports {
        input:
            kraken_summary_reports = krakenuniq_summary_reports
    }

    output {
        File krona_merged         = krona_merge.krona_report_html
        File metagenomics_summary = aggregate_metagenomics_reports.krakenuniq_aggregate_taxlevel_summary
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


