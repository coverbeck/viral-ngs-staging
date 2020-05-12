version 1.0



workflow filter_classified_bam_to_taxa {
    call metagenomics__filter_bam_to_taxa as filter_bam_to_taxa
    output {
        File    bam_filtered_to_taxa                        = filter_bam_to_taxa.bam_filtered_to_taxa
        Int     classified_taxonomic_filter_read_count_pre  = filter_bam_to_taxa.classified_taxonomic_filter_read_count_pre
        Int     classified_taxonomic_filter_read_count_post = filter_bam_to_taxa.classified_taxonomic_filter_read_count_post
        String  viral_classify_version                      = filter_bam_to_taxa.viralngs_version
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


