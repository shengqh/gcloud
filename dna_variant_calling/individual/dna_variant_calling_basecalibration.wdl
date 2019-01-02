# WORKFLOW DEFINITION
workflow dnaseq {
  String base_file_name
  
  File input_bam
  File input_bam_index

  File ref_fasta
  File ref_fasta_index
  File ref_dict

  File dbSNP_vcf
  File dbSNP_vcf_index
  Array[File] known_indels_sites_VCFs
  Array[File] known_indels_sites_indices

  Int preemptible_tries
  Int agg_preemptible_tries

  # Optional input to increase all disk sizes in case of outlier sample with strange size behavior
  Int? increase_disk_size

  # Some tasks need wiggle room, and we also need to add a small amount of disk to prevent getting a
  # Cromwell error from asking for 0 disk when the input is less than 1GB
  Int additional_disk = select_first([increase_disk_size, 20])

  # Mark Duplicates takes in as input readgroup bams and outputs a slightly smaller aggregated bam. Giving .25 as wiggleroom
  Float md_disk_multiplier = 2.25

  Float agg_bam_size = size(input_bam, "GB")

  Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB")
  Float dbsnp_size = size(dbSNP_vcf, "GB")

  call BaseRecalibrator {
    input:
      input_bam = input_bam,
      recalibration_report_filename = base_file_name + ".recal_data.csv",
      dbSNP_vcf = dbSNP_vcf,
      dbSNP_vcf_index = dbSNP_vcf_index,
      known_indels_sites_VCFs = known_indels_sites_VCFs,
      known_indels_sites_indices = known_indels_sites_indices,
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      disk_size = agg_bam_size + ref_size + dbsnp_size + additional_disk,
      preemptible_tries = agg_preemptible_tries
  }
  
  # Outputs that will be retained when execution is complete
  output {
    File recalibration_report = BaseRecalibrator.recalibration_report
  }
}

# TASK DEFINITIONS
task BaseRecalibrator {
  File input_bam
  String recalibration_report_filename
  File dbSNP_vcf
  File dbSNP_vcf_index
  Array[File] known_indels_sites_VCFs
  Array[File] known_indels_sites_indices
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  Float disk_size
  Int preemptible_tries

  command {
    /gatk/gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
      -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
      -Xloggc:gc_log.log -Xms4000m" \
      BaseRecalibrator \
      -R ${ref_fasta} \
      -I ${input_bam} \
      --use-original-qualities \
      -O ${recalibration_report_filename} \
      --known-sites ${dbSNP_vcf} \
      --known-sites ${sep=" --known-sites " known_indels_sites_VCFs}
  }
  runtime {
    preemptible: preemptible_tries
    memory: "6 GB"
    disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
    docker: "broadinstitute/gatk:4.0.11.0"
  }
  output {
    File recalibration_report = "${recalibration_report_filename}"
  }
}
