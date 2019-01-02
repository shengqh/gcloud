# WORKFLOW DEFINITION
workflow dnaseq {
  String base_file_name
  
  File input_bam

  Int preemptible_tries
  Int agg_preemptible_tries

  # Optional input to increase all disk sizes in case of outlier sample with strange size behavior
  Int? increase_disk_size

  # Some tasks need wiggle room, and we also need to add a small amount of disk to prevent getting a
  # Cromwell error from asking for 0 disk when the input is less than 1GB
  Int additional_disk = select_first([increase_disk_size, 20])

  # Mark Duplicates takes in as input readgroup bams and outputs a slightly smaller aggregated bam. Giving .25 as wiggleroom
  Float md_disk_multiplier = 2.25

  Float mapped_bam_size = size(input_bam, "GB")

  call MarkDuplicates {
    input:
      input_bam = input_bam,
      output_bam_basename = base_file_name + ".aligned.unsorted.duplicates_marked",
      metrics_filename = base_file_name + ".duplicate_metrics",
      # The merged bam will be smaller than the sum of the parts so we need to account for the unmerged inputs
      # and the merged output.
      disk_size = (md_disk_multiplier * mapped_bam_size) + additional_disk,
      preemptible_tries = agg_preemptible_tries
  }

  # Outputs that will be retained when execution is complete
  output {
    File duplicate_metrics = MarkDuplicates.duplicate_metrics
    File duplicate_bam = MarkDuplicates.output_bam
  }
}

# TASK DEFINITIONS

# Mark duplicate reads to avoid counting non-independent observations
task MarkDuplicates {
  File input_bam
  String output_bam_basename
  String metrics_filename
  Float disk_size
  Int preemptible_tries

  # The program default for READ_NAME_REGEX is appropriate in nearly every case.
  # Sometimes we wish to supply "null" in order to turn off optical duplicate detection
  # This can be desirable if you don't mind the estimated library size being wrong and optical duplicate detection is taking >7 days and failing
  String? read_name_regex

 # Task is assuming query-sorted input so that the Secondary and Supplementary reads get marked correctly
 # This works because the output of BWA is query-grouped and therefore, so is the output of MergeBamAlignment.
 # While query-grouped isn't actually query-sorted, it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname"
  command {
    java -Xmx40g -jar /usr/gitc/picard.jar \
      MarkDuplicates \
      INPUT=${input_bam} \
      OUTPUT=${output_bam_basename}.bam \
      METRICS_FILE=${metrics_filename} \
      VALIDATION_STRINGENCY=SILENT 
  }
  runtime {
    preemptible: preemptible_tries
    memory: "40 GB"
    disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.1-1540490856"
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
    File duplicate_metrics = "${metrics_filename}"
  }
}
