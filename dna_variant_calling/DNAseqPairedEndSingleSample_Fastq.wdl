## This WDL pipeline implements data pre-processing and initial variant calling (GVCF
## generation) according to the GATK Best Practices (June 2016) for germline SNP and
## Indel discovery in human whole-genome sequencing (WGS) data.
##
## Inputs file format: fastq 

# WORKFLOW DEFINITION
workflow dnaseq {

  File contamination_sites_ud
  File contamination_sites_bed
  File contamination_sites_mu
  File wgs_evaluation_interval_list

  String base_file_name
  String final_gvcf_base_name
  String unmapped_fastq_suffix
  
  File raw_fastq1
  File raw_fastq2

  File wgs_calling_interval_list
  Int haplotype_scatter_count
  Int break_bands_at_multiples_of

  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File ref_alt
  File ref_bwt
  File ref_sa
  File ref_amb
  File ref_ann
  File ref_pac

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
  # Germline single sample GVCFs shouldn't get bigger even when the input bam is bigger (after a certain size)
  Int GVCF_disk_size = select_first([increase_disk_size, 30])
  # Sometimes the output is larger than the input, or a task can spill to disk. In these cases we need to account for the
  # input (1) and the output (1.5) or the input(1), the output(1), and spillage (.5).
  Float bwa_disk_multiplier = 2.5
  # SortSam spills to disk a lot more because we are only store 300000 records in RAM now because its faster for our data
  # so it needs more disk space.  Also it spills to disk in an uncompressed format so we need to account for that with a
  # larger multiplier
  Float sort_sam_disk_multiplier = 3.25

  # Mark Duplicates takes in as input readgroup bams and outputs a slightly smaller aggregated bam. Giving .25 as wiggleroom
  Float md_disk_multiplier = 2.25

  # These values set the thresholds for what is considered outside the normal realm of "reasonable" data.
  Float max_duplication_in_reasonable_sample = 0.30
  Float max_chimerism_in_reasonable_sample = 0.15

  String bwa_commandline

  String recalibrated_bam_basename = base_file_name + ".aligned.duplicates_marked.recalibrated"

  Int compression_level = 2

  # Get the size of the standard reference files as well as the additional reference files needed for BWA
  Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB")
  Float bwa_ref_size = ref_size + size(ref_alt, "GB") + size(ref_amb, "GB") + size(ref_ann, "GB") + size(ref_bwt, "GB") + size(ref_pac, "GB") + size(ref_sa, "GB")
  Float dbsnp_size = size(dbSNP_vcf, "GB")

  # Align flowcell-level unmapped input bams in parallel

  Float unmapped_fastq_size = size(raw_fastq1, "GB")+size(raw_fastq2, "GB")

  # Map reads to reference
  call FastqBwaMem {
    input:
      input_fastq1 = raw_fastq1,
      input_fastq2 = raw_fastq2,
      bwa_commandline = bwa_commandline,
      output_bam_basename = base_file_name + ".aligned.unsorted",
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      ref_alt = ref_alt,
      ref_bwt = ref_bwt,
      ref_amb = ref_amb,
      ref_ann = ref_ann,
      ref_pac = ref_pac,
      ref_sa = ref_sa,
      # The merged bam can be bigger than only the aligned bam,
      # so account for the output size by multiplying the input size by 2.75.
      disk_size = unmapped_fastq_size + bwa_ref_size + (bwa_disk_multiplier * unmapped_fastq_size) + additional_disk,
      compression_level = compression_level,
      preemptible_tries = preemptible_tries
  }

  Float agg_bam_size = size(FastqBwaMem.output_bam, "GB")
  
  # Sort aggregated+deduped BAM file and fix tags
  call SortSam as SortSampleBam {
    input:
      input_bam = FastqBwaMem.output_bam,
      output_bam_basename = base_file_name + ".aligned.sorted",
      # This task spills to disk so we need space for the input bam, the output bam, and any spillage.
      disk_size = (sort_sam_disk_multiplier * agg_bam_size) + additional_disk,
      compression_level = compression_level,
      preemptible_tries = agg_preemptible_tries
  }

  Float mapped_bam_size = size(SortSampleBam.output_bam, "GB")
  
  call MarkDuplicates {
    input:
      input_bam = SortSampleBam.output_bam,
      output_bam_basename = base_file_name + ".aligned.sorted.duplicates_marked",
      metrics_filename = base_file_name + ".duplicate_metrics",
      disk_size = (md_disk_multiplier * mapped_bam_size) + additional_disk,
      compression_level = compression_level,
      preemptible_tries = agg_preemptible_tries
  }

  # Estimate level of cross-sample contamination
  call CheckContamination {
    input:
      input_bam = MarkDuplicates.output_bam,
      input_bam_index = MarkDuplicates.output_bam_index,
      contamination_sites_ud = contamination_sites_ud,
      contamination_sites_bed = contamination_sites_bed,
      contamination_sites_mu = contamination_sites_mu,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      output_prefix = base_file_name + ".preBqsr",
      disk_size = agg_bam_size + ref_size + additional_disk,
      preemptible_tries = agg_preemptible_tries,
      contamination_underestimation_factor = 0.75
  }

  # Create list of sequences for scatter-gather parallelization
  call CreateSequenceGroupingTSV {
    input:
      ref_dict = ref_dict,
      preemptible_tries = preemptible_tries
  }

  # We need disk to localize the sharded input and output due to the scatter for BQSR.
  # If we take the number we are scattering by and reduce by 3 we will have enough disk space
  # to account for the fact that the data is not split evenly.
  Int num_of_bqsr_scatters = length(CreateSequenceGroupingTSV.sequence_grouping)
  Int potential_bqsr_divisor = num_of_bqsr_scatters - 10
  Int bqsr_divisor = if potential_bqsr_divisor > 1 then potential_bqsr_divisor else 1

 # Perform Base Quality Score Recalibration (BQSR) on the sorted BAM in parallel
  scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping) {
    # Generate the recalibration model by interval
    call BaseRecalibrator {
      input:
        input_bam = MarkDuplicates.output_bam,
        input_bam_index = MarkDuplicates.output_bam_index,
        recalibration_report_filename = base_file_name + ".recal_data.csv",
        sequence_group_interval = subgroup,
        dbSNP_vcf = dbSNP_vcf,
        dbSNP_vcf_index = dbSNP_vcf_index,
        known_indels_sites_VCFs = known_indels_sites_VCFs,
        known_indels_sites_indices = known_indels_sites_indices,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        # We need disk to localize the sharded bam due to the scatter.
        disk_size = (agg_bam_size / bqsr_divisor) + ref_size + dbsnp_size + additional_disk,
        preemptible_tries = agg_preemptible_tries
    }
  }

  # Merge the recalibration reports resulting from by-interval recalibration
  # The reports are always the same size
  call GatherBqsrReports {
    input:
      input_bqsr_reports = BaseRecalibrator.recalibration_report,
      output_report_filename = base_file_name + ".recal_data.csv",
      disk_size = additional_disk,
      preemptible_tries = preemptible_tries
  }

  scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping_with_unmapped) {
    # Apply the recalibration model by interval
    call ApplyBQSR {
      input:
        input_bam = MarkDuplicates.output_bam,
        input_bam_index = MarkDuplicates.output_bam_index,
        output_bam_basename = recalibrated_bam_basename,
        recalibration_report = GatherBqsrReports.output_bqsr_report,
        sequence_group_interval = subgroup,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        # We need disk to localize the sharded bam and the sharded output due to the scatter.
        disk_size = ((agg_bam_size * 3) / bqsr_divisor) + ref_size + additional_disk,
        compression_level = compression_level,
        preemptible_tries = agg_preemptible_tries
    }
  }

  # Merge the recalibrated BAM files resulting from by-interval recalibration
  call GatherBamFiles {
    input:
      input_bams = ApplyBQSR.recalibrated_bam,
      output_bam_basename = base_file_name,
      # Multiply the input bam size by two to account for the input and output
      disk_size = (2 * agg_bam_size) + additional_disk,
      compression_level = compression_level,
      preemptible_tries = agg_preemptible_tries
  }

  #BQSR bins the qualities which makes a significantly smaller bam
  Float binned_qual_bam_size = size(GatherBamFiles.output_bam, "GB")

  # Break the calling interval_list into sub-intervals
  # Perform variant calling on the sub-intervals, and then gather the results
  call ScatterIntervalList {
    input:
      interval_list = wgs_calling_interval_list,
      scatter_count = haplotype_scatter_count,
      break_bands_at_multiples_of = break_bands_at_multiples_of
  }

  # We need disk to localize the sharded input and output due to the scatter for HaplotypeCaller.
  # If we take the number we are scattering by and reduce by 20 we will have enough disk space
  # to account for the fact that the data is quite uneven across the shards.
  Int potential_hc_divisor = ScatterIntervalList.interval_count - 20
  Int hc_divisor = if potential_hc_divisor > 1 then potential_hc_divisor else 1

  # Call variants in parallel over WGS calling intervals
  scatter (index in range(ScatterIntervalList.interval_count)) {
    # Generate GVCF by interval
    call HaplotypeCaller {
      input:
        contamination = CheckContamination.contamination,
        input_bam = GatherBamFiles.output_bam,
        input_bai = GatherBamFiles.output_bam_index,
        interval_list = ScatterIntervalList.out[index],
        gvcf_basename = base_file_name,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        # Divide the total output GVCF size and the input bam size to account for the smaller scattered input and output.
        disk_size = ((binned_qual_bam_size + GVCF_disk_size) / hc_divisor) + ref_size + additional_disk,
        preemptible_tries = agg_preemptible_tries
     }
  }

  # Combine by-interval GVCFs into a single sample GVCF file
  call MergeVCFs {
    input:
      input_vcfs = HaplotypeCaller.output_gvcf,
      input_vcfs_indexes = HaplotypeCaller.output_gvcf_index,
      output_vcf_name = final_gvcf_base_name + ".g.vcf.gz",
      disk_size = GVCF_disk_size,
      preemptible_tries = agg_preemptible_tries
  }

  Float gvcf_size = size(MergeVCFs.output_vcf, "GB")

  # Validate the GVCF output of HaplotypeCaller
  call ValidateGVCF {
    input:
      input_vcf = MergeVCFs.output_vcf,
      input_vcf_index = MergeVCFs.output_vcf_index,
      dbSNP_vcf = dbSNP_vcf,
      dbSNP_vcf_index = dbSNP_vcf_index,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      wgs_calling_interval_list = wgs_calling_interval_list,
      disk_size = gvcf_size + ref_size + dbsnp_size + additional_disk,
      preemptible_tries = agg_preemptible_tries
  }

  # QC the GVCF
  call CollectGvcfCallingMetrics {
    input:
      input_vcf = MergeVCFs.output_vcf,
      input_vcf_index = MergeVCFs.output_vcf_index,
      metrics_basename = final_gvcf_base_name,
      dbSNP_vcf = dbSNP_vcf,
      dbSNP_vcf_index = dbSNP_vcf_index,
      ref_dict = ref_dict,
      wgs_evaluation_interval_list = wgs_evaluation_interval_list,
      disk_size = gvcf_size + dbsnp_size + additional_disk,
      preemptible_tries = agg_preemptible_tries
  }

  # Outputs that will be retained when execution is complete
  output {
    File selfSM = CheckContamination.selfSM
    Float contamination = CheckContamination.contamination

    File duplicate_metrics = MarkDuplicates.duplicate_metrics
    File output_bqsr_reports = GatherBqsrReports.output_bqsr_report

    File output_vcf = MergeVCFs.output_vcf
    File output_vcf_index = MergeVCFs.output_vcf_index

    File gvcf_summary_metrics = CollectGvcfCallingMetrics.summary_metrics
    File gvcf_detail_metrics = CollectGvcfCallingMetrics.detail_metrics
  }
}

# TASK DEFINITIONS

# Read FASTQ and stream to BWA MEM for alignment
task FastqBwaMem {
  File input_fastq1
  File input_fastq2
  String bwa_commandline
  String output_bam_basename
  File ref_fasta
  File ref_fasta_index
  File ref_dict

  # This is the .alt file from bwa-kit (https://github.com/lh3/bwa/tree/master/bwakit),
  # listing the reference contigs that are "alternative".
  File ref_alt

  File ref_amb
  File ref_ann
  File ref_bwt
  File ref_pac
  File ref_sa
  Float disk_size
  Int compression_level
  Int preemptible_tries

  command <<<
    set -o pipefail
    set -e

    # set the bash variable needed for the command-line
    bash_ref_fasta=${ref_fasta}
    /usr/gitc/${bwa_commandline} ${input_fastq1} ${input_fastq2} | samtools view -bS -o ${output_bam_basename}.bam
    
  >>>
  runtime {
    preemptible: preemptible_tries
    memory: "40 GB"
    cpu: "16"
    disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.1-1540490856"
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
  }
}

task SortSam {
  File input_bam
  String output_bam_basename
  Int preemptible_tries
  Int compression_level
  Float disk_size

  command {
    java -Dsamjdk.compression_level=${compression_level} -Xms4000m -jar /usr/gitc/picard.jar \
      SortSam \
      INPUT=${input_bam} \
      OUTPUT=${output_bam_basename}.bam \
      SORT_ORDER="coordinate" \
      CREATE_INDEX=true \
      CREATE_MD5_FILE=true \
      MAX_RECORDS_IN_RAM=300000
  }
  runtime {
    disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
    cpu: "1"
    memory: "5000 MB"
    preemptible: preemptible_tries
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.1-1540490856"
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
    File output_bam_index = "${output_bam_basename}.bai"
    File output_bam_md5 = "${output_bam_basename}.bam.md5"
  }
}

# Mark duplicate reads to avoid counting non-independent observations
task MarkDuplicates {
  File input_bam
  String output_bam_basename
  String metrics_filename
  Float disk_size
  Int preemptible_tries
  Int compression_level

  # The program default for READ_NAME_REGEX is appropriate in nearly every case.
  # Sometimes we wish to supply "null" in order to turn off optical duplicate detection
  # This can be desirable if you don't mind the estimated library size being wrong and optical duplicate detection is taking >7 days and failing
  String? read_name_regex

 # Task is assuming query-sorted input so that the Secondary and Supplementary reads get marked correctly
 # This works because the output of BWA is query-grouped and therefore, so is the output of MergeBamAlignment.
 # While query-grouped isn't actually query-sorted, it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname"
  command {
    java -Dsamjdk.compression_level=${compression_level} -Xms4000m -jar /usr/gitc/picard.jar \
      MarkDuplicates \
      INPUT=${input_bam} \
      OUTPUT=${output_bam_basename}.bam \
      METRICS_FILE=${metrics_filename} \
      VALIDATION_STRINGENCY=SILENT \
      OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
      ASSUME_SORT_ORDER="coordinate" \
      CLEAR_DT="false" \
      ADD_PG_TAG_TO_READS=false
    if [[ -s ${output_bam_basename}.bam ]]; then
      samtools index ${output_bam_basename}.bam
    fi
  }
  runtime {
    preemptible: preemptible_tries
    cpu: "1"
    memory: "7 GB"
    disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.1-1540490856"
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
    File output_bam_index = "${output_bam_basename}.bam.bai"
    File duplicate_metrics = "${metrics_filename}"
  }
}

# Generate sets of intervals for scatter-gathering over chromosomes
task CreateSequenceGroupingTSV {
  File ref_dict
  Int preemptible_tries
  
  # Use python to create the Sequencing Groupings used for BQSR and PrintReads Scatter.
  # It outputs to stdout where it is parsed into a wdl Array[Array[String]]
  # e.g. [["1"], ["2"], ["3", "4"], ["5"], ["6", "7", "8"]]
  command <<<
    python <<CODE
    with open("${ref_dict}", "r") as ref_dict_file:
        sequence_tuple_list = []
        longest_sequence = 0
        for line in ref_dict_file:
            if line.startswith("@SQ"):
                line_split = line.split("\t")
                # (Sequence_Name, Sequence_Length)
                sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))
        longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]
    # We are adding this to the intervals because hg38 has contigs named with embedded colons and a bug in GATK strips off
    # the last element after a :, so we add this as a sacrificial element.
    hg38_protection_tag = ":1+"
    # initialize the tsv string with the first sequence
    tsv_string = sequence_tuple_list[0][0] + hg38_protection_tag
    temp_size = sequence_tuple_list[0][1]
    for sequence_tuple in sequence_tuple_list[1:]:
        if temp_size + sequence_tuple[1] <= longest_sequence:
            temp_size += sequence_tuple[1]
            tsv_string += "\t" + sequence_tuple[0] + hg38_protection_tag
        else:
            tsv_string += "\n" + sequence_tuple[0] + hg38_protection_tag
            temp_size = sequence_tuple[1]
    # add the unmapped sequences as a separate line to ensure that they are recalibrated as well
    with open("sequence_grouping.txt","w") as tsv_file:
      tsv_file.write(tsv_string)
      tsv_file.close()

    tsv_string += '\n' + "unmapped"

    with open("sequence_grouping_with_unmapped.txt","w") as tsv_file_with_unmapped:
      tsv_file_with_unmapped.write(tsv_string)
      tsv_file_with_unmapped.close()
    CODE
  >>>
  runtime {
    preemptible: preemptible_tries
    docker: "python:2.7"
    memory: "2 GB"
  }
  output {
    Array[Array[String]] sequence_grouping = read_tsv("sequence_grouping.txt")
    Array[Array[String]] sequence_grouping_with_unmapped = read_tsv("sequence_grouping_with_unmapped.txt")
  }
}

# Generate Base Quality Score Recalibration (BQSR) model
task BaseRecalibrator {
  File input_bam
  File input_bam_index
  String recalibration_report_filename
  Array[String] sequence_group_interval
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
      --known-sites ${sep=" --known-sites " known_indels_sites_VCFs} \
      -L ${sep=" -L " sequence_group_interval}
  }
  runtime {
    preemptible: preemptible_tries
    cpu: "1"
    memory: "6 GB"
    disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
  }
  output {
    File recalibration_report = "${recalibration_report_filename}"
  }
}

# Combine multiple recalibration tables from scattered BaseRecalibrator runs
task GatherBqsrReports {
  Array[File] input_bqsr_reports
  String output_report_filename
  Int disk_size
  Int preemptible_tries

  command {
    /gatk/gatk --java-options "-Xms3000m" \
      GatherBQSRReports \
      -I ${sep=' -I ' input_bqsr_reports} \
      -O ${output_report_filename}
    }
  runtime {
    preemptible: preemptible_tries
    memory: "3500 MB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_bqsr_report = "${output_report_filename}"
  }
}

# Apply Base Quality Score Recalibration (BQSR) model
task ApplyBQSR {
  File input_bam
  File input_bam_index
  String output_bam_basename
  Array[String] sequence_group_interval
  File recalibration_report
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  Float disk_size
  Int compression_level
  Int preemptible_tries 
  command {
    /gatk/gatk --java-options "-XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps \
      -XX:+PrintGCDetails -Xloggc:gc_log.log \
      -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Dsamjdk.compression_level=${compression_level} -Xms3000m" \
      ApplyBQSR \
      --create-output-bam-md5 \
      --add-output-sam-program-record \
      -R ${ref_fasta} \
      -I ${input_bam} \
      --use-original-qualities \
      -O ${output_bam_basename}.bam \
      -bqsr ${recalibration_report} \
      --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
      -L ${sep=" -L " sequence_group_interval}
  }
  runtime {
    preemptible: preemptible_tries
    cpu: "1"
    memory: "3500 MB"
    disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
  }
  output {
    File recalibrated_bam = "${output_bam_basename}.bam"
    File recalibrated_bam_index = "${output_bam_basename}.bai"
    File recalibrated_bam_checksum = "${output_bam_basename}.bam.md5"
  }
}


# Combine multiple recalibrated BAM files from scattered ApplyRecalibration runs
task GatherBamFiles {
  Array[File] input_bams
  String output_bam_basename
  Float disk_size
  Int compression_level
  Int preemptible_tries

  command {
    java -Dsamjdk.compression_level=${compression_level} -Xms2000m -jar /usr/gitc/picard.jar \
      GatherBamFiles \
      INPUT=${sep=' INPUT=' input_bams} \
      OUTPUT=${output_bam_basename}.bam \
      CREATE_INDEX=true \
      CREATE_MD5_FILE=true
    }
  runtime {
    preemptible: preemptible_tries
    memory: "3 GB"
    disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.1-1540490856"
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
    File output_bam_index = "${output_bam_basename}.bai"
    File output_bam_md5 = "${output_bam_basename}.bam.md5"
  }
}

# Notes on the contamination estimate:
# The contamination value is read from the FREEMIX field of the selfSM file output by verifyBamId
#
# In Zamboni production, this value is stored directly in METRICS.AGGREGATION_CONTAM
#
# Contamination is also stored in GVCF_CALLING and thereby passed to HAPLOTYPE_CALLER
# But first, it is divided by an underestimation factor thusly:
#   float(FREEMIX) / ContaminationUnderestimationFactor
#     where the denominator is hardcoded in Zamboni:
#     val ContaminationUnderestimationFactor = 0.75f
#
# Here, I am handling this by returning both the original selfSM file for reporting, and the adjusted
# contamination estimate for use in variant calling
task CheckContamination {
  File input_bam
  File input_bam_index
  File contamination_sites_ud
  File contamination_sites_bed
  File contamination_sites_mu
  File ref_fasta
  File ref_fasta_index
  String output_prefix
  Float disk_size
  Int preemptible_tries
  Float contamination_underestimation_factor

  command <<<
    set -e

    # creates a ${output_prefix}.selfSM file, a TSV file with 2 rows, 19 columns.
    # First row are the keys (e.g., SEQ_SM, RG, FREEMIX), second row are the associated values
    /usr/gitc/VerifyBamID \
    --Verbose \
    --NumPC 4 \
    --Output ${output_prefix} \
    --BamFile ${input_bam} \
    --Reference ${ref_fasta} \
    --UDPath ${contamination_sites_ud} \
    --MeanPath ${contamination_sites_mu} \
    --BedPath ${contamination_sites_bed} \
    1>/dev/null

    # used to read from the selfSM file and calculate contamination, which gets printed out
    python3 <<CODE
    import csv
    import sys
    with open('${output_prefix}.selfSM') as selfSM:
      reader = csv.DictReader(selfSM, delimiter='\t')
      i = 0
      for row in reader:
        if float(row["FREELK0"])==0 and float(row["FREELK1"])==0:
          # a zero value for the likelihoods implies no data. This usually indicates a problem rather than a real event.
          # if the bam isn't really empty, this is probably due to the use of a incompatible reference build between
          # vcf and bam.
          sys.stderr.write("Found zero likelihoods. Bam is either very-very shallow, or aligned to the wrong reference (relative to the vcf).")
          sys.exit(1)
        print(float(row["FREEMIX"])/${contamination_underestimation_factor})
        i = i + 1
        # there should be exactly one row, and if this isn't the case the format of the output is unexpectedly different
        # and the results are not reliable.
        if i != 1:
          sys.stderr.write("Found %d rows in .selfSM file. Was expecting exactly 1. This is an error"%(i))
          sys.exit(2)
    CODE
  >>>
  runtime {
    preemptible: preemptible_tries
    cpu: "1"
    memory: "2 GB"
    disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
    docker: "us.gcr.io/broad-gotc-prod/verify-bam-id:c8a66425c312e5f8be46ab0c41f8d7a1942b6e16-1500298351"
  }
  output {
    File selfSM = "${output_prefix}.selfSM"
    Float contamination = read_float(stdout())
  }
}

# This task calls picard's IntervalListTools to scatter the input interval list into scatter_count sub interval lists
# Note that the number of sub interval lists may not be exactly equal to scatter_count.  There may be slightly more or less.
# Thus we have the block of python to count the number of generated sub interval lists.
task ScatterIntervalList {
  File interval_list
  Int scatter_count
  Int break_bands_at_multiples_of

  command <<<
    set -e
    mkdir out
    java -Xms1g -jar /usr/gitc/picard.jar \
      IntervalListTools \
      SCATTER_COUNT=${scatter_count} \
      SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
      UNIQUE=true \
      SORT=true \
      BREAK_BANDS_AT_MULTIPLES_OF=${break_bands_at_multiples_of} \
      INPUT=${interval_list} \
      OUTPUT=out

    python3 <<CODE
    import glob, os
    # Works around a JES limitation where multiples files with the same name overwrite each other when globbed
    intervals = sorted(glob.glob("out/*/*.interval_list"))
    for i, interval in enumerate(intervals):
      (directory, filename) = os.path.split(interval)
      newName = os.path.join(directory, str(i + 1) + filename)
      os.rename(interval, newName)
    print(len(intervals))
    CODE
  >>>
  output {
    Array[File] out = glob("out/*/*.interval_list")
    Int interval_count = read_int(stdout())
  }
  runtime {
    memory: "2 GB"
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.1-1540490856"
  }
}

# Call variants on a single sample with HaplotypeCaller to produce a GVCF
task HaplotypeCaller {
  File input_bam
  File input_bai
  File interval_list
  String gvcf_basename
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  Float? contamination
  Float disk_size
  Int preemptible_tries

  # We use interval_padding 500 below to make sure that the HaplotypeCaller has context on both sides around
  # the interval because the assembly uses them.
  #
  # Using PrintReads is a temporary solution until we update HaploypeCaller to use GATK4. Once that is done,
  # HaplotypeCaller can stream the required intervals directly from the cloud.
  command {
    /gatk/gatk --java-options "-Xms2g" \
      PrintReads \
      -I ${input_bam} \
      --interval-padding 500 \
      -L ${interval_list} \
      -O local.sharded.bam \
    && \
    /gatk/gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms8000m" \
      HaplotypeCaller \
      -R ${ref_fasta} \
      -O ${gvcf_basename}.vcf.gz \
      -I local.sharded.bam \
      -L ${interval_list} \
      -ERC GVCF \
      --max-alternate-alleles 3 \
      --contamination ${default=0 contamination} \
      --read-filter OverclippedReadFilter
  }
  runtime {
    preemptible: preemptible_tries
    memory: "40 GB"
    cpu: "1"
    disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
  }
  output {
    File output_gvcf = "${gvcf_basename}.vcf.gz"
    File output_gvcf_index = "${gvcf_basename}.vcf.gz.tbi"
  }
}


# Combine multiple VCFs or GVCFs from scattered HaplotypeCaller runs
task MergeVCFs {
  Array[File] input_vcfs
  Array[File] input_vcfs_indexes
  String output_vcf_name
  Int disk_size
  Int preemptible_tries

  # Using MergeVcfs instead of GatherVcfs so we can create indices
  # See https://github.com/broadinstitute/picard/issues/789 for relevant GatherVcfs ticket
  command {
    java -Xms2000m -jar /usr/gitc/picard.jar \
      MergeVcfs \
      INPUT=${sep=' INPUT=' input_vcfs} \
      OUTPUT=${output_vcf_name}
  }
  runtime {
    preemptible: preemptible_tries
    memory: "3 GB"
    disks: "local-disk " + disk_size + " HDD"
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.1-1540490856"
  }
  output {
    File output_vcf = "${output_vcf_name}"
    File output_vcf_index = "${output_vcf_name}.tbi"
  }
}

# Validate a GVCF with -gvcf specific validation
task ValidateGVCF {
  File input_vcf
  File input_vcf_index
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File dbSNP_vcf
  File dbSNP_vcf_index
  File wgs_calling_interval_list
  Float disk_size
  Int preemptible_tries

  command {
    /gatk/gatk --java-options "-Xms3000m" \
      ValidateVariants \
      -V ${input_vcf} \
      -R ${ref_fasta} \
      -L ${wgs_calling_interval_list} \
      -gvcf \
      --validation-type-to-exclude ALLELES \
      --dbsnp ${dbSNP_vcf}
  }
  runtime {
    preemptible: preemptible_tries
    cpu: "1"
    memory: "3500 MB"
    disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
  }
}

# Collect variant calling metrics from GVCF output
task CollectGvcfCallingMetrics {
  File input_vcf
  File input_vcf_index
  String metrics_basename
  File dbSNP_vcf
  File dbSNP_vcf_index
  File ref_dict
  File wgs_evaluation_interval_list
  Float disk_size
  Int preemptible_tries

  command {
    java -Xms2000m -jar /usr/gitc/picard.jar \
      CollectVariantCallingMetrics \
      INPUT=${input_vcf} \
      OUTPUT=${metrics_basename} \
      DBSNP=${dbSNP_vcf} \
      SEQUENCE_DICTIONARY=${ref_dict} \
      TARGET_INTERVALS=${wgs_evaluation_interval_list} \
      GVCF_INPUT=true
  }
  runtime {
    preemptible: preemptible_tries
    cpu: "1"
    memory: "3 GB"
    disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.1-1540490856"
  }
  output {
    File summary_metrics = "${metrics_basename}.variant_calling_summary_metrics"
    File detail_metrics = "${metrics_basename}.variant_calling_detail_metrics"
  }
}
