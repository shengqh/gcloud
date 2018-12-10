## Quanhu Sheng, quanhu.sheng.1@vanderbilt.edu 08/2018
## Shenglai He, shenglai.he@vanderbilt.edu 07/2018
##
## This WDL pipeline implements data pre-processing and initial variant calling (GVCF
## generation) according to the GATK Best Practices (June 2016) for germline SNP and
## Indel discovery in human whole-genome sequencing (WGS) data.
##
## Inputs file format: fastq 

# WORKFLOW DEFINITION
workflow paired_dnaseq {

  File contamination_sites_ud
  File contamination_sites_bed
  File contamination_sites_mu
  File? fingerprint_genotypes_file
  File? haplotype_database_file
  File wgs_evaluation_interval_list
  File wgs_coverage_interval_list

  String sample_name
  String unmapped_fastq_suffix
  
  File raw_fastq1
  File raw_fastq2

  File wgs_calling_interval_list
  Int haplotype_scatter_count
  Int break_bands_at_multiples_of
  Int? read_length

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

  # ValidateSamFile runs out of memory in mate validation on crazy edge case data, so we want to skip the mate validation
  # in those cases.  These values set the thresholds for what is considered outside the normal realm of "reasonable" data.
  Float max_duplication_in_reasonable_sample = 0.30
  Float max_chimerism_in_reasonable_sample = 0.15

  String bwa_commandline

  String recalibrated_bam_basename = sample_name + ".aligned.duplicates_marked.recalibrated"

  Int compression_level = 2

  # Get the version of BWA to include in the PG record in the header of the BAM produced
  # by MergeBamAlignment.
  #call GetBwaVersion

  # Get the size of the standard reference files as well as the additional reference files needed for BWA
  Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB")
  Float bwa_ref_size = ref_size + size(ref_alt, "GB") + size(ref_amb, "GB") + size(ref_ann, "GB") + size(ref_bwt, "GB") + size(ref_pac, "GB") + size(ref_sa, "GB")
  Float dbsnp_size = size(dbSNP_vcf, "GB")

  # Align flowcell-level unmapped input bams in parallel

  Float unmapped_fastq_size = size(raw_fastq1, "GB")+size(raw_fastq2, "GB")


  String sub_strip_path = "gs://.*/"
  String sub_strip_unmapped = unmapped_fastq_suffix + "$"
    
  # Map reads to reference
  call FastqBwaMem {
    input:
      input_fastq1 = raw_fastq1,
      input_fastq2 = raw_fastq2,
      bwa_commandline = bwa_commandline,
      sample_name = sample_name,
      output_bam_basename = sample_name + ".aligned.unsorted",
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

  # Outputs that will be retained when execution is complete
  output {
    File output_sam = FastqBwaMem.output_sam
  }
}

# TASK DEFINITIONS

# Get version of BWA
task GetBwaVersion {
  command {
    # not setting set -o pipefail here because /bwa has a rc=1 and we dont want to allow rc=1 to succeed because
    # the sed may also fail with that error and that is something we actually want to fail on.
    /usr/gitc/bwa 2>&1 | \
    grep -e '^Version' | \
    sed 's/Version: //'
  }
  runtime {
    memory: "1 GB"
    docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
  }
  output {
    String version = read_string(stdout())
  }
}

# Read unmapped BAM, convert on-the-fly to FASTQ and stream to BWA MEM for alignment, then stream to MergeBamAlignment
task FastqBwaMem {
  File input_fastq1
  File input_fastq2
  String bwa_commandline
  String sample_name
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
    bash_sample_name=${sample_name}
    # if ref_alt has data in it,
    if [ -s ${ref_alt} ]; then
      
      /usr/gitc/${bwa_commandline} -R '@RG\tID:'$bash_sample_name'\tPU:illumina\tLB:'$bash_sample_name'\tSM:'$bash_sample_name'\tPL:illumina' -Y $bash_ref_fasta ${input_fastq1} ${input_fastq2} >${output_bam_basename}.sam \
      2> >(tee ${output_bam_basename}.bwa.stderr.log >&2)
      grep -m1 "read .* ALT contigs" ${output_bam_basename}.bwa.stderr.log | \
      grep -v "read 0 ALT contigs"

    # else ref_alt is empty or could not be found
    else
      exit 1;
    fi
  >>>
  runtime {
    preemptible: preemptible_tries
    memory: "14 GB"
    cpu: "16"
    disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
    docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
  }
  output {
    File output_sam = "${output_bam_basename}.sam"
    File bwa_stderr_log = "${output_bam_basename}.bwa.stderr.log"
  }
}

