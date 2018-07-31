## Shenglai He, shenglai.he@vanderbilt.edu 07/2018
##
## This WDL pipeline implements aligning of rna seq data using STAR and
## counts reads to genomic features 
##
## Inputs file format: fastq

# WORKFLOW DEFINITION
workflow rnaaligncounts {

  String sample_name
  String base_file_name
  String final_gvcf_base_name
  #Array[File] flowcell_unmapped_bams
  # String unmapped_bam_suffix
  String unmapped_fastq_suffix

  
  File raw_fastq1
  File raw_fastq2

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
  File chrLength
  File chrNameLength
  File chrName
  File chrStart
  File Genome
  File genomeParameters
  File SA
  File SAindex
  File exonGeTrInfo_tab
  File exonInfo_tab
  File geneInfo_tab
  File human_g1k_v37_dict
  File human_g1k_v37_fasta
  File human_g1k_v37_fasta_fai
  File human_g1k_v37_len
  File Log
  File sjdbInfo
  File sjdbList_fromGTF_out_tab
  File sjdbList_out_tab
  File transcriptInfo_tab
  File ref_gtf
  String ref_dir

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

  #String bwa_commandline

  String recalibrated_bam_basename = base_file_name + ".aligned.duplicates_marked.recalibrated"

  Int compression_level = 2


  # Get the size of the standard reference files as well as the additional reference files needed for STAR
  #Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB")
  #Float bwa_ref_size = ref_size + size(ref_alt, "GB") + size(ref_amb, "GB") + size(ref_ann, "GB") + size(ref_bwt, "GB") + size(ref_pac, "GB") + size(ref_sa, "GB")
  Float star_ref_size = 100.0 
  Float dbsnp_size = size(dbSNP_vcf, "GB")

  # Align flowcell-level unmapped input bams in parallel

  Float unmapped_fastq_size = size(raw_fastq1, "GB")+size(raw_fastq2, "GB")


    String sub_strip_path = "gs://.*/"
    String sub_strip_unmapped = unmapped_fastq_suffix + "$"
    String sub_sub = sub(sub(raw_fastq1, sub_strip_path, ""), sub_strip_unmapped, "")

    
    # Map reads to reference
    call FastqStar {
      input:
        input_fastq1 = raw_fastq1,
        input_fastq2 = raw_fastq2,
        sample_name = sample_name + "_",
        #output_bam_basename = sub_sub+".",
	chrLength = chrLength,
	chrNameLength = chrNameLength,
        chrName = chrName,
	chrStart = chrStart,
	Genome = Genome,
	genomeParameters = genomeParameters,
        SA = SA,
        SAindex = SAindex,
	exonGeTrInfo_tab = exonGeTrInfo_tab,
	exonInfo_tab = exonInfo_tab,
	geneInfo_tab = geneInfo_tab,
	human_g1k_v37_dict = human_g1k_v37_dict,
	human_g1k_v37_fasta = human_g1k_v37_fasta,
        human_g1k_v37_fasta_fai =  human_g1k_v37_fasta_fai,
        human_g1k_v37_len =  human_g1k_v37_len,
        Log = Log,
        sjdbInfo = sjdbInfo,
	sjdbList_fromGTF_out_tab = sjdbList_fromGTF_out_tab,
	sjdbList_out_tab =  sjdbList_out_tab,
        transcriptInfo_tab = transcriptInfo_tab,
	ref_dir = "/cromwell_root"+ref_dir,
        disk_size = unmapped_fastq_size + star_ref_size + (bwa_disk_multiplier * unmapped_fastq_size) + additional_disk,
        compression_level = compression_level,
        preemptible_tries = preemptible_tries
    }

    Float unsorted_bam_size = size(FastqStar.unsorted_bam, "GB")
    
    call Counts {
      input:
        ref_gtf = ref_gtf,
        input_bam = FastqStar.unsorted_bam,
        sample_name = sample_name,
        disk_size = unsorted_bam_size + additional_disk,
        preemptible_tries = preemptible_tries,
    }
    

  # Outputs that will be retained when execution is complete
  output {

   #  File unsorted_bam = FastqStar.unsorted_bam
     File output_counts = Counts.output_counts
     File sortedByCoordinate_bam =  FastqStar.sortedByCoordinate_bam
  }
}



task FastqStar {
  File input_fastq1
  File input_fastq2
  String sample_name
#  String output_bam_basename
  File chrLength
  File chrNameLength
  File chrName
  File chrStart
  File Genome
  File genomeParameters
  File SA
  File SAindex
  File exonGeTrInfo_tab
  File exonInfo_tab
  File geneInfo_tab
  File human_g1k_v37_dict
  File human_g1k_v37_fasta
  File human_g1k_v37_fasta_fai
  File human_g1k_v37_len
  File Log
  File sjdbInfo
  File sjdbList_fromGTF_out_tab
  File sjdbList_out_tab
  File transcriptInfo_tab
  String ref_dir 
		
 # File ref_fasta
 # File ref_fasta_index
 # File ref_dict

  # This is the .alt file from bwa-kit (https://github.com/lh3/bwa/tree/master/bwakit),
  # listing the reference contigs that are "alternative".
 # File ref_alt

 # File ref_amb
 # File ref_ann
 # File ref_bwt
 # File ref_pac
 # File ref_sa
  Float disk_size
  Int compression_level
  Int preemptible_tries

  command <<<
    set -o pipefail
    set -e

    # set the bash variable needed for the command-line
    # if ref_alt has data in it
    /STAR/bin/Linux_x86_64/STAR --twopassMode Basic \
    --outSAMprimaryFlag AllBestScore \
    --outSAMattrRGline ID:${sample_name} SM:${sample_name} LB:${sample_name} PL:ILLUMINA PU:ILLUMINA \
    --runThreadN 8 \
    --genomeDir ${ref_dir} \
    --readFilesIn  ${input_fastq1} ${input_fastq2} \
    --readFilesCommand zcat \
    --outFileNamePrefix ${sample_name} \
    --outSAMtype BAM SortedByCoordinate Unsorted

    # else ref_alt is empty or could not be found
    #else
    #  exit 1;
    #fi
  >>>
  runtime {
    preemptible: preemptible_tries
    memory: "40 GB"
    cpu: "16"
    disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
    docker: "gcr.io/gatk-test-206214/cqs_base:latest"
  }
  output {
    File  unsorted_bam = "${sample_name}Aligned.out.bam"
    File sortedByCoordinate_bam = "${sample_name}Aligned.sortedByCoord.out.bam"
   # File final_out = "${output_bam_basename}Log.final.out"
   # File log_out = "${output_bam_basename}Log.out"
   # File Log_progress_out = "${output_bam_basename}Log.progress.out"
   # File SJ_out_tab = "${output_bam_basename}SJ.out.tab"
   # File star_stderr_log = "${output_bam_basename}star.stderr.log"
  }
}

task Counts {
  File ref_gtf
  File input_bam
  String sample_name
  Float disk_size
  Int preemptible_tries

  command <<<
    set -o pipefail
    set -e

    /bin/featureCounts -g gene_id -t exon -p -T 8 -a ${ref_gtf} -o ${sample_name}.count ${input_bam}
     
  >>>
  runtime {
      preemptible: preemptible_tries
      memory: "8 GB"
      cpu: "8"
      disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
      docker: "gcr.io/gatk-test-206214/cqs_base:latest"
  }
  output{

   File output_counts = "${sample_name}.count"
   }
 }


