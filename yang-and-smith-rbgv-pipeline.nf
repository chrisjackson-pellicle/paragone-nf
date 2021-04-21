#!/usr/bin/env nextflow

///////////////////////////////////////////////////////////////////////////////////////////////////////////
////////// Nextflow Pipeline for processing HybPiper paralogs fasta files with Y&S pipeline  //////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

def helpMessage() {
    log.info """

    Usage:
    The typical command for running the pipeline is as follows:

    nextflow run yang-and-smith-rbgv-pipeline.nf \
    -c yang-and-smith-rbgv.config \
    --hybpiper_paralogs_directory <directory> \
    --external_outgroups_file <file>
    --external_outgroups <taxon1,taxon2,taxon3...>
    --internal_outgroups <taxon1,taxon2,taxon3...>
    -profile <profile>
    -resume
 
    Mandatory arguments:

      ################################################################################################################

      --hybpiper_paralogs_directory <directory>       Path to folder containing HybPiper paralog fasta files.

      ..and either

      --internal_outgroups <taxon1,taxon2,taxon3...>  A comma-separated list of taxa present in the paralog fasta files
                                                      to use as outgroups. Default is none

      ...or

      --external_outgroups_file <file>                File containing fasta sequences of outgroup sequences for each
                                                      gene

      ################################################################################################################

    Optional arguments:
    
      --internal_outgroups <taxon1,taxon2,taxon3...>  A comma-separated list of taxa present in the paralog fasta files
                                                      to use as outgroups. Default is none
      --external_outgroups_file <file>                File containing fasta sequences of outgroup sequences for each
                                                      gene
      --external_outgroups <taxon1,taxon2,taxon3...>  A comma-separated list of outgroup taxa to add from the
                                                      outgroups_file. Default is all
      -profile <profile>                              Configuration profile to use. Can use multiple (comma separated)
                                                      Available: standard (default), slurm
      --output_dir                                    Specify the name of the main output results folder.
                                                      Default is 'results'
      --pool <int>                                    Number of threads for the Python multiprocessing pool. Used in
                                                      e.g. alignments and tree-building steps. Default is 1, so
                                                      e.g. one alignment will be run at a time during alignment steps.
      --threads <int>                                 Number of threads per multiprocessing pool instance. Used for
                                                      programs that support multi-threading (e.g. mafft, IQTree).
                                                      Default is 1
      --no_supercontigs                               Use this flag if you are processing paralogs from a run of
                                                      HybPiper that used the --nosupercontigs flag. Mafft alignments
                                                      are re-aligned using clustal omega, which can do a better job in
                                                      these cases. Default is off
      --process_04_trim_tips_relative_cutoff <float>  When pruning long tips during the tree QC stage, provide a branch
                                                      length for the maximum imbalance between sister tips allowed.
                                                      Default is 0.2
      --process_04_trim_tips_absolute_cutoff <float>  When pruning long tips during the tree QC stage, provide a branch
                                                      length for the maximum allowed tip branch length. Default is 0.4.
      --process_06_branch_length_cutoff <float>       When pruning long internal branches (putative deep paralogs)
                                                      during the tree QC stage, provide a branch length for the maximum
                                                      allowed internal branch length. Default is 0.3
      --process_06_minimum_taxa <int>                 After the final tree-pruning step prior to paralogy resolution,
                                                      only retain trees with a minimum number of taxa remaining.
                                                      Default is 3
      --process_09_prune_paralog_MO_minimum_taxa <int>
                                                      For the MO method, only process trees with a minimum number of
                                                      taxa. Default is 2
      --process_10_prune_paralogs_RT_minimum_ingroup_taxa <int>
                                                      For the RT method, only process trees with a minumum number of
                                                      ingroup taxa. Default is 2
      --process_11_prune_paralogs_MI_relative_tip_cutoff <float>
                                                      Default is 0.2
      --process_11_prune_paralogs_MI_absolute_tip_cutoff <float>
                                                      Default is 0.4
      --process_11_prune_paralogs_MI_minimum_taxa <int>    
                                                      Default is 2


    """.stripIndent()
}

/*
Check that the input directory is provided, and either (at least) a file of external outgroups or list of internal
taxa to use:
*/
if (params.help || !params.hybpiper_paralogs_directory || (!params.external_outgroups_file && !params.internal_outgroups)) {
  helpMessage()
  exit 0
}


// Check that only allowed parameters are provided:
allowed_params = ["internal_outgroups", "external_outgroups_file", "hybpiper_paralogs_directory",
"external_outgroups", "threads", "pool", "nosupercontigs", "process_04_trim_tips_relative_cutoff",
"process_04_trim_tips_absolute_cutoff", "process_06_branch_length_cutoff", "process_06_minimum_taxa",
"process_09_prune_paralog_MO_minimum_taxa", "process_10_prune_paralogs_RT_minimum_ingroup_taxa",
"process_11_prune_paralogs_MI_relative_tip_cutoff", "process_11_prune_paralogs_MI_absolute_tip_cutoff",
"process_11_prune_paralogs_MI_minimum_taxa", "outdir", "help", "no_supercontigs"]

params.each { entry ->
  if (! allowed_params.contains(entry.key)) {
      println("The parameter <${entry.key}> is not known");
      exit 0;
  }
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////  Set up channels for input data  ///////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

// External outgroup sequences file:
Channel
  .fromPath("${params.external_outgroups_file}")
  .first()
  .set { outgroups_file_ch }

// HybPiper paralog fasta files:
Channel
  .fromPath("${params.hybpiper_paralogs_directory}", checkIfExists: true)
  .into { paralogs_ch_1; paralogs_ch_2 }


//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////// Run 01_check_outgroups_align_and_hmmclean.py  ///////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

process align_paralogs_01 {
  // echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
  file(external_outgroups_file) from outgroups_file_ch
  file(paralog_folder) from paralogs_ch_1

  output:
  stdout outgroup_coverage_ch
  file("outgroup_coverage_report.tsv")
  file("01_alignments")
  file("02_alignments_hmmcleaned") into (alignments_hmmcleaned_ch_1,
                                         alignments_hmmcleaned_ch_2,
                                         alignments_hmmcleaned_ch_3,
                                         alignments_hmmcleaned_ch_4)

  script:
  if (params.external_outgroups_file) {
    external_outgroups_file_string = "-external_outgroups_file ${external_outgroups_file}"
  } else {
    external_outgroups_file_string = ''
  }

  if (params.external_outgroups) {
    external_outgroups_list = params.external_outgroups?.tokenize(',')
    external_outgroups_string = ''

    for (outgroup in external_outgroups_list) {
       external_outgroup_string = "-external_outgroup ${outgroup} "
       external_outgroups_string = external_outgroups_string + external_outgroup_string
    }
  } else {
    external_outgroups_string = ''
  }

  if (params.internal_outgroups) {
    internal_outgroups_list = params.internal_outgroups?.tokenize(',')
    internal_outgroups_string = ''

    for (outgroup in internal_outgroups_list) {
      internal_outgroup_string = "-internal_outgroup ${outgroup} "
      internal_outgroups_string = internal_outgroups_string + internal_outgroup_string
    }
  } else {
    internal_outgroups_string = ''
  }

  if (!params.no_supercontigs) {
  """
  # python /Users/chrisjackson/PycharmProjects/Yang-and-Smith-RBGV-scripts/01_check_outgroups_align_and_hmmclean.py 
  python /Yang-and-Smith-RBGV-scripts/01_check_outgroups_align_and_hmmclean.py \
   ${paralog_folder} \
  ${external_outgroups_file_string} \
  ${external_outgroups_string} \
  ${internal_outgroups_string} \
  -pool ${params.pool} \
  -threads ${params.threads} \
  """
  } else if (params.no_supercontigs) {
  """
  # python /Users/chrisjackson/PycharmProjects/Yang-and-Smith-RBGV-scripts/01_check_outgroups_align_and_hmmclean.py 
  python /Yang-and-Smith-RBGV-scripts//01_check_outgroups_align_and_hmmclean.py \
  ${paralog_folder} \
  ${external_outgroups_file_string} \
  ${external_outgroups_string} \
  ${internal_outgroups_string} \
  -pool ${params.pool} \
  -threads ${params.threads} \
  -no_supercontigs
  """
  }
}

// Log outgroup coverage statistics to screen and `nextflow.log` file:    
outgroup_coverage_ch.subscribe { log.info '\n' + "$it" + '\n' }


//////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Run 03_alignment_to_tree.py  //////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

process alignment_to_tree_03 {
  //echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
  file(alignments_folder) from alignments_hmmcleaned_ch_1

  output:
  file("03_tree_files") into tree_folder_ch

  script:
  """ 
  # python /Users/chrisjackson/PycharmProjects/Yang-and-Smith-RBGV-scripts/03_alignment_to_tree.py
  python /Yang-and-Smith-RBGV-scripts/03_alignment_to_tree.py \
  ${alignments_folder} \
  -threads_pool ${params.pool} \
  -threads_iqtree ${params.threads}
  """
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Run 04_trim_tips.py  //////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

process trim_tips_04 {
  // echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
  file(trees_folder) from tree_folder_ch

  output:
  file("04_trim_tips") into trim_tips_ch

  script:
  """  
  # python /Users/chrisjackson/PycharmProjects/Yang-and-Smith-RBGV-scripts/04_trim_tips.py
  python /Yang-and-Smith-RBGV-scripts/04_trim_tips.py \
  ${trees_folder} \
  .treefile \
  ${params.process_04_trim_tips_relative_cutoff} \
  ${params.process_04_trim_tips_absolute_cutoff} \
  04_trim_tips
  echo "Finished running trim_tips_04"
  """
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////  Run 05_mask_tips_by_taxonID_transcripts.py  /////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

process mask_tips_05 {
  // echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
  file(trimmed_tips_folder) from trim_tips_ch
  file(alignments_folder) from alignments_hmmcleaned_ch_2

  output:
  file("05_masked_tips") into masked_tips_ch 

  script:
  """
  # python /Users/chrisjackson/PycharmProjects/Yang-and-Smith-RBGV-scripts/05_mask_tips_by_taxonID_transcripts.py
  python /Yang-and-Smith-RBGV-scripts/05_mask_tips_by_taxonID_transcripts.py \
  ${trimmed_tips_folder} \
  ${alignments_folder} \
  y \
  05_masked_tips
  """
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////  Run 06_cut_long_internal_branches.py  ///////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

process cut_long_internal_branches_06 {
  // echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
  file(masked_tips_folder) from masked_tips_ch

  output:
  file("06_cut_internal_branches") into cut_internal_branches_ch

  script:
  """
  mkdir 06_cut_internal_branches
  # python /Users/chrisjackson/PycharmProjects/Yang-and-Smith-RBGV-scripts/06_cut_long_internal_branches.py
  python /Yang-and-Smith-RBGV-scripts/06_cut_long_internal_branches.py \
  ${masked_tips_folder} \
  .mm \
  ${params.process_06_branch_length_cutoff} \
  ${params.process_06_minimum_taxa} \
  06_cut_internal_branches
  """
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Run 07_subset_fasta_from_tree.py  ///////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

process write_alignment_subset_07 {
  // echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'
                    
  input:
  file(cut_internal_branches_folder) from cut_internal_branches_ch
  file(alignment_folder) from alignments_hmmcleaned_ch_3

  output:
  file("07_selected_alignments") into selected_alignments_ch

  script:
  """
  mkdir 07_selected_alignments
  # python /Users/chrisjackson/PycharmProjects/Yang-and-Smith-RBGV-scripts/07_subset_fasta_from_tree.py
  python /Yang-and-Smith-RBGV-scripts/07_subset_fasta_from_tree.py \
  ${cut_internal_branches_folder} \
  .subtree \
  ${alignment_folder} \
  NotApplicable \
  07_selected_alignments
  """
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Run 08_mafft_alignment_and_iqtree.py  ///////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

process realign_and_iqtree_08 {
  //  echo  true
   label 'in_container'
   publishDir "${params.outdir}", mode: 'copy'

   input:
   file(selected_alignments_ch) from selected_alignments_ch
   file(hmm_cleaned_alignments) from alignments_hmmcleaned_ch_4
   file(external_outgroups_file) from outgroups_file_ch

   output:
   file("08_outgroups_added")
   file("09_realigned") into (realigned_with_outgroups_ch_1,realigned_with_outgroups_ch_2, realigned_with_outgroups_ch_3) 
   file("10_realigned_trees") into (realigned_trees_ch_1, realigned_trees_ch_2, realigned_trees_ch_3)
   file("in_and_outgroups_list.txt") into (in_out_list_ch_1, in_out_list_ch_2, in_out_list_ch_3)

   script:
   if (params.external_outgroups_file) {
    external_outgroups_file_string = "-external_outgroups_file ${external_outgroups_file}"
  } else {
    external_outgroups_file_string = ''
  }

   if (params.external_outgroups) {
    external_outgroups_list = params.external_outgroups?.tokenize(',')
    external_outgroups_string = ''

    for (outgroup in external_outgroups_list) {
       external_outgroup_string = "-external_outgroup ${outgroup} "
       external_outgroups_string = external_outgroups_string + external_outgroup_string
    }
  } else {
    external_outgroups_string = ''
  }

    if (params.internal_outgroups) {
    internal_outgroups_list = params.internal_outgroups?.tokenize(',')
    internal_outgroups_string = ''

    for (outgroup in internal_outgroups_list) {
      internal_outgroup_string = "-internal_outgroup ${outgroup} "
      internal_outgroups_string = internal_outgroups_string + internal_outgroup_string
    }
  } else {
    internal_outgroups_string = ''
  }
 
   if (!params.no_supercontigs) {
   """
   # python /Users/chrisjackson/PycharmProjects/Yang-and-Smith-RBGV-scripts/08_mafft_alignment_and_iqtree.py
   python /Yang-and-Smith-RBGV-scripts/08_mafft_alignment_and_iqtree.py \
   ${hmm_cleaned_alignments} \
   ${selected_alignments_ch} \
   ${external_outgroups_file_string} \
   ${external_outgroups_string} \
   ${internal_outgroups_string} \
   -threads_pool ${params.pool} \
   -threads_mafft ${params.threads}
   """
   } else {
   """
   # python /Users/chrisjackson/PycharmProjects/Yang-and-Smith-RBGV-scripts/08_mafft_alignment_and_iqtree.py
   python /Yang-and-Smith-RBGV-scripts/08_mafft_alignment_and_iqtree.py \
   ${hmm_cleaned_alignments} \
   ${selected_alignments_ch} \
   ${external_outgroups_file_string} \
   ${external_outgroups_string} \
   ${internal_outgroups_string} \
   -threads_pool ${params.pool} \
   -threads_mafft ${params.threads} \
   -no_supercontigs
   """
   }
 }


//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Run 09_prune_paralogs_MO.py  ////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

process prune_paralogs_MO_09 {
  // echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
  file(in_out_file) from in_out_list_ch_1
  file(realigned_trees_folder) from realigned_trees_ch_1

  output:
  file("11_prune_MO_trees") into MO_ch

  script:
  """
  mkdir 11_prune_MO_trees
  # python /Users/chrisjackson/PycharmProjects/Yang-and-Smith-RBGV-scripts/09_prune_paralogs_MO.py
  python /Yang-and-Smith-RBGV-scripts/09_prune_paralogs_MO.py \
  ${realigned_trees_folder} \
  .treefile \
  ${params.process_09_prune_paralog_MO_minimum_taxa} \
  11_prune_MO_trees \
  ${in_out_file}
  """
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Run 10_prune_paralogs_RT.py  ////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

process prune_paralogs_RT_10 {
  // echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
  file(in_out_file) from in_out_list_ch_2
  file(realigned_trees_folder) from realigned_trees_ch_2

  output:
  file("12_prune_RT_trees") into RT_ch

  script:
  """
  mkdir 12_prune_RT_trees
  # python /Users/chrisjackson/PycharmProjects/Yang-and-Smith-RBGV-scripts/10_prune_paralogs_RT.py
  python /Yang-and-Smith-RBGV-scripts/10_prune_paralogs_RT.py \
  ${realigned_trees_folder} \
  .treefile 12_prune_RT_trees \
  ${params.process_10_prune_paralogs_RT_minimum_ingroup_taxa} \
  ${in_out_file}
  """
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Run 11_prune_paralogs_MI.py  ////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

process prune_paralogs_MI_11 {
  // echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
  file(in_out_file) from in_out_list_ch_3
  file(realigned_trees_folder) from realigned_trees_ch_3

  output:
  file("13_prune_MI_trees") into MI_ch

  script:
  """
  mkdir 13_prune_MI_trees
  # python /Users/chrisjackson/PycharmProjects/Yang-and-Smith-RBGV-scripts/11_prune_paralogs_MI.py
  python /Yang-and-Smith-RBGV-scripts/11_prune_paralogs_MI.py \
  ${realigned_trees_folder} \
  .treefile \
  ${params.process_11_prune_paralogs_MI_relative_tip_cutoff} \
  ${params.process_11_prune_paralogs_MI_absolute_tip_cutoff} \
  ${params.process_11_prune_paralogs_MI_minimum_taxa} \
  13_prune_MI_trees
  """
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Run 07_subset_fasta_from_tree.py for MO /////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

process write_alignment_subset_MO_12 {
  // echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
  file(MO_folder) from MO_ch
  file(alignment_folder) from realigned_with_outgroups_ch_1

  output:
  file("14_selected_alignments_MO") into selected_alignments_MO_ch

  script:
  """
  mkdir 14_selected_alignments_MO
  # python /Users/chrisjackson/PycharmProjects/Yang-and-Smith-RBGV-scripts/07_subset_fasta_from_tree.py
  python /Yang-and-Smith-RBGV-scripts/07_subset_fasta_from_tree.py \
  ${MO_folder} \
  .tre \
  ${alignment_folder} \
  NotApplicable \
  14_selected_alignments_MO
  """
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Run 07_subset_fasta_from_tree.py for RT /////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

process write_alignment_subset_RT_13 {
  // echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
  file(RT_folder) from RT_ch
  file(alignment_folder) from realigned_with_outgroups_ch_2

  output:
  file("15_selected_alignments_RT") into selected_alignments_RT_ch

  script:
  """
  mkdir 15_selected_alignments_RT
  # python /Users/chrisjackson/PycharmProjects/Yang-and-Smith-RBGV-scripts/07_subset_fasta_from_tree.py
  python /Yang-and-Smith-RBGV-scripts/07_subset_fasta_from_tree.py \
  ${RT_folder} \
  .tre \
  ${alignment_folder} \
  NotApplicable \
  15_selected_alignments_RT
  """
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Run 07_subset_fasta_from_tree.py for MI /////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

process write_alignment_subset_MI_14 {
  // echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
  file(MI_folder) from MI_ch
  file(alignment_folder) from realigned_with_outgroups_ch_3

  output:
  file("16_selected_alignments_MI") into selected_alignments_MI_ch

  script:
  """
  mkdir 16_selected_alignments_MI
  # python /Users/chrisjackson/PycharmProjects/Yang-and-Smith-RBGV-scripts/07_subset_fasta_from_tree.py
  python /Yang-and-Smith-RBGV-scripts/07_subset_fasta_from_tree.py \
  ${MI_folder} \
  .tre \
  ${alignment_folder} \
  NotApplicable \
  16_selected_alignments_MI
  """
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Run 12_strip_names_and_mafft.py for MO //////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

process strip_names_and_realign_MO_15 {
  // echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy', pattern: '16_alignments_stripped_names', saveAs: { "17_alignments_stripped_names_MO"}
  publishDir "${params.outdir}", mode: 'copy', pattern: '17_alignments_stripped_names_realigned', saveAs: { "18_alignments_stripped_names_MO_realigned"}

  input:
  file(selected_alignments_MO) from selected_alignments_MO_ch

  output:
  file("16_alignments_stripped_names") 
  file("17_alignments_stripped_names_realigned")

  script:
  if (!params.no_supercontigs) {
  """
  # python /Users/chrisjackson/PycharmProjects/Yang-and-Smith-RBGV-scripts/12_strip_names_and_mafft.py
  python /Yang-and-Smith-RBGV-scripts/12_strip_names_and_mafft.py \
  ${selected_alignments_MO} \
  -threads_pool ${params.pool} \
  -threads_mafft ${params.threads} 
  """
  } else {
  """
  # python /Users/chrisjackson/PycharmProjects/Yang-and-Smith-RBGV-scripts/12_strip_names_and_mafft.py
  python /Yang-and-Smith-RBGV-scripts/12_strip_names_and_mafft.py \
  ${selected_alignments_MO} \
  -threads_pool ${params.pool} \
  -threads_mafft ${params.threads} \
  -no_supercontigs
  """
  }
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Run 12_strip_names_and_mafft.py for RT //////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

process strip_names_and_realign_RT_16 {
  // echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy', pattern: '16_alignments_stripped_names', saveAs: { "19_alignments_stripped_names_RT"}
  publishDir "${params.outdir}", mode: 'copy', pattern: '17_alignments_stripped_names_realigned', saveAs: { "20_alignments_stripped_names_RT_realigned"}

  input:
  file(selected_alignments_RT) from selected_alignments_RT_ch

  output:
  file("16_alignments_stripped_names") 
  file("17_alignments_stripped_names_realigned")

  script:
  if (!params.no_supercontigs) {
  """
  # python /Users/chrisjackson/PycharmProjects/Yang-and-Smith-RBGV-scripts/12_strip_names_and_mafft.py
  python /Yang-and-Smith-RBGV-scripts/12_strip_names_and_mafft.py \
  ${selected_alignments_RT} \
  -threads_pool ${params.pool} \
  -threads_mafft ${params.threads} 
  """
  } else {
  """
  # python /Users/chrisjackson/PycharmProjects/Yang-and-Smith-RBGV-scripts/12_strip_names_and_mafft.py
  python /Yang-and-Smith-RBGV-scripts/12_strip_names_and_mafft.py \
  ${selected_alignments_RT} \
  -threads_pool ${params.pool} \
  -threads_mafft ${params.threads} -no_supercontigs
  """
  }
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Run 12_strip_names_and_mafft.py for MI //////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

process strip_names_and_realign_MI_17 {
  // echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy', pattern: '16_alignments_stripped_names', saveAs: { "21_alignments_stripped_names_MI"}
  publishDir "${params.outdir}", mode: 'copy', pattern: '17_alignments_stripped_names_realigned', saveAs: { "22_alignments_stripped_names_MI_realigned"}

  input:
  file(selected_alignments_MI) from selected_alignments_MI_ch

  output:
  file("16_alignments_stripped_names")
  file("17_alignments_stripped_names_realigned")

  script:
  if (!params.no_supercontigs) {
  """
  # python /Users/chrisjackson/PycharmProjects/Yang-and-Smith-RBGV-scripts/12_strip_names_and_mafft.py
  python /Yang-and-Smith-RBGV-scripts/12_strip_names_and_mafft.py \
  ${selected_alignments_MI} \
  -threads_pool ${params.pool} \
  -threads_mafft ${params.threads} 
  """
  } else {
  """
  # python /Users/chrisjackson/PycharmProjects/Yang-and-Smith-RBGV-scripts/12_strip_names_and_mafft.py
  python /Yang-and-Smith-RBGV-scripts/12_strip_names_and_mafft.py \
  ${selected_alignments_MI} \
  -threads_pool ${params.pool} \
  -threads_mafft ${params.threads} \
  -no_supercontigs
  """
  }
}

//////////////////////////////////////  End of script ///////////////////////////////////////////////////////
