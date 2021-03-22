#!/usr/bin/env nextflow

///////////////////////////////////////////////////////////////////////////////////////////////////////////
////////// Nextflow Pipeline for processing HybPiper paralogs fasta files with Y&S pipeline  //////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////



def helpMessage() {
    log.info """

    Usage:
    The typical command for running the pipeline is as follows:

    nextflow run yang_and_smith_pipeline_v1_8.nf -c yang_and_smith.config --hybpiper_paralogs_directory <directory> --outgroups_file <file>
    --outgroups <outgroups, csv list>
    -profile <profile>

    Mandatory arguments:

      ##############################################################################################

      --hybpiper_paralogs_directory <directory>    Path to folder containing HybPiper paralog fasta files.
      --outgroups_file <file>                      File containing fasta sequences of outgroup sequences for each gene.
      --outgroups <taxon1,taxon2,taxon3...>        A comma-separated list of outgroup taxa to add, in order of 
                                                   preference.  


      ##############################################################################################

    Optional arguments:

      -profile <profile>                              Configuration profile to use. Can use multiple (comma separated)
                                                      Available: standard (default), slurm
      --output_dir                                    Specify the name of the main output results folder. Default is 'results'.
      --pool <int>                                    Number of threads for the Python multiprocessing pool. Used in e.g. alignments and tree-building steps. Default is 1, so e.g. one alignment will be run at a time during alignment steps.
      --threads <int>                                 Number of threads per multiprocessing pool instance. Used for programs that support multi-threading (e.g. mafft, IQTree). Default is 1.
      --no_supercontigs                               Use this flag if you are processing paralogs from a run of HybPiper that used the --nosupercontigs flag. Mafft alignments are re-aligned using clustal omega, which can do a better job in these cases. Default is off. 
      --process_02_trim_bad_ends_cutoff <int>         Number of bases either side of an internal gap that much match the reference before trimming stops.Default is 5.
      --process_02_trim_bad_ends_size <int>           Number of continuous internal gap positions for an internal gap to be investigated. Default is 15.
      --skip_process_02_trim_bad_ends                 Skips the step trimming the ends of internal gaps.
      --process_04_trim_tips_relative_cutoff <float>  When pruning long tips during the tree QC stage, provide a branch length for the maximum imbalance between sister tips allowed. Default is 0.2.
      --process_04_trim_tips_absolute_cutoff <float>  When pruning long tips during the tree QC stage, provide a branch length for the maximum allowed tip branch length. Default is 0.4.
      --process_06_branch_length_cutoff <float>       When pruning long internal branches (putative deep paralogs) during the tree QC stage, provide a branch length for the maximum allowed internal branch length. Default is 0.3
      --process_06_minimum_taxa <int>                 After the final tree-pruning step prior to paralogy resolution, only retain trees with a minimum number of taxa remaining. Default is 3. 
      --process_09_prune_paralog_MO_minimum_taxa <int>
                                                      For the MO method, only process trees with a minimum number of taxa. Default is 2.
      --process_10_prune_paralogs_RT_minimum_ingroup_taxa <int>
                                                      For the RT method, only process trees with a minumum number of ingroup taxa. Default is 2.  
      --process_11_prune_paralogs_MI_relative_tip_cutoff <float>
                                                      Default is 0.2
      --process_11_prune_paralogs_MI_absolute_tip_cutoff <float>
                                                      Default is 0.4
      --process_11_prune_paralogs_MI_minimum_taxa <int>    
                                                      Default is 2


    """.stripIndent()
}

// Check that input directories are provided
if (params.help || !params.hybpiper_paralogs_directory || !params.outgroups_file || !params.outgroups) {
  helpMessage()
  exit 0
}


// Check that only allowed paramters are provided:
allowed_params = ["outgroups_file", "hybpiper_paralogs_directory", "outgroups", "threads", "pool", "nosupercontigs", "process_02_trim_bad_ends_cutoff","process_02_trim_bad_ends_size", "skip_process_02_trim_bad_ends", "process_04_trim_tips_relative_cutoff", "process_04_trim_tips_absolute_cutoff", "process_06_branch_length_cutoff", "process_06_minimum_taxa", "process_09_prune_paralog_MO_minimum_taxa", "process_10_prune_paralogs_RT_minimum_ingroup_taxa", "process_11_prune_paralogs_MI_relative_tip_cutoff", "process_11_prune_paralogs_MI_absolute_tip_cutoff", "process_11_prune_paralogs_MI_minimum_taxa", "outdir", "help", "no_supercontigs"]

params.each { entry ->
  if (! allowed_params.contains(entry.key)) {
      println("The parameter <${entry.key}> is not known");
      exit 0;
  }
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  Set up channels for input data  ///////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

// Outgroup gene sequences file
Channel
  .fromPath("${params.outgroups_file}", checkIfExists: true)
  .first()
  .set { outgroups_file_ch }


// HybPiper paralog fasta files
Channel
  .fromPath("${params.hybpiper_paralogs_directory}", checkIfExists: true)
  .into { paralogs_ch_1; paralogs_ch_2 }


//////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Run 01_add_outgroup.py  ///////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

process align_paralogs_01 {
  // echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
  file(target_file) from outgroups_file_ch
  file(paralog_folder) from paralogs_ch_1

  output:
  // file("01_outgroup_added")
  file("01_alignments")
  file("03_alignments_hmmcleaned") into (alignments_hmmcleaned_ch_1,
                                         alignments_hmmcleaned_ch_2,
                                         alignments_hmmcleaned_ch_3,
                                         alignments_hmmcleaned_ch_4,
                                         alignments_hmmcleaned_ch_5,
                                         alignments_hmmcleaned_ch_6)

  script:
  outgroups_list = params.outgroups?.tokenize(',')
  outgroups_string = ''

  for (outgroup in outgroups_list) {
    outgroup_string = "-outgroup ${outgroup} "
    outgroups_string = outgroups_string + outgroup_string
  }

  println(outgroups_string)

  if (!params.no_supercontigs) {
  """
  python /opt/01_myScripts/02_YangSmith/01_add_outgroup.py ${target_file} ${paralog_folder} ${outgroups_string} -pool ${params.pool} -threads ${params.threads}
  """
  } else if (params.no_supercontigs) {
  """
  python /opt/01_myScripts/02_YangSmith/01_add_outgroup.py ${outgroups_file} ${paralog_folder} ${outgroups_string} -pool ${params.pool} -threads ${params.threads} -no_supercontigs
  """
  }
}
    

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////  Set up channels to allow process alexander_correct_bad_ends_02 to be skipped /////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////


// Used by process 'alexander_correct_bad_ends_02':
(process_02_noSkip_ch_1, process_02_skip_ch_1) = ( params.skip_process_02_trim_bad_ends
                 ? [Channel.empty(), alignments_hmmcleaned_ch_1]
                 : [alignments_hmmcleaned_ch_1, Channel.empty()] )

// Used by process 'mask_tips_05':
process_02_skip_ch_2 = ( params.skip_process_02_trim_bad_ends
                 ? alignments_hmmcleaned_ch_2 : Channel.empty() )

// Used by process 'write_alignment_subset_07':
process_02_skip_ch_3 = ( params.skip_process_02_trim_bad_ends
                 ? alignments_hmmcleaned_ch_3 : Channel.empty() )

// Used by process 'write_alignment_subset_MO_12':
process_02_skip_ch_4 = ( params.skip_process_02_trim_bad_ends
                 ? alignments_hmmcleaned_ch_4 : Channel.empty() )

// Used by process 'write_alignment_subset_RT_13':
process_02_skip_ch_5 = ( params.skip_process_02_trim_bad_ends
                 ? alignments_hmmcleaned_ch_5 : Channel.empty() )

// Used by process 'write_alignment_subset_MI_14':
process_02_skip_ch_6 = ( params.skip_process_02_trim_bad_ends
                 ? alignments_hmmcleaned_ch_6 : Channel.empty() )


//////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Run 02_alexander_correct_bad_ends.py  /////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

process alexander_correct_bad_ends_02 {
  //echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
  file(alignments_hmmcleaned) from process_02_noSkip_ch_1

  output:
  file("04_alignments_internalcut") into (alignments_internalcut_ch, alignments_internalcut_ch_2, alignments_internalcut_ch_3, alignments_internalcut_ch_4, alignments_internalcut_ch_5, alignments_internalcut_ch_6)

  script:
  targets_list = params.outgroups?.tokenize(',')
  targets_string = ''

  for (target in targets_list) {
    target_string = "-target ${target} "
    targets_string = targets_string + target_string
  }

  """
  python /opt/01_myScripts/02_YangSmith/02_alexander_correct_bad_ends.py -gene_folder ${alignments_hmmcleaned} \
-mingaplength ${params.process_02_trim_bad_ends_size} -fastaext .fasta -matchvalue ${params.process_02_trim_bad_ends_cutoff} -outputfolder 04_alignments_internalcut ${targets_string} -fastaext trimmed_hmm.fasta
  """
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Run 03_alignment_to_tree.py  //////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

process alignment_to_tree_03 {
  //echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
  file(alignments_folder) from alignments_internalcut_ch.mix(process_02_skip_ch_1)

  output:
  file("05_tree_files") into tree_folder_ch

  script:
  """ 
  python /opt/01_myScripts/02_YangSmith/03_alignment_to_tree.py ${alignments_folder} -threads_pool ${params.pool} -threads_iqtree ${params.threads} 
  """
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Run 04_trim_tips.py  //////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

process trim_tips_04 {
  //echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
  file(trees_folder) from tree_folder_ch


  output:
  file("06_trim_tips") into trim_tips_ch

  script:
  """  
  python /opt/01_myScripts/02_YangSmith/04_trim_tips.py ${trees_folder} .treefile ${params.process_04_trim_tips_relative_cutoff} ${params.process_04_trim_tips_absolute_cutoff} 06_trim_tips 
  echo "Finished running trim_tips_04"
  """
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////  Run 05_mask_tips_by_taxonID_transcripts.py  /////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

process mask_tips_05 {
  //echo true
  println("Running mask_tips_05")
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
  file(trimmed_tips_folder) from trim_tips_ch
  file(alignments_folder) from alignments_internalcut_ch_2.mix(process_02_skip_ch_2)

  output:
  file("07_masked_tips") into masked_tips_ch 

  script:
  println('test')
  """
  python /opt/01_myScripts/02_YangSmith/05_mask_tips_by_taxonID_transcripts.py ${trimmed_tips_folder} ${alignments_folder} y 07_masked_tips
  """
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////  Run 06_cut_long_internal_branches.py  ///////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

process cut_long_internal_branches_06 {
  //echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
  file(masked_tips_folder) from masked_tips_ch

  output:
  file("08_cut_internal_branches") into cut_internal_branches_ch

  script:
  """
  mkdir 08_cut_internal_branches
  python /opt/01_myScripts/02_YangSmith/06_cut_long_internal_branches.py ${masked_tips_folder} .mm ${params.process_06_branch_length_cutoff} ${params.process_06_minimum_taxa} 08_cut_internal_branches
  """
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Run 07_subset_fasta_from_tree.py  ///////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

process write_alignment_subset_07 {
  //echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
  file(cut_internal_branches_folder) from cut_internal_branches_ch
  file(alignment_folder) from alignments_internalcut_ch_3.mix(process_02_skip_ch_3)

  output:
  file("09_selected_alignments") into selected_alignments_ch

  script:
  """
  mkdir 09_selected_alignments
  python /opt/01_myScripts/02_YangSmith/07_subset_fasta_from_tree.py ${cut_internal_branches_folder} .subtree ${alignment_folder} NotApplicable 09_selected_alignments 
  """
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Run 08_mafft_alignment_and_iqtree.py  ///////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

process realign_and_iqtree_08 {
   // echo true
   label 'in_container'
   publishDir "${params.outdir}", mode: 'copy'

   input:
   file(selected_alignments_ch) from selected_alignments_ch
   file(target_file) from outgroups_file_ch

   output:
   file("10_realigned") 
   file("11_realigned_trees") into (realigned_trees_ch_1, realigned_trees_ch_2, realigned_trees_ch_3)
   file("in_and_outgroups_list.txt") into (in_out_list_ch, in_out_list_ch_2, in_out_list_ch_3)

   script:
   outgroups_list = params.outgroups?.tokenize(',')
   outgroups_string = ''

   for (outgroup in outgroups_list) {
     outgroup_string = "-outgroup ${outgroup} "
     outgroups_string = outgroups_string + outgroup_string
   }

   println(outgroups_string)

   if (!params.no_supercontigs) {
   """
   python /opt/01_myScripts/02_YangSmith/08_mafft_alignment_and_iqtree.py ${selected_alignments_ch} ${target_file} ${outgroups_string} -threads_pool ${params.pool} -threads_mafft ${params.threads}
   """
   } else {
   """
   python /opt/01_myScripts/02_YangSmith/08_mafft_alignment_and_iqtree.py ${selected_alignments_ch} ${target_file} ${outgroups_string} -threads_pool ${params.pool} -threads_mafft ${params.threads} -no_supercontigs
   """
   }
 }


//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Run 09_prune_paralogs_MO.py  ////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

process prune_paralogs_MO_09 {
  //echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
  file(in_out_file) from in_out_list_ch
  file(realigned_trees_folder) from realigned_trees_ch_1

  output:
  file("12_prune_MO_trees") into MO_ch

  script:
  """
  mkdir 12_prune_MO_trees
  python /opt/01_myScripts/02_YangSmith/09_prune_paralogs_MO.py ${realigned_trees_folder} .treefile ${params.process_09_prune_paralog_MO_minimum_taxa} 12_prune_MO_trees ${in_out_file} 
  """
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Run 10_prune_paralogs_RT.py  ////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

process prune_paralogs_RT_10 {
  //echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
  file(in_out_file) from in_out_list_ch_2
  file(realigned_trees_folder) from realigned_trees_ch_2

  output:
  file("13_prune_RT_trees") into RT_ch

  script:
  """
  mkdir 13_prune_RT_trees
  python /opt/01_myScripts/02_YangSmith/10_prune_paralogs_RT.py ${realigned_trees_folder} .treefile 13_prune_RT_trees ${params.process_10_prune_paralogs_RT_minimum_ingroup_taxa} ${in_out_file} 
  """
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Run 11_prune_paralogs_MI.py  ////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

process prune_paralogs_MI_11 {
  //echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
  file(in_out_file) from in_out_list_ch_3
  file(realigned_trees_folder) from realigned_trees_ch_3

  output:
  file("14_prune_MI_trees") into MI_ch

  script:
  """
  mkdir 14_prune_MI_trees
  python /opt/01_myScripts/02_YangSmith/11_prune_paralogs_MI.py ${realigned_trees_folder} .treefile ${params.process_11_prune_paralogs_MI_relative_tip_cutoff} ${params.process_11_prune_paralogs_MI_absolute_tip_cutoff} ${params.process_11_prune_paralogs_MI_minimum_taxa} 14_prune_MI_trees 
  """
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Run 07_subset_fasta_from_tree.py for MO /////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

process write_alignment_subset_MO_12 {
  //echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
  file(MO_folder) from MO_ch
  file(alignment_folder) from alignments_internalcut_ch_4.mix(process_02_skip_ch_4)

  output:
  file("15_selected_alignments_MO") into selected_alignments_MO_ch

  script:
  """
  mkdir 15_selected_alignments_MO
  python /opt/01_myScripts/02_YangSmith/07_subset_fasta_from_tree.py ${MO_folder} .tre ${alignment_folder} NotApplicable 15_selected_alignments_MO 
  """
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Run 07_subset_fasta_from_tree.py for RT /////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

process write_alignment_subset_RT_13 {
  //echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
  file(RT_folder) from RT_ch
  file(alignment_folder) from alignments_internalcut_ch_5.mix(process_02_skip_ch_5)

  output:
  file("16_selected_alignments_RT") into selected_alignments_RT_ch

  script:
  """
  mkdir 16_selected_alignments_RT
  python /opt/01_myScripts/02_YangSmith/07_subset_fasta_from_tree.py ${RT_folder} .tre ${alignment_folder} NotApplicable 16_selected_alignments_RT 
  """
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Run 07_subset_fasta_from_tree.py for MI /////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

process write_alignment_subset_MI_14 {
  //echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
  file(MI_folder) from MI_ch
  file(alignment_folder) from alignments_internalcut_ch_6.mix(process_02_skip_ch_6)

  output:
  file("17_selected_alignments_MI") into selected_alignments_MI_ch

  script:
  """
  mkdir 17_selected_alignments_MI
  python /opt/01_myScripts/02_YangSmith/07_subset_fasta_from_tree.py ${MI_folder} .tre ${alignment_folder} NotApplicable 17_selected_alignments_MI
  """
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Run 12_strip_names_and_mafft.py for MO //////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

process strip_names_and_realign_MO_15 {
  //echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy', pattern: '16_alignments_stripped_names', saveAs: { "18_alignments_stripped_names_MO"}
  publishDir "${params.outdir}", mode: 'copy', pattern: '17_alignments_stripped_names_realigned', saveAs: { "19_alignments_stripped_names_MO_realigned"}

  input:
  file(selected_alignments_MO) from selected_alignments_MO_ch

  output:
  file("16_alignments_stripped_names") 
  file("17_alignments_stripped_names_realigned")

  script:
  if (!params.no_supercontigs) {
  """
  python /opt/01_myScripts/02_YangSmith/12_strip_names_and_mafft.py ${selected_alignments_MO} -threads_pool ${params.pool} -threads_mafft ${params.threads} 
  """
  } else {
  """
  python /opt/01_myScripts/02_YangSmith/12_strip_names_and_mafft.py ${selected_alignments_MO} -threads_pool ${params.pool} -threads_mafft ${params.threads} -no_supercontigs
  """
  }
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Run 12_strip_names_and_mafft.py for RT //////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

process strip_names_and_realign_RT_16 {
  //echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy', pattern: '16_alignments_stripped_names', saveAs: { "20_alignments_stripped_names_RT"}
  publishDir "${params.outdir}", mode: 'copy', pattern: '17_alignments_stripped_names_realigned', saveAs: { "21_alignments_stripped_names_RT_realigned"}

  input:
  file(selected_alignments_RT) from selected_alignments_RT_ch

  output:
  file("16_alignments_stripped_names") 
  file("17_alignments_stripped_names_realigned")

  script:
  if (!params.no_supercontigs) {
  """
  python /opt/01_myScripts/02_YangSmith/12_strip_names_and_mafft.py ${selected_alignments_RT} -threads_pool ${params.pool} -threads_mafft ${params.threads} 
  """
  } else {
  """
  python /opt/01_myScripts/02_YangSmith/12_strip_names_and_mafft.py ${selected_alignments_RT} -threads_pool ${params.pool} -threads_mafft ${params.threads} -no_supercontigs
  """
  }
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Run 12_strip_names_and_mafft.py for MI //////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

process strip_names_and_realign_MI_17 {
  //echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy', pattern: '16_alignments_stripped_names', saveAs: { "22_alignments_stripped_names_MI"}
  publishDir "${params.outdir}", mode: 'copy', pattern: '17_alignments_stripped_names_realigned', saveAs: { "23_alignments_stripped_names_MI_realigned"}

  input:
  file(selected_alignments_MI) from selected_alignments_MI_ch

  output:
  file("16_alignments_stripped_names")
  file("17_alignments_stripped_names_realigned")

  script:
  if (!params.no_supercontigs) {
  """
  python /opt/01_myScripts/02_YangSmith/12_strip_names_and_mafft.py ${selected_alignments_MI} -threads_pool ${params.pool} -threads_mafft ${params.threads} 
  """
  } else {
  """
  python /opt/01_myScripts/02_YangSmith/12_strip_names_and_mafft.py ${selected_alignments_MI} -threads_pool ${params.pool} -threads_mafft ${params.threads} -no_supercontigs
  """
  }
}


//////////////////////////////////////  End of script ///////////////////////////////////////////////////////
