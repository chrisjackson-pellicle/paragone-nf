#!/usr/bin/env nextflow

///////////////////////////////////////////////////////////////////////////////////////////////////////////
////////// Nextflow Pipeline for processing HybPiper paralogs fasta files with Y&S pipeline  //////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

nextflow.enable.dsl=2

def helpMessage() {
    log.info """

    Usage:
    The typical command for running the pipeline is as follows:

    nextflow run yang-and-smith-rbgv-pipeline.nf \
    -c yang-and-smith-rbgv.config \
    --hybpiper_paralogs_directory <directory> \
    --external_outgroups_file <file> \
    --external_outgroups <taxon1,taxon2,taxon3...> \
    --internal_outgroups <taxon1,taxon2,taxon3...> \
    -profile <profile> \
    -resume
 
    Mandatory arguments:

      ############################################################################

      --hybpiper_paralogs_directory <directory>       
                                  Path to folder containing HybPiper paralog fasta 
                                  files.

      ..and either

      --internal_outgroups <taxon1,taxon2,taxon3...>  
                                  A comma-separated list of taxa present in the 
                                  paralog fasta files to use as outgroups. Default 
                                  is none

      ...and/or

      --external_outgroups_file <file>
                                  File containing fasta sequences of outgroup 
                                  sequences for each gene

      ############################################################################

    Optional arguments:
    
      --internal_outgroups <taxon1,taxon2,taxon3...>  
                                  A comma-separated list of taxa present in the 
                                  paralog fasta files to use as outgroups. Default 
                                  is none

      --external_outgroups_file <file>                
                                  File containing fasta sequences of outgroup 
                                  sequences for each gene

      --external_outgroups <taxon1,taxon2,taxon3...>  
                                  A comma-separated list of outgroup taxa to add 
                                  from the outgroups_file. Default is all

      -profile <profile>          Configuration profile to use. Can use multiple 
                                  (comma separated). Available: standard (default), 
                                  slurm

      --outdir                    Specify the name of the main output results 
                                  folder. Default is 'results'

      --pool <int>                Number of threads for the Python multiprocessing 
                                  pool. Used in e.g. alignments and tree-building 
                                  steps. Default is 1, so e.g. one alignment will 
                                  be run at a time during alignment steps.
      
      --threads <int>             Number of threads per multiprocessing pool 
                                  instance. Used for programs that support 
                                  multi-threading (e.g. mafft, IQTree). Default 
                                  is 1
      
      --no_supercontigs           Use this flag if you are processing paralogs 
                                  from a run of HybPiper that used the 
                                  --nosupercontigs flag. Mafft alignments are 
                                  re-aligned using clustal omega, which can do a 
                                  better job in these cases. Default is off
      
      --process_04_trim_tips_relative_cutoff <float>  
                                  When pruning long tips during the tree QC stage, 
                                  provide a branch length for the maximum imbalance 
                                  between sister tips allowed. Default is 0.2
      
      --process_04_trim_tips_absolute_cutoff <float>  
                                  When pruning long tips during the tree QC stage, 
                                  provide a branch length for the maximum allowed 
                                  tip branch length. Default is 0.4.
     
      --process_06_branch_length_cutoff <float>       
                                  When pruning long internal branches (putative 
                                  deep paralogs)during the tree QC stage, provide 
                                  a branch length for the maximum allowed internal 
                                  branch length. Default is 0.3
      
      --process_06_minimum_taxa <int>                 
                                  After the final tree-pruning step prior to 
                                  paralogy resolution, only retain trees with a 
                                  minimum number of taxa remaining. Default is 3
      
      --process_09_prune_paralog_MO_minimum_taxa <int>
                                  For the MO method, only process trees with a 
                                  minimum number of taxa. Default is 2
      
      --process_10_prune_paralogs_RT_minimum_ingroup_taxa <int>
                                  For the RT method, only process trees with a 
                                  minumum number of ingroup taxa. Default is 2
      
      --process_11_prune_paralogs_MI_relative_tip_cutoff <float>
                                  Default is 0.2
      
      --process_11_prune_paralogs_MI_absolute_tip_cutoff <float>
                                  Default is 0.4
      
      --process_11_prune_paralogs_MI_minimum_taxa <int>    
                                  Default is 2

      --iqtree_ufbootstraps       Generate bootstraps for trees using IQ-TREE's 
                                  ultrafasta bootstrap approximation (UFBoot) 
                                  via the options "-bb 1000 -bnni". Default is
                                  no bootstraps


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
"process_11_prune_paralogs_MI_minimum_taxa", "outdir", "help", "no_supercontigs", "iqtree_ufbootstraps"]

params.each { entry ->
  if (! allowed_params.contains(entry.key)) {
      println("The parameter <${entry.key}> is not known");
      exit 0;
  }
}



//////////////////////////////////////
//  Set up channels for input data  //
//////////////////////////////////////


// External outgroup sequences file:
Channel
  .fromPath("${params.external_outgroups_file}")
  .first()
  .set { outgroups_file_ch }

// HybPiper paralog fasta files:
Channel
  .fromPath("${params.hybpiper_paralogs_directory}", checkIfExists: true)
  .set { paralogs_ch }



/////////////////////////////
//  DEFINE DSL2 PROCESSES  //
/////////////////////////////


process ALIGN_PARALOGS_01 {
  /*
  Run script 01_check_outgroups_align_and_hmmclean.py
  */

  // echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
    path(external_outgroups_file)
    path(paralog_folder)

  output:
    stdout emit: outgroup_coverage_ch
    path("00_paralogs_gene_names_sanitised"), emit: paralogs_gene_names_sanitised_ch
    path("01_alignments")
    path("02_alignments_hmmcleaned"), emit: alignments_hmmcleaned_ch
    path("outgroup_coverage_report.tsv")
    path("*_sanitised.fasta"), emit: external_outgroups_sanitised_ch optional true
    

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
    python /Yang-and-Smith-RBGV-scripts/01_check_outgroups_align_and_hmmclean.py \
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


process ALIGNMENT_TO_TREE_02 {
  /*
  Run script 03_alignment_to_tree.py
  */

  //echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
    path(alignments_folder)

  output:
    path("03_tree_files")

  script:
    if (params.iqtree_ufbootstraps) {
      bootstraps_string = "-generate_bootstraps"
    } else {
      bootstraps_string = ''
    """ 
    python /Yang-and-Smith-RBGV-scripts/03_alignment_to_tree.py \
    ${alignments_folder} \
    -threads_pool ${params.pool} \
    -threads_iqtree ${params.threads} \
    ${bootstraps_string}
    """
}


process TRIM_TIPS_03 {
  /*
  Run script 04_trim_tips.py
  */

  // echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
    path(trees_folder)

  output:
    path("04_trim_tips")

  script:
    """  
    python /Yang-and-Smith-RBGV-scripts/04_trim_tips.py \
    ${trees_folder} \
    .treefile \
    ${params.process_04_trim_tips_relative_cutoff} \
    ${params.process_04_trim_tips_absolute_cutoff} \
    04_trim_tips
    echo "Finished running trim_tips_04"
    """
}


process MASK_TIPS_04 {
  /*
  Run script 05_mask_tips_by_taxonID_transcripts.py
  */


  // echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
    path(trimmed_tips_folder)
    path(alignments_folder)

  output:
    path("05_masked_tips")

  script:
    """
    python /Yang-and-Smith-RBGV-scripts/05_mask_tips_by_taxonID_transcripts.py \
    ${trimmed_tips_folder} \
    ${alignments_folder} \
    y \
    05_masked_tips
    """
}


process CUT_LONG_INTERNAL_BRANCHES_05 {
  /*
  Run script 06_cut_long_internal_branches.py
  */

  // echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
    path(masked_tips_folder)

  output:
    path("06_cut_internal_branches")

  script:
    """
    mkdir 06_cut_internal_branches
    python /Yang-and-Smith-RBGV-scripts/06_cut_long_internal_branches.py \
    ${masked_tips_folder} \
    .mm \
    ${params.process_06_branch_length_cutoff} \
    ${params.process_06_minimum_taxa} \
    06_cut_internal_branches
    """
}


process WRITE_ALIGNMENT_SUBSET_06 {
  /*
  Run script 07_subset_fasta_from_tree.py
  */

  // echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'
                    
  input:
    path(cut_internal_branches_folder)
    path(alignment_folder)

  output:
    path("07_selected_alignments")

  script:
    """
    mkdir 07_selected_alignments
    python /Yang-and-Smith-RBGV-scripts/07_subset_fasta_from_tree.py \
    ${cut_internal_branches_folder} \
    .subtree \
    ${alignment_folder} \
    07_selected_alignments \
    -from_cut_internal_branches
    """
}


process REALIGN_AND_IQTREE_07 {
  /*
  Run script 08_mafft_alignment_and_iqtree.py
  */

  // echo  true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
    path(selected_alignments_ch)
    path(hmm_cleaned_alignments)
    path(external_outgroups_file)

  output:
    path("08_outgroups_added")
    path("09_realigned"), emit: realigned_fasta_ch
    path("10_realigned_trees"), emit: realigned_trees_ch
    path("in_and_outgroups_list.txt"), emit: in_and_outgroups_list_ch

  script:
    if (params.iqtree_ufbootstraps) {
      bootstraps_string = "-generate_bootstraps"
    } else {
      bootstraps_string = ''

    if (params.no_supercontigs) {
      no_supercontigs_string = "-no_supercontigs"
    } else {
      no_supercontigs_string = ''

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

    """
    python /Yang-and-Smith-RBGV-scripts/08_mafft_alignment_and_iqtree.py \
    ${hmm_cleaned_alignments} \
    ${selected_alignments_ch} \
    ${external_outgroups_file_string} \
    ${external_outgroups_string} \
    ${internal_outgroups_string} \
    ${bootstraps_string} \
    ${no_supercontigs_string}
    -threads_pool ${params.pool} \
    -threads_mafft ${params.threads}

    """
}


process PRUNE_PARALOGS_MO_08 {
  /*
  Run script 09_prune_paralogs_MO.py
  */

  // echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
    path(in_out_file)
    path(realigned_trees_folder)

  output:
    path("11_prune_MO_trees")

  script:
    """
    mkdir 11_prune_MO_trees
    python /Yang-and-Smith-RBGV-scripts/09_prune_paralogs_MO.py \
    ${realigned_trees_folder} \
    .treefile \
    ${params.process_09_prune_paralog_MO_minimum_taxa} \
    11_prune_MO_trees \
    ${in_out_file}
    """
}


process PRUNE_PARALOGS_RT_09 {
  /*
  Run script 10_prune_paralogs_RT.py
  */

  // echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
  path(in_out_file)
  path(realigned_trees_folder)

  output:
  path("12_prune_RT_trees")

  script:
  """
  mkdir 12_prune_RT_trees
  python /Yang-and-Smith-RBGV-scripts/10_prune_paralogs_RT.py \
  ${realigned_trees_folder} \
  .treefile 12_prune_RT_trees \
  ${params.process_10_prune_paralogs_RT_minimum_ingroup_taxa} \
  ${in_out_file}
  """
}


process PRUNE_PARALOGS_MI_10 {
  /*
  Run script 11_prune_paralogs_MI.py
  */

  // echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
  path(in_out_file)
  path(realigned_trees_folder)

  output:
  path("13_prune_MI_trees")

  script:
  """
  mkdir 13_prune_MI_trees
  python /Yang-and-Smith-RBGV-scripts/11_prune_paralogs_MI.py \
  ${realigned_trees_folder} \
  .treefile \
  ${params.process_11_prune_paralogs_MI_relative_tip_cutoff} \
  ${params.process_11_prune_paralogs_MI_absolute_tip_cutoff} \
  ${params.process_11_prune_paralogs_MI_minimum_taxa} \
  13_prune_MI_trees
  """
}


process WRITE_ALIGNMENT_SUBSET_MO_11 {
  /*
  Run script 07_subset_fasta_from_tree.py
  */

  // echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
  path(MO_folder)
  path(alignment_folder)

  output:
  path("14_selected_alignments_MO")

  script:
  """
  mkdir 14_selected_alignments_MO
  python /Yang-and-Smith-RBGV-scripts/07_subset_fasta_from_tree.py \
  ${MO_folder} \
  .tre \
  ${alignment_folder} \
  14_selected_alignments_MO
  """
}


process WRITE_ALIGNMENT_SUBSET_RT_12 {
  /*
  Run script 07_subset_fasta_from_tree.py
  */

  // echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
  path(RT_folder)
  path(alignment_folder)

  output:
  path("15_selected_alignments_RT")

  script:
  """
  mkdir 15_selected_alignments_RT
  python /Yang-and-Smith-RBGV-scripts/07_subset_fasta_from_tree.py \
  ${RT_folder} \
  .tre \
  ${alignment_folder} \
  15_selected_alignments_RT
  """
}


process WRITE_ALIGNMENT_SUBSET_MI_13 {
  /*
  Run script 07_subset_fasta_from_tree.py
  */

  // echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
  path(MI_folder)
  path(alignment_folder)

  output:
  path("16_selected_alignments_MI")

  script:
  """
  mkdir 16_selected_alignments_MI
  python /Yang-and-Smith-RBGV-scripts/07_subset_fasta_from_tree.py \
  ${MI_folder} \
  .tre \
  ${alignment_folder} \
  16_selected_alignments_MI
  """
}


process STRIP_NAMES_AND_REALIGN_MO_14 {
  /*
  Run script 12_strip_names_and_mafft.py for MO
  */

  // echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy', pattern: '16_alignments_stripped_names', saveAs: { "17_alignments_stripped_names_MO"}
  publishDir "${params.outdir}", mode: 'copy', pattern: '17_alignments_stripped_names_realigned', saveAs: { "18_alignments_stripped_names_MO_realigned"}

  input:
  path(selected_alignments_MO)

  output:
  path("16_alignments_stripped_names") 
  path("17_alignments_stripped_names_realigned")

  script:
  if (!params.no_supercontigs) {
  """
  python /Yang-and-Smith-RBGV-scripts/12_strip_names_and_mafft.py \
  ${selected_alignments_MO} \
  -threads_pool ${params.pool} \
  -threads_mafft ${params.threads} 
  """
  } else {
  """
  python /Yang-and-Smith-RBGV-scripts/12_strip_names_and_mafft.py \
  ${selected_alignments_MO} \
  -threads_pool ${params.pool} \
  -threads_mafft ${params.threads} \
  -no_supercontigs
  """
  }
}


process STRIP_NAMES_AND_REALIGN_RT_15 {
  /*
  Run script 12_strip_names_and_mafft.py for RT
  */


  // echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy', pattern: '16_alignments_stripped_names', saveAs: { "19_alignments_stripped_names_RT"}
  publishDir "${params.outdir}", mode: 'copy', pattern: '17_alignments_stripped_names_realigned', saveAs: { "20_alignments_stripped_names_RT_realigned"}

  input:
  path(selected_alignments_RT)

  output:
  path("16_alignments_stripped_names") 
  path("17_alignments_stripped_names_realigned")

  script:
  if (!params.no_supercontigs) {
  """
  python /Yang-and-Smith-RBGV-scripts/12_strip_names_and_mafft.py \
  ${selected_alignments_RT} \
  -threads_pool ${params.pool} \
  -threads_mafft ${params.threads} 
  """
  } else {
  """
  python /Yang-and-Smith-RBGV-scripts/12_strip_names_and_mafft.py \
  ${selected_alignments_RT} \
  -threads_pool ${params.pool} \
  -threads_mafft ${params.threads} -no_supercontigs
  """
  }
}


process STRIP_NAMES_AND_REALIGN_MI_16 {
  /*
  Run script 12_strip_names_and_mafft.py for MI
  */

  // echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy', pattern: '16_alignments_stripped_names', saveAs: { "21_alignments_stripped_names_MI"}
  publishDir "${params.outdir}", mode: 'copy', pattern: '17_alignments_stripped_names_realigned', saveAs: { "22_alignments_stripped_names_MI_realigned"}

  input:
  path(selected_alignments_MI)

  output:
  path("16_alignments_stripped_names")
  path("17_alignments_stripped_names_realigned")

  script:
  if (!params.no_supercontigs) {
  """
  python /Yang-and-Smith-RBGV-scripts/12_strip_names_and_mafft.py \
  ${selected_alignments_MI} \
  -threads_pool ${params.pool} \
  -threads_mafft ${params.threads} 
  """
  } else {
  """
  python /Yang-and-Smith-RBGV-scripts/12_strip_names_and_mafft.py \
  ${selected_alignments_MI} \
  -threads_pool ${params.pool} \
  -threads_mafft ${params.threads} \
  -no_supercontigs
  """
  }
}




////////////////////////
//  Define workflows  //
////////////////////////


workflow {
  ALIGN_PARALOGS_01( outgroups_file_ch, paralogs_ch )
  ALIGN_PARALOGS_01.out.outgroup_coverage_ch.subscribe { log.info '\n' + "$it" + '\n' } // Log outgroup coverage statistics to screen and `nextflow.log` file
  ALIGNMENT_TO_TREE_02( ALIGN_PARALOGS_01.out.alignments_hmmcleaned_ch )
  TRIM_TIPS_03( ALIGNMENT_TO_TREE_02.out )
  MASK_TIPS_04( TRIM_TIPS_03.out, ALIGN_PARALOGS_01.out.alignments_hmmcleaned_ch )
  CUT_LONG_INTERNAL_BRANCHES_05( MASK_TIPS_04.out )
  WRITE_ALIGNMENT_SUBSET_06( CUT_LONG_INTERNAL_BRANCHES_05.out, ALIGN_PARALOGS_01.out.alignments_hmmcleaned_ch )

  if (params.external_outgroups_file) {
      REALIGN_AND_IQTREE_07( WRITE_ALIGNMENT_SUBSET_06.out, ALIGN_PARALOGS_01.out.alignments_hmmcleaned_ch, ALIGN_PARALOGS_01.out.external_outgroups_sanitised_ch )
  } else {
      REALIGN_AND_IQTREE_07( WRITE_ALIGNMENT_SUBSET_06.out, ALIGN_PARALOGS_01.out.alignments_hmmcleaned_ch, [] )
  }

  PRUNE_PARALOGS_MO_08( REALIGN_AND_IQTREE_07.out.in_and_outgroups_list_ch, REALIGN_AND_IQTREE_07.out.realigned_trees_ch )
  PRUNE_PARALOGS_RT_09( REALIGN_AND_IQTREE_07.out.in_and_outgroups_list_ch, REALIGN_AND_IQTREE_07.out.realigned_trees_ch )
  PRUNE_PARALOGS_MI_10( REALIGN_AND_IQTREE_07.out.in_and_outgroups_list_ch, REALIGN_AND_IQTREE_07.out.realigned_trees_ch )
  
  WRITE_ALIGNMENT_SUBSET_MO_11( PRUNE_PARALOGS_MO_08.out, REALIGN_AND_IQTREE_07.out.realigned_fasta_ch )
  WRITE_ALIGNMENT_SUBSET_RT_12( PRUNE_PARALOGS_RT_09.out, REALIGN_AND_IQTREE_07.out.realigned_fasta_ch )
  WRITE_ALIGNMENT_SUBSET_MI_13( PRUNE_PARALOGS_MI_10.out, REALIGN_AND_IQTREE_07.out.realigned_fasta_ch )
  
  STRIP_NAMES_AND_REALIGN_MO_14( WRITE_ALIGNMENT_SUBSET_MO_11.out )
  STRIP_NAMES_AND_REALIGN_RT_15( WRITE_ALIGNMENT_SUBSET_RT_12.out )
  STRIP_NAMES_AND_REALIGN_MI_16( WRITE_ALIGNMENT_SUBSET_MI_13.out )

}



/////////////////////
//  End of script  //
/////////////////////

