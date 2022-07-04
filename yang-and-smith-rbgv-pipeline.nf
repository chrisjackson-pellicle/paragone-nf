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

      --batch_size <int>          Number of paralog fasta files to align/generate 
                                  trees for in one batch. Default is 20

      -profile <profile>          Configuration profile to use. Can use multiple 
                                  (comma separated). Available: standard (default), 
                                  slurm

      --outdir                    Specify the name of the main output results 
                                  folder. Default is 'results'

      --pool <int>                Number of threads for the Python multiprocessing 
                                  pool. Used in e.g. alignments and tree-building 
                                  steps. Default is 1, so e.g. one alignment will 
                                  be run at a time during alignment steps
      
      --threads <int>             Number of threads per multiprocessing pool 
                                  instance. Used for programs that support 
                                  multi-threading (e.g. mafft, IQ-TREE). Default 
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
                                  tip branch length. Default is 0.4

      --mafft_algorithm           If using MAFFT (default), use this alignment 
                                  algorithm (e.g. linsi, einsi, etc.). Default is 
                                  to use the flag --auto 

      --use_muscle                Use MUSCLE to align sequences intead of MAFFT. 
                                  Default is MAFFT

      --use_fasttree              Use FastTreeMP to generate trees instead of 
                                  IQ-TREE. Default is IQ-TREE

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

      --bootstraps                Generate bootstraps for trees. For IQ-TREE, uses 
                                  ultrafast bootstrap approximation (UFBoot) 
                                  via the options '-bb 1000 -bnni'. For FastTree, 
                                  uses the SH test. Default is no bootstraps


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
"process_11_prune_paralogs_MI_minimum_taxa", "outdir", "help", "no_supercontigs", "bootstraps",
"use_fasttree", "use_muscle", "mafft_algorithm", "batch_size"]

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


process CHECK_AND_BATCH_PARALOGS_01 {
  /*
  Run script 01a_check_outgroups_and_batch.py
  */

  // echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
    path(external_outgroups_file)
    path(paralog_folder)

  output:
    stdout emit: outgroup_coverage_ch
    path("01_batch_folders/batch_*"), emit: batch_folders_ch
    path("outgroup_coverage_report.tsv")
    path("*sanitised.*"), emit: external_outgroups_sanitised_ch optional true
    

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
  

    """
    python /Yang-and-Smith-RBGV-scripts/01a_check_outgroups_and_batch.py \
    ${paralog_folder} \
    ${external_outgroups_file_string} \
    ${external_outgroups_string} \
    ${internal_outgroups_string} \
    -batch_size ${params.batch_size}
    """
}


process ALIGN_AND_HMMCLEAN_02 {
  /*
  02a_align_and_hmmclean.py
  */

  // echo true
  label 'in_container'
  publishDir "${params.outdir}/02_alignments", mode: 'copy', pattern: "batch_*_alignments"
  publishDir "${params.outdir}/03_alignments_hmmcleaned", mode: 'copy', pattern: "batch_*_alignments_hmmcleaned"

  input:
    path(alignments_folder)

  output:
    path("batch_*_alignments")
    path("batch_*_alignments_hmmcleaned"), emit: batch_alignments_hmmcleaned_ch

  script:
    if (params.no_supercontigs) {
      no_supercontigs_string = "-no_supercontigs"
    } else {
      no_supercontigs_string = ''
    }

    if (params.mafft_algorithm) {
      mafft_algorithm_string = "-mafft_algorithm ${params.mafft_algorithm}"
    } else {
      mafft_algorithm_string = ''
    }

    if (params.use_muscle) {
      muscle_string = "-use_muscle"
    } else {
      muscle_string = ''
    }

    if (params.no_supercontigs) {
    """ 
    python /Yang-and-Smith-RBGV-scripts/02a_align_and_hmmclean.py \
    ${alignments_folder} \
    ${no_supercontigs_string} \
    ${muscle_string} \
    ${mafft_algorithm_string} \
    -pool ${params.pool} \
    -threads ${params.threads}

    batch=(batch_*_alignments_clustal_hmmcleaned)
    batch_rename=\${batch/clustal_/}
    mv \$batch \$batch_rename
    """

    } else {
    """ 
    python /Yang-and-Smith-RBGV-scripts/02a_align_and_hmmclean.py \
    ${alignments_folder} \
    ${no_supercontigs_string} \
    ${muscle_string} \
    ${mafft_algorithm_string} \
    -pool ${params.pool} \
    -threads ${params.threads}
    """
    }
}


process ALIGNMENT_TO_TREE_03 {
  /*
  run script 03a_alignment_to_tree.py
  */

  // echo true
  label 'in_container'
  publishDir "${params.outdir}/04_tree_files", mode: 'copy', pattern: "batch_*_alignments_hmmcleaned_tree_files"

  input:
    path(alignments_folder)

  output:
    path("batch_*_alignments_hmmcleaned_tree_files"), emit: batch_treefiles_ch

  script:
    if (params.bootstraps) {
      bootstraps_string = "-generate_bootstraps"
    } else {
      bootstraps_string = ''
    }

    if (params.use_fasttree) {
      fasttree_string = "-use_fasttree"
    } else {
      fasttree_string = ''
    }

    """ 
    python /Yang-and-Smith-RBGV-scripts/03a_alignment_to_tree.py \
    ${alignments_folder} \
    -threads_pool ${params.pool} \
    -threads_iqtree ${params.threads} \
    ${bootstraps_string} \
    ${fasttree_string}
    """
}


process TRIM_TIPS_04 {
  /*
  Run script 04_trim_tips.py
  */

  // echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
    path(trees_folder)

  output:
    path("05_trim_tips")

  script:
    """ 
    mkdir all_trees_combined
    for trees_folder in ${trees_folder}
      do
        echo \${trees_folder}
        cp -r \${trees_folder}/* all_trees_combined
      done

    python /Yang-and-Smith-RBGV-scripts/04_trim_tips.py \
    all_trees_combined \
    .treefile \
    ${params.process_04_trim_tips_relative_cutoff} \
    ${params.process_04_trim_tips_absolute_cutoff} \
    05_trim_tips
    echo "Finished running trim_tips_04"

    """
}


process MASK_TIPS_05 {
  /*
  Run script 05_mask_tips_by_taxonID_transcripts.py
  */


  // echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy', pattern: "06_masked_tips"

  input:
    path(trimmed_tips_folder)
    path(alignments_folder)

  output:
    path("06_masked_tips"), emit: masked_tips_ch
    path("all_hmmclean_alignments_combined"), emit: all_hmmclean_alignments_combined_ch

  script:
    """
    mkdir all_hmmclean_alignments_combined
    for alignment_folder in ${alignments_folder}
      do
        echo \${alignment_folder}
        cp -r \${alignment_folder}/* all_hmmclean_alignments_combined
      done

    python /Yang-and-Smith-RBGV-scripts/05_mask_tips_by_taxonID_transcripts.py \
    ${trimmed_tips_folder} \
    all_hmmclean_alignments_combined \
    y \
    06_masked_tips
    """
}


process CUT_LONG_INTERNAL_BRANCHES_06 {
  /*
  Run script 06_cut_long_internal_branches.py
  */

  // echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
    path(masked_tips_folder)

  output:
    path("07_cut_internal_branches")

  script:
    """
    mkdir 07_cut_internal_branches
    python /Yang-and-Smith-RBGV-scripts/06_cut_long_internal_branches.py \
    ${masked_tips_folder} \
    .mm \
    ${params.process_06_branch_length_cutoff} \
    ${params.process_06_minimum_taxa} \
    07_cut_internal_branches
    """
}


process WRITE_ALIGNMENT_SUBSET_07 {
  /*
  Run script 07a_subset_fasta_from_tree_and_batch.py
  */

  // echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'
                    
  input:
    path(cut_internal_branches_folder)
    path(alignment_folder)

  output:
    path("08_selected_alignments_batch_folders/batch_*"), emit: subset_batch_folders_ch

  script:
    """
    mkdir 08_selected_alignments
    python /Yang-and-Smith-RBGV-scripts/07a_subset_fasta_from_tree_and_batch.py \
    ${cut_internal_branches_folder} \
    .subtree \
    ${alignment_folder} \
    08_selected_alignments \
    -from_cut_internal_branches \
    -batch_size ${params.batch_size}
    """
}


process REALIGN_AND_IQTREE_08 {
  /*
  Run script 08_mafft_alignment_and_iqtree.py
  */

  // echo  true
  label 'in_container'
  publishDir "${params.outdir}/09_outgroups_added", mode: 'copy', pattern: "batch_*_outgroups_added"
  publishDir "${params.outdir}/10_realigned", mode: 'copy', pattern: "batch_*_outgroups_added_alignments"
  publishDir "${params.outdir}/11_realigned_trees", mode: 'copy', pattern: "batch_*_outgroups_added_alignments_tree_files"

  input:
    path(selected_alignments_ch)
    path(hmm_cleaned_alignments)
    path(external_outgroups_file)

  output:
    path("batch_*_outgroups_added")
    path("batch_*_outgroups_added_alignments"), emit: realigned_fasta_ch
    path("batch_*_outgroups_added_alignments_tree_files"), emit: realigned_trees_ch
    path("in_and_outgroups_list*.txt"), emit: in_and_outgroups_list_ch


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

    if (params.no_supercontigs) {
      no_supercontigs_string = "-no_supercontigs"
    } else {
      no_supercontigs_string = ''
    }

    if (params.use_fasttree) {
      fasttree_string = "-use_fasttree"
    } else {
      fasttree_string = ''
    }

    if (params.use_muscle) {
      muscle_string = "-use_muscle"
    } else {
      muscle_string = ''
    }

    if (params.mafft_algorithm) {
      mafft_algorithm_string = "-mafft_algorithm ${params.mafft_algorithm}"
    } else {
      mafft_algorithm_string = ''
    }

    if (params.no_supercontigs) {
    """
    python /Yang-and-Smith-RBGV-scripts/08a_alignment_and_tree.py \
    ${hmm_cleaned_alignments} \
    ${selected_alignments_ch} \
    ${external_outgroups_file_string} \
    ${external_outgroups_string} \
    ${internal_outgroups_string} \
    ${no_supercontigs_string} \
    ${muscle_string} \
    ${mafft_algorithm_string} \
    ${fasttree_string} \
    -threads_pool ${params.pool} \
    -threads_mafft ${params.threads}

    batch_trees=(batch_*_outgroups_added_alignments_clustal_tree_files)
    batch_trees_rename=\${batch_trees/clustal_/}
    mv \$batch_trees \$batch_trees_rename


    batch_aln=(batch_*_outgroups_added_alignments)
    batch_aln_clustal=(batch_*_outgroups_added_alignments_clustal)
    rm -r \$batch_aln
    mv \$batch_aln_clustal \$batch_aln



    """
    } else {
    """
    python /Yang-and-Smith-RBGV-scripts/08a_alignment_and_tree.py \
    ${hmm_cleaned_alignments} \
    ${selected_alignments_ch} \
    ${external_outgroups_file_string} \
    ${external_outgroups_string} \
    ${internal_outgroups_string} \
    ${no_supercontigs_string} \
    ${muscle_string} \
    ${mafft_algorithm_string} \
    ${fasttree_string} \
    -threads_pool ${params.pool} \
    -threads_mafft ${params.threads}
    """
    }
  }


process RESOLVE_POLYTOMIES {
  /*
  Run script resolve_polytomies.py
  */

  label 'in_container'
  publishDir "${params.outdir}/11_realigned_trees", mode: 'copy', pattern: "treefiles_polytomies_resolved"
  // echo true

    input:
    path(realigned_trees_folder)

  output:
    path("treefiles_polytomies_resolved"), emit: trees_resolved_polytomies_ch

  script:
    """
    mkdir all_realigned_trees_combined
    for tree_folder in ${realigned_trees_folder}
      do
        echo \${tree_folder}
        cp -r \${tree_folder}/* all_realigned_trees_combined
      done

    python /Yang-and-Smith-RBGV-scripts/resolve_polytomies.py all_realigned_trees_combined
    """
}


process PRUNE_PARALOGS_MO_09 {
  /*
  Run script 09_prune_paralogs_MO.py
  */

  // echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
    path(in_out_file)
    path(realigned_resolved_trees_folder)

  output:
    path("12_prune_MO_trees")

  script:
    """
    mkdir 12_prune_MO_trees

    python /Yang-and-Smith-RBGV-scripts/09_prune_paralogs_MO.py \
    treefiles_polytomies_resolved \
    .treefile \
    ${params.process_09_prune_paralog_MO_minimum_taxa} \
    12_prune_MO_trees \
    ${in_out_file}
    """
}


process PRUNE_PARALOGS_RT_10 {
  /*
  Run script 10_prune_paralogs_RT.py
  */

  // echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
  path(in_out_file)
  path(realigned_resolved_trees_folder)

  output:
  path("13_prune_RT_trees")

  script:
  """
  mkdir 13_prune_RT_trees

  python /Yang-and-Smith-RBGV-scripts/10_prune_paralogs_RT.py \
  treefiles_polytomies_resolved \
  .treefile 13_prune_RT_trees \
  ${params.process_10_prune_paralogs_RT_minimum_ingroup_taxa} \
  ${in_out_file}
  """
}


process PRUNE_PARALOGS_MI_11 {
  /*
  Run script 11_prune_paralogs_MI.py
  */

  // echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
  path(in_out_file)
  path(realigned_resolved_trees_folder)

  output:
  path("14_prune_MI_trees")

  script:
  """
  mkdir 14_prune_MI_trees

  python /Yang-and-Smith-RBGV-scripts/11_prune_paralogs_MI.py \
  treefiles_polytomies_resolved \
  .treefile \
  ${params.process_11_prune_paralogs_MI_relative_tip_cutoff} \
  ${params.process_11_prune_paralogs_MI_absolute_tip_cutoff} \
  ${params.process_11_prune_paralogs_MI_minimum_taxa} \
  14_prune_MI_trees

  """
}


process WRITE_ALIGNMENT_SUBSET_MO_12 {
  /*
  Run script 07a_subset_fasta_from_tree_and_batch.py
  */

  // echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
  path(MO_folder)
  path(alignment_folder)

  output:
  path("15_selected_alignments_MO_batch_folders/batch_*"), emit: mo_subset_batch_folders_ch

  script:
  """
  mkdir all_realignments_combined
    for alignment in ${alignment_folder}
      do
        echo \${alignment}
        cp -r \${alignment}/* all_realignments_combined
      done


  mkdir 15_selected_alignments_MO

  python /Yang-and-Smith-RBGV-scripts/07a_subset_fasta_from_tree_and_batch.py \
  ${MO_folder} \
  .tre \
  all_realignments_combined \
  15_selected_alignments_MO \
  -batch_size ${params.batch_size}

  mv 08_selected_alignments_batch_folders 15_selected_alignments_MO_batch_folders

  """
}


process WRITE_ALIGNMENT_SUBSET_RT_13 {
  /*
  Run script 07a_subset_fasta_from_tree_and_batch.py
  */

  // echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
  path(RT_folder)
  path(alignment_folder)

  output:
  path("16_selected_alignments_RT_batch_folders/batch_*"), emit: rt_subset_batch_folders_ch

  script:
  """
  mkdir all_realignments_combined
    for alignment in ${alignment_folder}
      do
        echo \${alignment}
        cp -r \${alignment}/* all_realignments_combined
      done

  mkdir 16_selected_alignments_RT

  python /Yang-and-Smith-RBGV-scripts/07a_subset_fasta_from_tree_and_batch.py \
  ${RT_folder} \
  .tre \
  all_realignments_combined \
  16_selected_alignments_RT \
  -batch_size ${params.batch_size}

  mv 08_selected_alignments_batch_folders 16_selected_alignments_RT_batch_folders
  """
}


process WRITE_ALIGNMENT_SUBSET_MI_14 {
  /*
  Run script 07a_subset_fasta_from_tree_and_batch.py
  */

  // echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
  path(MI_folder)
  path(alignment_folder)

  output:
  path("17_selected_alignments_MI_batch_folders/batch_*"), emit: mi_subset_batch_folders_ch

  script:
  """
  mkdir all_realignments_combined
    for alignment in ${alignment_folder}
      do
        echo \${alignment}
        cp -r \${alignment}/* all_realignments_combined
      done

  mkdir 17_selected_alignments_MI

  python /Yang-and-Smith-RBGV-scripts/07a_subset_fasta_from_tree_and_batch.py \
  ${MI_folder} \
  .tre \
  all_realignments_combined \
  17_selected_alignments_MI \
  -batch_size ${params.batch_size}

  mv 08_selected_alignments_batch_folders 17_selected_alignments_MI_batch_folders
  """
}


process STRIP_NAMES_AND_REALIGN_MO_15 {
  /*
  Run script 12a_strip_names_and_align.py for MO
  */

  // echo true
  label 'in_container'
  publishDir "${params.outdir}/18_alignments_stripped_names_MO", mode: 'copy', pattern: "batch_*_stripped_names/*", saveAs: { filename -> file(filename).getName() }
  publishDir "${params.outdir}/19_alignments_stripped_names_MO_realigned", mode: 'copy', pattern: "batch_*_alignments/*", saveAs: { filename -> file(filename).getName() }


  input:
  path(selected_alignments_MO)

  output:
  path("batch_*_stripped_names/*"), emit: stripped_names_ch
  path("batch_*_alignments/*"), emit: stripped_names_aligned_ch

  script:
  if (params.no_supercontigs) {
      no_supercontigs_string = "-no_supercontigs"
    } else {
      no_supercontigs_string = ''
    }

  if (params.no_supercontigs) {
  """
  echo ${selected_alignments_MO}
  python /Yang-and-Smith-RBGV-scripts/12a_strip_names_and_align.py \
  ${selected_alignments_MO} \
  -threads_pool ${params.pool} \
  -threads_mafft ${params.threads} \
  ${no_supercontigs_string} 


  batch_aln=(batch_*_stripped_names_alignments)
  batch_aln_clustal=(batch_*_stripped_names_alignments_clustal)
  rm -r \$batch_aln
  mv \$batch_aln_clustal \$batch_aln
  """

  } else {
  """
  echo ${selected_alignments_MO}
  python /Yang-and-Smith-RBGV-scripts/12a_strip_names_and_align.py \
  ${selected_alignments_MO} \
  -threads_pool ${params.pool} \
  -threads_mafft ${params.threads} \
  ${no_supercontigs_string} 
  """

  }
}


process STRIP_NAMES_AND_REALIGN_RT_16 {
  /*
  Run script 12a_strip_names_and_align.py for RT
  */


  // echo true
  label 'in_container'
  publishDir "${params.outdir}/20_alignments_stripped_names_RT", mode: 'copy', pattern: "batch_*_stripped_names/*", saveAs: { filename -> file(filename).getName() }
  publishDir "${params.outdir}/21_alignments_stripped_names_RT_realigned", mode: 'copy', pattern: "batch_*_alignments/*", saveAs: { filename -> file(filename).getName() }

  input:
  path(selected_alignments_RT)

  output:
  path("batch_*_stripped_names/*"), emit: stripped_names_ch
  path("batch_*_alignments/*"), emit: stripped_names_aligned_ch

  script:
  if (params.no_supercontigs) {
      no_supercontigs_string = "-no_supercontigs"
    } else {
      no_supercontigs_string = ''
    }

  if (params.no_supercontigs) {
  """
  python /Yang-and-Smith-RBGV-scripts/12a_strip_names_and_align.py \
  ${selected_alignments_RT} \
  -threads_pool ${params.pool} \
  -threads_mafft ${params.threads} \
  ${no_supercontigs_string}

  batch_aln=(batch_*_stripped_names_alignments)
  batch_aln_clustal=(batch_*_stripped_names_alignments_clustal)
  rm -r \$batch_aln
  mv \$batch_aln_clustal \$batch_aln
  """

  } else {
  """
  python /Yang-and-Smith-RBGV-scripts/12a_strip_names_and_align.py \
  ${selected_alignments_RT} \
  -threads_pool ${params.pool} \
  -threads_mafft ${params.threads} \
  ${no_supercontigs_string}
  """
  }
}


process STRIP_NAMES_AND_REALIGN_MI_17 {
  /*
  Run script 12a_strip_names_and_align.py for MI
  */

  // echo true
  label 'in_container'
  publishDir "${params.outdir}/22_alignments_stripped_names_MI", mode: 'copy', pattern: "batch_*_stripped_names/*", saveAs: { filename -> file(filename).getName() }
  publishDir "${params.outdir}/23_alignments_stripped_names_MI_realigned", mode: 'copy', pattern: "batch_*_alignments/*"  , saveAs: { filename -> file(filename).getName() }

  input:
  path(selected_alignments_MI)

  output:
  path("batch_*_stripped_names/*"), emit: stripped_names_ch
  path("batch_*_alignments/*"), emit: stripped_names_aligned_ch

  script:
  if (params.no_supercontigs) {
      no_supercontigs_string = "-no_supercontigs"
    } else {
      no_supercontigs_string = ''
    }

  if (params.no_supercontigs) {
  """
  python /Yang-and-Smith-RBGV-scripts/12a_strip_names_and_align.py \
  ${selected_alignments_MI} \
  -threads_pool ${params.pool} \
  -threads_mafft ${params.threads} \
  ${no_supercontigs_string}

  batch_aln=(batch_*_stripped_names_alignments)
  batch_aln_clustal=(batch_*_stripped_names_alignments_clustal)
  rm -r \$batch_aln
  mv \$batch_aln_clustal \$batch_aln
  """

  } else {
  """
  python /Yang-and-Smith-RBGV-scripts/12a_strip_names_and_align.py \
  ${selected_alignments_MI} \
  -threads_pool ${params.pool} \
  -threads_mafft ${params.threads} \
  ${no_supercontigs_string}
  """
  }
}




////////////////////////
//  Define workflows  //
////////////////////////


workflow {
  CHECK_AND_BATCH_PARALOGS_01( outgroups_file_ch, paralogs_ch )
  CHECK_AND_BATCH_PARALOGS_01.out.outgroup_coverage_ch.subscribe { log.info '\n' + "$it" + '\n' } // Log outgroup coverage statistics to screen and `nextflow.log` file
  // CHECK_AND_BATCH_PARALOGS_01.out.batch_folders_ch.flatten()subscribe { log.info '\n' + "$it" + '\n' }

  ALIGN_AND_HMMCLEAN_02( CHECK_AND_BATCH_PARALOGS_01.out.batch_folders_ch.flatten() )
  ALIGNMENT_TO_TREE_03( ALIGN_AND_HMMCLEAN_02.out.batch_alignments_hmmcleaned_ch ) 
  TRIM_TIPS_04( ALIGNMENT_TO_TREE_03.out.batch_treefiles_ch.collect() )
  MASK_TIPS_05( TRIM_TIPS_04.out, 
                ALIGN_AND_HMMCLEAN_02.out.batch_alignments_hmmcleaned_ch.collect() )
  CUT_LONG_INTERNAL_BRANCHES_06( MASK_TIPS_05.out.masked_tips_ch )
  WRITE_ALIGNMENT_SUBSET_07( CUT_LONG_INTERNAL_BRANCHES_06.out, 
                             MASK_TIPS_05.out.all_hmmclean_alignments_combined_ch )

  if (params.external_outgroups_file) {
      REALIGN_AND_IQTREE_08( WRITE_ALIGNMENT_SUBSET_07.out.subset_batch_folders_ch.flatten(), 
                             MASK_TIPS_05.out.all_hmmclean_alignments_combined_ch, 
                             CHECK_AND_BATCH_PARALOGS_01.out.external_outgroups_sanitised_ch.first() ) // WARN: The operator `first` is useless when applied to a value channel which returns a single value by definition
  } else {
      REALIGN_AND_IQTREE_08( WRITE_ALIGNMENT_SUBSET_07.out.subset_batch_folders_ch.flatten(), 
                             MASK_TIPS_05.out.all_hmmclean_alignments_combined_ch.first(), 
                             [] )
  }

  RESOLVE_POLYTOMIES( REALIGN_AND_IQTREE_08.out.realigned_trees_ch.collect() )

  PRUNE_PARALOGS_MO_09( REALIGN_AND_IQTREE_08.out.in_and_outgroups_list_ch.first(), 
                        RESOLVE_POLYTOMIES.out.trees_resolved_polytomies_ch )
  PRUNE_PARALOGS_RT_10( REALIGN_AND_IQTREE_08.out.in_and_outgroups_list_ch.first(), 
                        RESOLVE_POLYTOMIES.out.trees_resolved_polytomies_ch )
  PRUNE_PARALOGS_MI_11( REALIGN_AND_IQTREE_08.out.in_and_outgroups_list_ch.first(), 
                        RESOLVE_POLYTOMIES.out.trees_resolved_polytomies_ch )
  
  WRITE_ALIGNMENT_SUBSET_MO_12( PRUNE_PARALOGS_MO_09.out, 
                                REALIGN_AND_IQTREE_08.out.realigned_fasta_ch.collect() )
  WRITE_ALIGNMENT_SUBSET_RT_13( PRUNE_PARALOGS_RT_10.out, 
                                REALIGN_AND_IQTREE_08.out.realigned_fasta_ch.collect() )
  WRITE_ALIGNMENT_SUBSET_MI_14( PRUNE_PARALOGS_MI_11.out, 
                                REALIGN_AND_IQTREE_08.out.realigned_fasta_ch.collect() )
  
  STRIP_NAMES_AND_REALIGN_MO_15( WRITE_ALIGNMENT_SUBSET_MO_12.out.mo_subset_batch_folders_ch.flatten() )
  STRIP_NAMES_AND_REALIGN_RT_16( WRITE_ALIGNMENT_SUBSET_RT_13.out.rt_subset_batch_folders_ch.flatten() )
  STRIP_NAMES_AND_REALIGN_MI_17( WRITE_ALIGNMENT_SUBSET_MI_14.out.mi_subset_batch_folders_ch.flatten() )

}

/////////////////////
//  End of script  //
/////////////////////

