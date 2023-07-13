#!/usr/bin/env nextflow

///////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Nextflow Pipeline for running ParaGone  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

nextflow.enable.dsl=2

/////////////////////////////////////////////////////////////////////////////////////////
// PRINT PARAGONE-NF NEXTFLOW SCRIPT VERSION
/////////////////////////////////////////////////////////////////////////////////////////

if( params.remove('version') ) {
    println('paragone-nf version 1.0.0, running ParaGone version 0.0.14c')
    exit 0
} 


/////////////////////////////////////////////////////////////////////////////////////////
// CHECK FOR UNRECOGNISED PARAMETERS
/////////////////////////////////////////////////////////////////////////////////////////

allowed_params = [
                  "help",
                  "outdir",
                  "gene_fasta_directory",
                  "gene_name_delimiter",
                  "gene_name_field_num",
                  "external_outgroups_file",
                  "external_outgroups",
                  "internal_outgroups",
                  "pool",
                  "threads",
                  "use_clustal",
                  "mafft_algorithm",
                  "mafft_adjustdirection",
                  "no_trimming",
                  "trimal_terminalonly_off",
                  "trimal_gapthreshold",
                  "trimal_simthreshold",
                  "trimal_cons",
                  "trimal_nogaps",
                  "trimal_noallgaps",
                  "trimal_gappyout",
                  "trimal_strict",
                  "trimal_strictplus",
                  "trimal_automated1",
                  "trimal_block",
                  "trimal_resoverlap",
                  "trimal_seqoverlap",
                  "trimal_w",
                  "trimal_gw",
                  "trimal_sw",
                  "no_cleaning",
                  "run_profiler",
                  "generate_bootstraps",
                  "use_fasttree",
                  "min_tips",
                  "treeshrink_q_value",
                  "cut_deep_paralogs_internal_branch_length_cutoff",
                  "mo",
                  "mi",
                  "rt",
                  "mo_algorithm_paragone",
                  "minimum_taxa",
                  "ignore_1to1_orthologs",
                  "debug",
                  "keep_intermediate_files"
                  ]

params.each { entry ->
  if (! allowed_params.contains(entry.key)) {
      println("The parameter <${entry.key}> is not known");
      exit 0;
  }
}


/////////////////////////////////////////////////////////////////////////////////////////
// CREATE THE HELP MESSAGE AND CHECK FOR MINIMUM PARAMETERS
/////////////////////////////////////////////////////////////////////////////////////////

def helpMessage() {
    log.info """

    Usage:
    The typical command for running the pipeline is as follows:

    nextflow run paragone.nf \
    -c paragone.config \
    --gene_fasta_directory <directory> \
    --external_outgroups_file <file> \
    --external_outgroups <taxon1,taxon2,taxon3...> \
    --internal_outgroups <taxon1,taxon2,taxon3...> \
    --mo \
    --mi \
    --rt \
    -profile <profile> \
    -resume

    Mandatory arguments:

      ############################################################################

      --gene_fasta_directory <directory>       
                                  Path to folder containing paralog fasta files for
                                  each gene (e.g. as output by HybPiper)

      ..and either

      --internal_outgroups <taxon1,taxon2,taxon3...>  
                                  A comma-separated list of taxa present in the 
                                  paralog fasta files to use as outgroups. Default 
                                  is none

      ...and/or

      --external_outgroups_file <file>
                                  File containing fasta sequences of outgroup 
                                  sequences for each gene

      ...and one or more of

      --mo                        Resolve paralogs using the Monophyletic Outgroup 
                                  (MO) algorithm

      --mi                        Resolve paralogs using the Maximum Inclusion (MI)
                                  algorithm

      --rt                        Resolve paralogs using the RooTed ingroups (RT)
                                  algorithm

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
                                  be run at a time during alignment steps
      
      --threads <int>             Number of threads per multiprocessing pool 
                                  instance. Used for programs that support 
                                  multi-threading (e.g. mafft, IQ-TREE). Default 
                                  is 1
      
      --use_clustal               If specified, alignments are performed using 
                                  Clustal Omega rather than MAFFT. If the parameter
                                  "--mafft_adjustdirection" is also provided, 
                                  alignments are performed with MAFFT first, followed 
                                  by realignment using with Clustal Omega

      --mafft_algorithm           If using MAFFT (default), use this alignment 
                                  algorithm (e.g. linsi, einsi, etc.). Default is 
                                  to use the flag --auto

      --mafft_adjustdirection     Allow MAFFT to generate reverse complement sequences, 
                                  as necessary, and align them together with the 
                                  remaining sequences. Note that the first sequence is 
                                  assumed to be in the correct orientation. If this 
                                  parameter is used, it is only applied during the 
                                  first alignment step of the pipeline

      --no_trimming               Do not trim alignments using Trimal

      --trimal_terminalonly_off   Consider all alignment positions when trimming using 
                                  Trimal, rather than only terminal columns

      --trimal_gapthreshold       1 - (fraction of sequences with a gap allowed) when 
                                  trimming alignments with Trimal. Range: [0 - 1]. 
                                  Default is 0.12

      --trimal_simthreshold       Trimal minimum average similarity allowed when 
                                  trimming alignments with Trimal. Range: [0 - 1]

      --trimal_cons               Minimum percentage of positions in the original 
                                  alignment to conserve when trimming alignments with 
                                  Trimal. Range: [0 - 100]

      --trimal_nogaps             Remove all positions with gaps in the alignment when 
                                  trimming

      --trimal_noallgaps          Remove columns composed only by gaps when trimming 
                                  alignments

      --trimal_gappyout           When trimming alignments with Trimal, use the automated 
                                  selection on "gappyout" mode. This method only uses 
                                  information based on gaps distribution

      --trimal_strict             When trimming alignments with Trimal, use the automated 
                                  selection on "strict" mode

      --trimal_strictplus         When trimming alignments with Trimal, use the automated 
                                  selection on "strictplus" mode

      --trimal_automated1         When trimming alignments with Trimal, use a heuristic 
                                  selection of the automatic method based on similarity 
                                  statistics

      --trimal_block              Minimum column block size to be kept in the trimmed 
                                  alignment. Available with manual and automatic (gappyout) 
                                  methods when trimming alignments with Trimal

      --trimal_resoverlap         Minimum overlap of a positions with other positions in 
                                  the column to be considered a "good position" when trimming 
                                  alignments with Trimal. Range: [0 - 1]

      --trimal_seqoverlap         Minimum percentage of "good positions" that a sequence 
                                  must have in order to be conserved when trimming alignments 
                                  with Trimal. Range: [0 - 100]

      --trimal_w                  (half) Window size, score of position i is the average of 
                                  the window (i - n) to (i + n), when trimming alignments with 
                                  Trimal

      --trimal_gw                 (half) Window size only applies to statistics/methods based 
                                  on gaps, when trimming alignments with Trimal

      --trimal_sw                 (half) Window size only applies to statistics/methods based 
                                  on 'similarity, when trimming alignments with Trimal

      --no_cleaning               Do not clean alignments using HmmCleaner.pl

      --run_profiler              If supplied, run the subcommand using cProfile. Saves a 
                                  *.csv file of results

      --generate_bootstraps       Generate bootstraps for trees. For IQ-TREE, uses 
                                  ultrafast bootstrap approximation (UFBoot) 
                                  via the options '-bb 1000 -bnni'. For FastTree, 
                                  uses the SH test. Default is no bootstraps

      --use_fasttree              Use FastTreeMP to generate trees instead of 
                                  IQ-TREE. Default is IQ-TREE

      --min_tips                  The minimum number of tips in a tree after trimming/masking 
                                  tips or pruning deep paralogs; if below this value, no output 
                                  tree is written. Default is 4

      --treeshrink_q_value        q value for TreeShrink; the quantile(s) to set threshold. 
                                  Default is 0.05

      --cut_deep_paralogs_internal_branch_length_cutoff
                                  When pruning long internal branches (putative 
                                  deep paralogs) during the tree QC stage, provide 
                                  a branch length for the maximum allowed internal 
                                  branch length. Default is 0.3

      --mo                        Run the Monophyletic Outgroups (MO) algorithm

      --mi                        Run the Maximum Inclusion (MI) algorithm

      --rt                        Run the RooTed ingroups (RT) algorithm

      --mo_algorithm_paragone     If pruning trees using the MO algorithm, use an updated 
                                  ParaGone implementation rather than the original Yang and
                                  Smith 2014 implementation

      --minimum_taxa              Minimum number of taxa required in pruned trees. Default 
                                  is 4

      --ignore_1to1_orthologs     Do not output 1to1 orthologs, i.e. trees with no paralogs

      --debug                     If supplied, log additional information when running the MO pruning algorithm. '
                                  This can make the log files much larger

      --keep_intermediate_files   N/A


    """.stripIndent()
}

/*
Check that the input directory is provided, and either (at least) a file of external outgroups or list of internal
taxa to use:
*/
if (params.help || !params.gene_fasta_directory || (!params.external_outgroups_file && !params.internal_outgroups) || 
  (!params.mo && !params.mi && !params.rt)) {
  helpMessage()
  println(
  """
  ERROR: Please provide the minimum parameters. The pipeline requires a directory of paralog fasta files 
  via the parameter "--gene_fasta_directory", and either a file of external outgroups via the parameter
  "--external_outgroups_file" or a comma-seperated list on ingroup taxa via the parameter "--internal_outgroups".
  At least one paralogy resoluton method should also be provided via "--mo", "--mi', "--rt".
  """
  )
  exit 0
}


/////////////////////////////////////////////////////////////////////////////////////////
// SET CHANNELS FOR PARALOG FOLDER AND EXTERNAL OUTGROUPS
/////////////////////////////////////////////////////////////////////////////////////////

// External outgroup sequences file:
Channel
  .fromPath("${params.external_outgroups_file}")
  .first()
  .set { outgroups_file_ch }

// Paralog fasta files:
Channel
  .fromPath("${params.gene_fasta_directory}", checkIfExists: true)
  .set { paralogs_ch }


/////////////////////////////////////////////////////////////////////////////////////////
//  CREATE COMMAND LISTS FOR EACH PARAGONE SUBCOMMAND
/////////////////////////////////////////////////////////////////////////////////////////

// Create subcommand lists:
def check_and_align_command_list = []
def alignment_to_tree_command_list = []
def qc_trees_and_extract_fasta_command_list = []
def align_selected_and_tree_command_list = []
def prune_paralogs_command_list = []
def final_alignments_command_list = []

// Parse internal outgroups, if provided:
if (params.internal_outgroups) {
  internal_outgroups_list = params.internal_outgroups?.toString().tokenize(',')
  internal_outgroups_string = ''
  for (outgroup in internal_outgroups_list) {
    internal_outgroup_string = "--internal_outgroup ${outgroup} "
    internal_outgroups_string = internal_outgroups_string + internal_outgroup_string
  }

  check_and_align_command_list << internal_outgroups_string
}

// Parse external outgroup file, if provided:
if (params.external_outgroups_file) {
  File external_outgroups_file = new File(params.external_outgroups_file)
  external_outgroups_file_basename = external_outgroups_file.getName()
  check_and_align_command_list << "--external_outgroups_file ${external_outgroups_file_basename}"
}

// Parse selected external outgroup taxa, if provided:
if (params.external_outgroups) {
  external_outgroups_list = params.external_outgroups?.toString().tokenize(',')
  external_outgroups_string = ''
  for (outgroup in external_outgroups_list) {
    external_outgroup_string = "--external_outgroup ${outgroup} "
    external_outgroups_string = external_outgroups_string + external_outgroup_string
  }

  check_and_align_command_list << external_outgroups_string
}

// Parse general parameters:
if (params.gene_name_delimiter) {
  check_and_align_command_list << "--gene_name_delimiter ${params.gene_name_delimiter}"
}

if (params.gene_name_field_num) {
  check_and_align_command_list << "--gene_name_field_num ${params.gene_name_field_num}"
}

if (params.pool) {
  check_and_align_command_list << "--pool ${params.pool}"
  alignment_to_tree_command_list << "--pool ${params.pool}"
  align_selected_and_tree_command_list << "--pool ${params.pool}"
  final_alignments_command_list << "--pool ${params.pool}"
}

if (params.threads) {
  check_and_align_command_list << "--threads ${params.threads}"
  alignment_to_tree_command_list << "--threads ${params.threads}"
  align_selected_and_tree_command_list << "--threads ${params.threads}"
  final_alignments_command_list << "--threads ${params.threads}"
}

if (params.use_clustal) {
  check_and_align_command_list << "--use_clustal"
  align_selected_and_tree_command_list << "--use_clustal"
  final_alignments_command_list << "--use_clustal"
}

if (params.mafft_algorithm) {
  check_and_align_command_list << "--mafft_algorithm ${params.mafft_algorithm}"
  align_selected_and_tree_command_list << "--mafft_algorithm ${params.mafft_algorithm}"
  final_alignments_command_list << "--mafft_algorithm ${params.mafft_algorithm}"
}

if (params.mafft_adjustdirection) {
  check_and_align_command_list << "--mafft_adjustdirection"
}

if (params.no_trimming) {
  check_and_align_command_list << "--no_trimming"
  align_selected_and_tree_command_list << "--no_trimming"
  final_alignments_command_list << "--no_trimming"
}

if (params.trimal_terminalonly_off) {
  check_and_align_command_list << "--trimal_terminalonly_off"
  align_selected_and_tree_command_list << "--trimal_terminalonly_off"
  final_alignments_command_list << "--trimal_terminalonly_off"
}

if (params.trimal_gapthreshold) {
  check_and_align_command_list << "--trimal_gapthreshold ${params.trimal_gapthreshold}"
  align_selected_and_tree_command_list << "--trimal_gapthreshold ${params.trimal_gapthreshold}"
  final_alignments_command_list << "--trimal_gapthreshold ${params.trimal_gapthreshold}"
}

if (params.trimal_simthreshold) {
  check_and_align_command_list << "--trimal_simthreshold ${params.trimal_simthreshold}"
  align_selected_and_tree_command_list << "--trimal_simthreshold ${params.trimal_simthreshold}"
  final_alignments_command_list << "--trimal_simthreshold ${params.trimal_simthreshold}"
}

if (params.trimal_cons) {
  check_and_align_command_list << "--trimal_cons ${params.trimal_cons}"
  align_selected_and_tree_command_list << "--trimal_cons ${params.trimal_cons}"
  final_alignments_command_list << "--trimal_cons ${params.trimal_cons}"
}

if (params.trimal_nogaps) {
  check_and_align_command_list << "--trimal_nogaps"
  align_selected_and_tree_command_list << "--trimal_nogaps"
  final_alignments_command_list << "--trimal_nogaps"
}

if (params.trimal_noallgaps) {
  check_and_align_command_list << "--trimal_noallgaps"
  align_selected_and_tree_command_list << "--trimal_noallgaps"
  final_alignments_command_list << "--trimal_noallgaps"
}

if (params.trimal_gappyout) {
  check_and_align_command_list << "--trimal_gappyout"
  align_selected_and_tree_command_list << "--trimal_gappyout"
  final_alignments_command_list << "--trimal_gappyout"
}

if (params.trimal_strict) {
  check_and_align_command_list << "--trimal_strict"
  align_selected_and_tree_command_list << "--trimal_strict"
  final_alignments_command_list << "--trimal_strict"
}

if (params.trimal_strictplus) {
  check_and_align_command_list << "--trimal_strictplus"
  align_selected_and_tree_command_list << "--trimal_strictplus"
  final_alignments_command_list << "--trimal_strictplus"
}

if (params.trimal_automated1) {
  check_and_align_command_list << "--trimal_automated1"
  align_selected_and_tree_command_list << "--trimal_automated1"
  final_alignments_command_list << "--trimal_automated1"
}

if (params.trimal_block) {
  check_and_align_command_list << "--trimal_block ${params.trimal_block}"
  align_selected_and_tree_command_list << "--trimal_block ${params.trimal_block}"
  final_alignments_command_list << "--trimal_block ${params.trimal_block}"
}

if (params.trimal_resoverlap) {
  check_and_align_command_list << "--trimal_resoverlap ${params.trimal_resoverlap}"
  align_selected_and_tree_command_list << "--trimal_resoverlap ${params.trimal_resoverlap}"
  final_alignments_command_list << "--trimal_resoverlap ${params.trimal_resoverlap}"
}

if (params.trimal_seqoverlap) {
  check_and_align_command_list << "--trimal_seqoverlap ${params.trimal_seqoverlap}"
  align_selected_and_tree_command_list << "--trimal_seqoverlap ${params.trimal_seqoverlap}"
  final_alignments_command_list << "--trimal_seqoverlap ${params.trimal_seqoverlap}"
}

if (params.trimal_w) {
  check_and_align_command_list << "--trimal_w ${params.trimal_w}"
  align_selected_and_tree_command_list << "--trimal_w ${params.trimal_w}"
  final_alignments_command_list << "--trimal_w ${params.trimal_w}"
}

if (params.trimal_gw) {
  check_and_align_command_list << "--trimal_gw ${params.trimal_gw}"
  align_selected_and_tree_command_list << "--trimal_gw ${params.trimal_gw}"
  final_alignments_command_list << "--trimal_gw ${params.trimal_gw}"
}

if (params.trimal_sw) {
  check_and_align_command_list << "--trimal_sw ${params.trimal_sw}"
  align_selected_and_tree_command_list << "--trimal_sw ${params.trimal_sw}"
  final_alignments_command_list << "--trimal_sw ${params.trimal_sw}"
}

if (params.no_cleaning) {
  check_and_align_command_list << "--no_cleaning"
}

if (params.generate_bootstraps) {
  alignment_to_tree_command_list << "--generate_bootstraps"
  align_selected_and_tree_command_list << "--generate_bootstraps"
}

if (params.use_fasttree) {
  alignment_to_tree_command_list << "--use_fasttree"
  align_selected_and_tree_command_list << "--use_fasttree"
}

if (params.min_tips) {
  qc_trees_and_extract_fasta_command_list << "--min_tips ${params.min_tips}"
}

if (params.treeshrink_q_value) {
  qc_trees_and_extract_fasta_command_list << "--treeshrink_q_value ${params.treeshrink_q_value}"
}

if (params.cut_deep_paralogs_internal_branch_length_cutoff) {
  qc_trees_and_extract_fasta_command_list << "--cut_deep_paralogs_internal_branch_length_cutoff ${params.cut_deep_paralogs_internal_branch_length_cutoff}"
}

if (params.mo) {
  prune_paralogs_command_list << "--mo"
  final_alignments_command_list << "--mo"
}

if (params.mi) {
  prune_paralogs_command_list << "--mi"
  final_alignments_command_list << "--mi"
}

if (params.rt) {
  prune_paralogs_command_list << "--rt"
  final_alignments_command_list << "--rt"
}

if (params.new_mo_algorithm) {
  prune_paralogs_command_list << "--new_mo_algorithm"
}

if (params.minimum_taxa) {
  prune_paralogs_command_list << "--minimum_taxa ${params.minimum_taxa}"
}

if (params.ignore_1to1_orthologs) {
  prune_paralogs_command_list << "--ignore_1to1_orthologs"
}

if (params.debug) {
  prune_paralogs_command_list << "--debug"
}

if (params.run_profiler) {
  check_and_align_command_list << "--run_profiler"
  alignment_to_tree_command_list << "--run_profiler"
  qc_trees_and_extract_fasta_command_list << "--run_profiler"
  align_selected_and_tree_command_list << "--run_profiler"
  prune_paralogs_command_list << "--run_profiler"
  final_alignments_command_list << "--run_profiler"
}


/////////////////////////////////////////////////////////////////////////////////////////
//  DEFINE THE PARAGONE WORKFLOW
/////////////////////////////////////////////////////////////////////////////////////////

workflow {
  CHECK_AND_ALIGN( paralogs_ch, outgroups_file_ch )


  // Get the correct alignment outdir depending on trim/clean settings:
  if (params.no_trimming && params.no_cleaning) {
    alignments_ch = CHECK_AND_ALIGN.out.alignments_ch
  } else if (params.no_cleaning && !params.no_trimming) {
    alignments_ch = CHECK_AND_ALIGN.out.alignments_trimmed_ch
  } else if (params.no_trimming && !params.no_cleaning) {
    alignments_ch = CHECK_AND_ALIGN.out.alignments_hmmcleaned_ch
  } else {
    alignments_ch = CHECK_AND_ALIGN.out.alignments_trimmed_hmmcleaned_ch
  }

  ALIGNMENT_TO_TREE( alignments_ch,
                     CHECK_AND_ALIGN.out.check_and_align_logs_and_reports_ch) 

  QC_TREES_AND_EXTRACT_FASTA( ALIGNMENT_TO_TREE.out.trees_pre_quality_control_ch,
                              alignments_ch,
                              ALIGNMENT_TO_TREE.out.alignment_to_tree_logs_and_reports_ch )

  if (params.external_outgroups_file) {
    ALIGN_SELECTED_AND_TREE( QC_TREES_AND_EXTRACT_FASTA.out.sequences_from_qc_trees_ch,
                             CHECK_AND_ALIGN.out.external_outgroups_sanitised_ch,
                             QC_TREES_AND_EXTRACT_FASTA.out.qc_trees_and_extract_fasta_log_and_reports_ch )
    } else {
      ALIGN_SELECTED_AND_TREE( QC_TREES_AND_EXTRACT_FASTA.out.sequences_from_qc_trees_ch,
                               [], // pass in empty list in place of sanitised external outgroup file
                               QC_TREES_AND_EXTRACT_FASTA.out.qc_trees_and_extract_fasta_log_and_reports_ch ) 
    }

  PRUNE_PARALOGS( ALIGN_SELECTED_AND_TREE.out.pre_paralog_resolution_trees_ch,
                  ALIGN_SELECTED_AND_TREE.out.align_selected_and_tree_logs_and_reports_ch )

  // Get input pruned trees depending on resolution algorithms supplied:
  if (params.mo) {
    mo_ch = PRUNE_PARALOGS.out.pruned_MO_ch
  } else:
    mo_ch = []

  if (params.mi) {
    mi_ch = PRUNE_PARALOGS.out.pruned_MI_ch
  } else:
    mi_ch = []

  if (params.rt) {
    rt_ch = PRUNE_PARALOGS.out.pruned_RT_ch
  } else:
    rt_ch = []

  // Get untrimmed or trimmed alignments depending on options supplied:
  if (!params.no_trimming) {
    pre_prune_alignments_ch = ALIGN_SELECTED_AND_TREE.out.pre_paralog_resolution_alignments_trimmed_ch
  } else:
    pre_prune_alignments_ch = ALIGN_SELECTED_AND_TREE.out.pre_paralog_resolution_alignments_ch

  FINAL_ALIGNMENTS( pre_prune_alignments_ch,
                    mo_ch,
                    mi_ch,
                    rt_ch,
                    PRUNE_PARALOGS.out.prune_paralogs_logs_and_reports_ch )

}


/////////////////////////////////////////////////////////////////////////////////////////
//  DEFINE DSL2 PROCESSES
/////////////////////////////////////////////////////////////////////////////////////////


process CHECK_AND_ALIGN {
  /*
  Run the `paragone check and align` command.
  */

  // echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
    path(paralog_folder)
    path(external_outgroups_file)
    

  output:
    stdout emit: outgroup_coverage_ch
    path ("00_logs_and_reports"), emit: check_and_align_logs_and_reports_ch
    path ("01_input_paralog_fasta_with_sanitised_filenames"), emit: paralogs_sanatised_filenames_ch
    path ("02_alignments"), emit: alignments_ch
    path ("03_alignments_trimmed"), emit: alignments_trimmed_ch
    path ("04_alignments_trimmed_hmmcleaned"), emit: alignments_trimmed_hmmcleaned_ch
    path("*sanitised.*"), emit: external_outgroups_sanitised_ch optional true
    

  script:
    check_and_align_command = "paragone check_and_align ${paralog_folder} " + check_and_align_command_list.join(' ')

    """
    echo "Executing command: ${check_and_align_command}"
    ${check_and_align_command}
    """
}


process ALIGNMENT_TO_TREE {
  /*
  Run the `paragone alignment_to_tree` command.
  */

  // echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
    path(alignments_folder)
    path(logs_and_reports)

  output:
    path ("00_logs_and_reports"), emit: alignment_to_tree_logs_and_reports_ch
    path ("05_trees_pre_quality_control"), emit: trees_pre_quality_control_ch

  script:
    alignment_to_tree_command = "paragone alignment_to_tree ${alignments_folder} " + alignment_to_tree_command_list.join(' ')
    
    """
    echo "Executing command: ${alignment_to_tree_command}"
    ${alignment_to_tree_command}
    """
}


process QC_TREES_AND_EXTRACT_FASTA {
  /*
  Run the `paragone qc_trees_and_extract_fasta` command.
  */

  // echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
    path(trees_pre_quality_control)
    path(alignments_folder)
    path(logs_and_reports)

  output:
    path ("00_logs_and_reports"), emit: qc_trees_and_extract_fasta_log_and_reports_ch
    path("06_trees_trimmed")
    path("07_trees_trimmed_masked")
    path("08_trees_trimmed_masked_cut")
    path("09_sequences_from_qc_trees"), emit: sequences_from_qc_trees_ch

  script:
    qc_trees_and_extract_fasta_command = "paragone qc_trees_and_extract_fasta ${alignments_folder} " + qc_trees_and_extract_fasta_command_list.join(' ')
    
    """
    echo "Executing command: ${qc_trees_and_extract_fasta_command}"
    ${qc_trees_and_extract_fasta_command}
    """
}


process ALIGN_SELECTED_AND_TREE {
  /*
  Run the `paragone align_selected_and_tree` command.
  */

  // echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
    path(alignments_from_qc_trees)
    path(external_outgroups_sanitised)
    path(logs_and_reports)

  output:
    path ("00_logs_and_reports"), emit: align_selected_and_tree_logs_and_reports_ch
    path("10_sequences_from_qc_outgroups_added")
    path("11_pre_paralog_resolution_alignments"), emit: pre_paralog_resolution_alignments_ch
    path("12_pre_paralog_resolution_alignments_trimmed"), emit: pre_paralog_resolution_alignments_trimmed_ch optional true
    path("13_pre_paralog_resolution_trees"), emit: pre_paralog_resolution_trees_ch

  script:
    align_selected_and_tree_command = "paragone align_selected_and_tree ${alignments_from_qc_trees} " + align_selected_and_tree_command_list.join(' ')
    
    """
    echo "Executing command: ${align_selected_and_tree_command}"
    ${align_selected_and_tree_command}
    """
}


process PRUNE_PARALOGS {
  /*
  Run the `paragone prune_paralogs` command.
  */

  // echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
    path(pre_paralog_resolution_trees)
    path(logs_and_reports)

  output:
    path ("00_logs_and_reports"), emit: prune_paralogs_logs_and_reports_ch
    path("14_pruned_MO"), emit: pruned_MO_ch optional true
    path("15_pruned_MI"), emit: pruned_MI_ch optional true
    path("16_pruned_RT"), emit: pruned_RT_ch optional true


  script:
    prune_paralogs_command = "paragone prune_paralogs " + prune_paralogs_command_list.join(' ')
    
    """
    echo "Executing command: ${prune_paralogs_command}"
    ${prune_paralogs_command}
    """
}


process FINAL_ALIGNMENTS {
  /*
  Run the `paragone final_alignments` command.
  */

  // echo true
  label 'in_container'
  publishDir "${params.outdir}", mode: 'copy'

  input:
    path(pre_prune_alignments)
    path(mo_trees)
    path(mi_trees)
    path(rt_trees)
    path(logs_and_reports)


  output:
    path ("00_logs_and_reports"), emit: final_alignments_logs_and_reports_ch
    path("17_selected_sequences_MO") optional true
    path("18_selected_sequences_MI") optional true
    path("19_selected_sequences_RT") optional true
    path("20_MO_stripped_names") optional true
    path("21_MI_stripped_names") optional true
    path("22_RT_stripped_names") optional true
    path("23_MO_final_alignments") optional true
    path("24_MI_final_alignments") optional true
    path("25_RT_final_alignments") optional true
    path("26_MO_final_alignments_trimmed") optional true
    path("27_MI_final_alignments_trimmed") optional true
    path("28_RT_final_alignments_trimmed") optional true

  script:
    final_alignments_command = "paragone final_alignments --keep_intermediate_files " + final_alignments_command_list.join(' ')
    
    """
    echo "Executing command: ${final_alignments_command}"
    ${final_alignments_command}
    """
}


/////////////////////////////////////////////////////////////////////////////////////////
// END OF SCRIPT  
/////////////////////////////////////////////////////////////////////////////////////////

