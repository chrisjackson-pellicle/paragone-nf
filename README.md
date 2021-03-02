# paralogy_resolution_tutorial
GAP tutorial for containerised Yang and Smith paralogy resolution pipeline 


https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4209138/

https://www.biorxiv.org/content/10.1101/2020.08.21.261925v2

https://bitbucket.org/dfmoralesb/target_enrichment_orthology/src/master/


This tuturial assumes that you have Singularity and Nextflow installed, and that you have the Y_and_S Singularity image downloaded. In addition, you should have the Nextflow pipeline script `xxx.nf` and its corresponding config file `xxx.config`.


## Input data

### Paralog sequences

If you have used the Nextflow pipeline `hybpiper_pipeline_v1_6_NO_INTRONERATE.nf` to run HybPiper, it will have produced the following subfolders in your main `results` folder:

 - `09_paralogs`
 - `10_paralogs_noChimeras`

See the tutorial[LINK] for a full description of the files in these output folders. Briefly, folder `09_paralogs` contains a fasta file for each gene in your HybPiper target file. Each fasta file contains the 'main' contig selected by HybPiper for each sample. Where HybPiper has detected putative paralog contigs, these sequences are also included; in such cases, the main contig has the fasta header suffix `.main`, whereas putative paralogs have the suffix `.0`, `.1` etc. Folder `10_paralogs_noChimeras` contains the same data, except putative chimeric contigs (see LINK) have been removed.        

***Tutorial step 1:*** Copy folder `09_paralogs` into your current working directory (i.e. the directory containing `xxx.nf` and `xxx.config`.

### Outgroup sequences

 - fasta file contianing outgroup sequences for each gene in your target file e.g.

HEAD OF FILE HERE.

Describe fasta header formatting - has to be same as HybPiper target file e.g. a suffix comprising a dash followed by a unique gene identifier.

***Tutorial step 2:*** Copy 


## Step 02: Running the pipeline

nextflow_20_04 run alex_YS_pipeline_v1_6.nf -c nextflow_alex_YS.config -profile slurm -resume --hybpiper_paralogs_directory 06_paralogs --target_file Angiosperms353_targetSequences.fasta --outgroups Ambtr --pool 4 --threads 4


Optional arguments:

      -profile <profile>                              Configuration profile to use. Can use multiple (comma separated)
                                                      Available: standard (default), slurm
      --pool <int>                                    Number of threads for Python multiprocessing pool. 
                                                      Default is 1.
      --threads <int>                                 Number of threads per pool instance. Default is 1.
      --skip_process_01_hmmcleaner                    Skip the HmmCleaner step of process_01 (useful when HybPiper has been run with -nosupercontigs). 
      --no_supercontigs
      --process_02_trim_bad_ends_cutoff <int>         Default is 5.
      --process_02_trim_bad_ends_size <int>           Default is 15
      --skip_process_02_trim_bad_ends
      --process_04_trim_tips_relative_cutoff          Default is 0.2
      --process_04_trim_tips_absolute_cutoff          Default is 0.4
      --process_06_branch_length_cutoff               Default is 0.3
      --process_06_minimum_taxa                       Default is 3 
      --process_09_prune_paralog_MO_minimum_taxa      Default is 2
      --process_10_prune_paralogs_RT_minimum_ingroup_taxa
                                                      Default is 2  
      --process_11_prune_paralogs_MI_relative_tip_cutoff
                                                      Default is 0.2
      --process_11_prune_paralogs_MI_absolute_tip_cutoff
                                                      Default is 0.4
      --process_11_prune_paralogs_MI_minimim_taxa    
                                                      Default is 2


## Step 03 Interpreting output data

ls -1 results
01_outgroup_added
02_alignments
03_alignments_hmmcleaned
04_alignments_internalcut
05_tree_files
06_trim_tips
07_masked_tips
08_cut_internal_branches
09_selected_alignments
10_realigned
11_realigned_trees
12_prune_MO_trees
13_prune_RT_trees
14_prune_MI_trees
15_selected_alignments_MO
16_selected_alignments_RT
17_selected_alignments_MI
18_alignments_stripped_names_MO
19_alignments_stripped_names_MO_realigned
20_alignments_stripped_names_RT
21_alignments_stripped_names_RT_realigned
22_alignments_stripped_names_MI
23_alignments_stripped_names_MI_realigned
in_and_outgroups_list.txt

### general notes

###################################


## Description

Bioinformatic sequence recovery for universal target-capture bait kits can be [substantially improved][12] by appropriate tailoring of target files to the group under study. To enable the best possible locus recovery from [Angiosperms353][10] capture data, we have developed an expanded target file (`mega353.fasta`) incorporating sequences from over 550 transcriptomes from the [1KP][9] project. To maximise computational efficiency we provide the script `filter_mega353.py`, which can be used to subsample the `mega353.fasta` file based on user-selected taxa or taxon groups. These groups can be defined using unique 1KP transcriptome codes, species, families, orders, or broader groups (e.g. Basal Eudicots, Monocots, etc). In addition, we  provide the script `BYO_transcriptome.py`, which can be used to add sequences from any transcriptome to any protein-coding nucleotide target file. These tailored and customised target files can be used directly in target-capture pipelines such as [HybPiper][8]. 

**Data files**
- `mega353.fasta` A target file for use with target enrichment datasets captured using the Angiosperms353 bait kit. 
- `filtering_options.csv` A comma-separated values file listing the options available for filtering the `mega353.fasta` file. This reference file can also be produced by the `filter_mega353.py` script (see below).

**Scripts**
- `filter_mega353.py` A script to filter the `mega353.fasta` target file.
- `BYO_transcriptome.py` A script to add sequences from any transcriptome dataset to any target file containing protein-coding sequences.

**Manuscript** 
- https://www.biorxiv.org/content/10.1101/2020.10.04.325571v1

## Dependencies

Dependencies for `filter_mega353.py`
- Python 3.7 or higher
- [BioPython][4] 1.76 or higher
- [pandas][11] 1.0.3 or higher


Please see the Wiki page [Installing dependencies][5] for further details.

## Installation

Assuming all dependencies are installed, either:

1. Download the NewTargets package directly from the repository home page and unzip it. Note that the `mega353.fasta` file in the unzipped package is also provided as a `.zip` file, and will need to be unzipped separately. 
2. Clone the repository using the command `git clone https://github.com/chrisjackson-pellicle/NewTargets.git`. Unzip the `mega353.zip` file.


[testing anchor link](#newtargets)

## Scripts

### filter_mega353.py

*Input*:
- `mega353.fasta`. The expanded Angiosperms353 target file.
- `select_file`. A text file containing a list of IDs; sequences from these samples will be retained in the filtered output target file.   



Please see the Wiki page [BYO_transcriptome][6] for further details.

[1]: https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate/ "Link to EXONERATE download page"
[2]: http://hmmer.org/ "Link to HMMER download page"
[3]: https://mafft.cbrc.jp/alignment/software/ "Link to MAFFT download page"
[4]: https://biopython.org/wiki/Download "Link to BioPython download page"
[5]: https://github.com/chrisjackson-pellicle/NewTargets/wiki/Installing-script-dependencies "Link to Installing dependencies Wiki page"
[6]: https://github.com/chrisjackson-pellicle/NewTargets/wiki/BYO_transcriptome.py:-adding-transcriptome-sequences-to-a-target-file "Link to BYO_transcriptome Wiki page"
[7]: https://github.com/chrisjackson-pellicle/NewTargets/wiki/filter_mega353.py:-filtering-the-mega353.fasta-target-file "Link to filter_mega353 Wiki page"
[8]: https://github.com/mossmatters/HybPiper/ "Link to the HybPiper GitHub repository"
[9]: https://sites.google.com/a/ualberta.ca/onekp/ "Link to the 1000 Plants website"
[10]: https://dx.doi.org/10.1093%2Fsysbio%2Fsyy086 "Link to Angiosperms353 manuscript"
[11]: https://pandas.pydata.org/pandas-docs/stable/getting_started/install.html "Link to pandas installation instructions"
[12]: https://www.biorxiv.org/content/10.1101/2020.10.04.325571v1 "Link to NewTargets manuscript"



