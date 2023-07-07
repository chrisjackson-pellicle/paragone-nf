# paragone-nf: a paralogy resolution pipeline

Current version: 1.0.0 (July 2023)

## Original orthology inference manuscript, documentation and scripts

This pipeline makes use of the **paralogy resolution** (also described as **orthology inference**) approaches described and implemented by Yang and Smith 2014 [here][4]. These approaches have since been adapted for target capture datasets as described in the manuscript [here][11]. The original documentation and scripts can be found [here][12].

## paragone-nf: containerised and pipelined using Singularity and Nextflow

The software [ParaGone](https://github.com/chrisjackson-pellicle/ParaGone) includes modified versions of the Yang and Smith scripts (with added functionality, logging and reporting), along with new scripts to facilitate seamless running of the pipeline. To further simplify the process, here I’ve provided a [Singularity][13] container based on the Linux distribution Ubuntu 22.04, containing Paragone along with all the dependencies: [IQTree][3], [Clustal Omega][6], [MAFFT][5], [BioPython][15], [HMMCleaner][2], [trimal][1], [FastTreeMP][23], [TreeShrink](https://github.com/uym2/TreeShrink). The container is called `hybpiper-paragone.sif`.

To run the paralogy resolution pipeline using this container, I’ve provided a [Nextflow][14] script that uses the software in the Singularity container. This pipeline runs all steps with a single command. The pipeline script is called `paragone.nf`. It comes with an associated config file called `paragone.config`. The only input required is a folder containing `.fasta` files for each of your target-capture loci, including paralogs, and an optional `.fasta` file containing outgroup sequences (used by some of the paralogy resolution methods, see below). The number of parallel processes running at any time, as well as computing resources given to each process (e.g. number of CPUs, amount of RAM etc) can be configured by the user by modifying the provided config file. The pipeline can be run directly on your local computer, or on an HPC system submitting jobs via a scheduler (e.g. SLURM, PBS, etc).

## Input data

### Paralog sequences

If you have used the Nextflow/Singularity pipeline `hybpiper-nf` to run HybPiper, as described [here][16], your main `results` folder will contain the following subfolders:

- `11_paralogs`
- `12_paralogs_noChimeras`

See the hybpiper-nf wiki entry [Output folders and files][7] for a full description of the files in these output folders. Briefly, folder `11_paralogs` contains a fasta file for each gene in your HybPiper target file. Each fasta file contains the 'main' contig selected by HybPiper for each sample. Where HybPiper has detected putative paralog contigs, these sequences are also included; in such cases, the main contig has the fasta header suffix `.main`, whereas putative paralogs have the suffix `.0`, `.1` etc. Folder `12_paralogs_noChimeras` contains the same data, except putative chimeric contigs (see [here][8] for an explanation) have been removed.

### Outgroup sequences

Some of the paralogy resolution methods used in this pipeline require an outgroup sequence for each of your genes. These outgroup sequences can be provided in two ways.

1) Designating one or more taxa in your HybPiper paralog files as outgroups, via the `--internal_outgroups <taxon1,taxon2,taxon3...>` option. For example, if your paralog `fasta` files contain sequences from the taxa `79686` and `79689`, you could designate these sequences as outgroups via `--internal_outgroups 79686,79689`.

2) Providing a fasta file (e.g. `outgroups.fasta`) containing 'external' outgroup sequences via the option `--external_outgroups_file outgroups.fasta`. The sequences in the file should have the same fasta header formatting and gene names as your HybPiper target file. For example, if you have used the Angiosperms353 target file for you HybPiper analysis, and you wish to use sequences from *Sesame* as your outgroup, your `outgroups.fasta` file might contain the following:

       >sesame-6995
       gtgggatatgaacaaaatccattgagcttgtattactgtta...
       >sesame-4757
       ctggtgcgtcgagcacttctcttgagcagcaacaatggcgg...
       >sesame-6933
       gaagtagatgctgtggtggtggaagcattcgacatatgcac...
    
       ...etc

Again, note that the gene identifier following the dash in the fasta headers (e.g. '6995' for header '>sesame-6995') needs to correspond to a gene identifier in your target file.

It's fine if your `outgroups.fasta` file contains additional sequences. When running the pipeline (see below) you can optionally provide one or more taxon names using the parameter `--external_outgroups <taxon1,taxon2,taxon3...>`, e.g. `--outgroups sesame`, and only these taxa will be included as outgroups. If this option isn't provided, all taxa/sequences in the `outgroups.fasta` file will be used. You can provide more than one outgroup taxon name using a comma-separated list, e.g. `--external_outgroups sesame,taxon2,taxon3` etc.

**NOTE:** at a minimum, you must provide either 'internal' ingroups via the `--internal_outgroups <taxon1,taxon2,taxon3...>` option, or a file of 'external' outgroup sequences via the `--external_outgroups_file outgroups.fasta` option.

## Running the pipeline on Linux

Please see the Wiki entry [Running on Linux][19].

## Running the pipeline on a Mac (macOS)

Please see the Wiki entry [Running on a Mac][17].

**NOTE:** Macs using the new Apple M1 chip are not yet supported.

## Running the pipeline on a PC (Windows)

Please see the Wiki entry [Running on a PC][18].

### Nextflow pipeline options and parameters

Example run command:

    nextflow run paragone.nf \
    -c paragone.config \
    -profile slurm_singularity \
    -resume \
    --gene_fasta_directory 11_paralogs \
    --external_outgroups_file outgroups.fasta \
    --outgroups sesame \
    --internal_outgroups <taxon1,taxon2,taxon3> \
    --mo --mi --rt 


See section [Pipeline parameters and options](#pipeline-parameters-and-options) for a full explanation of available parameters and flags. The required parameters are:

```
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
```

## Output folders and files

After running the pipeline, output can be found in the folder `results` (unless you have changed the name of the default output folder using the `--outdir <name>` parameter). This will consist of 28 subfolders. If you're just after the aligned `*.fasta` files for each of your target genes as output by each of the paralogy resolution methods, the three main output folders of interest are probably:

- `23_MO_final_alignments`
- `24_MI_final_alignments`
- `25_RT_final_alignments`

For a full explanation of output folders and files, please see the Wiki entry [Output folders and files][20].

## General information

For details on adapting the pipeline to run on local and HPC computing resources, see [here][22].

## Pipeline parameters and options

```
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
```

Please see the Wiki entry [Additional pipeline features and details][22] for further explanation of the parameters above, and general pipeline functionality.

## General notes

e.g.

- Manually reviewing trees from a preliminary run to select appropriate cut-off values for tree pruning
- Caveats of using these paralogy resolution approaches - relying largely on the fidelity of single-gene trees. Some loci will be better than others (i.e. short, low phylogenetic signal)

etc.

## Changelog

*07 July 2023*

- Add version to `paragone.nf` script;  v1.0.0
- Refactor `paragone.nf` for to use the Python package [`ParaGone`](https://github.com/chrisjackson-pellicle/ParaGone).
- New Singularity container with `ParaGone` installed.
- Added a `conda` and `conda_slurm` profile to the `paragone.config` file. This allows the pipeline to be run using conda packages rather than the Singularity container. The corresponding conda environment is created in the Nextflow `work` directory.

*28 November 2022*

- Change repository name from `Yang-and-Smith-paralogy-resolution` to `paragone-nf`.
- Update the required container name from `hybpiper-yang-and-smith-rbgv.sif` to `hybpiper-paragone.sif`.
- Arbitrarily resolve any polytomies produced when using FastTreeMP.
- If there are paralogs present in the outgroups, select a single representative sequence (i.e. most distant from a sample of ingroup taxa).

*09 November 2021* 

- Nextflow script updated to DSL2. **NOTE:** Nextflow version >= 21.04.1 is now required!
- Singularity *.def file updated to use Ubuntu 20.04, and to install Muscle and FastTreeMP.

[1]: http://trimal.cgenomics.org/ "Link to trimal website"
[2]: https://bmcecolevol.biomedcentral.com/articles/10.1186/s12862-019-1350-2 "Link to HmmCleaner manuscript"
[3]: http://www.iqtree.org/ "Link to IQtree website"
[4]: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4209138/ "Link to the Yang and Smith 2014 manuscript"
[5]: https://mafft.cbrc.jp/alignment/software/ "Link to mafft website"
[6]: https://www.ebi.ac.uk/seqdb/confluence/display/JDSAT/Clustal+Omega+Help+and+Documentation "Link to Clustal Omega website"
[7]: https://github.com/chrisjackson-pellicle/hybpiper-nf/wiki/Output-folders-and-files "Link to hybpiper-nf Wiki entry"
[8]: https://github.com/chrisjackson-pellicle/hybpiper-nf/wiki/Additional-pipeline-features-and-details#detection-of-putative-chimeric-contigs "Link to hybpiper-nf Wiki entry"
[11]: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8677558/ "Link to Yang 2022 manuscript"
[12]: https://bitbucket.org/dfmoralesb/target_enrichment_orthology/src/master/ "Link to Yang and Smith Bitbucket"
[13]:https://sylabs.io/docs/ "Link to Singularity website"
[14]:https://www.nextflow.io/ "Link to Nextflow website"
[15]:https://biopython.org/ "Link to BioPython website"
[16]: https://github.com/chrisjackson-pellicle/hybpiper-nf "Link to hybpiper-nf github"
[17]: https://github.com/chrisjackson-pellicle/paragone-nf/wiki/Running-on-a-Mac-(macOS)-with-Vagrant "Link to Running-on-a-Mac Wiki entry"
[18]: https://github.com/chrisjackson-pellicle/paragone-nf/wiki/Running-on-a-PC-(Windows)-with-Vagrant "Link to Running-on-a-PC Wiki entry"
[19]: https://github.com/chrisjackson-pellicle/paragone-nf/wiki/Running-on-Linux "Link to Running-on-Linux Wiki entry"
[20]: https://github.com/chrisjackson-pellicle/paragone-nf/wiki/Output-folders-and-files "Link to Output-folders-and-files Wiki entry"
[21]: https://github.com/chrisjackson-pellicle/paragone-nf/wiki/Additional-pipeline-features-and-details "Link to Additional-pipeline-features-and-details Wiki entry"
[22]: https://github.com/chrisjackson-pellicle/paragone-nf/wiki/Additional-pipeline-features-and-details#managing-computing-resources "Link to managing-computing-resources Wiki section"
[23]: http://www.microbesonline.org/fasttree/#OpenMP "Link to FastTree documentation"
[24]: https://www.drive5.com/muscle/downloads.htm "Link to MUSCLE download and documentation"
