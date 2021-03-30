# HybSeq Paralogy Resolution Tutorial

## Original orthology inference manuscript, documentation and scripts

This pipeline makes use of the **paralogy resolution** (also described as **orthology inference**) approaches described and implemented by Yang and Smith 2014 [here][4]. These approaches have since been adapted for target capture datasets as described in the bioRxiv manuscript [here][11]. The original documentation and scripts can be found [here][12].


## Yang_and_Smith-RBGV: containerised and pipelined using Singularity and Nextflow

To simplify running the Yang and Smith paralogy resolution methods on target capture data, I’ve provided a [Singularity][13] container containing the Linux distribution Ubuntu 18.04, containing the original scripts required to run the pipeline (including modifications, additions and bug fixes, see below for details), as well additonal new scripts and all the dependencies ([IQTree][3], [Clustal Omega][6], [mafft][5], [BioPython][15], [HMMCleaner][2], [trimal][1]). The container is called `yang_and_smith.sif`.

To run the paralogy resolution pipeline using this container, I’ve provided a [Nextflow][14] script that uses the software in the Singularity container. This pipeline runs all steps with a single command. The pipeline script is called `yang-and-smith-rbgv-pipeline.nf`. It comes with an associated config file called `yang-and-smith-rbgv.config`. The only input required is a folder containing `.fasta` files for each of your target-capture loci, including paralogs, and a `.fasta` file containing outgroup sequences (used by some of the paralogy resolution methods, see below). The number of parallel processes running at any time, as well as computing resources given to each process (e.g. number of CPUs, amount of RAM etc) can be configured by the user by modifying the provided config file. The pipeline can be run directly on your local computer, and on an HPC system submitting jobs via a scheduler (e.g. SLURM, PBS, etc).


## Input data

### Paralog sequences

If you have used the Nextflow pipeline `hybpiper-rbgv-pipeline.nf` to run HybPiper, as described [here][16], your main `results` folder will contain the following subfolders:

 - `11_paralogs`
 - `12_paralogs_noChimeras`

See the HybPiper-RBGV wiki entry [Output folders and files][7] for a full description of the files in these output folders. Briefly, folder `11_paralogs` contains a fasta file for each gene in your HybPiper target file. Each fasta file contains the 'main' contig selected by HybPiper for each sample. Where HybPiper has detected putative paralog contigs, these sequences are also included; in such cases, the main contig has the fasta header suffix `.main`, whereas putative paralogs have the suffix `.0`, `.1` etc. Folder `12_paralogs_noChimeras` contains the same data, except putative chimeric contigs (see [here][8] for an explanation) have been removed.        

***Tutorial step 1:***

    Copy folder `11_paralogs` into your current working directory (i.e. the directory containing `yang-and-smith-rbgv-pipeline.nf` and `yang-and-smith-rbgv.config`.

### Outgroup sequences

Some of the paralogy resolution methods used in this pipeline require an outgroup sequence for each of your genes. These outgroup sequences can be provided in a fasta file (e.g. `outgroups.fasta`), with the same fasta header formatting as your HybPiper target file. For example, if you have used the Angiosperms353 target file, and you wish to use sequences from Sesame as your outgroup, your `outgroups.fasta` file might contain the following:

    >sesame-6995
    gtgggatatgaacaaaatccattgagcttgtattactgtta...
    >sesame-4757
    ctggtgcgtcgagcacttctcttgagcagcaacaatggcgg...
    >sesame-6933
    gaagtagatgctgtggtggtggaagcattcgacatatgcac...
    
    ...etc
    
Again, note that the gene identifier following the dash in the fasta headers (e.g. '6995' for header '>sesame-6995') needs to correspond to a gene identifier in your target file. 

It's fine if your `outgroups.fasta` file contains additional sequences. When running the pipeline (see below) you'll provide one or more taxon names using the parameter `--outgroups <taxon_name>`, e.g. `--outgroups sesame`, and only these taxa will be included as outgroups. You can provide more than one outgroup taxon name using an comma-separated list, e.g. `--outgroups sesame,taxon2,taxon3` etc.

***Tutorial step 2:***

    Copy the fasta file containing outgroup sequences to your current working directory.


## Running the pipeline on Linux

Please see the Wiki entry [Running on Linux][19].

## Running the pipeline on a Mac (macOS)

Please see the Wiki entry [Running on a Mac][17].

## Running the pipeline on a PC (Windows)

Please see the Wiki entry [Running on a PC][18].



### Nextflow pipeline options and parameters

See section [Pipeline parameters and options](#pipeline-parameters-and-options) for a full explanation of available parameters and flags. The required parameters are:

    --hybpiper_paralogs_directory <directory>    Path to folder containing HybPiper paralog fasta files)
    --outgroups_file <file>                      File containing fasta sequences of target genes
    --outgroups <taxon1,taxon2,taxon3...>        A comma-separated list of outgroup taxa to add, in order of 
                                                 preference

***Tutorial step 3:***

    Run the pipeline using the command:
    
    nextflow run yang-and-smith-rbgv-pipeline.nf -c yang-and-smith-rbgv.config -profile slurm -resume --hybpiper_paralogs_directory 11_paralogs --outgroups_file outgroups.fasta --outgroups sesame

## Output folders and files

After running the pipeline, output can be found in the folder `results` (unless you have changed the name of the default output folder using the `--outdir <name>` parameter. This will consist of 23 subfolders. If you're just after the aligned .fasta files for each of your target genes as output by each of the paralogy resolution methods, the three main output folders of interest are probably:

- `19_alignments_stripped_names_MO_realigned`
- `21_alignments_stripped_names_RT_realigned`
- `23_alignments_stripped_names_MI_realigned`


For a full explanation of output folders and files, please see the Wiki entry [Output folders and files][20].



## Post-pipeline analyses

e.g.

- Concanated IQtree
- Astral

etc.

## General information

For details on adapting the pipeline to run on local and HPC computing resources, see [here][].

## Bug fixes and changes (WIP)

Please see the Wiki entry [Bug fixes and changes][].



## Issues still to deal with (WIP)

Please see the Wiki entry [Issues][].


## Pipeline parameters and options

Mandatory arguments:

      --hybpiper_paralogs_directory <directory>       Path to folder containing HybPiper paralog fasta files.
      --outgroups_file <file>                         File containing fasta sequences of outgroup sequences for each gene.
      --outgroups <taxon1,taxon2,taxon3...>           A comma-separated list of outgroup taxa to add, in order of 
                                                      preference.  
Optional arguments:

      -profile <profile>                              Configuration profile to use. Can use multiple (comma separated)
                                                      Available: standard (default), slurm
      --pool <int>                                    Number of threads for the Python multiprocessing pool. Used in e.g. alignments and tree-building steps. 
                                                      Default is 1, so e.g. one alignment will be run at a time during alignment steps.
      --threads <int>                                 Number of threads per multiprocessing pool instance. Used for programs that support multi-threading (e.g. mafft,
                                                      IQTree). Default is 1.
      --no_supercontigs                               Use this flag if you are processing paralogs from a run of HybPiper that used the --nosupercontigs flag. Mafft alignments with re-aligned using clustal omega, which can do a better job in these cases. Default is off.
      --process_02_trim_bad_ends_cutoff <int>         Number of bases either side of an internal gap that much match the reference before trimming stops.Default is 5.
      --process_02_trim_bad_ends_size <int>           Number of continuous internal gap positions for an internal gap to be investigated. Default is 15.
      --skip_process_02_trim_bad_ends                 Skips the step trimming the ends of internal gaps.
      --process_04_trim_tips_relative_cutoff <float>  When pruning long tips during the tree QC stage, provide a branch length for the maximum imbalance between sister tips allowed. Default is 0.2.
      --process_04_trim_tips_absolute_cutoff <float>  When pruning long tips during the tree QC stage, provide a branch length for the maximum allowed tip branch length. Default is 0.4.
      --process_06_branch_length_cutoff <float>       When pruning long internal branches (putative deep paralogs) during the tree QC stage, provide a branch length for the maximum allowed internal branch length. Default is 0.3.
      --process_06_minimum_taxa <int>                 After the final tree-pruning step prior to paralogy resolution, only retain trees with a minimum number of taxa remaining. Default is 3. 
      --process_09_prune_paralog_MO_minimum_taxa <int>
                                                      For the MO method, only process trees with a minimum number of taxa. Default is 2.
      --process_10_prune_paralogs_RT_minimum_ingroup_taxa <int>
                                                      For the RT method, only process trees with a minumum number of ingroup taxa. Default is 2. 
      --process_11_prune_paralogs_MI_relative_tip_cutoff <float>
                                                      Default is 0.2.
      --process_11_prune_paralogs_MI_absolute_tip_cutoff <float>
                                                      Default is 0.4.
      --process_11_prune_paralogs_MI_minimum_taxa <int>   
                                                      Default is 2.

Please see the Wiki entry [Additional pipeline features and details][] for further explanation of the parameters above, and general pipeline functionality.


## Managing computing resources

Please see [here][10] for instructions on how to customise the file `yang-and-smith-rbgv.config` for your particualar computing resources. Note that the Nextflow processes listed in the `yang-and-smith-rbgv.config` file will differ to those described in the instructions, but the same principles apply.  

### General notes

e.g.

- Manually reviewing trees from a preliminary run to select appropriate cut-off values for tree pruning
- Caveats of using these paralogy resolution approaches - relying largely on the fidelity of single-gene trees. Some loci will be better than others (i.e. short, low phylogenetic signal)

etc.


[1]: http://trimal.cgenomics.org/ "Link to trimal website"
[2]: https://bmcecolevol.biomedcentral.com/articles/10.1186/s12862-019-1350-2 "Link to HmmCleaner manuscript"
[3]: http://www.iqtree.org/ "Link to IQtree website"
[4]: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4209138/ "Link to the Yang and Smith 2014 manuscript"
[5]: https://mafft.cbrc.jp/alignment/software/ "Link to mafft website"
[6]: https://www.ebi.ac.uk/seqdb/confluence/display/JDSAT/Clustal+Omega+Help+and+Documentation "Link to Clustal Omega website"
[7]: https://github.com/chrisjackson-pellicle/HybPiper-RBGV/wiki/Output-folders-and-files "Link to HybPiper-RBGV Wiki entry"
[8]: https://github.com/chrisjackson-pellicle/HybPiper-RBGV/wiki/Additional-pipeline-features-and-details#detection-of-putative-chimeric-contigs "Link to HybPiper-RBGV Wiki entry"
[9]: https://github.com/chrisjackson-pellicle/HybPiper-RBGV/wiki/Running-on-a-Mac-(macOS)-with-Vagrant "Link to Vagrant macOS instructions"
[10]: https://github.com/chrisjackson-pellicle/HybPiper-RBGV/wiki/Additional-pipeline-features-and-details#managing-computing-resources "Link to customising Nexflow config instructions"
[11]: https://www.biorxiv.org/content/10.1101/2020.08.21.261925v2 "Link to Yang 2021 bioarchives manuscript"
[12]: https://bitbucket.org/dfmoralesb/target_enrichment_orthology/src/master/ "Link to Yang and Smith Bitbucket"
[13]:https://sylabs.io/docs/ "Link to Singularity website"
[14]:https://www.nextflow.io/ "Link to Nextflow website"
[15]:https://biopython.org/ "Link to BioPython website"
[16]: https://github.com/chrisjackson-pellicle/HybPiper-RBGV "Link to HybPiper-RBGV github"
[17]: https://github.com/chrisjackson-pellicle/paralogy_resolution_tutorial/wiki/Running-on-a-Mac-(macOS)-with-Vagrant "Link to Running-on-a-Mac Wiki entry"
[18]: https://github.com/chrisjackson-pellicle/paralogy_resolution_tutorial/wiki/Running-on-a-PC-(Windows)-with-Vagrant "Link to Running-on-a-PC Wiki entry"
[19]: https://github.com/chrisjackson-pellicle/paralogy_resolution_tutorial/wiki/Running-on-Linux "Link to Running-on-Linux Wiki entry"
[20]: https://github.com/chrisjackson-pellicle/paralogy_resolution_tutorial/wiki/Output-folders-and-files "Link to Output-folders-and-files Wiki entry"
