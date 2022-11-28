***Note: this repository was previous called `Yang-and-Smith-paralogy-resolution`. It was renamed to  `paragone-nf`  on 28 November 2022, to maintain a consistent naming scheme with the Python package  `ParaGone` (in development).***

---

# paragone-nf: a hyb-seq paralogy resolution pipeline

## Original orthology inference manuscript, documentation and scripts

This pipeline makes use of the **paralogy resolution** (also described as **orthology inference**) approaches described and implemented by Yang and Smith 2014 [here][4]. These approaches have since been adapted for target capture datasets as described in the bioRxiv manuscript [here][11]. The original documentation and scripts can be found [here][12].

## paragone-nf: containerised and pipelined using Singularity and Nextflow

To simplify running the Yang and Smith paralogy resolution methods on target capture data, I’ve provided a [Singularity][13] container based on the Linux distribution Ubuntu 20.04, containing the original scripts required to run the pipeline (including modifications, additions and bug fixes, see below for details), as well as additional new scripts and all the dependencies ([IQTree][3], [Clustal Omega][6], [MAFFT][5], [BioPython][15], [HMMCleaner][2], [trimal][1], [FastTreeMP][23], [MUSCLE][24]). The container is called `hybpiper-paragone.sif`.

To run the paralogy resolution pipeline using this container, I’ve provided a [Nextflow][14] script that uses the software in the Singularity container. This pipeline runs all steps with a single command. The pipeline script is called `paragone.nf`. It comes with an associated config file called `paragone.config`. The only input required is a folder containing `.fasta` files for each of your target-capture loci, including paralogs, and an optional `.fasta` file containing outgroup sequences (used by some of the paralogy resolution methods, see below). The number of parallel processes running at any time, as well as computing resources given to each process (e.g. number of CPUs, amount of RAM etc) can be user configured by modifying the provided config file. The pipeline can be run directly on your local computer, or on an HPC system submitting jobs via a scheduler (e.g. SLURM, PBS, etc).

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

    nextflow run paragone.nf -c paragone.config -profile slurm -resume --hybpiper_paralogs_directory 11_paralogs --external_outgroups_file outgroups.fasta --outgroups sesame --internal_outgroups taxon1,taxon2,taxon3>

See section [Pipeline parameters and options](#pipeline-parameters-and-options) for a full explanation of available parameters and flags. The required parameters are:

```
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
```

## Output folders and files

After running the pipeline, output can be found in the folder `results` (unless you have changed the name of the default output folder using the `--outdir <name>` parameter. This will consist of 23 subfolders. If you're just after the aligned .fasta files for each of your target genes as output by each of the paralogy resolution methods, the three main output folders of interest are probably:

- `18_alignments_stripped_names_MO_realigned`
- `20_alignments_stripped_names_RT_realigned`
- `22_alignments_stripped_names_MI_realigned`

For a full explanation of output folders and files, please see the Wiki entry [Output folders and files][20].

## General information

For details on adapting the pipeline to run on local and HPC computing resources, see [here][22].

## Pipeline parameters and options

```
    Usage:
    The typical command for running the pipeline is as follows:

    nextflow run paragone.nf \
    -c paragone.config \
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
                                  multi-threading (e.g. MAFFT, IQ-TREE). Default 
                                  is 1
      
      --no_stitched_contig        Use this flag if you are processing paralogs 
                                  from a run of HybPiper that used the 
                                  --no_stitched_contig flag. MAFFT alignments are 
                                  re-aligned using Clustal Omega, which can do a 
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

      --use_muscle                Use MUSCLE to align sequences instead of MAFFT. 
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
                                  minimum number of ingroup taxa. Default is 2
      
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
```

Please see the Wiki entry [Additional pipeline features and details][22] for further explanation of the parameters above, and general pipeline functionality.

## General notes

e.g.

- Manually reviewing trees from a preliminary run to select appropriate cut-off values for tree pruning
- Caveats of using these paralogy resolution approaches - relying largely on the fidelity of single-gene trees. Some loci will be better than others (i.e. short, low phylogenetic signal)

etc.

## Changelog

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
[11]: https://www.biorxiv.org/content/10.1101/2020.08.21.261925v2 "Link to Yang 2021 bioarchives manuscript"
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
