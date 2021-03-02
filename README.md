# paralogy_resolution_tutorial
GAP tutorial for containerised Yang and Smith paralogy resolution pipeline 



# NewTargets

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



