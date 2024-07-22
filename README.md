# How to run the program

## Authors

Emma Nguyen, Dayna Sais

## Description

We develop a program that takes miRNA sequencing reads as input files. These reads may belong to different replicates of different groups such as treatment, stages, timepoints, etc.

The outputs:

- Additional information added to existing output files of isomiR-SEA software, including: annotation, number of nt difference at each end, number of snps. Besides, we identified some inconsistencies in the way isomiR-SEA categorise isomiR type so we refine their categorisation.

- A set of visualisations that provide insights into the isomiR data:

  Graph 1: RPM and relative abundance of isomiRs vs miRNAs in different groups.

  Graph 2: RPM and unique tags of isomiR types (3p, 5p, both, canonical, others) in different groups.

  Graph 3: RPM and unique tags of isomiR types (addition/truncation + nt) in different groups.

  Graph 4: Proportion of templated vs untemplated at addition positions in different groups.

  Graph 5: Proportion of nucleotides (A, U, C, G) at addition positions in different groups.

## Install dependencies

Python: requests, pandas.
R: ggplot2, patchwork, dplyr, ggrepel.

## Run the program

1.  Put all outputs of isomiR-SEA in /0_isomiR-SEA_isomiRs folder: In this folder, there could be one or multiple subfolders, each represents a group (treatment, stage, time points, etc.) which contains one or multiple replicates.

    For example: We provided sample files in /data. There are 2 time points (14hr and 18hr). Each group has two replicates: 14hr includes 14hr_1.txt and 14hr_2.txt and 18hr includes 18hr_1.txt and 18hr_2.txt. The folder looks like this:

    ```
    /data/0_isomiR-SEA_isomiRs
            /14hr
                /14hr_1.txt
                /14hr_2.txt
            /18hr
                /18hr_1.txt
                /18hr_2.txt
    ```

2.  Set up genomic data: In order to get the precursor sequences extended from both sides of miRNAs, we need:

    - Genome sequence: (1) built-in from UCSC or (2) fasta file from users.
    - miRNA annotation: (1) gff file from miRBase or (2) an excel file from users which must have _chr_, _name_, _start_, _end_, _strand_ columns.

    Hence, there are four possible combinations. We provided sample file for case 1. Here are detailed setting up steps for mm10 in each case:

    **Case 1:** Built-in genome from UCSC and gff file from miRBase.

    - Put the gff file from miRBase in /data/0_genomic_data.

    ```
    /data/0_genomic_data
            /mmu.gff3
    ```

    - Modify /run.sh:

    ```bash
    path_coords_file=./data/0_genomic_data/mmu.gff3
    ...
    species='mm10'
    is_mirbase_gff=True
    is_built_in_genome=True
    ```

    **Case 2:** Built-in genome from UCSC and an excel file from users.

    - Put the gff file from miRBase in /data/0_genomic_data.

    **This must be an excel file with _chr_, _name_, _start_, _end_, _strand_ columns and coordinates are 1-based.**

    ```
    /data/0_genomic_data
            /mmu.xlsx
    ```

    - Modify /run.sh:

    ```bash
    path_coords_file=./data/0_genomic_data/mmu.xlsx
    ...
    species='mm10'
    is_mirbase_gff=False
    is_built_in_genome=True
    ```

    **Case 3:** Fasta file from users and gff file from miRBase.

    - Put the fasta genome from users and gff file from miRBase in /data/0_genomic_data

    ```
    /data/0_genomic_data
            /genome.fa
            /mmu.gff3
    ```

    - Modify /run.sh:

    ```bash
    path_genomic_file=./data/0_genomic_data/genome.fa
    path_coords_file=./data/0_genomic_data/mmu.gff3
    ...
    species=''
    is_mirbase_gff=True
    is_built_in_genome=False
    ```

    **Case 4:** Fasta file from users and an excel file from users.

    - Put the fasta genome and an excel file from users in /data/0_genomic_data

    ```
    /data/0_genomic_data
            /genome.fa
            /coords.xlsx
    ```

    - Modify /run.sh:

    ```bash
    path_genomic_file=./data/0_genomic_data/genome.fa
    path_coords_file=./data/0_genomic_data/coords.xlsx
    ...
    species=''
    is_mirbase_gff=False
    is_built_in_genome=False
    ```

3.  Run /run.sh file

    ```bash
    bash run.sh
    ```

**NOTE**: We are using API (= sending http request) to get the precursor sequence from built-in genome, so running _3_generate_precursor.py_ takes a while (15 mins). If it takes too long, just use the file 3_8_extended_precursor_seqs.csv by:

- Create /data/3_precursors folder
- Move 3_8_extended_precursor_seqs.csv to the /data/3_precursors folder
- Comment following lines in /run.sh like this:

```bash
# species='mm10'
# is_mirbase_gff=False
# is_built_in_genome=False
# python ./code/3_generate_precursor.py $path_summarised_output_folder $path_precursors_output_folder $path_genomic_file $path_coords_file $species $is_mirbase_gff $is_built_in_genome

```
