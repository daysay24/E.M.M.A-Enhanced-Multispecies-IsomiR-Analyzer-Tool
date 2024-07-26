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

1.  Put all outputs of isomiR-SEA in one folder: In this folder, there could be one or multiple subfolders, each represents a group (treatment, stage, time points, etc.) which contains one or multiple replicates.

    For example: We provided an example folder /test/mouse/isomiR-SEA_outputs. There are 2 time points (D0 and 18hr). Each group has two replicates: D0 includes D0_rpt1.txt and D0_rpt2.txt and 18hr includes 18hr_rpt1.txt and 18hr_rpt2.txt.

    ```
    ./test
        /mouse
            /isomiR-SEA_outputs
                /18hr
                    /18hr_rpt1.txt
                    /18hr_rpt2.txt
                /D0
                    /D0_rpt1.txt
                    /D0_rpt2.txt
    ```

    Then in the ./run.sh file, set path_raw_output_folder = <path to that folder>
    For example:

        ```
        path_raw_output_folder=./test/mouse/isomiR-SEA_outputs
        ```

2.  Set up genomic data: In order to get the precursor sequences extended from both sides of miRNAs, we need:

    - Genome sequence: (1) built-in from UCSC or (2) fasta file from users.
    - miRNA annotation: (1) gff file from miRBase or (2) an excel file from users which must have _chr_, _name_, _start_, _end_, _strand_ columns.

    Hence, there are four possible combinations. We provided sample files for Case 1 and Case 4 (download genome fasta for Schistosoma japonicum [here](https://parasite.wormbase.org/Schistosoma_japonicum_prjna724792/Info/Index), rename to _genome.fa_, save into ./test/s.jap/). Here are detailed setting up steps for each case:

    **Case 1:** Built-in genome from UCSC and gff file from miRBase.

    - Put the gff file from miRBase in a folder.

      For example:

      ```
      ./test
          /mouse
              /mmu.gff3
      ```

    - Modify /run.sh.

      For example:

      ```bash
      path_coords_file=./test/mouse/mmu.gff3
      ...
      species='mm10'
      is_mirbase_gff=True
      is_built_in_genome=True
      ```

    **Case 2:** Built-in genome from UCSC and an excel file from users.

    - Put the gff file from miRBase in a folder.

      **This must be an excel file with _chr_, _name_, _start_, _end_, _strand_ columns and coordinates are 1-based.**

      For example:

      ```
      ./test
          /mouse
              /mmu.xlsx
      ```

    - Modify ./run.sh.

      For example:

      ```bash
      path_coords_file=./test/mouse/mmu.xlsx
      ...
      species='mm10'
      is_mirbase_gff=False
      is_built_in_genome=True
      ```

    **Case 3:** Fasta file from users and gff file from miRBase.

    - Put the fasta genome from users and gff file from miRBase in a folder.

      For example:

      ```
      /test
        /mouse
            /genome.fa
            /mmu.gff3
      ```

    - Modify /run.sh.

      For example:

      ```bash
      path_genomic_file=./test/mouse/genome.fa
      path_coords_file=./test/mouse/mmu.gff3
      ...
      species=None
      is_mirbase_gff=True
      is_built_in_genome=False
      ```

    **Case 4:** Fasta file from users and an excel file from users.

    - Put the fasta genome and an excel file from users in a folder.

      For example:

      ```
      /test
        /s.jap
          /genome.fa
          /s.jap.xlsx
      ```

    - Modify ./run.sh.

      For example:

      ```bash
      path_genomic_file=./test/s.jap/genome.fa
      path_coords_file=./test/s.jap/s.jap.xlsx
      ...
      species=None
      is_mirbase_gff=False
      is_built_in_genome=False
      ```

3.  Run ./run.sh file

    ```bash
    bash run.sh
    ```

**NOTE**: We are using API (= sending http request) to get the precursor sequence from built-in genome, so running _3_generate_precursor.py_ takes a while (15 mins).
