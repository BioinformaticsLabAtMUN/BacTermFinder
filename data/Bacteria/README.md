# BacTermData: Experimentally Verified Bacteria Terminators

## Description of directories

In each directories, there are these files:

- The original file downloaded from the study additional files
- A `bedder.ipynb` which takes the downloaded file and does preprocessing on it and outputs a `.bed` file
- There is fasta file and gff files of the reference genome and the terminators of the study
- The naming convention is as follows
  - {BacteriaName}\_{Strain}-{Substrain}_{RefrenceGenome}_100bp
    - example: EscherichiaColi_K12-BW25113_CP009273.1_100bp
  - 10x_neg is the generated negative terminators that were randomly selected from the reference genome. See the paper and `.sh` files for more details.

## all_csv directory

There are CSV (Comma Separated Values) of terminators in that folder with some nucleotide frequency plots. 

## Technical step by step notes
1.  Found data from multiple sources

    1.  Term-seq paper citations 

    2.  Similar technologies like term-seq

    3.  Look into SRA - GEO databases for papers and data

2.  Downloaded their terminator files in different formats 

3.  For each data, create a bedder.ipynb in each study folder to create a bed file for TTS ( Transcription Termination Site)

    1.  Named the file {FullNameOfTheBacteria}_{strainCode-substrain}_{accessionNumber}.bed

4.  Used bed slope to widen the range to 100bp range (gene_file_creator_for_bed_slop.sh was needed to create a gene file for bed slope first and then bed_sloper_getfasta.sh)

    1.  {FullNameOfTheBacteria}_{strainCode-substrain}_{accessionNumber}_100bp.bed

5.  Used bed shuffle -excl  and get fasta to create neg files (randomly across the genome for 10x more data) neg_sample.sh

    1.  {FullNameOfTheBacteria}_{strainCode-substrain}_{accessionNumber}_100bp_100bp_neg.bed

    2.  {FullNameOfTheBacteria}_{strainCode-substrain}_{accessionNumber}_100bp_10x_neg.bed

6.  Used Esearch to get fasta files for each genome accession number (get_ref_genome.sh)

7.  Used get fasta -s to get the fasta file for positive and negative data (neg_sample.sh)

    1.  Same name with bed files but ends with .fasta instead of .bed

8.  Used underscore_remover.sh to change genome accession number from NC_XXXX to NCXXXX for easier managing file names (I don't know when i needed it)

9.  Used csv_generator.py to create CSV files from fasta files in each dir and have information about each sample as metadata. Information like study name, the position of the terminator start and end, reference genome, specie, strain, strand, ... . This file should be run multiple times with different inputs to get all_pos.csv and all_neg.csv

10. Deduplication (Details in the folder `Readme.md`)

11. With the vis.ipynb I created visualizations based on how log2fold of terminators versus random generation across the sequence. Which is like WebLogo but with lines instead of characters.

12. Used split_data.ipynb to combine all_pos and all_neg and also filter out the test files from training. I did a couple of EDA in that notebook, like duplication and whether we should delete them. The number of terminators per reference genome and per specie is also available in the ipynb.

    1.  Deleted some data that we deducted, it was low quality with the visualizations done in step 10

13. The train_data.csv was created and transformed into a special fasta format for iLearnPlus to create features. The transformation is done with data_format_transfer_ilearn.py

14. The df1x_new.fasta is generated and fed into an enhanced iLearnPlus script.

15. The enhanced iLearnPlus script is modified by me to handle parallel processing and batching of input data.