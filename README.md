# BacTermFinder: Bacteria-agnostic Comprehensive Terminator Finder using a CNN Ensemble

## Abstract 
Terminator is a region in the DNA that ends the transcription process. Finding bacterial terminators will lead to a better understanding of how bacteria's transcription works.  Currently, multiple tools are available for predicting bacterial terminators. However, most methods are specialized for certain bacteria or terminator types. In this work, we developed BacTermFinder, a tool that utilizes Convolutional Neural Networks (CNNs) with four different genomic representations trained on 41k bacterial terminators identified using RNA-seq technologies. Based on our results, BacTermFinder's recall score is  higher than that of the other four approaches we considered in our independent validation set of five different bacteria. Moreover, BacTermFinder's model identifies both types of terminators (intrinsic and factor-dependent) and even generalizes to archeal terminators. 

## How to run 
1. Create a virutal environment, I prefere miniconda,

`conda create -n bactermfinder python==3.9`

2. Then install the requirements via 

`pip install -r requirements.txt`

3. After that you can run 

`|        	         |      your fasta file    |  step size sliding window  | output name | Feature gen batch size |   > log.out`

`python genome_scan.py     YOUR_SEQ.fasta                 3 		                  out              10000              > log.out`

The results would be in a file called `outsequence.fasta_mean.csv`

BacTermfinder will encode sequence to 4 different encodings (binary, ENAC, NCP, PS2) and run CNNs for each of them. After that, the mean of the predictions will be stored in output file  `outsequence.fasta_mean.csv`.

The log of prediction will be in log.out text file. 

The sliding windows would be in the `out_sequence.fasta_sliding_windows.csv`

## Please cite it if you've used it!
Thank you very much for using our software, you can cite it as follows: 

BacTermFinder: A Comprehensive and General Bacterial Terminator Finder using a CNN Ensemble Seyed Mohammad Amin Taheri Ghahfarokhi, Lourdes Peña-Castillo. bioRxiv 2024.07.05.602086; doi: https://doi.org/10.1101/2024.07.05.602086 
