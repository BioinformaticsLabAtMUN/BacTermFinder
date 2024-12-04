# BacTermFinder: Bacteria-agnostic Comprehensive Terminator Finder using a CNN Ensemble

## Abstract 
Terminator is a region in the DNA that ends the transcription process. Finding bacterial terminators will lead to a better understanding of how bacteria's transcription works.  Currently, multiple tools are available for predicting bacterial terminators. However, most methods are specialized for certain bacteria or terminator types. In this work, we developed BacTermFinder, a tool that utilizes Convolutional Neural Networks (CNNs) with four different genomic representations trained on 41k bacterial terminators identified using RNA-seq technologies. Based on our results, BacTermFinder's recall score is  higher than that of the other four approaches we considered in our independent validation set of five different bacteria. Moreover, BacTermFinder's model identifies both types of terminators (intrinsic and factor-dependent) and even generalizes to archeal terminators. 

## How to run 
1. Create a virutal environment, we recommend [miniconda](https://docs.anaconda.com/miniconda/install/#),

```
conda create -n bactermfinder python==3.9
```

2. Then install the requirements via (you need to install numpy first for skbio package, then you can install other packages. If you are in Windows, you may need to install C++ compiler for skbio)

```
pip install numpy==1.23.0

pip install -r requirements.txt
```

3. After that you can run 

```bash
python genome_scan.py     YOUR_SEQ.fasta                 3 		                  out              10000              > log.out

|        	         |      your fasta file    |  step size sliding window  | output name | Feature gen batch size |   > log.out

```

The results would be in a file called `outsequence.fasta_mean.csv`

BacTermfinder will encode sequence to 4 different encodings (binary, ENAC, NCP, PS2) and run CNNs for each of them. After that, the mean of the predictions will be stored in output file  `outsequence.fasta_mean.csv`.

The log of prediction will be in log.out text file. 

The sliding windows would be in the `out_sequence.fasta_sliding_windows.csv`

## Threshold for different bacteria
We recommend to use the these threshold to classify terminators. High GC content bacteria genomes tend to have more factor-dependent termintors - which usually don't have strong motifs - and it's better to use less strict thresholds to find factor-dependent bacteria. One issue with less strict thresholds is that you could find false positive terminators. We got these numbers by maximizing the F-scores in our validation datasets.

|  Maximizing metric      | Threshold     |
| ----------------------- | ------------- |
| F2 - Less strict        |     0.13      |
| F1 - normal striction   |     0.3       |
| F0.5 - more strict      |     0.47      |


## Please cite it if you've used it!
Thank you very much for using our software, you can cite it as follows: 

BacTermFinder: A Comprehensive and General Bacterial Terminator Finder using a CNN Ensemble Seyed Mohammad Amin Taheri Ghahfarokhi, Lourdes Peña-Castillo. bioRxiv 2024.07.05.602086; doi: https://doi.org/10.1101/2024.07.05.602086 
