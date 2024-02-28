#!/bin/bash

# step 5 ( this file also has getfasta for neg bed files. the output are 3 files
# 1. 100bp bed file with 10x of the original bed file
# 2. 100bp bed file with 10x of the original bed file and shuffled (neg bed)
# 3. 100bp bed file with shuffled (neg fasta))

# for to go through each dir ( dir's name can have spaces )
SAVEIFS=$IFS
IFS=$(echo -en "\n\b")
for dir in $(ls -d */)
do
    echo "Processing $dir file..."
    # go through each dir (dir's name can have sapces so we use quotes)
    cd "$dir"
    # go through each .bed file
    for file in *.bed
    do
        # check that file name has ".bed"
        if [[ $file == *100bp.bed ]]
        then
            # get the name of the file
            refgenome=$(echo $file | cut -d "_" -f 3 | cut -d "." -f 1,2)
            # if refgenome has "NC" in it change it to "NC_"
            if [[ $refgenome == *NC* ]]
            then
                refgenome=$(echo $refgenome | sed 's/NC/NC_/')
            fi
            # get the file name of $ref_genome at the beginning and end with genome_length.txt and save it as $file_genome_length
            file_genome_length=${refgenome}-ref.fasta.fai
            file_neg_10x_100bp=${file%.bed}_10x_neg.bed # shuffles output (neg)
            file_10x_100bp=${file%.bed}_10x.bed # 10x of normal bed
            touch $file_10x_100bp
            for i in {1..10} # repeat original bed 10 times and append to new file
            do 
                cat $file >> "${file_10x_100bp}"
            done

            file_neg_100bp=${file%.bed}_100bp_neg.bed

            rm *5x* # remove 5x files

            bedtools shuffle -i $file -g $file_genome_length -excl $file -seed 42  > $file_neg_100bp 
            bedtools shuffle -i $file_10x_100bp -g $file_genome_length -excl $file_10x_100bp  -f 0.2 -seed 42  > $file_neg_10x_100bp 

            # file name without the ".bed"
            file_no_ext_fasta=${file_neg_100bp%.bed}.fasta
            file_no_ext_fasta_10x=${file_neg_10x_100bp%.bed}.fasta

            # get the file name of *genome_length.txt and save it as $file_genome_length
            file_genome_fasta=${refgenome}-ref.fasta
            # bed getfasta
            bedtools getfasta -s -fi $file_genome_fasta -bed $file_neg_100bp > $file_no_ext_fasta 
            bedtools getfasta -s -fi $file_genome_fasta -bed $file_neg_10x_100bp > $file_no_ext_fasta_10x 

        fi
    done
    # go back to the parent directory
    cd ..
done
IFS=$SAVEIFS