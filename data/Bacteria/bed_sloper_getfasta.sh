#step 4

SAVEIFS=$IFS
IFS=$(echo -en "\n\b")
for dir in $(ls -d */)
do
    # go through each dir
    cd $dir
    # go through each .bed file
    for file in *.bed
    do
        # check that file name has ".bed"
        if [[ $file == *\.[0-9].bed ]]
        then
            # append "100bp" to the $file
            # and save it as $file_100bp
            file_100bp=${file%.bed}_100bp.bed
            refgenome=$(echo $file | cut -d "_" -f 3 | cut -d "." -f 1,2)
            # if refgenome has "NC" in it change it to "NC_"
            if [[ $refgenome == *NC* ]]
            then
                refgenome=$(echo $refgenome | sed 's/NC/NC_/')
            fi
            # get the file name of $ref_genome at the beginning and end with genome_length.txt and save it as $file_genome_length
            file_genome_length=${refgenome}-ref.fasta-genome_length.txt # this should be .fai but that generates small trailing terminator which with this genome length file is not present
            # bed slop
            bedtools slop -i $file -g $file_genome_length -b 50 > $file_100bp
            # file name without the ".bed"
            file_no_ext_fasta=${file_100bp%.bed}.fasta
            # get the file name of *genome_length.txt and save it as $file_genome_length
            file_genome_fasta=${refgenome}-ref.fasta
            # bed getfasta
            bedtools getfasta -s -fi $file_genome_fasta -bed $file_100bp > $file_no_ext_fasta
        fi
    done
    # go back to the parent directory
    cd ..
done
IFS=$SAVEIFS