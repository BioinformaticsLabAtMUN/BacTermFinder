
# step 3, getting a file ready for bedtools getfasta

SAVEIFS=$IFS
IFS=$(echo -en "\n\b")
for dir in $(ls -d */)
do
    # go through each dir
    cd $dir
    # go through each .bed file
    for file in *.fasta
    do
        # check that file name has ".bed"
        if [[ $file == *-ref.fasta ]]
        then
            # get the file name and cut off with delimiter "-"
            name=${file%-ref.fasta}
            # get the length of the reference genome
            length=$(tail -n+2 $file | wc -m)
            # create a file and put name followed by length tab delimited
            echo -e "$name\t$length" > $file-genome_length.txt
        fi
    done
    # go back to the parent directory
    cd ..
done
IFS=$SAVEIFS