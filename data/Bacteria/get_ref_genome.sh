# This script will go through all of the .bed files in different dirs and download the reference genome from NCBI 

# Step 2

# for handling spaces in file names
SAVEIFS=$IFS
IFS=$(echo -en "\n\b")

# go through all of the .bed files in the same directory
for dir in $(ls -d */)
do
    # go through each dir
    cd $dir
    # go through each .bed file
    for file in *.bed
    do
        # if regex matches
        if [[ $file = *100bp.bed ]]
        then
            # get the reference genome
            # from the .bed file
            ref_genome=$(head -n 1 $file | cut -f 1 )
            # if the ref_genome only had 1 underscores
            ref_genome=$(echo $ref_genome | cut -d "_" -f -2)
            # get rid of trailing "_0" at the end of some ref_genom     
            ref_genome=$(echo $ref_genome | sed 's/_0$//')
            
            # esearch -db nuccore -query "$ref_genome" < /dev/null | efetch -format fasta > $ref_genome-ref.fasta & # commenting for now

            f=$(esearch -db assembly -query "$ref_genome" | esummary | xtract -pattern DocumentSummary -element FtpPath_RefSeq) && wget -O - $f/*gff.gz > $ref_genome.gff.gz &
            
            # sleep 1 second
            sleep 1
            
        fi

    done
    # go back to the parent directory
    cd ..
done
IFS=$SAVEIFS # restore $IFS