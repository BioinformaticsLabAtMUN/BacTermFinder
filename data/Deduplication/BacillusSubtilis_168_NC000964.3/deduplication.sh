# This script is used to deduplicate the regions from bed files

touch toKeep.bed

# sort
for i in 5 7 19 24 26 35
do
    namesa=$(ls $i* | grep -P "^$i+_\K.*\.bed")
    # append sorted to the end of the file name
    sorted_name=$(echo $namesa | sed 's/.bed/_sorted.bed/g')
    # sort the bed file
    sort -k 1,1 -k2,2n $namesa -o $sorted_name
done

for i in 5 7 19 24 26 35
do
    echo "Processing sample $i ##########################################"
    # intersect with the other bed files
    for j in 5 7 19 24 26 35
    do
    # if the file is not the same as the current file
        if [ $i -ne $j ]
        then
            echo "Intersecting with sample $j"
            # file name has the index i at the beggining followed by underscore and any character can happen and ending by .bed
            namesa=$(ls $i* | grep -P "^$i+_\K.*_sorted\.bed")
            namesb=$(ls $j* | grep -P "^$j+_\K.*_sorted\.bed")
            # intersect the files
            bedtools intersect -s -v -f 0.5 -r -b $namesa -a $namesb >> toKeep.bed
            bedtools intersect -s -v -f 0.5 -r -a $namesa -b $namesb >> toKeep.bed
            # intersect the files for perfectly aligned terminators
            bedtools intersect -s -wo -f 0.5 -r -b $namesa -a $namesb | awk -F "\t" '$13 == 101 {print}{next}' | cut -f1-6 >> toKeep.bed
            bedtools intersect -s -wo -f 0.5 -r -a $namesa -b $namesb | awk -F "\t" '$13 == 101 {print}{next}' | cut -f1-6 >> toKeep.bed
            # intersect the files for averaging between not perfectly aligned terminators
            bedtools intersect -s -wo -f 0.5 -r -b $namesa -a $namesb | awk -F "\t" '$13 != 101 {OFS="\t";print $1,int(($2+$8)/2), int(($3+$9)/2),$4,$5,$6} {next}'  >> toKeep.bed
            bedtools intersect -s -wo -f 0.5 -r -a $namesa -b $namesb | awk -F "\t" '$13 != 101 {OFS="\t";print $1,int(($2+$8)/2), int(($3+$9)/2),$4,$5,$6} {next}'  >> toKeep.bed
        fi
    done
    # merge the bed file
    # bedtools merge -s -d -51 -c 4,5,6  -o distinct -i toKeep.bed > merged_$i.bed
    # remove the temporary files
    # rm sorted_$i.bed
    # rm toKeep.bed
done


# bedtools merge -s -d -101 -c 4,5,6  -o distinct -i toKeep_sorted.bed | wc -l ( even with this, we had some terms that had more than 101 length)
# 7882

# bedtools merge -s -d -90 -c 4,5,6  -o distinct -i toKeep_sorted.bed | wc -l
#  bedtools merge -s -d -90 -c 4,5,6  -o distinct -i toKeep_sorted.bed | awk -F"\t" '{if ($3-$2 > 101) {OFS="\t"; print $1,int(($3-$2)/2)+$2-50,int(($3-$2)/2)+$2+51,$4,$5,$6} else {print}}' > toKeep_merged.bed 
# 4872

# bedtools merge -s -d -60 -c 4,5,6  -o distinct -i toKeep_sorted.bed | wc -l
# 3502