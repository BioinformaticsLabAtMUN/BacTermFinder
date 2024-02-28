# This script is used to deduplicate the regions from bed files

touch toKeep.bed

for i in 2 4 15 21 24 25 35
do
    echo "Processing sample $i ##########################################"
    # sort the bed file
    # sort -k 1,1 -k2,2n $i_*.bed > sorted_$i.bed
    # intersect with the other bed files
    for j in 2 4 15 21 24 25 35
    do
    # if the file is not the same as the current file
        if [ $i -ne $j ]
        then
            echo "Intersecting with sample $j"
            # file name has the index i at the beggining followed by underscore and any character can happen and ending by .bed
            namesa=$(ls $i* | grep -P "^$i+_\K.*\.bed")
            namesb=$(ls $j* | grep -P "^$j+_\K.*\.bed")
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
done

# bedtools merge -s -d -60 -c 4,5,6  -o distinct -i toKeep_sorted.bed | wc -l
# 3512

# bedtools merge -s -d -90 -c 4,5,6  -o distinct -i toKeep_sorted.bed | wc -l
# bedtools merge -s -d -90 -c 4,5,6  -o distinct -i toKeep_sorted.bed | awk -F"\t" '{if ($3-$2 > 101) {OFS="\t"; print $1,int(($3-$2)/2)+$2-50,int(($3-$2)/2)+$2+51,$4,$5,$6} else {print}}' > toKeep_merged.bed 
# 4139

# bedtools merge -s -d -101 -c 4,5,6  -o distinct -i toKeep_sorted.bed | wc -l
# 5992

# Number of genes in Ecoli is around 4,500 

# When i swaped a , b for -wo ones, the number of regions increased

# testing just merge them all
# for i in {102..122..1}; do echo -n $i ' '  ; bedtools merge -s -d -61 -c 4,5,6  -o distinct -i together_sorted.txt | awk -F"\t" '{if ($3-$2 < '$i') {print}}' | wc -l; done
# 102  2609
# 103  2934
# 104  3057
# 105  3122
# 106  3174
# 107  3209
# 108  3235
# 109  3253
# 110  3277
# 111  3293
# 112  3307
# 113  3312
# 114  3315
# 115  3316
# 116  3321
# 117  3327
# 118  3332
# 119  3343
# 120  3353
# 121  3359
# 122  3363