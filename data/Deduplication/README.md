# Steps for deduplication

## Step 0 
sort -k 1,1 -k2,2n 13_StreptomycesLividans_TK24_CP009124.1_100bp.bed > 13_sorted_StreptomycesLividans_TK24_CP009124.1_100bp.bed 

## step 1 & 2
bedtools intersect -s -v -f 0.5 -r -b 13_sorted_StreptomycesLividans_TK24_CP009124.1_100bp.bed -a 23_sorted_StreptomycesLividans_TK24_CP009124.1_100bp.bed > toKeep.bed
bedtools intersect -s -v -f 0.5 -r -a 13_sorted_StreptomycesLividans_TK24_CP009124.1_100bp.bed -b 23_sorted_StreptomycesLividans_TK24_CP009124.1_100bp.bed >> toKeep.bed

## step 3 for perfectly aligned terminators
bedtools intersect -s -wo -f 0.5 -r -b 13_sorted_StreptomycesLividans_TK24_CP009124.1_100bp.bed -a 23_sorted_StreptomycesLividans_TK24_CP009124.1_100bp.bed | awk -F "\t" '$13 == 101 {print}{next}' | cut -f1-6 >> toKeep.bed 

same above again switch a,bed

# step 4 for averaging between not perfectly aligned terminators
bedtools intersect -s -wo -f 0.5 -r -b 13_sorted_StreptomycesLividans_TK24_CP009124.1_100bp.bed -a 23_sorted_StreptomycesLividans_TK24_CP009124.1_100bp.bed | awk -F "\t" '$13 != 101 {OFS="\t";print $1,int(($2+$8)/2), int(($3+$9)/2),$4,$5,$6} {next}'  >> toKeep.bed 

same above again switch a,bed

# sort toKeep
sort -k 1,1 -k2,2n 

# Don't forget to merge
bedtools merge -s -d -90 -c 4,5,6  -o distinct -i toKeep_sorted.bed | awk -F"\t" '{if ($3-$2 > 101) {OFS="\t"; print $1,int(($3-$2)/2)+$2-50,int(($3-$2)/2)+$2+51,$4,$5,$6} else {print}}' > toKeep_merged.bed 

## What doesn't work

bedtools closest -s -d -b 13_sorted_StreptomycesLividans_TK24_CP009124.1_100bp.bed -a 23_sorted_StreptomycesLividans_TK24_CP009124.1_100bp.bed > closest.txt

awk and sort by numbers how many less than 49
awk  -F"\t" '{ if ($13 > 49) {print $1"\t"($2+$8)//2"\t"($3+$9)//2"\t"$10"\t"$11"\t"$12} else {pass}}' closest.txt | wc -l # 2198 (removes 325) ( ther are 650 pairs that their distance is less than 50) , 40 of them have a non-zero distance.

bedtools merge -s -d -51 -c 4,5,6  -o distinct -i together_sorted.txt | awk -F"\t" '{if ($3-$2 > 101) {print $3-$2}}' |sort -n| # 1986 (removes 537 duplication ) 68 of them have a non-zero distance.