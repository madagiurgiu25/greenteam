


cat gencode_lnc_hg38.bed | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$9"@"$10"@"$14"@"$15"\t0\t"$4}' > gencode_hg38_short.bed
cat gencode_lnc_mm10.bed | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$9"@"$10"@"$14"@"$15"\t0\t"$4}' > gencode_mm10_short.bed

giurgiu@bioclient2:~> cat noncode_hg38_short.bed| awk -F"\t" '{if ($6 != "") print $0}' > tmp
giurgiu@bioclient2:~> cat tmp > noncode_hg38_short.bed
giurgiu@bioclient2:~> cat noncode_mm10_short.bed| awk -F"\t" '{if ($6 != "") print $0}' > tmp
giurgiu@bioclient2:~> cat tmp > noncode_mm10_short.bed


cat noncode_lnc_hg38.bed | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$9"@"$10"@"$14"@"$15"\t0\t"$4}' > noncode_hg38_short.bed
cat noncode_lnc_mm10.bed | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$9"@"$10"@"$14"@"$15"\t0\t"$4}' > noncode_mm10_short.bed

/home/proj/biosoft/software/bedtools-2.17.0/bin/bedtools sort -i lncrnadb_hg38_short.bed >  lncrnadb_hg38_short_sorted.bed
/home/proj/biosoft/software/bedtools-2.17.0/bin/bedtools sort -i lncrnadb_mm10_short.bed >  lncrnadb_mm10_short_sorted.bed
/home/proj/biosoft/software/bedtools-2.17.0/bin/bedtools sort -i lncipedia_hg38_short.bed >  lncipedia_hg38_short_sorted.bed
/home/proj/biosoft/software/bedtools-2.17.0/bin/bedtools sort -i noncode_hg38_short.bed >  noncode_hg38_short_sorted.bed
/home/proj/biosoft/software/bedtools-2.17.0/bin/bedtools sort -i noncode_mm10_short.bed >  noncode_mm10_short_sorted.bed
/home/proj/biosoft/software/bedtools-2.17.0/bin/bedtools sort -i gencode_hg38_short.bed >  gencode_hg38_short_sorted.bed
/home/proj/biosoft/software/bedtools-2.17.0/bin/bedtools sort -i gencode_mm10_short.bed >  gencode_mm10_short_sorted.bed


############ find overlaps

/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/Noncoding/softwares/bedtools2/bin/bedtools intersect -wa -wb \
-a gencode_hg38_short_sorted.bed \
-b noncode_hg38_short_sorted.bed lncipedia_hg38_short_sorted.bed lncrnadb_hg38_short_sorted.bed \
-names noncode lncipedia lncrnadb \
-s -sorted -f 0.9 -r | uniq > overlapExons/gencode_hg38_overlappAll.bed


/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/Noncoding/softwares/bedtools2/bin/bedtools intersect -wa -wb \
-a noncode_hg38_short_sorted.bed \
-b lncipedia_hg38_short_sorted.bed lncrnadb_hg38_short_sorted.bed gencode_hg38_short_sorted.bed \
-names lncipedia lncrnadb gencode \
-s -sorted -f 0.9 -r | uniq > overlapExons/noncode_hg38_overlappAll.bed


/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/Noncoding/softwares/bedtools2/bin/bedtools intersect -wa -wb \
-a lncipedia_hg38_short_sorted.bed \
-b lncrnadb_hg38_short_sorted.bed gencode_hg38_short_sorted.bed noncode_hg38_short_sorted.bed \
-names lncrnadb gencode noncode \
-s -sorted -f 0.9 -r | uniq > overlapExons/lncipedia_hg38_overlappAll.bed

/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/Noncoding/softwares/bedtools2/bin/bedtools intersect -wa -wb \
-a lncrnadb_hg38_short_sorted.bed \
-b gencode_hg38_short_sorted.bed noncode_hg38_short_sorted.bed lncipedia_hg38_short_sorted.bed \
-names gencode noncode lncipedia \
-s -sorted -f 0.9 -r | uniq > overlapExons/lncrnadb_hg38_overlappAll.bed

cat overlapExons/gencode_hg38_overlappAll.bed | awk -F"\t" '{print $4"\t"$7"\t"$11}' | sort --parallel=4 -k1,1 -k3,3 | uniq -c > overlapExons/gencode_hg38_overlappAll_exons.bed

cat overlapExons/noncode_hg38_overlappAll.bed | awk -F"\t" '{print $4"\t"$7"\t"$11}' | sort --parallel=4 -k1,1 -k3,3 | uniq -c > overlapExons/noncode_hg38_overlappAll_exons.bed

cat overlapExons/lncrnadb_hg38_overlappAll.bed | awk -F"\t" '{print $4"\t"$7"\t"$11}' | sort --parallel=4 -k1,1 -k3,3 | uniq -c > overlapExons/lncrnadb_hg38_overlappAll_exons.bed


cat overlapExons/lncipedia_hg38_overlappAll.bed | awk -F"\t" '{print $4"\t"$7"\t"$11}' | sort --parallel=4 -k1,1 -k3,3 | uniq -c > overlapExons/lncipedia_hg38_overlappAll_exons.bed

tail -n+2 gencode_hg38_short.bed | awk '{print $4}' | sort | uniq -c > overlapExons/gencode_hg38_short_exons.bed
tail -n+2 noncode_hg38_short.bed | awk '{print $4}' | sort | uniq -c > overlapExons/noncode_hg38_short_exons.bed
tail -n+2 lncipedia_hg38_short.bed | awk '{print $4}' | sort | uniq -c > overlapExons/lncipedia_hg38_short_exons.bed
tail -n+2 lncrnadb_hg38_short.bed | awk '{print $4}' | sort | uniq -c > overlapExons/lncrnadb_hg38_short_exons.bed

