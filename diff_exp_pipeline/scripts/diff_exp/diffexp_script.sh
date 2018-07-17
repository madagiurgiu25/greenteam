
#bam="/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/daten/STAR"
#sample="M103"
#subset_size=4
#data="/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/daten/Ballgown_Supergenes_exclude20"
#wd="/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/daten/DiffExp"
#printGenes="TRUE"
#file="/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/daten/Ballgown_Supergenes_exclude20/test.txt"


sample=$1
subset_size=$2
data=$3
bam=$4
wd=$5
#printGenes=$6
#plots=$7
#transcriptExpression=$8
file=$6
file_genes=$7
anno=$8

ballgown="/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/Noncoding/greenteam/scripts/diffexp/ballgown_run.R"
r="/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/Noncoding/softwares/R-3.5.0/bin/Rscript"
cuffdiff="/home/proj/biosoft/software/cufflinks-2.2.1.Linux_x86_64/cuffdiff"

ballgownPath="${wd}/${sample}/Ballgown"
mkdir -p $ballgownPath

############################## Run Ballgown
echo "Run ballgown ... ${sample} ... ${subset_size} ...."
#"$r" "$ballgown" "$sample" "$subset_size" "$data" "$bam" "$ballgownPath" "TRUE" "TRUE" "TRUE" "$file" "$file_genes"
#newsize=$((subset_size-1))
#echo "Run ballgown subsets... ${sample} ... ${newsize} ...."
#"$r" "$ballgown" "$sample" "${newsize}" "$data" "$bam" "$ballgownPath" "FALSE" "FALSE" "TRUE" "$file" "$file_genes"
#
############################### Run Cuffidff

echo "Run Cuffdiff .... ${sample} ..."
condition1=`cat "${ballgownPath}/bams.txt" | awk '{print $1}'`
condition2=`cat "${ballgownPath}/bams.txt" | awk '{print $2}'`

pathCuffdiff="${wd}/${sample}/Cuffdiff"
mkdir -p "$pathCuffdiff"


condition1=`ls "${star_root}${sample}"*/star_sorted.bam | paste -sd "," -`
condition2=`ls "${star_root}${sample}"*/star_sorted.bam | paste -sd "," -`


$cuffdiff -o $pathCuffdiff -p 5 -c 5 -v --FDR 0.05 \
--library-type fr-firststrand \
"$anno" \
"$condition1" "$condition2"


#condition1=`ls "${star_root}"H103A*/star_sorted.bam | paste -sd "," -`
#condition2=`ls "${star_root}"H103B*/star_sorted.bam | paste -sd "," -`

#$cuffdiff -o "$pathCuffdiff" -p 5 -c 5 --FDR 0.05 -q \
#"$anno" \
#"$condition1" "$condition2"
