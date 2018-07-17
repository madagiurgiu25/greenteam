 #!/bin/sh

path=$1
condition1=$2
condition2=$3
anno=$4

####################################################################################################
################################################################################################
######################### Diff exp using Cufflinks #############################################
################################################################################################
cuffdiff=/home/proj/biosoft/software/cufflinks-2.2.1.Linux_x86_64/cuffdiff

echo "$condition1"
echo "$condition2"

$cuffdiff -o $path -p 5 -c 5 -v --FDR 0.05 \
--library-type fr-firststrand \
"$anno" \
"$condition1" "$condition2"







