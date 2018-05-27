{\rtf1\ansi\ansicpg1252\cocoartf1404\cocoasubrtf470
{\fonttbl\f0\fnil\fcharset0 Menlo-Regular;}
{\colortbl;\red255\green255\blue255;}
\paperw11900\paperh16840\margl1440\margr1440\vieww23260\viewh8400\viewkind0
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural\partightenfactor0

\f0\fs22 \cf0 \CocoaLigature0 #!/bin/bash\
\
cat count_exons_gencode_lnc_hg38.txt | awk -F" " '\{print $1"\\t"$2\}' > count_exons_gencode_lnc_hg38_formatted.txt \
cat count_exons_lncipedia_lnc_hg38.txt | awk -F" " '\{print $1"\\t"$2\}' > count_exons_lncipedia_lnc_hg38_formatted.txt \
cat count_exons_lncrnadb_lnc_hg38.txt | awk -F" " '\{print $1"\\t"$2\}' > count_exons_lncrnadb_lnc_hg38_formatted.txt \
cat count_exons_noncode_lnc_hg38.txt | awk -F" " '\{print $1"\\t"$2\}' > count_exons_noncode_lnc_hg38_formatted.txt \
cat count_exons_gencode_lnc_mm10.txt | awk -F" " '\{print $1"\\t"$2\}' > count_exons_gencode_lnc_mm10_formatted.txt \
cat count_exons_lncrnadb_lnc_mm10.txt | awk -F" " '\{print $1"\\t"$2\}' > count_exons_lncrnadb_lnc_mm10_formatted.txt \
cat count_exons_noncode_lnc_mm10.txt | awk -F" " '\{print $1"\\t"$2\}' > count_exons_noncode_lnc_mm10_formatted.txt }