#!/bin/bash

file=$1

### format anno file for NOISeq diff exp
cat "$file" | grep -P "gene\t" | awk '$1 ~ /^chr/ {split($1,a,"chr"); split($10,b,";"); if($7 == "+" ) {print a[2]"\t"$4"\t"$5"\t"b[1];} else print a[2]"\t"$5"\t"$4"\t"b[1];}' | awk -F"\"" '{print $1"\t"$2"\t"$3"\t"$5}' > "$file"_noiseq

### format anno file for cuffdiff
grep -v -P "gene\t" "$file" > "$file"_nogenes
