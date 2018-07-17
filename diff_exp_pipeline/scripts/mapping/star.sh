 #!/bin/sh

genome="$1"
file1="$2"
file2="$3"
anno="$4"
path="$5"
nrThreads="$6"

/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/Noncoding/softwares/STAR-2.6.0c/bin/Linux_x86_64/STAR --runThreadN "$nrThreads" \
	--genomeDir "$genome" \
	--readFilesIn <(zcat "$file1") <(zcat "$file2") \
	--outBAMsortingThreadN 6 \
	--outSAMtype BAM SortedByCoordinate \
	--outReadsUnmapped Fastx \
	--sjdbGTFfile "$anno" \
	--outFileNamePrefix  "${path}"

