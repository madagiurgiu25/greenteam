#!/bin/bash

declare -A config # init array
config=( # set default values in config array
    [name]=""
    [nameMapping]=""
    [nameQuant]=""
    [samples]=""
    [size]=0
    [genelist]=""
    [path]=""
    [pathFastq]=""
    [genomeDir]=""
    [anno]=""
    [jsonAnno]=""
    [command]=""
    [nrThreads]=10
    [stringTie]="/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/Noncoding/softwares/stringtie-1.3.4d.Linux_x86_64/stringtie"
    [rawCounts]="/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/daten/expression/pipeline/scripts/quantification/rawCounts.py"
    [formatAnno]="/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/daten/expression/pipeline/scripts/diff_exp/format_anno.sh"
    [runCuffDiff]="/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/daten/expression/pipeline/scripts/diff_exp/runCuffDiff.sh"
    [runBallgown]="/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/daten/expression/pipeline/scripts/diff_exp/runBallgown.R"
    [runNOISeq]="/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/daten/expression/pipeline/scripts/diff_exp/runNOISeq.sh"
    [starScript]="/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/daten/expression/pipeline/scripts/mapping/star.sh"
    [samtools]="/home/proj/biosoft/software/samtools-1.2/samtools"
    [bedCov]="/home/proj/biosoft/software/bedtools-2.17.0/bin/genomeCoverageBed"
    [python]="/home/g/giurgiu/python/Python-3.6.0/python"
    [r]="/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/Noncoding/softwares/R-3.5.0/bin/Rscript"
)

conf_file=$1

while read line
do
    if echo $line | grep -F = &>/dev/null
    then
        varname=$(echo "$line" | cut -d '=' -f 1)
        config[$varname]=$(echo "$line" | cut -d '=' -f 2-)
    fi
done < $conf_file

star_root="${config[path]}STAR/${config[nameMapping]}"
stringtie_root="${config[path]}StringTie/${config[nameQuant]}/"
ballgown_input_root="${config[path]}Ballgown/${config[nameQuant]}"
ballgown_root="${config[path]}DiffExp/${config[name]}/Ballgown"
cuffdiff_root="${config[path]}DiffExp/${config[name]}/Cuffdiff"
noiseq_root="${config[path]}DiffExp/${config[name]}/NOISeq"

> output.log
echo "#########################################################################"  | tee -a output.log
echo "############################## START ####################################"  | tee -a output.log
echo "#########################################################################"  | tee -a output.log
echo "# "
echo "# GENERAL SETTINGS:" | tee -a output.log
echo "# Name: ${config[name]}" | tee -a output.log
echo "# Samples: "${config[samples]}  | tee -a output.log
cat ${config[samples]} | tee -a output.log
echo "" | tee -a output.log
echo "# Replicates count: ${config[size]}" | tee -a output.log
echo "# Working dir: "${config[path]}  | tee -a output.log
echo "# GenomeDir: "${config[genomeDir]}  | tee -a output.log
echo "# GeneList: "${config[genelist]}  | tee -a output.log
echo "# Command: "${config[command]}  | tee -a output.log
echo "# Log file: output.log" | tee -a output.log

case ${config[command]} in
        mapping)

            echo "#########################################################################"  | tee -a output.log
            echo "################################ start STAR #############################"  | tee -a output.log
            echo "#########################################################################"  | tee -a output.log

            echo "STAR root: $star_root" | tee -a output.log
            mkdir -p $star_root

            while IFS='' read -r line || [[ -n "$line" ]];
            do
	            file1=`echo "$line" | awk -F"\t" '{print $1}'`
	            file2=`echo "$line" | awk -F"\t" '{print $2}'`
	            sample_name=`echo "$line" | awk -F"\t" '{print $3}'`

	            mkdir -p "${star_root}/${sample_name}"
	            echo "map ... $file1 $file2 "  | tee -a output.log
                ${config[starScript]} \
                    ${config[genomeDir]} \
                    ${config[pathFastq]}${file1} \
                    ${config[pathFastq]}${file2} \
                    ${config[anno]} \
                    "${star_root}/${sample_name}/mappedstar_" \
                    ${config[nrThreads]}

                #########################################################################
                ################################ sort samtools ##########################
                #########################################################################
	            echo "sort ..."  | tee -a output.log

	            cd "${star_root}"
	            ${config[samtools]} sort -@ 8 "${star_root}/${sample_name}/mappedstar_Aligned.sortedByCoord.out.bam" "${star_root}/${sample_name}/star_sorted"

	            echo "create index ... "  | tee -a output.log
	            ${config[samtools]} index "${star_root}/${sample_name}/star_sorted.bam"

	            rm "${star_root}/${sample_name}/mappedstar_Aligned.sortedByCoord.out.bam"

	            ############### BAM coverage #######################
	            echo "run coverage ..." | tee -a output.log
	            ${config[bedCov]} -dz -split \
	                -ibam "${star_root}/${sample_name}/star_sorted.bam" > "${star_root}/${sample_name}/coverage.bed"

	            ############### Fragment length ####################
	            echo "fragment len ... " | tee -a output.log
	            ${config[samtools]} view "${star_root}/${sample_name}/star_sorted.bam" | awk '{print $9}' | sort | uniq -c > "${star_root}/${sample_name}/fragmentLen.txt"

            done < ${config[samples]}

            ;;


        quantification)

            echo "#########################################################################"  | tee -a output.log
            echo "################################ start StringTie ########################"  | tee -a output.log
            echo "#########################################################################"  | tee -a output.log

            mkdir -p "$stringtie_root"
            mkdir -p "$ballgown_input_root"
            > "${ballgown_input_root}gene_exp.tab"

            echo "# StringTie root: $stringtie_root"  | tee -a output.log
            echo "# Ballgown quant root: $ballgown_input_root"  | tee -a output.log
            echo "#" | tee -a output.log

            while IFS='' read -r line || [[ -n "$line" ]];
            do

                file1=`echo "$line" | awk -F"\t" '{print $1}'`
	        file2=`echo "$line" | awk -F"\t" '{print $2}'`
	        sample_name=`echo "$line" | awk -F"\t" '{print $3}'`

	        pathStringTie="${stringtie_root}${sample_name}/"
	        pathBallgown="${ballgown_input_root}/${sample_name}/"

	        mkdir -p "$pathStringTie"
	        mkdir -p "$pathBallgown"

                echo "# Sample: $sample_name" | tee -a output.log
                echo "# StringTie path: $pathStringTie" | tee -a output.log
                echo "# Ballgown path: $pathBallgown" | tee -a output.log
	        echo "# stringtie transcript assembly..." | tee -a output.log
		bamfile="${star_root}${sample_name}/star_sorted.bam"

	            ${config[stringTie]} \
	                "${bamfile}" \
	                -e -b "${pathBallgown}" \
	                -G "${config[anno]}" \
	                -p 8 \
	                -A "${pathStringTie}${sample_name}_geneabundance.tab" \
	                -o "${pathStringTie}${sample_name}_star.gtf"

                echo "# compute raw counts ..." | tee -a output.log
                echo "# Output: ${pathStringTie}${sample_name}_geneabundance.tab_rawcounts" | tee -a output.log
                ${config[python]} -u ${config[rawCounts]} "${pathStringTie}${sample_name}_geneabundance.tab" ${config[jsonAnno]} | tee -a output.log

                > "${pathStringTie}${sample_name}_rawcounts.tab"
                > "${pathStringTie}${sample_name}_tmp.tab"
                > "${pathStringTie}${sample_name}_fpkm.tab"
                echo "Gene_ID	RAW.$sample_name" >> "${pathStringTie}${sample_name}_rawcounts.tab"
                tail -n+2 "${pathStringTie}${sample_name}_geneabundance.tab_rawcounts"  | awk -F"\t" '{print $1"\t"$10}' >> "${pathStringTie}${sample_name}_rawcounts.tab"

                echo "Gene_ID	TPM.$sample_name" >> "${pathStringTie}${sample_name}_tmp.tab"
                tail -n+2 "${pathStringTie}${sample_name}_geneabundance.tab_rawcounts"  | awk -F"\t" '{print $1"\t"$9}' >> "${pathStringTie}${sample_name}_tmp.tab"

                echo "Gene_ID	FPKM.$sample_name" >> "${pathStringTie}${sample_name}_fpkm.tab"
                tail -n+2 "${pathStringTie}${sample_name}_geneabundance.tab_rawcounts"  | awk -F"\t" '{print $1"\t"$8}' >> "${pathStringTie}${sample_name}_fpkm.tab"

                echo "# ----------------------------------------------------" | tee -a output.log
                echo "# " | tee -a output.log

            done < ${config[samples]}

            ;;


        diffexp)

            echo "#########################################################################" | tee -a output.log
            echo "################################ start DiffExp ##########################" | tee -a output.log
            echo "#########################################################################" | tee -a output.log

            # format anno file for cuffidff and noiseq
            # output: *.gft_nogenes
            # output: *.gtf_noiseq
            bash ${config[formatAnno]} ${config[anno]}

            mkdir -p "$ballgown_root/"
            mkdir -p "$cuffdiff_root/"
            mkdir -p "$noiseq_root/"

            ############################## Run Ballgown
            echo "Run ballgown ... ${config[name]} ... ${config[size]} ...." | tee -a output.log
            ${config[r]} ${config[runBallgown]} ${config[name]} ${config[size]} $ballgown_input_root $star_root $ballgown_root "TRUE" "TRUE" "TRUE" ${config[samples]} ${config[genelist]}
            #newsize=$((config[size]-1))
            #echo "Run ballgown subsets... ${config[samples]}  ... ${newsize} ...." | tee -a output.log
            #${config[r]} ${config[runBallgown]} ${config[name]} $newsize $ballgown_input_root $star_root $ballgown_root "TRUE" "FALSE" "FALSE" ${config[samples]} ${config[genelist]}


            ############################### Run Cuffidff
            echo "Run Cuffdiff .... ${config[name]}  ..." | tee -a output.log
            condition1=`cat "${ballgown_root}/bams.txt" | awk '{print $1}'`
            condition2=`cat "${ballgown_root}/bams.txt" | awk '{print $2}'`

            bash ${config[runCuffDiff]} $cuffdiff_root $condition1 $condition2 "${config[anno]}_nogenes"


            ;;

        *)
            echo $"Usage: $0 {command=mapping|quantification|diffexp}" | tee -a output.log
            printf '%s\n' "${config[@]}" | tee -a output.log
            exit 1
esac



#if [[ $mapping || $quantification ]]
#then
#
#
#
#    while IFS='' read -r line || [[ -n "$line" ]];
#    do
#	    file1=`echo "$line" | awk -F"\t" '{print $1}'`
#	    file2=`echo "$line" | awk -F"\t" '{print $2}'`
#	    name=`echo "$line" | awk -F"\t" '{print $3}'`
#
#	    #########################################################################
#	    ################################ start STAR #############################
#	    #########################################################################
#	    mkdir -p "${star_root}${name}"
#	        "$starScript" \
#	        "$genomeDir" \
#	        "${path}${file1}" \
#	        "${path}${file2}" \
#	        "$anno" \
#	"${star_root}${name}/mappedstar_" \
#	"$nrThreads"
#
#	#########################################################################
#        ################################ sort samtools ##########################
#        #########################################################################
#	echo "$file1 $file2"
#	echo "sort ..."
#
#	cd "${star_root}${name}/"
#	/home/proj/biosoft/software/samtools-1.2/samtools sort -@ 8 "${star_root}${name}/mappedstar_Aligned.sortedByCoord.out.bam" "${star_root}${name}/star_sorted"
#
#	echo "create_ index"
#	/home/proj/biosoft/software/samtools-1.2/samtools index "${star_root}${name}/star_sorted.bam"
#
#	rm "${star_root}${name}/mappedstar_Aligned.sortedByCoord.out.bam"
#
#	 ############### BAM coverage #######################
#	 echo "run coverage"
#	 /home/proj/biosoft/software/bedtools-2.17.0/bin/genomeCoverageBed -dz -split \
#	 -ibam "${star_root}${name}/star_sorted.bam" > "${star_root}${name}/coverage.bed"
#
#	############### Fragment length ####################
#	echo "fragment len"
#	/home/proj/biosoft/software/samtools-1.2/samtools view "${star_root}${name}/star_sorted.bam" | awk '{print $9}' | sort | uniq -c > "${star_root}${name}/fragmentLen.txt"
#
#
#	#################################################################################################
#        ################################ start StringTie - Transcript assembly and abundace estimations##
#        #################################################################################################
#
#	pathStringTie="${stringtie_root}${name}/"
#	pathBallgown="${ballgown_root}${name}/"
#
#	mkdir -p "$pathStringTie"
#	mkdir -p "$pathBallgown"
#
#	echo "stringtie transcript assembly"
#	bamfile="${star_root}${name}/star_sorted.bam"
#
#	$stringTie \
#	"${bamfile}" \
#	-e -b "${pathBallgown}" \
#	-G "$anno" \
#	-p 8 \
#	-A "${pathStringTie}${name}_geneabundance.tab" \
#	-o "${pathStringTie}${name}_star.gtf"
#
#done < "$samples"
#fi

############################################ merge transcripts ######################################
#####################################################################################################
#
##/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/Noncoding/softwares/stringtie-1.3.4d.Linux_x86_64/stringtie \
##        --merge -p 8 \
##        -G "$ianno" \
##        -o "$path"StringTie/stringtie_merged.gtf "$path"StringTie/mergeList.txt
#
#
#################################################################################################
########################## Diff exp using Cufflinks #############################################
#################################################################################################
#
##pathCuffdiff="${cuffdiff_root}H103/"
#
##mkdir -p $pathCuffdiff
#
##condition1=`ls "${star_root}"H103A*/star_sorted.bam | paste -sd "," -`
##condition2=`ls "${star_root}"H103B*/star_sorted.bam | paste -sd "," -`
#
#
##$cuffdiff -o $pathCuffdiff -p 5 -c 5 -v --FDR 0.05 \
##--library-type fr-firststrand \
##"$anno" \
##"$condition1" "$condition2"
#
#
#
#

