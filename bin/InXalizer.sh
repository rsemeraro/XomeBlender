#!/bin/bash

######################################################################
#             														 #
#              ##      ##											 #
#               ##    ##         								     #
# ## ###    ##   ##  ##        #      ##      ##  ###### ##### ####  #
# ## ####   ##    ####        ###     ##      ##     ##  ##    ##  # #
# ## ## ##  ##    ####       ## ##    ##      ##    ##   ####  ## #  #
# ## ##  ## ##   ##  ##     #######   ##      ##   ##    ##    ##    #
# ## ##   ####  ##    ##   ##     ##  ##      ##  ##     ##    # #   #
# ## ##    ### ##      ## ##       ## ####### ## ####### ##### #  #  #
#                               									 #
######################################################################

######################
#                    #
#  Define functions  #
#                    #
######################

function display_usage()
{
	COLUMNS=$(tput cols)
	printf "%*s\n" $((($COLUMNS/2)+1)) "~~~~~~~~~~~~~~~~~";printf "%*s\n" $((($COLUMNS/2)+1)) "|               |";printf "%*s\n" $((($COLUMNS/2)+1)) "|   INXALIZER   |";printf "%*s\n" $((($COLUMNS/2)+1)) "|               |";printf "%*s\n" $((($COLUMNS/2)+1)) "~~~~~~~~~~~~~~~~~";printf "%*s\n" $((($COLUMNS/2)+30)) "______________________________________________________________________________";printf "%*s\n" $((($COLUMNS/2)+29)) "InXalizer generates the data useful for Xome-Blender running.       "; printf "%*s\n" $((($COLUMNS/2)+30)) "------------------------------------------------------------------------------"; printf "%*s\n" $((($COLUMNS/2)+30)) "Usage: ./InXalizer.sh -f file -l my_label -r ref.fa -scn 4 -vn 1000 -sa Linear"; printf "%*s\n" $((($COLUMNS/2)+30)) "______________________________________________________________________________"; printf '\e[1;39m%*s\n\e[m' $((($COLUMNS/2)-2)) "Arguments:";printf "%*s" $((($COLUMNS/3)+1)) "-f,--file"; printf "%*s\n" $((($COLUMNS/3)-10)) "Input file (bam format)";printf "%*s" $((($COLUMNS/3)+2)) "-l,--label"; printf "%*s\n" $((($COLUMNS/3)-11)) "A prefix for the files to be generated"; printf "%*s" $((($COLUMNS/3)+14)) "-scn,--subclone_number"; printf "%*s\n" $((($COLUMNS/3)-23)) "The number of subclones (tumoral) to generate"; printf "%*s" $((($COLUMNS/3)+13)) "-vn,--variants_number"; printf "%*s\n" $((($COLUMNS/3)-22)) "The number of somatic variants"; printf "%*s" $((($COLUMNS/3)+20)) "-sa,--subclonal_architecture"; printf "%*s\n" $((($COLUMNS/3)-29)) "The kind of subclonal architecture"; printf '\e[8;39m%*s\e[m' $((($COLUMNS/3)+20)) "-sa,--subclonal_architecture"; printf "%*s\n" $((($COLUMNS/3)-29)) "to be reproduced (Linear or Branched)"; printf "%*s" $((($COLUMNS/3)+6)) "-r,--reference"; printf "%*s\n" $((($COLUMNS/3)-15)) "The reference file (fasta format)"; printf "%*s" $((($COLUMNS/3))) "-b,--bed"; printf "%*s\n" $((($COLUMNS/3)-9)) "The target file in .bed format (not required)";printf "%*s" $((($COLUMNS/3))) "-c,--cnv"; printf "%*s\n" $((($COLUMNS/3)-9)) "Generates a CNV file containing random generated coordinates"; printf '\e[8;39m%*s\e[m' $((($COLUMNS/3))) "-c,--cnv"; printf "%*s\n" $((($COLUMNS/3)-9)) "for duplication or deletion events. [#Events,EventSize]"; printf "%*s" $((($COLUMNS/3)+4)) "-t,--threads"; printf "%*s\n" $((($COLUMNS/3)-13)) "The number of threads used for this analysis";printf "%*s" $((($COLUMNS/3)+3)) "-o,--output"; printf "%*s\n" $((($COLUMNS/3)-12)) "Output directory. If omitted, generates a";printf '\e[8;39m%*s\e[m' $((($COLUMNS/3)+3)) "-o,--output"; printf "%*s\n" $((($COLUMNS/3)-12)) "results directory in the current position"; printf "%*s\n"; printf "%*s\n" $((($COLUMNS/2)+30)) "_______________________________~~~~~~~~~~~~~~~~_______________________________"; printf "%*s\n" $((($COLUMNS/2)-1)) "|   CNV LIST   |"; printf "%*s\n" $((($COLUMNS/2)-1)) "~~~~~~~~~~~~~~~~"; printf "%*s\n" $((($COLUMNS/2)+30)) "______________________________________________________________________________"; printf "%*s\n" $((($COLUMNS/2)+15)) "Usage: ./InXalizer.sh -r ref.fa -cl cnv_list.txt"; printf "%*s\n" $((($COLUMNS/2)+30)) "------------------------------------------------------------------------------"; printf '\e[1;39m%*s\n\e[m' $((($COLUMNS/2)-2)) "Arguments:"; printf "%*s" $((($COLUMNS/3)+6)) "-cl,--cnv_list"; printf "%*s\n" $((($COLUMNS/3)-15)) "A file containing the information"; printf '\e[8;39m%*s\e[m' $((($COLUMNS/3)+6)) "-cl,--cnv_list"; printf "%*s\n" $((($COLUMNS/3)-15)) "to generate multiple cnv file"; printf "%*s" $((($COLUMNS/3)+6)) "-r,--reference"; printf "%*s\n" $((($COLUMNS/3)-15)) "The reference file (fasta format)"; printf "%*s" $((($COLUMNS/3)+4)) "-t,--threads"; printf "%*s\n" $((($COLUMNS/3)-13)) "The number of threads used for this analysis";printf "%*s" $((($COLUMNS/3)+3)) "-o,--output"; printf "%*s\n" $((($COLUMNS/3)-12)) "Output directory. If omitted, generates a";printf '\e[8;39m%*s\e[m' $((($COLUMNS/3)+3)) "-o,--output"; printf "%*s\n" $((($COLUMNS/3)-12)) "results directory in the current position"; printf "%*s\n" $((($COLUMNS/2)+30)) "______________________________________________________________________________"; printf "%*s\n" $((($COLUMNS/2)+30)) "Each line of the list file must contain a set of IDs for each subclone (no con"; printf "%*s\n" $((($COLUMNS/2)+30)) "-trol necessary), a label for the cnv file to be generated, the number and the"; printf "%*s\n" $((($COLUMNS/2)+30)) "size of duplication or deletion events (Number of events, Size of each event) "; printf "%*s\n" $((($COLUMNS/2)-6)) "and the path to a bed file (not required)."; printf "%*s\n" $((($COLUMNS/2)+30)) "------------------------------------------------------------------------------"; printf "%*s\n" $((($COLUMNS/2)+30)) "Xome-Blender. Written by Roberto Semeraro, Department of Clinical and Sperime-"; printf "%*s\n" $((($COLUMNS/2)+30)) "ntal Medicine, University of Florence. For bug reports or suggestion write to "; printf "%*s\n" $((($COLUMNS/2)-25)) "robe.semeraro@gmail.com"; printf "%*s\n"
}

function dependencies_check()
{
	type samtools >/dev/null 2>&1 || { echo -e '\e[31mIt seems that samtools is not installed on your computer, please install it.\e[0m'; exit; }
	if [[ $(wc -l < $FilesFolder"/deps.txt") -lt 3 ]] ; then
		if [[ -z $(grep "picard.jar" "$FilesFolder/deps.txt") ]] ; then
			echo -e '\e[31mIt seems that Picard is not installed on your computer.\e[0m' && exit
		elif [[ -z $(grep "GenomeAnalysisTK.jar" "$FilesFolder/deps.txt") ]] ; then
			echo -e '\e[31mIt seems that GATK is not installed on your computer.\e[0m' && exit
		elif [[ -z $(grep "vcftools" "$FilesFolder/deps.txt") ]] ; then
			echo -e '\e[31mIt seems that vcftools is not installed on your computer.\e[0m' && exit
		fi
	fi
	PicardExeLine=$(head -n 1 "$FilesFolder/deps.txt")
	GATKExeLine=$(sed -n '2p' "$FilesFolder/deps.txt")
	PicardExe=$(echo "java -jar $PicardExeLine")
	GATKExe=$(echo "java -jar $GATKExeLine")
	Inx_counter=$(echo "$FilesFolder/Scripts/inx_counter")
	Vars_counter=$(echo "$FilesFolder/Scripts/count_vars")
	Randlines=$(echo "$FilesFolder/Scripts/randlines")
}

gold_gen_wg()
{
	SamActiveProcess=$(ps aux | grep "samtools depth ${SAMBED}" | awk '{print $2}' | wc -l)
	if [[ -z $GATKActiveProcess ]]; then
		GATKActiveProcess=0
	fi
		if [[ ${THREADS} -gt 1 ]] ; then
		REMCPUs=$((${THREADS}-$SamActiveProcess-$GATKActiveProcess))
		GATKCpus=$REMCPUs
	else
		GATKCpus=1
	fi
	paste <(date +"%d/%m/%y %T") <(echo "GATK" "SAMPLE="${LABEL} "GATKCPus="$GATKCpus "REMCPUs="$REMCPUs "COVCPUs="$COVCPUs "- start") -d " " >> $WorkingFolder"/"${LABEL}.log	
	eval $GATKExe -T HaplotypeCaller \
	-R ${REFERENCE} \
	-I ${FILE} \
	-XL "ExcludeList.intervals" \
	-o ${LABEL}".vcf" \
	-rf BadCigar \
	-nct $GATKCpus \
	--max_alternate_alleles 1 \
	--doNotRunPhysicalPhasing \
	-stand_call_conf 1000 \
	-stand_emit_conf 1000 \
	-out_mode EMIT_VARIANTS_ONLY --variant_index_type LINEAR --variant_index_parameter 128000 &>/dev/null
	paste <(date +"%d/%m/%y %T") <(echo "GATK" "SAMPLE="${LABEL} "GATKCPus="$GATKCpus "REMCPUs="$REMCPUs "COVCPUs="$COVCPUs "- finish") -d " " >> $WorkingFolder"/"${LABEL}.log
}

gold_gen_wes()
{
	SamActiveProcess=$(ps aux | grep "samtools depth ${SAMBED}" | awk '{print $2}' | wc -l)
	if [[ -z $GATKActiveProcess ]]; then
		GATKActiveProcess=0
	fi
		if [[ ${THREADS} -gt 1 ]] ; then
		REMCPUs=$((${THREADS}-$SamActiveProcess-$GATKActiveProcess))
		GATKCpus=$REMCPUs
	else
		GATKCpus=1
	fi
	paste <(date +"%d/%m/%y %T") <(echo "GATK" "SAMPLE="${LABEL} "GATKCPus="$GATKCpus "REMCPUs="$REMCPUs "COVCPUs="$COVCPUs "- start") -d " " >> $WorkingFolder"/"${LABEL}.log
	eval $GATKExe -T HaplotypeCaller \
	-R ${REFERENCE} \
	-I ${FILE} \
	-L ${BED} \
	-XL "ExcludeList.intervals" \
	-o ${LABEL}".vcf" \
	-rf BadCigar \
	-nct $GATKCpus \
	--max_alternate_alleles 1 \
	--doNotRunPhysicalPhasing \
	-stand_call_conf 1000 \
	-stand_emit_conf 1000 \
	-out_mode EMIT_VARIANTS_ONLY --variant_index_type LINEAR --variant_index_parameter 128000 &>/dev/null
	paste <(date +"%d/%m/%y %T") <(echo "GATK" "SAMPLE="${LABEL} "GATKCPus="$GATKCpus "REMCPUs="$REMCPUs "COVCPUs="$COVCPUs "- finish") -d " " >> $WorkingFolder"/"${LABEL}.log		
}

val_extract()
{
	arr=($(samtools depth ${SAMBED} -r $c ${FILE} | awk '{print $3}'))
	IFS=$'\n' sortedarr=($(sort -n <<<"${arr[*]}"))
	MedianElement=$((${#sortedarr[@]}/2))
	Median=$(echo ${sortedarr[$MedianElement]})
	echo $Median >> $CovRun\.Medians
	paste <(date +"%d/%m/%y %T") <(echo "SAMTOOLS" "SAMPLE="${LABEL} "CHR="$c "FLAG="$Flag "COVCPUs="$COVCPUs "- finish") -d " " >> $WorkingFolder"/"${LABEL}.log
}

Cov_func()
{
	Flag=0 # Tiene conto delle cpu usate nei processi di calcolo coverage
	progress_cov &
	BarPid=$!
	for c in $Chromosomes
	do
	Flag=$(($Flag+1))
	val_extract &
	paste <(date +"%d/%m/%y %T") <(echo "SAMTOOLS" "SAMPLE="${LABEL} "CHR="$c "FLAG="$Flag "COVCPUs="$COVCPUs "- start") -d " " >> $WorkingFolder"/"${LABEL}.log
	if [[ $Flag -eq $COVCPUs ]]; then
		SamPid=`ps -ef | grep -w "samtools depth ${SAMBED} -r $c ${FILE}" | sed '/grep/d' | awk '{print $3}'`
		while [[ -z $SamPid ]]; do sleep 0.1; SamPid=`ps -ef | grep -w "samtools depth ${SAMBED} -r $c ${FILE}" | sed '/grep/d' | awk '{print $3}'`; done
		SamWaitStrng=$( echo "while [ -e /proc/$SamPid ]; do sleep 0.1; done" )
		eval "$SamWaitStrng"
		Flag=0
	fi
	done
	while [[ ! -e $CovRun\.Medians  || `wc -l $CovRun\.Medians | awk '{print $1}'` -lt 6 ]]; do sleep 0.1; done
        awk '{ total += $1; count++ } END { printf "%.0f\n", total/count }' $CovRun\.Medians > ${LABEL}.cov
        rm $CovRun\.Medians
	echo -ne "\r\033[K      [====================>] - Finished." && kill $BarPid 
	wait $BarPid 2>/dev/null

}

cov_calc()
{
	if [[ ! -e "${LABEL}.cov" ]] ; then
		echo -e " 1) - "'\e[37mCalculation of coverage\e[0m'
		CovRun=.C$RANDOM
		Cov_func
	fi
}

prepare()
{
	if [[ ! -e ${LABEL}".vcf" ]] ; then
		echo -e "\n 2) - "'\e[37mCreation of variant file\e[0m'
		progress &
		BarPid2=$!
		$GS # Generate gVCF
		
		grep "0/1" ${LABEL}".vcf" > ${LABEL}".eth.vcf"
		grep "#" ${LABEL}".vcf" > ${LABEL}".head"
		rm ${LABEL}".vcf"
		cat ${LABEL}".head" ${LABEL}".eth.vcf" > ${LABEL}".vcf"
		rm ${LABEL}".head"
		rm ${LABEL}".eth.vcf"
	
		FileVCFIn=$WorkingFolder"/"${LABEL}".vcf"
		R --slave --args $FileVCFIn,${LABEL},${SA},${VN},${SCN},$WorkingFolder < $FilesFolder/Scripts/SubCloneGenerator.R 2>/dev/null
		rm ExcludeList.intervals
		echo -ne "\r\033[K      [======================>] - Finished." && kill $BarPid2 
		wait $BarPid2 2>/dev/null
	fi
}

intervals()
{
	
	PathForR=$WorkingFolder

	INT_IDs=${IDs//,/-}
	CHECK_CHR=`echo $PidChr | sed 's/-r//' | tr -d " \t\n\r"`

	if [[ ! -z ${BED} ]] ; then
		R --slave --args $PathForR,${LABEL},$DelSize,$Events,$CHECK_CHR,$FilesFolder,$INT_IDs,${REFERENCE},${BED} < $FilesFolder/Scripts/Intervals_generator.R  2>/dev/null
	else
		R --slave --args $PathForR,${LABEL},$DelSize,$Events,$CHECK_CHR,$FilesFolder,$INT_IDs,${REFERENCE} < $FilesFolder/Scripts/Intervals_generator.R  2>/dev/null
	fi
}

genovars()
{
	cd $SubDir
	paste <(date +"%d/%m/%y %T") <(echo "PURIFICATION" "SAMPLE="${ID[$f]} "HCCPUS="${HCCPUS} "- start") -d " " >> $WorkingFolder"/"${LABEL}.log
	if [[ -e ${ID[$f]}".remove" ]] ; then
		if [[ ! -e ${ID[$f]}".addreads.1.bam" ]] ; then
			VarFile=`echo ${ID[$f]}".remove"`
			VarSize=`wc -l $VarFile | awk '{print $1}'`
			while read CMD; do
        			case "$CMD" in \#*) continue ;; esac
        			VariantGeneration $CMD
			done < "$VarFile"
			FilesToRemove=0
			IFS=' ' read -ra CountToRemove <<< "$FilesToRemove"
			while [ ${#CountToRemove[@]} -lt $VarSize ]; do FilesToRemove=`ls -a  | grep ${ID[$f]} | grep refreads | tr '\n' ' '` ; IFS=' ' read -ra CountToRemove <<< "$FilesToRemove"; sleep 0.1; done
			cat "."${ID[$f]}".varreads" | cut -f 4 > ${ID[$f]}".varreadsname"
			$PicardExe FilterSamReads I= $FILE FILTER=excludeReadList QUIET=TRUE VALIDATION_STRINGENCY=LENIENT READ_LIST_FILE= ${ID[$f]}".varreadsname" WRITE_READS_FILES=false O= ${ID[$f]}".1.bam" 2>/dev/null &
			RefToRemove=$FilesToRemove
			FilesToRemove=0
			IFS=' ' read -ra CountToRemove <<< "$FilesToRemove"
			$Vars_counter "."${ID[$f]}".varreads" | while read CMD; do arrIN=(${CMD// / }); $Randlines .${ID[$f]}.${arrIN[0]}.refreads ${arrIN[1]} > ${ID[$f]}.${arrIN[0]}.addreads; done
			while [ ${#CountToRemove[@]} -lt $VarSize ]; do FilesToRemove=`ls -1 *addreads | tr '\n' ' '` ; IFS=' ' read -ra CountToRemove <<< "$FilesToRemove"; sleep 0.1; done
			rm $RefToRemove ".temp.read" 2> /dev/null
			cat $FilesToRemove > ${ID[$f]}".addreads" 2> /dev/null
			cat ${ID[$f]}".addreads" 2> /dev/null | cut -f 4- > ${ID[$f]}".addreads.sam" 2> /dev/null
			wait
			rm $FilesToRemove ${ID[$f]}".addreads" 2> /dev/null
			samtools view -bT $REFERENCE ${ID[$f]}".addreads.sam" > ${ID[$f]}".addreads.bam"
			samtools view -H $FILE > "header"
			samtools reheader "header" ${ID[$f]}".addreads.bam" > ${ID[$f]}".addreads.0.bam"
			samtools sort  -T "/tmp/"${ID[$f]}".aln.sorted" -o ${ID[$f]}".addreads.1.bam" ${ID[$f]}".addreads.0.bam"
			wait
			rm "header" 2> /dev/null
			rm ${ID[$f]}".varreadsname" 2> /dev/null
			rm ${ID[$f]}".addreads.sam" 2> /dev/null
			rm ${ID[$f]}".addreads.bam" 2> /dev/null
			rm "."${ID[$f]}".varreads" 2> /dev/null
			rm ${ID[$f]}".addreads.0.bam" 2> /dev/null
		fi	
		while [ $(( $(date +%s) - $(stat -c %Y ${ID[$f]}".1.bam") )) -lt 20 ]; do sleep 1; done
		samtools merge ${ID[$f]}".bam" ${ID[$f]}".1.bam" ${ID[$f]}".addreads.1.bam"
		rm ${ID[$f]}".1.bam" 2> /dev/null
		rm ${ID[$f]}".addreads.1.bam" 2> /dev/null
		rm $VarFile 2> /dev/null
		mv ${ID[$f]}".bam" $WorkingFolder 2> /dev/null
		if [[ -e ${ID[$f]}".variants" ]] ; then
			if [[ ${ID[$f]} != "${LABEL}_Control" ]] ; then
				mv ${ID[$f]}".variants" $WorkingFolder 2> /dev/null
				mv ${ID[$f]}".vcf" $WorkingFolder 2> /dev/null
			fi
		fi
	elif [[ ! -e ${ID[$f]}".remove" ]]; then
		cp $FILE . 2> /dev/null
		mv *".bam" ${ID[$f]}".bam" 2> /dev/null
		mv ${ID[$f]}".variants" $WorkingFolder 2> /dev/null
		mv ${ID[$f]}".vcf" $WorkingFolder 2> /dev/null
		wait
		mv ${ID[$f]}".bam" $WorkingFolder 2> /dev/null
	fi
	paste <(date +"%d/%m/%y %T") <(echo "PURIFICATION" "SAMPLE="${ID[$f]} "HCCPUS="${HCCPUS} "- finish") -d " " >> $WorkingFolder"/"${LABEL}.log
	cd $WorkingFolder
	rmdir $SubDir 2> /dev/null
}

VariantGeneration(){
        VcfChr=$(echo $1)
        VcfPos=$(echo $2)
        VcfRef=$(echo $4)
        VcfAlt=$(echo $5)

        ### VarType Check ###

        if [[ ${#VcfRef} -gt ${#VcfAlt} ]] ; then
                VarType="D"
                Elgth=$((${#VcfRef} - 1))
        elif [[ ${#VcfRef} -lt ${#VcfAlt} ]] ; then
                VarType="I"
                Elgth=0
        elif [[ ${#VcfRef} -eq ${#VcfAlt} ]] ; then
                VarType="S"
                Elgth=0
        fi

        ### Length Calculation ###

        TrueLength=$(($Elgth + $VcfPos))

        ### Read discovery ###

        samtools view -F 0x04 $FILE $VcfChr:$VcfPos-$TrueLength |
                while IFS= read -r line
                do
			echo -e "${ID[$f]}\t$VcfPos\t$VarType\t$line" > ".temp.read"
                        $Inx_counter ".temp.read" $VcfRef 2>/dev/null |:
                done
}

startgraph()
{
	PriID=`echo " - SAMPLES = "${IDs//,/ }`
	Samplength=`echo ${#PriID}`
	LineLength=$(($Samplength+1))
	HalfLength=$(($LineLength/2))
	BorderLines=`for ((x = 0; x < $LineLength ; x++)); do   printf %s =; done`
	printf '\e[1;39m%*s\n\e[m' $(($HalfLength+26)) " _           _    _          _                     "
	printf '\e[1;39m%*s\n\e[m' $(($HalfLength+26)) "| |  _      \ \  / /        | | _   ___   ___  _ _ "
	printf '\e[1;39m%*s\n\e[m' $(($HalfLength+26)) "| | | |___   \ \/ /   __ _  | |(_) \  _| / _ \| '_/"
	printf '\e[1;39m%*s\n\e[m' $(($HalfLength+26)) "| | |  _  \  / /\ \  / _' | | || | _\ \ |  __/| |  "
	printf '\e[1;39m%*s\n\e[m' $(($HalfLength+26)) "|_| |_| |_| /_/  \_\ \__,_| |_||_||____| \___||_|  "
	echo $BorderLines
	echo -e " - "'\e[31mSAMPLES\e[0m' = ${IDs//,/ }
	echo $BorderLines
	printf '\e[1;39m%*s\n\e[m' $((($LineLength/2)+5)) "InXalizing:"
	printf "%*s\n"
}

function progress () {
    s=0.5;
    f=0.25;
    echo -ne "\r";
    while true; do
           sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [                      ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [=>                    ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [==>                   ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [===>                  ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [====>                 ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [=====>                ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [======>               ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [=======>              ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [========>             ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [=========>            ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [==========>           ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [===========>          ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [============>         ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [=============>        ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [==============>       ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [===============>      ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [================>     ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [=================>    ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [==================>   ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [===================>  ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [====================> ] Elapsed: ${s} secs.";
           sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [=====================>] Elapsed: ${s} secs.";
    done;
}

function progress_cov () {
    s=0.5;
    f=0.25;
    echo -ne "\r";
    while true; do
           sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [                     ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [=>                   ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [==>                  ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [===>                 ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [====>                ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [=====>               ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [======>              ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [=======>             ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [========>            ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [=========>           ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [==========>          ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [===========>         ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [============>        ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [=============>       ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [==============>      ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [===============>     ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [================>    ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [=================>   ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [==================>  ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [===================> ] Elapsed: ${s} secs.";
           sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [====================>] Elapsed: ${s} secs.";
    done;	
}

function progress_somatic () {
    s=0.5;
    f=0.25;
    echo -ne "\r";
    while true; do
           sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [               ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [=>             ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [==>            ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [===>           ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [====>          ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [=====>         ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [======>        ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [=======>       ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [========>      ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [=========>     ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [==========>    ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [===========>   ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [============>  ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [=============> ] Elapsed: ${s} secs.";
           sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [==============>] Elapsed: ${s} secs.";
    done;	
}

function progress () {
    s=0.5;
    f=0.25;
    echo -ne "\r";
    while true; do
           sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [                      ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [=>                    ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [==>                   ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [===>                  ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [====>                 ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [=====>                ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [======>               ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [=======>              ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [========>             ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [=========>            ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [==========>           ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [===========>          ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [============>         ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [=============>        ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [==============>       ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [===============>      ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [================>     ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [=================>    ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [==================>   ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [===================>  ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [====================> ] Elapsed: ${s} secs.";
           sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [=====================>] Elapsed: ${s} secs.";
    done;
}

function progress_cnv () {
    s=0.5;
    f=0.25;
    echo -ne "\r";
    while true; do
           sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [                         ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [=>                       ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [==>                      ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [===>                     ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [====>                    ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [=====>                   ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [======>                  ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [=======>                 ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [========>                ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [=========>               ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [==========>              ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [===========>             ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [============>            ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [=============>           ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [==============>          ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [===============>         ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [================>        ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [=================>       ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [==================>      ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [===================>     ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [====================>    ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [=====================>   ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [======================>  ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [=======================> ] Elapsed: ${s} secs.";
           sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [========================>] Elapsed: ${s} secs.";
    done;
}

######################
### Prepare Module ###
######################

if [  $# -le 1 ]
	then
	display_usage
	exit 1
fi

if [[ ( $1 == "--help") ||  $1 == "-h" ]]
	then
	display_usage
	exit 0
fi

while [[ $# > 1 ]]
do
	key="$1"
	shift
	case $key in
	-f|--file)
	FILE="$1"
	shift
	;;
	-l|--label)
	LABEL="$1"
	shift
	;;
	-scn|--subclone_number)
	SCN="$1"
	shift
	;;
	-vn|--variants_number)
	VN="$1"
	shift
	;;
	-sa|--subclonal_architecture)
	SA="$1"
	shift
	;;
	-r|--reference)
	REFERENCE="$1"
	shift
	;;
	-c|--cnv)
	CNV="$1"
	shift
	;;
	-cl|--cnvlist)
	CNVLIST="$1"
	shift
	;;
	-b|--bed)
	BED="$1"
	shift
	;;
	-t|--threads)
	THREADS="$1"
	shift
	;;
	-o|--output)
	OUTPUT="$1"
	shift
	;;
	esac
done

######################
# Enviroment setting #
######################

if [[ -z ${REFERENCE} ]] ; then
	echo -e '\e[31mNo reference selected!\e[0m' && exit
else
	if [[ ! -e "${REFERENCE}.fai" || ! -e ""${REFERENCE%.*}".dict" ]] ; then
		echo -e '\e[31mThe fasta index or the dictionary are missing!\e[0m' && exit
	fi
fi

if [[ -z ${THREADS} ]] ; then
	THREADS=1
elif [[ ${THREADS} =~ [[:alpha:]] || ${THREADS} -lt 0 ]] ; then
	echo "Number of threads it's not an integer!" && exit
fi

SOURCE="${BASH_SOURCE[0]}"
TARGET="$(readlink "$SOURCE")"
if [[ -z $TARGET ]]; then
	TARGET=$SOURCE
fi

FilesFolder=`echo $(dirname "${TARGET%/*}")`

WorkingFolder=`pwd`

######################
# Dependencies check #
######################

dependencies_check

#####################
# Chr Style Control #
#####################

ChrStyle=`grep ">" ${REFERENCE} | awk '{print $1'} | head -n 1 | sed 's/^.//'`
if [[ "$ChrStyle" == *"chr"* ]]; then
	Chromosomes="chr1 chr2 chr3 chr4 chr5 chr6"
	PidChr="-r chr1"
else
	Chromosomes="1 2 3 4 5 6"
	PidChr="-r 1"
fi

### CNV list ###

if [[ ! -z ${CNVLIST} ]]; then
	COLUMNS=$(tput cols)
	printf '\e[37m%*s\n\e[m' $((($COLUMNS/3)+2)) "#####################"; printf '\e[37m%*s\n\e[m' $((($COLUMNS/3)+2)) " # CNV file creation #"; printf '\e[37m%*s\n\e[m' $((($COLUMNS/3)+2)) " #####################"

	CNVCounter=0

	array=()
	getArray() {
	i=0
	while read line
		do
			array[i]=$line
			i=$(($i + 1))
		done < $1
	}
	getArray ${CNVLIST}
	for e in "${array[@]}"
	do
		#IDs=`echo $e | awk '{print $1}'`
		LABEL=`echo $e | awk '{print $1}'`
		CNV=`echo $e | awk '{print $2}'`
                ### Modifico qui il 15/09/16 ###
		BED=`echo $e | awk '{print $3}'`
		### Label Check ###
		if [[ ${LABEL} == *"-"* ]]; then
	  		echo -e '\e[31mHyphens are not allowed for labels. Please, replace it.\e[0m' && exit
		fi
		
		### IDs Gen ###
		ControlLabel=`echo ${LABEL}"_Control"`
		SubClLabels=""
		for i in $(seq 1 $SCN)
		do
			SubClLabels=$SubClLabels,${LABEL}"_Subclone_"$i
		done
		IDs=$ControlLabel$SubClLabels
		IFS=',' read -ra ID <<< "$IDs"

		Fields=`echo ${CNV} | tr -cd , | wc -c`
		if [[ $Fields -ne 1 ]] ; then
			echo -e '\e[31mBad formatting, fix it. It must be 44,1000000 (EventsNumber,EventSize).\e[0m' && exit
		else
			IFS=',' read -ra Features <<< "$CNV"
			Events=`echo ${Features[0]}`
			if [[ $Events -lt 2 ]]; then
				echo -e '\e[31mThe number of CNV events must be at least 2.\e[0m' && exit
			fi
			DelSize=`echo ${Features[1]}`
			if [[ $(($DelSize*$Events)) -gt 500000000 ]] ; then
				echo -e '\e[31mThe size of CNV events is too big.\e[0m' && exit
			fi
		fi
		CNVCounter=$(($CNVCounter+1))
		echo -e " - "'\e[31mSAMPLES\e[0m' = ${IDs//,/ }
		echo -e "   - "'\e[30mNumber of events =\e[0m' $Events  "   - "'\e[30mSize of event =\e[0m' $DelSize
		#IDs=`echo CONTROL,$IDs`
		intervals & # Function
		if [[ "$CNVCounter" -eq ${THREADS} ]]; then
			wait
			CNVCounter=0
		fi
	done
	wait
	if [[ ! -z ${OUTPUT} ]] ; then
		if [[ ! -d ${OUTPUT} ]]; then
			mkdir ${OUTPUT}
		fi
		mv *"_CNV.txt" ${OUTPUT}
	fi
	echo -e "     "'\e[37mDone\e[0m'
	exit
fi

### Input control ###

if [[ -z ${FILE} ]] ; then
	echo -e '\e[31mNo input file selected!\e[0m' && exit
else
	if [[ ! -e "${FILE}.bai" ]] ; then
		echo -e '\e[31mThe bam index file is missing!\e[0m' && exit
	fi
fi

FileSize=`stat --printf="%s" $FILE`

if [[ -z ${SCN} ]] ; then
	echo -e '\e[31mYou have to set the number of clones to be generated.\e[0m' && exit
fi

if [[ ${SCN} =~ [[:alpha:]] || ${SCN} -lt 1 || ${SCN} -gt 5 ]]; then
	echo -e '\e[31mNumber of clones it is not an integer or it is not comprised between 1 and 5.\e[0m' && exit
fi

if [[ -z ${VN} ]] ; then
	echo -e '\e[31mPlease, set the variants number!\e[0m' && exit
fi

if [[ ${VN} =~ [[:alpha:]] || ${VN} -gt 5000 ]]; then
	echo -e '\e[31mNumber of variants it is not an integer or it is greater than 5000.\e[0m' && exit
fi

if [[ ${VN} -eq 1 ]]; then
	VN=2
fi

if [[ -z ${SA} ]] ; then
	echo -e '\e[31mNo subclonal architecture selected!\e[0m' && exit
fi

if [[ ${SA} != "Linear" && ${SA} != "Branched" ]] ; then
	echo -e '\e[31mPlease, select Linear or Branched architecture.\e[0m' && exit
fi

if [[ ! -z ${BED} ]] ; then
	GS="gold_gen_wes"
	SAMBED="-b ${BED}"
else
	echo -e '\e[31mNo bed file selected! You using genomic data!\e[0m'
	GS="gold_gen_wg"
	SAMBED=""
fi

if [[ -z ${LABEL} ]] ; then
	echo -e '\e[31mNo label defined!\e[0m' && exit
fi

### Label Check ###

if [[ ${LABEL} == *"-"* ]]; then
	echo -e '\e[31mHyphens are not allowed for labels. Please, replace it.\e[0m' && exit
fi

if [[ ! -z ${CNV} ]] ; then
	Fields=`echo ${CNV} | tr -cd , | wc -c`
	if [[ $Fields -ne 1 ]] ; then
		echo -e '\e[31mBad formatting, fix it. It must be 44,1000000 (EventsNumber,EventSize).\e[0m' && exit
	else
		IFS=',' read -ra Features <<< "$CNV"
		Events=`echo ${Features[0]}`
		if [[ $Events -lt 2 ]]; then
			echo -e '\e[31mThe number of CNV events must be at least 2.\e[0m' && exit
		fi
		DelSize=`echo ${Features[1]}`
		if [[ $(($DelSize*$Events)) -gt 500000000 ]] ; then
			echo -e '\e[31mThe size of CNV events is too big.\e[0m' && exit
		fi
	fi
fi

### Cpu distribution ###

HCCPUS=1
if [[ ${THREADS} -gt 2 ]] ; then
 	CPUGROUPS=$(echo $(perl -E "say (int(${THREADS}/2) )")) # Determina il numero di cpu a gruppi di 2 (dato che 2 cpu servono per la funzione di rimozione varianti)
 	REMAININGCPUs=$(echo $(perl -E "say ${THREADS}-(2*$CPUGROUPS)")) # Se il numero di cpu non è multiplo di 2 calcola quante cpu restano
 	if [[ $REMAININGCPUs -ne 0 ]]; then
 		HCCPUS=$(echo $(perl -E "say $CPUGROUPS+1")) # Se ci sono cpu che avanzano, aumenta il numero di cpu utilizzabili ma modifica il numero di quelle utilizzate per il calcolo del coverage
 	else
 		HCCPUS=$(echo $CPUGROUPS)
 	fi
fi

COVCPUs=6

### IDs Gen ###

ControlLabel=`echo ${LABEL}"_Control"`
SubClLabels=""
for i in $(seq 1 $SCN)
do
	SubClLabels=$SubClLabels,${LABEL}"_Subclone_"$i
done
IDs=$ControlLabel$SubClLabels
IFS=',' read -ra ID <<< "$IDs"

#########
# Start #
#########

startgraph

End=$(echo $(perl -E "say ${#ID[@]}-1"))

if [[ $THREADS -lt $COVCPUs ]]; then
	COVCPUs=$(echo $THREADS)
fi

if [[ ! -e "ExcludeList.intervals" ]]; then
	samtools view -H $FILE | grep SN | awk '{print $2}' | sed -r 's/^.{3}//' | awk -v chchek="chr" '{if ($1 ~ /^ *chr/ ) print substr($1,4), chchek; else print $1}' | awk '{ if ($1 > 22 && $1 !~ /NC_007605/ && $1 !~ /hs37d5/) print $0}' | awk '{ if ($2 == "chr") print $2$1; else print $1}' > ExcludeList.intervals
fi

cov_calc

prepare

CloneCount=0
BarCount=0

for f in $(seq 0 $End)
do
	SubDir=`echo ${ID[$f]}`
	if [[ (! -d $SubDir) && (! -e $SubDir".bam") ]] ; then
		BarCount=$(($BarCount+1))
	fi
done

if [[ $BarCount -ge 1 ]]; then
	echo -e "\n 3) - "'\e[37mSample generation\e[0m'
	progress_somatic &
	BarPid3=$!
else
	echo -e "\n 3) - "'\e[37mSample generation\e[0m'
fi

cov_calc

ConditionCounter=0
for f in $(seq 0 $End)
do
	SubDir=`echo ${ID[$f]}`
	if [[ ! -d $SubDir ]]; then
		mkdir $SubDir
		Bait=`ls -1p | grep -v / | grep $SubDir | tr '\n' ' '`
		cp $Bait $SubDir
	fi
	CloneCount=$(($CloneCount+1))
	genovars &
	if [[ -e $SubDir"/"$SubDir".variants" ]] ; then
		if [[ -e $SubDir".variants" ]] ; then
			rm $Bait
		fi
	elif [[ -e $SubDir"/"$SubDir".remove" ]] ; then
		if [[ -e $SubDir".remove" ]] ; then
			rm $Bait
		fi
	fi
	if [[ $CloneCount -eq ${HCCPUS} ]]; then
		if [[ $SubDir != ${ID[1]} ]] ; then
			Pido=`ps -ef | grep "READ_LIST_FILE= $SubDir.varreadsname WRITE_READS_FILES=false O= $SubDir.1.bam" | sed '/grep/d' | awk '{print $3}'`
			while [[ -z $Pido ]]; do sleep 0.1; Pido=`ps -ef | grep "READ_LIST_FILE= $SubDir.varreadsname WRITE_READS_FILES=false O= $SubDir.1.bam" | sed '/grep/d' | awk '{print $3}'`; done
			WaitStrng=$( echo "while [ -e /proc/$Pido ]; do sleep 0.1; done" )
			eval "$WaitStrng"
			CloneCount=0
		fi
	fi
done

while [[ `ls | grep -c .bam$` -lt ${#ID[@]} ]]; do
	sleep 0.1
done
if [[ $BarCount -ge 1 ]]; then
	echo -ne "\r\033[K      [==============>] - Finished." && kill $BarPid3
	wait $BarPid3 2>/dev/null
fi

if [[ ! -z ${CNV} ]] ; then
	echo -e "\n 4) - "'\e[37mCreation of CNV events file\e[0m'
	progress_cnv &
	BarPid4=$!
	intervals
	echo -ne "\r\033[K      [========================>] - Finished." && kill $BarPid4
	wait $BarPid4 2>/dev/null
	echo -e "\n      "'\e[37mDone\e[0m'
fi

if [[ ! -z ${OUTPUT} ]] ; then
	if [[ ! -d ${OUTPUT} ]]; then
		mkdir ${OUTPUT}
	fi
	mv ${LABEL}".cov" ${OUTPUT}
	mv ${LABEL}".vcf" ${OUTPUT}
	for i in $(seq 1 $End)
	do
		mv ${ID[$i]}".vcf" ${OUTPUT}
	done
	mv *.bam ${OUTPUT}
	mv *.variants ${OUTPUT}
	if [[ ! -z ${CNV} ]] ; then
		mv ${LABEL}"_CNV.txt" ${OUTPUT}
	fi
fi
