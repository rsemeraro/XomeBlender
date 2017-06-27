#!/bin/bash

############################################################
#                                                          #
#     ##    ##     #####    ###        ###  #########      #
#      ##  ##    ##    ##   ####      ####  ##             #
#       ####    ##      ##  ## ##    ## ##  ######         #
#       ####    ##      ##  ##   ####   ##  ######         #
#      ##  ##    ##    ##   ##          ##  ##             #
#     ##    ##    #####     ##          ##  #########      #
#                                                          #
#  ####    ##     #####  ####    ##  ## #    #####  ####   #
#  ##  #   ##     ##     ## #    ##  ##  #   ##     ##  #  #
#  ## #    ##     ##     ##  #   ##  ##   #  ##     ## #   #
#  ###     ##     ####   ##   #  ##  ##   #  ####   ##     #
#  ## #    ##     ##     ##    # ##  ##   #  ##     # #    #
#  ##  #   ##     ##     ##     ###  ##  #   ##     #  #   #
#  #####   #####  ###### ##     ###  ####    #####  #   #  #
#                                                          #
############################################################

######################
#                    #
#  Define functions  #
#                    #
######################

function display_usage()
{
	COLUMNS=$(tput cols)
	printf "%*s\n" $(($COLUMNS/2)) "~~~~~~~~~~~~~~~~~~";printf "%*s\n" $(($COLUMNS/2)) "|                |";printf "%*s\n" $(($COLUMNS/2)) "|  XOME-BLENDER  |";printf "%*s\n" $(($COLUMNS/2)) "|                |";printf "%*s\n" $(($COLUMNS/2)) "~~~~~~~~~~~~~~~~~~";printf "%*s\n" $((($COLUMNS/2)+29)) "____________________________________________________________________________";printf "%*s\n" $((($COLUMNS/2)+29)) "Xome-Blender is a tool for the generation of .bam files made up of different";printf "%*s\n" $((($COLUMNS/2)+29)) "samples. For example, you can simulate the subclonal architecture of tumoral";printf "%*s\n" $((($COLUMNS/2)-20)) "cells or the contamination."; printf "%*s\n" $((($COLUMNS/2)+29)) "----------------------------------------------------------------------------"; printf "%*s\n" $((($COLUMNS/2)+27)) " Usage: ./Xome-Blender.sh -f f1,f2,f3 -la label -sc 80 -p 10,40,50 -fc 100"; printf "%*s\n" $((($COLUMNS/2)+29)) "____________________________________________________________________________"; printf '\e[1;39m%*s\n\e[m' $((($COLUMNS/2)-4)) "Arguments:";printf "%*s" $((($COLUMNS/3)+3)) "-f,--files"; printf "%*s\n" $((($COLUMNS/3)-13)) "Input files comma separated";printf "%*s" $((($COLUMNS/3)+4)) "-la,--label"; printf "%*s\n" $((($COLUMNS/3)-14)) "The same label used previously in InXalizer"; printf '\e[8;39m%*s\e[m' $((($COLUMNS/3)+4)) "-la,--label";printf "%*s\n" $((($COLUMNS/3)-14)) "to generate the files to be merged"; printf "%*s" $((($COLUMNS/3)+16)) "-sc,--starting_coverage"; printf "%*s\n" $((($COLUMNS/3)-26)) "The coverage of the input bam (file.cov)";printf "%*s" $((($COLUMNS/3)+9)) "-p,--percentages"; printf "%*s\n" $((($COLUMNS/3)-19)) "The desidered percentage of each sample";printf "%*s" $((($COLUMNS/3)+13)) "-fc,--final_coverage"; printf "%*s\n" $((($COLUMNS/3)-23)) "The desidered coverage for the output bam";printf "%*s" $((($COLUMNS/3)+6)) "-v,--variants"; printf "%*s\n" $((($COLUMNS/3)-16)) "The vcf files generated with InXalizer"; printf "%*s" $((($COLUMNS/3)+3)) "-cnv,--cnv"; printf "%*s\n" $((($COLUMNS/3)-13)) "A CNV file genertaed with InXalizer"; printf "%*s" $((($COLUMNS/3)+2)) "-l,--list"; printf "%*s\n" $((($COLUMNS/3)-12)) "Automated Mode: Get a tab separated file" ;printf '\e[8;39m%*s\e[m' $((($COLUMNS/3)+2)) "-l,--list";printf "%*s\n" $((($COLUMNS/3)-11)) "containing different anlaysis. Every row "; printf '\e[8;39m%*s\e[m' $((($COLUMNS/3)+2)) "-l,--list";printf "%*s\n" $((($COLUMNS/3)-12)) "must contain all the options above"; printf "%*s" $((($COLUMNS/3)+5)) "-t,--threads"; printf "%*s\n" $((($COLUMNS/3)-15)) "The number of threads used for this analysis"; printf "%*s" $((($COLUMNS/3)+4)) "-o,--output"; printf "%*s\n" $((($COLUMNS/3)-14)) "Output directory. If omitted, generates a";printf '\e[8;39m%*s\e[m' $((($COLUMNS/3)+4)) "-o,--output"; printf "%*s\n" $((($COLUMNS/3)-14)) "results directory in the current position"; printf "%*s\n" $((($COLUMNS/2)+29)) "____________________________________________________________________________"; printf "%*s\n" $((($COLUMNS/2)+20)) "NB: The order of samples must be respected in each option!"; printf "%*s\n" $((($COLUMNS/2)+29)) "----------------------------------------------------------------------------"; printf "%*s\n" $((($COLUMNS/2)+29)) "Xome-Blender. Written by Roberto Semeraro, Department of Clinical and Speri-"; printf "%*s\n" $((($COLUMNS/2)+29)) "mental Medicine, University of Florence. For bug reports or suggestion write"; printf "%*s\n" $((($COLUMNS/2)-21)) "to robe.semeraro@gmail.com"
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
	vcftoolsExe=$(sed -n '3p' "$FilesFolder/deps.txt")
	PicardExe=$(echo "java -jar $PicardExeLine")
	XomeCounter=$(echo "$FilesFolder/Scripts/xome_counter")

}

startgraph()
{
	PriID=`echo " - SAMPLES = "${IDs//,/ }`
	Samplength=`echo ${#PriID}`
	LineLength=$(($Samplength+1))
	HalfLength=$(($LineLength/2))
	BorderLines=`for ((x = 0; x < $LineLength ; x++)); do   printf %s =; done`
	printf '\e[1;39m%*s\n\e[m' $(($HalfLength+19)) "   _    _                            "
	printf '\e[1;39m%*s\n\e[m' $(($HalfLength+19)) "  \ \  / /  ___   _           ___    "
	printf '\e[1;39m%*s\n\e[m' $(($HalfLength+19)) "   \ \/ /  / _ \ | |_______  / _ \   "
	printf '\e[1;39m%*s\n\e[m' $(($HalfLength+19)) "   / /\ \ | (_) ||  _   _  \|  __/   "
	printf '\e[1;39m%*s\n\e[m' $(($HalfLength+19)) "  /_/  \_\ \___/ |_| |_| |_| \___|   "
	printf '\e[1;39m%*s\n\e[m' $(($HalfLength+19)) " ___  _                 _            "
	printf '\e[1;39m%*s\n\e[m' $(($HalfLength+19)) "| . \| |  ___  _       | |  ___  _ _ "
	printf '\e[1;39m%*s\n\e[m' $(($HalfLength+19)) "|  _/| | / _ \| |___  _| | / _ \| '_/"
	printf '\e[1;39m%*s\n\e[m' $(($HalfLength+19)) "| . \| ||  __/|  _  \/ . ||  __/| |  "
	printf '\e[1;39m%*s\n\e[m' $(($HalfLength+19)) "|___/|_| \___||_| |_|\___| \___||_|  "
	echo $BorderLines
	echo -e " - "'\e[31mSAMPLES\e[0m' = ${IDs//,/ }
	echo -e " - "'\e[31mPERCENTAGES\e[0m' = ${PERCENTAGES//,/% }'%'
	echo $BorderLines
	printf '\e[1;39m%*s\n\e[m' $((($LineLength/2)+5)) "Blending:"
	printf "%*s\n"
}

function subsampling()
{
	if [[ ! -f ${FILE[$e]} ]]; then
		echo -e '\e[31mFile '"${FILE[$e]}"' not found!\e[0m' && exit
	fi
	New_percentage=$(echo $(printf "%.0f\n" $(perl -E "say ${FINALCOVERAGE}/${STARTINGCOVERAGE}*${PERCENTAGE[$e]}")))
	if [[ $New_percentage -eq 0 ]]; then
		New_percentage=1
	fi
	if [[ ! -f "$InvisibleDir/"${ID[$e]}"_"${PERCENTAGE[$e]}"%.bam" ]]; then 
		if [[ ${#New_percentage} == 1 ]]; then
			MCSBait=$(echo "samtools view -s $(( ( RANDOM % 100 )  + 1 )).0$New_percentage -b ${FILE[$e]} > "$InvisibleDir/"${ID[$e]}"_"${PERCENTAGE[$e]}"%.bam" 2>/dev/null |:")
			MCBait=$(echo $MCSBait | awk '{split($0,a,/ >/);$0=a[1]}1')
			eval $MCSBait &
			if [[ "$Flag" -eq ${THREADS} ]]; then
                        	SamPid=`ps -ef | grep -w "$MCBait" | sed '/grep/d' | awk '{print $3}'`
                        	echo $SamPid
                        	while [[ -z $SamPid ]]; do sleep 0.1; SamPid=`ps -ef | grep -w "$MCBait" | sed '/grep/d' | awk '{print $3}'`; done
                        	SamWaitStrng=$( echo "while [ -e /proc/$SamPid ]; do sleep 0.1; done" )
                        	eval "$SamWaitStrng"
                        	Flag=0
                        fi
		else
			MCSBait=$(echo "samtools view -s $(( ( RANDOM % 100 )  + 1 )).$New_percentage -b ${FILE[$e]} > "$InvisibleDir/"${ID[$e]}"_"${PERCENTAGE[$e]}"%.bam" 2>/dev/null |:")
			MCBait=$(echo $MCSBait | awk '{split($0,a,/ >/);$0=a[1]}1')
			eval $MCSBait &
			if [[ "$Flag" -eq ${THREADS} ]]; then
                                SamPid=`ps -ef | grep -w "$MCBait" | sed '/grep/d' | awk '{print $3}'`
                                while [[ -z $SamPid ]]; do sleep 0.1; SamPid=`ps -ef | grep -w "$MCBait" | sed '/grep/d' | awk '{print $3}'`; done
                                SamWaitStrng=$( echo "while [ -e /proc/$SamPid ]; do sleep 0.1; done" )
                                eval "$SamWaitStrng"
                                Flag=0
                        fi
		fi
	fi
}

function progress_prep () {
    s=0.5;
    f=0.25;
    echo -ne "\r";
    while true; do
           sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [                  ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [=>                ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [==>               ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [===>              ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [====>             ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [=====>            ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [======>           ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [=======>          ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [========>         ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [=========>        ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [==========>       ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [===========>      ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [============>     ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [=============>    ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [==============>   ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [===============>  ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [================> ] Elapsed: ${s} secs.";
           sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [=================>] Elapsed: ${s} secs.";
    done;	
}

function progress_final () {
    s=0.5;
    f=0.25;
    echo -ne "\r";
    while true; do
	   sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [                 ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [=>               ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [==>              ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [===>             ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [====>            ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [=====>           ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [======>          ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [=======>         ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [========>        ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [=========>       ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [==========>      ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [===========>     ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [============>    ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [=============>   ] Elapsed: ${s} secs." \
	&& sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [==============>  ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [===============> ] Elapsed: ${s} secs.";
           sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [================>] Elapsed: ${s} secs.";
    done;	
}

function CNV()
{
	if [[ $INT_SAMP == *"-"* ]]; then
		IFS='-' read -ra INT_SAMPLES <<< "$INT_SAMP"
		for t in "${INT_SAMPLES[@]}"
		do
			INT_SAMP=`echo $t`
			INT_FILE=`ls *.bam | grep "$INT_SAMP""_"`
			if [[ ! -e "${INT_FILE}.bai" ]] ; then
				echo -e '\e[31mThe bam index file is missing!\e[0m' && exit
			fi
			if [[ ! -f $INT_SAMP"."$INT_CHR"."$INT_START"."$INT_STOP".sorted.NonEthR.50.bam" ]]; then
				AddIntervals & # For many samples
			fi
		done
	else

		INT_FILE=`ls *.bam | grep "$INT_SAMP""_"`
		if [[ ! -e "${INT_FILE}.bai" ]] ; then
			echo -e '\e[31mThe bam index file is missing!\e[0m' && exit
		fi
		if [[ ! -f $INT_SAMP"."$INT_CHR"."$INT_START"."$INT_STOP".sorted.NonEthR.50.bam" ]]; then
			AddIntervals # For one sample
		fi
	fi
	wait
}

function progress_cnv () {
    s=0.5;
    f=0.25;
    echo -ne "\r";
    while true; do
           sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [        ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [=>      ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [==>     ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [===>    ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [====>   ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [=====>  ] Elapsed: ${s} secs." \
        && sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [======> ] Elapsed: ${s} secs.";
           sleep $f && s=`echo ${s} + ${f} + ${f} | bc` && echo -ne "\r      [=======>] Elapsed: ${s} secs.";
    done;	
}

function AddIntervals()
{
	samtools view -b -h $INT_FILE $INT_CHR":"$INT_START"-"$INT_STOP > $INT_SAMP"."$INT_CHR"."$INT_START"."$INT_STOP".bam"
	wait
	samtools index $INT_SAMP"."$INT_CHR"."$INT_START"."$INT_STOP".bam"
	wait
	
	VariantFile=`echo $INT_SAMP".vcf"`
	VARFILE=`printf -- '%s\n' "${VARIANT[@]}" | grep "$VariantFile"`
	if [[ ! -f "$VARFILE" ]]; then
		echo -e '\e[31mFile '"$VARFILE"' not found!\e[0m' && exit
	fi

	VarLength=`$vcftoolsExe --vcf "$VARFILE" --chr $INT_CHR --from-bp $INT_START --to-bp $INT_STOP --recode --stdout -c | grep -v "#" | wc -l`

	if [[ $VarLength -gt 0 ]]; then
		$vcftoolsExe --vcf "$VARFILE" --chr $INT_CHR --from-bp $INT_START --to-bp $INT_STOP --recode --stdout -c | stdbuf -oL grep -v "#" | stdbuf -oL awk 'BEGIN {FS = "\t"} ; {print $2}' | while IFS= read -r line; do samtools view -h $INT_SAMP"."$INT_CHR"."$INT_START"."$INT_STOP".bam" $INT_CHR":"$line"-"$line | samtools fillmd -e - $INT_REF 2>/dev/null | grep -v "^@" | awk -v pos=$line 'BEGIN {OFS = FS = "\t" } ; {system(" '$XomeCounter' "$4" "pos" "$10" "$6" "$1)}' >> "."$INT_SAMP"."$INT_CHR"."$INT_START"."$INT_STOP".varreads" 2>/dev/null; done ; wait ; sort -u "."$INT_SAMP"."$INT_CHR"."$INT_START"."$INT_STOP".varreads" > "."$INT_SAMP"."$INT_CHR"."$INT_START"."$INT_STOP".svarreads" # Extract and sort reads name of reads containing the alternative allele
	
		$vcftoolsExe --vcf "$VARFILE" --chr $INT_CHR --from-bp $INT_START --to-bp $INT_STOP --recode --stdout -c | stdbuf -oL grep -v "#" | stdbuf -oL awk 'BEGIN {FS = "\t"} ; {print $2}' | while IFS= read -r line; do samtools view -h $INT_SAMP"."$INT_CHR"."$INT_START"."$INT_STOP".bam" $INT_CHR":"$line"-"$line | samtools fillmd -e - $INT_REF 2>/dev/null | grep -v "^@" | awk '{print $1}' | sort -u >> "."$INT_SAMP"."$INT_CHR"."$INT_START"."$INT_STOP".allvarreads"; done # Extract all reads name

		rm "."$INT_SAMP"."$INT_CHR"."$INT_START"."$INT_STOP".varreads"
	
		$PicardExe FilterSamReads I= $INT_SAMP"."$INT_CHR"."$INT_START"."$INT_STOP".bam" FILTER=includeReadList QUIET=TRUE READ_LIST_FILE= "."$INT_SAMP"."$INT_CHR"."$INT_START"."$INT_STOP".svarreads" WRITE_READS_FILES=false O= $INT_SAMP"."$INT_CHR"."$INT_START"."$INT_STOP".EthR.bam" 2>/dev/null # Generate Ethero Variant Reads File
		wait
		samtools sort -T /tmp/$INT_SAMP"."$INT_CHR"."$INT_START"."$INT_STOP".sorted.EthR.bam" -o $INT_SAMP"."$INT_CHR"."$INT_START"."$INT_STOP".sorted.EthR.bam" $INT_SAMP"."$INT_CHR"."$INT_START"."$INT_STOP".EthR.bam" # Sort Ethero Variant Reads File
		wait
		rm "."$INT_SAMP"."$INT_CHR"."$INT_START"."$INT_STOP".svarreads"
		rm $INT_SAMP"."$INT_CHR"."$INT_START"."$INT_STOP".EthR.bam"

		$PicardExe FilterSamReads I= $INT_SAMP"."$INT_CHR"."$INT_START"."$INT_STOP".bam" FILTER=excludeReadList QUIET=TRUE READ_LIST_FILE= "."$INT_SAMP"."$INT_CHR"."$INT_START"."$INT_STOP".allvarreads" WRITE_READS_FILES=false O= $INT_SAMP"."$INT_CHR"."$INT_START"."$INT_STOP".NonEthR.bam"  2>/dev/null # Generate Non Ethero Reads File
		rm "."$INT_SAMP"."$INT_CHR"."$INT_START"."$INT_STOP".allvarreads"
		wait
		samtools sort -T /tmp/$INT_SAMP"."$INT_CHR"."$INT_START"."$INT_STOP".sorted.NonEthR.bam" -o $INT_SAMP"."$INT_CHR"."$INT_START"."$INT_STOP".sorted.NonEthR.bam" $INT_SAMP"."$INT_CHR"."$INT_START"."$INT_STOP".NonEthR.bam"
		rm $INT_SAMP"."$INT_CHR"."$INT_START"."$INT_STOP".NonEthR.bam" # Remove non ethero reads full bam	
		wait
		samtools index $INT_SAMP"."$INT_CHR"."$INT_START"."$INT_STOP".sorted.EthR.bam" 
		samtools index $INT_SAMP"."$INT_CHR"."$INT_START"."$INT_STOP".sorted.NonEthR.bam"
	else
		mv $INT_SAMP"."$INT_CHR"."$INT_START"."$INT_STOP".bam" $INT_SAMP"."$INT_CHR"."$INT_START"."$INT_STOP".sorted.NonEthR.bam"
		mv $INT_SAMP"."$INT_CHR"."$INT_START"."$INT_STOP".bam.bai" $INT_SAMP"."$INT_CHR"."$INT_START"."$INT_STOP".sorted.NonEthR.bam.bai"
	fi
	samtools view -s $(( ( RANDOM % 100 )  + 1 )).50 -b $INT_SAMP"."$INT_CHR"."$INT_START"."$INT_STOP".sorted.NonEthR.bam"  > $INT_SAMP"."$INT_CHR"."$INT_START"."$INT_STOP".sorted.NonEthR.50.bam" # SubSample the 50% of non ethero reads
	samtools index $INT_SAMP"."$INT_CHR"."$INT_START"."$INT_STOP".sorted.NonEthR.50.bam"
	wait
	if [[ $INT_TYPE == "Del" ]]; then
		if [[ $VarLength -gt 0 ]]; then
			samtools view $INT_SAMP"."$INT_CHR"."$INT_START"."$INT_STOP".bam" | awk '{print $1}' >> "."$INT_SAMP".remove"
		else
			samtools view $INT_SAMP"."$INT_CHR"."$INT_START"."$INT_STOP".sorted.NonEthR.bam" | awk '{print $1}' >> "."$INT_SAMP".remove"
		fi
	fi
	rm $INT_SAMP"."$INT_CHR"."$INT_START"."$INT_STOP".sorted.NonEthR.bam"  # Remove sorted non ethero reads full bam
	rm $INT_SAMP"."$INT_CHR"."$INT_START"."$INT_STOP".sorted.NonEthR.bam.bai"
	if [[ $VarLength -gt 0 ]]; then
		rm $INT_SAMP"."$INT_CHR"."$INT_START"."$INT_STOP".bam" # Remove bam region
		rm $INT_SAMP"."$INT_CHR"."$INT_START"."$INT_STOP".bam.bai"
	fi
}

function DelGenerator()
{
	sort -u "."$INT_SAMP".remove" > "."$INT_SAMP".sorted.remove"
	rm "."$INT_SAMP".remove"
	INT_FILE=`ls *.bam | grep "$INT_SAMP""_"`
	eval $PicardExe FilterSamReads I= $INT_FILE FILTER=excludeReadList QUIET=TRUE READ_LIST_FILE= "."$INT_SAMP".sorted.remove" WRITE_READS_FILES=false VALIDATION_STRINGENCY=LENIENT O= $INT_FILE"_2" 2>/dev/null
	echo "Finished!" > $INT_SAMP"DelFuncEnd"
	rm $INT_FILE
	mv $INT_FILE"_2" $INT_FILE
	rm "."$INT_SAMP".sorted.remove"
}

function MultipleMixing()
{
	if [[ -z ${FINALCOVERAGE} ]] ; then
		echo -e '\e[31mFinal coverage is missing!\e[0m' && exit
	fi

	if [[ -z ${FILES} ]] ; then
		echo -e '\e[31mNo files selected!\e[0m' && exit
	fi

	if [[ -z ${LABEL} ]] ; then
		echo -e '\e[31mNo label assigned!\e[0m' && exit
	fi

	if [[ -z ${PERCENTAGES} ]] ; then
		echo -e '\e[31mNo percentages selected!\e[0m' && exit
	fi

	if [[ -z ${STARTINGCOVERAGE} ]] ; then
		echo -e '\e[31mNo starting coverage selected!\e[0m' && exit
	fi

	if [[ -z ${VARIANTS} ]] ; then
		echo -e '\e[31mNo variants file selected!\e[0m' && exit
	fi

	BarGrep=""

	IFS=',' read -ra FILE <<< "$FILES"
	IFS=',' read -ra PERCENTAGE <<< "${PERCENTAGES}"
	IFS=',' read -ra VARIANT <<< "${VARIANTS}"

	if [[ ${#VARIANT[@]} != $((${#PERCENTAGE[@]}-1)) ]]; then
		echo -e '\e[31mWrong number of variants file!\e[0m' && exit
	fi
	
	MaxPercentage=`IFS=$'\n'; echo "${PERCENTAGE[*]}" | sort -nr | head -n1`

	### Label Check ###

	if [[ ${LABEL} == *"-"* ]]; then
  		echo -e '\e[31mHyphens are not allowed for labels. Please, replace it.\e[0m' && exit
	fi

	### IDs Gen ###

	ControlLabel=`echo ${LABEL}"_Control"`
	SubClLabels=""
	for i in $(seq 1 $((${#PERCENTAGE[@]}-1)))
	do
		SubClLabels=$SubClLabels,${LABEL}"_Subclone_"$i
	done

	IDs=$ControlLabel$SubClLabels
	IFS=',' read -ra ID <<< "$IDs"

	if [[ ${#PERCENTAGE[@]} != ${#FILE[@]} || ${#PERCENTAGE[@]} != ${#ID[@]} ]]; then
		echo -e '\e[31mThe number of values per input is different, fix it!\e[0m' && exit
	fi

	TotalPercent=$(echo $(perl -E "say ${PERCENTAGES//,/ + }"))
	if [[ $TotalPercent != 100 ]]; then 
		echo -e '\e[31mThe sum of your percentages is not equal to 100, fix it!\e[0m' && exit
	fi
	End=$(echo $(perl -E "say ${#PERCENTAGE[@]}-1"))
	
	for i in $(seq 0 $((${#PERCENTAGE[@]}-1)))
	do
		New_coverage=$(echo $(printf "%.0f\n" $(perl -E "say ${FINALCOVERAGE}*(${PERCENTAGE[$i]}/100)")))
		if [[ "$New_coverage" -eq 0 ]]; then
			New_coverage=1
		fi
		if [[ "$New_coverage" -gt "${STARTINGCOVERAGE}" ]]; then
			echo -e '\e[31mI cannot reach the desidered coverage for the output with this inputs!\e[0m' && exit
		fi
	done

	startgraph

	WorkDir=$(pwd)

	cd $WorkDir

	InvisibleDir=.NewCoverageBam$RANDOM

	if [[ ! -d $InvisibleDir ]]; then
		mkdir $InvisibleDir
	fi

	### calculate percentages ###

	echo -e " 1) - "'\e[37mSubclone preparation\e[0m'
	progress_prep &
	BarPid=$!	
	
	Flag=0

	for e in $(seq 0 $End)
	do
		Flag=$(($Flag+1))
		subsampling
	done
	
	GrepSubClonesIndexes=`ls | grep -c "%.bam$"`
	while [[ $GrepSubClonesIndexes -lt ${#ID[@]} ]]; do
		cd $InvisibleDir
		GrepSubClonesIndexes=`ls | grep -c "%.bam$"`
		cd $WorkDir
		sleep  0.1
	done

	while [ -z $BarGrep ]; do
		cd $InvisibleDir
		BarGrep=`ls *%.bam 2>/dev/null | grep "$MaxPercentage%" | head -n 1`
		cd $WorkDir
		sleep  0.1
	done
	
	MaxSample=${BarGrep%*_*}	

	while [ ! -f $InvisibleDir"/"$BarGrep ]; do sleep 0.1;done	

	while [ $(( $(date +%s) - $(stat -c %Y $InvisibleDir"/"$BarGrep) )) -lt 20 ]; do sleep 1; done
	if [[ ! -z ${CNVS} ]]; then
		cd $InvisibleDir
		IndexingVec=`ls -1 *%.bam | tr "\n" " "`
		for i in $IndexingVec
		do
			if [[ ! -e $i.bai || $i.bai -ot $i  ]]; then
				samtools index $i 
			fi
		done
		cd $WorkDir
	fi	
	if [[ ! -z ${CNVS} ]]; then
		while [ ! -f $InvisibleDir"/"$BarGrep".bai" ]; do sleep 0.1;done
		while [ $(( $(date +%s) - $(stat -c %Y $InvisibleDir"/"$BarGrep".bai") )) -lt 20 ]; do sleep 1; done		
		echo -ne "\r\033[K      [=================>] - Finished." && kill $BarPid 
		wait $BarPid 2>/dev/null
	else
		while [ ! -f $InvisibleDir"/"$BarGrep ]; do sleep 0.1;done
		while [ $(( $(date +%s) - $(stat -c %Y $InvisibleDir"/"$BarGrep) )) -lt 20 ]; do sleep 1; done		
		echo -ne "\r\033[K      [=================>] - Finished." && kill $BarPid 
		wait $BarPid 2>/dev/null
	fi

	### Add CNV ###

	if [[ ! -z ${CNVS} ]]; then
		if [[ ! -f ${CNVS} ]]; then
			echo -e '\e[31mFile '"${CNVS}"' not found!\e[0m' && exit
		fi	
		if [[ $THREADS -ge 3 ]]; then
			CNVTHREADS=$(echo $(perl -E "say (int(${THREADS}/3) )"))
		else
			CNVTHREADS=1
		fi
		echo -e "\n 2) - "'\e[37mCNV adding\e[0m'
		progress_cnv &
		BarPid2=$!
		CNVFlag=0
		### Deleted sample Number ###
		DelOccurences=`grep Del ${CNVS} | cut -f5 | grep -o -n "-" | cut -d : -f 1 | uniq -c | awk '{print $1}' | awk -v idx=1 'NR==1 || $idx>max{max=$idx} END{print max}'`
		DeletedSamples=$(($DelOccurences+1))
		array=()
		getArray() {
		i=0
		while read line
			do
				array[i]=$line
				i=$(($i + 1))
			done < $1
		}
		Couple=0
		getArray ${CNVS}
		for e in "${array[@]}"
		do
			INT_CHR=`echo $e | awk '{print $1}'`
			INT_START=`echo $e | awk '{print $2}'`
			INT_STOP=`echo $e | awk '{print $3}'`
			INT_TYPE=`echo $e | awk '{print $4}'`
			INT_SAMP=`echo $e | awk '{print $5}'`
			INT_REF=`echo $e | awk '{print $6}'`
			if [[ $INT_SAMP == *"-"* ]]; then
				Couple=$(($Couple+1))
				needle="-"
				number_of_occurrences=$(grep -o "$needle" <<< "$INT_SAMP" | wc -l)
				CNVev=$((1+$number_of_occurrences))
			else
				CNVev=1
			fi
			CNVFlag=$(($CNVFlag+$CNVev))
			if [[ pwd == $WorkDir ]]; then
				cd $InvisibleDir
			else
				cd $WorkDir
				cd $InvisibleDir
			fi
			CNV & # Function
			Pido=`echo $!`
			if [[ "$CNVFlag" -ge ${CNVTHREADS} ]]; then
				if [[ $CNVFlag -gt ${CNVTHREADS} ]]; then
					CNVFlag=`echo ${CNVTHREADS}`
				fi
				WaitStrng=$( echo "while [ -e /proc/$Pido ]; do sleep 0.1; done" )
				eval "$WaitStrng"
				CNVFlag=$(($CNVFlag-$CNVev))
			fi
			cd $WorkDir
		done

		TotalIntervals=$((${#array[@]}+$Couple))

		cd $InvisibleDir
		GrepIntervals=`ls | grep -c "NonEthR.50.bam.bai"`
		while [[ $GrepIntervals -lt $TotalIntervals ]]; do
			GrepIntervals=`ls | grep -c "NonEthR.50.bam.bai"`
			sleep  0.1
		done

		DelFlag=0
		
		for k in $(seq 1 $((${#ID[@]}-1)))
		do
			DelFlag=$(($DelFlag+1))
			INT_SAMP=${ID[$k]}
			DelGenerator &
			if [[ $DelFlag -eq ${THREADS} ]]; then
				Delcheck=`ls -1 | grep -w $INT_SAMP"DelFuncEnd"`
				while [[ ! -f $Delcheck ]]; do sleep 0.1; Delcheck=`ls -1 | grep -w $INT_SAMP"DelFuncEnd"`; done
				DelFlag=0
			fi
		done
		cd $WorkDir
		DeletingProcesses=`cd $InvisibleDir; ls -1 | grep -c "DelFuncEnd"; cd $WorkDir`
		while [[ $DeletingProcesses -lt $DeletedSamples ]]; do sleep 0.1; DeletingProcesses=`cd $InvisibleDir; ls -1 | grep -c "DelFuncEnd"; cd $WorkDir`; done
		echo -ne "\r\033[K      [========>] - Finished." && kill $BarPid2
		wait $BarPid2 2>/dev/null
	fi
	
	### create good read group ###

	echo -e "\n 3) - "'\e[37mSample finalization\e[0m'
	progress_final &
	BarPid3=$!
	IDVec=$(echo ${IDs//,/_})
	PercVec=$(echo ${PERCENTAGES//,/_})
	samtools view -H "$InvisibleDir/"${ID[0]}"_"${PERCENTAGE[0]}"%.bam" | grep "@RG" > "$InvisibleDir/.temp.rg"
	
	sed -e 's/SM[^[:space:]]*/SM:'"$IDVec"'_'"$PercVec"'/g' "$InvisibleDir/.temp.rg" > "$InvisibleDir/.rg.txt"
	
	### merge new files ###

	cd $InvisibleDir
	rm .temp.rg
	
	ls -1 *.bam > Merge.input

	MergeInputPercs=`grep -c "%" Merge.input`
	while [[ $MergeInputPercs -lt ${#PERCENTAGE[@]} ]]; do sleep 0.1 ; ls -1 *.bam > Merge.input ; MergeInputPercs=`grep -c "%" Merge.input`; done

	samtools merge -h ".rg.txt" -cp  $IDVec"_"$PercVec".unsrt.bam" -b Merge.input

	### sort final file ###

	samtools sort -T /tmp/$IDVec"_"$PercVec".bam" -o $IDVec"_"$PercVec".bam" $IDVec"_"$PercVec".unsrt.bam" -@ ${THREADS} 2>/dev/null

	rm $IDVec"_"$PercVec".unsrt.bam"
	rm Merge.input

	samtools index $IDVec"_"$PercVec".bam"

	cd $WorkDir

	### remove single files ###

	for z in $(seq 0 $End)
	do
		rm "$InvisibleDir/"${ID[$z]}"_"${PERCENTAGE[$z]}"%.bam"
		if [[ ! -z ${CNVS} ]]; then
			rm "$InvisibleDir/"${ID[$z]}"_"${PERCENTAGE[$z]}"%.bam.bai"
		fi
	done

	if [[ -z ${OUTPUT} ]] ; then
		if [[ ! -d XomeBlender_results ]]; then
			mkdir XomeBlender_results
		fi
		cp "$InvisibleDir/"$IDVec"_"$PercVec".bam" "XomeBlender_results"
		cp "$InvisibleDir/"$IDVec"_"$PercVec".bam.bai" "XomeBlender_results"

	### delete single files ###

		if [[ ! -e "XomeBlender_results/"$IDVec"_"$PercVec".bam" ]]; then
			echo "File does not exist!"
		else
			rm "$InvisibleDir/"$IDVec"_"$PercVec".bam"
			rm "$InvisibleDir/"$IDVec"_"$PercVec".bam.bai"
			rm "$InvisibleDir/.rg.txt"
		fi
		if [[ ! -z ${CNVS} ]]; then
			cd $InvisibleDir
			rm *".bam"
			rm *".bai"
			rm "out.log"
			rm *"DelFuncEnd"
			cd $WorkDir
		fi 
	else
		if [[ ! -d ${OUTPUT} ]]; then
			mkdir ${OUTPUT}
		fi
		cp "$InvisibleDir/"$IDVec"_"$PercVec".bam" ${OUTPUT}
		cp "$InvisibleDir/"$IDVec"_"$PercVec".bam.bai" ${OUTPUT}

	### delete single files ###

		if [[ ! -e ${OUTPUT}"/"$IDVec"_"$PercVec".bam" ]]; then
			echo "File does not exist!"
		else
			rm "$InvisibleDir/"$IDVec"_"$PercVec".bam"
			rm "$InvisibleDir/"$IDVec"_"$PercVec".bam.bai"
			rm "$InvisibleDir/.rg.txt"
		fi
		if [[ ! -z ${CNVS} ]]; then
			cd $InvisibleDir
                        rm *".bam"
                        rm *".bai"
			rm *"DelFuncEnd"
			rm "out.log"
			cd $WorkDir
                fi
	fi
if [ `find $InvisibleDir -prune -empty` ]
then
	rm -rf $InvisibleDir
fi
echo -ne "\r\033[K      [================>] - Finished." && kill $BarPid3
wait $BarPid3 2>/dev/null
}

##########################
# Multiple Contamination #
##########################

#######################
#                     #
#  Parameters Config  #
#                     #
#######################

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
	-f|--files)
	FILES="$1"
	shift
	;;
	-la|--label)
	LABEL="$1"
	shift
	;;
	-p|--percentages)
	PERCENTAGES="$1"
	shift
	;;
	-sc|--starting_coverage)
	STARTINGCOVERAGE="$1"
	shift
	;;
	-fc|--final_coverage)
	FINALCOVERAGE="$1"
	shift
	;;
	-v|--variant_file)
	VARIANTS="$1"
	shift
	;;
	-cnv|--cnv)
	CNVS="$1"
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
	-l|--list)
	LIST="$1"
	shift
	;;
	esac
done

### Threads setting ###

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

######################
# Dependencies check #
######################

dependencies_check

if [[ -z ${LIST} ]] ; then
	COLUMNS=$(tput cols)
	MultipleMixing # Start analysis
else
	COLUMNS=$(tput cols)

	### generate array ####

	array=()
	getArray() {
		i=0
		while read line
		do
			array[i]=$line
			i=$(($i + 1))
		done < $1
	}
	getArray ${LIST}
	for e in "${array[@]}"
	do
		FILES=`echo $e | awk '{print $1}'`
		LABEL=`echo $e | awk '{print $2}'`
		STARTINGCOVERAGE=`echo $e | awk '{print $3}'`
		PERCENTAGES=`echo $e | awk '{print $4}'`
		FINALCOVERAGE=`echo $e | awk '{print $5}'`
		VARIANTS=`echo $e | awk '{print $6}'`
		CNVS=`echo $e | awk '{print $7}'`
		MultipleMixing # Start analysis
	done
fi

echo -e "\n      "'\e[37mDone\e[0m'
