#!/bin/bash

######################################################
#                                                    #
#  ##    ##     #####    ###        ###  #########   #
#   ##  ##    ##    ##   ####      ####  ##          # 
#    ####    ##      ##  ## ##    ## ##  ######      #
#    ####    ##      ##  ##   ####   ##  ######      # 
#   ##  ##    ##    ##   ##          ##  ##          #
#  ##    ##    #####     ##          ##  #########   #
#                                                    #   
######################################################

######################
#                    #
#  Define functions  #
#                    #
######################

function display_usage()
{
	COLUMNS=$(tput cols)
	printf "%*s\n" $(($COLUMNS/2)) "~~~~~~~~~~~~~~~~~~";printf "%*s\n" $(($COLUMNS/2)) "|                |";printf "%*s\n" $(($COLUMNS/2)) "|  XOME-BLENDER  |";printf "%*s\n" $(($COLUMNS/2)) "|                |";printf "%*s\n" $(($COLUMNS/2)) "~~~~~~~~~~~~~~~~~~";printf "%*s\n" $((($COLUMNS/2)+29)) "____________________________________________________________________________";printf "%*s\n" $((($COLUMNS/2)+29)) "Xome-Blender is a tool for the generation of .bam files made up of different";printf "%*s\n" $((($COLUMNS/2)+29)) "samples. For example, you can simulate the subclonal architecture of tumoral";printf "%*s\n" $((($COLUMNS/2)-20)) "cells or the contamination."; printf "%*s\n" $((($COLUMNS/2)+29)) "----------------------------------------------------------------------------"; printf "%*s\n" $((($COLUMNS/2)+22)) "Usage: $0 [-f file1,file2,file3] [-i id1,id2,id3]"; printf "%*s\n" $((($COLUMNS/2)+29)) "____________________________________________________________________________"; printf '\e[1;39m%*s\n\e[m' $((($COLUMNS/2)-4)) "Arguments:";printf "%*s" $((($COLUMNS/3)+3)) "-f,--files"; printf "%*s\n" $((($COLUMNS/3)-13)) "Input files comma separated";printf "%*s" $((($COLUMNS/3))) "-i,--id"; printf "%*s\n" $((($COLUMNS/3)-10)) "Samples ID's; same order of files";printf "%*s" $((($COLUMNS/3)+7)) "-c,--coverages"; printf "%*s\n" $((($COLUMNS/3)-17)) "The native coverage of each bam";printf "%*s" $((($COLUMNS/3)+9)) "-p,--percentages"; printf "%*s\n" $((($COLUMNS/3)-19)) "The desidered percentage of each sample";printf "%*s" $((($COLUMNS/3)+12)) "-tc,--totalcoverage"; printf "%*s\n" $((($COLUMNS/3)-22)) "The desidered coverage for the output bam";printf "%*s" $((($COLUMNS/3)+6)) "-v,--variants"; printf "%*s\n" $((($COLUMNS/3)-16)) "The variants file generated with InXalizer";printf "%*s" $((($COLUMNS/3)+2)) "-l,--list"; printf "%*s\n" $((($COLUMNS/3)-12)) "Automated Mode: Get a tab separated file" ;printf '\e[8;39m%*s\e[m' $((($COLUMNS/3)+2)) "-l,--list";printf "%*s\n" $((($COLUMNS/3)-11)) "containing different anlaysis. Every row "; printf '\e[8;39m%*s\e[m' $((($COLUMNS/3)+2)) "-l,--list";printf "%*s\n" $((($COLUMNS/3)-12)) "must contain all the options above";printf "%*s" $((($COLUMNS/3)+4)) "-o,--output"; printf "%*s\n" $((($COLUMNS/3)-14)) "Output directory. If omitted, generates a";printf '\e[8;39m%*s\e[m' $((($COLUMNS/3)+4)) "-o,--output"; printf "%*s\n" $((($COLUMNS/3)-14)) "results directory in the current position"; printf "%*s\n" $((($COLUMNS/2)+29)) "____________________________________________________________________________"; printf "%*s\n" $((($COLUMNS/2)+20)) "NB: The order of samples must be respected in each option!"; printf "%*s\n" $((($COLUMNS/2)+29)) "----------------------------------------------------------------------------"; printf "%*s\n" $((($COLUMNS/2)+29)) "Xome-Blender. Written by Roberto Semeraro, Department of Clinical and Speri-"; printf "%*s\n" $((($COLUMNS/2)+29)) "mental Medicine, University of Florence. For bug reports or suggestion write"; printf "%*s\n" $((($COLUMNS/2)-19)) "to roberto.semeraro@unifi.it"
}

function MultipleMixing()
{
	if [[ -z ${TOTALCOVERAGE} ]] ; then
		echo -e '\e[31mTotal coverage is missing!\e[0m' && exit
	fi

	if [[ -z ${FILES} ]] ; then
		echo -e '\e[31mNo files selected!\e[0m' && exit
	fi

	if [[ -z ${IDs} ]] ; then
		echo -e '\e[31mNo IDs assigned!\e[0m' && exit
	fi

	if [[ -z ${PERCENTAGES} ]] ; then
		echo -e '\e[31mNo percentages selected!\e[0m' && exit
	fi

	if [[ -z ${COVERAGES} ]] ; then
		echo -e '\e[31mNo coverages selected!\e[0m' && exit
	fi

	if [[ -z ${VARIANTS} ]] ; then
		echo -e '\e[31mNo variants file selected!\e[0m' && exit
	fi

	IFS=',' read -ra FILE <<< "$FILES"
	IFS=',' read -ra ID <<< "$IDs"
	IFS=',' read -ra PERCENTAGE <<< "${PERCENTAGES}"
	IFS=',' read -ra COVERAGE <<< "${COVERAGES}"

	if [[ ${#PERCENTAGE[@]} != ${#FILE[@]} || ${#PERCENTAGE[@]} != ${#COVERAGE[@]} || ${#PERCENTAGE[@]} != ${#ID[@]} ]]; then
		echo -e '\e[31mThe number of values per input is different, fix it!\e[0m' && exit
	fi
	TotalPercent=$(echo $(perl -E "say ${PERCENTAGES//,/ + }"))
	if [[ $TotalPercent != 100 ]]; then 
		echo -e '\e[31mThe sum of your percentages is not equal to 100, try again!\e[0m' && exit
	fi
	End=$(echo $(perl -E "say ${#COVERAGE[@]}-1"))
	for i in $(seq 0 $End)
	do
		New_coverage=$(echo $(printf "%.0f\n" $(perl -E "say ${TOTALCOVERAGE}*(${PERCENTAGE[$i]}/100)")))
		if [[ "$New_coverage" -gt "${COVERAGE[$i]}" ]]; then
			echo -e '\e[31mI cannot reach the desidered coverage for the output with this inputs!\e[0m' && exit
		fi
	done

	echo -e " - "'\e[31mSAMPLES\e[0m' = ${IDs//,/ }
	echo -e " - "'\e[31mPERCENTAGES\e[0m' = ${PERCENTAGES//,/% }'%'

	WorkDir=$(pwd)

	cd $WorkDir

	InvisibleDir=.NewCoverageBam$RANDOM

	if [[ ! -d $InvisibleDir ]]; then
		mkdir $InvisibleDir
	fi

	### calculate percentages ###

	for e in $(seq 0 $End)
	do
		New_percentage=$(echo $(printf "%.0f\n" $(perl -E "say ${TOTALCOVERAGE}/${COVERAGE[$e]}*${PERCENTAGE[$e]}")))
		if [[ ${#New_percentage} == 1 ]]; then
			samtools view -s $(( ( RANDOM % 100 )  + 1 )).0$New_percentage -b ${FILE[$e]} > "$InvisibleDir/"${ID[$e]}"_"${PERCENTAGE[$e]}"%.bam"
		else
			samtools view -s $(( ( RANDOM % 100 )  + 1 )).$New_percentage -b ${FILE[$e]} > "$InvisibleDir/"${ID[$e]}"_"${PERCENTAGE[$e]}"%.bam"
		fi
	done

	wait

	### create good read group ###

	IDVec=$(echo ${IDs//,/_})
	PercVec=$(echo ${PERCENTAGES//,/_})		
	NumberOfFileds=$(samtools view -H "$InvisibleDir/"${ID[0]}"_"${PERCENTAGE[0]}"%.bam" | grep "@RG" | awk 'BEGIN {FS="\t"} { print NF }')
	SMPositionField=$(samtools view -H "$InvisibleDir/"${ID[0]}"_"${PERCENTAGE[0]}"%.bam" | grep "@RG" | awk 'BEGIN {FS="\t"} { print $0 }'| samtools view -H "$InvisibleDir/"${ID[0]}"_"${PERCENTAGE[0]}"%.bam" | grep "@RG" | awk 'BEGIN {FS="\t"} { print $0 }' | tr $'\t' $'\n' | grep -n "SM:")
	SMPos=${SMPositionField:0:1}
	Stop=$(perl -E "say $SMPos-1")
	Subs=$(perl -E "say $SMPos+1")

	if [[ "$SMPos" = "$NumberOfFileds" ]]; then
		{ samtools view -H "$InvisibleDir/"${ID[0]}"_"${PERCENTAGE[0]}"%.bam" | grep "@RG" | awk -v f=1 -v t=$Stop '{for(i=f;i<=t;i++) printf("%s\t%s",$i,(i==t)?"\n":OFS)}'; echo "SM:"$IDVec"_"$PercVec; } | tr "\n" " " > "$InvisibleDir/.rg.txt"
	else
		{ samtools view -H "$InvisibleDir/"${ID[0]}"_"${PERCENTAGE[0]}"%.bam" | grep "@RG" | awk -v f=1 -v t=$Stop '{for(i=f;i<=t;i++) printf("%s\t%s",$i,(i==t)?"\n":OFS)}'; { echo "SM:"$IDVec"_"$PercVec; samtools view -H "$InvisibleDir/"${ID[0]}"_"${PERCENTAGE[0]}"%.bam" | grep "@RG" | awk 'BEGIN {FS="\t"} {print $Subs}'; } | tr "\n" "\t"; } | tr "\n" " " > "$InvisibleDir/.rg.txt"
	fi

	### merge new files ###

	cd $InvisibleDir

	MergeFiles=$(ls)
	samtools merge -rh ".rg.txt" $IDVec"_"$PercVec".bam" $MergeFiles
	samtools index $IDVec"_"$PercVec".bam"

	wait

	WorkingFolderForR=$WorkDir"/"$InvisibleDir

	R --slave --args $IDVec,${VARIANTS},$WorkingFolderForR < $FilesFolder/Scripts/Xome_VCF_generator.R 2>/dev/null

	cd $WorkDir

	### remove single files ###

	for z in $(seq 0 $End)
	do
		rm "$InvisibleDir/"${ID[$z]}"_"${PERCENTAGE[$z]}"%.bam"
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
	fi
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
	-i|--id)
	IDs="$1"
	shift
	;;
	-p|--percentages)
	PERCENTAGES="$1"
	shift
	;;
	-c|--coverages)
	COVERAGES="$1"
	shift
	;;
	-tc|--totalcoverage)
	TOTALCOVERAGE="$1"
	shift
	;;
	-v|--variants)
	VARIANTS="$1"
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

FilesFolder=`( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )`

if [[ -z ${LIST} ]] ; then
	COLUMNS=$(tput cols)
	printf '\e[37m%*s\n\e[m' $((($COLUMNS/3)+2)) "##########################"; printf '\e[37m%*s\n\e[m' $((($COLUMNS/3)+2)) " # Multiple Contamination #"; printf '\e[37m%*s\n\e[m' $((($COLUMNS/3)+2)) " ##########################"
	echo -e '\e[37m ___________\n|           |\n| User Mode |\n|___________|\n\e[0m'
	MultipleMixing # Start analysis
else
	COLUMNS=$(tput cols)
	printf '\e[37m%*s\n\e[m' $((($COLUMNS/3)+2)) "##########################"; printf '\e[37m%*s\n\e[m' $((($COLUMNS/3)+2)) " # Multiple Contamination #"; printf '\e[37m%*s\n\e[m' $((($COLUMNS/3)+2)) " ##########################"
	echo -e '\e[37m ________________\n|                |\n| Automated Mode |\n|________________|\n\e[0m'

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
		FILES=`echo $e | gawk '{print $1}'`
		IDs=`echo $e | gawk '{print $2}'`
		COVERAGES=`echo $e | gawk '{print $3}'`
		PERCENTAGES=`echo $e | gawk '{print $4}'`
		TOTALCOVERAGE=`echo $e | gawk '{print $5}'`
		MultipleMixing # Start analysis
	done
fi

### remove traces ###	

if [ `find $InvisibleDir -prune -empty` ]
then
	rm -rf $InvisibleDir
fi
echo -e '\e[31m   Done!\e[0m'
