#!/bin/bash

#################################
#                               #
#  ##  ###    ##  ##  ##    ##  #
#  ##  ####   ##  ##   ##  ##   #
#  ##  ## ##  ##  ##    ####    #
#  ##  ##  ## ##  ##    ####    #
#  ##  ##   ####  ##   ##  ##   #
#  ##  ##    ###  ##  ##    ##  #
#                               #
#################################

######################
#                    #
#  Define functions  #
#                    #
######################

function display_usage()
{
	COLUMNS=$(tput cols)
	printf "%*s\n" $(($COLUMNS/2)) "~~~~~~~~~~~~~~~~~~";printf "%*s\n" $(($COLUMNS/2)) "|                |";printf "%*s\n" $(($COLUMNS/2)) "|   INIXALIZER   |";printf "%*s\n" $(($COLUMNS/2)) "|                |";printf "%*s\n" $(($COLUMNS/2)) "~~~~~~~~~~~~~~~~~~";printf "%*s\n" $((($COLUMNS/2)+29)) "____________________________________________________________________________";printf "%*s\n" $((($COLUMNS/2)+29)) "InXalizer generates the data useful for Xome-Blender running.       "; printf "%*s\n" $((($COLUMNS/2)+29)) "----------------------------------------------------------------------------"; printf "%*s\n" $((($COLUMNS/2)+26)) "Usage: $0 [-f file1,file2] [-i id1,id2] -l my_label -r ref.fa"; printf "%*s\n" $((($COLUMNS/2)+29)) "____________________________________________________________________________"; printf '\e[1;39m%*s\n\e[m' $((($COLUMNS/2)-4)) "Arguments:";printf "%*s" $((($COLUMNS/3)+3)) "-f,--files"; printf "%*s\n" $((($COLUMNS/3)-13)) "Input files comma separated";printf "%*s" $((($COLUMNS/3))) "-i,--id"; printf "%*s\n" $((($COLUMNS/3)-10)) "Samples ID's; same order of files";printf "%*s" $((($COLUMNS/3)+3)) "-l,--label"; printf "%*s\n" $((($COLUMNS/3)-13)) "The label for the variants file";printf "%*s" $((($COLUMNS/3)+7)) "-r,--reference"; printf "%*s\n" $((($COLUMNS/3)-17)) "The reference file";printf "%*s" $((($COLUMNS/3)+1)) "-b,--bed"; printf "%*s\n" $((($COLUMNS/3)-11)) "The target file in .bed format (not required)";printf "%*s" $((($COLUMNS/3)+5)) "-t,--threads"; printf "%*s\n" $((($COLUMNS/3)-15)) "The number of threads used for this analysis";printf "%*s" $((($COLUMNS/3)+4)) "-o,--output"; printf "%*s\n" $((($COLUMNS/3)-14)) "Output directory. If omitted, generates a";printf '\e[8;39m%*s\e[m' $((($COLUMNS/3)+4)) "-o,--output"; printf "%*s\n" $((($COLUMNS/3)-14)) "results directory in the current position"; printf "%*s\n" $((($COLUMNS/2)+29)) "____________________________________________________________________________"; printf "%*s\n" $((($COLUMNS/2)+20)) "NB: The order of samples must be respected in each option!"; printf "%*s\n" $((($COLUMNS/2)+29)) "----------------------------------------------------------------------------"; printf "%*s\n" $((($COLUMNS/2)+29)) "Xome-Blender. Written by Roberto Semeraro, Department of Clinical and Speri-"; printf "%*s\n" $((($COLUMNS/2)+29)) "mental Medicine, University of Florence. For bug reports or suggestion write"; printf "%*s\n" $((($COLUMNS/2)-19)) "to roberto.semeraro@unifi.it"
}

gold_gen_wg()
{
	SamActiveProcess=$(ps aux | grep "samtools depth ${SAMBED}" | awk '{print $2}' | wc -l)
	GATKActiveProcess=`ps -ef | grep "HaplotypeCaller -R ${REFERENCE}" |  sed 's/^.*nct/nct/' | awk '{print $2}' | head -n -1 | awk '{sum+=$1};END{print sum}'`
	if [[ -z $GATKActiveProcess ]]; then
		GATKActiveProcess=0
	fi
	REMCPUs=$((${THREADS}-$SamActiveProcess-$GATKActiveProcess-1))
	if [[ $Samp -le $REMCPUs ]]; then
		GATKCpus=$REMCPUs
		paste <(date +"%d/%m/%y %T") <(echo "GATK" "SAMPLE="${ID[$f]} "GATKCPus="$GATKCpus "REMCPUs="$REMCPUs "Samp="$Samp "COVCPUs="$COVCPUs "HCFlag="$HCFlag "- start") -d " " >> $WorkingFolder"/"${LABEL}.log
	fi	
	GenomeAnalysisTK -T HaplotypeCaller \
	-R ${REFERENCE} \
	-I ${FILE[$f]} \
	-o ${ID[$f]}".g.vcf" \
	-rf BadCigar \
	-nct $GATKCpus \
	--max_alternate_alleles 1 \
	--doNotRunPhysicalPhasing \
	--emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 &>/dev/null
	paste <(date +"%d/%m/%y %T") <(echo "GATK" "SAMPLE="${ID[$f]} "GATKCPus="$GATKCpus "REMCPUs="$REMCPUs "Samp="$Samp "COVCPUs="$COVCPUs "HCFlag="$HCFlag "- finish") -d " " >> $WorkingFolder"/"${LABEL}.log
}

gold_gen_wes()
{  
	SamActiveProcess=$(ps aux | grep "samtools depth ${SAMBED}" | awk '{print $2}' | wc -l)
	GATKActiveProcess=`ps -ef | grep "HaplotypeCaller -R ${REFERENCE}" |  sed 's/^.*nct/nct/' | awk '{print $2}' | head -n -1 | awk '{sum+=$1};END{print sum}'`
	if [[ -z $GATKActiveProcess ]]; then
		GATKActiveProcess=0
	fi
	REMCPUs=$((${THREADS}-$SamActiveProcess-$GATKActiveProcess-1))
	if [[ $Samp -le $REMCPUs ]]; then
		GATKCpus=$REMCPUs
		paste <(date +"%d/%m/%y %T") <(echo "GATK" "SAMPLE="${ID[$f]} "GATKCPus="$GATKCpus "REMCPUs="$REMCPUs "Samp="$Samp "COVCPUs="$COVCPUs "HCFlag="$HCFlag "- start") -d " " >> $WorkingFolder"/"${LABEL}.log
	fi
	GenomeAnalysisTK -T HaplotypeCaller \
	-R ${REFERENCE} \
	-I ${FILE[$f]} \
	-L ${BED} \
	-o ${ID[$f]}".g.vcf" \
	-rf BadCigar \
	-nct $GATKCpus \
	--max_alternate_alleles 1 \
	--doNotRunPhysicalPhasing \
	-stand_call_conf 30 \
	-stand_emit_conf 10 \
	--emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 &>/dev/null
	paste <(date +"%d/%m/%y %T") <(echo "GATK" "SAMPLE="${ID[$f]} "GATKCPus="$GATKCpus "REMCPUs="$REMCPUs "Samp="$Samp "COVCPUs="$COVCPUs "HCFlag="$HCFlag "- finish") -d " " >> $WorkingFolder"/"${LABEL}.log		
}

vcf_gen()
{
	AllSampBack=$(echo ${ID[@]})
	N_AllSampBack=$(echo $AllSampBack | wc -w)
	THNUM=$((${THREADS}-$HCFlag))
	SB_ID_VEC=$(printf ' --variant \n%.0s' `seq 1 $N_AllSampBack | awk 'BEGIN { ORS = " " } { print }'`)
	FORMATTED_ID=$(echo ${ID[@]} | awk '{OFS=RS;$1=$1}1')
	FORMAT_COL=$(printf '.g.vcf\n%.0s' `seq 1 $N_AllSampBack | awk 'BEGIN { ORS = " " } { print }'`)
	IDX_COL=$(printf '*\n%.0s' `seq 1 $N_AllSampBack | awk 'BEGIN { ORS = " " } { print }'`)
	VARIANTS=$(paste <(echo "$SB_ID_VEC") <(echo "$FORMATTED_ID") <(echo "$FORMAT_COL") -d '')
	VARIANTSFORRM=$(paste <(echo "$FORMATTED_ID") <(echo "$FORMAT_COL") <(echo "$IDX_COL") -d '' | awk 'BEGIN { ORS = " " } { print }')
	VARIANTSROW=$(echo "$VARIANTS" | awk 'BEGIN { ORS = " " } { print }')
	GATKLINE=$(echo "GenomeAnalysisTK -T GenotypeGVCFs -R ${REFERENCE} --max_alternate_alleles 1 -nt $THNUM -stand_emit_conf 10 -stand_call_conf 30 $VARIANTSROW -o ${LABEL}.vcf &>/dev/null; wait")

	eval "$GATKLINE" # Generates Variant File

	PathForR=$WorkingFolder"/"$InvisibleDir"/"${LABEL}".vcf"
	R --slave --args $PathForR < $FilesFolder/Scripts/GS_Cleaner.R 2>/dev/null
	wait
	if [[ -s ${LABEL}.vcf ]]; then
		cd $WorkingFolder
		mv $InvisibleDir"/"${LABEL}".vcf" "."
		rm $InvisibleDir"/"${LABEL}".vcf.idx"
		wait
		if [[ -e ${LABEL}.vcf ]]; then
			cd $InvisibleDir
			RMVar=$(echo rm $VARIANTSFORRM)
			$RMVar
		fi
	fi
}

cov_calc()
{
	val_extract()
	{
		samtools depth ${SAMBED} -r $c ${FILE[$f]} | tee >(printf "%.0f\n" $(perl -E "say ($(wc -l | awk '{print $1}')/2)") > $CovRun\.Flngt$c) > /dev/null  >(awk '{print $3}' | sort -n > $CovRun\.srtChr$c) ; while [[ ! -s $CovRun\.Flngt$c ]] ;do wait; done; declare Lg_$c=`cat $CovRun\.Flngt$c` ; while [[ `perl -E "say ($(wc -l $CovRun\.srtChr$c | awk '{print $1}'))"` -lt "$((($(eval echo \$Lg_$c)) * 2-1))" ]] ;do wait; done ;sed -n $(echo "$(eval echo \$Lg_$c)"p) $CovRun\.srtChr$c >> $CovRun\.Means; if grep -q $(sed -n $(echo "$(eval echo \$Lg_$c)"p) $CovRun\.srtChr$c) $CovRun\.Means; then rm $CovRun\.srtChr$c $CovRun\.Flngt$c; fi
		paste <(date +"%d/%m/%y %T") <(echo "SAMTOOLS" "SAMPLE="${ID[$f]} "CHR="$c "FLAG="$Flag "COVCPUs="$COVCPUs "HCCPUS="$HCCPUS "HCFlag="$HCFlag "SAMP="$Samp "- finish") -d " " >> $WorkingFolder"/"${LABEL}.log
	}
	Flag=0 # Tiene conto delle cpu usate nei processi di calcolo coverage
	for c in $Chromosomes
	do
		if [[ $COVCPUs -ne 6 && $HCFlag -ne ${HCCPUS} ]]; then
			COVCPUs=6
		fi
		Flag=$(($Flag+1))
		val_extract &
		paste <(date +"%d/%m/%y %T") <(echo "SAMTOOLS" "SAMPLE="${ID[$f]} "CHR="$c "FLAG="$Flag "COVCPUs="$COVCPUs "HCCPUS="$HCCPUS "HCFlag="$HCFlag "SAMP="$Samp "- start") -d " " >> $WorkingFolder"/"${LABEL}.log
		if [[ $Flag -eq $COVCPUs ]]; then
			wait             
			Flag=0
		fi
	done
	wait
	if [[ `perl -E "say ($(wc -l $CovRun\.Means | awk '{print $1}'))"` == 6 ]]; then paste <(awk '{ total += $1; count++ } END { printf "%.0f\n", total/count }' $CovRun\.Means) <(echo ${ID[$f]}) >> CoverageStats.txt; else echo "Something went wrong!" && exit; fi
	wait
	rm $CovRun\.Means
}

prepare()
{
	echo -e "   - "${ID[$f]}
	echo -e "     "'\e[37mCalculation of coverage and preparation for variants file...\e[0m' && cov_calc 
	$GS # Generate gVCF
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
	-f|--files)
	FILES="$1"
	shift
	;;
	-i|--id)
	IDs="$1"
	shift
	;;
	-l|--label)
	LABEL="$1"
	shift
	;;
	-r|--reference)
	REFERENCE="$1"
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

if [[ -z ${FILES} ]] ; then
	echo -e '\e[31mNo files selected!\e[0m' && exit
fi

if [[ -z ${IDs} ]] ; then
	echo -e '\e[31mNo IDs assigned!\e[0m' && exit
fi

if [[ ! -z ${BED} ]] ; then
	GS="gold_gen_wes"
	SAMBED="-b ${BED}"
else
	echo -e '\e[31mNo bed file selected! You using genomic data!\e[0m'
	GS="gold_gen_wg"
	SAMBED=""
fi

if [[ -z ${REFERENCE} ]] ; then
	echo -e '\e[31mNo reference selected!\e[0m' && exit
fi

if [[ -z ${LABEL} ]] ; then
	echo -e '\e[31mNo label defined!\e[0m' && exit
fi

if [[ -z ${THREADS} ]] ; then
	THREADS=1
elif [[ ${THREADS} =~ [[:alpha:]] || ${THREADS} -lt 0 ]] ; then
	echo "Number of threads it's not an integer!" && exit
fi

IFS=',' read -ra FILE <<< "$FILES"
IFS=',' read -ra ID <<< "$IDs"

if [[ ${#ID[@]} != ${#FILE[@]} ]]; then
	echo -e '\e[31mThe number of files and IDs is different, fix it!\e[0m' && exit
fi

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

######################
# Enviroment setting #
######################

FilesFolder=`( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )`

WorkingFolder=`pwd`

echo -e " - "'\e[31mSAMPLES\e[0m' = ${IDs//,/ }

End=$(echo $(perl -E "say ${#FILE[@]}-1"))

HCCPUS=1
if [[ ${THREADS} -gt 6 ]] ; then
	CPUGROUPS=$(echo $(perl -E "say (int(${THREADS}/6) )")) # Determina il numero di cpu a gruppi di 6 (dato che 6 cpu servono per il calcolo del coverage)
	REMAININGCPUs=$(echo $(perl -E "say ${THREADS}-(6*$CPUGROUPS)")) # Se il numero di cpu non è multiplo di 6 calcola quante cpu restano
	if [[ $REMAININGCPUs -ne 0 ]]; then
		HCCPUS=$(echo $(perl -E "say $CPUGROUPS+1")) # Se ci sono cpu che avanzano, aumenta il numero di cpu utilizzabili ma modifica il numero di quelle utilizzate per il calcolo del coverage
	else
		HCCPUS=$(echo $CPUGROUPS)
	fi
fi

InvisibleDir=.Run$RANDOM

if [[ ! -d $InvisibleDir ]]; then
	mkdir $InvisibleDir
fi

cd $InvisibleDir

GATKCpus=1
HCFlag=0 # Tiene conto dei processi primari in corsa
Samp=${#ID[@]}
COVCPUs=6

#########
# Start #
#########

for f in $(seq 0 $End)
do
	CovRun=.C$RANDOM
	HCFlag=$(($HCFlag+1))
	if [[ "$HCFlag" -eq 1 ]]; then
		StopFile=${FILE[$f]}
	fi

	if [[ "$HCFlag" -eq ${HCCPUS} && $REMAININGCPUs -ne 0 ]]; then
		COVCPUs=$(echo $REMAININGCPUs)
	fi

	prepare &
	Samp=$(($Samp-1))
	if [[ $HCFlag -eq ${HCCPUS} ]]; then
		Pido=$(ps aux | grep "samtools depth ${SAMBED} $PidChr ${FILE[0]}" | awk '{print $2}' | head -n 1)
		while [[ -z $Pido ]]; do sleep 0.1; done
		WaitStrng=$( echo "while [ -e /proc/$Pido ]; do sleep 0.1; done" )
		eval "$WaitStrng"
		HCFlag=0
	fi
done
wait
echo -e "     "'\e[37mDone\e[0m'
wait
echo -e "   - "'\e[30mVariants File Creation\e[0m' && vcf_gen
wait
echo -e "     "'\e[37mDone\e[0m'
cd $WorkingFolder

mv $InvisibleDir"/"CoverageStats.txt "."

if [[ ! -z ${OUTPUT} ]] ; then
	if [[ ! -d ${OUTPUT} ]]; then
		mkdir ${OUTPUT}
	fi
	mv CoverageStats.txt ${OUTPUT}
	mv ${LABEL}".vcf" ${OUTPUT}
fi
rmdir $InvisibleDir
