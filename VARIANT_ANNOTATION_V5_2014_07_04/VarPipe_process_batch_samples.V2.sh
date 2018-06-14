#!/bin/bash
error_message() { # $1 error msg $2 error code
	echo ""
	echo "########################################"
	echo "Script annotate variants in BCF files"
	echo "########################################"
	echo ""
	echo "ERROR ....."
	echo -e $1
	echo ""
	echo "Mandatory parameters:"
        echo -e "\t-S\tSamplesplan                  (each line contains a sample to annotate and to merge with the others,"
        echo -e "\t  \t                              the first col is the sample name, only the first column is mandatory)"
	echo -e "\t-N\tSample Name                  (usually starting with IonXpress)"
 	echo -e "\t-I\tInput directory              (where sample folders are stored)"
	echo -e "\t-O\tOutput directory             (where annotated variants are stored - if it does not exist, then it is created)"
	echo -e "\t-F\tbcftools directory           (binary directory)"
	echo -e "\t-B\tUse BCF_cleaned or BCF dir   (values: BCF or BCF_cleaned)"
	echo -e "\t-D\tSize of genomic interval     (used to split each BCF file in equal sized intervals - each interval is annotated in parallel"
	echo -e "\t  \t                              integer >0  if no split is desidered than set to 999999999"
        echo -e "\t  \t                              Proposed size: for a trio Constit + Tumor DNA and T RNA Whole Genome --> size=10000000"
	echo -e "\t  \t                              This is roughly equivalent to 3 300 000 000 / 10 000 000 = 330 jobs that are run in parallel)"
	echo -e "\t-V\tBcftools version             (values: 0.1.19 / 0.2.0)"
	echo -e "\t-C\tConfirm launch               (values: yes / no)"
	echo -e "\t-T\tChromInfo table file"
	echo -e "\t-P\tool folder path"
	echo -e "\t-A\tANNOVAR db folder path"
	echo ""
        if [ ! -z "$2" ]; then exit $2; fi
}

while getopts "S:N:I:O:F:B:D:V:C:T:P:A:" optionName; do
case "$optionName" in

S) samplesplan="$OPTARG";;                      ##  path to sample plan
N) sampleName="$OPTARG";;			##  sample name
I) inpDIR="$OPTARG";;				##  samples directory
O) outDIR="$OPTARG";;				##  output directory
F) bcfDir="$OPTARG";;				##  bcftools directory
B) BCF="$OPTARG";;				##  BCF_cleaned or BCF
D) SIZE="$OPTARG";;				##  Size of interval to annotate
V) VER="$OPTARG";;				##  Bcftools version
C) conf="$OPTARG";;                             ##  confirm option
T) chrSizeTab="$OPTARG";;                       ##  UCSC chromosomes info table file
P) DNASEQ_VARIANTS="$OPTARG";;                  ##  VARIANT_ANNOTATION_V5_2014_07_04 tool folder
A) ANNODB="$OPTARG";;                        ##  ANNOVAR DB path 
esac
done

# test samplesplan
if [[ -z "$samplesplan" ]]; then error_message "Missing Samplesplan (-S)" 97; fi
TEST_ABSOLUTE_PATH=`echo ${samplesplan} | awk '{if ($0 ~ "^/") {print "TRUE"} else {print "FALSE"}}'`
if [[ "${TEST_ABSOLUTE_PATH}" == "FALSE" ]]; then
 error_message "\nSampleplan path (${samplesplan}) is not an absolute path!" 97
fi
if [[ ! -e "$samplesplan" ]]; then error_message "Samplesplan $samplesplan does not exist" 97; fi

# sample name exists ?
if [[ -z "$sampleName" ]]; then
 error_message "\nsampleName not provided (${sampleName})!" 92
fi

# test input dir
if [[ -z "${inpDIR}" ]]; then error_message "Missing input dir (-I)" 99; fi
TEST_ABSOLUTE_PATH=`echo ${inpDIR} | awk '{if ($0 ~ "^/") {print "TRUE"} else {print "FALSE"}}'`
if [[ "${TEST_ABSOLUTE_PATH}" == "FALSE" ]]; then
 error_message "\nInput directory path (${inpDIR}) is not an absolute path!" 99
fi
if [[ ! -e "${inpDIR}" ]]; then error_message "Input directory (${inpDIR}) does not exist." 99; fi

# test out dir
if [[ -z "$outDIR" ]]; then error_message "Missing output directory (-O)" 98; fi
TEST_ABSOLUTE_PATH=`echo ${outDIR} | awk '{if ($0 ~ "^/") {print "TRUE"} else {print "FALSE"}}'`
if [[ "${TEST_ABSOLUTE_PATH}" == "FALSE" ]]; then
 error_message "\nOutput directory path (${outDIR}) is not an absolute path!" 98
fi

# test use of deeply removed BCF or not (BCF)
if [[ -z "$BCF" ]]; then error_message "Missing indicator of BCF_cleaned or BCF directory" 96; fi
if [[ "$BCF" != "BCF_cleaned" && "$BCF" != "BCF" ]]; then error_message "BCF dir indicator not correct ($BCF)" 96; fi

# Size of genomic interval to annotate
if [[ -z "$SIZE" ]]; then
 error_message "\nSize not provided (${SIZE})!" 95
fi

# Version of bcftools
if [[ "${VER}" != "0.1.19" && "${VER}" != "0.2.0" ]]; then
 error_message "\nBCFTOOL version not correct or not provided (${VER})!" 94
fi

# Confirm not provided
if [[ "$conf" != "yes" && "$conf" != "no" ]]; then
 error_message "\nConfirmation not provided or wrong (${conf})!" 93
fi

###  CHROMOSOME SIZE TABLE
#chrSizeTab=/home/carlucchesi/BRIO/annotation/UCSC_genome_Hg19_Spikein_2013_08_22/UCSC_hg19_chr1_22_X_Y_chromInfo_table.txt
#chrSizeTab=/home/carlucchesi/BRIO/annotation/UCSC_genome_Hg19_Spikein_2013_08_22/UCSC_hg19_chrY_chromInfo_table.txt

mkdir -p ${outDIR}
nsamples=`awk 'NR>1 && $1 !~ "#" {s++} END {print s}' ${samplesplan}`

BCF_cleaned_string() {   # $1 the file of sample names 
 ## returns the chain of BCF full path names (separated by :) for each chromosome
	i=0
	bcf=""
	for ARG in $*; do
		if [[ $i -ne 0 ]]; then bcf="${bcf}:"; fi                              # but the first time 
		bcf="${bcf}${inpDIR}/${ARG}/BCF_cleaned/${ARG}_${c}_mpileup.bcf"  #  all times
		i=1
	done
	echo ${bcf}
}

BCF_string() {    # $1 the file of sample names
 ## returns the chain of BCF full path names (separated by :) for each chromosome
        i=0
        bcf=""
        for ARG in $*; do
                if [[ $i -ne 0 ]]; then bcf="${bcf}:"; fi                           #  but the first time
                bcf="${bcf}${inpDIR}/${ARG}/BCF/${ARG}_${c}_mpileup.bcf"       #  all times
                i=1
        done
        echo ${bcf}
}


# process BCF by sample
#echo $inpDIR
nbJobsLaunched=0
if [ "$conf" == "yes" ]; then
  echo "Perform variants annotation by chr:region for the UNION of the ${nsamples} samples in samplesplan"
  while read line; do
    c=`echo $line | cut -f1 -d " "`
    cS=`echo $line | cut -f2 -d " "`
    # echo -e "\n$c size:$cS"
    export c=${c}
    start=1
    let "stop=SIZE - 1"
    limitReached=0
    while [ $limitReached  -eq 0 ]; do
          if [ $stop -gt $cS ]; then 
            stop=$cS; 
            limitReached=1; 
          fi
          # echo ${c}:${start}-${stop}
          chr=${c}:${start}-${stop}
          let "start=stop + 1"
          let "stop=start + SIZE - 1"
          ##  create sample list from the first column of the samplesplan
          samples=`awk 'NR>1 && $1 !~ "#" {print $1}' ${samplesplan}`
          # echo ${samples}
          if [ "$BCF" == "BCF_cleaned" ]; then bcfs=`BCF_cleaned_string ${samples}`; fi
          if [ "$BCF" == "BCF" ]; then bcfs=`BCF_string ${samples}`; fi
#        if [[ "TRUE" == "FALSE" ]]; then
          echo -e "\n ... launching ${chr}"
          mkdir -p ${outDIR}/$chr
          echo ${bcfs} > ${outDIR}/${chr}/samplesPlan.txt
		  # echo "${DNASEQ_VARIANTS}/allelicExpression_launcher.V2.sh -a ${ANNODB} -c ${chr} -o ${outDIR}/${chr} -i ${outDIR}/${chr}/samplesPlan.txt -V ${VER} -r yes -P ${DNASEQ_VARIANTS}" #debug line
          bash ${DNASEQ_VARIANTS}/allelicExpression_launcher.V2.sh -a ${ANNODB} -n ${sampleName} -b ${bcfDir} -c ${chr} -o ${outDIR}/${chr} -i ${outDIR}/${chr}/samplesPlan.txt -V ${VER} -r yes -P ${DNASEQ_VARIANTS}
#        fi              
    let "nbJobsLaunched = nbJobsLaunched + 1"
    done
  done < $chrSizeTab
else
  error_message "\nConfirmation is NO - Analysis not launched!" 94
fi

echo -e "\n\nNb of jobs launched: $nbJobsLaunched"

exit 0

