# Script to look for variants in aligned BAM file
#!/bin/bash
function error_message {
	echo ""
	echo "Script annotate variants in BCF files"
	echo ""
	echo "Mandatory parameters:"
	echo -e "\t-i\tSamples plan (BCFs separated by ':')"
	echo -e "\t-o\tOutput directory (must exists)"
	echo -e "\t-c\tChromosome name (chr([0-9]{1,2}|X|Y|All))"
	echo -e "\t-a\tAnnovar database directory"
	echo -e "\t-b\tbcftools binary directory"
	echo -e "\t-V\tBCFTOOLS version (0.1.19 or 0.2.0)"
	echo ""
	echo "Optional parameters:"
	echo -e "\t-r\tRun commands: yes (default: no)"
	echo -e "\t-h\tHG version (default: hg19)"
	echo -e "\t-C\tRemove temporary files (default: yes)"
	echo -e "\t-P\ttool folder path"
	exit 99
}

while getopts "i:o:n:c:a:b:V:h:C:r:P:" optionName; do
case "$optionName" in

i) BCFS="$OPTARG";;                             ##  e.g. $somePath/sample1.bcf:$somePath2/sample2.bcf
o) OUTDIR="$OPTARG";;				##  e.g. Variants
n) barcodeID="$OPTARG";;			##  barcodeID
c) CHR="$OPTARG";;				##  e.g. chr22
a) ANNODB="$OPTARG";;				##  Path to Annovar DB
b) bcfDir="$OPTARG";;				## bcftools binary directory
V) VER="$OPTARG";;				##  BCFTOOLS version (0.1.19 / 0.2.0)
h) HG="$OPTARG";;				##  e.g. hg19
C) CLEAN="$OPTARG";;				##  e.g. yes/no
r) RUN="$OPTARG";;                              ##  yes or no (use no to check that everything is OK, then use yes to process)
P) DNASEQ_VARIANTS="$OPTARG";;                  ##  VARIANT_ANNOTATION_V5_2014_07_04 tool folder

esac
done

#echo "ANNOVAR ${ANNODB}" #debug line
#echo "BCFS ${BCFS}" #debug line
#echo "OUTDIR ${OUTDIR}" #debug line
#echo "CHR ${CHR}" #debug line
#echo "DNASEQ_VARIANTS ${DNASEQ_VARIANTS}" #debug line
# Check mandatory parameters
if [[ -z "$BCFS" ]]; then echo "Missing input BCFs (-b)."; error_message; fi
if [[ -z "$barcodeID" ]]; then echo "Missing barcodeID (-n)"; error_message; fi
if [[ ! -e "$BCFS" ]]; then echo "Samples plan does not exist"; error_message;fi
if [[ -z "$CHR" ]]; then echo "Missing chromosome name (-c)."; error_message; fi
if [[ -z "$OUTDIR" ]]; then echo "Missing output directory (-o)."; error_message; fi
if [[ ! -e "$OUTDIR" ]]; then echo "Output directory does not exist."; error_message; fi
if [[ -z "$ANNODB" ]]; then echo "Missing Annovar database (-a)"; error_message; fi
if [[ -z "$bcfDir" ]]; then echo "Missing bcftools directory (-b)"; error_message; fi
if [[ -z "$VER" ]]; then echo "Missing BCFTOOLS version (-v)"; error_message; fi

# Init optional parameters to default values when arguments are missing
if [[ -z "$HG" ]]; then HG="hg19"; fi
if [[ -z "$CLEAN" ]]; then CLEAN="yes"; fi

CHR_NB=`echo ${CHR} | awk '{split($0,part,":"); print part[1]}'`
CHR_COORD=`echo ${CHR} | awk '{split($0,part,":"); print part[2]}'`
sample_nb=`echo ${barcodeID} | awk '{split($0,part,"_"); print part[2]}'`

query="qsub $qsubMail \
	-N X${sample_nb}_${CHR_NB}_${CHR_COORD} \
	-e `hostname`:$OUTDIR/Annot_$CHR.ER -o `hostname`:$OUTDIR/Annot_$CHR.OU \
	-v workdir=$OUTDIR,annoDB=$ANNODB,hg=$HG,chr=$CHR,BCFs=$BCFS,BCFDIR=$bcfDir,clean=$CLEAN,ver=$VER,dnaPath=$DNASEQ_VARIANTS \
	$DNASEQ_VARIANTS/allelicExpression.V2.sh"
	# echo -e "\t"$query
	if [[ $RUN == "yes" ]]; then eval $query && sleep 1; fi
exit 0
