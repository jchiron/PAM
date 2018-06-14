
#!/bin/bash
# This script CONVERTs a VCF file generated in IonXpress format TO a VCF file compatible with mpileup BCF format

# Print time
exec_time() {
	echo "[$(date +'%T')] $1"
}

error_message() {    # $1 error msg $2 error code
        echo ""
        echo "#########################################"
        echo "# SCRIPT to convert IonXpress VCF files #"
	echo "#  into samtools compatible BCF files   #"
        echo "#########################################"
        echo ""
        echo "ERROR ....."
        echo $1
        echo ""
        echo "Mandatory parameters:"
        echo "\t-i\tFull path IonXpress VCF file name (path begins with /, file name without .vcf suffix)"
        echo "\t-o\tFull path BCF output file         (path begins with /, file name without .bcf suffix," 
        echo "\t\t                                    file name can be the same as VCF file name,"
        echo "\t\t                                    if path does not exist it is created)"
        echo "\t-U\tFull path UCSC reference genome   (path begins with /, file name without .fa  suffix)"
        echo "\t-C\tConfirm                           (yes / no)"
	echo ""
        echo "Optional parameters:"
        echo "\t-1\tTorrentServer VCF code for the STRAND FORWARD REF ALLELE COUNT  in  INFO field (default: SRF)"
        echo "\t-2\tTorrentServer VCF code for the STRAND REVERSE REF ALLELE COUNT  in  INFO field (default: SRR)"
        echo "\t-3\tTorrentServer VCF code for the STRAND FORWARD ALT ALLELE COUNT  in  INFO field (default: SAF)"
        echo "\t-4\tTorrentServer VCF code for the STRAND REVERSE ALT ALLELE COUNT  in  INFO field (default: SAR)"
        echo ""
        if [ ! -z "$2" ]; then exit $2; fi
}

while getopts "i:o:U:C:1:2:3:4:" optionName; do
case "$optionName" in

i) input="$OPTARG";;                              ## Full path IonXpress VCF file name (without .vcf)
o) output="$OPTARG";;                             ## Full path BCF output file name (without .bcf)
U) UCSCgenome="$OPTARG";;                         ## UCSC reference genome (w/o .fa)
C) CONF="$OPTARG";;                               ## Confirm yes / no
1) SRF="$OPTARG";;                                ## code for the STRAND FORWARD REF ALLELE COUNT      INFO field (default: SRF)
2) SRR="$OPTARG";;                                ## code for the STRAND REVERSE REF ALLELE COUNT      INFO field (default: SRR)
3) SAF="$OPTARG";;                                ## code for the STRAND FORWARD ALT ALLELE COUNT      INFO field (default: SAF)
4) SAR="$OPTARG";;                                ## code for the STRAND REVERSE ALT ALLELE COUNT      INFO field (default: SAR)

esac
done

#################  VCF FILE ################
if [ -z "${input}" ]; then
        error_message "Missing full path VCF input file"  #99
fi
TEST_ABSOLUTE_PATH=`echo ${input} | awk '{if ($0 ~ "^/") {print "TRUE"} else {print "FALSE"}}'`
if [ "${TEST_ABSOLUTE_PATH}" == "FALSE" ]; then
 error_message "\nInput path of VCF file (${input}) is not an absolute path!" #99
fi
input_suffix=`echo "$input"|awk -F . '{print $NF}'`
if [ "${input_suffix}" == "vcf" ]; then
 error_message "\nDo not provide .vcf suffix to your VCF input file name!" 99
fi
if [ ! -e "${input}.vcf" ]; then
        error_message "Input VCF file ${input}.vcf does not exist!" 99
fi

##############  BCF output file ###################"
if [ -z "$output" ]; then
        error_message "Missing full path BCF output file" 98
fi
TEST_ABSOLUTE_PATH=`echo ${output} | awk '{if ($0 ~ "^/") {print "TRUE"} else {print "FALSE"}}'`
if [ "${TEST_ABSOLUTE_PATH}" == "FALSE" ]; then
 error_message "\nOutput BCF path (${output}) is not an absolute path!" 98
fi
output_suffix=`echo "$output"|awk -F . '{print $NF}'`
if [ "${output_suffix}" == "bcf" ]; then
 error_message "\nDo not provide .bcf suffix to your BCF output file name!" 98
fi
outdir=`dirname ${output}`
if [ ! -e "${outdir}" ]; then mkdir -p $outdir; fi     ##  create directory for BCF files


###########   UCSC genome #############
if [ -z "${UCSCgenome}" ]; then
        error_message "Missing UCSC genome"  97
fi
TEST_ABSOLUTE_PATH=`echo ${UCSCgenome} | awk '{if ($0 ~ "^/") {print "TRUE"} else {print "FALSE"}}'`
if [ "${TEST_ABSOLUTE_PATH}" == "FALSE" ]; then
 error_message "\nUCSC path (${UCSCgenome}) is not an absolute path!" 97
fi
if [ ! -e "${UCSCgenome}.fa" ]; then
        error_message "File UCSC genome ${UCSCgenome} does nor exist!" 97
fi

############   CONFIRMATION to LAUNCH
if [[ "${CONF}" != "yes" ]] && [[ "${CONF}" != "no" ]]; then
 error_message "Confirmation (-C) not provided or wrong (${CONF})!" #96
fi


##  optional read depth fields
if [ -z "$SRF" ]; then SRF="SRF"; echo "WARNING! STRAND FORWARD REF ALLELE COUNT = $SRF"; fi    ## STRAND FORWARD REF ALLELE COUNT
if [ -z "$SRR" ]; then SRR="SRR"; echo "WARNING! STRAND REVERSE REF ALLELE COUNT = $SRR"; fi    ## STRAND REVERSE REF ALLELE COUNT
if [ -z "$SAF" ]; then SAF="SAF"; echo "WARNING! STRAND FORWARD ALT ALLELE COUNT = $SAF"; fi    ## STRAND FORWARD ALT ALLELE COUNT
if [ -z "$SAR" ]; then SAR="SAR"; echo "WARNING! STRAND REVERSE ALT ALLELE COUNT = $SAR"; fi    ## STRAND REVERSE ALT ALLELE COUNT

chromFile="${UCSCgenome}_Chromosomes.txt"
if [ ! -e "${chromFile}" ]; then
        error_message "UCSC genome chromosome index file ${chromFile} does not exist!" 95
fi

if [ "$CONF" == "yes" ]; then
#### write vcf files for TVC alleles called in alleles_IonXpress_nnn.xls file #########
exec_time "Write vcf line for all alleles found in alleles_IonXpress_nnn.xls."

#variantCaller sample output dir
SAMPLE_DIR=$(dirname ${input})

#sample IonXpress ID 
IonXpress=$(basename ${SAMPLE_DIR})

#sample alleles xls file
TSV_file=${SAMPLE_DIR}/alleles_${IonXpress}.xls

#temp vcf file
temp_vcf=${outdir}/${IonXpress}_temp.vcf
alleles_vcf=${outdir}/alleles_${IonXpress}.vcf

awk 'BEGIN {OFS="\t"};
$1 ~ "##" {print}' $workdir/${input}.vcf > ${outdir}/alleles_${IonXpress}.vcf

echo '##INFO=<ID=DP4,Number=4,Type=Integer,Description="# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">' >> ${outdir}/alleles_${IonXpress}.vcf

awk 'BEGIN {OFS="\t"};
$1 ~ "#CHROM" {print}' $workdir/${input}.vcf >> ${outdir}/alleles_${IonXpress}.vcf

#awk 'BEGIN {OFS="\t"};
#$1 ~ "#" {print}
#$1 !~ "#"{
#status=match($4,",",matchVal)
#if(status){i=0;for val in matchVal{print $1,$2,$3,matchVal[$i]}}}' $workdir/${input}.vcf > ${outdir}/alleles_${IonXpress}.vcf


#awk 'BEGIN {OFS="\t"};$1 !~ "Chrom"{DPF=$21;DPR=$23;SRF=DPF-$26;SRR=DPR-$27;status=match($12,"---",matchVal);if(status){known="."}else{known=matchVal[0]};print $1,$2,known,$3,$4,$8,"PASS","AF="$7",SRF="SRF",SRR="SRR",SAF="$26",SAR="$27,"AF:SRF:SRR:SAF:SAR",$7":"SRF":"SRR":"$26":"$27}' ${TSV_file} >> ${outdir}/alleles_${IonXpress}.vcf

awk 'BEGIN {OFS="\t"};$1 !~ "Chrom"{DPF=$23;DPR=$25;SRF=DPF-$28;SRR=DPR-$29;
#verifify deletion/insertion retreive vcf pos, ref and alt
status=match($5,/(Heterozygous|Homozygous)/,matchVal);if(status){
status=match($4,"-",matchVal);if(status){
POS=$2-1;
if(POS==$17){
ANCHOR=substr($18,1,1);
REF=ANCHOR$3;
ALT=ANCHOR;}
else{
dist=$2-$17;
ANCHOR=substr($18,dist,1);
REF=ANCHOR$3;
ALT=ANCHOR;}
}
else{
status=match($3,"-",matchVal);if(status){
POS=$2-1;
if(POS==$17){
ANCHOR=substr($19,1,1);
REF=ANCHOR;
ALT=ANCHOR$4;}
else{
dist=$2-$17;
ANCHOR=substr($18,dist,1);
REF=ANCHOR;
ALT=ANCHOR$4;}
}
else{
POS=$2;REF=$3;ALT=$4
}
};
status=match($13,"---",matchVal);if(status){known="."}else{known=$13};status=match($7,"100",matchVal);if(status){gt="1/1"}else{gt="0/1"};print $1,POS,known,REF,ALT,$8,"PASS","AF="$7",SRF="SRF",SRR="SRR",SAF="$28",SAR="$29,"GT:AF:SRF:SRR:SAF:SAR",gt":"$7":"SRF":"SRR":"$28":"$29}}' ${TSV_file} >> ${outdir}/alleles_${IonXpress}.vcf



exec_time "Convert IonXpress VCF file."
# For each chromosome:
for i in {1..22} X Y; do
#	- get INFO field 
#	- split into sub-fields; look for SAF, SAR, SRF, SRR fields  
#	- add DP4 field at the end of INFO field   DP4=SRF,SRR,SAF,SAR
#	- convert VCF to BCF
#	- index BCF
	awk -v i=$i -v SAF=$SAF -v SAR=$SAR -v SRF=$SRF -v SRR=$SRR 'BEGIN {OFS="\t";err=0;DP4="";DP4SAF="";DP4SAR="";DP4SRF="";DP4SRR=""}; 
               
                # Print header (start with #)
                $1 ~ "#" {print}
                $1 !~"#" && $1 == "chr"i { 
		  ### In INFO column ($8), get SAF, SAR, SRF, SRR using regex and add DP4 tag and value
                  # SAF
		  status=match($8,"SAF=[0-9]+",matchVal);
		  if (status) { split(matchVal[0],TAG,"="); DP4SAF=TAG[2] } else {err=1}
                  # SAR
                  status=match($8,"SAR=[0-9]+",matchVal);
                  if (status) { split(matchVal[0],TAG,"="); DP4SAR=TAG[2] } else {err=1}
                  # SRF
                  status=match($8,"SRF=[0-9]+",matchVal);
                  if (status) { split(matchVal[0],TAG,"="); DP4SRF=TAG[2] } else {err=1}
                  # SRR
                  status=match($8,"SRR=[0-9]+",matchVal);
                  if (status) { split(matchVal[0],TAG,"="); DP4SRR=TAG[2] } else {err=1}
                  #
                  if (!err) {
                    DP4="DP4="DP4SRF","DP4SRR","DP4SAF","DP4SAR;
		    #  Add DP4 at the end of INFO
                    $8=$8";"DP4;
		    # If alt base is X, set it as ref. X means ref when quality
		    if ($5=="X") {$5=".";}
		    print;
                  } 
                 err=0;DP4="";DP4SAF="";DP4SAR="";DP4SRF="";DP4SRR=""
	        }' ${outdir}/alleles_${IonXpress}.vcf | bcftools sort - > ${outdir}/alleles_${IonXpress}.chr${i}.vcf
	bgzip -c ${outdir}/alleles_${IonXpress}.chr${i}.vcf > ${outdir}/alleles_${IonXpress}.chr${i}.vcf.gz 
	bcftools index -t ${outdir}/alleles_${IonXpress}.chr${i}.vcf.gz 
	bcftools view -O b -r chr${i} ${outdir}/alleles_${IonXpress}.chr${i}.vcf.gz > $workdir/${output}_chr${i}_mpileup.bcf 
	bcftools index $workdir/${output}_chr${i}_mpileup.bcf
done
else
 error_message "Confirmation was ${CONF} - process not launched!" 94
fi

exec_time "End." && exit 0

#$workdir/${input}.vcf | \
#$outdir/alleles_${IonXpress}.vcf | \ 
