#!/bin/bash
############
# Routines #
############
function exec_time {
        # Print time and message given as arg
         echo "[$(date +'%T')] $1"
}

function error_message { # $1 error msg  $2 error code
	echo ""
	echo "########################################################"
	echo "Script to concatenate variant tables of each chromosomes"
	echo "########################################################"
	echo ""
	echo "ERROR ....."
	echo -e $1
	echo ""
	echo "Mandatory parameters:"
	echo -e "\t-V\tFull path of Variant directory  (where chromosome variant tables are stored - must begin with /)"
	echo -e "\t-A\tAggregate variants     (yes or no)"
	echo -e "\t-S\tSplit by chromosome    (yes or no)"
	echo -e "\t-Z\tCompress .gz           (yes or no)"
	echo -e "\t-C\tConfirmation to launch (yes or no)"
	echo ""
        if [ ! -z "$2" ]; then exit $2; fi
}

while getopts "V:A:S:Z:C:" optionName; do
case "$optionName" in

V) workdir="$OPTARG";;				##  samples directory
A) aggregate="$OPTARG";;				##  confirmation to aggregate
S) split="$OPTARG";;				##  confirmation to split
Z) compress="$OPTARG";;				##  confirmation to compress
C) conf="$OPTARG";;				##  confirmation to launch

esac
done

# Check mandatory parameters
if [[ -z "${workdir}" ]]; then error_message "Missing Variant dir (-V)" 99; fi
TEST_ABSOLUTE_PATH=`echo ${workdir} | awk '{if ($0 ~ "^/") {print "TRUE"} else {print "FALSE"}}'`
if [[ "${TEST_ABSOLUTE_PATH}" == "FALSE" ]]; then
 error_message "\nVariant directory path (${workdir}) is not an absolute path!" 99
fi
if [[ ! -e "${workdir}" ]]; then
 error_message "\nVariant directory path (${workdir}) does not exist!" 99
fi

if [[ "${aggregate}" != "yes" && "${aggregate}" != "no" ]]; then
 error_message "\nConfirmation to aggregate is absent or incorrect (${aggregate})!" 98
fi

if [[ "${split}" != "yes" && "${split}" != "no" ]]; then
 error_message "\nConfirmation to split is absent or incorrect (${split})!" 98
fi

if [[ "${compress}" != "yes" && "${compress}" != "no" ]]; then
 error_message "\nConfirmation to compress is absent or incorrect (${compress})!" 97
fi

if [[ "${conf}" != "yes" && "${conf}" != "no" ]]; then
 error_message "\nConfirmation is absent or incorrect (${conf})!" 96
fi


if [ "$conf" == "yes" ]; then

    ###  get list of Union sample files for which the Annotation has been created
    exec_time "Get list of variant files ..."
    clf=""    ##  list of Union_sample chromosome files
    for f in $(ls -d ${workdir}/*/); do echo ${f%%/}; done | xargs -n 1 basename | \
    awk '{split($0,rec,":"); split(rec[2],pos,"-"); 
          chr=substr(rec[1],1,3); c=substr(rec[1],4); 
          if (c=="X") {c=23}; if (c=="Y") {c=24}; print chr, c, pos[1], pos[2]}' | \
    sort -k2n -k3n | awk '{if ($2==23) {$2="X"}; if ($2==24) {$2="Y"}; print $1$2":"$3"-"$4 }' > ${workdir}/part_file.txt
    while read l ; do
      f=${workdir}/${l}/Union_samples_chr*.txt
      #echo $f
      if [ -e $f ]; then
        clf=${clf}" "${f}
      fi
    done < ${workdir}/part_file.txt
    exec_time "List of variant files created ..."
#  echo -e "\nList of folders to aggregate:\n${clf}"

    ###  get header from CHROMOSOME 1
    exec_time "Get header of each part files ..."
    header=${workdir}/header.txt
    head -n1 `echo  ${clf} | cut -d " " -f 1`  > ${header}
    exec_time "Header created ..."

  if [ "$aggregate" == "yes" ]; then
    exec_time "Start aggregating part files ..."

    ###  consolidate only existing Union-files
    cat ${header} > ${workdir}/Union_samples_chrAll.txt
    prev_chr=""
    for chr_file in ${clf}; do
      ### split file by chromosomes
      if [ "$split" == "yes" ]; then
          ## get chr of chr_file
          cur_chr=`basename $chr_file | cut -d ":" -f 1`
          if [ "${prev_chr}" != "${cur_chr}" ]; then 
             cat ${header} > ${workdir}/Union_samples_chr${cur_chr}.txt
             prev_chr=${cur_chr}
             exec_time "start saving ${cur_chr} ..."
          fi 
          grep -Ev '"CHROM"[[:space:]]"START"[[:space:]]"STOP"' ${chr_file} >> ${workdir}/Union_samples_chr${cur_chr}.txt
      fi
          grep -Ev '"CHROM"[[:space:]]"START"[[:space:]]"STOP"' ${chr_file} >> ${workdir}/Union_samples_chrAll.txt
          exec_time "file ${chr_file} aggregated ..."
    done
  fi 
  if [ "$split" == "yes" ]; then
    rm  ${workdir}/Union_samples_chrAll.txt
    exec_time "Union_samples_chrAll.txt removed"
  fi

  ### compress the complete file if it still exists
  if [[ "$compress" == "yes" && -e ${workdir}/Union_samples_chrAll.txt ]]; then
    gzip ${workdir}/Union_samples_chrAll.txt
    exec_time "Union_samples_chrAll.txt compressed"
  fi

else
  error_message "\nConfirmation is ${conf} - cannot launch aggregation!" 95
fi

exec_time "Done!"

exit 0

