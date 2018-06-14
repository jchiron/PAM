#!/bin/bash
############
# Routines #
############
function exec_time {
        # Print time and message given as arg
         echo -e "[$(date +'%T')] $1"
}

function error_message { # $1 error msg  $2 error code
        echo ""
        echo "########################################################"
        echo "Script to check the status of annotation jobs           "
        echo "########################################################"
        echo ""
        echo "ERROR ....."
        echo -e $1
        echo ""
        echo "Mandatory parameters:"
        echo -e "\t-V\tFull path of Variant directory  (where chromosome variant tables are stored - must begin with /)"
        echo -e "\t-N\tSample name (usually starting with IonXpress...)"
        echo ""
        if [ ! -z "$2" ]; then exit $2; fi
}

while getopts "V:N:" optionName; do
case "$optionName" in

V) workdir="$OPTARG";;                          ##  samples directory
N) sampleName="$OPTARG";;                       ##  sample name

esac
done

# Check mandatory parameters
if [[ -z "${workdir}" ]]; then error_message "Missing Variant dir (-V)" 99; fi
if [[ -z "${sampleName}" ]]; then error_message "Missing sample name value (-N)" 99; fi

TEST_ABSOLUTE_PATH=`echo ${workdir} | awk '{if ($0 ~ "^/") {print "TRUE"} else {print "FALSE"}}'`
if [[ "${TEST_ABSOLUTE_PATH}" == "FALSE" ]]; then
 error_message "\nVariant directory path (${workdir}) is not an absolute path!" 99
fi
if [[ ! -e "${workdir}" ]]; then
 error_message "\nVariant directory path (${workdir}) does not exist!" 99
fi

sample_nb=`echo ${sampleName} | awk '{split($0,part,"_"); print part[2]}'`
user=`whoami`
NSEC=20                  # number of seconds to wait for chr by chr annotation completion

nbParts=`ls -d1 ${workdir}/*/ | wc -l`
# nbPartsDoneCorrect=`ls -d1 ${workdir}/*/Union_samples_* | wc -l`
nbPartsRunning=`qstat -u ${user} | grep "X${sample_nb}_chr" | wc -l`

while [ "${nbPartsRunning}" != "0" ]
	do
	nbPartsRunning=`qstat -u ${user} | grep "X${sample_nb}_chr" | wc -l`
	nbPartsDoneNoVariants=`grep "All AVINPUT files are empty" ${workdir}/*/*OU | wc -l`
	if [ "${nbPartsRunning}" == "0" ]; then
		nbPartsDoneCorrect=`ls -d1 ${workdir}/*/Union_samples_* | wc -l`
	else
		nbPartsDoneCorrect='...'
	fi

	let "nbPartsDone=nbParts - nbPartsRunning"
	let "pctPartsDone=100*(nbPartsDone)/nbParts"
	if [ "${nbPartsDoneCorrect}" != "..." ]; then
		let "nbPartsDoneOtherErrors=nbPartsDone - (nbPartsDoneCorrect + nbPartsDoneNoVariants)"
	fi

	printf '\r%s' "[$(date +'%T')] Number of Parts : ${nbParts}	finished running: ${nbPartsDone}	completed correctly (s 0): ${nbPartsDoneCorrect}	completed with no variants detected (s 101): ${nbPartsDoneNoVariants}	completed with other errors (s <> 101) : ${nbPartsDoneOtherErrors}	% completed : ${pctPartsDone}	still running: ${nbPartsRunning}"
	echo `qstat -u ${user} | grep "X${sample_nb}_chr"`
	sleep ${NSEC}
done

if [ "${nbPartsRunning}" == "0" ]; then exit 0; fi
