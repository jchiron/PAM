#!/bin/bash
#AUTORUNDISABLE
VERSION="0.1"

#Now to have your plugin do some work or call other programs
OUTFILE=${RESULTS_DIR}/${PLUGINNAME}.html    
                                  
#echo "<html><body>raw data dir: ${RAW_DATA_DIR} <br/>" > $OUTFILE
#echo "analysis dir: ${ANALYSIS_DIR} <br/>" >> $OUTFILE
#echo "results dir: ${RESULTS_DIR} <br/>" >> $OUTFILE
#echo "analysis name: ${TSP_ANALYSIS_NAME} <br/>" >> $OUTFILE
#echo "path to BAM file: ${TSP_FILEPATH_BAM} <br/>" >> $OUTFILE
#echo "path to barcode list: ${TSP_FILEPATH_BARCODE_TXT} <br/>" >> $OUTFILE
#echo "project name: ${TSP_PROJECT} <br/>" >> $OUTFILE
#echo "run ID: ${TSP_RUNID} <br/>" >> $OUTFILE
#echo "PGM experiment ID: ${TSP_RUN_NAME} <br/>" >> $OUTFILE
#echo "Sample: ${TSP_SAMPLE} <br/>" >> $OUTFILE
#echo "vcf file: ${ANALYSIS_DIR}/plugin_out/variantCaller_out/${TSP_RUN_NAME}.vcf.zip <br/>" >> $OUTFILE

echo "<html><body>" > $OUTFILE
echo python ${DIRNAME}/variant_annotation_plugin.py \
--install-dir ${DIRNAME} \
--output-dir ${RESULTS_DIR} \
--report-dir ${ANALYSIS_DIR} \
--sample-list ${TSP_FILEPATH_BARCODE_TXT} \
--genome-fasta ${TSP_FILEPATH_GENOME_FASTA} \
--output-url ${TSP_URLPATH_PLUGIN_DIR} >> $OUTFILE

#set variable environment defining output annotation format to generate
if [[ "${PLUGINCONFIG__VA_DIAG}" != "True" ]]
then
export PLUGINCONFIG__VA_DIAG="False"
fi
if [[ "${PLUGINCONFIG__VA_SAFIR02}" != "True" ]]
then
export PLUGINCONFIG__VA_SAFIR02="False"
fi
echo "<br/> Diag ${PLUGINCONFIG__VA_DIAG}" >> $OUTFILE
echo "<br/> Safir02 ${PLUGINCONFIG__VA_SAFIR02}" >> $OUTFILE
 
echo "</body></html>" >> $OUTFILE

#echo ${TSP_URLPATH_PLUGIN_DIR}
#echo ${TSP_FILEPATH_PLUGIN_DIR}

${DIRNAME}/variant_annotation_plugin.py \
--install-dir  ${DIRNAME} \
--output-dir   ${RESULTS_DIR} \
--report-dir   ${ANALYSIS_DIR} \
--sample-list  ${TSP_FILEPATH_BARCODE_TXT} \
--genome-fasta ${TSP_FILEPATH_GENOME_FASTA} \
--output-url ${TSP_URLPATH_PLUGIN_DIR}