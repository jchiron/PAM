#!/usr/bin/env python
# Copyright (C) 2013 Ion Torrent Systems, Inc. All Rights Reserved

import sys
import os
import subprocess
import json
import time
import traceback
import re
import glob
import ast
from optparse import OptionParser
from django.conf import settings
from django.template.loader import render_to_string


def printtime(message, *args):
	if args:
		message = message % args
	print "[ " + time.strftime('%X') + " ] " + message
	sys.stdout.flush()
	sys.stderr.flush()


def plugin_main():

	# script parameters
	parser = OptionParser()
	parser.add_option('-s', '--sample-name', help='Sample name', dest='s_name')
	parser.add_option('-i', '--sample-id', help='Sample barcode', dest='s_id')
	parser.add_option('-v', '--sample-vcf', help='Sample .vcf file from variantCaller', dest='s_vcf')
	parser.add_option('-d', '--install-dir', help='Directory containing plugin files', dest='install_dir') # mettre /results/plugins/Check_contamination/variantAnnotation_custom
	parser.add_option('-o', '--output-dir', help='Directory for results files', dest='output_dir') # mettre .../Check_contamination_out.*/variantAnnotation/
	parser.add_option('-f', '--genome-fasta', help='Reference genome fasta file', dest='genome_fasta')
	(options, args) = parser.parse_args()

	SAMPLE_NAME = options.s_name # sample name to annotate
	SAMPLE_ID = options.s_id  # sample barcode to annotate
	SAMPLE_VCF = options.s_vcf  # sample .vcf file 
	DIRNAME = options.install_dir # home directory for the plugin files
	TSP_FILEPATH_PLUGIN_DIR = options.output_dir # target plugin results directory
	PGM_GENOME_FASTA= options.genome_fasta # genome .fasta

	######################################################################################

	# retreive genome version to read genome parameter file
	regexp = ".*/(.*).fasta$"
	genome_id = re.match(regexp, PGM_GENOME_FASTA)
	if genome_id.group(1) is not None:
		if re.match("hg19", genome_id.group(1)):
			GENOME_PARAMETER = DIRNAME +"/conf/hg19.conf"

	# read genome parameter file and construct genome UCSC files path
	genome_parameter_file = open(GENOME_PARAMETER,'r')
	for line in genome_parameter_file:
		line = line.rstrip('\n')
		if line.startswith('UCSC_FASTA'):
			UCSC_GENOME_FASTA = line.split(',')[1]
		if line.startswith('UCSC_CHR_INFO'):
			UCSC_CHR_TAB_INFO = line.split(',')[1]

	
	DNASEQ_VARIANTS="/results/tools/VARIANT_ANNOTATION_V5_2014_07_04"
	IonXpress_ANNOTATION="/results/tools/IonXpress_VCF_conversion_2014_07_04"
	ANNODB="/results/tools/annovar/humandb"
	BCFDIR="/usr/bin"
	SAMDIR="/usr/bin"
	SRF="SRF"
	SRR="SRR"
	SAF="SAF"
	SAR="SAR"
	
	SAMPLE_ID_CHR = TSP_FILEPATH_PLUGIN_DIR + "/" + SAMPLE_ID + "/CHR"
	SAMPLE_ID_BCF = TSP_FILEPATH_PLUGIN_DIR + "/" + SAMPLE_ID + "/BCF"
	
	# creation sous dossier sample / CHR et BCF
	cmd = "mkdir -p " + SAMPLE_ID_CHR
	os.system(cmd)
	cmd = "mkdir -p " + SAMPLE_ID_BCF
	os.system(cmd)

	sample_id_txt = os.path.join(TSP_FILEPATH_PLUGIN_DIR,SAMPLE_ID + '.txt')
	sample_plan=open(sample_id_txt,'w')
	sample_plan.write("Header \n" + SAMPLE_ID)
	sample_plan.close()

	#########################
	# PIPELINE D'ANNOTATION #
	#########################
	
	print "Starting annotation pipeline :"
	#1# convert_TorSer_VCF_generate_BCF.sh
	print "Execute convert_TorSer_VCF_generate_BCF.sh"
	cmd = "bash /results/plugins/Check_contamination/variantAnnotation_custom/convert_TorSer_VCF_generate_BCF.sh -i " + SAMPLE_VCF +\
	    " -o " + SAMPLE_ID_BCF + "/" + SAMPLE_ID +\
	    " -b " + SAMPLE_ID +\
	    " -U " + UCSC_GENOME_FASTA +\
	    " -C yes " +\
	    " -1 " + SRF +\
	    " -2 " + SRR +\
	    " -3 " + SAF +\
	    " -4 " + SAR +"\n"
	    
	os.system(cmd)

	#2# VarPipe_process_batch_samples.V2.sh
	print "Execute VarPipe_process_batch_samples.V2.sh"
	cmd = "bash /results/tools/VARIANT_ANNOTATION_V5_2014_07_04/VarPipe_process_batch_samples.V2.sh -S " + sample_id_txt +\
	    " -N " + SAMPLE_ID +\
	    " -I " + TSP_FILEPATH_PLUGIN_DIR +\
	    " -O " + SAMPLE_ID_CHR +\
	    " -F " + SAMPLE_ID_BCF +\
	    " -B BCF" +\
	    " -D 999999999" +\
	    " -V 0.1.19" +\
	    " -C yes " +\
	    " -T " + UCSC_CHR_TAB_INFO +\
	    " -P /results/tools/VARIANT_ANNOTATION_V5_2014_07_04" +\
	    " -A " + ANNODB +"\n"

	os.system(cmd)
	   
	#3# check_variant_annotation_queue.sh
	print "Execute check_variant_annotation_queue.sh"
	cmd = "bash /results/tools/VARIANT_ANNOTATION_V5_2014_07_04/check_variant_annotation_queue.sh -V " + SAMPLE_ID_CHR +\
	    " -N " + SAMPLE_ID +"\n"
	    
	os.system(cmd)

	#4# MergeVar_all_CHRs.V2.sh
	print "Execute MergeVar_all_CHRs.V2.sh"
	cmd = "bash /results/tools/VARIANT_ANNOTATION_V5_2014_07_04/MergeVar_all_CHRs.V2.sh -V " + SAMPLE_ID_CHR +\
	    " -A yes" +\
	    " -S no" +\
	    " -Z no" +\
	    " -C yes" +"\n"
	    
	os.system(cmd)

	cmd = "cp " + SAMPLE_ID_CHR + "/Union_samples_chrAll.txt " + SAMPLE_ID_CHR + "/" + SAMPLE_ID +"_annot.tsv"
	os.system(cmd)
	
	# diag formating routine
	print "Execute format_diag.py"
	cmd = "python /results/plugins/variantAnnotation/scripts/format_diag.py --tsv-file " + SAMPLE_ID_CHR + "/" + SAMPLE_ID + "_annot.tsv" +\
	    " --sample-id " + SAMPLE_ID +\
	    " --sample-name " + SAMPLE_NAME +\
	    " --plugin-path " + TSP_FILEPATH_PLUGIN_DIR +"\n"

        os.system(cmd)

	return 0

if __name__ == "__main__":
	
	exit(plugin_main())
