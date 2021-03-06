#!/usr/bin/env python

import sys
import os
import subprocess
import json
import time
import traceback
import re
from optparse import OptionParser
from django.conf import settings
from django.template.loader import render_to_string


# critical environment variables:
DIRNAME                     = '' # home directory for the plugin files
TSP_FILEPATH_PLUGIN_DIR      = ''
startplugin_json            = {}

# File names generated by the plugin
BASENAME_VARIANTS_TSV      = 'Union_samples_chrAll.txt'
#BASENAME_PARAMETERS_JSON    = 'local_parameters.json'
HTML_BLOCK                  = 'variantAnnotation_block.html'    # Sample report page block


def plugin_main():
        
    #global PLUGIN_DEV_SKIP_VARIANT_CALLING    
    global DIRNAME
    global TSP_FILEPATH_PLUGIN_DIR
    global startplugin_json

    global UCSC_GENOME_FASTA

    # script parameters
    parser = OptionParser()
    parser.add_option('-a', '--annotation-dir', help='Directory containing annotation results', dest='annotation_dir')
    parser.add_option('-s', '--sample-id', help='IonXpress sample name', dest='sample_id')
    
    (options, args) = parser.parse_args()

    DIRNAME                     = os.environ['DIRNAME'] # home directory for the plugin files
    TSP_FILEPATH_PLUGIN_DIR     = os.environ['TSP_FILEPATH_PLUGIN_DIR'] # target plugin results directory
    ANNOTATION_DIR              = options.annotation_dir # annotation results folder containing IonXpress sub-folders
    SAMPLE_ID                   = options.sample_id # sample IonXpress name

    
    if ANNOTATION_DIR is not None and SAMPLE_ID is not None :

        #rename file
        cmd = 'cp '+ ANNOTATION_DIR + '/' + SAMPLE_ID + '/CHR/' + BASENAME_VARIANTS_TSV + ' ' + ANNOTATION_DIR + '/' + SAMPLE_ID + '/CHR/' + SAMPLE_ID + '_all_annot.tsv'

        print cmd +'\n' #debug line
        os.system(cmd)

        # TODO: filter all annotation file (restrict columns and transcript id)
        #xls file

        #render html block
        tsv_file_path = ANNOTATION_DIR + '/' + SAMPLE_ID + '/CHR/' + SAMPLE_ID + '_all_annot.tsv'
        # TODO xls_file_path
        xls_file_path = ANNOTATION_DIR + '/' + SAMPLE_ID + '/CHR/' + SAMPLE_ID + '_annot.xls'

        render_context = {
            'sample_id' : SAMPLE_ID,
            'variants_tsv_link' : tsv_file_path,
            'variants_xls_link' : xls_file_path
            }

        html_block_path = ANNOTATION_DIR + '/' + SAMPLE_ID
        out = open(html_block_path + '/' + HTML_BLOCK, 'w'
    else:
        print 'ERROR: missing argument\n'
