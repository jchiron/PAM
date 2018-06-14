#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import subprocess
import json
import time
import traceback
import re
import glob
import ast
import csv
import operator
from optparse import OptionParser
from django.conf import settings
from django.template.loader import render_to_string

def printtime(message, *args):
    if args:
        message = message % args
    print "[ " + time.strftime('%X') + " ] " + message
    sys.stdout.flush()
    sys.stderr.flush()


def run_command(command,description):
    printtime(' ')
    printtime('Task    : ' + description)
    printtime('Command : ' + command)
    printtime(' ')
    return subprocess.call(command,shell=True)


def execute_output(cmd):
    try:
        process = subprocess.Popen(cmd, stderr=subprocess.STDOUT, stdout=subprocess.PIPE, shell=True)
        return process.communicate()[0]
    except:
        traceback.print_exc()
        return ''


def plugin_main():
    print("Starting conversion")
    # script parameters
    parser = OptionParser()
    parser.add_option('-t', '--tsv-file', help='Directory containing plugin files', dest='tsv_path')
    parser.add_option('-i', '--sample-id', help='Sample ID', dest='sample_id')
    parser.add_option('-s', '--sample-name', help='Sample name', dest='sample_name')
    parser.add_option('-p', '--plugin-path', help='Plugin output path', dest='plugin_path')
    # parser.add_option('-b', '--bed-file', help='Directory for results files', dest='bed_path')
    (options, args) = parser.parse_args()

    ANALYSIS_DIR                = os.environ['ANALYSIS_DIR']

###### Reprendre le sample name pour nommer les fichiers de sortie pour Safir 02

    file_tsv_path = options.tsv_path
    sample_id = options.sample_id
    sample_name = options.sample_name
    plugin_path = options.plugin_path
    # file_bed_option = options.bed_path

    non_decimal = re.compile(r'[^\d]+')
    columns="restrict"
    # print(file_tsv_path)
    # print(file_bed_option)
    # file_tsv_path = "/media/sf_k.tran/Documents/Tools/pgm-plugins/_custom/annotation_html_file"
    # file_bed_path = "C:\\Users\\k.tran\\Documents\\Bergonié\\Projects\\Safir02\\BedFiles\\IAD43772_Designed_v2.bed"
    # file_tsv_path = "C:\\Users\\k.tran\\Documents\\Tools\\pgm-plugins\\_custom\\annotation_html_file\\IonXpress_018_annot.tsv"

    # file_output_path = file_tsv_path.rsplit(".", 1)[0] + '_NGS_VC_TS.' + file_tsv_path.rsplit(".", 1)[1]
    file_output_path = plugin_path + '/' + sample_id + '/CHR/' + sample_name + '_NGS_VC_TS.' + file_tsv_path.rsplit(".", 1)[1]
    file_temp_path = file_tsv_path.rsplit(".", 1)[0] + '.tmp'

    # print("output path : " + file_output_path)
    # print("tsv path : " + file_tsv_path)
    # print("sample id : " + sample_id)
    # print("sample name : " + sample_name)
    # print("plugin path : " + plugin_path)
    
    columns_final = "\
    Concl.Bio\t\
    Comm.Bio\t\
    Var.Class\t\
    GENE\t\
    p.\t\
    exon\t\
    IonXpress_ALT_RATIO\t\
    IonXpress_DP4_TOT\t\
    Ref.NM\t\
    TRANSCRIPTS\t\
    c.\t\
    Codon\t\
    CHROM\t\
    START\t\
    STOP\t\
    Strand\t\
    GENE_REGION\t\
    TYPE\t\
    REF\t\
    ALT\t\
    IonXpress_DP4_ALT\t\
    Variant_Class\t\
    ESP Freq\t\
    COSMIC\t\
    dbSNP\t\
    1000G_ALL_AF\t\
    1000G_ALL_PROBS_AA:AB:BB\t\
    1000G_EUR_AF\t\
    1000G_EUR_PROBS_AA:AB:BB\t\
    1000G_AMR_AF\t\
    1000G_AMR_PROBS_AA:AB:BB\t\
    1000G_ASN_AF\t\
    1000G_ASN_PROBS_AA:AB:BB\t\
    1000G_AFR_AF\t\
    1000G_AFR_PROBS_AA:AB:BB\t\
    SIFT_SCORE\t\
    SIFT_INTERPRETATION\t\
    POLYPHEN2_HDIV_SCORE\t\
    POLYPHEN2_HDIV_INTERPRETATION\t\
    POLYPHEN2_HVAR_SCORE\t\
    POLYPHEN2_HVAR_INTERPRETATION\t\
    Amplicon\
    ".replace(" ", "").split("\t")

    columns_renamed = "\
    Concl.Bio\t\
    Comm.Bio\t\
    Var.Class\t\
    Gene\t\
    p.\t\
    exon\t\
    Var.freq.\t\
    Pos.Cov.\t\
    Ref.NM\t\
    TRANSCRIPTS\t\
    c.\t\
    Codon\t\
    Chr.\t\
    Start_Position\t\
    Stop_Position\t\
    Strand\t\
    Region\t\
    Type\t\
    Ref.seq\t\
    Var.seq\t\
    Var.Cov.\t\
    Variant_Class\t\
    ESP Freq\t\
    COSMIC\t\
    dbSNP\t\
    1000G_ALL_AF\t\
    1000G_ALL_PROBS_AA:AB:BB\t\
    1000G_EUR_AF\t\
    1000G_EUR_PROBS_AA:AB:BB\t\
    1000G_AMR_AF\t\
    1000G_AMR_PROBS_AA:AB:BB\t\
    1000G_ASN_AF\t\
    1000G_ASN_PROBS_AA:AB:BB\t\
    1000G_AFR_AF\t\
    1000G_AFR_PROBS_AA:AB:BB\t\
    SIFT_SCORE\t\
    SIFT_INTERPRETATION\t\
    POLYPHEN2_HDIV_SCORE\t\
    POLYPHEN2_HDIV_INTERPRETATION\t\
    POLYPHEN2_HVAR_SCORE\t\
    POLYPHEN2_HVAR_INTERPRETATION\t\
    Amplicon\
    ".replace(" ", "").split("\t")

    columns_to_remove = "\
    UCSC_LINK\t\
    QUICKGO_LINK\
    ".replace(" ", "").split("\t")

    ########## functions ##########
    # Returns a new list. "lst" is not modified.
    def delete_by_indices(lst, indices):
        s = set(indices)
        return [ value for (i, value) in enumerate(lst) if i not in s ]

    ######################################################################################
    #read startplugin_json file, extract bed file informations
    try:
        json_file = open(os.path.join(ANALYSIS_DIR,'ion_params_00.json'), 'r')
        ion_params_json = json.load(json_file,parse_float=str)
        json_file.close()
    except:
        printtime('ERROR: Failed to load and parse ion_params_00.json')
        return 1
    file_bed_path = str(ion_params_json['experimentAnalysisSettings'].get('targetRegionBedFile'))
    print("Using the following bed file : " + file_bed_path)

    ########## add missing columns ##########
    reader = csv.reader(open(file_tsv_path), delimiter='\t')
    columns_current=next(reader)
    name2index = dict((name, index) for index, name in enumerate(columns_current))
    columns_to_remove_index = [name2index[name] for name in columns_to_remove]
    # using List Comprehensions instead of substracting sets, preserve list order
    columns_to_keep = [x for x in columns_current if x not in columns_to_remove]
    columns_to_add = [x for x in columns_final if x not in columns_to_keep]
    columns_names_all = columns_to_add + columns_to_keep
    # open output file
    file_temp = open(file_temp_path, 'w')
    writer = csv.writer(file_temp, delimiter='\t')
    # write header
    writer.writerow(columns_names_all)
    # write file content
    for row in reader:
        col_to_keep = delete_by_indices(row, columns_to_remove_index)
        rows_all = ['']*len(columns_to_add) + col_to_keep
        writer.writerow(rows_all)
    file_temp.close()

    ########## reorder columns ##########
    f = open(file_temp_path)
    reader = csv.reader(f, delimiter='\t')
    columns_current=next(reader)
    # define the index of each column, based on name or based on index
    name2index = dict((name, index) for index, name in enumerate(columns_current))
    index2name = dict((index, name) for index, name in enumerate(columns_current))
    columns_custom_index = [name2index[name] for name in columns_final]
    # define the remaining columns
    l3 = [x for x in list(name2index.keys()) if x not in columns_final]
    # print(l3)
    if (columns == "all"):
        remaining_index=[name2index[name] for name in l3]
        columns_final_index=columns_custom_index+remaining_index
    else:
        columns_final_index=columns_custom_index
    # print(columns_final_index)
    reorderfunc = operator.itemgetter(*columns_final_index)

    # list - define column header
    columns_final_order=[]
    for idx, val in enumerate(columns_final_index):
        columns_final_order.append(index2name[val])
    # print(columns_final_order)

    # redefine name2index and index2name
    name2index = dict((name, index) for index, name in enumerate(columns_final_order))
    index2name = dict((index, name) for index, name in enumerate(columns_final_order))

    # open output file
    file_output = open(file_output_path, 'w')
    writer = csv.writer(file_output, delimiter='\t')
    # write header
    del(columns_renamed[name2index["TRANSCRIPTS"]])
    # print(columns_renamed)
    writer.writerow(columns_renamed)

    # Get NM data
    with open(file_bed_path) as f_bed:
        next(f_bed)
        last_gene_id = ""
        dict_gene_nm = {}
        for line in f_bed:
            line_parts=line.strip().split("\t")[-1].split(";")
            GENE_id = line_parts[0].split("=")[-1]
            if (GENE_id != last_gene_id):
                NM_id = "NM_" + line_parts[-1].split("_")[-1].split(".")[0]
                dict_gene_nm[GENE_id] = NM_id
        # print(dict_gene_nm)

    # write file content
    for row in reader:
        # get reordered content
        row_content=list(reorderfunc(row))

        # fill / modify data
        row_content[name2index["CHROM"]] = non_decimal.sub('', row_content[name2index["CHROM"]])
        ref_nm = dict_gene_nm.get(row_content[name2index["GENE"]])
        row_content[name2index["Ref.NM"]] = ref_nm
        if(row_content[name2index["TRANSCRIPTS"]] != 'NA'):
            if(row_content[name2index["TRANSCRIPTS"]].endswith(',')):
                transcripts = row_content[name2index["TRANSCRIPTS"]][:-1]
            for transcript_data in transcripts.split(","):
                nm_data = transcript_data.split(":")
                if(nm_data[0] == ref_nm):
                    if(len(nm_data)>=2):
                        row_content[name2index["exon"]] = non_decimal.sub('', nm_data[1])
                    if(len(nm_data)>=3):
                        row_content[name2index["c."]] = nm_data[2]
                    if(len(nm_data)>=4):
                        row_content[name2index["p."]] = nm_data[3]
                        row_content[name2index["Codon"]] = nm_data[3][2:]
        if(row_content[name2index["IonXpress_ALT_RATIO"]] != 'NA'):
            row_content[name2index["IonXpress_ALT_RATIO"]] = int(float(row_content[name2index["IonXpress_ALT_RATIO"]])*100)
        # convert to tuple then write in file
        del(row_content[name2index["TRANSCRIPTS"]])
        row_content = tuple(row_content)
        # print(row_content)
        writer.writerow(row_content)
    file_output.close()
    f.close()
    # delete the temporary file created
    os.remove(file_temp_path)
    return 0

if __name__ == "__main__":
    exit(plugin_main())
