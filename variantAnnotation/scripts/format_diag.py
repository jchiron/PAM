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
import urllib2

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
    (options, args) = parser.parse_args()

    ANALYSIS_DIR                = os.environ['ANALYSIS_DIR']

###### Reprendre le sample name pour nommer les fichiers de sortie pour Diag

    file_tsv_path = options.tsv_path
    sample_id = options.sample_id
    sample_name = options.sample_name
    plugin_path = options.plugin_path
   
    non_decimal = re.compile(r'[^\d]+')
    columns="restrict"
  
    file_output_path = plugin_path + '/' + sample_id + '/CHR/' + sample_name + '_NGS_Diag_VC_TS.' + file_tsv_path.rsplit(".", 1)[1]
    file_temp_path = file_tsv_path.rsplit(".", 1)[0] + '.tmp'

    # print("output path : " + file_output_path)
    # print("tsv path : " + file_tsv_path)
    # print("sample id : " + sample_id)
    # print("sample name : " + sample_name)
    # print("plugin path : " + plugin_path)

    ALT_RATIO = ""
    DP4_ALT = ""
    DP4_TOT = ""
    reader = csv.reader(open(file_tsv_path), delimiter='\t')	
    columns_name = next(reader)
    #print columns_name
    for name in columns_name:
        #print "name : " + name
	if re.match(".*_ALT_RATIO",name):
		ALT_RATIO = name
  		print "ALT_RATION = " + ALT_RATIO
	if re.match(".*_DP4_TOT",name):
		DP4_TOT = name
		print "ADP4_TOT = " + DP4_TOT
        if re.match(".*_DP4_ALT",name):
		DP4_ALT = name
		print "DP4_ALT = " + DP4_ALT
		
    columns_final = "\
    Comm.Bio\t\
    CHROM\t\
    Ref. NM\t\
    TRANSCRIPTS\t\
    GENE\t\
    exon\t\
    c. (Annovar)\t\
    p. (Annovar)\t\
    c. (Mutalyzer)\t\
    p. (Mutalyzer)\t\
    %s\t\
    %s\t\
    %s\t\
    GENE_REGION\t\
    TYPE\t\
    START\t\
    REF\t\
    ALT\t\
    COSMIC\t\
    ESP Freq\t\
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
    CLINVAR\t\
    MUTATIONTASTER_SCORE\t\
    MUTATIONTASTER_INTERPRETATION\t\
    SIFT_SCORE\t\
    SIFT_INTERPRETATION\t\
    POLYPHEN2_HDIV_SCORE\t\
    POLYPHEN2_HDIV_INTERPRETATION\t\
    POLYPHEN2_HVAR_SCORE\t\
    POLYPHEN2_HVAR_INTERPRETATION\t\
    Amplicon\
    " % (ALT_RATIO,DP4_ALT,DP4_TOT)
    columns_final = columns_final.replace(" ", "").split("\t")

    #print columns_final

    columns_renamed = "\
    Comm.Bio\t\
    Chr\t\
    Ref.NM\t\
    TRANSCRIPTS\t\
    Gene\t\
    exon\t\
    c. (Annovar)\t\
    p. (Annovar)\t\
    c. (Mutalyzer)\t\
    p. (Mutalyzer)\t\
    Var.freq.\t\
    Var.Cov.\t\
    Pos.Cov.\t\
    Region\t\
    Type\t\
    Start_Position\t\
    Ref.seq\t\
    Var.seq\t\
    COSMIC\t\
    ESP Freq\t\
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
    CLINVAR\t\
    MUTATIONTASTER_SCORE\t\
    MUTATIONTASTER_INTERPRETATION\t\
    SIFT_SCORE\t\
    SIFT_INTERPRETATION\t\
    POLYPHEN2_HDIV_SCORE\t\
    POLYPHEN2_HDIV_INTERPRETATION\t\
    POLYPHEN2_HVAR_SCORE\t\
    POLYPHEN2_HVAR_INTERPRETATION\t\
    Amplicon\
    ".replace(" ", "").split("\t")

    columns_to_remove = "\
    STOP\t\
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
    ## debugging when target bed changed from review plan
    try :
    	api_url = os.environ['RUNINFO__API_URL']
    	pk =  os.environ['RUNINFO__PK']
    	url = api_url + '/v1/pluginresult/?format=json&plugin__name=variantCaller&result=%s' % pk
    	f = urllib2.urlopen(url)
    	d = json.loads(f.read())
    	#for variantcaller in d['objects']:
        #for barcodes in d['objects'][0]['store']['barcodes']:
            #variantcallerTargetBedPath = variantcaller['store']['targets_bed']
            #variantcallerTargetBedPath = barcodes['targets_bed']
            #break
        variantcallerTargetBedPath = d['objects'][0]['store']['barcodes'][sample_id]['targets_bed']
        #2 samples names with 2 barcodes (different panels)
        #if sample_id=='IonXpress_078' or sample_id=='IonXpress_082':
        #    variantcallerTargetBedPath = '/results/uploads/BED/25/hg19/unmerged/detail/IAD118795_231_Designed.with_NM.bed'
    except:
    	print "TARGET BED NOT FOUND WITH REST API"
    	variantcallerTargetBedPath = str(ion_params_json['experimentAnalysisSettings'].get('targetRegionBedFile'))
    file_bed_path = variantcallerTargetBedPath
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
    #for oncomine run : target file used in IR must not be modified use file containing GENE_ID\tNM to retreive transcrit for diag annotation
    if "oncomine" in variantcallerTargetBedPath.lower():
        #MET NM_001127500 PGM882 included
        #oncomine_NM_file = open("/media/PGM/doc_internes/Oncomine/oncomine_transcrits_refNM.tsv",'r')
        #MET NM_000245 start 17nov17
        oncomine_NM_file = open("/media/PGM/doc_internes/Oncomine/oncomine_transcrits_refNM_MET_000245.tsv",'r')
        dict_gene_nm = {}
        for line in oncomine_NM_file:
            line_parts=line.strip().split("\t")
            GENE_id = line_parts[0]
            if len(line_parts)==2: #some genes without default NM attributed at this time
                NM_id = line_parts[1]
            else:
                NM_id = ""
            dict_gene_nm[GENE_id] = NM_id
        oncomine_NM_file.close()
    elif "ocav3" in variantcallerTargetBedPath.lower():
        OCA_NM_file = open("/media/PGM/bioinfo/RefSeq_labo_310518.txt",'r')
        dict_gene_nm = {}
        for line in OCA_NM_file:
            line_parts=line.strip().split("\t")
            GENE_id = line_parts[0]
            print 'gene ' + GENE_id
            if len(line_parts)==2: #some genes without default NM attributed at this time
                NM_id = line_parts[1]
                print 'NM ' + NM_id
            else:
                NM_id = ""
            dict_gene_nm[GENE_id] = NM_id
        OCA_NM_file.close()
    else:        
        with open(file_bed_path) as f_bed:
            next(f_bed)
            last_gene_id = ""
            dict_gene_nm = {}
            for line in f_bed:
                line_parts=line.strip().split("\t")[-1].split(";")
                GENE_id = line_parts[0].split("=")[-1]
                #GENE_id = line_parts[0].split("=")[-1].split("_")[0]
                if (GENE_id != last_gene_id):
                    NM_id = "NM_" + line_parts[-1].split("_")[-1].split(".")[0]
                    dict_gene_nm[GENE_id] = NM_id
                    #print(dict_gene_nm)

    # open false positives file
    #fp_file = open(os.environ['DIRNAME'] + "/scripts/FalsePositives.csv", 'r')
    #fp_file = open("/media/PGM/bioinfo/known_mutations_info.tsv",'r')
    fp_file = open("/results/plugins/variantAnnotation/conf/known_mutations_info.tsv",'r')
    fp_reader = csv.reader(fp_file,delimiter='\t')

    # write file content
    for row in reader:
        # get reordered content
        row_content=list(reorderfunc(row))

        # fill / modify data
        # Keep only exonic & non synonymous rows
        # if (row_content[name2index["Region"]] == 'exonic' && row_content[name2index["Type"]] != 'synonymous'):
        #row_content[name2index["CHROM"]] = non_decimal.sub('', row_content[name2index["CHROM"]])
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
                    	row_content[name2index["c.(Annovar)"]] = nm_data[2]
                    	
                        # requete MUTALYZER pour obtenir p. plus précis qu'annovar et en hgvs
                        if re.match("c\.[0-9]+_[0-9]+[ATCG]+",nm_data[2]): # nomenclature hgvs : rajout 'delins' quand absent, pour que requete mutalyzer marche
                    		s = re.split("([ATCG]+)",nm_data[2])
                    		nm_data[2] = s[0] + "delins" + s[1]
                    		row_content[name2index["Comm.Bio"]] = row_content[name2index["Comm.Bio"]] + " Note : c.(Annovar) converted to %s for Mutalyzer request. " % nm_data[2]
                        try:
                        	succes = True
                        	url_mutalyzer = 'https://mutalyzer.nl/json/runMutalyzer?variant=%s:%s' % (ref_nm, nm_data[2])
    				f = urllib2.urlopen(url_mutalyzer)
    				d = json.loads(f.read())
    				if d['transcriptDescriptions']:
    					nuc_description = d['transcriptDescriptions'][0]
    					nuc_hgvs = nuc_description.split(":")[-1]
    					row_content[name2index["c.(Mutalyzer)"]] = nuc_hgvs
    				else:
    					succes = False
    					print ref_nm + ':' + nm_data[2] + ' : Failed to get c. from Mutalyzer.'
    				if d['proteinDescriptions']:
    					prot_description = d['proteinDescriptions'][0]
    					prot_hgvs = prot_description.split(":")[-1]#.replace("(","").replace(")","")
    					row_content[name2index["p.(Mutalyzer)"]] = prot_hgvs
    				else:
    					succes = False
    					print ref_nm + ':' + nm_data[2] + ' : Failed to get p. from Mutalyzer.'
    				if not succes:
    					row_content[name2index["Comm.Bio"]] = row_content[name2index["Comm.Bio"]] + " Mutalyzer : "
	    				for message in d['messages']:
	    					if message['errorcode'] != 'WNOVER':
	    						row_content[name2index["Comm.Bio"]] = row_content[name2index["Comm.Bio"]] + message['message']
	    				#link = " See https://mutalyzer.nl/name-checker?description=%s%%3A%s" % (ref_nm, nm_data[2])
	    				#row_content[name2index["Comm.Bio"]] = row_content[name2index["Comm.Bio"]] + link
    			except: 
    				print "Unexpected error:", sys.exc_info()[0]
    				print 'Failed to get data from Mutalyzer'
    				
    			# recherche si faux-positif connu 
                        for fp_row in fp_reader:
                            #print fp_row
                            #print "FP NM " + fp_row[1]
                            #print "FP c. " + fp_row[4]
                            #print "com. bio. " + fp_row[0]
                            if fp_row[1] == ref_nm and (fp_row[4] == row_content[name2index["c.(Annovar)"]] or fp_row[4] == row_content[name2index["c.(Mutalyzer)"]]):
                                row_content[name2index["Comm.Bio"]] = row_content[name2index["Comm.Bio"]] + fp_row[0]
    			fp_file.seek(0)
    				
                    if(len(nm_data)>=4):
                        row_content[name2index["p.(Annovar)"]] = nm_data[3]
                        #row_content[name2index["Codon"]] = nm_data[3][2:]

			
        #if(row_content[name2index["IonXpress_ALT_RATIO"]] != 'NA'):
        #    row_content[name2index["IonXpress_ALT_RATIO"]] = int(float(row_content[name2index["IonXpress_ALT_RATIO"]])*100)
        if(row_content[name2index[ALT_RATIO]] != 'NA'):
            row_content[name2index[ALT_RATIO]] = int(float(row_content[name2index[ALT_RATIO]])*100)

        #correct BRCA1 exon number
        if(row_content[name2index["GENE"]]== "BRCA1"):
            if (row_content[name2index["exon"]]!=''):
                if(int(row_content[name2index["exon"]])>3):
                    #print row_content[name2index["exon"]]
                    row_content[name2index["exon"]] = int(row_content[name2index["exon"]]) + 1
        #correct APC exon number
        if(row_content[name2index["GENE"]]== "APC"):
            if (row_content[name2index["exon"]]!=''):
                row_content[name2index["exon"]] = int(row_content[name2index["exon"]]) - 1


        # convert to tuple then write in file
        del(row_content[name2index["TRANSCRIPTS"]])
        row_content = tuple(row_content)
        # print(row_content)
        writer.writerow(row_content)
    fp_file.close()
    file_output.close()
    f.close()

    if (re.search(r'acrometrix',sample_name,re.I)==None) and (re.search(r'TRU-Q',sample_name,re.I)==None)and (re.search(r'temoin-positif-TRUQ',sample_name,re.I)==None):
        sample_name = sample_name.replace(" ","\ ")
        file_output_path = file_output_path.replace(" ","\ ")

        #add to sort mutations by type and frequencies in categories
        file_sort_path = plugin_path + '/' + sample_id + '/CHR/' + sample_name + '_NGS_Diag_sort.tsv'
    
    #print file_output_path
    #print file_sort_path
    

        #header
        cmd = "grep Comm.Bio " + file_output_path +" >" + file_sort_path
        os.system(cmd)

        #category 1 : exonic frameshift, missense, nonsense and sort decreasing frequency
        cmd = "grep -v -E \"splicing|synonymous|intronic|UTR|ncRNA|Comm.Bio\" " + file_output_path + '|sort -t \'\t\' -rn -k10,10 >> ' + file_sort_path
        os.system(cmd)
        cmd = "echo >>" + file_sort_path
        os.system(cmd)

        #category 2 : splicing and sort decreasing frequency
        cmd = "grep splicing " + file_output_path + "|sort -t \'\t\' -rn -k10,10 >> " + file_sort_path
        os.system(cmd)
        cmd = "echo >>" + file_sort_path
        os.system(cmd)

        #category 3 : exonic synonymous and sort decreasing frequency
        cmd = "grep synonymous " + file_output_path + "|sort -t \'\t\' -rn -k10,10 >> " + file_sort_path
        os.system(cmd)
        cmd = "echo >>" + file_sort_path
        os.system(cmd)
    
        #category 4 : intronic and sort decreasing frequency
        cmd = "grep intronic " + file_output_path + "|sort -t \'\t\' -rn -k10,10 >> " + file_sort_path
        os.system(cmd)
        cmd = "echo >>" + file_sort_path
        os.system(cmd)

        #category 5 : UTR and sort decreasing frequency
        cmd = "grep UTR " + file_output_path + "|sort -t \'\t\' -rn -k10,10 >> " + file_sort_path
        os.system(cmd)
        cmd = "echo >>" + file_sort_path
        os.system(cmd)

        #category 6 : ncRNA and sort decreasing frequency
        cmd = "grep ncRNA " + file_output_path + "|sort -t \'\t\' -rn -k10,10 >> " + file_sort_path
        os.system(cmd)

        #replace file_output NGS_diag_VC.tsv by the sorted one
        cmd = "cp " + file_sort_path + " " + file_output_path
        os.system(cmd)

        #delete temp files created
        os.remove(file_temp_path)
        os.remove(file_sort_path)
    return 0

if __name__ == "__main__":
    exit(plugin_main())
