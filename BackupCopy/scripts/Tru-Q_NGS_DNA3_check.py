#!/usr/bin/env python
import sys
import json
import os
import csv
import xlrd
import xlwt
from xlutils.copy import copy
from xlutils.styles import Styles
from ion.plugin import *
from subprocess import *
import urllib2
import re

def getPluginPath(run_pk,pluginName):
	pluginPath = False
	#try:
	api_url = os.environ['RUNINFO__API_URL']
	#print api_url
	#pk =  os.environ['RUNINFO__PK']
	#api_url = 'http://23TJQW1/rundb/api/v1/pluginresult/?format=json&plugin__name=%s&result=%s' % (pluginName,str(run_pk))
	api_url = api_url + '/v1/pluginresult/?format=json&plugin__name=%s&result=%s' % (pluginName,str(run_pk))
	#print api_url
	f = urllib2.urlopen(api_url)
	d = json.loads(f.read())
		#print d
	for plugin in d['objects']: # already sorted from latest to oldest
			#print plugin
		if plugin['state'] == 'Completed':
			pluginPath = plugin['path']
			pluginPath = pluginPath.split('/')[-1]
			break
	if not pluginPath:
		print 'WARNING! No completed instance found for plugin "%s"' % pluginName
	#except:
	#	print 'ERROR!  Failed to get plugin path'
	
	print 'INFO: using %s path: %s' % (pluginName,pluginPath)
	return pluginPath
	
def representsInt(s): # pour eviter avertissement "nombre ecrit en texte" sous excel
	try: 
		s = int(s)
		return s
	except ValueError:
		return s

################################################################################
HD_names = ["Temoin_TRUQ_DNA3","Temoin-positif-TRU-Q","temoin-positif-TRU-Q","TEMOIN-TRUQ","TEMOIN-POSITIF-TRU-Q","Temoin TRU-Q Horizon","TEMOIN-TRU-Q","Temoin-TRUQ-Q","Temoin_TRUQ_DNA3_EQ","Temoin_HD730"]
xls_path = "/media/PGM/results_S5/suivi_temion_positif_horizon_Tru-Q_NGS_DNA3.xls"

run_dir = sys.argv[1]
run_pk = sys.argv[2]
os.chdir(run_dir)

variantCallerPath = getPluginPath(run_pk,"variantCaller")
#variantCallerPath = "variantCaller_out.2100"
variantAnnotationPath = getPluginPath(run_pk,"variantAnnotation")
#variantAnnotationPath = "variantAnnotation_out.2101"

try:
	with open('plugin_out/%s/barcodes.json' % variantCallerPath, 'r') as fh:
		json_dat = json.load(fh)
except:
	print 'Error reading plugin json.'
	
barcode = False
for ionxpress in json_dat.keys():
	#if json_dat[ionxpress]['sample'] in HD_names:
        for ctl_name in HD_names :
		if ctl_name in json_dat[ionxpress]['sample'] :
		#if re.search(ctl_name,json_dat[ionxpress]['sample'],re.I):
			sample = json_dat[ionxpress]['sample']
			#print sample
			barcode = ionxpress
		
                        if barcode:
                                va_path = "plugin_out/" + variantAnnotationPath + "/%s/CHR/%s_NGS_Diag_VC_TS.tsv" % (barcode,sample)
                        else:
                                print "ERROR :  Tru-Q_NGS_DNA3 sample not found. aborting."
                                exit()

                        # xls styles

                        headerStyle = xlwt.easyxf('alignment: horiz centre, wrap on; font: name Calibri, bold on, height 220; borders: left thin, top thin, bottom thin, right thin')
                        generalStyle = xlwt.easyxf('alignment: horiz left; font: name Calibri, height 220; borders: left thin, top thin, bottom thin, right thin')

                        # i/o
                        va_file = open(va_path,"r")	
                        va_reader = csv.reader(va_file,delimiter="\t")

                        xls_reader = xlrd.open_workbook(xls_path,formatting_info=True)
                        xls_sheet = xls_reader.sheet_by_index(0)
                        s = Styles(xls_reader)

                        final_xls = copy(xls_reader)
                        final_sheet = final_xls.get_sheet(0)

                        # header 1
                        header1 = xls_sheet.row(0)
                        if run_dir.split('/')[-1] == "":
                                run_name = run_dir.split('/')[-2].replace('Auto_user_','')
                        else:
                                run_name = run_dir.split('/')[-1].replace('Auto_user_','')
                        #print run_name
                        final_sheet.row(0).write(len(header1),run_name+' '+sample,headerStyle)

                        # header 2
                        header2 = xls_sheet.row(1)
                        final_sheet.row(1).write(len(header2),"Var.freq",generalStyle)

                        # variants lines
                        for i in range(xls_sheet.nrows-2):
                                variant_line = xls_sheet.row(i+2)
                                variant2check = (xls_sheet.cell(i+2,5).value,xls_sheet.cell(i+2,6).value)
	
                                va_file.seek(0)
                                va_reader.next()
                                for va_line in va_reader:
                                        #print va_line
                                        if va_line:
                                                va_variant = (va_line[7],va_line[8])
                                                va_variant_freq = "?"
                                                if variant2check == va_variant:
                                                        va_variant_freq = va_line[9]
                                                        break
                                final_sheet.row(i+2).write(len(variant_line),representsInt(va_variant_freq),generalStyle)

                        cmd = Popen(["rm","-f",xls_path], stdout=PIPE)
                        out, err = cmd.communicate()
                        print '-Removing previous xls : \nOUT: %s\nERR: %s'%(out, err)
                        
                        final_xls.save(xls_path) # overwrite previous xls
