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

envDict = dict(os.environ)

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
HD_names = ["TEMOIN_BRCA_SOMATIC","temoin_brca_somatic","TEMOIN-BRCA-SOMATIC","temoin-brca-somatic","TEMOINBRCA_SOMATIC"]
xls_path = "/media/PGM/results_S5/suivi_temoin_brca_somatique_multiplex.xls"

run_dir = sys.argv[1]
run_pk = sys.argv[2]
os.chdir(run_dir)

variantCallerPath = getPluginPath(run_pk,"variantCaller")
#variantCallerPath = "variantCaller_out.2100"
variantAnnotationPath = getPluginPath(run_pk,"variantAnnotation")
#variantAnnotationPath = "variantAnnotation_out.2101"
if run_dir.split('/')[-1] == "":
	run_name = run_dir.split('/')[-2].replace('Auto_user_','')
else:
	run_name = run_dir.split('/')[-1].replace('Auto_user_','')
#print run_name

# get the project name.
try:
	project = envDict['TSP_PROJECT']
except:
	project = 'unknown'
			
# get the run name.
try:
	runName = envDict['TSP_ANALYSIS_NAME']
except:
	runName = 'unknown'
backupPath = project + '/' + runName

temoin_xls_path = '/media/PGM/results_S5/%s/temoin_brca_somatic_multiplex.xls' % backupPath
print "fichier temoin du run : %s" %temoin_xls_path

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
	print "ERROR : TEMOIN_BRCA_SOMATIC sample not found. aborting."
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

temoin_xls_workbook = xlwt.Workbook()
temoin_xls_sheet = temoin_xls_workbook.add_sheet('variants')

# header 1
header1 = xls_sheet.row(0)

final_sheet.row(0).write(len(header1),run_name,headerStyle)
temoin_xls_sheet.row(0).write(0,header1[0].value,headerStyle)
temoin_xls_sheet.row(0).write(1,header1[1].value,headerStyle)
temoin_xls_sheet.row(0).write(2,header1[2].value,headerStyle)
temoin_xls_sheet.row(0).write(3,header1[3].value,headerStyle)
temoin_xls_sheet.row(0).write(4,header1[4].value,headerStyle)
temoin_xls_sheet.row(0).write(5,header1[5].value,headerStyle)
temoin_xls_sheet.row(0).write(6,header1[6].value,headerStyle)
temoin_xls_sheet.row(0).write(7,header1[7].value,headerStyle)
temoin_xls_sheet.row(0).write(8,run_name,headerStyle)

# header 2
header2 = xls_sheet.row(1)
final_sheet.row(1).write(len(header2),"Var.freq",generalStyle)
temoin_xls_sheet.row(1).write(0,header2[0].value,generalStyle)
temoin_xls_sheet.row(1).write(1,header2[1].value,generalStyle)
temoin_xls_sheet.row(1).write(2,header2[2].value,generalStyle)
temoin_xls_sheet.row(1).write(3,header2[3].value,generalStyle)
temoin_xls_sheet.row(1).write(4,header2[4].value,generalStyle)
temoin_xls_sheet.row(1).write(5,header2[5].value,generalStyle)
temoin_xls_sheet.row(1).write(6,header2[6].value,generalStyle)
temoin_xls_sheet.row(1).write(7,header2[7].value,generalStyle)
temoin_xls_sheet.row(1).write(8,"Var. freq",generalStyle)

# variants lines
for i in range(xls_sheet.nrows-2):
	variant_line = xls_sheet.row(i+2)
	variant2check = (str(xls_sheet.cell(i+2,5).value),str(xls_sheet.cell(i+2,6).value))
	#print 'V2C>'
	#print variant2check
	#print '<V2C'

	va_file.seek(0)
	va_reader.next()
	for va_line in va_reader:
		#print va_line
		if va_line:
			if variant2check == ('c.7397T>C','p.Val2466Ala'):
				#print 'mut to modif'
				variant2check = ('c.7397T>C','p.V2466A')
				#print 'V2C modif >'
				#print variant2check
				#print'< V2C modif'
				va_variant = (va_line[5],va_line[6])
				
			else:
				va_variant = (va_line[7],va_line[8])
			#print va_variant
			va_variant_freq = "?"
			if variant2check == va_variant:
				va_variant_freq = va_line[9]
				break
	final_sheet.row(i+2).write(len(variant_line),representsInt(va_variant_freq),generalStyle)
	temoin_xls_sheet.row(i+2).write(0,variant_line[0].value,generalStyle)
	temoin_xls_sheet.row(i+2).write(1,variant_line[1].value,generalStyle)
	temoin_xls_sheet.row(i+2).write(2,variant_line[2].value,generalStyle)
	temoin_xls_sheet.row(i+2).write(3,variant_line[3].value,generalStyle)
	temoin_xls_sheet.row(i+2).write(4,variant_line[4].value,generalStyle)
	temoin_xls_sheet.row(i+2).write(5,variant_line[5].value,generalStyle)
	temoin_xls_sheet.row(i+2).write(6,variant_line[6].value,generalStyle)
	temoin_xls_sheet.row(i+2).write(7,variant_line[7].value,generalStyle)
	temoin_xls_sheet.row(i+2).write(8,representsInt(va_variant_freq),generalStyle)

temoin_xls_workbook.save(temoin_xls_path)

cmd = Popen(["rm","-f",xls_path], stdout=PIPE)
out, err = cmd.communicate()
print '-Removing previous xls : \nOUT: %s\nERR: %s'%(out, err)

final_xls.save(xls_path) # overwrite previous xls
