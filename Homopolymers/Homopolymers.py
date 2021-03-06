#!/usr/bin/env python
import os
import json
import csv
from ion.plugin import *
from subprocess import *
import glob
import urllib2
import numpy
import matplotlib
matplotlib.use('Agg')
from pylab import *
import xlwt
import xlrd
import xlsxwriter


class Homopolymers(IonPlugin):
	"""Plugin object searching ins/del mutations in homopolymers regions 6n+ in BRCA1 BRCA2"""
	version = "0.1"
	envDict = dict(os.environ)
	json_dat = {}
	barcodeNames = []
	sampleNameLookup = {} # dictionary allowing us to get the sample name associated with a particular barcode name

	def launch(self):
		try:
			with open('startplugin.json', 'r') as fh:
				self.json_dat = json.load(fh)
		except:
			print 'Error reading plugin json.'

		try:
			htmlOut = open('%s/plugin_out/Homopolymers_out/Homopolymers_block.html'%self.envDict['ANALYSIS_DIR'], 'a+')
		except:
			htmlOut = open('Homopolymers_block.html', 'w')
		htmlOut.write('<html><body>\n')
	
		if os.path.isfile(self.envDict['TSP_FILEPATH_BARCODE_TXT']):
			self.isBarcodedRun = True
			
		if not self.isBarcodedRun:
			print "This plugin only work for barcoded run for now. \n Exiting..."
			sys.exit(0)
		
	
		# set defaults for user options
		variantCallerPath = False
		variantCallerInstance = False
		result_dir = False
		
		try:
			result_dir=self.json_dat['runinfo']['results_dir']
			
		except:
			print 'Warning: plugin does not appear to be configured'
		
		try:
			api_url = self.json_dat['runinfo']['api_url'] + '/v1/pluginresult/?format=json&plugin__name=variantCaller&result=' + str(self.json_dat['runinfo']['pk'])
			api_key = self.json_dat['runinfo'].get('api_key', None)
			if api_key is not None:
				api_url = api_url + '&api_key=%s' % api_key
				print 'Using API key: %s' % api_key
			else:
				print 'No API key available'
			f = urllib2.urlopen(api_url)
			d = json.loads(f.read())
	
			for plugin in d['objects']: # they are already sorted from newest to oldest, so take the first of this list which is 'completed'
				if plugin['state'] == 'Completed':
					variantCallerPath = plugin['path']
					variantCallerInstance = variantCallerPath.split('.')[-1]
					break
			if not variantCallerPath:
				print 'WARNING! No completed instance found for plugin variantCallerPath'
		except:
			print 'ERROR!  Failed to get variantCallerPath'
	
		# Get bam filenames.
		with open(os.path.join(self.json_dat['runinfo']['basecaller_dir'], 'datasets_basecaller.json'), 'r') as f:
			json_basecaller = json.load(f)

#		samples = json.loads(self.json_dat['plan']['barcodedSamples'])
		if isinstance(self.json_dat['plan']['barcodedSamples'],dict):
			samples = self.json_dat['plan']['barcodedSamples']
		else:
			samples = json.loads(self.json_dat['plan']['barcodedSamples'])
		bamPaths = []
		bams = []
	
		try:
			reference_path = self.envDict['TSP_FILEPATH_GENOME_FASTA']
		except:
			reference_path = ''
			
		for datum in json_basecaller['datasets']:
			if reference_path != '':
				tempPath = os.path.join(self.json_dat['runinfo']['alignment_dir'], datum['file_prefix']+'.bam')
			else:
				tempPath = os.path.join(self.json_dat['runinfo']['basecaller_dir'], datum['file_prefix']+'.basecaller.bam')

			if os.path.exists(tempPath):
				bamPaths.append(tempPath)
				if datum['dataset_name'][datum['dataset_name'].rfind('/')+1:] != 'No_barcode_match' and '/' in datum['dataset_name']:
					bams.append(datum['dataset_name'])

		# get the list of 'valid' barcodes or samples (could be either depending on whether user altered names with run planning
		# and sort of hacky, but extract this from the BAM file names we just got above
		for bamFileName in bamPaths:
			barcodeName = bamFileName.split('/')[-1] # get the last part, just the name with no path (probably can use os method here too)
			barcodeName = barcodeName.split('_rawlib')[0] # get just the barcode part of the name
			# find a possible matching sample name
			for sampleItemName in samples:
				sampleItem = samples[sampleItemName]
				if barcodeName in sampleItem['barcodes']:
					self.sampleNameLookup[barcodeName] = sampleItemName
			if barcodeName in self.sampleNameLookup.keys():
				sampleName = self.sampleNameLookup[barcodeName]
			else:
				sampleName = ''
				self.sampleNameLookup[barcodeName] = '' # makes it much easier later to do the lookup with no conditional tests
			# MGD note: I considered setting blank sample names to the barcode name instead, but might not be what customer intended
			#print 'BARCODE FOUND: %s SAMPLE ID: %s' % (barcodeName, sampleName)
			self.barcodeNames.append(barcodeName)
		
		#not to analyse BLANC_PCR sample
		for key in self.sampleNameLookup.keys():
			if self.sampleNameLookup[key] == "BLANC_PCR":
				print "To remove dict value "+key+" : "+self.sampleNameLookup[key]
				del self.sampleNameLookup[key]
                        #elif self.sampleNameLookup[key] =="NOMINE_LW305":
                        #        print "Failed sample : To remove dict value "+key+" : "+self.sampleNameLookup[key]
                        #        del self.sampleNameLookup[key]
                        #elif self.sampleNameLookup[key] =="ROY_LW650":
                        #        print "Failed sample : To remove dict value "+key+" : "+self.sampleNameLookup[key]
                        #        del self.sampleNameLookup[key]

		data = {}
		
		# WORKBOOK STYLES
		#HomoReport_path = "%s/%s.homopolymersReport.xls" % (self.json_dat['runinfo']['results_dir'],run_id)
		HomoReport_path = "%s/homopolymersReport.xlsx" % (self.json_dat['runinfo']['results_dir'])
		
		HomoReport = xlsxwriter.Workbook(HomoReport_path)

		# headerStyle_string = "font:name Calibri, height 220, bold on;borders"
		headerStyle = HomoReport.add_format({'font_name': 'Calibri', 'bold': True, 'border': 1})
		generalStyle = HomoReport.add_format({'font_name': 'Calibri', 'border': 1})
		#generalStyle = HomoReport.add_format('font: name Calibri, height 220')
		#greenBGStyle = HomoReport.add_format('font: name Calibri, height 220; pattern: pattern solid, fore-colour 0x2A')
		boldStyle = HomoReport.add_format({'font_name': 'Calibri', 'bold': True})
		#redBGStyle = HomoReport.add_format('font: name Calibri, height 220; pattern: pattern solid, fore-colour 0x1D')

	

	
		#boucle pour toutes les mutations
		conf_file = "/results/plugins/Homopolymers/conf/min6_homopolymers_positions.tsv"
		confReader = csv.reader(open(conf_file,"rb"),delimiter="\t")
		for row in confReader:
			if row[0].startswith('#'):	# pass headers
				continue
			else:
				chrom = row[0]
				positions = row[1]
				gene = row[2]
				hgvsc = row[3]
		
				#add mutSheet
				sheet_id = gene +" "+hgvsc 
				mutSheet = HomoReport.add_worksheet(sheet_id)
				line = 0
				col = 0
				mutSheet.write(line, col,sheet_id,generalStyle)
				line=line+1

				header = "Sample\tIns\tDel\tDepth\tReal depth\tFreq Ins\tFreq Del\tFlow Dup"
				header_fields = header.split('\t')
				for i in range(len(header_fields)):
					mutSheet.write(line,i,header_fields[i],headerStyle)
				#line=line+1

				header = "Sample\tFlow Ins\tFlow Del\tFlow Depth\tFlow Real depth\tFlow Freq Ins\tFlow Freq Del\tFlow Dup"
				header_fields = header.split('\t')
				for i in range(len(header_fields)):
					mutSheet.write(line,i+13,header_fields[i],headerStyle)
				line=line+1

                                print 'CHROMOSOME %s' % chrom
                                print 'POSITION %s' % positions

				for vcf in glob.glob('%s/plugin_out/variantCaller_out.%s/IonXpress_*/TSVC_variants.vcf' % (self.envDict['ANALYSIS_DIR'],variantCallerInstance)):
					barcode = vcf.split('/')[-2]
					
					
					try :
						sample = self.sampleNameLookup[barcode]
					except KeyError:
						continue
					print "Sample : " + sample
					
					
					data[sample] = {}

                                        		
					vcfReader = csv.reader(open(vcf,"rb"),delimiter="\t")
					typ = False
					for row in vcfReader:
						if row[0].startswith('#'):	# pass headers
							continue
                                                #elif "COSM34210" in row[2].split(';'):
                                                elif chrom in row[0] :
                                                        if ',' in positions:
                                                                for position in positions.split(','):
                                                                        if position in row[1]:
                                                                                if len(row[2].split(';')) > 1 :
                                                                                        print "WARNING : more than one hotspot found at dupG position."
                                                                                        print "-> Check .VCF to see if mutation type are style correct (ins =/= complex)"
                                                                                        #htmlOut.write('<br>Warning : see plugin log</b>')
                                                                                infoField = row[7]
                                                                                print "mut position found in TSVC_variants.vcf"
                                                                                print "%s : %s : %s : %s > %s" %(row[0],row[1],row[2],row[3],row[4])
                                                                                data[sample]["alt"] = row[4].split('=')[-1]
                                                                                if row[6]=="PASS":
                                                                                        data[sample]["filter"] = "No"
                                                                                        data[sample]["reason"] = ""
                                                                                else:
                                                                                        data[sample]["filter"] = "Yes"
                                                                                        print "INFO Field = " + row[7]
                                                                                taglist = infoField.split(";")
                                                                                for tag in taglist:
                                                                                        if tag.startswith("AO="):
                                                                                                ao = tag.split('=')[-1].split(',')
                                                                                        if tag.startswith("DP="):
                                                                                                dp = tag.split('=')[-1]
                                                                                                data[sample]["DP"] = dp
                                                                                        if tag.startswith("TYPE="):
                                                                                                typ = tag.split('=')[-1].split(',')
                                                                                        if tag.startswith("RBI="):
                                                                                                data[sample]["rbi"] = tag.split('=')[-1]
                                                                                        if tag.startswith("FAO="):
                                                                                                fao = tag.split('=')[-1].split(',')
                                                                                                #print "VCF FAO: "
                                                                                                #print fao
                                                                                        if tag.startswith("FDP="):
                                                                                                fdp = tag.split('=')[-1]
                                                                                                data[sample]["FDP"] = fdp
                                                                                                #print "VCF FDP: "
                                                                                                #print fdp
                                                                                        if data[sample]["filter"]=="Yes":
                                                                                                if tag.startswith("FR="):
                                                                                                        data[sample]["reason"] = tag.split('=')[-1]
                                                                                break

                                                        elif positions in row[1]:
                                                                if len(row[2].split(';')) > 1 :
                                                                        print "WARNING : more than one hotspot found at dupG position."
                                                                        print "-> Check .VCF to see if mutation type are style correct (ins =/= complex)"
                                                                        #htmlOut.write('<br>Warning : see plugin log</b>')
                                                                infoField = row[7]
                                                                print "mut position found in TSVC_variants.vcf"
                                                                print "%s : %s : %s : %s > %s" %(row[0],row[1],row[2],row[3],row[4])
                                                                data[sample]["alt"] = row[4].split('=')[-1]
                                                                if row[6]=="PASS":
                                                                        data[sample]["filter"] = "No"
                                                                        data[sample]["reason"] = ""
                                                                else:
                                                                        data[sample]["filter"] = "Yes"
                                                                        print "INFO Field = " + row[7]
                                                                taglist = infoField.split(";")
                                                                for tag in taglist:
                                                                        if tag.startswith("AO="):
                                                                                ao = tag.split('=')[-1].split(',')
                                                                        if tag.startswith("DP="):
                                                                                dp = tag.split('=')[-1]
                                                                                data[sample]["DP"] = dp
                                                                        if tag.startswith("TYPE="):
                                                                                typ = tag.split('=')[-1].split(',')
                                                                        if tag.startswith("RBI="):
                                                                                data[sample]["rbi"] = tag.split('=')[-1]
                                                                        if tag.startswith("FAO="):
                                                                                fao = tag.split('=')[-1].split(',')
                                                                                #print "VCF FAO: "
                                                                                #print fao
                                                                        if tag.startswith("FDP="):
                                                                                fdp = tag.split('=')[-1]
                                                                                data[sample]["FDP"] = fdp
                                                                                #print "VCF FDP: "
                                                                                #print fdp
                                                                        if data[sample]["filter"]=="Yes":
                                                                                if tag.startswith("FR="):
                                                                                        data[sample]["reason"] = tag.split('=')[-1]
                                        ao_ins = 0
                                        ao_del = 0
                                        fao_ins = 0
                                        fao_del = 0

                                        #mut not found in TSVC_variants.vcf look for it in small_variant_filtered.vcf
                                        if not typ:
                                                data[sample]["DP"] = 0
                                                data[sample]["FDP"] = 0
                                                print "WARNING : no dupA found in TSVC_variants.vcf for %s" % sample
                                                vcf = "%s/plugin_out/variantCaller_out.%s/%s/small_variants_filtered.vcf" % (self.envDict['ANALYSIS_DIR'],variantCallerInstance,barcode)
                                                print "vcf filtered : %s" % vcf
                                                vcfReader = csv.reader(open(vcf,"rb"),delimiter="\t")
                                                typ = False
                                                for row in vcfReader:
                                                        if row[0].startswith('#'):	# pass headers
                                                                continue
                                                                #elif "COSM34210" in row[2].split(';'):
                                                        elif chrom in row[0] :
                                                                if ',' in positions:
                                                                        for position in positions.split(','):
                                                                                if position in row[1]:
                                                                                        if len(row[2].split(';')) > 1 :
                                                                                                print "WARNING : more than one hotspot found at dupG position."
                                                                                                print "-> Check .VCF to see if mutation type are style correct (ins =/= complex)"
                                                                                                #htmlOut.write('<br>Warning : see plugin log</b>')
                                                                                        infoField = row[7]
                                                                                        print "mut position found in filtered variants"
                                                                                        print "%s : %s : %s : %s > %s" %(row[0],row[1],row[2],row[3],row[4])
                                                                                        data[sample]["alt"] = row[4].split('=')[-1]
                                                                                        data[sample]["filter"] = "Yes"
                                                                                        print "INFO Field = " + row[7]
                                                                                        taglist = infoField.split(";")
                                                                                        for tag in taglist:
                                                                                                if tag.startswith("AO="):
                                                                                                        ao = tag.split('=')[-1].split(',')
                                                                                                if tag.startswith("DP="):
                                                                                                        dp = tag.split('=')[-1]
                                                                                                        data[sample]["DP"] = dp
                                                                                                if tag.startswith("TYPE="):
                                                                                                        typ = tag.split('=')[-1].split(',')
                                                                                                if tag.startswith("RBI="):
                                                                                                        data[sample]["rbi"] = tag.split('=')[-1]
                                                                                                if tag.startswith("FAO="):
                                                                                                        fao = tag.split('=')[-1].split(',')
                                                                                                        #print "FILTERED FAO :"
                                                                                                        #print fao
                                                                                                if tag.startswith("FDP="):
                                                                                                        fdp = tag.split('=')[-1]
                                                                                                        data[sample]["FDP"] = fdp
                                                                                                        #print "FILTERED FDP: "
                                                                                                        #print fdp
                                                                                                if tag.startswith("FR="):
                                                                                                        FR = tag.split('=')[-1]
                                                                                                        data[sample]["reason"] = FR
                                                                                        break
                                                                elif positions in row[1]:
                                                                        if len(row[2].split(';')) > 1 :
                                                                                print "WARNING : more than one hotspot found at dupG position."
                                                                                print "-> Check .VCF to see if mutation type are style correct (ins =/= complex)"
                                                                                #htmlOut.write('<br>Warning : see plugin log</b>')
                                                                        infoField = row[7]
                                                                        print "mut position found in filtered variants"
                                                                        print "%s : %s : %s : %s > %s" %(row[0],row[1],row[2],row[3],row[4])
                                                                        data[sample]["alt"] = row[4].split('=')[-1]
                                                                        data[sample]["filter"] = "Yes"
                                                                        print "INFO Field = " + row[7]
                                                                        taglist = infoField.split(";")
                                                                        for tag in taglist:
                                                                                if tag.startswith("AO="):
                                                                                        ao = tag.split('=')[-1].split(',')
                                                                                if tag.startswith("DP="):
                                                                                        dp = tag.split('=')[-1]
                                                                                        data[sample]["DP"] = dp
                                                                                if tag.startswith("TYPE="):
                                                                                        typ = tag.split('=')[-1].split(',')
                                                                                if tag.startswith("RBI="):
                                                                                        data[sample]["rbi"] = tag.split('=')[-1]
                                                                                if tag.startswith("FAO="):
                                                                                        fao = tag.split('=')[-1].split(',')
                                                                                        #print "FILTERED FAO :"
                                                                                        #print fao
                                                                                if tag.startswith("FDP="):
                                                                                        fdp = tag.split('=')[-1]
                                                                                        data[sample]["FDP"] = fdp
                                                                                        #print "FILTERED FDP: "
                                                                                        #print fdp
                                                                                if tag.startswith("FR="):
                                                                                        FR = tag.split('=')[-1]
                                                                                        data[sample]["reason"] = FR



                                                if not typ:
                                                        data[sample]["DP"] = 0
                                                        data[sample]["AO_ins"] = 0
                                                        data[sample]["AO_del"] = 0
                                                        data[sample]["FDP"] = 0
                                                        data[sample]["FAO_ins"] = 0
                                                        data[sample]["FAO_del"] = 0
                                                        data[sample]["filter"] = "Not candidate"
                                                        data[sample]["reason"] = "Not candidate"
                                                        print "WARNING : no dupA found neither in TSVC_variants.vcf neither in small_variant_filtered.vcf for %s. Default values 0 for ins and del." % sample
                                                        #htmlOut.write('<br>Warning : see plugin log</b>')
                                                else:
                                                        for i in range (len(typ)):
                                                                if typ[i] == "ins":	#or typ[i] == "complex" ?
                                                                        ao_ins = ao_ins + int(ao[i])
                                                                        fao_ins = fao_ins + int(fao[i])
                                                                elif typ[i] == "del":
                                                                        ao_del = ao_del + int(ao[i])
                                                                        fao_del = fao_del + int(fao[i])
                                                        data[sample]["AO_ins"] = ao_ins
                                                        data[sample]["AO_del"] = ao_del
                                                        data[sample]["FAO_ins"] = fao_ins
                                                        data[sample]["FAO_del"] = fao_del

                                        #dup found in TSVC_variant.vcf
                                        else:
                                                for i in range (len(typ)):
                                                        if typ[i] == "ins":	#or typ[i] == "complex" ?
                                                                ao_ins = ao_ins + int(ao[i])
                                                                fao_ins = fao_ins + int(fao[i])
                                                        elif typ[i] == "del":
                                                                ao_del = ao_del + int(ao[i])
                                                                fao_del = fao_del + int(fao[i])
                                                data[sample]["AO_ins"] = ao_ins
                                                data[sample]["AO_del"] = ao_del
                                                data[sample]["FAO_ins"] = fao_ins
                                                data[sample]["FAO_del"] = fao_del
                                                #data[sample]["filter"] = "No"
                                                #data[sample]["reason"] = ""
			
                                #variant found in vcf or small_variant_filtered_vcf but not allele type "del" or "ins" defined, set to "0" default value
                                if not ao_del :
                                        ao_del = 0
                                        data[sample]["AO_del"] = ao_del
                                if not ao_ins :
                                        ao_ins = 0
                                        data[sample]["AO_ins"] = ao_ins
                                if not fao_del :
                                        fao_del = 0
                                        data[sample]["FAO_del"] = fao_del
                                if not fao_ins :
                                        fao_ins = 0
                                        data[sample]["FAO_ins"] = fao_ins
				### CALCUL DE LA DETECTION DE DUP ###
				#list_ao_ins = []
				#list_ao_del = []
				#list_dp = []
				#list_dp_real = []
				list_freq_ins = []
				list_freq_fao_ins = []
				#list_freq_del = []
		
				sample_list = self.sampleNameLookup.values()
				if '' in sample_list:
					sample_list.remove('')
				for sample in sample_list:
					#list_ao_ins.append(data[sample]["AO_ins"])
					#list_ao_del.append(data[sample]["AO_del"])
					#list_dp.append(data[sample]["DP"])
					
					#raw data
					dp_real = data[sample]["AO_del"] + int(data[sample]["DP"])
					data[sample]["DP_real"] = dp_real
					
					#list_dp_real.append(dp_real)
					
					if dp_real is 0 :
						freq_ins = 0
						data[sample]["FREQ_ins"] = freq_ins
						freq_del = 0
						data[sample]["FREQ_del"] = freq_del
					else:
						freq_ins = numpy.divide(float(data[sample]["AO_ins"]),float(dp_real))*100
						data[sample]["FREQ_ins"] = freq_ins
						list_freq_ins.append(freq_ins)
						freq_del = numpy.divide(float(data[sample]["AO_del"]),float(dp_real))*100
						data[sample]["FREQ_del"] = freq_del
						#list_freq_del.append(freq_del)
					
					#recalibrated data
					fdp_real = data[sample]["FAO_del"] + int(data[sample]["FDP"])
					data[sample]["FDP_real"] = fdp_real
					#print "FDP REAL: "
					#print fdp_real
					
					if fdp_real is 0 :
						freq_Fins = 0
						data[sample]["FREQ_FAO_ins"] = freq_Fins
						freq_Fdel = 0
						data[sample]["FREQ_FAO_del"] = freq_Fdel
					else:
						freq_fao_ins = numpy.divide(float(data[sample]["FAO_ins"]),float(fdp_real))*100
						data[sample]["FREQ_FAO_ins"] = freq_fao_ins
						list_freq_fao_ins.append(freq_fao_ins)
						freq_fao_del = numpy.divide(float(data[sample]["FAO_del"]),float(fdp_real))*100
						data[sample]["FREQ_FAO_del"] = freq_fao_del
						#list_freq_del.append(freq_del)


				#raw data
				moy_freq_ins = numpy.mean(list_freq_ins)
				ecart_type_freq_ins = numpy.std(list_freq_ins,ddof=1)
				threshold = moy_freq_ins + 2*ecart_type_freq_ins

				#recalibrated data
				moy_freq_fao_ins = numpy.mean(list_freq_fao_ins)
				ecart_type_freq_fao_ins = numpy.std(list_freq_fao_ins,ddof=1)
				fao_threshold = moy_freq_fao_ins + 2*ecart_type_freq_fao_ins
							
				for sample in sample_list:
					#raw data
					if data[sample]["FREQ_ins"] > threshold:
						data[sample]["Dup"] = "Probable"
					else:
						data[sample]["Dup"] = "Non"
					#recalibrated data
					if data[sample]["FREQ_FAO_ins"] > fao_threshold:
						data[sample]["FAO_Dup"] = "Probable"
					else:
						data[sample]["FAO_Dup"] = "Non"
		
				filename = "dup.csv"
				csvfile = open(result_dir+"/"+filename,"wb")
				results = csv.writer(csvfile, delimiter =";")

				fao_filename = "recalibrated_dup.csv"
				fao_csvfile = open(result_dir+"/"+fao_filename,"wb")
				fao_results = csv.writer(fao_csvfile, delimiter =";")
				#header
				results.writerow(['Dup','mean = %s' % numpy.around(moy_freq_ins,1),'standard_deviation = %s' % numpy.around(ecart_type_freq_ins,1), 'threshold = %s' % numpy.around(threshold,1)])
				results.writerow(['Sample','Ins','Del','Depth','Real depth','Freq Ins','Freq Del','Dup'])

				fao_results.writerow(['FAO_Dup','FAO_mean = %s' % numpy.around(moy_freq_fao_ins,1),'standard_deviation = %s' % numpy.around(ecart_type_freq_fao_ins,1), 'threshold = %s' % numpy.around(fao_threshold,1)])
				fao_results.writerow(['Sample','FAO_Ins','FAO_Del','FDepth','FReal depth','Freq FAO Ins','Freq FAO Del','FAO_Dup'])
				for sample in sample_list:	
					results.writerow([sample,data[sample]["AO_ins"],data[sample]["AO_del"],data[sample]["DP"],data[sample]["DP_real"],numpy.around(data[sample]["FREQ_ins"],1),numpy.around(data[sample]["FREQ_del"],1),data[sample]["Dup"]])
					#write xls file
					#raw data
					mutSheet.write(line,0,sample,generalStyle)
					mutSheet.write(line,1,data[sample]["AO_ins"],generalStyle)
					mutSheet.write(line,2,data[sample]["AO_del"],generalStyle)
					mutSheet.write(line,3,data[sample]["DP"],generalStyle)
					mutSheet.write(line,4,data[sample]["DP_real"],generalStyle)
					mutSheet.write(line,5,numpy.around(data[sample]["FREQ_ins"],1),generalStyle)
					mutSheet.write(line,6,numpy.around(data[sample]["FREQ_del"],1),generalStyle)
					mutSheet.write(line,7,data[sample]["Dup"],generalStyle)
					#recalibrated data
					mutSheet.write(line,13,sample,generalStyle)
					mutSheet.write(line,14,data[sample]["FAO_ins"],generalStyle)
					mutSheet.write(line,15,data[sample]["FAO_del"],generalStyle)
					mutSheet.write(line,16,data[sample]["FDP"],generalStyle)
					mutSheet.write(line,17,data[sample]["FDP_real"],generalStyle)
					mutSheet.write(line,18,numpy.around(data[sample]["FREQ_FAO_ins"],1),generalStyle)
					mutSheet.write(line,19,numpy.around(data[sample]["FREQ_FAO_del"],1),generalStyle)
					mutSheet.write(line,20,data[sample]["FAO_Dup"],generalStyle)
					line=line+1
					
				csvfile.close()
				fao_csvfile.close()
		
				### RAW DATA GRAPH ###
				fig = figure()
		
				N = len(sample_list)
				ind = numpy.arange(N)
				width = 0.2
				
				tab_freq_ins = []
				tab_freq_del = []
				for sample in sample_list:
					tab_freq_ins.append(data[sample]["FREQ_ins"])
					tab_freq_del.append(data[sample]["FREQ_del"])
			
				ax = fig.add_subplot(111)
				rects1 = ax.bar(ind, tab_freq_ins, width, facecolor='#9999ff', label='d')
				rects2 = bar(ind+width, tab_freq_del, width, facecolor='#ff9999')
				ax.set_xlabel('Patients')
				ax.set_ylabel('Frequence')
		
				xticks(ind+width,sample_list,fontsize=6)
				yticks(fontsize=6)

				ax.legend( (rects1[0], rects2[0]), ('freq ins', 'freq del'))
		
				# horizontal line indicating the threshold
				ax.plot([0., N], [threshold, threshold], "k--")
		
				# ax.set_xlim([0,max(tab_freq_del)])
				#max_y = max([max(tab_freq_ins),max(tab_freq_del)])
				max_y = 100
				ax.set_ylim([0,max_y])	#augmentation un peu de l'ordonne, trop courte

				fig.autofmt_xdate()
		
				fig_name = gene +"_"+hgvsc+".png"
				fig.savefig(fig_name,figsize=(10,50),dpi=200)

				### RECALIBRATED GRAPH ###
				fig = figure()
		
				N = len(sample_list)
				ind = numpy.arange(N)
				width = 0.2
				
				tab_freq_ins = []
				tab_freq_del = []
				for sample in sample_list:
					tab_freq_ins.append(data[sample]["FREQ_FAO_ins"])
					tab_freq_del.append(data[sample]["FREQ_FAO_del"])
			
				ax = fig.add_subplot(111)
				rects1 = ax.bar(ind, tab_freq_ins, width, facecolor='#9999ff', label='d')
				rects2 = bar(ind+width, tab_freq_del, width, facecolor='#ff9999')
				ax.set_xlabel('Patients')
				ax.set_ylabel('Frequence')
		
				xticks(ind+width,sample_list,fontsize=6)
				yticks(fontsize=6)

				ax.legend( (rects1[0], rects2[0]), ('freq ins', 'freq del'))
		
				# horizontal line indicating the threshold
				ax.plot([0., N], [fao_threshold, fao_threshold], "k--")
		
				# ax.set_xlim([0,max(tab_freq_del)])
				#max_y = max([max(tab_freq_ins),max(tab_freq_del)])
				max_y = 100
				ax.set_ylim([0,max_y])	#augmentation un peu de l'ordonne, trop courte

				fig.autofmt_xdate()
		
				fig_name_flow = gene +"_"+hgvsc+"flowdata.png"
				fig.savefig(fig_name_flow,figsize=(50,10),dpi=200)


				
				reader = csv.reader(open("%s/dup.csv" % self.envDict['RESULTS_DIR'],"r"), delimiter =";")
				rownum = 0
				firstline = reader.next()
				moy = firstline[1]
				ec = firstline[2]
				thr = firstline[3]

				freader = csv.reader(open("%s/recalibrated_dup.csv" % self.envDict['RESULTS_DIR'],"r"), delimiter =";")
				frownum = 0
				firstline = freader.next()
				fmoy = firstline[1]
				fec = firstline[2]
				fthr = firstline[3]

				
				# add fig in xls sheet
				line = line+2
				#mutSheet.insert_image(line,0,fig_name,{'x_scale':0.5, 'y_scale':0.5})
                                mutSheet.insert_image(line,0,fig_name)
				#mutSheet.insert_image(line,13,fig_name_flow,{'x_scale':0.5, 'y_scale':0.5})
                                mutSheet.insert_image(line,13,fig_name_flow)
		

				#add stats data to xls sheet
				mutSheet.write(1,9,moy,boldStyle)
				mutSheet.write(2,9,ec,boldStyle)
				mutSheet.write(3,9,thr,boldStyle)

				mutSheet.write(1,22,fmoy,boldStyle)
				mutSheet.write(2,22,fec,boldStyle)
				mutSheet.write(3,22,fthr,boldStyle)
						
		
				#alt, rbi and filter info to xls sheet
				line = line+30
				mutSheet.write(line,0,"Sample",headerStyle)
				mutSheet.write(line,1,"ALT",headerStyle)
				mutSheet.write(line,2,"RBI",headerStyle)
				mutSheet.write(line,3,"Filter",headerStyle)
				mutSheet.write(line,4,"Reason",headerStyle)
				line = line+1
				for sample in sample_list:	
					#write xls file
					if "alt" in data[sample].keys():
						alt = '%s' % data[sample]["alt"]
					else:
						alt = "NA"
						data[sample]["alt"] = "NA";
					if "rbi" in data[sample].keys():
						rbi = '%s' % data[sample]["rbi"]
					else:
						rbi = "NA"
						data[sample]["rbi"] = "NA"
					#print rbi
					mutSheet.write(line,0,sample,generalStyle)
					mutSheet.write(line,1,data[sample]["alt"],generalStyle)
					mutSheet.write(line,2,data[sample]["rbi"],generalStyle)
					mutSheet.write(line,3,data[sample]["filter"],generalStyle)
					mutSheet.write(line,4,data[sample]["reason"],generalStyle)
					line=line+1
				
				#print "CHROMOSOME %s" % chrom
				#print "POSITION %s" % position
				if chrom == "chr13" and position == "32910661":
					print "dupA temoin KQ99"
					### HTML ###
					htmlOut.write('<table width="80%" border="1" align="left">')
					for row in reader :
						#write header row. assumes first row in csv contains header	
						if rownum == 0 :
							htmlOut.write('<tr>') # write <tr> tag
							for column in row:
								htmlOut.write('<th>' + column + '</th>')
							htmlOut.write('</tr>')
						#write all other rows   
						else:
							htmlOut.write('<tr>')    
							for column in row:
								htmlOut.write('<td>' + column + '</td>')
							htmlOut.write('</tr>')
			     
						#increment row count    
						rownum += 1
					htmlOut.write('</table>')
					htmlOut.write('<br><b>%s</b>' % moy)
					htmlOut.write('<br><b>%s</b>' % ec)
					htmlOut.write('<br><b>%s</b>' % thr)
					htmlOut.write('<br><img src="BRCA2_c.2175dupA.png" style="width: 50%;max-height: 50%" align="left" />')
				
		homopolymersReport = "homopolymersReport.xlsx"
		htmlOut.write('''<li><br/><a href="%s">%s</a><br></li> '''  % (homopolymersReport, "Telecharger le rapport d'homopomlymeres"))
							
		htmlOut.close()
		HomoReport.close()
		return True

'''		
		###############################
		####### HotCount script #######
		###############################
		
		# generer les fastq necessaires pour le script
		FastqCreatorCp = Popen(['cp', '%s/FastqCreator.py' % self.envDict['DIRNAME'], self.envDict['RESULTS_DIR']], stdout=PIPE)
		out, err = FastqCreatorCp.communicate()
		print 'exec: FastqCreator.py'
		FastqCmd = Popen(['python', '%s/FastqCreator.py' % self.envDict['RESULTS_DIR']], stdout=PIPE)
		FastqOut, FastqErr = FastqCmd.communicate()
			
		# lancement script
		os.chdir(self.envDict['RESULTS_DIR'])
		#HotCountCmd1 = Popen(["bash","%s/HotCount/do_it.sh" % self.envDict["DIRNAME"],"%s/HotCount/asxl1.txt" % self.envDict["DIRNAME"],"%s/*.fastq" % self.envDict["RESULTS_DIR"],">","%s/comptage-asxl1.tsv" % self.envDict["RESULTS_DIR"]], stdout=PIPE, shell=True)
		HotCountCmd1 = Popen("bash %s/HotCount/do_it.sh %s/HotCount/asxl1.txt *.fastq > comptage-asxl1.tsv" % (self.envDict["DIRNAME"],self.envDict["DIRNAME"]), stdout=PIPE, shell=True)
		HotCountCmd1Out, HotCountCmd1Err = HotCountCmd1.communicate()
		
		#HotCountCmd2 = Popen(["Rscript", "%s/HotCount/do_it.R" % self.envDict["DIRNAME"], "%s/comptage-asxl1.tsv" % self.envDict["RESULTS_DIR"],"ALL","dupG","delG"], stdout=PIPE, shell=True)
		HotCountCmd2 = Popen("Rscript %s/HotCount/do_it.R comptage-asxl1.tsv ALL dupG delG > HotCount_results.txt" % self.envDict["DIRNAME"], stdout=PIPE, shell=True)
		HotCountCmd2Out, HotCountCmd2Err = HotCountCmd2.communicate()
		
		# nettoyer les fastq
		FastqDel = Popen("rm *.fastq", stdout=PIPE, shell=True)
		FastqDelOut, FastqDelErr = FastqDel.communicate()

		comptage = open("comptage-asxl1.tsv")
		comptagelines = comptage.readlines()
		txt = open("HotCount_results.txt")
		txtlines = txt.readlines()
		htmlOut.write('<BR CLEAR="all"><hr>')
		htmlOut.write('<b> comptage-asxl1.tsv :</b>')
		for line in comptagelines:
			htmlOut.write('<br>%s' % line)
		htmlOut.write('<br><hr>')
		htmlOut.write('<b> resultats :</b>')
		for line in txtlines:
			htmlOut.write('<br>%s' % line)
'''
		#htmlOut.close()
		#return True


if __name__ == "__main__":
    PluginCLI(Homopolymers())
