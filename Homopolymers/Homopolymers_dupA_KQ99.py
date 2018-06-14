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
			htmlOut = open('%s/plugin_out/Homopolymers_out/Homopolymrs_block.html'%self.envDict['ANALYSIS_DIR'], 'a+')
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
		
		data = {}
		for vcf in glob.glob('%s/plugin_out/variantCaller_out.%s/IonXpress_*/TSVC_variants.vcf' % (self.envDict['ANALYSIS_DIR'],variantCallerInstance)):
			barcode = vcf.split('/')[-2]
			sample = self.sampleNameLookup[barcode]
			print "Sample : " + sample
			data[sample] = {}
		
			vcfReader = csv.reader(open(vcf,"rb"),delimiter="\t")
			typ = False
			for row in vcfReader:
				if row[0].startswith('#'):	# pass headers
					continue
				#elif "COSM34210" in row[2].split(';'):
				elif "chr13" in row[0] and "32910661" in row[1]:
					if len(row[2].split(';')) > 1 :
						print "WARNING : more than one hotspot found at dupG position."
						print "-> Check .VCF to see if mutation type are style correct (ins =/= complex)"
						htmlOut.write('<br>Warning : see plugin log</b>')
					infoField = row[7]
					print "dupA position found"
					print "%s : %s : %s : %s > %s" %(row[0],row[1],row[2],row[3],row[4])
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
					
			ao_ins = 0
			ao_del = 0

			#dupA not found in TSVC_variants.vcf look for it in small_variant_filtered.vcf
			if not typ:
				data[sample]["DP"] = 0
				print "WARNING : no dupA found in TSVC_variants.vcf for %s" % sample
				vcf = "%s/plugin_out/variantCaller_out.%s/%s/small_variants_filtered.vcf" % (self.envDict['ANALYSIS_DIR'],variantCallerInstance,barcode)
				print "vcf filtered : %s" % vcf
				vcfReader = csv.reader(open(vcf,"rb"),delimiter="\t")
				typ = False
				for row in vcfReader:
					if row[0].startswith('#'):	# pass headers
						continue
				#elif "COSM34210" in row[2].split(';'):
					elif "chr13" in row[0] and "32910661" in row[1]:
						if len(row[2].split(';')) > 1 :
							print "WARNING : more than one hotspot found at dupG position."
							print "-> Check .VCF to see if mutation type are style correct (ins =/= complex)"
							htmlOut.write('<br>Warning : see plugin log</b>')
						infoField = row[7]
						print "dupA position found"
						print "%s : %s : %s : %s > %s" %(row[0],row[1],row[2],row[3],row[4])
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
				if not typ:
					data[sample]["DP"] = 0
					print "WARNING : no dupA found neither in TSVC_variants.vcf neither in small_variant_filtered.vcf for %s" % sample
					htmlOut.write('<br>Warning : see plugin log</b>')
				else:
					for i in range (len(typ)):
						if typ[i] == "ins":	#or typ[i] == "complex" ?
							ao_ins = ao_ins + int(ao[i])
						elif typ[i] == "del":
							ao_del = ao_del + int(ao[i])

			#dupA found in TSVC_variant.vcf
			else:
				for i in range (len(typ)):
					if typ[i] == "ins":	#or typ[i] == "complex" ?
						ao_ins = ao_ins + int(ao[i])
					elif typ[i] == "del":
						ao_del = ao_del + int(ao[i])
			data[sample]["AO_ins"] = ao_ins
			data[sample]["AO_del"] = ao_del
			
		### CALCUL DE LA DETECTION DE DUPA ###
		#list_ao_ins = []
		#list_ao_del = []
		#list_dp = []
		#list_dp_real = []
		list_freq_ins = []
		#list_freq_del = []
		
		sample_list = self.sampleNameLookup.values()
		if '' in sample_list:
			sample_list.remove('')
		for sample in sample_list:
			#list_ao_ins.append(data[sample]["AO_ins"])
			#list_ao_del.append(data[sample]["AO_del"])
			#list_dp.append(data[sample]["DP"])
			
			dp_real = data[sample]["AO_del"] + int(data[sample]["DP"])
			data[sample]["DP_real"] = dp_real
			#list_dp_real.append(dp_real)
			
			freq_ins = numpy.divide(float(data[sample]["AO_ins"]),float(dp_real))*100
			data[sample]["FREQ_ins"] = freq_ins
			list_freq_ins.append(freq_ins)
			freq_del = numpy.divide(float(data[sample]["AO_del"]),float(dp_real))*100
			data[sample]["FREQ_del"] = freq_del
			#list_freq_del.append(freq_del)
			
		moy_freq_ins = numpy.mean(list_freq_ins)
		ecart_type_freq_ins = numpy.std(list_freq_ins,ddof=1)
		threshold = moy_freq_ins + 2*ecart_type_freq_ins
		
		for sample in sample_list:
			if data[sample]["FREQ_ins"] > threshold:
				data[sample]["DupA"] = "Probable"
			else:
				data[sample]["DupA"] = "Non"
		
		filename = "dupA.csv"
		csvfile = open(result_dir+"/"+filename,"wb")
		results = csv.writer(csvfile, delimiter =";")
		#header
		results.writerow(['DupA','mean = %s' % numpy.around(moy_freq_ins,1),'standard_deviation = %s' % numpy.around(ecart_type_freq_ins,1), 'threshold = %s' % numpy.around(threshold,1)])
		results.writerow(['Sample','Ins','Del','Depth','Real depth','Freq Ins','Freq Del','DupG'])
		for sample in sample_list:	
			results.writerow([sample,data[sample]["AO_ins"],data[sample]["AO_del"],data[sample]["DP"],data[sample]["DP_real"],numpy.around(data[sample]["FREQ_ins"],1),numpy.around(data[sample]["FREQ_del"],1),data[sample]["DupA"]])
		csvfile.close()	
		
		### GRAPH ###
		fig = figure()
		
		N = len(sample_list)
		ind = numpy.arange(N)
		width = 0.4
		
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
		
		xticks(ind+width,sample_list)
 
		ax.legend( (rects1[0], rects2[0]), ('freq ins', 'freq del') )
		
		# horizontal line indicating the threshold
		ax.plot([0., N], [threshold, threshold], "k--")
		
		#ax.set_xlim([0,max(tab_freq_del)])
		max_y = max([max(tab_freq_ins),max(tab_freq_del)])
		ax.set_ylim([0,max_y*1.5])	#augmentation un peu de l'ordonne, trop courte

		fig.autofmt_xdate()
		
		fig.savefig("freq_ins_del_histogram.png",dpi=200)

		### HTML ###
		reader = csv.reader(open("%s/dupA.csv" % self.envDict['RESULTS_DIR'],"r"), delimiter =";")
		rownum = 0
		firstline = reader.next()
		moy = firstline[1]
		ec = firstline[2]
		thr = firstline[3]
		htmlOut.write('<table width="80%" border="1" align="left">')
		for row in reader :
			# write header row. assumes first row in csv contains header	
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
		htmlOut.write('<br><img src="freq_ins_del_histogram.png" style="width: 50%;max-height: 50%" align="left" />')

		htmlOut.close()
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
