#!/usr/bin/env python
import os
import json
import csv
import re
import pysam
import operator
import urllib2
import glob
from ion.plugin import *
from subprocess import *


class Check_contamination(IonPlugin):
    """Plugin object to check contamination on control sample"""
    version = "0.1"
    envDict = dict(os.environ)
    json_dat = {}
    amplicons = {}
    barcodeNames = []
    sampleNameLookup = {} # dictionary allowing us to get the sample name associated with a particular barcode name
    

    def launch(self, data=None):
        try:
		with open('startplugin.json', 'r') as fh:
			self.json_dat = json.load(fh)
	except:
		print 'Error reading plugin json.'

	try:
		htmlOut = open('%s/plugin_out/Check_contamination_out/Check_contamination_block.html'%self.envDict['ANALYSIS_DIR'], 'a+')
	except:
		htmlOut = open('Check_contamination_block.html', 'w')
	htmlOut.write('<html><body>\n')
	
	if os.path.isfile(self.envDict['TSP_FILEPATH_BARCODE_TXT']):
		self.isBarcodedRun = True
			
	if not self.isBarcodedRun:
		print "This plugin only work for barcoded run for now. \n Exiting..."
		sys.exit(0)
		
	# set defaults for user options
	targetbed = False
	sample = False
	barcode = False
	read_num = False
	read_len = False
	
	# Parse pluginconfig json.
	try:
		targetbed = self.json_dat['pluginconfig']['targetbed']
		sample = self.json_dat['pluginconfig']['sample']
		barcode = self.json_dat['pluginconfig']['barcode']
		read_num = int(self.json_dat['pluginconfig']['read_num'])
		read_len = int(self.json_dat['pluginconfig']['read_len'])
		
	except:
		print 'Warning: plugin does not appear to be configured'
	
	htmlOut.write('<b> Settings :</b> <br/>\n')
	htmlOut.write('-Sample : %s<br/>\n' % sample)
	htmlOut.write('-Barcode : %s<br/>\n' % barcode)
	htmlOut.write('-Read length : %s<br/>\n' % read_len)
	htmlOut.write('-Read number : %s<br/>\n' % read_num)
	
	# ouverture du bam a lire
	bampath = os.environ['ANALYSIS_DIR'] + "/" + str(barcode) + "_rawlib.bam"
	bamfile = pysam.AlignmentFile(bampath,"rb")
	sample_name = sample.replace(" ","_")
	output_name = "%s_filtered.bam" % sample_name
	bamfile_filtered = pysam.AlignmentFile(output_name, "wb", template=bamfile)
	
	# ouverture du target bed, recuperation des infos de chaque amplicon
	bedfile = open(targetbed,"r")
	target_reader = csv.reader(bedfile, delimiter = '\t')
	target_reader.next()
	for row in target_reader:
		s = row[7].split(';')
		s[0] = s[0].split('=')[-1]
		s[1] = s[1].split('=')[-1]
		self.amplicons[row[3]] = {}
		self.amplicons[row[3]]['chr'] = row[0]
		self.amplicons[row[3]]['startpos'] = row[1]
		self.amplicons[row[3]]['endpos'] = row[2]
		self.amplicons[row[3]]['gene_id'] = s[0]
		self.amplicons[row[3]]['transcrit'] = s[1]
		self.amplicons[row[3]]['contamination count'] = 0
		self.amplicons[row[3]]['reads_to_analyze'] = []
	bedfile.close()

	# parcours des reads du bam, comptage de chaque read resultant d'une possible contamination
	
	#for amplicon in self.amplicons.keys():
	#	for read in bamfile.fetch(self.amplicons[amplicon]['chr'],int(self.amplicons[amplicon]['startpos']),int(self.amplicons[amplicon]['endpos'])):
	#		###if (len(read.seq) >= read_len):	read.seq comprend aussi les soft-clipped, ce n'est pas bon
	#		if (len(read.query_alignment_sequence) >= read_len): 	#aligned portion of the read, exclude soft-clipped bases
	#			self.amplicons[amplicon]['contamination count'] = self.amplicons[amplicon]['contamination count'] + 1
	#			bamfile_filtered.write(read)
	
	# Creation bam avec reads > len
	for read in bamfile.fetch():
		if (len(read.query_alignment_sequence) >= read_len): 	#aligned portion of the read, exclude soft-clipped bases
			bamfile_filtered.write(read)
	
	# Remise en odre du bam et creation du BAI
	print "samtools sort..."
	cmd = Popen([["samtools","sort","-o",output_name,"%s_filtered.bam" % sample_name], stdout=PIPE)
	out, err = cmd.communicate()
	print 'OUT: %s\nERR: %s'%(out, err)
	
	print "samtools index..."
	cmd = Popen(["samtools","index","%s_filtered.bam" % sample_name], stdout=PIPE)
	out, err = cmd.communicate()
	print 'OUT: %s\nERR: %s'%(out, err)
	
	# API recuperation parametres
	
	print "Gathering params..."
	ref = self.envDict['TSP_FILEPATH_GENOME_FASTA']
	target_unmerged = self.envDict['PLAN__BEDFILE']
	target_merge = False
	hot_unmerged = False
	hot_vcf = False
	param = False
	
	try:
		api_url = self.json_dat['runinfo']['api_url'] + '/v1/pluginresult/?format=json&plugin__name=%s&result=%s' % ("variantCaller",str(self.json_dat['runinfo']['pk']))
		f = urllib2.urlopen(api_url)
		d = json.loads(f.read())
		for plugin in d['objects']: # they are already sorted from newest to oldest, so take the first of this list which is 'completed'
			if plugin['state'] == 'Completed':
				param = plugin['store']['files'][0]['server_path']
				print "param found"
				target_merge = plugin['config']['meta']['targetregions_merge']
				print "target_merge found"
				hot_unmerged = plugin['config']['meta']['targetloci']
				print "hot_unmerged found"
				break
	except:
		print 'ERROR!  Failed to get vc parameters'
	print "Gathering params... done"
	
	# Calcule nouvelle couverture par amplicon (script plugin coverageAnalysis)
	print "amplicon coverage for filtered bam..."
	cmd = Popen(["/results/plugins/coverageAnalysis/run_coverage_analysis.sh","-B",target_unmerged,ref,"%s_filtered.bam" % sample_name], stdout=PIPE)
	out, err = cmd.communicate()
	print 'OUT: %s\nERR: %s'%(out, err)
	
	#####################################################################
	### LANCEMENT VARIANTCALLER POUR LE BAM CONTENANT LES READS > LEN ###
	#####################################################################

	print "preparing hotspot"
	
	#call("mkdir","variantCalling/")
	os.makedirs("variantCalling")
	
	cmd = Popen(["tvcutils","prepare_hotspots","-b",hot_unmerged,"-o","variantCalling/hotspot_prepared.vcf","-r",ref,"-a"], stdout=PIPE)
	out, err = cmd.communicate()
	print 'OUT: %s\nERR: %s'%(out, err)
	
	
	print "Launching VariantCaller"
	cmd = Popen(["python","/results/plugins/variantCaller/variant_caller_pipeline.py","-i","%s_filtered.bam" % sample_name,"-r",ref,"-b",target_merge,"-s","variantCalling/hotspot_prepared.vcf","-p",param,"-o","variantCalling"], stdout=PIPE)
	out, err = cmd.communicate()
	print 'OUT: %s\nERR: %s'%(out, err)
	
	print "Generate variants table..."
	cmd = Popen(["python","/results/plugins/variantCaller/scripts/generate_variant_tables.py","--input-vcf","variantCalling/TSVC_variants.vcf","--region-bed",target_unmerged,"--hotspots","--output-xls","variantCalling/variants.xls","--alleles2-xls","variantCalling/alleles.xls","--summary-json","variantCalling/variant_summary.json"], stdout=PIPE)
	out, err = cmd.communicate()
	#,"--scatter-png","scatter.png",
	print 'OUT: %s\nERR: %s'%(out, err)
	
	bamfile.close()
	bamfile_filtered.close()
		
	#################################################################
	### RAJOUT COLONNES COUVERTURE DE CHAQUE AMPLICON PAR PATIENT ###
	#################################################################
	
	try:
		reference_path = self.envDict['TSP_FILEPATH_GENOME_FASTA']
	except:
		reference_path = ''
	samples = json.loads(self.json_dat['plan']['barcodedSamples'])
	# Get bam filenames.
	with open(os.path.join(self.json_dat['runinfo']['basecaller_dir'], 'datasets_basecaller.json'), 'r') as f:
		json_basecaller = json.load(f)
	bamPaths = []
	bams = []
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
		self.barcodeNames.append(barcodeName)
	self.sampleNameLookup[''] = '' # allows us to easily handle case where barcode might not have been found
	
	coverageAnalysisPath = False
	try:
		api_url = self.json_dat['runinfo']['api_url'] + '/v1/pluginresult/?format=json&plugin__name=%s&result=%s' % ("coverageAnalysis",str(self.json_dat['runinfo']['pk']))
		f = urllib2.urlopen(api_url)
		d = json.loads(f.read())
		for plugin in d['objects']: # they are already sorted from newest to oldest, so take the first of this list which is 'completed'
			if plugin['state'] == 'Completed':
				coverageAnalysisPath = plugin['path']
				break
	except:
		print 'ERROR!  Failed to get variant caller path'
			
	if coverageAnalysisPath:
		for covFile in glob.glob('%s/plugin_out/%s/IonXpress_*/*.amplicon.cov.xls' % (self.envDict['ANALYSIS_DIR'],coverageAnalysisPath.split('/')[-1])):
			barcodeName = covFile.split("/")[-2]
			sampleName = self.sampleNameLookup[barcodeName]
			# ouverture du fichier de couverture par amplicon de chaque patient, recuperation des infos
			coveragefile = open(covFile,"r")
			reader = csv.reader(coveragefile, delimiter = '\t')
			reader.next()
			for row in reader:
				amplicon = row[3]
				totalReads = row[9]
				self.amplicons[amplicon][sampleName] = totalReads

	else:
		print 'WARNING! No completed instance found for coverageAnalysis plugin. Coverage for each sample will not be displayed in final result'
	
	##############################
	### ecriture des resultats ###
	##############################
	
	# creation fichier resultat		
	#self.OUTPUT_DIR = os.path.join(self.envDict['ANALYSIS_DIR'], 'plugin_out', 'downloads')
	#if not os.path.isdir(self.OUTPUT_DIR):
	#	cmd = Popen(['mkdir', self.OUTPUT_DIR], stdout=PIPE, env=self.envDict)
	#	out, err = cmd.communicate()
	#	print 'OUT: %s\nERR: %s'%(out, err)
	
	sample_tested = self.sampleNameLookup[barcode]
	filename = "Check_contamination_%s.csv" % sample_tested
	#csvfile = open(self.OUTPUT_DIR+"/"+filename,"wb")
	csvfile = open(filename,"wb")
	results = csv.writer(csvfile, delimiter =";")
	
	htmlOut.write('<br/><b> RESULTS :</b> <br/>\n')
	samples_list = self.sampleNameLookup.values()
	try:
		samples_list.remove('')
	except:
		pass
	header = ["Amplicon","Chr","Start","End","Gene_id", "Transcript", "Number of reads > %s nucleotides" % read_len, "Contamination"]
	for sample in samples_list:
		header.append("Total_Reads_"+sample) 
	results.writerow(header)
	
	contam = False
	for amplicon in self.amplicons.keys():
		if self.amplicons[amplicon]['contamination count'] >= read_num :
			contam = True
			htmlOut.write('Possible contamination for this amplicon : %s (%s, %s, %s, %s, %s) <br/>\n' % (amplicon,self.amplicons[amplicon]['chr'],self.amplicons[amplicon]['startpos'],self.amplicons[amplicon]['endpos'],self.amplicons[amplicon]['gene_id'],self.amplicons[amplicon]['transcrit']))
			htmlOut.write('&#8627; %s reads found over %s nucletoides <br/>' % (self.amplicons[amplicon]['contamination count'], read_len))
			result_line = [amplicon,self.amplicons[amplicon]['chr'],self.amplicons[amplicon]['startpos'],self.amplicons[amplicon]['endpos'],self.amplicons[amplicon]['gene_id'],self.amplicons[amplicon]['transcrit'],self.amplicons[amplicon]['contamination count'],"Yes"]
		else:
			result_line = [amplicon,self.amplicons[amplicon]['chr'],self.amplicons[amplicon]['startpos'],self.amplicons[amplicon]['endpos'],self.amplicons[amplicon]['gene_id'],self.amplicons[amplicon]['transcrit'],self.amplicons[amplicon]['contamination count'],"No"]
		for sample in samples_list:
			result_line.append(self.amplicons[amplicon][sample])
		results.writerow(result_line)
	if not contam:
		htmlOut.write('No contamination. <br/>\n')
		
	csvfile.close()
	
	# classement par odre decroissant de nombre de reads > read_len dans le fichier csv de resultat
	#data = csv.reader(open(self.OUTPUT_DIR+"/"+filename),delimiter=";")
	data = csv.reader(open(filename),delimiter=";")
	header = data.next()
	sortedlist = sorted(data, key=lambda x: int(x[6]), reverse=True)
	#with open((self.OUTPUT_DIR+"/"+filename),"wb") as f:
	with open(filename,"wb") as f:
		filewriter = csv.writer(f,delimiter=";")
		filewriter.writerow(header)
		for row in sortedlist:
			filewriter.writerow(row)
	f.close()
	#htmlOut.write('<br/><a href="../%s">%s</a><br>'%("downloads/"+filename, filename))
	htmlOut.write('<br/><a href="../%s">%s</a><br>'%(self.envDict['RESULTS_DIR'].split('/')[-1]+'/'+filename, filename))
	htmlOut.write('<br><b>VariantCaller for %s, using reads only > %s nucleotides : </b>' % (sample_tested,read_len))
	htmlOut.write('<br/><a href="../%s">%s</a><br>'%(self.envDict['RESULTS_DIR'].split('/')[-1]+'/variantCalling/alleles.xls', '%s_filtered.variantcaller.xls' % sample_tested))
	htmlOut.close()
	
        return True

if __name__ == "__main__":
    PluginCLI(Check_contamination())
