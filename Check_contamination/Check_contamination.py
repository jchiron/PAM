#!/usr/bin/env python
from ion.plugin import *
from subprocess import *
import os
import json
import csv
import xlwt
import pysam
import urllib2
import glob

class Check_contamination(IonPlugin):
    """Plugin object to check control sample contamination"""
    version = "0.1"
    envDict = dict(os.environ)
    json_dat = {}
    amplicons = {}
    barcodeNames = []
    sampleNameLookup = {} # dictionary allowing us to get the sample name associated with a particular barcode name
    
    def representsInt(self, s): # pour eviter avertissement "nombre ecrit en texte" sous excel
		try: 
			s = int(s)
			return s
		except ValueError:
			return s

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
	sample = False
	barcode = False
	read_len = False
	
	# Parse pluginconfig json.
	try:
		sample = self.json_dat['pluginconfig']['sample']
		barcode = self.json_dat['pluginconfig']['barcode']
		read_len = int(self.json_dat['pluginconfig']['read_len'])
		
	except:
		print 'Warning: plugin does not appear to be configured'
	
	htmlOut.write('<b> Settings :</b> <br/>\n')
	htmlOut.write('-Sample : %s<br/>\n' % sample)
	htmlOut.write('-Barcode : %s<br/>\n' % barcode)
	htmlOut.write('-Read length : %s<br/>\n' % read_len)
	
	# API recuperation parametres
	print "Gathering params..."
	ref = self.envDict['TSP_FILEPATH_GENOME_FASTA']
	target_unmerged = False
	target_merge = False
	hot_unmerged = False
#	hot_vcf = False
	param = False
	
	# NOTE : gathering param here is designed for variantcaller launched with only ONE targetbed, ONE hotspot, ONE parameter file for ALL samples
	try:
		api_url = self.json_dat['runinfo']['api_url'] + '/v1/pluginresult/?format=json&plugin__name=%s&result=%s' % ("variantCaller",str(self.json_dat['runinfo']['pk']))
		f = urllib2.urlopen(api_url)
		d = json.loads(f.read())
		for plugin in d['objects']: # they are already sorted from newest to oldest, so take the first of this list which is 'completed'
			if plugin['state'] == 'Completed':
				for data_file in plugin['store']['files']:
					if data_file['type'] == "parameters_json":	# get the same parameters as the last variantcaller
						param = data_file['server_path']
                                                param = param.replace('/_parameters.json','/local_parameters.json')
						print "  - json parameters found : " + param
						break
				for data_file in plugin['store']['files']:
					if data_file['type'] == "hotspots_bed":
						hot_unmerged = data_file['server_path']
						print "  - hot_unmerged found : " + hot_unmerged
						break
				break
	except:
		print 'ERROR!  Failed to get VC parameters (json parameters / hotspot bed file)'
	try:
		api_url = self.json_dat['runinfo']['api_url'] + '/v1/results/' + str(self.json_dat['runinfo']['pk'])
		f = urllib2.urlopen(api_url)
		d = json.loads(f.read())
		eas = d['eas']
		
		api_url = self.json_dat['runinfo']['net_location'] + eas
		f = urllib2.urlopen(api_url)
		d = json.loads(f.read())
		target_unmerged = d['targetRegionBedFile']	
		print "  - target_unmerge found : " + target_unmerged
		target_merge = target_unmerged.replace('unmerged/detail','merged/plain')
		print "  - target_merge should be : " + target_merge
	except:
		print 'ERROR!  Failed to get target bed path (unmerged / merged)'
		
	print "...Gathering params DONE"
	
#	param = '/results/plugins/Check_contamination/test_local_parameters.json'
#	target_unmerged = '/results/uploads/BED/78/hg19/unmerged/detail/IAD72953_231_Designed.with_NM.bed'
#	target_merge = target_unmerged.replace('unmerged/detail','merged/plain')
#	hot_unmerged = '/results/plugins/Check_contamination/test_hotspot.bed'
	
	# ouverture du bam a lire
	bampath = os.environ['ANALYSIS_DIR'] + "/" + str(barcode) + "_rawlib.bam"
	#bamfile = pysam.AlignmentFile(bampath,"rb") #pysam v0.8.1
        bamfile = pysam.Samfile(bampath,"rb") #pysam v0.7

	sample_name = sample.replace(" ","_")
	output_name = "%s_filtered.bam" % sample_name
	#bamfile_filtered = pysam.AlignmentFile(output_name, "wb", template=bamfile) #pysam v0.8.1
	bamfile_filtered = pysam.Samfile(output_name, "wb", template=bamfile) #pysam v0.7
	
	# ouverture du target bed, recuperation des infos de chaque amplicon
	bedfile = open(target_unmerged,"r")
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
	bedfile.close()
	
	# Creation bam avec reads > len
	for read in bamfile.fetch():
		#if (len(read.query_alignment_sequence) >= read_len): 	#pysam v0.8.1 # aligned portion of the read, exclude soft-clipped bases 
            if (len(read.query) >= read_len): 	#pysam v0.7 # aligned portion of the read, exclude soft-clipped bases 
                bamfile_filtered.write(read)
			
	bamfile.close()
	bamfile_filtered.close()
	
	# Remise en odre du bam et creation du BAI
	print "samtools sort..."
	cmd = Popen(["samtools","sort", output_name, "%s_filtered.sorted" % sample_name], stdout=PIPE)
	out, err = cmd.communicate()
	print 'OUT: %s\nERR: %s'%(out, err)
	
	print "samtools index..."
	cmd = Popen(["samtools","index","%s_filtered.sorted.bam" % sample_name], stdout=PIPE)
	out, err = cmd.communicate()
	print 'OUT: %s\nERR: %s'%(out, err)
	
	# Calcule nouvelle couverture par amplicon (script plugin coverageAnalysis)
	print "amplicon coverage for filtered bam..."
	os.makedirs("coverage")
	cmd = Popen(["/results/plugins/coverageAnalysis/run_coverage_analysis.sh","-ag","-D","coverage","-B",target_unmerged,ref,"%s_filtered.sorted.bam" % sample_name], stdout=PIPE)
	out, err = cmd.communicate()
	print 'OUT: %s\nERR: %s'%(out, err)
	
	for covFile in glob.glob('%s/coverage/*.amplicon.cov.xls' % (self.envDict['RESULTS_DIR'])):
		# ouverture du fichier de couverture par amplicon apres filtre des reads < len, recuperation des total reads
		coveragefile = open(covFile,"r")
		reader = csv.reader(coveragefile, delimiter = '\t')
		reader.next()
		for row in reader:
			amplicon = row[3]
			totalReads = row[9]
			self.amplicons[amplicon]["contamination count"] = int(totalReads)
	
	#####################################################################
	### LANCEMENT VARIANTCALLER POUR LE BAM CONTENANT LES READS > LEN ###
	#####################################################################

	os.makedirs("variantCalling")
	
	#print "preparing hotspot"
	#cmd = Popen(["tvcutils","prepare_hotspots","-b",hot_unmerged,"-o","variantCalling/hotspot_prepared.vcf","-r",ref,"--left-alignment","on","--allow-block-substitutions","on"], stdout=PIPE)
	#out, err = cmd.communicate()
	#print 'OUT: %s\nERR: %s'%(out, err)
	
	print "Launching VariantCaller"
	#cmd = Popen(["python","/results/plugins/variantCaller/bin/variant_caller_pipeline.py","-i","%s_filtered.sorted.bam" % sample_name,"-r",ref,"-b",target_merge,"-s","variantCalling/hotspot_prepared.vcf","-p",param,"-o","variantCalling","--primer-trim-bed",target_unmerged,"--postprocessed-bam","%s_filtered_processed.sorted.bam" % sample_name,"--error-motifs","/results/plugins/variantCaller/share/TVC/sse/motifset.txt"], stdout=PIPE)
        cmd = Popen(["python","/results/plugins/variantCaller/bin/variant_caller_pipeline.py","-i","%s_filtered.sorted.bam" % sample_name,"-r",ref,"-b",target_merge,"-p",param,"-o","variantCalling","--primer-trim-bed",target_unmerged,"--postprocessed-bam","%s_filtered_processed.sorted.bam" % sample_name,"--error-motifs","/results/plugins/variantCaller/share/TVC/sse/motifset.txt"], stdout=PIPE)
	out, err = cmd.communicate()
	print 'OUT: %s\nERR: %s'%(out, err)
	
	print "Generate variants table..."
	#cmd = Popen(["python","/results/plugins/variantCaller/scripts/generate_variant_tables.py","--input-vcf","variantCalling/TSVC_variants.vcf","--region-bed",target_unmerged,"--hotspots","--output-xls","variantCalling/variants.xls","--alleles2-xls","variantCalling/alleles.xls","--summary-json","variantCalling/variant_summary.json"], stdout=PIPE)
        cmd = Popen(["python","/results/plugins/variantCaller/scripts/generate_variant_tables.py","--input-vcf","variantCalling/TSVC_variants.vcf","--region-bed",target_unmerged,"--output-xls","variantCalling/variants.xls","--alleles2-xls","variantCalling/alleles.xls","--summary-json","variantCalling/variant_summary.json"], stdout=PIPE)
	out, err = cmd.communicate()
	#,"--scatter-png","scatter.png",
	print 'OUT: %s\nERR: %s'%(out, err)
	
	##################
	### ANNOTATION ###
	##################
	
	os.makedirs("variantAnnotation")
	print "Launching variantAnnotation"
	
	#cmd = Popen(["cp",hot_unmerged,self.envDict['RESULTS_DIR']+"/variantCalling/hotspot.bed"], stdout=PIPE) # pour etape generate vcf de l'annotation 
	#out, err = cmd.communicate()
	#print 'OUT: %s\nERR: %s'%(out, err)
	
	cmd = Popen(["python","/results/plugins/Check_contamination/variantAnnotation_custom/variant_annotation_plugin.py","--sample-name",sample,"--sample-id",barcode,"--sample-vcf",self.envDict['RESULTS_DIR']+"/variantCalling/TSVC_variants","--install-dir","/results/plugins/Check_contamination/variantAnnotation_custom","--output-dir",self.envDict['RESULTS_DIR']+"/variantAnnotation","--genome-fasta",ref], stdout=PIPE)
	out, err = cmd.communicate()
	print 'OUT: %s\nERR: %s'%(out, err)

	#################################################################
	### RAJOUT COLONNES COUVERTURE DE CHAQUE AMPLICON PAR PATIENT ###
	#################################################################
	
	try:
		reference_path = self.envDict['TSP_FILEPATH_GENOME_FASTA']
	except:
		reference_path = ''
		
	if isinstance(self.json_dat['plan']['barcodedSamples'],dict):
		samples = self.json_dat['plan']['barcodedSamples']
	else:
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
                                print 'coverageAnalysis completed found : ' + coverageAnalysisPath
				break
	except:
		print 'ERROR!  Failed to get coverage analysis path'
			
	if coverageAnalysisPath:
		for covFile in glob.glob('%s/plugin_out/%s/Ion*/*.amplicon.cov.xls' % (self.envDict['ANALYSIS_DIR'],coverageAnalysisPath.split('/')[-1])):
			barcodeName = covFile.split("/")[-2]
			sampleName = self.sampleNameLookup[barcodeName]
                        if sampleName != "BLANC_PCR":
                            print 'fichier couverture par amplicon : '+ covFile
			# ouverture du fichier de couverture par amplicon de chaque patient, recuperation des infos
                            coveragefile = open(covFile,"r")
                            reader = csv.reader(coveragefile, delimiter = '\t')
                            reader.next()
                            for row in reader:
				amplicon = row[3]
                                totalReads = row[9]
                                if amplicon in self.amplicons:
                                    self.amplicons[amplicon][sampleName] = totalReads

	else:
		print 'WARNING! No completed instance found for coverageAnalysis plugin. Coverage for each sample will not be displayed in final result'
	
	##############################
	### ecriture des resultats ###
	##############################
	
	# WORKBOOK STYLES
		
	headerFont = xlwt.Font()
	headerFont.name = 'Calibri'
	headerFont.height = 220
	headerFont.bold = 1
	headerBorders = xlwt.Borders()
	headerBorders.left = xlwt.Borders.MEDIUM
	headerBorders.right = xlwt.Borders.MEDIUM
	headerBorders.top = xlwt.Borders.MEDIUM
	headerBorders.bottom = xlwt.Borders.MEDIUM
	headerStyle = xlwt.XFStyle()
	headerStyle.font = headerFont
	headerStyle.borders = headerBorders
	generalFont = xlwt.Font()
	generalFont.name = 'Calibri'
	generalFont.height = 220
	generalStyle = xlwt.XFStyle()
	generalStyle.font = generalFont
	redBG = xlwt.Pattern()
	redBG.pattern = redBG.SOLID_PATTERN
	redBG.pattern_fore_colour = 0x1D #0x2D rose
	redBGStyle = xlwt.XFStyle()
	redBGStyle.pattern = redBG
	redBGStyle.font = generalFont
	
	#remove bug patient '' when present
	samples_list = self.sampleNameLookup.values()
	try:
		samples_list.remove('')
	except:
		pass
	
	results = xlwt.Workbook()
	ampliconSheet = results.add_sheet("Amplicons Contamination")
	variantAnnotationSheet = results.add_sheet("Annotation (filtered)")
	variantCallerSheet = results.add_sheet("variantCaller (filtered)")
	
	#Amplicon coverage Sheet	
	header = ["Amplicon","Chr","Start","End","Gene_id", "Transcript", "%s reads > %s nucleotides" % (sample,read_len)]
	for s in samples_list:
            if s != "BLANC_PCR":
		header.append(s+"_Total_Reads") 
	headerline = ampliconSheet.row(0)
	for i in range(len(header)):
		headerline.write(i,header[i],headerStyle)
	
	to_write = []
	for amplicon in self.amplicons.keys():
		line = [amplicon,self.amplicons[amplicon]['chr'],self.amplicons[amplicon]['startpos'],self.amplicons[amplicon]['endpos'],self.amplicons[amplicon]['gene_id'],self.amplicons[amplicon]['transcrit'],self.amplicons[amplicon]['contamination count']]
		for s in samples_list:
                    if s != "BLANC_PCR":
                        if s in self.amplicons[amplicon]:
                            line.append(self.amplicons[amplicon][s])
		to_write.append(line)
	to_write.sort(key=lambda x: x[6]) # classement par ordre decroissant de nombre de reads > read_len
	to_write.reverse()
	
	l=1
	for line in to_write:
		row = ampliconSheet.row(l)
		for i in range(len(line)):
                    if line[6]>1 :#au moins 2 reads de plus de 100b pour BLANC_PCR sur l'amplicon
			if i < 7:
				row.write(i,self.representsInt(line[i]),generalStyle)
			else:
				if int(line[i]) < 5*(int(line[6])):
					row.write(i,self.representsInt(line[i]),redBGStyle)
				else:
					row.write(i,self.representsInt(line[i]),generalStyle)
		l=l+1
	
	#VariantCaller Sheet
	vc_filtered = csv.reader(open(self.envDict['RESULTS_DIR']+"/variantCalling/alleles.xls"),delimiter="\t")
	headerline = variantCallerSheet.row(0)
	header = vc_filtered.next()
	for i in range(len(header)):
		headerline.write(i,header[i],headerStyle)
	l=1
	for vc_line in vc_filtered:
		if vc_line[4] == "Homozygous" or vc_line[4] == "Heterozygous":
			row = variantCallerSheet.row(l)
			for i in range(len(vc_line)):
				row.write(i,self.representsInt(vc_line[i]),generalStyle)
			l=l+1
	
	#VariantAnnotation Sheet
	try:
		if " " in sample: # correction bug quand nom patient contient espace, le fichier tsv est creer avec la partie du nom patient avant l'espace..
			cutname_sample = sample.split(" ")[0]
			va_filtered = csv.reader(open(self.envDict['RESULTS_DIR']+"/variantAnnotation/%s/CHR/%s_NGS_Diag_VC_TS.tsv" % (barcode,cutname_sample)),delimiter="\t")
		else:
			va_filtered = csv.reader(open(self.envDict['RESULTS_DIR']+"/variantAnnotation/%s/CHR/%s_NGS_Diag_VC_TS.tsv" % (barcode,sample)),delimiter="\t")
		headerline = variantAnnotationSheet.row(0)
		header = va_filtered.next()
		for i in range(len(header)):
			headerline.write(i,header[i],headerStyle)
		l=1
		for va_line in va_filtered:
			row = variantAnnotationSheet.row(l)
			for i in range(len(va_line)):
				if i == 17 and va_line[17] != "NA": # test presence COSMIC
					row.write(i,self.representsInt(va_line[i]),redBGStyle)
				else:
					row.write(i,self.representsInt(va_line[i]),generalStyle)
			l=l+1
	except:
		print "Warning : annotation not found (could mean absence of variants in variantcaller). VariantAnnotation Sheet is empty"
		variantAnnotationSheet.row(0).write(0,"Annotation file not found (could mean absence of variants in variantCaller results) ")
	
	
	filename = "Check_contamination_%s.xls" % sample_name
	results.save("%s/%s" % (self.envDict['RESULTS_DIR'],filename))
	
	htmlOut.write('<br/><b> RESULTS :</b> <br/>\n')
	htmlOut.write('<br/><a href="../%s">%s</a><br>'%(self.envDict['RESULTS_DIR'].split('/')[-1]+'/'+filename, filename))
	htmlOut.close()
	
        return True

if __name__ == "__main__":
    PluginCLI(Check_contamination())
