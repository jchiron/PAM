#!/usr/bin/env python

from ion.plugin import *
import sys
import os
import json
import urllib2
import math
import matplotlib
matplotlib.use('Agg')
from pylab import *
import numpy as np
from fnmatch import fnmatch
import zipfile
import glob
import csv
import re

class plotCoverage(IonPlugin):
    """Plugin object which compare two variant caller output"""
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
		htmlOut = open('%s/plugin_out/plotCoverage_out/plotCoverage.html'%self.envDict['ANALYSIS_DIR'], 'a+')
	except:
		htmlOut = open('plotCoverage.html', 'w')
	#htmlOut.write('<html><body>\n')
	
	if os.path.isfile(self.envDict['TSP_FILEPATH_BARCODE_TXT']):
		self.isBarcodedRun = True
			
	if not self.isBarcodedRun:
		print "This plugin only work for barcoded run for now. \n Exiting..."
		sys.exit(0)
		
	
	# set options
	coverage_path = False
	result_dir = False
        variantCaller_path = False
	target_path = {}
	
	# Parse pluginconfig json.
	try:
		result_dir=self.json_dat['runinfo']['results_dir']
		
	except:
		print 'Warning: plugin does not appear to be configured'
	
	try:
		#api_url = self.json_dat['runinfo']['api_url'] + '/v1/pluginresult/?format=json&plugin__name=coverageAnalysis&result=' + str(self.json_dat['runinfo']['pk'])
                api_url = os.getenv('RUNINFO__API_URL',
                           'http://localhost/rundb/api') + '/v1/pluginresult/?format=json&plugin__name=coverageAnalysis&result=' + str(self.json_dat['runinfo']['pk'])
                #api_url = '129.21.10.3/v1/pluginresult/?format=json&plugin__name=coverageAnalysis&result=' + str(self.json_dat['runinfo']['pk'])
                #print '!!!!!!!!!!API URL:' + api_url
                
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
				coverage_path = plugin['path']
				#coverage_path = coverage_path.split('/')[-1]
				break
		if not coverage_path:
			print 'WARNING! No completed instance found for plugin coverageAnalysis'
	
		#api_url = self.json_dat['runinfo']['api_url'] + '/v1/pluginresult/?format=json&plugin__name=variantCaller&result=' + str(self.json_dat['runinfo']['pk'])
		api_url = os.getenv('RUNINFO__API_URL',
                           'http://localhost/rundb/api') + '/v1/pluginresult/?format=json&plugin__name=variantCaller&result=' + str(self.json_dat['runinfo']['pk'])
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
				variantCaller_path = plugin['path']
                                print 'variantCaller path:' + variantCaller_path
				break
		if not coverage_path:
			print 'WARNING! No completed instance found for plugin coverageAnalysis'
	
                #get target for specific barcode
                for barcode in d['objects'][0]['store']['barcodes']:
                    target = d['objects'][0]['store']['barcodes'][barcode]['targets_bed']
                    #print 'BARCODE : ' + barcode + ' TARGET : ' + target
                    target_path[barcode] = target
                    
	except:
            print 'ERROR!  Failed to get coverage analysis path or target bed path or variantCaller path'

        
	
	# Get bam filenames.
	with open(os.path.join(self.json_dat['runinfo']['basecaller_dir'], 'datasets_basecaller.json'), 'r') as f:
		json_basecaller = json.load(f)
        
        # test if/else for version 4.6
        barcodeSamples_fromJson = self.json_dat['plan']['barcodedSamples']
        if isinstance(barcodeSamples_fromJson,dict):
            samples = barcodeSamples_fromJson
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
			print 'BARCODE FOUND: %s SAMPLE ID: %s' % (barcodeName, sampleName)
			self.barcodeNames.append(barcodeName)
		else:
			print 'NO SAMPLE NAME associated with barcode %s => maybe contamination' % (barcodeName)			
			sampleName = ''
			self.sampleNameLookup[barcodeName] = 'None' # makes it much easier later to do the lookup with no conditional tests
                # MGD note: I considered setting blank sample names to the barcode name instead, but might not be what customer intended
                print 'BARCODE FOUND: %s SAMPLE ID: %s' % (barcodeName, sampleName)
                self.barcodeNames.append(barcodeName)
		
	########### PROCESS ################

	#find samples amplicon.cov.xls files in IonXpress subdirectories in coverageAnalysis output dir

	cov_pattern = "*.amplicon.cov.xls"
	#list of containing amplicon.cov.xls files for all samples
	cov_files = []

	for dirPath, subdirName, fileName in os.walk(coverage_path):
	    #print('Found directory: %s' % dirPath)
	    for fname in fileName:
		if fnmatch(fname,cov_pattern):
		    if not fnmatch (fname,"link.amplicon.cov.xls"):
			#print os.path.join(dirPath,fname)
			cov_files.append(os.path.join(dirPath,fname))
	#samples id list
	samples_list = {}
        amplicon_chr = {}

	#chromosomes list
	chromosomes = ''

	#read amplicon_cov file to find amplicon with less than 500 reads mapped
	#hash key= "amplicon id" value="sample name(total reads mapped on amplicon)"
	amplicon_sample_cov = {}

        amplicon_cov = {}

        print "nb amplicon.cov.xls files :" + str(len(cov_files)) + "\n"
	if len(cov_files) >=1:
            
	    #hash key="sample name" value="string tab separated of total reads per amplicon"
	    amplicon_cov_for_sample = {}

            #amplicon id list
            amplicon_list = {}
                        
            targets_list = []
            for val in target_path.values():
                if val in targets_list:
                    continue
                else:
                    targets_list.append(val)

            for target in targets_list:
                #hash key="amplicon id" value="string tab separated of total reads per amplicon"
                amplicon_cov[target] = {}
                amplicon_cov_for_sample[target] = {}
                samples_list[target] = []

            
                if "NGHS-102X_BRCA1BRCA2" in target: #sur les panels brca pour paola
                    ##top block html file
                    TSP_FILEPATH_PLUGIN_DIR = os.environ['TSP_FILEPATH_PLUGIN_DIR']
                    block_html_file = TSP_FILEPATH_PLUGIN_DIR + "/variantAnnotation_block.html"
                    block_html=open(block_html_file,'w')
                    block_content = "<!DOCTYPE html>\n" \
                                    "<html lang=\"en\">\n" \
                                    "<head>\n" \
                                    "<base target=\"_parent\"/>\n" \
                                    "<link rel=\"stylesheet\" media=\"all\" href=\"/site_media/resources/bootstrap/css/bootstrap.min.css\">\n" \
                                    "<link href=\"/site_media/resources/kendo/styles/kendo.common.min.css\" rel=\"stylesheet\">\n" \
                                    "<link href=\"/site_media/resources/less/kendo.tb.min.css\" rel=\"stylesheet\">\n" \
                                    "<link type=\"text/css\" rel=\"stylesheet\" href=\"/site_media/resources/styles/tb-styles.min.css\">\n" \
                                    "</head>\n\n" \
                                    "<body>\n" \
                                    "<b>Percentage of CDS coverage with a depth of 300x calculated on rawlib_processed.bam </b><br/>\n"\
                                    "<div class=\"k-widget k-grid\">\n" \
                                    "<table class=\"table-striped\">\n" \
                                    "    <thead class=\"k-grid-header\">\n" \
                                    "      <tr>\n" \
                                    "        <th><span class=\"help\" title=\"Sample name\">Sample name</span></th>\n" \
                                    "        <th><span class=\"help\" title=\"BRCA1\">BRCA1</span></th>\n" \
                                    "        <th><span class=\"help\" title=\"BRCA2\">BRCA2</span></th>\n" \
                                    "      </tr>\n" \
                                    "    </thead>\n" \
                                    "    \n"

                    block_html.write(block_content)
                    block_html.close()
            
	    #parse sample amplicon cov file 
	    for file_path in cov_files:
		#print file_path
		dir_path,sample_file = os.path.split(file_path)
		#print dir_path
		#print sample_file
		
		sample_id=os.path.basename(dir_path)
                if sample_id in target_path.keys():
                    target = target_path[sample_id]
                else:
                    continue
		if sample_id in self.sampleNameLookup:
			sample_name = self.sampleNameLookup[sample_id]+'_'+sample_id
			#samples_list.append(sample_name)
        	        # !!! TODO: decommenter >>>>
			#control = 'BLANC_PCR'
                        #control_pcr = 'BLANC_'
                        #control_pcr1 = 'Blanc_'
                        safir = 'SAFIR'
                        #diff_panel = 'ROBOT'
                        #diff_panel = 'repassage'
                        sample1 ='BOURDENS_'
                        sample2 = 'VIRELEGOUX_LN398'

			if re.search(safir,sample_name,re.I):
# re.search(control,sample_name,re.I) or re.search(safir,sample_name,re.I) or re.search(control_pcr,sample_name,re.I) or re.search(control_pcr1,sample_name,re.I):
# or re.search(diff_panel,sample_name,re.I):
# or re.search(sample1,sample_name,re.I) or re.search(sample2,sample_name,re.I) :
# or re.search(diff_panel,sample_name,re.I):
			    #print 'H2O control not analyzed'
			    continue
                        #if not re.search(diff_panel,sample_name,re.I) :
                        #if not re.search(sample1,sample_name,re.I) and not re.search(sample2,sample_name,re.I):
                        #if not re.search(sample1,sample_name,re.I):
                        #if not re.search('IonXpress_086',sample_id,re.I) and not re.search('IonXpress_084',sample_id,re.I):
                        #    print "not analysed"
                        #    continue
			#add sample to the list of samples processed
                        #if re.search('IonXpress_086',sample_id,re.I) or re.search('IonXpress_084',sample_id,re.I):
                        #    print "not analysed"
                        #    continue
			samples_list[target].append(sample_name)
			# <<<< fin decommenter 
	
			print sample_name
#process less than 300X coord
			print('find position less than 300X detph calculated on rawlib_processed.bam')
			prefix = sample_file.split(".")[0]
			#print prefix
			#bam_path = dir_path +"/"+prefix+".bam"
                        #bam_path = variantCaller_path+"/"+sample_id+"/"+sample_id+"_rawlib_processed.bam"
			bam_path = variantCaller_path+"/"+sample_id+"/"+sample_id+"_rawlib.realigned_processed.bam"
                        #print bam_path

			#>>>>>>>>>>>> decommenter et tester une fois le fichier panel recuperer
			sample_cov_path = result_dir+"/"+sample_name.replace(" ","-")+"_cov.tsv"
			cmd = "bedtools coverage -d -abam %s -b %s > %s" % (bam_path, target_path[sample_id], sample_cov_path)
			os.system(cmd)

			cov_300x_path = result_dir+"/"+sample_name.replace(" ","-")+"_cov_max_300x.bed"
			cov_300x_file = open(cov_300x_path,'w')
			cov300x_writer = csv.writer(cov_300x_file, delimiter = "\t")

			chr_pos_final = result_dir+"/"+sample_name.replace(" ","-")+"_chr_pos_300x_max.tsv"
                        chr_pos_path = result_dir+"/"+sample_name.replace(" ","-")+"_chr_pos_300x_max.tmp"
			chr_pos_file =open(chr_pos_path,'w')
			chr_pos_writer = csv.writer(chr_pos_file, delimiter = "\t", lineterminator="\n")

			sample_cov_file = open(sample_cov_path,'r')
			for line in sample_cov_file:
			    line = line.rstrip('\n')
			    chrom = line.split('\t')[0]
			    start = line.split('\t')[1]
			    end = line.split('\t')[2]
			    name = line.split('\t')[3]
			    pos_in_region = line.split('\t')[8]
			    depth = line.split('\t')[9]
			    depth = int(depth)
	
			    if depth<300:
				cov300x_writer.writerow(line.split('\t'))
				pos_in_chrom = int(start)+int(pos_in_region)
				chr_pos_writer.writerow([chrom,str(pos_in_chrom),str(pos_in_chrom)])
			cov_300x_file.close()
			chr_pos_file.close()
			sample_cov_file.close()
                      
                        #remove duplicates in chr_pos_file due to overlaping amplicons
                        os.system('sort -u %s > %s' % (chr_pos_path,chr_pos_final))

                        print "target :" + target_path[sample_id]
                        if "NGHS-102X_BRCA1BRCA2" in target_path[sample_id]: #sur les panels brca pour paola
                            TSP_FILEPATH_PLUGIN_DIR = os.environ['TSP_FILEPATH_PLUGIN_DIR']
                            block_html_file = TSP_FILEPATH_PLUGIN_DIR + "/variantAnnotation_block.html"
                           
                            #> debut stats BRCA
                            print 'stats BRCA'
                            #htmlOut.write('<b>Percentage of CDS coverage with a depth 300x</b><br/>')
                            #htmlOut.write('Sample Name<TAB/>BRCA1<TAB/>BRCA2<br/>')
                            bc1 = "/results/plugins/plotCoverage/brca_files/BRCA1_coding_ensembl75.sort.bed"
                            bc2 = "/results/plugins/plotCoverage/brca_files/BRCA2_coding_ensembl75.sort.bed"

                            bc1_cds_size = 5592
                            bc2_cds_size = 10257
                        
                            #fichier bed des bases bc1 couvertes a max 300x
                            bc1_300x = result_dir+"/"+sample_name.replace(" ","-")+"_cov_less300x_BC1.tsv"
                            cmd = "bedtools intersect -a %s -b %s > %s" % (chr_pos_final, bc1, bc1_300x)
                            os.system(cmd)
                            line_count=sum(1 for line in open(bc1_300x))
                            pc_bc1_300x = 100-(line_count*100/float(bc1_cds_size))
                            pc_bc1_300x_float = "%.2f" % pc_bc1_300x
                            print 'BRCA1 :' + str(pc_bc1_300x_float) + '% covered at minimum 300x' 
                        
                            #fichier bed des bases bc2 couvertes a max 300x
                            bc2_300x = result_dir+"/"+sample_name.replace(" ","-")+"_cov_less300x_BC2.tsv"
                            cmd = "bedtools intersect -a %s -b %s > %s" % (chr_pos_final, bc2, bc2_300x)
                            os.system(cmd)
                            line_count=sum(1 for line in open(bc2_300x))
                            pc_bc2_300x = 100-(line_count*100/float(bc2_cds_size))
                            pc_bc2_300x_float = "%.2f" % pc_bc2_300x
                            print 'BRCA2 :' + str(pc_bc2_300x_float) + '% covered at minimum 300x'

                            #htmlOut.write('%s<TAB/>%s<TAB/>%s<br/>' % (sample_name,pc_bc1_300x,pc_bc2_300x))
                            block_html=open(block_html_file,'a')
                           
                            block_content = "<tr>\n" \
                                            "<th><span class=\"help\" title=\"Sample name\">%s</span></th>\n" \
                                            "<th><span class=\"help\" title=\"BRCA1\">%s</span></th>\n" \
                                            "<th><span class=\"help\" title=\"BRCA2\">%s</span></th>\n" \
                                            "</tr>\n" \
                                            " \n" % (sample_name,pc_bc1_300x_float,pc_bc2_300x_float)
                                            
                            block_html.write(block_content)
                            block_html.close()                            
                            #<fin stats BRCA

		#write region with less thna 300X depth in bed file
			region_300X_path = result_dir+"/"+sample_name+"_cov_less_than_300x.bed"
			region_300X_file = open(region_300X_path,'w')
			region_300X_writer = csv.writer(region_300X_file, delimiter = "\t")

			chr_pos_infile = open(chr_pos_final,'r')
        		chr_infile = csv.reader(chr_pos_infile,delimiter='\t')

			chrom = ''
			start = ''
			end = ''
			for line in chr_infile:
				if chrom=='' and start=='':
					chrom = line[0]
					start = line[1]
					end = line[1]
					previous_pos = start
					
				else:
					if line[0]==chrom:
						if int(line[1])==(int(previous_pos)+1) :
							end = line[1]
							previous_pos = line[1]
						else:
							region_300X_writer.writerow([chrom,start,end])
							chrom = line[0]
							start = line[1]
							end = line[1]
							previous_pos = line[1]
					else:
						region_300X_writer.writerow([chrom,start,end])
						chrom = line[0]
						start = line[1]
						end = line[1]
						previous_pos = line[1]

			#last region
			region_300X_writer.writerow([chrom,start,end])

			region_300X_file.close()
			chr_pos_infile.close()
		
			#delete temp files
			cmd = "rm " + sample_cov_path
			os.system(cmd)			
			
			cmd = "rm " + cov_300x_path
			os.system(cmd)
		
			cmd = "rm " + chr_pos_path
			os.system(cmd)
                        
                        cmd = "rm " + chr_pos_final
			os.system(cmd)
                        
			#<<<<<<<<<<<<< fin decommenter

                        '''
			#process less than 100X coord
			print('find position less than 100X detph')
			prefix = sample_file.split(".")[0]
			#print prefix
			bam_path = dir_path +"/"+prefix+".bam"
			#print bam_path

			#>>>>>>>>>>>> decommenter et tester une fois le fichier panel recuperer
			sample_cov_path = result_dir+"/"+sample_name.replace(" ","-")+"_cov.tsv"
			cmd = "bedtools coverage -d -abam %s -b %s > %s" % (bam_path, target_path, sample_cov_path)
			os.system(cmd)

			cov_100x_path = result_dir+"/"+sample_name+"_cov_less_than_100x.bed"
			cov_100x_file = open(cov_100x_path,'w')
			cov100x_writer = csv.writer(cov_100x_file, delimiter = "\t")
			chr_pos_path = result_dir+"/"+sample_name+"_chr_pos_100x_max.tsv"
			chr_pos_file =open(chr_pos_path,'w')
			chr_pos_writer = csv.writer(chr_pos_file, delimiter = "\t")

			sample_cov_file = open(sample_cov_path,'r')
			for line in sample_cov_file:
			    line = line.rstrip('\n')
			    chrom = line.split('\t')[0]
			    start = line.split('\t')[1]
			    end = line.split('\t')[2]
			    name = line.split('\t')[3]
			    pos_in_region = line.split('\t')[8]
			    depth = line.split('\t')[9]
			    depth = int(depth)
	
			    if depth<100:
				cov100x_writer.writerow(line.split('\t'))
				pos_in_chrom = int(start)+int(pos_in_region)
				chr_pos_writer.writerow([chrom,str(pos_in_chrom)])
			cov_100x_file.close()
			chr_pos_file.close()
			sample_cov_file.close()
			#<<<<<<<<<<<<< fin decommenter
                        '''
                        
                        amplicon_list[target] = ''

			#find amplicon with less than 500 reads mapped on
		
			cov_file = open(file_path,'r')
			for line in cov_file:
			    if line.startswith('contig_id'):
				next;
			    else:
				line = line.rstrip('\n')
		    
				amplicon_id = line.split('\t')[3]
				total_reads = line.split('\t')[9]
				if int(total_reads)<500:
				    if amplicon_id in amplicon_sample_cov.keys():
					if sample_name in amplicon_sample_cov[amplicon_id].keys():
					    amplicon_sample_cov[amplicon_id][sample_name] = total_reads
					else:
					    amplicon_sample_cov[amplicon_id][sample_name] = total_reads
				    else:
					amplicon_sample_cov[amplicon_id] = {}
					if sample_name in amplicon_sample_cov[amplicon_id].keys():
					    amplicon_sample_cov[amplicon_id][sample_name] = total_reads
					else:
					    amplicon_sample_cov[amplicon_id][sample_name] = total_reads
				else:
				    if amplicon_id in amplicon_sample_cov.keys():
					if sample_name in amplicon_sample_cov[amplicon_id].keys():
					    amplicon_sample_cov[amplicon_id][sample_name] = " "
					else:
					    amplicon_sample_cov[amplicon_id][sample_name] = " "
				    else:
					amplicon_sample_cov[amplicon_id] = {}
					if sample_name in amplicon_sample_cov[amplicon_id].keys():
					    amplicon_sample_cov[amplicon_id][sample_name] = " "
					else:
					    amplicon_sample_cov[amplicon_id][sample_name] = " "
			cov_file.close()

			#amplicon_cov.xls file is sorted by total reads in amplicon, reorder file by genomic position to plot values
			file_path = file_path.replace(" ","\ ")
			sample_name = sample_name.replace(" ","\ ")
			cmd = "sort +0.3n -1 +1n -2 +2n -3 " + file_path + " > " + result_dir + "/" + sample_name + "_amplicon.cov.sorted.tsv"
			#print cmd
			os.system(cmd)

			###############################
			#find amplicons cov for sample

			chromosomes = ''
                        amplicon_chr[target] = {}
			previous_chr = None
			sample_name = sample_name.replace("\ "," ")
			sample_sorted = sample_name + "_amplicon.cov.sorted.tsv"
			sorted_file = open(os.path.join(result_dir,sample_sorted),'r')
			for line in sorted_file:
                            if line.startswith('contig_id'):
				next;
			    else:
				line = line.rstrip('\n')
			    
				amplicon_id = line.split('\t')[3]
				total_reads = line.split('\t')[9]

                                if (total_reads=='0'): #debug: pas de reads sur l'amplicon mais log(10)retourne NA
                                    total_reads='1'
				total_reads = log10(float(total_reads))
	
				total_reads = str(total_reads)

				#retrieve chromosomes list
				chr_id = line.split('\t')[0]
                                amplicon_chr[target][amplicon_id] = chr_id
				if (previous_chr is None):
				    chromosomes = chr_id
				    previous_chr = chr_id
				else:
				    if (chr_id == previous_chr):
					next
				    else:
					chromosomes = chromosomes + "\t" + chr_id
					previous_chr = chr_id
	
                                if sample_name in amplicon_cov_for_sample[target].keys():
				    amplicon_cov_for_sample[target][sample_name] = amplicon_cov_for_sample[target][sample_name] + "\t" + total_reads
				else:
				    amplicon_cov_for_sample[target][sample_name] = total_reads
				    #print "dict " + amplicon_cov_for_sample[sample_name]

				if amplicon_list[target] is not None:
				    amplicon_list[target] = str(amplicon_list[target]) + "\t" + amplicon_id
				else:
				    amplicon_list[target] = amplicon_id

			    
				#hash key=amplicon value=cov
				if amplicon_id in amplicon_cov[target].keys():
				    amplicon_cov[target][amplicon_id] = amplicon_cov[target][amplicon_id] + "\t" + total_reads
				else:
				    amplicon_cov[target][amplicon_id] = total_reads
	
			sorted_file.close()
			#print chromosomes
			#print amplicon_list

		######################
		#write result file 
	
		#tab file containing header=region_id and lines=samples amplicon total reads for which that have less than 500 reads
		uncov_path = result_dir + "/amplicon.uncov.tsv"
		uncov_file = open(uncov_path,'w')

		#print header = amplicon id
		for sample in samples_list[target]:
		    uncov_file.write("\t" + sample)
		uncov_file.write("\n")
		#values for amplicon
		for amplicon_id in amplicon_sample_cov.keys():
		    uncov_file.write(amplicon_id)
		    for sample in samples_list[target]:
                        if sample in amplicon_sample_cov[amplicon_id]:
                            uncov_file.write("\t" + amplicon_sample_cov[amplicon_id][sample])
                        else:
                            uncov_file.write("\t ")
		    uncov_file.write("\n")

            for sample in target_path.keys():
                if "NGHS-102X_BRCA1BRCA2" in target_path[sample_id]: #sur les panels brca pour paola
                    block_html=open(block_html_file,'a')
                    block_content = "</table>\n</div></body></html>\n" 
                    block_html.write(block_content)
                    block_html.close()

        ### htmlOut ###

        htmlOut.write('''
        <?xml version="1.0" encoding="iso-8859-1"?>
        <!DOCTYPE HTML>
        <html>
        <!-- <base target="_parent"/> -->
        <META HTTP-EQUIV="Content-Type" CONTENT="text/html; charset=utf-8">
        <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        
        <link rel="stylesheet" media="all" href="/site_media/resources/bootstrap/css/bootstrap.min.css" />
        <link href="/site_media/resources/kendo/styles/kendo.common.min.css" rel="stylesheet" />
        <link href="/site_media/resources/less/kendo.tb.min.css" rel="stylesheet" />
        <link type="text/css" rel="stylesheet" href="/site_media/resources/styles/tb-layout.css" />
        <link type="text/css" rel="stylesheet" href="/site_media/resources/styles/tb-styles.min.css" />

        <link rel="stylesheet" type="text/css" href="/site_media/stylesheet.css"/>
        <link rel="stylesheet" type="text/css" href="/site_media/resources/styles/print.css" media="print" />
        <link rel="stylesheet" type="text/css" href="/site_media/resources/styles/report.css" media="screen" />

        <script type="text/javascript" src="/site_media/resources/jquery/jquery-1.8.2.min.js"></script>
        <script type="text/javascript" src="/site_media/resources/scripts/kendo.custom.min.js"></script>

        <style type="text/css">
        body {background:white}
        .help {cursor:help; border-bottom: 1px dotted #A9A9A9}
        </style>

        </head>
        </head>

        <title>Torrent Coverage Plot Report</title>
        <body>

        <div class="container-fluid">
        
        <h1><center>Coverage Plot Report</center></h1>
        <h3><center>%s<center></h3>
        ''' % self.envDict['TSP_ANALYSIS_NAME'])

        

        nb_target = 0
        for target in targets_list:

            htmlOut.write('''
            <br/>
            <h2><center>Library %s</center></h2>
            ''' % target)

            nb_target = nb_target+1
            #list of files containing amplicon log10(total reads) : one file with all amplicons and 4 others (subset of amplicon to have less data in 1 plot to interpret easier)
            files_list = []

            #tab file containing header=region_id and lines=samples amplicon total reads for all amplicons and all samples
            all_cov_path = result_dir + "/all_amplicon_samples-target" + str(nb_target) + ".cov.tsv"
            all_cov_file = open(all_cov_path,'w')

            #print header = amplicon id
            all_cov_file.write("sample" + amplicon_list[target] + "\n")

            #print samples data
            for sample in samples_list[target]:
                #print sample
                all_cov_file.write(sample + "\t" + amplicon_cov_for_sample[target][sample] + "\n")
            all_cov_file.close()

            files_list.append(all_cov_path)

            #convert column to row to make file usable by biologist
            all_infile = open(all_cov_path,'r')
            csv_all_infile = csv.reader(all_infile,delimiter='\t')
            lines = []
        
            for line in csv_all_infile:
                lines.append(line)
                #print line        
            all_infile.close()

            transposed = [[lines[j][i] for j in range(len(lines))] for i in range(len(lines[0]))]
        
            log10_cov_file =  result_dir + "/all_samples_amplicons_log10.target" + str(nb_target) + ".cov.tsv"
            log10_out_file = open(log10_cov_file,'w')
       
            for i in range(len(transposed)):
                for j in range(len(transposed[i])):
                    if j==0:
                        #log10_out_file.write(amplicon_chr[transposed[i][j])
                        log10_out_file.write(transposed[i][j])
                    else:
                        line = "\t"+transposed[i][j]
                        line = line.replace(".",",")
                        log10_out_file.write(line)
                log10_out_file.write("\n")
            log10_out_file.close()
            #cmd = "sed -i 's/\./,/g' "+log10_cov_file
            #os.system(cmd)

            ##generate amplicon and subset files
            if len(samples_list[target])>10:
		#subset data files all amplicov cov for 8 samples max
		chunks = [samples_list[target][i:i+8] for i in xrange(0,len(samples_list[target]),8)]
		sample_subset = 0
		for chunk in chunks:
		    sample_subset = sample_subset + 1
		    sample_subset_path = result_dir+"/all_amplicon-sample_subset"+str(sample_subset)+ "-target" + str(nb_target) + ".tsv"
		    sample_subset_file = open(sample_subset_path,'w')
		    sample_subset_file.write("sample" + amplicon_list[target] + "\n")
		    for sample in chunk:
			sample_subset_file.write(sample + "\t" + amplicon_cov_for_sample[target][sample] + "\n")
		    sample_subset_file.close()
		    files_list.append(sample_subset_path)


            #4 subset graphs to view less amplicon info in one graph for all samples
            cmd  = "head -n 1 "+ all_cov_path +"| wc -w"
            nb_columns = os.popen(cmd).read()
            #print "COLONNES = " + nb_columns
            parts = 4
            nb_col_subset = int((int(nb_columns)-1)/parts)
            #print "SUBSETS = " + str(nb_col_subset) + "\n"

            start = 1
            end = start + nb_col_subset - 1
            for i in range(1,5):
                if (end>nb_columns):
                    end = nb_columns
                cmd = "cut -f 1,"+str(start)+"-"+str(end)+" "+all_cov_path+" > "+result_dir+"/amplicon_subset"+str(i)+"-all_samples"+"-target"+str(nb_target)+".tsv"
                #print cmd + "\n"
                os.system(cmd)
                subset_file = result_dir+"/amplicon_subset"+str(i)+"-all_samples"+"-target"+str(nb_target)+".tsv"
                files_list.append(subset_file)
	    
                if len(samples_list[target])>10:
		    #subset files containing 8 samples max and amplicons_number/4 (mandatory by IH)
		    sample_subset = 0
		    for chunk in chunks:
			sample_subset = sample_subset + 1
			sample_subset_path = result_dir+"/all_amplicon-sample_subset"+str(sample_subset)+"-target"+str(nb_target)+".tsv"
			cmd = "cut -f 1,"+str(start)+"-"+str(end)+" "+sample_subset_path+" > "+result_dir+"/amplicon_S"+str(i)+"-samples_S"+str(sample_subset)+"-target"+str(nb_target)+".tsv"
			#print cmd + "\n"
			os.system(cmd)
			sub_file_path = result_dir+"/amplicon_S"+str(i)+"-samples_S"+str(sample_subset)+"-target"+str(nb_target)+".tsv"
			files_list.append(sub_file_path)
	    
                start = start + nb_col_subset
                end = end + nb_col_subset
	    

            #add chromosome info to log10 cov file
            log10_cov_infile = open(log10_cov_file,'r')
            #csv_all_infile = csv.reader(log10_cov_infile,delimiter='\t')
            #lines = []
        
            log10_cov_file_chr = result_dir + "/all_samples_amplicons_log10_with_chr.target"+str(nb_target)+".cov.tsv"
            log10_out_file_chr = open(log10_cov_file_chr,'w')

            for line in log10_cov_infile:
                line = line.rstrip('\n')
                amplicon = line.split('\t')[0]
                if amplicon=="sample":
                    log10_out_file_chr.write(" \t"+line+"\n")
                else:
                    #cmd = "sed -i 's/,/\./g' "+log10_cov_file
                    #os.system(cmd)
                    log10_out_file_chr.write(amplicon_chr[target][amplicon]+"\t"+line+"\n")
                #lines.append(line)
                #print line        
            log10_cov_infile.close()
            log10_out_file_chr.close()

            cmd = "mv "+log10_cov_file_chr+" "+ log10_cov_file
            os.system(cmd)
        
        
       

            #############################
            #create amplicons cov graphes
            cpt = 1
            for f in files_list:
                print f + "\n"
	    
                #figure(figsize=(150,50)
                
                fname = f.split("/")[-1]
                prefix = fname.split(".")[0]
	    
                #print "FILE NAME = "+fname+"\n"
                #print "PREFIX = "+prefix+"\n"

                cov_file = open(f,'r')
                
                my_colors = ['blue','green','tomato','grey','sienna','orange','palevioletred','darkorchid','lime','magenta','dodgerblue','lightcoral','lightblue','darkkhaki','cyan','yellow','plum','peru','steelblue','mediumspringgreen','red','yellowgreen','maroon','gold','black','darkseagreen','burlywood','deeppink','slategrey','greenyellow','lightpink','fuchsia','blueviolet','navy','peachpuff']
                nb = 0	    
                for line in cov_file:
                    if line.startswith('sample'):
                        x_names = line.split('\t')
                        x_names.pop(0) #remove first element of the list (header line, remove "sample" label and retreive amplicon id for x_axis)
                        x_len = len(x_names)
                        x_data = range(len(x_names))
                        x = array(x_data)

                    else:
                        y_data = line.split('\t')
                        y_legend = y_data[0] #first field of the line = sample name for legend
                        y_data.pop(0)
                        y2_data =[]
                        #fig=figure(1,figsize=(50,20))
                        fig=figure(1,figsize=(50,10))
                        for n in y_data:
                            y2_data.append(n)
		
                        y2 = array(y2_data)
                        my_color = my_colors[nb]
                        plot(x,y2,label=y_legend,color=my_color,linewidth=3.5)
                        nb = nb + 1 
                cov_file.close()
                width = 200
                height = 5000
                xticks(x_data,x_names,rotation=90,fontsize=20)
                yticks(fontsize=14)
                ylim(ymin=0,ymax=5)
                xlabel("amplicons",fontsize=20)
                ylabel("total reads (log10)",fontsize=20)
                axhline(y=3,color='black',linestyle='--',label='1000')
                axhline(y=2.7,color='b',linestyle='--',label='500')
                axhline(y=2.48,color='orange',linestyle='--',label='300')
                axhline(y=2,color='r',linestyle='--',label='100')
                #lgd = legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=4, mode="expand", borderaxespad=0.)
                legend(loc='center left',bbox_to_anchor=(1,0.5)) #legend(fontsize=20,loc='center left',bbox_to_anchor=(1,0.5))
                #autoscale(enable=True)
                #tight_layout(rect=[0,0,0.9,1])
	    
                #show()
		
                #!!!!!!TODO retreive prefixe and create file name
                #output file graph
                file_name = prefix+'_plot.png'
                img_path = os.path.join(result_dir,file_name)
                #fig.savefig(img_path,dpi=50,bbox_extra_artists=([lgd,x_ticks]), bbox_inches='tight')
                fig.savefig(img_path,dpi=50)
                #savefig(img_path)
                htmlOut.write('<h4><b>Figure %s:</b> %s <br/></h4>' % (cpt,file_name.replace(".png","")))
                htmlOut.write('<img src="%s" /> ' % file_name)
                cpt = cpt+1
        
                close()
	
        ### create zip ###
        zipname = "plotCoverage_figures.zip"
        f = zipfile.ZipFile(zipname,'w',zipfile.ZIP_DEFLATED)
        for filename in glob.glob("*.png"):
            f.write(filename)
        f.close()
	
        htmlOut.write('''<li>
        <br/><a href="%s">%s</a><br>
        </li> '''  % (zipname, "Telecharger toutes les figures dans un ZIP"))
	
	

        return True


if __name__ == "__main__":
    PluginCLI(plotCoverage())
