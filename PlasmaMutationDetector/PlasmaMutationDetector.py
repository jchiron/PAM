#!/usr/bin/env python
# part of this code is from FileExporter plugin
import json
import os
import sys
import glob
from ion.plugin import *
from subprocess import *
from django.utils.datastructures import SortedDict
import urllib2
import httplib2
import requests
import re
import xlwt
import xlrd

class PlasmaMutationDetector(IonPlugin):
	""" Plugin object which launch R module PlasmaMutationDetector on remote host IB58B."""
	
	version = "0.1"
	
	envDict = dict(os.environ)
	json_dat = {}
	renameString = ""
	barcodeNames = []
	sampleNameLookup = {} # dictionary allowing us to get the sample name associated with a particular barcode name
	isBarcodedRun = False
	runName = 'unknown'
	project = 'unknown'
	runPath = ''
        projects = {}
        folder = '/media/PGM/results_S5/ADN_circulant_colon_lung/'
	
	doneSamples=[]
	
	# paths to the most recent, completed instance of the BackupCopy plugin
	#BackupCopy_Path = False
        	
	def getBarcodeNameFromFileName(self, fileName):
		for testBarcode in self.barcodeNames:
			testBarcode2 = testBarcode + '_'
			if testBarcode2 in fileName:
				barcodeName = testBarcode
				return barcodeName
		

	def getPluginPath(self, pluginName):
		pluginPath = False
		try:
			api_url = self.json_dat['runinfo']['api_url'] + '/v1/pluginresult/?format=json&plugin__name=%s&result=%s' % (pluginName,str(self.json_dat['runinfo']['pk']))
			api_key = self.json_dat['runinfo'].get('api_key', None)
			if api_key is not None:
				api_url = api_url + '&api_key=%s' % api_key
				#print 'Using API key: %s' % api_key
			else:
				print 'No API key available'
			f = urllib2.urlopen(api_url)
			d = json.loads(f.read())
			
			for plugin in d['objects']: # they are already sorted from newest to oldest, so take the first of this list which is 'completed'
				if plugin['state'] == 'Completed':
					pluginPath = plugin['path']
					pluginPath = pluginPath.split('/')[-1]
					break
			if not pluginPath:
				print 'WARNING! No completed instance found for plugin "%s"' % pluginName
		except:
			print 'ERROR!  Failed to get variant caller path'
		
		print 'INFO: using %s path: %s' % (pluginName,pluginPath)
		return pluginPath
		
       
	def launch(self):
		try:
			with open('startplugin.json', 'r') as fh:
				self.json_dat = json.load(fh)
		except:
			print 'Error reading plugin json.'
		
		try:
			htmlOut = open('%s/plugin_out/PlasmaMutationDetector_out/PlasmaMutationDetector_block.html'%self.envDict['ANALYSIS_DIR'], 'a+')
		except:
			htmlOut = open('PlasmaMutatioDetector_block.html', 'w')
		htmlOut.write('<html><body>\n')

		if os.path.isfile(self.envDict['TSP_FILEPATH_BARCODE_TXT']):
			self.isBarcodedRun = True
			
		if not self.isBarcodedRun:
			print "This plugin only work for barcoded run. \n Exiting..."
			sys.exit(0)

                # get the project name.
		#try:
		#	self.project = self.envDict['TSP_PROJECT']
		#except:
		#	self.project = 'unknown'
			
		# get the run name.
		try:
			self.runName = self.envDict['TSP_ANALYSIS_NAME']
		except:
			self.runName = 'unknown'

		self.runPath = self.project + '/' + self.runName

		runNumber = re.search('(?<=S5XL-0151-)(?P<Number>[0-9]+).*',self.runName)
		runNumber = runNumber.group('Number')
		print 'run name '+self.runName
		print 'run path '+self.runPath
		print 'run number '+runNumber

		runRemote = 'S5XL-'+runNumber
		#create PGM directory on remote host IB58B
		cmd = 'echo ionadmin | sudo -u ionadmin -S rsync -r -t --rsync-path=\"mkdir -p /u02/MiSeq/test/error_rate_ctDNA/'+runRemote+' && rsync\" /media/PGM/results_S5/BPER/touch.txt genetique@129.10.20.36:/u02/MiSeq/test/error_rate_ctDNA/'+runRemote+'/.'
		#cmd = 'echo ionadmin | sudo -u ionadmin -S rsync -r -t /media/PGM/results/ADN_circulant_colonlung/touch.txt genetique@ib58b:/u02/MiSeq/test/error_rate_ctDNA/'+runRemote+'/.'
		print 'create remote dir '+cmd

		os.system(cmd)


                # Get bam filenames.
                with open(os.path.join(self.json_dat['runinfo']['basecaller_dir'], 'datasets_basecaller.json'), 'r') as f:
                        json_basecaller = json.load(f)

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
                        #print 'adding BAM: %s' % tempPath
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
                                if sampleName not in self.doneSamples :
                                        self.barcodeNames.append(barcodeName)
                                        self.doneSamples.append(sampleName)
                                else :
                                        print 'PATIENT %s already analysed'%sampleName
                        else:
                                sampleName = 'None'
                                self.sampleNameLookup[barcodeName] = 'None' # makes it much easier later to do the lookup with no conditional tests
                                # MGD note: I considered setting blank sample names to the barcode name instead, but might not be what customer intended

                json_file = open(os.path.join(self.envDict['ANALYSIS_DIR'],'ion_params_00.json'), 'r')
                ion_params_json = json.load(json_file,parse_float=str)

                barcodedSamples_fromJson = ion_params_json['experimentAnalysisSettings'].get('barcodedSamples')
                if isinstance(barcodedSamples_fromJson,dict):
                        barcodedSamples = barcodedSamples_fromJson
                else:
                        barcodedSamples=json.loads(ion_params_json['experimentAnalysisSettings'].get('barcodedSamples'))

                target_path = {}
                for key in barcodedSamples :
                        #target = d['objects'][0]['store']['barcodes'][barcode]['targets_bed']
                        #target_path[barcode] = target
                        #print 'BARCODE : ' + barcode + ' TARGET : ' + target
                        all_barcodes= barcodedSamples.get(key).get('barcodes')
                        all_barcodes_str = ','.join(all_barcodes)
                        for barcode in all_barcodes :
                                bc_info = barcodedSamples.get(key).get('barcodeSampleInfo').get(barcode)
                                nucType=bc_info['nucleotideType']

                                if nucType == 'DNA' :
                                        #target = d['objects'][0]['store']['barcodes'][barcode]['targets_bed']
                                        target = bc_info['targetRegionBedFile']
                                        target_path[barcode] = target
                                        if target == '/results/uploads/BED/5/hg19/unmerged/detail/ColonLungV2.20140523.designed.with_NM.bed':
                                                self.projects[barcode] = 'ADN_circulant_colon_lung'
                                        elif target == '/results/uploads/BED/2/hg19/unmerged/detail/IAD118795_231_Designed_MET_000245.with_NM.bed':
                                                self.projects[barcode] = 'colon_lung_routine'
                                        elif target == '/results/uploads/BED/7/hg19/unmerged/detail/Oncomine_Focus.20150316.designed.plusPurpose.bed':
                                                self.projects[barcode] = 'essais_precoces'
                                        elif target == '/results/uploads/BED/6/hg19/unmerged/detail/IAD112641_173_designed.with_NM.bed':
                                                self.projects[barcode] = 'Safir02'
                                elif nucType =='RNA' :
                                        target = ''
                                        self.projects[barcode] = 'essais_precoces'

                                project = self.projects[barcode]
                                print 'BARCODE : ' + barcode + ' TARGET : ' + target + ' PROJECT : ' + project

                #todo check BackupCopy completed
                for barcode in self.projects.keys():
                        if (self.projects[barcode] == 'ADN_circulant_colon_lung'):
                                
                                #todo check BackupCopy completed
                                backupDir = '/media/PGM/results_S5/'+self.projects[barcode]+'/'+self.runName
                                print 'backup '+backupDir

                                print "BARCODE: %s" % barcode
                                #copy processed bam for all samples in runRemote directory
                                for name in glob.glob('%s/*/*processed*' % backupDir):
                                        print "bam file: %s" % name
                                        if barcode in name:
                                                print "bam file %s used for %s" % (name,barcode)
                                                cmd = 'echo ionadmin | sudo -u ionadmin -S rsync -r -t '+name+' genetique@129.10.20.36:/u02/MiSeq/test/error_rate_ctDNA/'+runRemote+'/.'
                                                print cmd
                                                os.system(cmd)

                #launch GATK on remote host IB58B
                cmd = 'echo ionadmin | sudo -u ionadmin -S rsync -r -t --rsync-path=\"cd /u02/MiSeq/test/error_rate_ctDNA/'+runRemote+'/ && bash /home/J.dupiot-chiron/Tools/PlasmaMutationDetector/GATK_indelrealign_baqr.sh /u02/MiSeq/test/error_rate_ctDNA/'+runRemote+'/ >/dev/null 2>&1 && rsync\" /media/PGM/results_S5/BPER/ADNtc_colon_lung_v2/touch.txt genetique@129.10.20.36:/u02/MiSeq/test/error_rate_ctDNA/'+runRemote+'/.'
                print cmd
                os.system(cmd)

		cmd = 'echo ionadmin | sudo -u ionadmin -S rsync -r -t --rsync-path=\"mkdir /u02/MiSeq/test/error_rate_ctDNA/'+runRemote+'/rBAM && rsync\" /media/PGM/results_S5/BPER/ADNtc_colon_lung_v2/touch.txt genetique@129.10.20.36:/u02/MiSeq/test/error_rate_ctDNA/'+runRemote+'/.'
		print cmd
		os.system(cmd)

		cmd = 'echo ionadmin | sudo -u ionadmin -S rsync -r -t --rsync-path=\"mv /u02/MiSeq/test/error_rate_ctDNA/'+runRemote+'/recal_* /u02/MiSeq/test/error_rate_ctDNA/'+runRemote+'/rBAM/. && rsync\" /media/PGM/results_S5/ADNtc_colon_lung_v2/touch.txt genetique@129.10.20.36:/u02/MiSeq/test/error_rate_ctDNA/'+runRemote+'/.'
		print cmd
		os.system(cmd)

		#launch R module PlasmaMutationDetector on remote host IB58B
		print 'launch R on remote host'
		cmd = 'echo ionadmin | sudo -u ionadmin -S rsync -r -t --rsync-path=\"Rscript --vanilla /u02/MiSeq/test/error_rate_ctDNA/PlasmaMutationDetector_on_remote.R '+runRemote+' >/dev/null 2>&1 && rsync\" /media/PGM/results_S5/ADNtc_colon_lung_v2/touch.txt genetique@129.10.20.36:/u02/MiSeq/test/error_rate_ctDNA/'+runRemote+'/.' 
		print cmd
		os.system(cmd)

		print 'copy result dirctory'
		#copy result on PGM directory
		cmd = 'echo ionadmin | sudo -u ionadmin -S rsync -r -t genetique@129.10.20.36:/u02/MiSeq/test/error_rate_ctDNA/'+runRemote+' /media/PGM/results_S5/BPER/ADNtc_colon_lung_v2/.' 
		print cmd
		os.system(cmd)



if __name__ == "__main__":
	PluginCLI(PlasmaMutationDetector())

