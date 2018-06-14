#!/usr/bin/env python
# part of this code is from FileExporter plugin
import json
import os
import sys
import glob
import re
from ion.plugin import *
from subprocess import *
from django.utils.datastructures import SortedDict
import urllib2
import requests

class BackupCopy(IonPlugin):
	""" Plugin object which copy all desired results file from the run into a folder."""
	""" A hierarchical filesystem is created (project->run->patients) to store the data, and files name are uniformly renamed"""
	version = "0.1"
	
	envDict = dict(os.environ)
	json_dat = {}
	renameString = ""
	barcodeNames = []
	sampleNameLookup = {} # dictionary allowing us to get the sample name associated with a particular barcode name
	isBarcodedRun = False
        projects = {}
        target = {}
        sample_type = {}
	runName = 'unknown'
	project = 'unknown'
	runPath = ''
	
	# paths to the most recent, completed instance of the plugins
	variantCallerPath = False
	variantAnnotationPath = False
	coverageAnalysisPath = False
	compareVcfPath = False
	plotCoveragePath = False
	HomopolymersPath = False
	check_contaminationPath = False
	finalReportPath = False
	run_Safir = False

	def getBarcodeNameFromFileName(self, fileName):
		for testBarcode in self.barcodeNames:
			testBarcode2 = testBarcode + '_'
			if testBarcode2 in fileName:
				barcodeName = testBarcode
		return barcodeName
		
	def createOneFolderByPatient(self, bamFileList):
		for fileName in bamFileList:
			fileName = fileName.replace('./', '')
			barcodeName = self.getBarcodeNameFromFileName(fileName)
			patientName = self.sampleNameLookup[barcodeName]
                        #print 'PATIENT NAME : ' + patientName + 'BARCODE NAME : ' + barcodeName
                        #for bc in self.projects.keys():
                        #        print 'BC : ' + bc + ' PROJECT : ' + self.projects[bc]
                        #if self.porjects[barcodeName]:
                        project = self.projects[barcodeName]
                        #else :
                        #        target = d['objects'][0]['store']['barcodes'][barcode]['targets_bed']
                                
                        print 'PATIENT NAME : ' + patientName + ' BARCODE NAME : ' + barcodeName + ' PROJECT : ' + project 
                        print 'RUN PATH : ' + self.runPath + ' RUN NAME : ' + self.runName
			print 'CREATING FOLDER : %s -'%(patientName)
			cmd = Popen(['%s/backup_pgm.sh' % (self.envDict['DIRNAME']), 'mkdir', project + '/' + self.runName + '/' + patientName], stdout=PIPE)
			out, err = cmd.communicate()
			print 'OUT: %s\nERR: %s'%(out, err)
			
	def getPluginPath(self, pluginName):
		pluginPath = False
                target_path = {} #hash key=barcode value=target to get project for specific barcode
		try:
			api_url = self.json_dat['runinfo']['api_url'] + '/v1/pluginresult/?format=json&plugin__name=%s&result=%s' % (pluginName,str(self.json_dat['runinfo']['pk']))
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
					pluginPath = plugin['path']
					pluginPath = pluginPath.split('/')[-1]
					break
			if not pluginPath:
				print 'WARNING! No completed instance found for plugin "%s"' % pluginName
		except:
			print 'ERROR!  Failed to get variant caller path'
		
		print 'INFO: using %s path: %s' % (pluginName,pluginPath)
		return pluginPath
		

	# Method to rename and copy the .bam and .bai files.
	def copyBam(self, bamFileList):
		for fileName in bamFileList:
			fileName = fileName.replace('./', '')
			barcodeName = self.getBarcodeNameFromFileName(fileName)
                        project = self.projects[barcodeName]

			# build our new filename
			finalName = self.renameString.replace('@BARINFO@', barcodeName)
			finalName = finalName.replace('@SAMPLEID@', self.sampleNameLookup[barcodeName])
                        if barcodeName in self.sample_type.keys() :#for barcodes in plan sample_type is defined 
                                if self.sample_type[barcodeName] == 'DNA':
                                        destName = finalName + '.bam'
                                elif self.sample_type[barcodeName] ==  'RNA':
                                        destName = finalName + '_RNA_rawlib.bam'
                        else: #no sample type defined for barcode, it is to manage barcodes wtih no name associated found in sequencing (barcode switch or contamination)
                                self.sample_type[barcodeName] = 'DNA'
                                destName = finalName + '.bam'

                        

			# And, copy.
			print 'COPYING BAM: %s --> %s'%(fileName, destName)
			cmdBam = Popen(['%s/backup_pgm.sh' % (self.envDict['DIRNAME']), 'cp', fileName, '%s/%s/%s/%s' % (project,self.runName, self.sampleNameLookup[barcodeName], destName)], stdout=PIPE)
			outBam, errBam = cmdBam.communicate()
			print 'OUT: %s\nERR: %s'%(outBam, errBam)
			
			fileName = fileName.replace('.bam', '.bam.bai')
			destName = destName.replace('.bam', '.bam.bai')
			
			print 'COPYING BAI: %s --> %s'%(fileName, destName)
			cmdBai = Popen(['%s/backup_pgm.sh' % (self.envDict['DIRNAME']), 'cp', fileName, '%s/%s/%s/%s' % (project,self.runName, self.sampleNameLookup[barcodeName], destName)], stdout=PIPE)
			outBai, errBai = cmdBai.communicate()
			print 'OUT: %s\nERR: %s'%(outBai, errBai)

			# 02/10/2015 add processed bam if variantCaller exists (in variantCaller plugin folder)
			if self.variantCallerPath:	# if == False, no instance of variantCaller found
                                if self.sample_type[barcodeName] == 'DNA' : #variantCaller done for DNA samples not for RNA, so only processed bam for DNA samples
                                        instance = self.variantCallerPath.split(".")[-1]
                                        srcNameProcessedBam = '%s/plugin_out/%s/%s/%s_rawlib_processed.bam' % (self.envDict['ANALYSIS_DIR'], self.variantCallerPath, barcodeName, barcodeName)
                                        srcNameProcessedBai = '%s/plugin_out/%s/%s/%s_rawlib_processed.bam.bai' % (self.envDict['ANALYSIS_DIR'], self.variantCallerPath, barcodeName, barcodeName)
				
                                        srcNameProcessedBam = '%s/plugin_out/%s/%s/%s_rawlib.realigned_processed.bam' % (self.envDict['ANALYSIS_DIR'], self.variantCallerPath, barcodeName, barcodeName)
                                        srcNameProcessedBai = '%s/plugin_out/%s/%s/%s_rawlib.realigned_processed.bam.bai' % (self.envDict['ANALYSIS_DIR'], self.variantCallerPath, barcodeName, barcodeName)
                                
                                        destNameProcessedBam = finalName + '.processed_%s.bam' % instance
                                        destNameProcessedBai = finalName + '.processed_%s.bam.bai' % instance
                                        # And, copy.
                                        print 'COPYING PROCESSED BAM : %s --> %s'%(srcNameProcessedBam, destNameProcessedBam)
                                        cmd = Popen(['%s/backup_pgm.sh' % (self.envDict['DIRNAME']), 'cp', srcNameProcessedBam, '%s/%s/%s/%s' % (project,self.runName, self.sampleNameLookup[barcodeName], destNameProcessedBam)], stdout=PIPE)
                                        out, err = cmd.communicate()
                                        print 'OUT: %s\nERR: %s'%(out, err)
                                        
                                        print 'COPYING PROCESSED BAI : %s --> %s'%(srcNameProcessedBai, destNameProcessedBai)
                                        cmd = Popen(['%s/backup_pgm.sh' % (self.envDict['DIRNAME']), 'cp', srcNameProcessedBai, '%s/%s/%s/%s' % (project,self.runName, self.sampleNameLookup[barcodeName], destNameProcessedBai)], stdout=PIPE)
                                        out, err = cmd.communicate()
                                        print 'OUT: %s\nERR: %s'%(out, err)
                                #elif self.sample_type[barcodeName] == 'RNA' : #variantCaller not done for RNA samples, get rawlib.bam from run folder
                                #        srcNameRawlibBam = '%s/%s_rawlib.bam' % (self.envDict['ANALYSIS_DIR'], barcodeName)
                                #        destNameRawlibBam = 

	# Method to rename and copy the .vcf and .xls files.
	def copyVcf(self, bamFileList):
		instance = self.variantCallerPath.split(".")[-1]
		for fileName in bamFileList:
			fileName = fileName.replace('./', '')
			barcodeName = self.getBarcodeNameFromFileName(fileName)
                        project = self.projects[barcodeName]

			# build our new filename and move this file
			finalName = self.renameString.replace('@BARINFO@', barcodeName)
			finalName = finalName.replace('@SAMPLEID@', self.sampleNameLookup[barcodeName])

			destNameVcf = finalName + '.variantcaller.%s.vcf' % instance
			destNameXls = finalName + '.variantcaller.%s.xls' % instance
                        destNameFiltered = finalName + '.variantcaller.%s.nocall.vcf' % instance

			srcNameVcf = '%s/plugin_out/%s/%s/TSVC_variants.vcf' % (self.envDict['ANALYSIS_DIR'], self.variantCallerPath, barcodeName)
			srcNameXls = '%s/plugin_out/%s/%s/alleles_%s.xls' % (self.envDict['ANALYSIS_DIR'], self.variantCallerPath, barcodeName, barcodeName)
                        srcNameFiltered = '%s/plugin_out/%s/%s/small_variants_filtered.vcf' % (self.envDict['ANALYSIS_DIR'], self.variantCallerPath, barcodeName)

			#test added to take account of mixed panels on same run (VaraintCaller results for all samples are not in same folder) or RNA samples haven't VC results
			if os.path.isfile(srcNameVcf) and os.access(srcNameVcf, os.R_OK):
				# And, copy.
				print 'COPYING VCF FILE : %s --> %s'%(srcNameVcf, destNameVcf)
				cmdVcf = Popen(['%s/backup_pgm.sh' % (self.envDict['DIRNAME']), 'cp', srcNameVcf, '%s/%s/%s/%s' % (project,self.runName, self.sampleNameLookup[barcodeName], destNameVcf)], stdout=PIPE)
				outVcf, errVcf = cmdVcf.communicate()
				print 'OUT: %s\nERR: %s'%(outVcf, errVcf)
			
				print 'COPYING XLS FILE : %s --> %s'%(srcNameXls, destNameXls)
				cmdXls = Popen(['%s/backup_pgm.sh' % (self.envDict['DIRNAME']), 'cp', srcNameXls, '%s/%s/%s/%s' % (project,self.runName, self.sampleNameLookup[barcodeName], destNameXls)], stdout=PIPE)
				outXls, errXls = cmdXls.communicate()
				print 'OUT: %s\nERR: %s'%(outXls, errXls)
                                
                                print 'COPYING NOCALL VCF FILE : %s --> %s'%(srcNameFiltered, destNameFiltered)
                                cmdNocall = Popen(['%s/backup_pgm.sh' % (self.envDict['DIRNAME']), 'cp', srcNameFiltered, '%s/%s/%s/%s' % (project,self.runName, self.sampleNameLookup[barcodeName], destNameFiltered)], stdout=PIPE)
                                outNocall, errNocall = cmdNocall.communicate()
                                print 'OUT: %s\nERR: %s'%(outNocall, errNocall)
			else: 
				print 'OUT: %s does not exist\nNo results found for this sample in this variantCaller instance.' % (srcNameXls)
                
                for project in self.projects.values():
                        #also copy vc parameter settings
                        if project == 'ADN_circulant_colon_lung':
                                srcNameJson = '%s/plugin_out/%s/ADN circulant colon lung_parameters.json' % (self.envDict['ANALYSIS_DIR'], self.variantCallerPath)
                                destNameJson = 'variantcaller_parameters.%s.json' % instance
                        elif project == 'colon_lung_routine':
                                srcNameJson = '%s/plugin_out/%s/Colon lung MET14_parameters.json' % (self.envDict['ANALYSIS_DIR'], self.variantCallerPath)
                                destNameJson = 'variantcaller_parameters.%s.json' % instance
                        elif project == 'essais_precoces':
                                srcNameJson = '%s/plugin_out/%s/Oncomine focus_parameters.json' % (self.envDict['ANALYSIS_DIR'], self.variantCallerPath)
                                destNameJson = 'variantcaller_parameters.%s.json' % instance
                        elif project == 'Safir02':
                                srcNameJson = '%s/plugin_out/%s/Safir02_parameters.json' % (self.envDict['ANALYSIS_DIR'], self.variantCallerPath)
                                destNameJson = 'variantcaller_parameters.%s.json' % instance
                        elif project == 'BRCA_qiagen':
                                srcNameJson = '%s/plugin_out/%s/BRCA_qiagen_TSS-5.6_parameters.json' % (self.envDict['ANALYSIS_DIR'], self.variantCallerPath)
                                destNameJson = 'variantcaller_parameters.%s.json' % instance
                        else :
                                #project not yet defined maby first run
                                print 'NO JSON PARAMETERS FILE FIND FOR PROJECT %' % project
                                
                        print 'COPYING JSON PARAMETERS FILE : %s --> %s'%(srcNameJson, destNameJson)
                        cmdJson = Popen(['%s/backup_pgm.sh' % (self.envDict['DIRNAME']), 'cp', srcNameJson, '%s/%s/%s' % (project,self.runName, destNameJson)], stdout=PIPE)
                        outJson, errJson = cmdJson.communicate()
                        print 'OUT: %s\nERR: %s'%(outJson, errJson)
			
			
	def copyTsv(self):
		print '!!!!!!!!!!!!!!!!! IN copyTSV def :'
		#print '%s/plugin_out/%s/IonXpress_*/CHR/*_NGS_Diag_VC_TS.tsv' % (self.envDict['ANALYSIS_DIR'],self.variantAnnotationPath)
                api_url = self.json_dat['runinfo']['api_url'] + '/v1/results/%s/?format=json' % str(self.json_dat['runinfo']['pk'])
                print 'API URL: ' + api_url
                f = urllib2.urlopen(api_url)
                d = json.loads(f.read())
                for project in d['projects']:
                        num = int(project.split('/')[-2])
                        print "Project number found :"
                        print num
                        if num == 5: #project = Safir2
                                self.run_Safir = True

                #if self.run_Safir :
                for fileName in glob.glob('%s/plugin_out/%s/*/CHR/*_NGS_VC.tsv' % (self.envDict['ANALYSIS_DIR'],self.variantAnnotationPath)):
                        print 'VA PATH: %s' %(self.variantAnnotationPath)
                        srcName = fileName
                        print 'TSV : srcName'
                        barcodeName = fileName.split("/")[-3]
                        project = self.projects[barcodeName]
                        # build our new filename and move this file                                        
                        finalName = self.renameString.replace('@BARINFO@', barcodeName)
                        finalName = finalName.replace('@SAMPLEID@', self.sampleNameLookup[barcodeName])
                        instance = self.variantAnnotationPath.split(".")[-1]
                        destName = finalName + '.annotation.%s.tsv' % instance
                        
                        #test added to take account of mixed panels on same run (VaraintAnnotation results for all samples are not in same folder) or RNA samples haven't VA results                                  
                        if os.path.isfile(srcName) and os.access(srcName, os.R_OK):

                                # And, copy.                                                              
                                print 'COPYING TSV: %s --> %s'%(srcName, destName)
                                cmd = Popen(['%s/backup_pgm.sh' % (self.envDict['DIRNAME']), 'cp', srcName, '%s/%s/%s/%s' % (project,self.runName, self.sampleNameLookup[barcodeName], destName)], stdout=PIPE)
                                out, err = cmd.communicate()
                                print 'OUT: %s\nERR: %s'%(out, err)
                        else:
                                print 'OUT: %s does not exist\nNo results found for this sample in this variantAnnotation instance.' % (srcName)


                #else:
                for fileName in glob.glob('%s/plugin_out/%s/*/CHR/*_NGS_Diag_VC_TS.tsv' % (self.envDict['ANALYSIS_DIR'],self.variantAnnotationPath)):
                        print 'VA PATH: %s' %(self.variantAnnotationPath)
                        srcName = fileName
                        print 'TSV : srcName'
                        barcodeName = fileName.split("/")[-3]
                        project = self.projects[barcodeName]
                        # build our new filename and move this file
                        finalName = self.renameString.replace('@BARINFO@', barcodeName)
                        finalName = finalName.replace('@SAMPLEID@', self.sampleNameLookup[barcodeName])
                        instance = self.variantAnnotationPath.split(".")[-1]
                        destName = finalName + '.annotation.%s.tsv' % instance
			
                        #test added to take account of mixed panels on same run (VaraintAnnotation results for all samples are not in same folder) or RNA samples haven't VA results
                        if os.path.isfile(srcName) and os.access(srcName, os.R_OK):
																     
                                # And, copy.
                                print 'COPYING TSV: %s --> %s'%(srcName, destName)
                                cmd = Popen(['%s/backup_pgm.sh' % (self.envDict['DIRNAME']), 'cp', srcName, '%s/%s/%s/%s' % (project,self.runName, self.sampleNameLookup[barcodeName], destName)], stdout=PIPE)
                                out, err = cmd.communicate()
                                print 'OUT: %s\nERR: %s'%(out, err)
                        else:
                                print 'OUT: %s does not exist\nNo results found for this sample in this variantAnnotation instance.' % (srcName)


	def copyCov(self):	
		for fileName in glob.glob('%s/plugin_out/%s/*/*.amplicon.cov.xls' % (self.envDict['ANALYSIS_DIR'],self.coverageAnalysisPath)):
			srcName = fileName
			barcodeName = fileName.split("/")[-2]
                        project = self.projects[barcodeName]
			# build our new filename and move this file
			finalName = self.renameString.replace('@BARINFO@', barcodeName)
			finalName = finalName.replace('@SAMPLEID@', self.sampleNameLookup[barcodeName])
			instance = self.coverageAnalysisPath.split(".")[-1]
			destName = finalName + '.amplicon.cov.%s.xls' % instance

			# And, copy.
			print 'COPYING COVERAGE ANALYSIS FILE : %s --> %s'%(srcName, destName)
			cmd = Popen(['%s/backup_pgm.sh' % (self.envDict['DIRNAME']), 'cp', srcName, '%s/%s/%s/%s' % (project,self.runName, self.sampleNameLookup[barcodeName], destName)], stdout=PIPE)
			out, err = cmd.communicate()
			print 'OUT: %s\nERR: %s'%(out, err)
			
		for fileName in glob.glob('%s/plugin_out/%s/*.bc_summary.xls' % (self.envDict['ANALYSIS_DIR'],self.coverageAnalysisPath)):
			srcName = fileName
			analysisName = fileName.split(".")[0]
			# build our new filename and move this file
			runName = self.envDict['ANALYSIS_DIR'].split("/")[-1]
                        #finalName = self.renameString.replace(runName, analysisName)
			instance = self.coverageAnalysisPath.split(".")[-1]
			destName = runName + '.%s.bc_summary.xls' %instance

                        for project in self.projects.values():
                                # And, copy.
                                print 'COPYING COVERAGE ANALYSIS FILE : %s --> %s'%(srcName, destName)
                                cmd = Popen(['%s/backup_pgm.sh' % (self.envDict['DIRNAME']), 'cp', srcName, '%s/%s/%s' % (project,self.runName, destName)], stdout=PIPE)
                                out, err = cmd.communicate()
                                print 'OUT: %s\nERR: %s'%(out, err)


			
	def copyCom(self):
		for fileName in glob.glob('%s/plugin_out/downloads/Compare_VCF_*_*.csv' % self.envDict['ANALYSIS_DIR']):
			srcName = fileName
			destName = fileName.split("/")[-1]

			# And, copy.
			print 'COPYING COVERAGE ANALYSIS FILE : %s --> %s'%(srcName, destName)
			cmd = Popen(['%s/backup_pgm.sh' % (self.envDict['DIRNAME']), 'cp', srcName, '%s/%s' % (self.runName, destName)], stdout=PIPE)
			out, err = cmd.communicate()
			print 'OUT: %s\nERR: %s'%(out, err)
			
	def copyPlot(self):
                for project in self.projects.values():
                        # Create a folder for plots
                        print 'CREATING PLOT FOLDER'
                        cmd = Popen(['%s/backup_pgm.sh' % (self.envDict['DIRNAME']), 'mkdir', project + '/' + self.runName  + '/_plotCoverage'], stdout=PIPE)
                        out, err = cmd.communicate()
                        print 'OUT: %s\nERR: %s'%(out, err)
		
		fileNames = []
		for fileName in glob.glob('%s/plugin_out/%s/*.png' % (self.envDict['ANALYSIS_DIR'],self.plotCoveragePath)):
			fileNames.append(fileName)
		for fileName in glob.glob('%s/plugin_out/%s/*_cov_less_than_300x.bed' % (self.envDict['ANALYSIS_DIR'],self.plotCoveragePath)):
			fileNames.append(fileName)
		for fileName in glob.glob('%s/plugin_out/%s/all_samples_amplicons_log10.cov.tsv' % (self.envDict['ANALYSIS_DIR'],self.plotCoveragePath)):
			fileNames.append(fileName)
		for fileName in glob.glob('%s/plugin_out/%s/amplicon.uncov.tsv' % (self.envDict['ANALYSIS_DIR'],self.plotCoveragePath)):
			fileNames.append(fileName)
		
		for fileName in fileNames:
			srcName = fileName
			destName = fileName.split("/")[-1]

                        for project in self.projects.values():
                                # And, copy.
                                print 'COPYING PLOT COVERAGE FILE : %s --> %s'%(srcName, destName)
                                cmd = Popen(['%s/backup_pgm.sh' % (self.envDict['DIRNAME']), 'cp', srcName, '%s/%s/_plotCoverage/%s' % (project, self.runName, destName)], stdout=PIPE)
                                out, err = cmd.communicate()
                                print 'OUT: %s\nERR: %s'%(out, err)
			
	def copyHomopoly(self):
		# Create a folder for Homopolymers
		print 'CREATING Homopolymers FOLDER'
		cmd = Popen(['%s/backup_pgm.sh' % (self.envDict['DIRNAME']), 'mkdir', self.runPath + '/_Homopolymers'], stdout=PIPE)
		out, err = cmd.communicate()
		print 'OUT: %s\nERR: %s'%(out, err)
		
		#srcName1 = "%s/plugin_out/%s/BRCA2_c.2175dupA.png" % (self.envDict['ANALYSIS_DIR'],self.HomopolymersPath)
		#destName1 = srcName1.split("/")[-1]
		srcName2 = "%s/plugin_out/%s/homopolymersReport.xlsx" % (self.envDict['ANALYSIS_DIR'],self.HomopolymersPath)
		#destName2 = srcName2.split("/")[-1]

                try:
                        self.runName = self.envDict['TSP_ANALYSIS_NAME']
                except:
                        self.runName = 'unknown'
                runNumber = re.search('(?<=S5XL-0151-)(?P<Number>[0-9]+).*',self.runName)
                runNumber = runNumber.group('Number')
                destName2 = 'homopolymersReport_S5XL-0151-'+runNumber+'.xlsx'

		# And, copy.
		print 'COPYING HOMOPOLYMERS FILES :'
		#cmd = Popen(['%s/backup_pgm.sh' % (self.envDict['DIRNAME']), 'cp', srcName1, '%s/_Homopolymers/%s' % (self.runPath, destName1)], stdout=PIPE)
		#out, err = cmd.communicate()
		#print 'OUT: %s\nERR: %s'%(out, err)
		cmd = Popen(['%s/backup_pgm.sh' % (self.envDict['DIRNAME']), 'cp', srcName2, '%s/_Homopolymers/%s' % (self.runPath, destName2)], stdout=PIPE)
		out, err = cmd.communicate()
		print 'OUT: %s\nERR: %s'%(out, err)
		
	def copyCheck(self):
                for project in self.projects.values():
                        # Create a folder for Check_contamination
                        print 'CREATING Check_contamination FOLDER'
                        cmd = Popen(['%s/backup_pgm.sh' % (self.envDict['DIRNAME']), 'mkdir', project + '/' + self.runName + '/_Check-contamination'], stdout=PIPE)
                        out, err = cmd.communicate()
                        print 'OUT: %s\nERR: %s'%(out, err)
		
		fileNames = []
		for fileName in glob.glob("%s/plugin_out/*/Check_contamination_*.xls" % (self.envDict['ANALYSIS_DIR'])):
                        #print 'CHECK CONTA FILE : %s' % fileName
			fileNames.append(fileName)

                for project in self.projects.values():
                        # And, copy.
                        print 'COPYING CHECK CONTAMINATION FILES :'
                        for fileName in fileNames:
                                cmd = Popen(['%s/backup_pgm.sh' % (self.envDict['DIRNAME']), 'cp', fileName, '%s/%s/_Check-contamination/' %(project,self.runName)], stdout=PIPE)
                                out, err = cmd.communicate()
                                print 'OUT: %s\nERR: %s'%(out, err)
			
	def copyFin(self):	
		for fileName in glob.glob('%s/plugin_out/%s/*/*.finalReport.xls' % (self.envDict['ANALYSIS_DIR'],self.finalReportPath)):
			srcName = fileName
			barcodeName = fileName.split("/")[-2]
                        project = self.projects[barcodeName]
			try :
				sampl_name = self.sampleNameLookup[barcodeName]
			# build our new filename and move this file
				instance = self.finalReportPath.split(".")[-1]
				destName = sampl_name + "_" + barcodeName + '.finalReport.%s.xls' % instance

				if os.path.isfile(srcName) and os.access(srcName, os.R_OK):
					# And, copy.
					print 'COPYING FINALREPORT : %s --> %s'%(srcName, destName)
					cmd = Popen(['%s/backup_pgm.sh' % (self.envDict['DIRNAME']), 'cp', srcName, '%s/%s/%s/%s' % (project,self.runName, self.sampleNameLookup[barcodeName], destName)], stdout=PIPE)
					out, err = cmd.communicate()
					print 'OUT: %s\nERR: %s'%(out, err)
				else:
					print 'OUT: %s does not exist\nNo results found for this sample in this FinalReport instance.' % (srcName)
			except:
				print 'No finalReport found for sample %s' % (barcodeName)

			

	def copyFastqSff(self, suffix):
		for fileName in glob.glob('*.%s' % suffix):
			print 'DEBUG: checking %s file: %s' % (suffix, fileName)
			barcodeName = self.getBarcodeNameFromFileName(fileName)
                        project = self.projects[barcodeName]
			finalName = self.renameString.replace('@BARINFO@', barcodeName)
			finalName = finalName.replace('@SAMPLEID@', self.sampleNameLookup[barcodeName])
			destName = finalName + '.' + suffix

			print 'COPYING %s FILE : %s to %s' % (suffix, fileName, destName)
			cmd = Popen(['%s/backup_pgm.sh' % (self.envDict['DIRNAME']), 'cp', fileName, '%s/%s/%s/%s' % (project,self.runName, self.sampleNameLookup[barcodeName], destName)], stdout=PIPE)
			out, err = cmd.communicate()
			print 'OUT: %s\nERR: %s'%(out, err)

	def launch(self):
		try:
			with open('startplugin.json', 'r') as fh:
				self.json_dat = json.load(fh)
		except:
			print 'Error reading plugin json.'
		
		try:
			htmlOut = open('%s/plugin_out/BackupCopy_out/BackupCopy_block.html'%self.envDict['ANALYSIS_DIR'], 'a+')
		except:
			htmlOut = open('BackupCopy_block.html', 'w')
		htmlOut.write('<html><body>\n')

		if os.path.isfile(self.envDict['TSP_FILEPATH_BARCODE_TXT']):
			self.isBarcodedRun = True
			
		if not self.isBarcodedRun:
			print "This plugin only work for barcoded run. \n Exiting..."
			sys.exit(0)

		# set defaults for user options
		bamCreate = False
		sffCreate = False
		fastqCreate = False
		vcCreate = False
		tsvCreate = False
		covCreate = False
		pdfCreate = False
		logCreate = False
		#comCreate = False
		plotCreate = False
		homoCreate = False
		cheCreate = False
		finCreate = False

		# Parse pluginconfig json.
		try:
			delim = self.json_dat['pluginconfig']['delimiter_select']
			selections = self.json_dat['pluginconfig']['select_dialog']
			try:
				temp = self.json_dat['pluginconfig']['bamCreate']
				if (temp == 'on'):
					bamCreate = True
			except:
				print 'Logged: no BAM creation.'
			try:
				temp = self.json_dat['pluginconfig']['sffCreate']
				if (temp == 'on'):
					sffCreate = True
			except:
				print 'Logged: no SFF creation.'
			try:
				temp = self.json_dat['pluginconfig']['fastqCreate']
				if (temp == 'on'):
					fastqCreate = True
			except:
				print 'Logged: no FASTQ creation.'
			try:
				temp = self.json_dat['pluginconfig']['vcCreate']
				if (temp == 'on'):
					vcCreate = True
			except:
				print 'Logged: no VC linking.'
			try:
				temp = self.json_dat['pluginconfig']['tsvCreate']
				if (temp == 'on'):
					tsvCreate = True
			except:
				print 'Logged: no TSV creation.'
			try:
				temp = self.json_dat['pluginconfig']['covCreate']
				if (temp == 'on'):
					covCreate = True
			except:
				print 'Logged: no COV creation.'
			try:
				temp = self.json_dat['pluginconfig']['pdfCreate']
				if (temp == 'on'):
					pdfCreate = True
			except:
				print 'Logged: no PDF creation.'
			try:
				temp = self.json_dat['pluginconfig']['logCreate']
				if (temp == 'on'):
					logCreate = True
			except:
				print 'Logged: no LOG creation.'
			#try:
			#	temp = self.json_dat['pluginconfig']['comCreate']
			#	if (temp == 'on'):
			#		comCreate = True
			#except:
			#	print 'Logged: no COM creation.'
			try:
				temp = self.json_dat['pluginconfig']['plotCreate']
				if (temp == 'on'):
					plotCreate = True
			except:
				print 'Logged: no PLOT creation.'
			try:
				temp = self.json_dat['pluginconfig']['homoCreate']
				if (temp == 'on'):
					homoCreate = True
			except:
				print 'Logged: no HOMOPOLYMERS creation.'
			try:
				temp = self.json_dat['pluginconfig']['cheCreate']
				if (temp == 'on'):
					cheCreate = True
			except:
				print 'Logged: no CHECKCONTA creation.'
			try:
				temp = self.json_dat['pluginconfig']['finCreate']
				if (temp == 'on'):
					finCreate = True
			except:
				print 'Logged: no FINALREPORT creation.'
		except:
			print 'Warning: plugin does not appear to be configured, will default to run name'
			delim = '.'
			selections = ['TSP_RUN_NAME']
			fastqCreate = True
		try:
			self.runlevel = self.json_dat['runplugin']['runlevel']
		except:
			self.runlevel = ""
			print 'No run level detected.'

		# DEBUG: Print barcoded sampleID data.
		#test if/else for version 4.6
		barcodedSamples_fromJson = self.json_dat['plan']['barcodedSamples']
		if isinstance(barcodedSamples_fromJson,dict):
			samples = barcodedSamples_fromJson
		else:
			samples = json.loads(self.json_dat['plan']['barcodedSamples'])

		#print 'SAMPLEID DATA: %s'%samples
		#print 'TYPE: %s' % type(samples)
		
		htmlOut.write('<b>Copy BAM / BAI?</b> %s<br/>\n'%bamCreate)
		htmlOut.write('<b>Create SFF?</b> %s<br/>\n'%sffCreate)
		htmlOut.write('<b>Create FASTQ?</b> %s<br/>\n'%fastqCreate)
		htmlOut.write('<b>Copy variant caller file(s)?</b> %s<br/>\n'%vcCreate)
		htmlOut.write('<b>Copy variant annotation file(s)?</b> %s<br/>\n'%tsvCreate)
		htmlOut.write('<b>Copy coverage analysis files(?)</b> %s<br/>\n'%covCreate)
		#htmlOut.write('<b>Copy compare vcf files(?)</b> %s<br/>\n'%comCreate)
		htmlOut.write('<b>Copy plot coverage files(?)</b> %s<br/>\n'%plotCreate)
		htmlOut.write('<b>Copy homopolymers files(?)</b> %s<br/>\n'%homoCreate)
		htmlOut.write('<b>Copy check_contamination files(?)</b> %s<br/>\n'%cheCreate)
		htmlOut.write('<b>Copy finalReport(?)</b> %s<br/>\n'%finCreate)
		htmlOut.write('<b>Copy PDF report?</b> %s<br/>\n'%pdfCreate)
		htmlOut.write('<b>Copy log files(?)</b> %s<br/>\n'%logCreate)

		# Remove empty values.
		if not isinstance(selections, unicode):
			selections[:] = [entry for entry in selections if entry != '']
		elif selections != u'':
			selections = [selections]
		else:
			print 'Warning: No options selected, will use default TSP_RUN_NAME'
			selections = ['TSP_RUN_NAME']
		
		# Get appropriate values.
		for i in range(len(selections)):
			# Use an arbitrary value that nobody will ever use otherwise, so they're easy to replace.
			# '@' is an invalid character, right? Maybe not, actually...
			if (selections[i] == 'OPT_BARCODE'):
				selections[i] = '@BARINFO@'
			elif (selections[i] == 'TSP_SAMPLE'):
				selections[i] = '@SAMPLEID@'
			else:
				selections[i] = self.envDict[selections[i]]
				selections[i] = selections[i].replace('\\', '') # no idea why, but some chips look like \314R\ with backslashes?
		try:
			reference_path = self.envDict['TSP_FILEPATH_GENOME_FASTA']
		except:
			reference_path = ''

		# won't make sense to create vcf links if no reference was specified, so don't waste the time
		if reference_path == '':
			vcCreate = False

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
		#self.runPath = self.project + '/' + self.runName
		
                # get target and project for barcode
                target_path = {}
                #projects = {}
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
                if not variantCaller_path:
                        print 'WARNING! No completed instance found for plugin coverageAnalysis'

                                

                try:
                        json_file = open(os.path.join(self.envDict['ANALYSIS_DIR'],'ion_params_00.json'), 'r')
                        ion_params_json = json.load(json_file,parse_float=str)
                        json_file.close()
                except:
                        print('ERROR: Failed to load and parse ion_params_00.json')
                        return 1
                #~ barcodedSamples=ast.literal_eval(ion_params_json['experimentAnalysisSettings'].get('barcodedSamples'))
                barcodedSamples_fromJson = ion_params_json['experimentAnalysisSettings'].get('barcodedSamples')
                if isinstance(barcodedSamples_fromJson,dict):
                        barcodedSamples = barcodedSamples_fromJson
                else:
                        barcodedSamples=json.loads(ion_params_json['experimentAnalysisSettings'].get('barcodedSamples'))

                #for barcode in d['objects'][0]['store']['barcodes']:
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
                                        self.sample_type[barcode] = 'DNA'
                                        if target == '/results/uploads/BED/5/hg19/unmerged/detail/ColonLungV2.20140523.designed.with_NM.bed':
                                                self.projects[barcode] = 'ADN_circulant_colon_lung'
                                        elif target == '/results/uploads/BED/2/hg19/unmerged/detail/IAD118795_231_Designed_MET_000245.with_NM.bed':
                                                self.projects[barcode] = 'colon_lung_routine'
                                        elif target == '/results/uploads/BED/7/hg19/unmerged/detail/Oncomine_Focus.20150316.designed.plusPurpose.bed':
                                                self.projects[barcode] = 'essais_precoces'
                                        elif target == '/results/uploads/BED/6/hg19/unmerged/detail/IAD112641_173_designed.with_NM.bed':
                                                self.projects[barcode] = 'Safir02'
                                        elif target == '/results/uploads/BED/18/hg19/unmerged/detail/NGHS-102X_BRCA1BRCA2.amplicons.with_NM.bed':
                                                self.projects[barcode] = 'BRCA_qiagen'
                                elif nucType =='RNA' :
                                        self.sample_type[barcode] = 'RNA'
                                        target = ''
                                        self.projects[barcode] = 'essais_precoces'

                                project = self.projects[barcode]
                                print 'BARCODE : ' + barcode + ' TARGET : ' + target + ' PROJECT : ' + project

                        # creating folders of projects
                        print 'CREATING FOLDER : %s -'%(project)
                        cmd = Popen([self.envDict['DIRNAME']+'/backup_pgm.sh', 'mkdir',project], stdout=PIPE)
                        out, err = cmd.communicate()
                        print 'OUT: %s\nERR: %s'%(out, err)
		
                        self.runPath = project + '/' + self.runName
                        # creating folders of the run in project
                        print 'CREATING FOLDER : %s -'%(self.runName)
                        cmd = Popen([self.envDict['DIRNAME']+'/backup_pgm.sh', 'mkdir', self.runPath], stdout=PIPE)
                        out, err = cmd.communicate()
                        print 'OUT: %s\nERR: %s'%(out, err)		
		
		# get the actual path to the plugins
		if vcCreate:
			self.variantCallerPath = self.getPluginPath("variantCaller")
			htmlOut.write('&#8594; VariantCaller instance : <b>%s</b><br/>\n'%self.variantCallerPath)
			if not self.variantCallerPath:
				vcCreate = False
		if tsvCreate:
			self.variantAnnotationPath = self.getPluginPath("variantAnnotation")
			htmlOut.write('&#8594; VariantAnnotation instance : <b>%s</b><br/>\n'%self.variantAnnotationPath)
			if not self.variantAnnotationPath:
				tsvCreate = False
		if covCreate:
			self.coverageAnalysisPath = self.getPluginPath("coverageAnalysis")
			htmlOut.write('&#8594; coverageAnalysis instance : <b>%s</b><br/>\n'%self.coverageAnalysisPath)
			if not self.coverageAnalysisPath:
				covCreate = False	
		#if comCreate: # for this option, all compare_vcf files will be copied so the getting the last instance is not necessary. 
		#	self.compareVcfPath = self.getPluginPath("Compare_VCF")# though, using getPluginPath() also allow us to check if one instance exist (if not not, the download folder may not exist).
		#	if not self.compareVcfPath:
		#		comCreate = False
		if plotCreate:
			self.plotCoveragePath = self.getPluginPath("plotCoverage")
			htmlOut.write('&#8594; plotCoverage instance : <b>%s</b><br/>\n'%self.plotCoveragePath)
			if not self.plotCoveragePath:
				plotCreate = False	
		if homoCreate:
			self.HomopolymersPath = self.getPluginPath("Homopolymers")
			htmlOut.write('&#8594; Homopolymers instance : <b>%s</b><br/>\n'%self.HomopolymersPath)
			if not self.HomopolymersPath:
				homoCreate = False
		if cheCreate:
			self.check_contaminationPath = self.getPluginPath("Check_contamination")
			htmlOut.write('&#8594; Check_contamination instance : <b>%s</b><br/>\n'%self.check_contaminationPath)
			if not self.check_contaminationPath:
				cheCreate = False
		if finCreate:
			self.finalReportPath = self.getPluginPath("FinalReport")
			htmlOut.write('&#8594; finalReport instance : <b>%s</b><br/>\n'%self.finalReportPath)
			if not self.finalReportPath:
				finCreate = False			
		
		# Get bam filenames.
		with open(os.path.join(self.json_dat['runinfo']['basecaller_dir'], 'datasets_basecaller.json'), 'r') as f:
			json_basecaller = json.load(f)
		bamPaths = []
		bams = []
		for datum in json_basecaller['datasets']:
			#if reference_path != '':
                        tempPath = os.path.join(self.json_dat['runinfo']['alignment_dir'], datum['file_prefix']+'.bam')
			if os.path.exists(tempPath):
                                print 'adding BAM: %s' % tempPath
                                bamPaths.append(tempPath)
                                if datum['dataset_name'][datum['dataset_name'].rfind('/')+1:] != 'No_barcode_match' and '/' in datum['dataset_name']:
                                        bams.append(datum['dataset_name'])
                        else:
				tempPath = os.path.join(self.json_dat['runinfo']['basecaller_dir'], datum['file_prefix']+'.basecaller.bam')
                                if os.path.exists(tempPath) and not re.search('nomatch',tempPath) :
                                        print 'adding BAM: %s' % tempPath
                                        bamPaths.append(tempPath)
                                        if datum['dataset_name'][datum['dataset_name'].rfind('/')+1:] != 'No_barcode_match' and '/' in datum['dataset_name']:
                                                bams.append(datum['dataset_name'])

		# get the list of 'valid' barcodes or samples (could be either depending on whether user altered names with run planning
		# and sort of hacky, but extract this from the BAM file names we just got above
		for bamFileName in bamPaths:
                        print 'BAM PATH : ' + bamFileName
			barcodeName = bamFileName.split('/')[-1] # get the last part, just the name with no path (probably can use os method here too)
			barcodeName = barcodeName.split('_rawlib')[0] # get just the barcode part of the name
			# find a possible matching sample name
                        print 'barcode = ' + barcodeName
			for sampleItemName in samples:
				sampleItem = samples[sampleItemName]
                                #print 'Item = ' + str(sampleItem)
				if barcodeName in sampleItem['barcodes']:
					self.sampleNameLookup[barcodeName] = sampleItemName
                                        print barcodeName + ' = ' + sampleItemName 
			if barcodeName in self.sampleNameLookup.keys():
				sampleName = self.sampleNameLookup[barcodeName]
                                print sampleName +' => '+ barcodeName
			else:
				sampleName = 'None'
				self.sampleNameLookup[barcodeName] = 'None' # makes it much easier later to do the lookup with no conditional tests
                                target = d['objects'][0]['store']['barcodes'][barcode]['targets_bed']
                                if target == '/results/uploads/BED/5/hg19/unmerged/detail/ColonLungV2.20140523.designed.with_NM.bed':
                                        self.projects[barcodeName] = 'ADN_circulant_colon_lung'
                                elif target == '/results/uploads/BED/2/hg19/unmerged/detail/IAD118795_231_Designed_MET_000245.with_NM.bed':
                                        self.projects[barcodeName] = 'colon_lung_routine'
                                elif target == '/results/uploads/BED/7/hg19/unmerged/detail/Oncomine_Focus.20150316.designed.plusPurpose.bed':
                                        self.projects[barcodeName] = 'essais_precoces'
                                elif target == '/results/uploads/BED/6/hg19/unmerged/detail/IAD112641_173_designed.with_NM.bed':
                                        self.projects[barcodeName] = 'Safir02'

			# MGD note: I considered setting blank sample names to the barcode name instead, but might not be what customer intended
			print 'BARCODE FOUND: %s SAMPLE ID: %s POJECT: %s' % (barcodeName, sampleName, self.projects[barcodeName])
			self.barcodeNames.append(barcodeName)
		self.sampleNameLookup[''] = '' # allows us to easily handle case where barcode might not have been found

		# log basic info for debug purposes
		print 'PLUGINCONFIG:'
		print '----------------------------------------------'
		print 'DELIMETER: "%s"'%delim
		htmlOut.write('<b>DELIMITER:</b> "%s"<br/>\n<b>SELECTIONS:</b><br/>\n'%delim)
		print 'SELECTIONS:'
		for sel in selections:
			print '  %s'%sel
			htmlOut.write('\t%s<br/>\n'%sel)
		print '----------------------------------------------'
		
		# Produce string to rename to.
		self.renameString = ""
		# Avoid delimiting anywhere but in between the arguments.
		firstSelectionDone = False
		for sel in selections:
			if sel != '':
				if firstSelectionDone:
					self.renameString += delim
				self.renameString += sel
				firstSelectionDone = True
		print 'BASE RENAME STRING: %s' % self.renameString

		# Create patients folders
		print 'CREATING PATIENTS FOLDERS...'
		self.createOneFolderByPatient(bamPaths)
		
		# Copy .bam and .bai file(s) if requested.
		if (bamCreate):
			self.copyBam(bamPaths)

		# Create fastq file(s) if requested.
		if (fastqCreate):
			# create FASTQ file(s)
			fromDir = '%s/FastqCreator.py'%self.envDict['DIRNAME']
			toDir = self.envDict['TSP_FILEPATH_PLUGIN_DIR']
			print 'cp: from %s\n to %s\n'%(fromDir, toDir)
			fqCmd = Popen(['cp', fromDir, toDir], stdout=PIPE, env=self.envDict)
			fqOut, fqErr = fqCmd.communicate()
			print 'exec: %s/FastqCreator.py'%toDir
			FastqCmd = Popen(['python', 'FastqCreator.py'], stdout=PIPE, env=self.envDict)
			FastqOut, FastqErr = FastqCmd.communicate()
			print 'mv: fastq -> specified format.'

			self.copyFastqSff('fastq')

		# Create sff file(s) if requested.
		if (sffCreate):
			fromDir = '%s/SFFCreator.py'%self.envDict['DIRNAME']
			toDir = self.envDict['TSP_FILEPATH_PLUGIN_DIR']
			print 'cp: from %s\n to %s\n'%(fromDir, toDir)
			sfCmd = Popen(['cp', fromDir, toDir], stdout=PIPE, env=self.envDict)
			sfOut, sfErr = sfCmd.communicate()
			print 'exec: %s/SFFCreator.py'%toDir
			SFFCmd = Popen(['python', 'SFFCreator.py'], stdout=PIPE, env=self.envDict)
			SFFOut, SFFErr = SFFCmd.communicate()
			print 'mv: sff -> specified format.'
			# Copy fastq and sff file(s).
			self.copyFastqSff('sff')

		# Copy Vcf and Variant Caller tabular file(s) (xls) if requested.
		if (vcCreate):
			self.copyVcf(bamPaths)
			
		# Copy Variant Annotation file(s) if requested.
		if (tsvCreate):
			print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! COPY VARIANT ANNOTATION TRUE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
			self.copyTsv()

		# Copy Coverage Analysis file(s) if requested.
		if (covCreate):
			self.copyCov()
			
		# Copy Compare_VCF file(s) if requested.
		#if (comCreate):
		#	self.copyCom()
		
		# Copy plot coverage file(s) if requested.
		if (plotCreate):
			self.copyPlot()
			
		# Copy homopolymers file(s) if requested.
		if (homoCreate):
			self.copyHomopoly()
			
		# Copy check_contamination file(s) if requested.
		if (cheCreate):
			self.copyCheck()
			
		# Copy finalReport file(s) if requested.
		if (finCreate):
			self.copyFin()
			
		# Copy pdf report is requested
		if (pdfCreate):
			report = '%s/report.pdf' %self.envDict['REPORT_ROOT_DIR']
			if not os.path.exists(report):
				pdf_url = self.json_dat['runinfo']['net_location'] + "/rundb/getPDF/%s.pdf" % str(self.json_dat['runinfo']['pk'])
				f = urllib2.urlopen(pdf_url)
                        for project in self.projects.values():
                                print 'COPYING PDF FILES: %s --> %s'%(report, project+'/'+self.runName + 'report.pdf')
                                cmd = Popen(['%s/backup_pgm.sh' % (self.envDict['DIRNAME']), 'cp', report, project +'/'+self.runName + '/report.pdf'], stdout=PIPE)
                                out, err = cmd.communicate()
                                print 'OUT: %s\nERR: %s'%(out, err)
                        plugins = '%s/plugins.pdf' %self.envDict['REPORT_ROOT_DIR']
                        if not os.path.exists(plugins):
                                pdf_url = self.json_dat['runinfo']['net_location'] + "/rundb/getPDF/%s.pdf" % str(self.json_dat['runinfo']['pk'])
                                f = urllib2.urlopen(pdf_url)
                        for project in self.projects.values():
                                print 'COPYING PDF FILES: %s --> %s'%(plugins, project+'/'+self.runName + 'plugins.pdf')
                                cmd = Popen(['%s/backup_pgm.sh' % (self.envDict['DIRNAME']), 'cp', plugins, project +'/'+self.runName + '/plugins.pdf'], stdout=PIPE)
                                out, err = cmd.communicate()
                                print 'OUT: %s\nERR: %s'%(out, err)

			
		# Copy log files if requested
		if (logCreate):
			logs = '%s/pgm_logs.zip' %self.envDict['REPORT_ROOT_DIR']
                        for project in self.projects.values():
                                print 'COPYING LOG FILES: %s --> %s'%(logs, project+'/'+self.runName + 'pgm_logs.zip')
                                cmd = Popen(['%s/backup_pgm.sh' % (self.envDict['DIRNAME']), 'cp', logs, project +'/'+self.runName + '/pgm_logs.zip'], stdout=PIPE)
                                out, err = cmd.communicate()
                                print 'OUT: %s\nERR: %s'%(out, err)
		
		htmlOut.close()	
		#check temoins positifs
		run_dir=str(os.environ['ANALYSIS_DIR'])
		run_pk=str(self.json_dat['runinfo']['pk'])
		cmd_check = "python /results/plugins/BackupCopy/scripts/Tru-Q_NGS_DNA3_check.py " + run_dir + " " + run_pk
		print cmd_check
		os.system(cmd_check)
		
                cmd_check = "python /results/plugins/BackupCopy/scripts/HD730_check.py " + run_dir + " " + run_pk
                print cmd_check
                os.system(cmd_check)
                cmd_check = "python /results/plugins/BackupCopy/scripts/HD779_check.py " + run_dir + " " + run_pk
                print cmd_check
		os.system(cmd_check)
                cmd_check = "python /results/plugins/BackupCopy/scripts/HD798_check.py " + run_dir + " " + run_pk
                print cmd_check
                os.system(cmd_check)
                cmd_check = "python /results/plugins/BackupCopy/scripts/TEMOIN_BRCA_SOMATIC_check.py " + run_dir + " " + run_pk
                print cmd_check
                os.system(cmd_check)
                cmd_check = "python /results/plugins/BackupCopy/scripts/Seraseq_check.py " + run_dir + " " + run_pk
                print cmd_check
                os.system(cmd_check)

		return True

if __name__ == "__main__":
	PluginCLI(BackupCopy())
