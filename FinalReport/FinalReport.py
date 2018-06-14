 #!/usr/bin/env python
import os
import json
import csv
import urllib2
from ion.plugin import *
from subprocess import *
import xlwt
import xlrd

class FinalReport(IonPlugin):
	"""Plugin object to check contamination on control sample"""
	version = "0.1"
	json_dat = {}
	minX = 300
        projects = {}
        project = 'unknown'
        sample_type = {}
	#run_CLv2 = False
	#run_CHP2 = False
	#run_BRCA = False
	#run_ADNtc = False
        #run_safir = False
	
	def getPluginPath(self, pluginName):
		pluginPath = False
		try:
			print self.json_dat['runinfo']['api_url']
			api_url = self.json_dat['runinfo']['api_url'] + '/v1/pluginresult/?format=json&plugin__name=%s&result=%s' % (pluginName,str(self.json_dat['runinfo']['pk']))
			api_key = self.json_dat['runinfo'].get('api_key', None)
			if api_key is not None:
				api_url = api_url + '&api_key=%s' % api_key
				print 'Using API key: %s' % api_key
			else:
				print 'No API key available'
			f = urllib2.urlopen(api_url)
			d = json.loads(f.read())
		
			for plugin in d['objects']: # they are already sorted from newest to oldest, so take the first 'completed' one
				if plugin['state'] == 'Completed':
					pluginPath = plugin['path']
					break
			if not pluginPath:
				print 'WARNING! No completed instance found for plugin "%s"' % pluginName
		except:
			print 'ERROR!  Failed to get variant caller path'
	
		print 'INFO: using %s path: %s' % (pluginName,pluginPath)
		return pluginPath
		
	#def getProject(self):		
#		try:
#			api_url = self.json_dat['runinfo']['api_url'] + '/v1/results/%s/?format=json' % str(self.json_dat['runinfo']['pk'])
#			print 'API URL: ' + api_url
 #                       f = urllib2.urlopen(api_url)
#			d = json.loads(f.read())
#			for project in d['projects']:
#				num = int(project.split('/')[-2])
#				print "Project number found :"
#				print num
#				if num == 4: #project = essais_precoces
#					self.run_CLv2 = False
#					self.run_CHP2 = True
#					self.run_BRCA = False
#					self.run_ADNtc = False
 #                                       self.run_safir = False
#				if num == 2: #project = colun_lung_plus_routine
#					self.run_CLv2 = True
#					self.run_CHP2 = False
#					self.run_BRCA = False
#					self.run_ADNtc = False
 #                                       self.run_safir = False
#				if num == 6: #project = BRCA_qiagen
#					self.run_CLv2 = False
#					self.run_CHP2 = False
#					self.run_BRCA = True
#					self.run_ADNtc = False
 #                                       self.run_safir = False
#				if num == 3: #project = ADN_circulant_colonlung
#					self.run_CLv2 = False
#					self.run_CHP2 = False
##					self.run_BRCA = False
#					self.run_ADNtc = True
 #                                       self.run_safir = False
  #                              if num == 5: #project = Safir2                             
   #                                    self.run_CLv2 = False
    #                                   self.run_CHP2 = False
     #                                  self.run_BRCA = False
      #                                 self.run_ADNtc = False
       #                                self.run_safir = True
	#	except:
	#		print 'ERROR!  Failed to get project number'
    
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
				
		sample2barcode = {}
		if isinstance(self.json_dat['plan']['barcodedSamples'],dict):
			samples = self.json_dat['plan']['barcodedSamples']
		else:
			samples = json.loads(self.json_dat['plan']['barcodedSamples'])
		for key in samples.keys():
			#sample2barcode[key] = samples[key]["barcodes"][0]
                        for bc in samples[key]["barcodes"]:
                                bc_info=samples.get(key).get('barcodeSampleInfo').get(bc)
                                nucType=bc_info['nucleotideType']
                                
                                if nucType == 'DNA' :
                                        barcode=bc
                                        #sample_name[bc]=key
                                        sample2barcode[key] = bc
                                        self.sample_type[barcode] = 'DNA'
                                        target = bc_info['targetRegionBedFile']
                                        #target_path[barcode] = target
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
                                        barcode=bc
                                        self.sample_type[barcode] = 'RNA'
                                        target = ''
                                        self.projects[barcode] = 'essais_precoces'

                                project = self.projects[barcode]
                                print 'BARCODE : ' + barcode + ' TARGET : ' + target + ' PROJECT : ' + project


		variantAnnotationPath = self.getPluginPath('variantAnnotation')
		coverageAnalysisPath = self.getPluginPath('coverageAnalysis')
		plotCoveragePath = self.getPluginPath('plotCoverage')
		variantCallerPath = self.getPluginPath('variantCaller')
		
		#in case of heterogen run and final report + backup not done for the first target sample group : enter the path of plugin output of interest and execute FinalReport then Backup_Copy
		#variantAnnotationPath = "/results/analysis/output/Home/Auto_user_PGM-404-Ion_Chef_Colon_Lung_AmpliSeq_v2_318_484_516/plugin_out/variantAnnotation_out.2980"
		#coverageAnalysisPath = self.getPluginPath('coverageAnalysis')
		#plotCoveragePath = "/results/analysis/output/Home/Auto_user_PGM-404-Ion_Chef_Colon_Lung_AmpliSeq_v2_318_484_516/plugin_out/plotCoverage_out.2983"
		#variantCallerPath = "/results/analysis/output/Home/Auto_user_PGM-404-Ion_Chef_Colon_Lung_AmpliSeq_v2_318_484_516/plugin_out/variantCaller_out.2979"
		
		# WORKBOOK STYLES
		
		# headerStyle_string = "font:name Calibri, height 220, bold on;borders"
		headerStyle = xlwt.easyxf('font: name Calibri, bold on, height 220; borders: left thin, top thin, bottom thin, right thin')
		generalStyle = xlwt.easyxf('font: name Calibri, height 220')
		greenBGStyle = xlwt.easyxf('font: name Calibri, height 220; pattern: pattern solid, fore-colour 0x2A')
		redBGStyle = xlwt.easyxf('font: name Calibri, height 220; pattern: pattern solid, fore-colour 0x1D')
		
		#self.getProject()
		
		#if self.run_CLv2:
	#		print "\n - Run type = Colon Lung Plus Routine\n"
			#hotspots_xls = xlrd.open_workbook("%s/database_hotspot_sensitivity_MAJ_decembre2015.xls" % self.json_dat['runinfo']['plugin']['path'])
	#	elif self.run_CHP2:
	#		print "\n - Run type = Essais Precoces\n"
			#hotspots_xls = xlrd.open_workbook("%s/hotspots_details_HEMATO.xls" % self.json_dat['runinfo']['plugin']['path'])
	#	elif self.run_BRCA:
	#		print "\n- Run type = BRCA qiagen\n"
	#	elif self.run_ADNtc:
	#		print "\n- Run type = ADNtc Colon Lung\n"
			#hotspots_xls = xlrd.open_workbook("%s/database_hotspot_sensitivity_MAJ_decembre2015.xls" % self.json_dat['runinfo']['plugin']['path'])
         #       elif self.run_safir:
          #              print "\n- Run type = Safir2\n"
 

		#h_sheet = hotspots_xls.sheet_by_index(0)

		target_path={}
                #api_url = self.json_dat['runinfo']['api_url'] + '/v1/pluginresult/?format=json&plugin__name=variantCaller&result=' + str(self.json_dat['runinfo']['pk'])
                api_url = os.getenv('RUNINFO__API_URL',
                           'http://localhost/rundb/api') + '/v1/pluginresult/?format=json&plugin__name=variantCaller&result=' + str(self.json_dat['runinfo']['pk'])
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
                                variantCaller_path = plugin['path']
                                #print 'variantCaller path:' + variantCaller_path
                                break
                if not variantCaller_path:
                        print 'WARNING! No completed instance found for plugin variantCaller'

                #get target for specific barcode
                for barcode in d['objects'][0]['store']['barcodes']:
                    target = d['objects'][0]['store']['barcodes'][barcode]['targets_bed']
                    #print 'BARCODE : ' + barcode + ' TARGET : ' + target
                    target_path[barcode] = target


		for sample in sample2barcode.keys():
			barcode = sample2barcode[sample]
                        project = self.projects[barcode]
                        #if sample=='VIRELEGOUX_LN398':
                        #        barcode='IonXpress_086'
                        #elif sample=='MAHAGNE_LK891':
                        #        barcode='IonSelect-15'
                        #elif sample=='none':
                        #        barcode='IonSelect-10'
                        #print 'sample: ' + sample + 'index: ' + barcode
                        if project=='Safir02':
                                va_filePath = variantAnnotationPath + "/%s/CHR/%s_NGS_VC.tsv" % (barcode,sample)
                        else:
                                va_filePath = variantAnnotationPath + "/%s/CHR/%s_NGS_Diag_VC_TS.tsv" % (barcode,sample)
			covFilePath = coverageAnalysisPath + "/%s/%s_%s_%s.amplicon.cov.xls" % (barcode,barcode,self.json_dat['expmeta']['run_name'],self.json_dat['expmeta']['results_name'])
                        #print 'covFilePath: '+ covFilePath
			baseCovFilePath = plotCoveragePath + "/%s_%s_cov_less_than_300x.bed" % (sample,barcode)
			vc_filePath = variantCallerPath + "/%s/alleles.xls" % barcode

			finalReport = xlwt.Workbook(encoding='utf-8')
			annotationSheet = finalReport.add_sheet("Annotation", cell_overwrite_ok=True)
			coverageSheet = finalReport.add_sheet("Amplicon Coverage")
			baseCoverageSheet = finalReport.add_sheet("Base Coverage -300X")
			variantSheet = finalReport.add_sheet("VariantCaller")

			########################
			## FEUILLE ANNOTATION ##
			########################
			
			#hotspot_cp = []
			#hotspot_cp_sensi = []
			#if self.run_LAM:
			#	for row in range(h_sheet.nrows-1):
			#		hotspot_cp.append((h_sheet.cell(row+1,5).value,h_sheet.cell(row+1,7).value,h_sheet.cell(row+1,8).value))
			#if self.run_SBT:
			#	for row in range(h_sheet.nrows-1):
			#		hotspot_cp_sensi.append((h_sheet.cell(row+1,5).value,h_sheet.cell(row+1,7).value,h_sheet.cell(row+1,8).value,h_sheet.cell(row+1,10).value))
					
			try:	# note : making a try is for H2O, when 0 mutations found = no annotation file produced = bug
				if os.path.isfile(va_filePath) and os.access(va_filePath, os.R_OK):
					va_file_ok = 1
					if va_file_ok == 1 :
						if os.stat(va_filePath).st_size>0:
							va_file = open(va_filePath,"r")
							va_file_reader = csv.reader(va_file,delimiter="\t")
					
							header = va_file_reader.next()
				#if self.run_SBT:
				#	header.insert(14,"Sensitivity")
							for i in range(len(header)):
								annotationSheet.row(0).write(i,header[i],headerStyle)
							l=1
							for va_line in va_file_reader:
								if va_line : # for empty lines (exist empty lines in diag format for PAM after sorting variants)
									row = annotationSheet.row(l)
							# run SBT : rajout donnees sensibilite si presente
							# if self.run_SBT:
							#	va_gene_cp = (va_line[3],va_line[7],va_line[8])
							#	sensi = "not found"
							#	for h in hotspot_cp_sensi:
							#		if va_gene_cp == (h[0],h[1],h[2]):
							#			sensi = h[3]
							#	va_line.insert(14,sensi)
									if (va_line[12] != "intronic" and va_line[12] != "UTR3" and va_line[12] != "UTR5") and (va_line[13] != "synonymous"):					
										for i in range(len(va_line)):
											cell_val = self.representsInt(va_line[i])
											row.write(i,cell_val,greenBGStyle)
									elif (va_line[12] == "intronic" and va_line[20] == "NA"): #intronic without MAF
										for i in range(len(va_line)):
											cell_val = self.representsInt(va_line[i])
											row.write(i,cell_val,greenBGStyle)
									elif (va_line[17] != "NA"): #WITH COSMIC ID
										for i in range(len(va_line)):
											cell_val = self.representsInt(va_line[i])
											row.write(i,cell_val,greenBGStyle)
									else:
										for i in range(len(va_line)):
											cell_val = self.representsInt(va_line[i])
											row.write(i,cell_val,generalStyle)
					
							# modif Lisa run hemato : hotspot en rouge
							# if self.run_LAM:
							#	va_gene_cp = (va_line[3],va_line[5],va_line[6])
							#	for h in hotspot_cp:
							#		if va_gene_cp == h:
							#			row.write(5,va_line[5],redBGStyle)
							#			row.write(6,va_line[6],redBGStyle)
								l=l+1

					# ajout "Amplicons <300X: "
					#row = annotationSheet.row(l+2)
					#row.write(0,"Amplicons < 300X: ", generalStyle)
				
							va_file.close()
					else:
						print "WARNING: annotation file found for %s but reading problem." % sample



				else:
					print "WARNING : Annotation file not found for %s. Annotation sheet is empty." % sample
					va_file_ok = 0
			except:
				print "WARNING : Annotation file not found for %s. Annotation sheet is empty." % sample
				
			
		
			######################
			## FEUILLE COVERAGE ##
			######################
			try:
				covFile = open(covFilePath,"r")
				covFile_reader = csv.reader(covFile,delimiter="\t")

				header = covFile_reader.next()
				header[4] = "gene_id"
				for i in range(len(header)):
					coverageSheet.row(0).write(i,header[i],headerStyle)
				l=1
				red_ampl = []
				for cov_line in covFile_reader:
					cov_line[4] = cov_line[4].split(";")[0].split("GENE_ID=")[-1] #colonne attribute, garder seulement gene id
					row = coverageSheet.row(l)
					if int(cov_line[9]) < self.minX:
						red_ampl.append(cov_line[3])
						for i in range(len(cov_line)):
							cell_val = self.representsInt(cov_line[i])
							row.write(i,cell_val,redBGStyle)
					else:
						for i in range(len(cov_line)):
							cell_val = self.representsInt(cov_line[i])
							row.write(i,cell_val,generalStyle)
					l=l+1
				covFile.close()
			except:
				print "WARNING : Coverage file not found for %s. Coverage sheet is empty." % sample

			###########################
			## FEUILLE BASE COVERAGE ##
			###########################

			try:	#note : this is for H2O, cov_less_than_300x do not exist

				bcovFile = open(baseCovFilePath,"r") 
				bcovFile_reader = csv.reader(bcovFile,delimiter="\t")
				bcov_file_ok = 1

				if bcov_file_ok == 1 :
				# creation dic amplicon par chr
					design = {}
                                        #target for barcode
					targetbed = target_path[barcode]
                                        targetbedFile = open(targetbed,"r")
					targetbed_reader = csv.reader(targetbedFile,delimiter="\t")
					targetbed_reader.next()
					for targetbed_line in targetbed_reader:
						chrom = targetbed_line[0]
						start = int(targetbed_line[1])
						stop = int(targetbed_line[2])
						ampl = targetbed_line[3]
						gene_id = targetbed_line[7].split(";")[0].split("GENE_ID=")[-1]
						if not chrom in design:
							design[chrom] = []
						design[chrom].append([ampl,start,stop,gene_id])

					to_write = []

					for bcov_line in bcovFile_reader:
					#print bcov_line
						if bcov_line[0] == '':	# test ligne vide = fichier vide, pas de base >300x
							continue
						chrom = bcov_line[0]
						start = int(bcov_line[1])
						stop = int(bcov_line[2])

					# recherche du/des amplicons de la region -300x
						first_ampl = "NF"
						last_ampl = "NF"
						for ampl in design[chrom]:
							ampl_name =  ampl[0]
							ampl_start = ampl[1]
							ampl_stop = ampl[2]
							if (first_ampl == "NF") and (start >= ampl_start and start <= ampl_stop):
								first_ampl = ampl_name
							if (last_ampl == "NF") and (stop >= ampl_start and stop <= ampl_stop):
								last_ampl = ampl_name
								break

					#recherche des hotspots de la region -300x
#					hots = []
#					for row in range(h_sheet.nrows-1):
						# tester si start ou stop du hotspot se trouve dans region
#						if (h_sheet.cell(row+1,1).value >= start and h_sheet.cell(row+1,1).value <= stop) or (h_sheet.cell(row+1,2).value >= start and h_sheet.cell(row+1,2).value <= stop):
							#hots.append([h_sheet.cell(row+1,3).value,h_sheet.cell(row+1,1).value,h_sheet.cell(row+1,2).value,h_sheet.cell(row+1,9).value])
#							hots.append(h_sheet.cell(row+1,8).value) # uniquement p. pour l'instant

						if first_ampl == "NF" or last_ampl == "NF":
							print "error : amplicon not found for : " + chrom + " : " + str(start) + " - " + str(stop)
							to_write.append([bcov_line[0],start,stop,"not found"])
						else:
							started = False
							list_ampl = []
							list_region = []
							line = [chrom,start,stop,"",""]
							for ampl in design[chrom]:
								if first_ampl == ampl[0]:
									started = True
								if last_ampl == ampl[0]:
									list_ampl.append(ampl[0])
									line[3] = list_ampl[0]
									if len(list_ampl)>1:
										for i in range(len(list_ampl)-1):
											line[3] = line[3] + ", " + list_ampl[i+1]
									if ampl[3] not in list_region:
										list_region.append(ampl[3])
									line[4] = list_region[0]
									if len(list_region)>1:
										for i in range(len(list_region)-1):
											line[4] = line[4] + ", " + list_region[i+1]
								
								# lines below commented 
								# h=""
								#if hots:
								#	if len(hots) > 1:
								#		h = hots[0]
								#		for i in range(len(hots)-1):
								#			h = h + ", " + hots[i+1]
								#	else:
								#		h = hots[0]
								#line.append(h)
									to_write.append(line)
									break
								elif started:
									list_ampl.append(ampl[0])
									if ampl[3] not in list_region:
										list_region.append(ampl[3])

					red_line = []
					white_line = []
					for line in to_write:
						red = False
						for red_a in red_ampl:
							if red_a in line[3].split(","):
								red = True
						if red:
							red_line.append(line)
						else:
							white_line.append(line)

					headerline = baseCoverageSheet.row(0)
				#header = ["Chr","Start","End","Amplicon","Region","Hotspots"]
					header = ["Chr","Start","End","Amplicon","Region"]
					for i in range(len(header)):
						headerline.write(i,header[i],headerStyle)
					l=1
					for line in red_line:
						row = baseCoverageSheet.row(l)
						for i in range(len(line)):
							row.write(i,line[i],redBGStyle)
						l=l+1
					for line in white_line:
						row = baseCoverageSheet.row(l)
						for i in range(len(line)):
							row.write(i,line[i],generalStyle)
						l=l+1

				bcovFile.close()

			except:
				print "WARNING : Base coverage file not found for %s. Base Coverage sheet is skiped." % sample

			
			
			

			###########################
			## FEUILLE VARIANTCALLER ##
			###########################
			
			try:	#note : if unsufficient reads, variantcaller don't produce results
				vc_file = open(vc_filePath,"r")
				vc_file_reader = csv.reader(vc_file,delimiter="\t")
			
				header = vc_file_reader.next()
				for i in range(len(header)):
					variantSheet.row(0).write(i,header[i],headerStyle)
				l=1
				for vc_line in vc_file_reader:
					row = variantSheet.row(l)
					for i in range(len(vc_line)):
						cell_val = self.representsInt(vc_line[i])
						row.write(i,cell_val,generalStyle)
					l=l+1
				vc_file.close()
			except:
				print "WARNING : alleles.xls file not found for %s. variantCaller sheet is skiped." % sample
			
			######
			
			cmd = Popen(['mkdir', self.json_dat['runinfo']['results_dir']+"/"+barcode], stdout=PIPE)
			out, err = cmd.communicate()
			finalReport.save("%s/%s/%s.finalReport.xls" % (self.json_dat['runinfo']['results_dir'],barcode,sample))

		return True

if __name__ == "__main__":
    PluginCLI(FinalReport())
