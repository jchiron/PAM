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

class BackupIR(IonPlugin):
	""" Plugin object which copy all desired IonReporter results file from the run into a folder."""
	""" A hierarchical filesystem is created (project->run->patients) to store the data, and files name are uniformly renamed"""
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
        folder = '/media/PGM/results/IonReporter'
	IR_samples=[]
	doneSamples=[]
	
	# paths to the most recent, completed instance of the IonReporter uploader plugin
	#IRuploaderPath = False
        	
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
			try:
				patientName = self.sampleNameLookup[barcodeName]
				
				if not re.search("_RNA_",patientName,re.I): # 1 result per patient (DNA +RNA) under same folder named as DNA sample
					for s_ir in self.IR_samples:
						if patientName in s_ir:
							sampleName = s_ir
							print 'CREATING FOLDER : %s -'%(s_ir)
					# cmd = 'mkdir '+ self.folder + '/' + self.runPath + '/' + patientName
					# os.system(cmd)
					# print 'cmd %s' % cmd

							cmd = Popen(['mkdir', self.folder+'/'+self.runPath+'/'+s_ir], stdout=PIPE)
					# cmd = Popen(['%s/backup_pgm.sh' % (self.envDict['DIRNAME']), 'mkdir', self.runPath + '/' + patientName], stdout=PIPE)
							out, err = cmd.communicate()
					#print 'OUT: %s\nERR: %s'%(out, err)

					# get unfiltered zip, current tsv files paths and dowload them
							#self.getSampleIRdata(patientName)
							self.getSampleIRdata(s_ir)
			except:
				'No name found for barcode'+ str(barcodeName)
                        
			
        def getSampleIRdata(self, sampleName):
            #get json file for sample
            jsonFile = '%s/%s/%s/%s.json' % (self.folder,self.runPath, sampleName, sampleName)
            #print 'Destination json file : %s' % jsonFile

	    #troubleshoot
            #cmd = 'curl -k -H "Authorization:YzdkYzM5MmFlZTg4NDJhMzBkODQ3ZWIwMzgyZmJkMmI0NzVhZDNmZjQzMTkzYjY1NDYyMDg5MzAwMGRiNGRkZg" "https://129.21.10.4/webservices_42/rest/api/analysis?format=json&name=%s &type=sample" > %s '%(sampleName,jsonFile)
	    #print 'cmd %s' % cmd

            print 'Get IR json for sample %s ' % sampleName 
	    #print self.IR_samples
	    #for s_ir in self.IR_samples:
		#    print s_ir
		 #   if sampleName in s_ir:
		#	    print 'sample %s = ir_id %s'% (sample_name,s_ir)
		#	    sampleName = s_ir
		    #else :
			#  print 'sample %s not in ir_id %s'% (sample_name,s_ir)  
            output = open(jsonFile,"w")
	    #!!!!!!!!!! TODO : changer car pas toujours v1, ex ROUBY_KZ492_v4 !!!!!!!!!
            url = "https://129.21.10.4/webservices_42/rest/api/analysis?format=json&name="+sampleName+"&type=sample"
	    #print 'URL %s'%url
            cmd = Popen(['curl', '-k', '-H', "Authorization:NThiZGY1MzE0MTA5NWM4MTMzOGIyOWNkOGQzYjg1MjY0MjBmOTg1NjE5ODNmNjM5YjQzOTQxMjM5NzI3MzQ3Mg", url],stdout=output)
            out, err = cmd.communicate()
            print 'OUT: %s\nERR: %s'%(out, err)
            output.close()

            #get unfiltered variants  zip file for sample
            with open(jsonFile, 'r') as f:
		    #json_list = json.load(f)
		    json_IR = json.load(f,parse_float=str)
	    f.close()
	    #print type(json_IR)
	    	    
	    status = False
	    analysis_name = False
	    #for data in json_IR:
		    #if isinstance(data['status'],dict):
		#	    status = data['status']
		 #   else:
	    #status = str(json_IR.get('status'))
	
	    status = json_IR[0]['status']
	    #print 'status :' + status
	    if status == 'SUCCESSFUL':
		    print 'Get IR unfiltered.zip for sample %s ' % sampleName
		    # get zip url
		    zip_url = json_IR[0]['data_links']['unfiltered_variants']
		    zip_url.replace('\\','')
		    print 'unfiltered_variants url :'+zip_url

		    # upload IR unfiltered.zip for sample to torrent server
		    uploaded_path = '%s/%s/%s/%s_filtered.zip' % (self.folder,self.runPath, sampleName, sampleName)
		    #to troubleshoot
		    cmd = 'curl -k -H "Authorization:YzdkYzM5MmFlZTg4NDJhMzBkODQ3ZWIwMzgyZmJkMmI0NzVhZDNmZjQzMTkzYjY1NDYyMDg5MzAwMGRiNGRkZg" "%s" -o %s '%(zip_url,uploaded_path)
		    print 'cmd %s' % cmd
                    cmd = Popen(['curl', '-k', '-H', "Authorization:YzdkYzM5MmFlZTg4NDJhMzBkODQ3ZWIwMzgyZmJkMmI0NzVhZDNmZjQzMTkzYjY1NDYyMDg5MzAwMGRiNGRkZg", '-H', "Content-Type:application/x-www-form-urlencoded", zip_url, '-o', uploaded_path])
                    
                    
		    #cmd = Popen(['curl', '-k', '-H', "Authorization:NThiZGY1MzE0MTA5NWM4MTMzOGIyOWNkOGQzYjg1MjY0MjBmOTg1NjE5ODNmNjM5YjQzOTQxMjM5NzI3MzQ3Mg", zip_url, '-o', uploaded_path])
		    out, err = cmd.communicate()
		    #print 'OUT: %s\nERR: %s'%(out, err)
		    
		    #extract zip filtered variants to get IR vcf file
		    tmp_dir = '%s/%s/%s/temp' % (self.folder,self.runPath,sampleName)
		    os.mkdir(tmp_dir)
		    cmd = 'unzip -d %s %s' % (tmp_dir,uploaded_path)
		    os.system(cmd)

		    for f in os.listdir(tmp_dir) :
			    #print 'FILE %s\n' % f
			    if f.endswith('.zip'):
				    #print 'ZIP FILE %s\n' % f
				    cmd = 'unzip -d %s/%s/%s %s/%s' % (self.folder,self.runPath,sampleName,tmp_dir,f)
				    #print 'UNZIP CMD %s\n' % cmd
				    os.system(cmd)

		    cmd = 'rm -r %s' % tmp_dir
		    os.system(cmd)

		    variantsDir = '%s/%s/%s/Variants/' % (self.folder,self.runPath,sampleName)
		    for dirPath,dirNames,f in os.walk(variantsDir) :
			    #print 'FILENAME %s\n' % f
			    for filename in f :
				    if filename.endswith('.vcf'):
                                            #print 'VCF FILE %s\n' % filename
                                            vcfFile = os.path.join(dirPath,filename)
                                            #print 'vcfFile PATH %s\n' % vcfFile
					    cmd = "sed -i \'s/\"//g\' %s" % vcfFile
					    os.system(cmd)
				    elif filename.endswith('full.tsv'):
					    #print 'FULL FILE %s\n' % filename
					    fullTsv = os.path.join(dirPath,filename)
					    #print 'FULL PATH %s\n' % fullTsv
					    cmd = "sed -i \'s/\"//g\' %s" % fullTsv
                                            os.system(cmd)
				    elif filename.endswith('oncomine.tsv'):
                                            #print 'ONCOMINE FILE %s\n' % filename
                                            oncomineTsv = os.path.join(dirPath,filename)
                                            #print 'ONCOMINE PATH %s\n' % oncomineTsv
					    cmd = "sed -i \'s/\"//g\' %s" % oncomineTsv
                                            #print cmd
					    os.system(cmd)

		    ######
		    # OKR report creation
		    ######

		    #????? TODO change to backupCopy floder path
		    print 'Get OKR report.pdf for sample %s ' % sampleName
		    report_path = '%s/%s/%s/%s_OKR_report.pdf' % (self.folder,self.runPath, sampleName, sampleName)
		    #print 'OKR report ouput path : %s\n' % report_path
		    #cmd = Popen(['curl', '--user', "democonfig:123456", "--request", "POST", "--location", "http://129.21.10.4:8088/api/okr/report/file", "-F","file=",uploaded_apth, "-F","options=\"{\"filterPresetID\":102,\"reportTemplateID\":101,\"reportFormat\":\"PDF\"}\"", '-o', report_path])
		    cmd = "curl --user democonfig:123456 --request POST --location \"http://129.21.10.4:8088/api/okr/report/file\" -F file=@%s -F options=\"{\\\"filterPresetID\\\":102,\\\"reportTemplateID\\\":101,\\\"reportFormat\\\":\\\"PDF\\\",\\\"fields\\\":[{\\\"name\\\":\\\"SampleID\\\",\\\"value\\\":\\\"%s\\\"}]}\" -o %s"% (vcfFile,sampleName,report_path)
		    #print 'OKR CMD: %s\n'%cmd
		    #out,err = cmd.communicate()
		    os.system(cmd)
		    #print 'OKR OUT: %s\nOKR ERR: %s'%(out, err)

		    output_path = '%s/%s/%s/' % (self.folder,self.runPath, sampleName)
                    print 'OUTPUT_PATH: %s -' % output_path
		    
		    tsv_head = '%s%s.head.tsv' % (output_path, sampleName)
		    print 'tsv_head path : '+ tsv_head
		    tmp_ARN = '%s%s.ARN.tsv' % (output_path, sampleName)
		    print 'tmp_ARN path : '+ tmp_ARN
		    tmp_SNV = '%s%s.SNV.tsv' % (output_path, sampleName)
		    print 'tmp_SNV path : '+ tmp_SNV
		    tmp_CNV = '%s%s.CNV.tsv' % (output_path, sampleName)
		    print 'tmp_CNV path : '+ tmp_CNV

		    #header
		    cmd = 'grep \#\#fileDate ' + vcfFile + '>' + tsv_head
		    #print cmd
		    os.system(cmd)
		    cmd = 'grep \#\#FusionSampleQC ' + vcfFile + '>>' + tsv_head
                    os.system(cmd)
		    cmd = 'grep \#\#FusionSampleOverallCall ' + vcfFile + '>>' + tsv_head
                    os.system(cmd)
		    cmd = 'grep \#\#source ' + vcfFile + '>>' + tsv_head
                    os.system(cmd)
		    cmd = 'grep \#\#reference ' + vcfFile + '>>' + tsv_head
                    os.system(cmd)
		    cmd = 'grep \#\#sampleGender ' + vcfFile + '>>' + tsv_head
                    os.system(cmd)
		    cmd = 'grep \#\#sampleDiseaseType ' + vcfFile + '>>' + tsv_head
                    os.system(cmd)
		    cmd = 'grep \#\#AssumedGender ' + vcfFile + '>>' + tsv_head
                    os.system(cmd)
		    cmd = 'grep \#\#CellularityAsAFraction ' + vcfFile + '>>' + tsv_head
		    os.system(cmd)
		    cmd = 'grep \#\#DiploidCopyNumber ' + vcfFile + '>>' + tsv_head
                    os.system(cmd)
		    cmd = 'grep \#\#TotalMappedFusionPanelReads ' + vcfFile + '>>' + tsv_head
                    os.system(cmd)
		    cmd = 'grep \#\#gc_mad ' + vcfFile + '>>' + tsv_head
                    os.system(cmd)
		    cmd = 'grep \#\#ilen_mad ' + vcfFile + '>>' + tsv_head
                    os.system(cmd)
		    cmd = 'grep \#\#mapd ' + vcfFile + '>>' + tsv_head
                    os.system(cmd)
		    cmd = 'grep \#\#pc2_mad ' + vcfFile + '>>' + tsv_head
                    os.system(cmd)
		    cmd = 'grep \#\#pc_mad ' + vcfFile + '>>' + tsv_head
                    os.system(cmd)
		    cmd = 'grep \#\#pcw ' + vcfFile + '>>' + tsv_head
                    os.system(cmd)
		    cmd = 'grep \#\#percent_aligned_reads ' + vcfFile + '>>' + tsv_head
                    os.system(cmd)
		    cmd = 'grep \#\#percent_non_zero_amplicons ' + vcfFile + '>>' + tsv_head
                    os.system(cmd)
		    cmd = 'grep \#\#tmap ' + vcfFile + '>>' + tsv_head
                    os.system(cmd)
		    cmd = 'grep \#\#total_read_count ' + vcfFile + '>>' + tsv_head
                    os.system(cmd)
		    cmd = 'grep \#\#fileUTCtime ' + vcfFile + '>>' + tsv_head
                    os.system(cmd)
		    cmd = 'grep \#\#IonReporterExportVersion ' + vcfFile + '>>' + tsv_head
                    os.system(cmd)
		    cmd = 'grep \#\#OncomineVariantAnnotationToolVersion ' + vcfFile + '>>' + tsv_head
                    os.system(cmd)
		    cmd = 'grep \#\#annotationSources ' + vcfFile + '>>' + tsv_head
                    os.system(cmd)
		    
		    #oncomine.tsv header
		    cmd = "grep rowtype " +oncomineTsv + '>' +tmp_ARN
                    os.system(cmd)
		    #RNA control
		    cmd = 'grep ExprControl ' + oncomineTsv + '>>' + tmp_ARN
		    #print cmd
		    os.system(cmd)
		    
		    #RNA imbalance present
		    cmd = "awk 'BEGIN{FS=\"\t\"};{ if ($2==\"5p3pAssays\" && $3==\"POS\") print $0 }' "+oncomineTsv + '>>' + tmp_ARN
		    #print cmd
		    os.system(cmd)

		    #RNA FUSION present
		    cmd = "awk 'BEGIN{FS=\"\t\"};{ if ($2==\"Fusion\" && $3==\"POS\") print $0 }' "+oncomineTsv + '>>' +tmp_ARN
		    #print cmd
		    os.system(cmd)

		    #RNA imbalance NoCall
		    cmd = "awk 'BEGIN{FS=\"\t\"};{ if ($2==\"5p3pAssays\" && $16==\"NOCALL\") print $0 }' "+oncomineTsv + '>>' +tmp_ARN
		    #print cmd
		    os.system(cmd)

		    #RNA imbalance absent
		    cmd = "awk 'BEGIN{FS=\"\t\"};{ if ($2==\"5p3pAssays\" && $16==\"FAIL\") print $0 }' "+oncomineTsv + '>>' +tmp_ARN
		    #print cmd
		    os.system(cmd)

		    #RNA FUSION absent
		    cmd = "awk 'BEGIN{FS=\"\t\"};{ if ($2==\"Fusion\" && $3==\"NEG\") print $0 }' "+oncomineTsv + '>>' +tmp_ARN
		    #print cmd
		    os.system(cmd)

		    #oncomine.tsv header
		    cmd = "grep rowtype " +oncomineTsv + '>' +tmp_SNV
                    os.system(cmd)
		    #SNV present
		    #category 1 : exonic frameshift, missense, nonsense
		    cmd = "awk 'BEGIN{FS=\"\t\"};{ if (($2==\"snp\" || $2==\"del\" || $2==\"ins\" || $2==\"mnp\" || $2==\"complex\") && $3==\"POS\") print $0 }' "+oncomineTsv + '| grep -v -E \"splicesite|synonymous|intronic|utr\"' + ">>" + tmp_SNV
		    #print cmd
		    os.system(cmd)
		    #category 2 : splicing
		    cmd = "awk 'BEGIN{FS=\"\t\"};{ if (($2==\"snp\" || $2==\"del\" || $2==\"ins\" || $2==\"mnp\" || $2==\"complex\") && $3==\"POS\") print $0 }' "+oncomineTsv + '| grep \"splicesite\"' + ">>" + tmp_SNV
		    #print cmd
		    os.system(cmd)
		    #category 3 : exonic synonymous
		    cmd = "awk 'BEGIN{FS=\"\t\"};{ if (($2==\"snp\" || $2==\"del\" || $2==\"ins\" || $2==\"mnp\" || $2==\"complex\") && $3==\"POS\") print $0 }' "+oncomineTsv + '| grep \"synonymous\"' + ">>" + tmp_SNV
		    #print cmd
		    os.system(cmd)
		    #category 4 : intronic
		    cmd = "awk 'BEGIN{FS=\"\t\"};{ if (($2==\"snp\" || $2==\"del\" || $2==\"ins\" || $2==\"mnp\" || $2==\"complex\") && $3==\"POS\") print $0 }' "+oncomineTsv + '| grep \"intronic\"' + ">>" + tmp_SNV
		    #print cmd
		    os.system(cmd)
		    #category 5 : UTR
		    cmd = "awk 'BEGIN{FS=\"\t\"};{ if (($2==\"snp\" || $2==\"del\" || $2==\"ins\" || $2==\"mnp\" || $2==\"complex\") && $3==\"POS\") print $0 }' "+oncomineTsv + '| grep \"utr\"' + ">>" + tmp_SNV
		    #print cmd
		    os.system(cmd)

		    #SNV NoCall
		    cmd = "awk 'BEGIN{FS=\"\t\"};{ if (($2==\"snp\" || $2==\"del\" || $2==\"ins\" || $2==\"mnp\" || $2==\"complex\") && $3==\"NEG\" && $16==\"NOCALL\") print $0 }' "+oncomineTsv + '>>' +tmp_SNV
		    #print cmd
		    os.system(cmd)

		    #full.tsv header
		    cmd = "grep locus " +fullTsv + '>' +tmp_CNV
		    os.system(cmd)
		    #CNV
		    cmd = "awk 'BEGIN{FS=\"\t\"};{ if ($2==\"CNV\") print $0 }' "+fullTsv + '>>' +tmp_CNV
		    #print cmd
		    os.system(cmd)
		    
		    #tsv_head.close()
		    #tmp_ARN.close()
		    #tmp_SNV.close()
		    #tmp_CNV.close()
		    #write XLS file
		    # WORKBOOK STYLES
		    print '!!!!!!!!!!!!!! CREATING EXCEL FILE !!!!!!!!!!!!!!!!!!'
		    #headerStyle_string = "font:name Calibri, height 220, bold on;borders"
		    #headerStyle = xlwt.easyxf('font: name Calibri, bold on, height 220; borders: left thin, top thin, bottom thin, right thin')
		    generalStyle = xlwt.easyxf('font: name Calibri, height 220')
		    greenBGStyle = xlwt.easyxf('font: name Calibri, height 220; pattern: pattern solid, fore-colour 0x2A')
		    redBGStyle = xlwt.easyxf('font: name Calibri, height 220; pattern: pattern solid, fore-colour 0x1D')

		    finalReport = xlwt.Workbook()
		    arnSheet = finalReport.add_sheet("FUSION")
		    cnvSheet = finalReport.add_sheet("CNV")
		    snvSheet = finalReport.add_sheet("ADN variants")
		    
		    #define header dict fields needed are differents for all sheets
		    #header =''
		    #header_dict=[]
		    #try:	
		    #if os.path.isfile(tsv_head) and os.access(tsv_head, os.R_OK):
		#	    print "ACCESS HEADER file OK"
			    #if os.stat(tsv_head).st_size>0:
	
		    try:
			    head_file = open(tsv_head,'r')
			    print 'OK READ TSV_HEAD FILE'
		    #head_reader = csv.reader(head_file,delimiter="\t")
			    #header = head_reader.next()
			    l = 0
			    l_arn = 0
			    l_cnv = 0
			    l_snv = 0
			    for header_line in head_file:
				    #print header_line
				    fields = header_line.split('\t')
			    #if fields[0]!="Locus":
				    #print 'diff Locus %s'%fields[0]
				    arnSheet.write(l,0,fields[0])
				    cnvSheet.write(l,0,fields[0])
				    snvSheet.write(l,0,fields[0])
				    l = l +1
				    #header = header + header_line
			    #else:
			    #print 'Locus line'
			    #arn header line
			    arnSheet.write(l,0,'Genes',generalStyle)
			    arnSheet.write(l,1,'Locus',generalStyle)
			    arnSheet.write(l,2,'Type',generalStyle)
			    arnSheet.write(l,3,'Variant ID',generalStyle)
			    arnSheet.write(l,4,'Variant Name',generalStyle)
			    arnSheet.write(l,5,'COSMIC/NCBI',generalStyle)
			    arnSheet.write(l,6,'Read Counts',generalStyle)
			    arnSheet.write(l,7,'Detection',generalStyle)
			    arnSheet.write(l,8,"3\'/5\' imbalance",generalStyle)
			    arnSheet.write(l,9,'Read counts Per Million',generalStyle)
			    l_arn = l
			    #cnv header line
			    cnvSheet.write(l,0,'Genes',generalStyle)
			    cnvSheet.write(l,1,'Locus',generalStyle)
			    cnvSheet.write(l,2,'Type',generalStyle)
			    cnvSheet.write(l,3,'CytoBand',generalStyle)
			    cnvSheet.write(l,4,'Length',generalStyle)
			    cnvSheet.write(l,5,'Copy Number',generalStyle)
			    cnvSheet.write(l,6,'CNV confidence',generalStyle)
			    cnvSheet.write(l,7,'Precision',generalStyle)
			    cnvSheet.write(l,8,'DGV',generalStyle)
			    l_cnv = l
			    #snv header line
			    snvSheet.write(l,0,'Gene',generalStyle)
			    snvSheet.write(l,1,'Transcrit',generalStyle)
			    snvSheet.write(l,2,'Locus',generalStyle)
			    snvSheet.write(l,3,'Type',generalStyle)
			    snvSheet.write(l,4,'Nocall Reason',generalStyle)
			    snvSheet.write(l,5,'Coding',generalStyle)
			    snvSheet.write(l,6,'Amino Acid Change',generalStyle)
			    snvSheet.write(l,7,'Reference',generalStyle)
			    snvSheet.write(l,8,'Genotype',generalStyle)
			    snvSheet.write(l,9,'Allele Frequency',generalStyle)
			    snvSheet.write(l,10,'Coverage',generalStyle)
			    snvSheet.write(l,11,'Allele coverage',generalStyle)
			    snvSheet.write(l,12,'Ref+/Ref-/Var+/Var-',generalStyle)
			    snvSheet.write(l,13,'Info',generalStyle)
			    snvSheet.write(l,14,'Variant ID',generalStyle)
			    snvSheet.write(l,15,'Location',generalStyle)
			    snvSheet.write(l,16,'Variant Effect',generalStyle)
			    snvSheet.write(l,17,'p-value',generalStyle)
			    snvSheet.write(l,18,'Oncomine Variant Class',generalStyle)
			    snvSheet.write(l,19,'Oncomine Gene Class',generalStyle)
			    snvSheet.write(l,20,'Variant ID',generalStyle)
			    snvSheet.write(l,21,'dbSNP',generalStyle)
			    snvSheet.write(l,22,'dbSNP MAF',generalStyle)
			    snvSheet.write(l,23,'5000Exomes MAF',generalStyle)
			    snvSheet.write(l,24,'drugbank',generalStyle)
			    snvSheet.write(l,25,'Sift',generalStyle)
			    snvSheet.write(l,26,'Grantham',generalStyle)
			    l_snv = l
		
			    
				    #for i in range(len(fields)):
				#	    print fields[i]
					    #header_dict[i]=fields[i]
					    #print header_dict[i]
			    head_file.close()
		    except :
			    print "WARNING: no header file found"
		    #ARN sheet
		    try :
			    arn_file = open(tmp_ARN,'r')
			    print 'OK READ TMP_ARN FILE'
			    l_arn = l_arn +1
			    tag_fusion = 0
			    tag_head = 0
			    for arn_line in arn_file:
				    #print arn_line
				    if tag_head == 1:
					    fields = arn_line.split('\t')
					    #ExprControl
					    if fields[1]=='ExprControl':
						    arnSheet.write(l_arn,0,fields[gene],generalStyle)
						    locus = fields[chrom]+':'+fields[pos]
						    arnSheet.write(l_arn,1,locus,generalStyle)
						    arnSheet.write(l_arn,2,'EXPR_CONTROL',generalStyle)
						    arnSheet.write(l_arn,3,fields[v_id],generalStyle)
						    arnSheet.write(l_arn,4,'',generalStyle)
						    arnSheet.write(l_arn,5,fields[v_name],generalStyle)
						    arnSheet.write(l_arn,6,fields[count],generalStyle)
						    arnSheet.write(l_arn,7,fields[fr],generalStyle)
						    arnSheet.write(l_arn,8,'',generalStyle)
						    arnSheet.write(l_arn,9,fields[rpm],generalStyle)
						    l_arn = l_arn + 1
					    #imbalanced
					    elif fields[1]=='5p3pAssays':
						    arnSheet.write(l_arn,0,fields[gene],generalStyle)
						    locus = fields[chrom]+':'+fields[pos]+','+fields[chrom]+':'+fields[R_pos]
						    arnSheet.write(l_arn,1,locus,generalStyle)
						    arnSheet.write(l_arn,2,'ASSAYS_5P_3P',generalStyle)
						    arnSheet.write(l_arn,3,fields[v_id],generalStyle)
						    arnSheet.write(l_arn,4,'',generalStyle)
						    arnSheet.write(l_arn,5,fields[v_name],generalStyle)
						    arnSheet.write(l_arn,6,fields[count],generalStyle)
						    arnSheet.write(l_arn,7,fields[fr],generalStyle)
						    arnSheet.write(l_arn,8,fields[imbalanced],generalStyle)
						    arnSheet.write(l_arn,9,fields[rpm],generalStyle)
						    l_arn = l_arn + 1
					    #known fusion
					    elif fields[1]=='Fusion':
						    if tag_fusion==0 :
							    L_gene = fields[gene]
							    L_exon = fields[exon]
							    L_chrom = fields[chrom]
							    L_pos = fields[pos]
							    tag_fusion = 1
						    #left = fields
						    elif tag_fusion==1 :
							    tag_fusion = 0
							    genes = L_gene+'('+L_exon+')-'+fields[gene]+'('+fields[exon]+')'
							    arnSheet.write(l_arn,0,genes,generalStyle)
							    locus = L_chrom+':'+L_pos+'-'+fields[chrom]+':'+fields[pos]
							    arnSheet.write(l_arn,1,locus,generalStyle)
							    arnSheet.write(l_arn,2,'FUSION',generalStyle)
							    var_id = fields[v_id].split('_')[0]
							    arnSheet.write(l_arn,3,var_id,generalStyle)
							    arnSheet.write(l_arn,4,'',generalStyle)
							    arnSheet.write(l_arn,5,fields[v_name],generalStyle)
							    arnSheet.write(l_arn,6,fields[count],generalStyle)
							    arnSheet.write(l_arn,7,fields[fr],generalStyle)
							    arnSheet.write(l_arn,8,'',generalStyle)
							    arnSheet.write(l_arn,9,fields[rpm],generalStyle)
							    l_arn = l_arn + 1
				    else:
					    fields = arn_line.split('\t')
					    gene = fields.index('FUNC1.gene')
					    chrom = fields.index('CHROM')
					    pos = fields.index('POS')
					    v_id = fields.index('ID')
					    v_name = fields.index('INFO.1.ANNOTATION')
					    count = fields.index('INFO...READ_COUNT')
					    fr = fields.index('FILTER')
					    rpm = fields.index('INFO...RPM')
					    imbalanced = fields.index('INFO.1.5P_3P_ASSAYS')
					    R_pos = fields.index('INFO.1.THREE_PRIME_POS')
					    exon = fields.index('FUNC1.exon')
					    tag_head = 1
					    
			    arn_file.close()	
		    except :
			    print "!!!!!!!WARNING: could not open tmp_ARN file"
		    #CNV sheet
		    try :
			    cnv_file = open(tmp_CNV,'r')
			    print 'OK READ TMP_CNV FILE'
			    l_cnv = l_cnv +1
			    tag_head = 0
			    for cnv_line in cnv_file:
				#print cnv_line
				    if tag_head == 1:
					    fields = cnv_line.split('\t')
					    cnvSheet.write(l_cnv,0,fields[gene],generalStyle)
					    cnvSheet.write(l_cnv,1,fields[locus],generalStyle)
					    cnvSheet.write(l_cnv,2,'CNV',generalStyle)
					    cytoband = fields[iscn].split('x')[0]
					    cnvSheet.write(l_cnv,3,cytoband,generalStyle)
					    cnvSheet.write(l_cnv,4,fields[length],generalStyle)
					    copy_number = fields[iscn].split('x')[-1]
					    cnvSheet.write(l_cnv,5,copy_number,generalStyle)
					    cnvSheet.write(l_cnv,6,fields[confidence],generalStyle)
					    cnvSheet.write(l_cnv,7,fields[precision],generalStyle)
					    cnvSheet.write(l_cnv,8,fields[dgv],generalStyle)
					    l_cnv = l_cnv + 1
				    else :
					    fields = cnv_line.split('\t')
					    gene = fields.index('gene')
					    locus = fields.index('# locus')
					    iscn = fields.index('iscn')
					    length = fields.index('length')
					    confidence = fields.index('confidence')
					    precision = fields.index('precision')
					    dgv = fields.index('dgv')
					    tag_head = 1
			    cnv_file.close()	
		    except :
			    print "!!!!!!!WARNING: could not open tmp_CNV file"
		    #SNV sheet
		    try :
			    snv_file = open(tmp_SNV,'r')
			    print 'OK READ TMP_SNV FILE'
			    l_snv = l_snv +1
			    tag_head = 0
			    for snv_line in snv_file :
				    if tag_head == 1 :
					    #print snv_line
					    fields = snv_line.split('\t')
                                            snvSheet.write(l_snv,0,fields[gene],generalStyle)
                                            snvSheet.write(l_snv,1,fields[transcript],generalStyle)
                                            locus = fields[chrom]+':'+fields[pos]
                                            snvSheet.write(l_snv,2,locus,generalStyle)
                                            if fields[call] == 'NOCALL':
                                                    row_type = 'NOCALL '+ fields[rowtype]
                                                    snvSheet.write(l_snv,3,row_type,generalStyle)
                                                    snvSheet.write(l_snv,4,fields[nocall_reason],generalStyle)
                                            else :
                                                    snvSheet.write(l_snv,3,fields[rowtype],generalStyle)
                                                    snvSheet.write(l_snv,4,'',generalStyle)
					    snvSheet.write(l_snv,5,fields[coding],generalStyle)
                                            snvSheet.write(l_snv,6,fields[aa_coding],generalStyle)
                                            snvSheet.write(l_snv,7,fields[ref],generalStyle)
                                            snvSheet.write(l_snv,8,fields[geno],generalStyle)
                                            snvSheet.write(l_snv,9,fields[allele_freq],generalStyle)
                                            snvSheet.write(l_snv,10,fields[cov],generalStyle)
                                            snvSheet.write(l_snv,11,fields[allele_cov],generalStyle)
                                            strand_ref_alt = fields[ref_forward]+'/'+fields[ref_reverse]+','+fields[alt_forward]+'/'+fields[alt_reverse]
                                            snvSheet.write(l_snv,12,strand_ref_alt,generalStyle)
                                            if fields[hs]=='True':
                                                    info = 'Hotspot'
                                            else:
                                                    info = ''
					    snvSheet.write(l_snv,13,info,generalStyle)
                                            snvSheet.write(l_snv,14,fields[v_id],generalStyle)
                                            snvSheet.write(l_snv,15,fields[location],generalStyle)
                                            snvSheet.write(l_snv,16,fields[function],generalStyle)
                                            snvSheet.write(l_snv,17,fields[qd],generalStyle)
                                            try:
						    snvSheet.write(l_snv,18,fields[var_class],generalStyle)
					    except:
						    snvSheet.write(l_snv,18,'',generalStyle)
					    try:
						    snvSheet.write(l_snv,19,fields[gene_class],generalStyle)
					    except:
						    snvSheet.write(l_snv,19,'',generalStyle)
                                            snvSheet.write(l_snv,20,fields[v_id],generalStyle)
                                            snvSheet.write(l_snv,21,'',generalStyle)
                                            snvSheet.write(l_snv,22,'',generalStyle)
                                            snvSheet.write(l_snv,23,'',generalStyle)
                                            snvSheet.write(l_snv,24,'',generalStyle)
                                            snvSheet.write(l_snv,25,'',generalStyle)
                                            snvSheet.write(l_snv,26,'',generalStyle)
                                            
					    l_snv = l_snv +1
					    #print 'SNV line number : %s ' % l_snv
				    else :
					    fields = snv_line.split('\t')
					    for field in fields :
						    print field
                                            gene = fields.index('FUNC1.gene')
					    print 'gene index %s' %gene
					    chrom = fields.index('CHROM')
					    print 'chrom index %s' %chrom
                                            pos = fields.index('POS')
					    print 'pos index %s'%pos
                                            v_id = fields.index('ID')
					    print ' v_id index %s' %v_id
                                            v_name = fields.index('INFO.1.ANNOTATION')
					    print ' v_name index %s' %v_name
					    exon = fields.index('FUNC1.exon')
					    print ' exon index %s' %exon
					    transcript = fields.index('FUNC1.transcript')
					    print 'transcript index %s' %transcript
					    rowtype = fields.index('rowtype')
					    print 'rowtype index %s' %rowtype
					    call = fields.index('call')
					    print 'call index %s' %call
					    nocall_reason = fields.index('INFO...FR')
					    print 'nocall_reason index %s' %nocall_reason
					    coding = fields.index('FUNC1.coding')
					    print 'coding index %s' %coding
					    aa_coding = fields.index('FUNC1.protein')
					    print 'aa_coding index %s' %aa_coding
					    ref = fields.index('REF')
					    print 'ref index %s' %ref
					    geno = fields.index('FORMAT.1.GT')
					    print 'geno index %s' %geno
					    allele_freq = fields.index('INFO.A.AF')
					    print 'allele_freq index %s' %allele_freq
					    cov = fields.index('INFO.1.FDP')
					    print 'cov index %s' %cov
					    allele_cov = fields.index('INFO.A.FAO')
					    print 'allele_cov index %s' %allele_cov
					    ref_forward = fields.index('INFO.1.FSRF')
					    print 'ref_forward index %s' %ref_forward
					    ref_reverse = fields.index('INFO.1.FSRR')
					    print 'ref_reverse index %s' %ref_reverse
					    alt_forward = fields.index('INFO.A.FSAF')
					    print 'alt_forward index %s' %alt_forward
					    alt_reverse = fields.index('INFO.A.FSAR')
					    print 'alt_reverse index %s' %alt_reverse
					    hs = fields.index('INFO.0.HS')
					    print 'hs index %s' %hs
					    location = fields.index('FUNC1.location')
					    print 'location index %s' %location
					    function = fields.index('FUNC1.function')
					    print 'function index %s' %function
					    qd = fields.index('INFO.1.QD')
					    print 'qd index %s' %qd
					    try:
						    var_class = fields.index('FUNC1.oncomineVariantClass')
					    except:
						    var_class = ''
					    #print 'var_class index %s' %var_class
					    try:
						    gene_class = fields.index('FUNC1.oncomineGeneClass')
					    except:
						    gene_class = ''
					    #print 'gene_class index %s' %gene_class
					    tag_head = 1
			    snv_file.close()	
		    except :
			    print "!!!!!!!WARNING: could not open tmp_SNV file"
		    #save excel file
		    print '!!!!!!!!!! EXCEL OUT SAVING !!!!!!!!!!!!!!!'
		    print 'output excel %s%s.xls' % (output_path,sampleName)
		    finalReport.save("%s%s.xls" % (output_path, sampleName))    
	    else :
		    print 'No IR uploader successfull instance found for sample '+sampleName

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
			htmlOut = open('%s/plugin_out/BackupIR_out/BackupIR_block.html'%self.envDict['ANALYSIS_DIR'], 'a+')
		except:
			htmlOut = open('BackupIR_block.html', 'w')
		htmlOut.write('<html><body>\n')

		if os.path.isfile(self.envDict['TSP_FILEPATH_BARCODE_TXT']):
			self.isBarcodedRun = True
			
		if not self.isBarcodedRun:
			print "This plugin only work for barcoded run. \n Exiting..."
			sys.exit(0)

                # get the project name.
		try:
			self.project = self.envDict['TSP_PROJECT']
		except:
			self.project = 'unknown'
			
		# get the run name.
		try:
			self.runName = self.envDict['TSP_ANALYSIS_NAME']
		except:
			self.runName = 'unknown'
		self.runPath = self.project + '/' + self.runName

                # creating folder of the project
		print 'CREATING FOLDER : %s -'%(self.project)
		cmd = Popen(['mkdir', self.folder+'/'+self.project], stdout=PIPE)
		out, err = cmd.communicate()
		#print 'OUT: %s\nERR: %s'%(out, err)

                # creating folder of the run
		print 'CREATING FOLDER : %s -'%(self.runName)
		cmd = Popen(['mkdir', self.folder+'/'+self.runPath], stdout=PIPE)
		out, err = cmd.communicate()
		#print 'OUT: %s\nERR: %s'%(out, err)

                #get IonReporter uploader path
                self.IRuploaderPath = self.getPluginPath("IonReporterUploader")
                htmlOut.write('&#8594; IonReporterUploader instance : <b>%s</b><br/>\n'%self.IRuploaderPath)
                if not self.IRuploaderPath:
                    IRuploader = False

		summary_path = '%s/plugin_out/%s/summary.html'%(self.envDict['ANALYSIS_DIR'],self.IRuploaderPath)
		print "SUMMARY "+summary_path
		summary = open(summary_path,'r')
		#IR_samples=[]
		for line in summary:
			#print line
			if 'Valid&nbsp;Samples' in line:
				line = line.strip().split("<br>")[0]
				sample_ir = line.split("&nbsp;")[-1]
				#print "SAMPLE IR " + sample_ir
				sample_ir = sample_ir.replace('RNA_','')
				#print "SAMPLE IR replace " + sample_ir
				if sample_ir not in self.IR_samples :
					self.IR_samples.append(sample_ir)
				#else:
				#	print sample_ir + 'already in list'
	

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
				sampleName = ''
				#self.sampleNameLookup[barcodeName] = '' # makes it much easier later to do the lookup with no conditional tests
			# MGD note: I considered setting blank sample names to the barcode name instead, but might not be what customer intended
		# Create patients folders
		print 'CREATING PATIENTS FOLDERS...'
		self.createOneFolderByPatient(bamPaths)	
		self.sampleNameLookup[''] = '' # allows us to easily handle case where barcode might not have been found

               

if __name__ == "__main__":
	PluginCLI(BackupIR())

