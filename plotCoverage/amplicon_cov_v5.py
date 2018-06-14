#!/usr/bin/env python

import sys
import os
import subprocess
import json
import time
import traceback
import re
import glob
import ast
import csv
import operator
import math
from optparse import OptionParser
import matplotlib
matplotlib.use('Agg')
 # tester mettre plus loin
from fnmatch import fnmatch
#from django.conf import settings
#from django.template.loader import render_to_string


parser = OptionParser()
parser.add_option('-r', '--result-dir', help='coverageAnalysis plugin output dir', dest='result_dir')
parser.add_option('-c', '--coverage-path', help='coverageAnalysis plugin output dir', dest='coverage_path')
parser.add_option('-t', '--target-path', help='coverageAnalysis plugin output dir', dest='target_path')
 
(options, args) = parser.parse_args()

result_dir  = options.result_dir #coverageAnalysis result dir path
coverage_path  = options.coverage_path #coverageAnalysis result dir path
target_path  = options.target_path #coverageAnalysis result dir path
print coverage_path

################################################################
#find samples amplicon.cov.xls files in IonXpress subdirectories in coverageAnalysis output dir

cov_pattern = "*.amplicon.cov.xls"
#bam_pattern = "*.bam"
#list of cnotaining amplicon.cov.xls files for all samples
cov_files = []
#list of cnotaining bam files for all samples
#bam_files = []

for dirPath, subdirName, fileName in os.walk(coverage_path):
    #print('Found directory: %s' % dirPath)
    for fname in fileName:
        if fnmatch(fname,cov_pattern):
            if not fnmatch (fname,"link.amplicon.cov.xls"):
                #print os.path.join(dirPath,fname)
                cov_files.append(os.path.join(dirPath,fname))
        #if fnmatch(fname,bam_pattern):
            #if not fnmatch (fname,".bam.bai"):
                #print os.path.join(dirPath,fname)
                #bam_files.append(os.path.join(dirPath,fname))
 
#if len(bam_files) >1:
    #parse sample bam file 
    #for file_path in bam_files:
        #print file_path
        #dir_path,sample_file = os.path.split(file_path)
        #print dir_path
        #print sample_file

        #file_extension="(.*).amplicon.cov.xls"
        #fields = re.match(file_extension,sample_file)
        #sample_id=os.path.basename(dir_path)
        
        #print sample_id


#samples id list
samples_list = []

#chromosomes list
chromosomes = ''

#read amplicon_cov file to find amplicon with less than 500 reads mapped
#hash key= "amplicon id" value="sample name(total reads mapped on amplicon)"
amplicon_sample_cov = {}

if len(cov_files) >1:
    #hash key="sample name" value="string tab separated of total reads per amplicon"
    amplicon_cov_for_sample = {}

    #hash key="amplicon id" value="string tab separated of total reads per amplicon"
    amplicon_cov = {}
    
    #parse sample amplicon cov file 
    for file_path in cov_files:
        #print file_path
        dir_path,sample_file = os.path.split(file_path)
        #print dir_path
        #print sample_file

        #file_extension="(.*).amplicon.cov.xls"
        #fields = re.match(file_extension,sample_file)
        sample_id=os.path.basename(dir_path)
        
        print sample_id

        #process less than 100X coord
        print('find position less than 100X detph')
        prefix = sample_file.split(".")[0]
        #print prefix
        bam_path = dir_path +"/"+prefix+".bam"
        #print bam_path

        """#>>>>>>>>>>>> decommenter et tester une fois le fichier panel recuperer
        # !!!! TODO : retrieve panel file path
        sample_cov = result_dir+"/"+sample_id+"_cov.out"
        cmd = "bedtools coverage -d -abam bam_path -b target_path > sample_cov"

        cov_100x_path = result_dir+"/"+sample_id+"_cov_less_than_100x.out"
        cov_100x_file = open(cov_100x_path,'w')
        chr_pos_path = result_dir+"/"+sample_id+"_chr_pos_100x_max.out"
        chr_pos_file =open(chr_pos_path,'w')
        
        sample_cov_file = open(sample_cov,'r')
        for line in sample_cov_file:
            line = line.rstrip('\n')
            chrom = line.split('\t')[0]
            start = line.split('\t')[1]
            end = line.split('\t')[2]
            name = line.split('\t')[3]
            pos_in_region = line.split('\t')[8]
            detph = line.split('\t')[9]

            if depth<100:
                cov_100X_file.write(line+"\n")
                pos_in_chrom = int(start)+int(pos_in_region)
                chr_pos_file.write(chrom+"\n"+pos_in_chrom+"\n")
        """#<<<<<<<<<<<<< fin decommenter


        # !!! TODO: recuperer le sample_name associe au sample_id et l'utiliser dans la suite du script a la place de sample_id
           
        '''
        #read startplugin_json file, extract names associated with barcodes
        
        sample_name = {}

        try:
        json_file = open(os.path.join(ANALYSIS_DIR,'ion_params_00.json'), 'r')
        ion_params_json = json.load(json_file,parse_float=str)
        json_file.close()
        except:
        printtime('ERROR: Failed to load and parse ion_params_00.json')
        return 1
    
        #put results in a new dictionnary
        barcodedSamples=ast.literal_eval(ion_params_json['experimentAnalysisSettings'].get('barcodedSamples'))
        #go through the dictionnary and revert it to use the barcode as key, assuming there is only one barcode per samplesname
        for key in barcodedSamples:
        barcode=barcodedSamples.get(key).get('barcodes')[0]
        sample_name[barcode]=key
        '''

        # !!! TODO: decommenter >>>>
        #control = 'BLANC PCR'
        #if control in sample_name :
        #    print 'H2O control not analyzed'
        #    continue
        #add sample to the list of samples processed
        #samples_list.append(sample_name)
        # <<<< fin decommenter 

        sample_name = sample_id
        samples_list.append(sample_name)

        #amplicon id list
        amplicon_list = ''

        #find amplicon with less than 500 reads mapped on
        #sample_name = sample_name.replace(" ","\ ")
        #file_path = file_path.replace(" ","\ ")
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
        
        #chromosomes = ''
        #cmd = "awk '{print $1}' " + dir_path + "/" + sample_name + "_amplicon.cov.sorted.tsv | sort -u"
        #print cmd
        #chromosomes = os.popen(cmd).read()
        #print  "CHROM:" + chromosomes
        
        ###############################
        #find amplicons cov for sample
        
        chromosomes = ''
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

                total_reads = log10(float(total_reads))

                total_reads = str(total_reads)

                #retrieve chromosomes list
                chr_id = line.split('\t')[0]
                if (previous_chr is None):
                    chromosomes = chr_id
                    previous_chr = chr_id
                else:
                    if (chr_id == previous_chr):
                        next
                    else:
                        chromosomes = chromosomes + "\t" + chr_id
                        previous_chr = chr_id

                if sample_name in amplicon_cov_for_sample.keys():
                    amplicon_cov_for_sample[sample_name] = amplicon_cov_for_sample[sample_name] + "\t" + total_reads
                else:
                    amplicon_cov_for_sample[sample_name] = total_reads
                    print "dict " + amplicon_cov_for_sample[sample_name]

                if amplicon_list is not None:
                    amplicon_list = amplicon_list + "\t" + amplicon_id
                else:
                    amplicon_list = amplicon_id

                    
                #hash key=amplicon value=cov
                if amplicon_id in amplicon_cov.keys():
                    amplicon_cov[amplicon_id] = amplicon_cov[amplicon_id] + "\t" + total_reads
                else:
                    amplicon_cov[amplicon_id] = total_reads

        sorted_file.close()
        #print chromosomes
        #print amplicon_list

######################
#write result file 

#tab file containing header=region_id and lines=samples amplicon total reads for which that have less than 500 reads
uncov_path = result_dir + "/amplicon.uncov.tsv"
uncov_file = open(uncov_path,'w')

#print header = amplicon id
for sample in samples_list:
    uncov_file.write("\t" + sample)
uncov_file.write("\n")
#values for amplicon
for amplicon_id in amplicon_sample_cov.keys():
    uncov_file.write(amplicon_id)
    for sample in samples_list:
        uncov_file.write("\t" + amplicon_sample_cov[amplicon_id][sample])
    uncov_file.write("\n")


#list of files containing amplicon log10(total reads) : one file with all amplicons and 4 others (subset of amplicon to have less data in 1 plot to interpret easier)
files_list = []

#tab file containing header=region_id and lines=samples amplicon total reads for all amplicons and all samples
all_cov_path = result_dir + "/all_amplicon_samples.cov.tsv"
all_cov_file = open(all_cov_path,'w')

#print header = amplicon id
all_cov_file.write("sample" + amplicon_list + "\n")

#print samples data
for sample in samples_list:
    #print sample
    all_cov_file.write(sample + "\t" + amplicon_cov_for_sample[sample] + "\n")
all_cov_file.close()

files_list.append(all_cov_path)


#subset data files all amplicov cov for 10 samples max
chunks = [samples_list[i:i+10] for i in xrange(0,len(samples_list),10)]
sample_subset = 0
for chunk in chunks:
    sample_subset = sample_subset + 1
    sample_subset_path = result_dir+"/all_amplicon-sample_subset"+str(sample_subset)+".tsv"
    sample_subset_file = open(sample_subset_path,'w')
    sample_subset_file.write("sample" + amplicon_list + "\n")
    for sample in chunk:
        sample_subset_file.write(sample + "\t" + amplicon_cov_for_sample[sample] + "\n")
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
    cmd = "cut -f 1,"+str(start)+"-"+str(end)+" "+all_cov_path+" > "+result_dir+"/amplicon_subset"+str(i)+"-all_samples.tsv"
    #print cmd + "\n"
    os.system(cmd)
    subset_file = result_dir+"/amplicon_subset"+str(i)+"-all_samples.tsv"
    files_list.append(subset_file)
    
    #subset files containing 10 samples max and amplicons_number/4 (mandatory by IH)
    sample_subset = 0
    for chunk in chunks:
        sample_subset = sample_subset + 1
        sample_subset_path = result_dir+"/all_amplicon.sample_subset"+str(sample_subset)+".tsv"
        cmd = "cut -f 1,"+str(start)+"-"+str(end)+" "+sample_subset_path+" > "+result_dir+"/amplicon_S"+str(i)+"-samples_S"+str(sample_subset)+".tsv"
        #print cmd + "\n"
        os.system(cmd)
        sub_file_path = result_dir+"/amplicon_S"+str(i)+"-samples_S"+str(sample_subset)+".tsv"
        files_list.append(sub_file_path)
    
    start = start + nb_col_subset
    end = end + nb_col_subset
    

#############################
#create amplicons cov graphes

for f in files_list:
    print f + "\n"
    
    #figure(figsize=(150,50)

    fname = f.split("/")[-1]
    prefix = fname.split(".")[0]
    
    print "FILE NAME = "+fname+"\n"
    print "PREFIX = "+prefix+"\n"

    cov_file = open(f,'r')

    for line in cov_file:
        #print line + "\n"
        if line.startswith('sample'):
            #print line + "\n"
            x_names = line.split('\t')
            x_names.pop(0) #remove first element of the list (header line, remove "sample" label and retreive amplicon id for x_axis)
            #print "X"
            #print len(x_data)
            x_len = len(x_names)
        
            x_data = range(len(x_names))
            x = array(x_data)

            #x_data1 = range(0,20)
            #x1 = array(x_data1)

        else:
            #if nb_amplicon<20:
            y_data = line.split('\t')
            y_legend = y_data[0] #first field of the line = sample name for legend
            y_data.pop(0)
            #print "Y"
            #print len(y_data)
            
            #y = array(y_data)
            #plot(x,y,label=y_legend,linewidth=2)
            
            y2_data =[]
            #nb_amplicon = 0
            #nb_file = 0
            fig=figure(1,figsize=(50,20))
            rcParams['axes.color_cycle']=['b','green','r','c','maroon','gold','k','sienna','orange','palevioletred','grey','darkkhaki','darkorchid','lime','magenta','dodgerblue','lightcoral','lightblue','yellow','plum','peru','steelblue','mediumspringgreen','yellowgreen','tomato','darkseagreen','burlywood','deeppink','slategrey','greenyellow']
            #rcParams['lines.linestyle_cycle']=['-','--','.','-.']
            for n in y_data:
                #if nb_amplicon<20:
                #n = int(n)
                #log_n = math.log10(n)
                #y2_data.append(log_n)
                #y2_data.append(int(n)/1000+nb_ech)
                y2_data.append(n)

                #def transform(n): return n+nb_ech
                #y2_data = map(transform(y_data))
                #for i in y2_data:
                #    i = i
                
            y2 = array(y2_data)
            plot(x,y2,label=y_legend,linewidth=2)
                #data = concatenate(x_names,0)
                #boxplot(x)
            #if nb_amplicon==20:
            #    nb_file = nb_file + 1
            #    nb_amplicon = 0
             #   width = 1000
            #    height = 500
    cov_file.close()
    width = 200
    height = 5000
    xticks(x_data,x_names,rotation=90,fontsize=20)
    yticks(fontsize=20)
    ylim(ymax=5)
    xlabel("amplicons",fontsize=22)
    ylabel("total reads (log10)",fontsize=22)
    axhline(y=3,color='black',linestyle='--',label='1000')
    axhline(y=2.7,color='b',linestyle='--',label='500')
    axhline(y=2.48,color='orange',linestyle='--',label='300')
    axhline(y=2,color='r',linestyle='--',label='100')
    legend(fontsize=20,loc='center left',bbox_to_anchor=(1,0.5)) #legend(loc='center left',bbox_to_anchor=(1,0.5))
    autoscale(enable=True)
    tight_layout(rect=[0,0,0.9,1])
    
    #show()
                

    #!!!!!!TODO retreive prefixe and create file name
    #output file graph
    file_name = prefix+'_plot.png'
    img_path = os.path.join(result_dir,file_name)
    fig.savefig(img_path,dpi=200)
    #savefig(img_path)
    close()

