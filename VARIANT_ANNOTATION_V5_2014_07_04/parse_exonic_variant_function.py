#!/usr/bin/env python

import sys
import csv
import os

filename = sys.argv[1]
pwd = sys.argv[2]
output = open(pwd+'/Exonic_variant_function_temp.txt',"w")
exonic_variant_function = open(pwd+'/'+filename,"r")

reader = csv.reader(exonic_variant_function,delimiter="\t")

for row in reader:
	protcons = row[1]
	if " " in protcons:
		s = protcons.split(" ")
		protcons1 = s[0]
		protcons2 = ""
		for word in s[1:]:
			protcons2 = protcons2 + word
		#for i in range(len(s[1:])):			# A TESTER SI CAS 3 MOTS APPARAIT, ex:"frameshift block substitution"
			#if i==0:
				#protcons2 = word
			#elif i>0:
				#protcons2 = protcons2 + " " + word
	else:
		protcons1 = protcons
		protcons2 = "-"
	output.write("%s:%s:%s:%s:%s\t%s\t%s\t%s\n" % (row[3],row[4],row[5],row[6],row[7],protcons1,protcons2,row[2]))

exonic_variant_function.close()
output.close()

sys.exit()
