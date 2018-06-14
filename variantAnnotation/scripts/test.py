import sys
import os

file_output_path = "/results/analysis/output/Home/PGM194_avec_BC_OK_278/plugin_out/variantAnnotation_out.1212/diag/MORISSONNEAU-kl353_NGS_Diag_VC_TS.tsv"

#add to sort mutations by type and frequencies in categories
#file_sort_tmp = plugin_path + '/' + sample_id + '/CHR/' + sample_name + '_NGS_Diag_sort.tmp'
file_sort_path = "/results/plugins/variantAnnotation/scripts/MORISSONNEAU-kl353_NGS_Diag_sort.tsv"

#header
cmd = "grep Comm.Bio " + file_output_path +" >" + file_sort_path
os.system(cmd)
#category 1 : exonic frameshift, missense, nonsense and sort decreasing frequency
cmd = "grep -v -E \"splicing|synonymous|intronic|UTR|ncRNA|Comm.Bio\" " + file_output_path + '|sort -t \'\t\' -rn -k10,10 >> ' + file_sort_path
os.system(cmd)
cmd = "echo >>" + file_sort_path
os.system(cmd)

#category 2 : splicing and sort decreasing frequency
cmd = "grep splicing " + file_output_path + "|sort -t \'\t\' -rn -k10,10 >> " + file_sort_path
os.system(cmd)
cmd = "echo >>" + file_sort_path
os.system(cmd)

#category 3 : exonic synonymous and sort decreasing frequency
cmd = "grep synonymous " + file_output_path + "|sort -t \'\t\' -rn -k10,10 >> " + file_sort_path
os.system(cmd)
cmd = "echo >>" + file_sort_path
os.system(cmd)

#category 4 : intronic and sort decreasing frequency
cmd = "grep intronic " + file_output_path + "|sort -t \'\t\' -rn -k10,10 >> " + file_sort_path
os.system(cmd)
cmd = "echo >>" + file_sort_path
os.system(cmd)

#category 5 : UTR and sort decreasing frequency
cmd = "grep UTR " + file_output_path + "|sort -t \'\t\' -rn -k10,10 >> " + file_sort_path
os.system(cmd)
cmd = "echo >>" + file_sort_path
os.system(cmd)

#category 6 : ncRNA and sort decreasing frequency
cmd = "grep ncRNA " + file_output_path + "|sort -t \'\t\' -rn -k10,10 >> " + file_sort_path
os.system(cmd)
    
