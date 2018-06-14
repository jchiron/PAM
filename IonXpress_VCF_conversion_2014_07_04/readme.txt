##  C.Lucchesi Jul.4, 2014
##
##  BASH script to convert LifeTech Torrent Server .vcf files into standard samtools / bcftools .bcf files
##    and to annotate the variant lines
##
Launching the script with no parameters will provide help on how to launch the script

sh ./Convert_Annotate_TorSer_VCF.V1.sh

#########################################
## SCRIPT to convert IonXpress VCF files #
##  into samtools compatible BCF files   #
##  and annotate variants                #
##########################################

Mandatory parameters:
        -P      Full path IonXpress Samples Folder (path begins with /)
        -L      Full path IonXpress Samplesplan    (path begins with /, 1st column is sample name - mandatory)
        -A      Full path Annotation Folder        (path begins with /, if it does not exist it is created
        -U      Full path UCSC reference genome    (path begins with /, file name without .fa  suffix)
        -X      Visualize results with less -S     (yes / no)
        -C      Confirm                            (yes / no)

Optional parameters:
        -1      TorrentServer VCF code for the STRAND FORWARD REF ALLELE COUNT  in  INFO field (default: SRF)
        -2      TorrentServer VCF code for the STRAND REVERSE REF ALLELE COUNT  in  INFO field (default: SRR)
        -3      TorrentServer VCF code for the STRAND FORWARD ALT ALLELE COUNT  in  INFO field (default: SAF)
        -4      TorrentServer VCF code for the STRAND REVERSE ALT ALLELE COUNT  in  INFO field (default: SAR)



#####################  FOLDER STRUCTURE  #############################""
A "sample" read-me study folder is provided in this same directory. 
As you can see from the parameters, the annotation process IS DRIVEN BY a SAMPLESPLAN file that can contain 1 or more samples.
The file has 1 header line (skipped) and 1 line per sample. Only the 1st field is processed (the sample name)
If the sample name begins with #, then the line is skipped.
If more than 1 sample is specified then a single annotation file is produced that merges the variants from each VCF sample file.

WARNING !!!!!!!!!!!!!
Sample names MUST not contain "_" characters. Otherwise the sample name is the set of characters till the first "_"

The structure of the study is as follows
./samples                                    ##  is the folder where samples are stored 
            /samplesplan_test.txt            ##  contains 1 line for each sample  to convert and annotate
            /IonXpress005                              ##  sample folder - one for each sample, folder name is sample name
                         /IonXpress005.vcf     ##  VCF in TorrentServer format - VCF name is the sample name
                         /BCF                  ##  converted VCF to BCF files are registered here - BCF names are prefixed by sample names

            /IonXpress006                              ##  sample folder - one for each sample, folder name is sample name
                         /IonXpress006.vcf     ##  VCF in TorrentServer format
                         /BCF                  ##  converted VCF to BCF files are registered here

              
            /ANNOTATE_5_6                              ##  folder where annotations are produced for the pool of sample 5 and 6
                         /Union_samples_chrAll.txt ##  the tab delimited annotation file 

#####################   COMMAND LINE and EXAMPLE DATASET ##############
Launch as follows:

sh ./Convert_Annotate_TorSer_VCF.V1.sh -P `pwd`/samples                      \
                                       -L `pwd`/samples/samplesplan_test.txt \
                                       -A `pwd`/samples/ANNOTATE_5_6         \
                                       -U ${UCSC_genome}                     \
                                       -X yes                                \
                                       -C yes

The variable ${UCSC_genome} must point to the full path name of the ref fasta Hgxx without .fa

#######################  THE LOGIC BEHIND  ########################
LifeTech Torrent Server uses specific fields for the depth of ref and alt alleles in the VCF outfile. The following has been found in the 
Life Tech documentation.
STRAND FORWARD REF ALLELE COUNT is SRF in INFO fields
STRAND REVERSE REF ALLELE COUNT is SRR in INFO fields
STRAND FORWARD ALT ALLELE COUNT is SAF in INFO fields
STRAND REVERSE ALT ALLELE COUNT is SAR in INFO fields

The bash script scans the INFO field of each VCF position looking for those fields; it then adds the DP4 filed at the end of the INFO field as:
                                DP4=SRF,SRR,SAF,SAR
Then it splits the VCF by chromosome and creates BCF files for each chromosome
BCF files are indexed
This is done for each sample in the samplesplan file

Then, the BCF of the samples in samplesplan are annotated chromosome by chromosome. 
An annotated variant file is created, with as many lines as the variants that are present in at least 1 sample BCF.
The name of the annotation file is
Union_samples_chrAll.txt

#####################   ENVIRONMENT VARIABLES  ######################
The variable ${UCSC_genome} must point to the full path name of the ref fasta
Hgxx without .fa
The variable ${DNASEQ_VARIANTS} must point to the folder where the Variant
Annotation scripts are located
The variable ${IonXpress_ANNOTATION} must point to the folder where the main
scripts are located


