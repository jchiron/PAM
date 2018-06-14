# tlesluyes, clucchesi 2014_01_13		
# clucchesi, ktran 2014_02_09		
# clucchesi 2014_03_25		
		
This version of the SNV / InDel VARIANT pipeline performs:		
Union of BCF files and construction of variant table for a set of samples defined in a samplesplan.
Two scripts must be called for that:

sh ./VarPipe_process_batch_samples.V2.sh

########################################
#Script annotate variants in BCF files
#########################################
Mandatory parameters:
        -S      Samplesplan                  (each line contains a sample to annotate and to merge with the others,
                                              the first col is the sample name, only the first column is mandatory)
        -I      Input directory              (where sample folders are stored)
        -O      Output directory             (where annotated variants are stored - if it does not exist, then it is created)
        -B      Use BCF_cleaned or BCF dir   (values: BCF or BCF_cleaned)
        -D      Size of genomic interval     (used to split each BCF file in equal sized intervals - each interval is annotated in parallel
                                              integer >0  if no split is desidered than set to 999999999
                                              Proposed size: for a trio Constit + Tumor DNA and T RNA Whole Genome --> size=10000000
                                              This is roughly equivalent to 3 300 000 000 / 10 000 000 = 330 jobs that are run in parallel)
        -V      Bcftools version             (values: 0.1.19 / 0.2.0)
        -C      Confirm launch               (values: yes / no)

and 

sh ./MergeVar_all_CHRs.V2.sh

########################################################
#Script to concatenate variant tables of each chromosomes
#########################################################
Mandatory parameters:
        -V      Full path of Variant directory  (where chromosome variant tables are stored - must begin with /)
        -A      Aggregate variants     (yes or no)
        -S      Split by chromosome    (yes or no)
        -Z      Compress .gz           (yes or no)
        -C      Confirmation to launch (yes or no)


One optional script

sh ./check_variant_annotation_queue.sh

########################################################
#Script to check the status of annotation jobs
#########################################################
Mandatory parameters:
        -V      Full path of Variant directory  (where chromosome variant tables are stored - must begin with /)

allows to test whether all jobs have finished (return code 0) or not (return code 1)
 


#########################################
The main entry point of the pipeline is the script "VarPipe_process_batch_samples.sh", that produces an annotated variant file for each chromosome of all the samples in the batch; 
the execution is driven by the samplesplan, a tab delimited text file of the form:

sample  stranded reference       proc    mem    
IB111P2  yes     transcriptome   12      40      
#IB111P4 yes     transcriptome   12      40      

The first line is the header; the first col is the sample name, only the first column is mandatory - when the the first column is dashed the sample is skept.

The script "MergeVar_all_CHRs.sh" is used to merge the variants of all chromosomes in a single file 
Mandatory parameters:
        -V      Variant directory  (where chromosome variant tables are stored)

##########################  ENV VARIABLES  #####################
To use this pipeline your .bashrc must include a kind of:
 		
export qsubMail="-m abe -M c.lucchesi@bordeaux.unicancer.fr"    ##  TORQUE mailing

#######  other modules available via PATH
export BRIO_ROOT=/home/carlucchesi/BRIO                    #  at CBiB: /home/clucchesi/BRIO

####### modules to add   OR  PATHs to modify
module load gcc
module add gcc/4.6.0
module load torque
module add samtools/0.1.19
module add lftp/4.3.8
module add qualimap/0.7.1

####### BRIO R library (enables installation of adds-on packages
export R_LIBS=${BRIO_ROOT}/tools/BRIO_R_library_for_locally_installed_packages
export picardPATH=${BRIO_ROOT}/tools/picard-tools-1.99

PATH=$PATH:${BRIO_ROOT}/tools/bowtie2-2.1.0
PATH=$PATH:${BRIO_ROOT}/tools/HTSeq-0.5.4p5/build/scripts-2.6
PATH=$PATH:${BRIO_ROOT}/tools/cufflinks-2.1.1.Linux_x86_64
PATH=$PATH:${BRIO_ROOT}/tools/FastQC
PATH=$PATH:${BRIO_ROOT}/tools/sickle-master
PATH=$PATH:${BRIO_ROOT}/tools/bedtools-2.17.0/bin
PATH=$PATH:${BRIO_ROOT}/tools/seqtk-master
PATH=$PATH:${BRIO_ROOT}/tools/SeqPrep-master
PATH=$PATH:${BRIO_ROOT}/tools/R-3.0.2/bin

PATH=$PATH:${BRIO_ROOT}/tools/annovar_2013_10_31
PATH=$PATH:${BRIO_ROOT}/tools/cutadapt-1.3/bin
PATH=$PATH:${BRIO_ROOT}/tools/tophat-2.0.10.Linux_x86_64

#######  SET variables for analysis folders
export ANNOVARDB=${BRIO_ROOT}/annotation/Annovar_DB_2013_05_20
export UCSC_genome=${BRIO_ROOT}/annotation/UCSC_genome_Hg19_Spikein_2013_08_22/G_Hg19_Spikein

export UCSC_transc=${BRIO_ROOT}/annotation/UCSC_transcriptome_Hg19_NEW_2014_01_20/T_Hg19
export UCSC_annot=${BRIO_ROOT}/annotation/GTF_2013_08_22/UCSC_Hg19_curated_NEW.gtf

#######  PIPELINES
export RNASEQ_VARIANTS=${BRIO_ROOT}/tools/VARIANT_ANNOTATION_V4_2014_03_25



