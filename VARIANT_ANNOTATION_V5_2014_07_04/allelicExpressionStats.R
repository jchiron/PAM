############
# Routines #
############

# Merge two DF
mergeDF=function(DF1,DF2) {
  DF1=merge(DF1,DF2,by="row.names",sort=F)
  rownames(DF1)=DF1[,1]
  DF1[,1]=NULL
  return(DF1)
}

# Add one column
addCol=function(DF,vector,colName) {
  DF=cbind(DF,vector)
  colnames(DF)[ncol(DF)]=colName
  return(DF)
}

# Show progression using index
progressInfo=function(index,max) {
  if (index%%5000==0) {
    print(paste("   ",index,"/",max))
  }
}

# Get index and name for each sample
getSampleInfo=function(colNames,firstColSuffix) {
  return(gsub(firstColSuffix,"",colNames[grep(firstColSuffix,colNames)]))
}

# Determine genotype of sample using statistics
genotypeCall=function(DF,sampleName,DP4TotCol,AltTotCol,
                      expectedMinProb,expectedMaxProb,confInt) {
  est=c()
  intMin=c()
  intMax=c()
  genotype=c()
  nRow=nrow(DF)
  print("  Determining genotypes")
  # For each line of DF
  for (i in 1:nRow) {
    # If record exists at this location, process it
    if (DF[i,DP4TotCol]>0) {
      # Get the number of trials
      binomTot=DF[i,DP4TotCol]
      # Get the expected minimum number of successes
      binomMin=round(binomTot*expectedMinProb)
      # Compute the binomial
      resMin=binom.test(x=binomMin,n=binomTot,p=expectedMinProb,
                     alternative="less",conf.level=confInt)
      # Store confidence intervals
      intMin=c(intMin,resMin$conf.int[2])
      intMax=c(intMax,1-resMin$conf.int[2])
      # Binomial is symetric so do not compute it in highiest side
      #binomMax=round(binomTot*expectedMaxProb)
      #resMax=binom.test(x=binomMax,n=binomTot,p=expectedMaxProb,
      #               alternative="greater",conf.level=confInt)
      #intMax=c(intMax,resMax$conf.int[1])
      AltRatio=DF[i,AltTotCol]/DF[i,DP4TotCol]
      est=c(est,AltRatio)
      if (AltRatio<resMin$conf.int[2]) { genotype=c(genotype,"AA") }
      else if (AltRatio>1-resMin$conf.int[2]) { genotype=c(genotype,"BB") }
      else { genotype=c(genotype,"AB") }
      
    } else {
      # Case there are no records, put NA
      est=c(est,NA)
      intMin=c(intMin,NA)
      intMax=c(intMax,NA)
      genotype=c(genotype,NA)
    }
    # Show progression
    #progressInfo(i,nRow)
  }
  # Init a DF to store output
  out=data.frame(est=round(est,2),intMin=round(intMin,2),
                 intMax=round(intMax,2),genotype=genotype)
  # Modify colnames
  colnames(out)=c(paste(sampleName,"_ALT_RATIO",sep=""),
                  paste(sampleName,"_CONF_INT_MIN",sep=""),
                  paste(sampleName,"_CONF_INT_MAX",sep=""),
                  paste(sampleName,"_GENOTYPE_(THRESH=",expectedMinProb,")",sep=""))
  rownames(out)=rownames(DF)
  return(out)
}

# Detect strand biais
strandBiais=function(DF,sampleName,RefFwCol,AltFwCol,DP4TotCol,pvalThresh,
                     expectedProb,confInt,adjustMethod,altHyp) {
  est=c()
  pvals=c()
  print("  Strand biais analysis")
  nRow=nrow(DF)
  # For each row in DF
  for (i in 1:nRow) {
    # Get the forward total
    totFW=DF[i,RefFwCol]+DF[i,AltFwCol]
    # Get the total
    tot=DF[i,DP4TotCol]
    nRow=nrow(DF)
    # If there is a record at this location
    if (tot>0) {
      # Compute the binomial
      res=binom.test(x=totFW,n=tot,p=expectedProb,
                     alternative=altHyp,conf.level=confInt)
      # Store metrics
      est=c(est,res$estimate)
      pvals=c(pvals,res$p.value)
    } else {
      # No record, put NA
      est=c(est,NA)
      pvals=c(pvals,NA)
    }
    # Show progression
    #progressInfo(i,nRow)
  }
  # Adjust p-values
  pvals=p.adjust(p=pvals,method=adjustMethod)
  strandbiais=c()
  # For each p-value
  for (i in 1:length(pvals)) {
    # If NA, set strand biais to NA
    if (is.na(pvals[i])) {
      strandbiais=c(strandbiais,NA)
    } else if (pvals[i]>pvalThresh) {
      # Else if pval is not significant -> no strand biais
      strandbiais=c(strandbiais,"No")
    } else {
      # Else -> strand biais
      strandbiais=c(strandbiais,"Yes")
    }
  }
  # Init a DF for output
  out=data.frame(est=round(est,2),strandbiais=strandbiais)
  # Modify colnames
  colnames(out)=c(paste(sampleName,"_FW_RATIO",sep=""),
                  paste(sampleName,"_STRAND_BIAIS",sep=""))
  rownames(out)=rownames(DF)
  return(out)
}

# Analyse ref and alt bases to associate mutation to a defined model
baseModification=function(DF) {
  # Set all mutations to undefined
  DF["BASE_MODIFICATION"]="UNDEFINED"
  DF["BASE_SUBSTITUTION"]="UNDEFINED"
  # Set indels
  indelIndices=which(DF[,4]=="-" | DF[,5]=="-")
  DF[indelIndices,"BASE_MODIFICATION"]="INDEL"
  DF[indelIndices,"BASE_SUBSTITUTION"]="INDEL"
  # Set each combination
  comb1Indices=which((DF[,4]=="C" & DF[,5]=="A") | (DF[,4]=="G" & DF[,5]=="T"))
  DF[comb1Indices,"BASE_MODIFICATION"]="C:G>A:T"
  DF[comb1Indices,"BASE_SUBSTITUTION"]="Transversion"
  comb2Indices=which((DF[,4]=="C" & DF[,5]=="G") | (DF[,4]=="G" & DF[,5]=="C"))
  DF[comb2Indices,"BASE_MODIFICATION"]="C:G>G:C"
  DF[comb2Indices,"BASE_SUBSTITUTION"]="Transversion"
  comb3Indices=which((DF[,4]=="C" & DF[,5]=="T") | (DF[,4]=="G" & DF[,5]=="A"))
  DF[comb3Indices,"BASE_MODIFICATION"]="C:G>T:A"
  DF[comb3Indices,"BASE_SUBSTITUTION"]="Transition"
  comb4Indices=which((DF[,4]=="T" & DF[,5]=="A") | (DF[,4]=="A" & DF[,5]=="T"))
  DF[comb4Indices,"BASE_MODIFICATION"]="T:A>A:T"
  DF[comb4Indices,"BASE_SUBSTITUTION"]="Transversion"
  comb5Indices=which((DF[,4]=="T" & DF[,5]=="C") | (DF[,4]=="A" & DF[,5]=="G"))
  DF[comb5Indices,"BASE_MODIFICATION"]="T:A>C:G"
  DF[comb5Indices,"BASE_SUBSTITUTION"]="Transition"
  comb6Indices=which((DF[,4]=="T" & DF[,5]=="G") | (DF[,4]=="A" & DF[,5]=="C"))
  DF[comb6Indices,"BASE_MODIFICATION"]="T:A>G:C"
  DF[comb6Indices,"BASE_SUBSTITUTION"]="Transversion"
  # Test is there are some undefined mutations left
  undefined=grep("UNDEFINED",DF[,"BASE_MODIFICATION"])
  for (i in undefined) {
    # If ref and alt base are weird (GAAAA -> GAAA), set INDEL
    if (nchar(as.character(DF[i,4]))>1 & nchar(as.character(DF[i,5]))>1) {
      DF[i,"BASE_MODIFICATION"]="INDEL"
      DF[i,"BASE_SUBSTITUTION"]="INDEL"
    }
  }
  # If there are still some weird base mutations, print them
  undefined=grep("UNDEFINED",DF[,"BASE_MODIFICATION"])
  if (length(undefined)>0) {
    print(paste("Warning: undefined base modifications"))
    print(DF[undefined,])
  }
  return(DF)
}

# Round LJB2 scores
roundLJB2=function(DF) {
  LJB2=c("SIFT_SCORE","POLYPHEN2_HDIV_SCORE","POLYPHEN2_HVAR_SCORE","LRT_SCORE",
         "MUTATIONTASTER_SCORE","MUTATIONASSESSOR_SCORE","FATHMM_SCORE",
         "GERP++_SCORE","PHYLOP_SCORE","SIPHY_SCORE")
  for (i in LJB2) {
    DF[,i]=round(DF[,i],2)
  }
  return(DF)
}

########
# Main #
########

# First paramaters
args=commandArgs(trailingOnly = TRUE)
if (length(args)!=2) {
  print("You must provide 2 parameters : input and output file.")
  quit("yes")
}
print("Loading data and initiazing parameters")
input=args[1]
output=args[2]
## CHANGE CLucchesi  ---> add sep="\t"
#data=read.table(input,header=T,comment.char="",check.names=F)
# data=read.table(file=input, sep="\t", header=T,comment.char="",check.names=F)
data=read.table(file=input, sep="\t", header= TRUE, comment.char="", check.names=FALSE, colClasses=c("REF"="character","ALT"="character"), as.is = TRUE, quote="\"")
rownames(data)=paste(data[,1],data[,2],data[,3],data[,4],data[,5],sep=":")

# Binomial parameters
pvalThresh=0.05
expectedProb=0.5
expectedMinProb=0.05
expectedMaxProb=0.95
confInt=0.95
adjustMethod="BH"
altHyp="two.sided"

# Input suffixes of interest
firstColSuffix="_DP4_REF_FW"
DP4TotSuffix="_DP4_TOT"
AltTotSuffix="_DP4_ALT"
RefFwSuffix="_DP4_REF_FW"
RefBwSuffix="_DP4_REF_BW"
AltFwSuffix="_DP4_ALT_FW"
AltBwSuffix="_DP4_ALT_BW"
firstInfoColSuffix="GENE"

# Get sample names
sampleNames=getSampleInfo(firstColSuffix=firstColSuffix,colNames=colnames(data))

# Init out DF with some data
out=data[,1:5]

# For each sample
for (sampleName in sampleNames) {
  print(paste("Processing:",sampleName))
  # Create names for columns of Interest
  DP4TotCol=paste(sampleName,DP4TotSuffix,sep="")
  AltTotCol=paste(sampleName,AltTotSuffix,sep="")
  RefFwCol=paste(sampleName,RefFwSuffix,sep="")
  RefBwCol=paste(sampleName,RefBwSuffix,sep="")
  AltFwCol=paste(sampleName,AltFwSuffix,sep="")
  AltBwCol=paste(sampleName,AltBwSuffix,sep="")

  # Add some DP4 metrics
  out=addCol(DF=out,vector=data[,paste(sampleName,AltTotSuffix,sep="")],
             colName=paste(sampleName,AltTotSuffix,sep=""))
  out=addCol(DF=out,vector=data[,paste(sampleName,DP4TotSuffix,sep="")],
             colName=paste(sampleName,DP4TotSuffix,sep=""))
  
  # Compute genotype call
  res=genotypeCall(DF=data,sampleName=sampleName,
                   DP4TotCol=DP4TotCol,AltTotCol=AltTotCol,
                   expectedMinProb=expectedMinProb,
                   expectedMaxProb=expectedMaxProb,confInt=confInt)
  # Add result in output
  out=mergeDF(DF1=out,DF2=res)
  
  # Compute strand biais
  res=strandBiais(DF=data,sampleName=sampleName,
                  RefFwCol=RefFwCol,
                  AltFwCol=AltFwCol,
                  DP4TotCol=DP4TotCol,
                  pvalThresh=pvalThresh,expectedProb=expectedProb,
                  confInt=confInt,adjustMethod=adjustMethod,
                  altHyp=altHyp)
  # Add result in output
  out=mergeDF(DF1=out,DF2=res)
}
print("Formating output data")
# Concatenate genotypes and create a new column
genotypes=c()
genotypeIndices=grep("_GENOTYPE_",colnames(out))
for (i in seq(1,length(genotypeIndices),1)) {
  if (i==1)
    genotypes=out[,genotypeIndices[i]]
  else
    genotypes=paste(genotypes,out[,genotypeIndices[i]],sep=":")
} 
out=addCol(DF=out,vector=genotypes,colName="ALL_GENOTYPES")
# Set indices to firstInfoCol (before 1000G) and lastInfoCol (after1000G)
firstInfoColIndex=grep(paste("^",firstInfoColSuffix,"$",sep=""),colnames(data))
transcriptsCol=grep("^TRANSCRIPTS$",colnames(data))
cosmicCol=grep("^COSMIC$",colnames(data))
out=mergeDF(DF1=out,DF2=data[,firstInfoColIndex:transcriptsCol])
for (i in c("ALL","EUR","AMR","ASN","AFR")) {
  AAprob=paste("1000G_",i,"_AA_PROB",sep="")
  ABprob=paste("1000G_",i,"_AB_PROB",sep="")
  BBprob=paste("1000G_",i,"_BB_PROB",sep="")
  g1000Name=paste("1000G_",i,"_AF",sep="")
  g1000Probs=paste("1000G_",i,"_PROBS_AA:AB:BB",sep="")
  out[,g1000Name]=data[,substr(g1000Name,1,nchar(g1000Name)-3)]
  out[,g1000Probs]=paste(round((1-out[,g1000Name])^2,2),":",
    round(2*out[,g1000Name]*(1-out[,g1000Name]),2),":",
    round(out[,g1000Name]^2,2),sep="")
}

lastInfoColIndex=ncol(data)
# Add info
out=mergeDF(DF1=out,DF2=data[,cosmicCol:lastInfoColIndex])
# Add base mutation information : transition/transversion/indel
out=baseModification(out)
out=roundLJB2(out)
# Sort and output final DF
out=out[with(out,order(CHROM,START,STOP)),]
write.table(x=out,file=output,sep="\t",col.names=colnames(out),row.names=F)

#source("Plots/Substitutions.R")
#source("Plots/Mutations.R")
#source("Plots/SomePlots.R")

#pdf("SomeMetrics_All.pdf")
#plotMetrics(out)
#dev.off()

#pdf("SomeMetrics_253mut.pdf")
#plotMetrics(interest)
#dev.off()

#pdf("SomeMetrics_1231mut.pdf")
#plotMetrics(interest2)
#dev.off()

#pdf("Filtering_steps.pdf")
#par(mar=c(5.1,4.1,4.1,5.1))
#interest=plotFilteringSteps(out,sampleNames,"T")
#interest2=plotFilteringSteps(out,sampleNames,"N")
#par(mar=c(5.1,4.1,4.1,2.1))
#dev.off()
