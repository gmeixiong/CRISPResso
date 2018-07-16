##Pulling the names of all fastq files in the directed folder, as well as reference file excel sheet
######CHANGE THIS VARIABLE##########
datafolder <- "Sortedfastq"
#######DID YOU CHANGE THE VARIABLE?#######

######CHANGE THIS VARIABLE##########
referenceseq <- "CRISPResso pipeline 96wp1sorted sample sheet.xlsx"
#######DID YOU CHANGE THE VARIABLE?#######

##Where are the trimmomatic PE adapters?###
adapterloc <- "/home/ubuntu/genome/TruSeq3-PE-2.fa"
#######################################

##What is the read lenght used in the experiment?###
readleng <- "250"
#######################################

##What is average fragment length (before adapters)?###
fragleng <- "300"
##################################################

##What is the standard deviation of the fragment legnths?###
###If unknown, assume 10% of the average for now#######
stdevfrag <- "45"
#######################################


#Incorporation of libraries for script
library(S4Vectors)
library(readxl)
library(plyr)
library(stringr)
library(gdata)
library(knitr)
library(tidyverse)
library(Rsamtools)
library(rtracklayer)
library(Biostrings)
library(readxl)

foldertxt <- paste0(datafolder,"/")
all_fq_fnames <- dir(datafolder, "R._001.fastq.gz$", recursive=TRUE, full=TRUE)

##Creating a dataframe with columns for all fastqs, forward fastqs, and reverse fastqs
fq_df <- data_frame(all_fq_fnames)
##Naming fq_fwd and fq_rev columns for incorporation of both forward and reverse reads
fq_df$fq_fwd <- NA
fq_df$fq_rev <- NA


##For loop iterating through all of the fastq files, pulling out any reads that are designated as "forward" with the _R1_ designation
##If a read is found, it is placed in the second column of fq_df - all forward reads will be placed here
for(i in 1:nrow(fq_df)){
  for(j in 1:nrow(fq_df)){
    if(grepl("_R1_001.fastq.gz$",fq_df[j,1]) == TRUE){
      fq_df[j,2] <- fq_df[j,1]
    }
  } 
}
##For loop iterating through all of the fastq files pulling out any reads that are designated as "reverse" through subbing R2 with R1, and matching every other part of the fastq filename
##If a mate pair is found, it is placed adjacent to the forward read(column2) in the third column.
for(i in 1:nrow(fq_df)){
  pe_check <- fq_df[i,1]
  pe_check <- gsub("R1", "NA", pe_check)
  pe_check <- gsub("R2", "R1", pe_check)
  for(j in 1:nrow(fq_df)){
    if(grepl(pe_check,fq_df[j,1]) == TRUE){
      fq_df[j,3] <- fq_df[i,1]
    }
  } 
}

#Selecting for only the forward and reverse reads, discarding the un-sorted list. Removes any cases that contain an NA
fq_df <- subset(fq_df, select=c(2,3))
fq_df <- fq_df[complete.cases(fq_df), ]

###############################################################################################
###############################################################################################
##Creating reference amplicons from imported excel file
refimport <- read_excel(paste0(datafolder, "/", referenceseq))

fq_df$reference <- NA
fq_df$HDR <- NA
fq_df$guideseq <- NA

for(i in 1:nrow(refimport)){
  iterate <- as.character(refimport$Sample_ID[i])
  for(j in 1:nrow(fq_df)){
    if(grepl(iterate,fq_df[j,1]) == TRUE){
      fq_df$reference[j] <- refimport$`WT amplicon`[i]
      fq_df$HDR[j] <- refimport$`HDR amplicon`[i]
      fq_df$guideseq[j] <- refimport$`Guide Sequence`[i]
    }
  } 
}


###############################################################################################
###############################################################################################
###############################################################################################
#Selecting fq_fnames only for fwd reads to re-name the files after merging
fq_namesonly <- c(NA,NA)
for(i in 1:nrow(fq_df)){
  fq_namesonly[i] <- gsub(foldertxt, "", fq_df$fq_fwd[i])
  fq_namesonly[i] <- sub("^(.*)[_].*", "\\1", fq_namesonly[i])
}


##Create Folder Name for output data from fastq name
dir.create(paste0(datafolder,"_CRISPResso_Out"))
dataout <- paste0(datafolder,"_CRISPResso_Out")

##Merging PE reads using FLASH prior to adapter trimming and mapping
##flash -r 151 -f 227 -s 21 -o "NAME" -d "DIRECTORY" -z input1 input2

for(i in 1:nrow(fq_df)){
  mergecmd <- paste0("flash -r ", readleng, " -f ", fragleng, " -s ", stdevfrag, 
                     " -o ", fq_namesonly[i],
                     " -z ",
                     fq_df[i,1], " ", fq_df[i,2])
  print(mergecmd)
  message(mergecmd, "\n"); system(mergecmd)
}

##Creating variable names for quality trimmed reads. Qual_P is the named designation 
fq_df$combined <- gsub("_001.fastq.gz$",".extendedFrags.fastq.gz",basename(fq_df$fq_fwd))
fq_df$qual_P <- gsub("_R1_001.fastq.gz$","_R1_qual_P.fastq.gz",basename(fq_df$fq_fwd))

##Adding the folder designation for quality trimmed reads
##Creating folder name for filter
dir.create(paste0(datafolder,"_filter"))
fq_df$qual_P <- paste0(datafolder,"_filter/",fq_df$qual_P)


###Removing adapters from reads, then passing them through quality filter, saving in new folder
##See trimmomatic for settings. A sliding window searches along for quality score at a minimum of 
##20 PHRED33 score, also trimming any adapter sequences.  
###############SET THESE VARIABLES###########################
leadingqual <- 3
trailingqual <- 3
slidewindsiz <- 4
slidewindqual <- 20
minlength <- 50
############################################################


for(i in 1:nrow(fq_df)){
  qualcmd <- paste0("trimmomatic SE ", fq_df$combined[i], " ", fq_df$qual_P[i], " ",
                    "ILLUMINACLIP:",adapterloc,":2:30:10", " ",
                    "LEADING:", as.character(leadingqual), " ", "TRAILING:", as.character(trailingqual), " ", 
                    "SLIDINGWINDOW:", as.character(slidewindsiz), ":", as.character(slidewindqual), " ", 
                    "MINLEN:", as.character(minlength))
  print(qualcmd)
  message(qualcmd, "\n"); system(qualcmd)
}

##Clearing intermediate files
remcmd <- paste0("find . -name '*extendedFrags.fastq.gz' -type f -delete")
print(remcmd)
message(remcmd, "\n"); system(remcmd)
remcmd <- paste0("find . -name '*notCombined_[12].fastq.gz' -type f -delete")
print(remcmd)
message(remcmd, "\n"); system(remcmd)
remcmd <- paste0("find . -name '*.hist' -type f -delete")
print(remcmd)
message(remcmd, "\n"); system(remcmd)
remcmd <- paste0("find . -name '*.histogram' -type f -delete")
print(remcmd)
message(remcmd, "\n"); system(remcmd)

##Maping sorting and indexing bam files from sucessfully paired reads using bwa mem
##Only mapping paired reads, fwd, rev, and then placing them in the bm_fname folder location.

for(i in 1:nrow(fq_df)) {
   cmd <- paste0("CRISPResso -r1 ", fq_df$qual_P[i], " -a ", fq_df$reference[i], " -e ", fq_df$HDR[i], " -g ", 
                 fq_df$guideseq[i], " -o ", dataout)
   print(cmd)
   message(cmd, "\n"); system(cmd)
}
