#
###### Get Salmon Gene counts on R to build host counts matrix ######
#

#use package tximportData to make counts table with all samples 
BiocManager::install("tximportData") 
library(tximportData)

#identify the directory where the Salmon output data are locate
dir = ("/Users/andreaunzueta/Dropbox/NSF_PRFB_GirguisLab/Bioinformatics/Metatranscriptomics_analysis/OysterTranscripts/Salmon_counts/Salmon_output2")
list.files(dir)
#define samples
samples <- read.table(file.path(dir, "samples.txt"), header = TRUE) # this text file contains sample names
samples$sample

#path to each quantification file
files = file.path(dir, samples$sample, "quant.sf")
names(files) = c("e_6h_62A","e_12h_64A","c_9hr_43","neg1",
                 "e_6h_27A","e_12h_29A","c_9hr_48","neg2",
                 "e_6h_32A","e_12h_34A","c_9hr_53","neg3",
                 "e_6h_37","e_12h_39A","c_9hr_58", "neg4",
                 "e_9h_63B","c_6hr_42","c_12hr_44",
                 "e_9h_33B","c_6hr_47","c_12hr_49",
                 "e_9h_28B","c_6hr_52","c_12hr_54A","e_9h_38B",
                 "c_6hr_57","c_12hr_59A")

metadata$Sample_ID 
all(file.exists(files))

#upload salmon output files (containing counts for each sample)
txi.tx <- tximport(files, type = "salmon", txOut = TRUE) #the abundace matrix is already normalized to TPM (this was part of Salmon specifications)

dim(txi.tx$abundance)

#remove rows that have a sum of 0
# Identify rows with a sum of zero
non_zero_rows <- rowSums(txi.tx$abundance) > 0
table(non_zero_rows) #59212 transcripts are non zero across samples AKA rows, 7239 transcripts have row sums of zero 

# Subset the matrix to keep only rows with non-zero sums!!
txi.tx$abundance <- txi.tx$abundance[non_zero_rows, ]

# View the updated matrix to confirm
dim(txi.tx$abundance) # Check dimensions to confirm rows were removed, dimentions should be 59212 genes and 28 samples

#save counts matrix to use for downstream analyses and upload on github
write.csv(x = txi.tx$counts, file = "oyster_counts_matrix.csv", row.names = TRUE)


#
######  Load Annotations of oytser transcripts #####
#

#NCBI downloaded annotations
Oyster_transcript_annotations_NCBI = data.frame(rtracklayer::import("/Users/andreaunzueta/Dropbox/NSF_PRFB_GirguisLab/Bioinformatics/Metatranscriptomics_analysis/OysterTranscripts/NCBI\Annotations/cvirginica_clean.gtf"))

#emapper annotatations (created in the emmaper annotations step fo the raw reads to raw counts pipeline)
Oyster_transcript_annotations_emapper = read.table("/Users/andreaunzueta/Dropbox/NSF_PRFB_GirguisLab/Bioinformatics/Metatranscriptomics_analysis/OysterTranscripts/eggNOG_annotations/oyster_transcriptome_anotations.emapper.annotations.tsv",sep = '\t',header = TRUE)

