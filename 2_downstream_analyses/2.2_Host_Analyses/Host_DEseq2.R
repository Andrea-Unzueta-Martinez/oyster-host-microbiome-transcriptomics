###
###### Differential gene expression with DEseq2 #####
###

library(DESeq2)

#very important to have metadata row names match column names of counts matrix exactly 
metadata # remove the two mock community samples 
#remove mocks communities
metadata_nomocks = metadata %>% filter(Sample_ID != c("mock1","mock2"))
#set sample names as row names
rownames(metadata_nomocks) = metadata_nomocks$Sample_ID
#add column with interaction term
metadata_nomocks$Interaction = paste(metadata_nomocks$Treatment,metadata_nomocks$Timepoint, sep = "_")

rownames(metadata_nomocks)
colnames(txi.tx$counts)

#order metadata row names in the order columns show up in the counts table
sampleTable = metadata_nomocks[match(names(files), rownames(metadata_nomocks)),]

#use salmon counts imported with Tximport
dds2 = DESeqDataSetFromTximport(txi.tx, sampleTable, ~ Interaction)

###quality filter dds2 object 

#remove control samples from dds2
colnames(dds2) #figure out what col number each of the neagatives is
dds2_filt = dds2[,-c(4,8,12,16)] #remove negaties
colnames(dds2_filt)
colData(dds2_filt)$Interaction

#re-define levels, to exclude "procedural_control_NA"
dds2_filt$Interaction = factor(dds2_filt$Interaction, levels = c("control_6 hour","control_9 hour","control_12 hour","exposure_6 hour","exposure_9 hour","exposure_12 hour"))

#remove low abundance transcripts that have < 5 reads
dds2_filt # 66451 features AKA isophorms
dds2_filt = dds2_filt[rowSums(counts(dds2_filt)) > 5 ] #55005 isophorms to go into DEseq2

colSums(counts(dds2_filt))
rowSums(counts(dds2_filt)) 

#names of genes that went into deseq analysis
transcrips_ids_deseq2_in = rownames(dds2_filt)

#note which levels are reference
dds2_filt$Treatment = relevel(factor(dds2_filt$Treatment), ref = "control")
dds2_filt$Timepoint = relevel(factor(dds2_filt$Timepoint), ref = "9 hour")

#now try differential abundance test
dds2_test = DESeq2::DESeq(dds2_filt)

colData(dds2_filt)

#to extract normalized counts from the deseq object:
counts(dds2_test,normalized = TRUE)

#extract results
#res2 = DESeq2::results(dds2_test, name ="Interaction_exposure_9.hour_vs_control_9.hour")
resultsNames(dds2_test)
res2 = DESeq2::results(dds2_test, contrast = c("Interaction","exposure_9 hour","control_9 hour"))
#sort results by p value
res2 = res2[order(res2$padj),]
#examine results summary
summary(res2)
plotMA(res2)
#p = 0.05
res05 <- results(dds1, alpha=0.05)
summary(res05)

#how many p-values were less than 0.05?
sum(res2$padj < 0.05, na.rm=TRUE) #791
#how many adjusted p-values were less than 0.1?
sum(res2$padj < 0.1, na.rm=TRUE) #952

#extract comparisons of interest:
#control vs exposure at each time point independently
control_9hour_vs_exposure_9hour = results(dds2_test, contrast = c("Interaction","exposure_9 hour","control_9 hour"))
control_9hour_vs_exposure_9hour_resLFC = lfcShrink(dds2_test,contrast = c("Interaction","exposure_9 hour","control_9 hour"), type="normal")

control_6hour_vs_exposure_6hour = results(dds2_test, contrast = c("Interaction","control_6 hour","exposure_6 hour"))
control_6hour_vs_exposure_6hour_resLFC = lfcShrink(dds2_test,contrast = c("Interaction","exposure_6 hour","control_6 hour"), type="normal")

control_12hour_vs_exposure_12hour = results(dds2_test, contrast = c("Interaction","control_12 hour","exposure_12 hour"))
control_12hour_vs_exposure_12hour_resLFC = lfcShrink(dds2_test,contrast = c("Interaction","exposure_12 hour","control_12 hour"), type="normal")


