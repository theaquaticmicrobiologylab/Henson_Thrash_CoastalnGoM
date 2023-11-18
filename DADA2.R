
##Code for analyzing sequencing data from Microbial ecology of coastal northern Gulf of Mexico water
#Installing DADA2 via https://github.com/benjjneb/dada2
#Citation :https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4927377/#:~:text=The%20Divisive%20Amplicon%20Denoising%20Algorithm,errors%20without%20constructing%20OTUs5.&text=Sample%20composition%20is%20inferred%20by,applicable%20to%20any%20genetic%20locus.
# Will only need to do once
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", force= TRUE)

##Workflow
library(Matrix)
library(dada2); packageVersion("dada2")


#should be path to files
path<- ""

#Pulling out fastq file names
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#Inspecting the quality of your F and R reads 
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

#Filtering  and trimming reads based on the quality from above. 
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
FiltReads <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                    maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                    compress=TRUE, multithread=TRUE)
#Check that you didn't lose large chuncks of reads (sanity check)
head(FiltReads)

#Calculate error rates from reads
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

#Unique Sequences/sequence variants
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errF, multithread=TRUE)

#Merge reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

#Creating a sequence table and removing reads < 250 bp and >256 bp
seqtab <- makeSequenceTable(mergers)
table(nchar(getSequences(seqtab)))
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:256]
table(nchar(getSequences(seqtab2)))

#removing Chimeras 
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)

#Tracking reads through the process.
getN <- function(x) sum(getUniques(x))
track <- cbind(FiltReads, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#Taxonomy (will need to have silva datasbases in your directory)
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr_v138_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "silva_species_assignment_v138.fa.gz")

#Wrote out tax and seq tables to dataframes and saved as .csv files.
seqTab_output<-as.data.frame(seqtab.nochim)
seqTab_output<-t(seqTab_output)
taxaTab_output<-as.data.frame(taxa)
write.csv(taxaTab_output, "coastalAMP_DADA2_taxtab.csv")
write.csv(seqTab_output, "coastalAMP_DAD2_counttab.csv")
#Manually removed "Chloroplast, Mitochondria, Eukaryotes, and NAs"




