#Phyloseq code used to analyze "Microbial ecology of coastal northern Gulf of Mexico waters"

rm(list=ls())

require(ggplot2)
require(grid)
require(plyr)
require(lubridate)
require(phyloseq)
require(devtools)
require(knitr)
require(reshape2)
require(vegan)
require(data.table)
require(ggsignif)


#Set working directory
setwd("~") 


# Gets in the Nutrient Data
NUT<- read.csv("NUT.csv", colClasses = "character")



###############################################################################################

# IMPORTANT!!!
# Sample data and OTU table need to have matching sample names in the same order!!!!!!!!!!!!!!!
###############################################################################################

# Read in just the OTU for each sample
OTU<- read.csv("OTU.csv")


# Create matrix from data frame because phyloseq is prissy
VOL<- data.matrix(OTU, rownames.force = NA)


# Name the Columns based on OTU numer
rownames(VOL)<- paste0("OTU", 1:nrow(VOL))


# Read in Taxonomy Table
TAX<- read.csv("TAX.csv", colClasses = "character")


# Create another matrix from data frame and label first column as OTU number
LSU<- as.matrix(TAX, rownames.force = NA)
rownames(LSU)<- paste0("OTU", 1:nrow(LSU))


# Load phyloseq
library("phyloseq")

#Tell phyloseq what is what otu and tax tables need to be matrix
OTU = otu_table(VOL, taxa_are_rows = TRUE)
TAX = tax_table(LSU)



# Read the data into phyloseq
physeq = phyloseq(OTU, TAX)
physeq

# Make nutrient data the sample data
sampledata = sample_data(NUT)

# Add row names thats correspond to other phyloseq data
rownames(sampledata) = sample_names(physeq)

# MAKE SURE THE SAMPLE NAMES MATCH EACH OTHER
sampledata

# Nutrient data can be in data.frame
SAM = sample_data(sampledata)

# Read the nutrient data into phyloseq also
ALL = phyloseq(OTU, TAX, SAM)
ALL
OTUs<-cbind(otu_table(ALL), tax_table(ALL))
write.csv(OTUs, "ALL_otus.csv")
set.seed(200)

alpha_meas = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")
p<-plot_richness(ALL, x = "Filter", color="Filter", measures = alpha_meas)
p + geom_boxplot(data=p$data, aes(x=Filter, color=NULL)) + geom_signif(comparisons = list(c("Pre", "Ster")), 
                                                                       map_signif_level=TRUE) +theme_bw()


ALL_rarefy<- rarefy_even_depth(ALL, rngseed = 200)
ALL_rarefy
OTUs_rarefied<-cbind(tax_table(ALL_rarefy), otu_table(ALL_rarefy))
write.csv(OTUs, "ALL_rarefy_otus.csv")

#NMDS
iDist<- phyloseq::distance(ALL_rarefy, method= 'bray')
MDS2<- ordinate(ALL_rarefy, 'NMDS', distance = iDist) 
MDS2
p5<- plot_ordination(ALL_rarefy, MDS2, color="Water_Type", shape = "Filter", label = 'Site_Specific')
p6<-p5 + scale_color_manual(values = c("steelblue", "orange2", "red", "purple", "green")
) +
  geom_point(aes(color = Water_Type), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5) + theme_bw() 
p6
p6 +geom_segment(data= NMDS,
                 aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2))


#envfit for NUT data
NUT_nmds<-as.data.frame((sample_data(ALL_rarefy)))
NUT_nmds<- as.data.frame(t(t(NUT_nmds[,c(6,7,9,10)])))
NUT_nmds[]<-lapply(NUT_nmds, as.numeric)
NUT_scaled<-scale(NUT_nmds)
nmds.envfit<-envfit(MDS2, NUT_scaled, permu=999, na.rm = TRUE)
nmds.envfit

NMDS<-read.csv("NMDS.csv", row.names = 1 )

#Getting Only Sterviex Samples
OTUSTER<-subset_samples(ALL_rarefy, Filter == "Ster")
write.csv(otu_table(OTUSTER), "OTUSTER.csv")
OTUSTER_RA<-transform_sample_counts(OTUSTER, function(x) {x/sum(x)} )
filt<-21/24
OTUSTER_sal<-filter_taxa(OTUSTER_RA, function(x) sum(x>= 0.001) >= (0.9*length(x)), TRUE)
tax_table(OTUSTER_sal)
write.csv(otu_table(OTUSTER_RA), "OTUSTER_RA.csv")
tax.tab <- data.frame(tax_table(OTUSTER_RA))
ModifyTax <- function(x,ind){
  #   xth row in the dataframe
  #   ind taxonomy level to change
  if(is.na(x[ind])){
    nonNa <- which(!is.na(x[-ind]))
    maxNonNa <- max(nonNa)
    x[ind] <- paste(x[maxNonNa],".",x[ind])
  }else{x[ind] <- x[ind]}
}

#replace the genus taxonomy with the highest known taxonomy
tax_table(OTUSTER_RA)[,6] <- apply(tax.tab,1,ModifyTax,ind=6)
write.csv(tax_table(OTUSTER_RA), "OTUSTER_Taxmod.csv")

##Repeat for Prefilter
OTUPre<-subset_samples(ALL_rarefy, Filter == "Pre")
write.csv(otu_table(OTUPre), "OTUPre.csv")
OTUPre_RA<-transform_sample_counts(OTUPre, function(x) {x/sum(x)} )
write.csv(otu_table(OTUPre_RA), "OTUPre_RA.csv")
OTUPre_sal<-filter_taxa(OTUPre_RA, function(x) sum(x>= 0.001) >= (0.9*length(x)), TRUE)
tax_table(OTUPre_sal)
tax.tab <- data.frame(tax_table(OTUPre_RA))

ModifyTax <- function(x,ind){
  #   xth row in the dataframe
  #   ind taxonomy level to change
  if(is.na(x[ind])){
    nonNa <- which(!is.na(x[-ind]))
    maxNonNa <- max(nonNa)
    x[ind] <- paste(x[maxNonNa],".",x[ind])
  }else{x[ind] <- x[ind]}
}

#   replace the genus taxonomy with the highest known taxonomy
tax_table(OTUPre_RA)[,6] <- apply(tax.tab,1,ModifyTax,ind=6)
write.csv(tax_table(OTUPre_RA), "OTUPre_Taxmod.csv")


#plot Ster Top50 (median) data
data<-read.csv("OTUSTER_RA.csv")
melted<-melt(data)
ggplot(melted, aes(variable, value)) +geom_boxplot(aes(), outlier.shape=NA) + geom_point(aes(color=Salinity), size =1) +
   scale_color_manual(values = c("steelblue", "orange2", "red", "purple", "green")
  ) +
  geom_point(aes(color = Salinity), alpha = 0.7, size = 1.5) +
  geom_point(colour = "grey90", size = 0.25) + theme_bw() + theme(axis.text.x = element_text(angle = 90))

#plot Pre Top50 (median) data
data<-read.csv("OTUPre_RA_top50.csv")
melted<-melt(data, value.name = "RA")
ggplot(melted, aes(variable, RA)) +geom_boxplot(aes(), outlier.shape=NA) + geom_point(aes(color=Salinity)) +
  scale_color_manual(values = c("steelblue", "orange2", "red", "purple", "green")
  ) +
  geom_point(aes(color = Salinity), alpha = 0.7, size = 1.5) +
  geom_point(colour = "grey90", size = 0.25) + theme_bw() + theme(axis.text.x = element_text(angle = 90))

# Site Map
require(ggspatial)
require(ggplot2)
require(sf)
require(tigris)
require(dplyr)
require(nhdplusTools)


LAmap<-read.csv("LatLong.csv", header=T)
coast<-st_read("../Downloads/louisiana_coastline/")
rivers<-st_read("../Downloads/USA_Rivers_and_Streams-shp/")
louisiana_r <- rivers %>%
    filter(State == "LA")
louisiana <- states(cb = TRUE) %>%
  filter(STUSPS == "LA")                      
ggplot() + geom_sf(data=louisiana)  + theme_bw()  +   annotation_scale() +  annotation_north_arrow(location = "tr", which_north = "true",) + 
                        geom_sf(data=basin, alpha=0, linetype = "dashed", size = 1.5)
ggplot() + geom_sf(data=louisiana)  + theme_bw()  +   annotation_scale() +  
                        annotation_north_arrow(location = "tr", which_north = "true",width = 0.2) + geom_sf(data=louisiana_r, color="light blue", alpha=0.75) + 
                        geom_sf(data=basin, alpha=0, linetype = "dashed", size = 1.5) + coord_sf(xlim = c(min_lon, max_lon), ylim = c(min_lat, max_lat)) +
                        geom_jitter(data=LAmap, aes(Longitude, Latitude, color=Descript, size=Salinity),alpha=0.5, stroke=0.5, width=0.2) + 
                        scale_color_manual(values = c("steelblue","red", "orange2" ))+scale_size_continuous(range = c(3, 6))
