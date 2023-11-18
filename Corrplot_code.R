# Code was developed from Dr. Pat Schloss, 
# https://riffomonas.org/code_club/2021-06-03-filter-by-significance
# https://www.youtube.com/watch?v=lhst1oc9mKQ&list=PLmNrK_nkqBpIIRdQTS2aOs5OD7vVMKWAi&index=36
# https://www.youtube.com/watch?v=-DfCitMqdUc&list=PLmNrK_nkqBpIIRdQTS2aOs5OD7vVMKWAi&index=37
#DOI: https://doi.org/10.1128/mra.01310-22


rm(list = ls())

library(tidyverse)
library(readxl)
library(ggtext)
library(glue)
library(broom)

cbbPalette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

metadata <- read_excel(path="/Path/to/file/Metadata.csv") %>%
  rename(sample_id = Sample)

##This OTU file should be Sample Names in column 1 and OTU# in row 1. ##
otu_counts <- read_csv("Path/to/file/OTU.csv") %>%
  select(Group, starts_with("OTU")) %>%
  rename(sample_id = Group) %>%
  pivot_longer(-sample_id, names_to="otu", values_to = "count")

nseqs_per_sample <- otu_counts %>%
  group_by(sample_id) %>%
  summarize(N = sum(count), .groups="drop") %>%
  count(N) %>%
  pull(N)

nseqs_per_sample 
stopifnot(length(nseqs_per_sample) == 1)

lod <- 100* 1/nseqs_per_sample

taxonomy <- read_tsv("Path/to/file/Tax.txt") %>%
  select("OTU", "Taxonomy") %>%
  rename_all(tolower) %>%
  mutate(taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)", ""),
         taxonomy = str_replace(taxonomy, ";$", "")) %>%
  separate(taxonomy,
           into=c("kingdom", "phylum", "class", "order", "family", "genus"),
           sep=";") %>%
  mutate(pretty_otu = str_replace(string=otu,
                                  pattern="tu0*",
                                  replacement = "TU "),
         order = str_replace(string=order,
                             pattern="(.*)",
                             replacement="*\\1*"),
         order = str_replace(string=order,
                             pattern="\\*(.*)_unclassified\\*",
                             replacement="Unclassified *\\1*"),
         taxon = glue("{order} ({pretty_otu})"),
         taxon = str_replace_all(taxon, "_", " ")) %>%
  select(otu, taxon,phylum)

#Moving to working on OTU data 
otu_rel_abund <- inner_join(metadata, otu_counts, by="sample_id") %>%
  inner_join(., taxonomy, by="otu") %>%
  group_by(sample_id) %>%
  mutate(rel_abund = 100*count / sum(count)) %>%
  ungroup() %>%
  select(-count) 

##For Ster (0.2-2.7)
ster_otu_rel_abund<-  otu_rel_abund %>%
  filter(str_detect(sample_id, "Ster"))
nicelooking<-ster_otu_rel_abund %>%
  group_by(taxon) %>%
  summarize(mean=mean(rel_abund),median=median(rel_abund),freq = sum(rel_abund != 0) / sum(rel_abund >= 0), maxRA=max(rel_abund),.groups="drop")
ster_taxon_pool <- ster_otu_rel_abund %>%
  group_by(taxon) %>%
  summarize(median=median(rel_abund),.groups="drop") 

ster_otu_habitat_rel_abund <- inner_join(ster_otu_rel_abund, ster_taxon_pool, by="taxon") %>%
  group_by(sample_id, Salinity, taxon, phylum) %>%
  summarize(rel_abund = sum(rel_abund),
            median = max(median),
            .groups="drop") %>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, median, .desc=FALSE)

#Cor test is used to do a two.sided spearman rank correlation of ASV RA and Salinity
experiment_significance <- ster_otu_habitat_rel_abund %>%
  nest(data = -taxon) %>%
  mutate(experiment_tests = map(.x=data, 
                                ~cor.test(~ rel_abund + Salinity, alternative="two.sided", method="spearman", data=.x) %>%
                                  tidy())) %>%
  unnest(experiment_tests) %>%
  mutate(p.experiment = p.adjust(p.value, method="BH")) %>%
  select(taxon, p.experiment, estimate)
tax_phylum<-taxonomy %>%
  group_by(taxon) %>%
  select(taxon,phylum)
ggplot_Salster_corr<-inner_join(experiment_significance, nicelooking, by="taxon")
ggplot_Salster_corr<-inner_join(ggplot_Salster_corr, tax_phylum, by="taxon")
##Plot data, I added taxa label 
ggplot(ggplot_Salster_corr, aes(estimate, freq)) + geom_jitter(aes(color=phylum, size=median)) + theme_bw() + 
  scale_color_manual(values =cbbPalette)
#For Density plot
ggplot(plot, aes(estimate)) + geom_density()
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##want to seperate two filters to plot seperately since they are distinct
pre_otu_rel_abund<-  otu_rel_abund %>%
  filter(str_detect(sample_id, "Pre"))

pre_nicelooking<-pre_otu_rel_abund %>%
  group_by(taxon) %>%
  summarize(mean=mean(rel_abund),median=median(rel_abund),freq = sum(rel_abund != 0) / sum(rel_abund >= 0), maxRA=max(rel_abund),.groups="drop")
pre_taxon_pool <- pre_otu_rel_abund %>%
  group_by(taxon) %>%
  summarize(median=median(rel_abund),.groups="drop") 

pre_otu_habitat_rel_abund <- inner_join(pre_otu_rel_abund, pre_taxon_pool, by="taxon") %>%
  group_by(sample_id, Salinity, taxon, phylum) %>%
  summarize(rel_abund = sum(rel_abund),
            median = max(median),
            .groups="drop") %>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, median, .desc=FALSE))


experiment_significance <- pre_otu_habitat_rel_abund %>%
  nest(data = -taxon) %>%
  mutate(experiment_tests = map(.x=data, 
                                ~cor.test(~ rel_abund + Salinity, alternative="two.sided", method="spearman", data=.x) %>%
                                  tidy())) %>%
  unnest(experiment_tests) %>%
  mutate(p.experiment = p.adjust(p.value, method="BH")) %>%
  select(taxon, p.experiment, estimate)
tax_phylum<-taxonomy %>%
  group_by(taxon) %>%
  select(taxon,phylum)
ggplot_Salpre_corr<-inner_join(experiment_significance, pre_nicelooking, by="taxon")
ggplot_Salpre_corr<-inner_join(ggplot_Salpre_corr, tax_phylum, by="taxon")
ggplot(ggplot_Salpre_corr, aes(estimate)) +geom_density()+ theme_bw()
ggplot(plot, aes(estimate, freq)) + geom_jitter(aes(color=phylum, size=median)) + theme_bw() + 
  scale_color_manual(values =cbbPalette) 
#Plot density based on rho
ggplot(plot, aes(estimate)) + geom_density()+ theme_bw()+theme(legend.position = "none", axis.text.y = element_blank())
