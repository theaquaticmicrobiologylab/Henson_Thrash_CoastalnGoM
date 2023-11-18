# Code based on Filter by significance exercise using diversity data
# Code is from Dr. Pat Schloss, 
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

#needed function from below
get_max_abund <- function(x) {
  
  x %>%
    group_by(habitat) %>%
    summarize(third_q = quantile(rel_abund, prob=0.75), .groups = "drop") %>%
    summarize(max_quartile = max(third_q)) %>%
    pull(max_quartile)
}


metadata <- read_excel(path="/Path/to/file/Metadata.csv") %>%
  rename(sample_id = Sample)

##This OTU file should be Sample Names in column 1 and OTU# in row 1. 
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
         genus = str_replace(string=genus,
                             pattern="(.*)",
                             replacement="*\\1*"),
         genus = str_replace(string=genus,
                             pattern="\\*(.*)_unclassified\\*",
                             replacement="Unclassified *\\1*"),
         taxon = glue("{genus} ({pretty_otu})"),
         taxon = str_replace_all(taxon, "_", " ")) %>%
  select(otu, taxon)

##The mutate function should test your different salinity groupings
otu_rel_abund <- inner_join(metadata, otu_counts, by="sample_id") %>%
  inner_join(., taxonomy, by="otu") %>%
  group_by(sample_id) %>%
  mutate(rel_abund = 100*count / sum(count)) %>%
  ungroup() %>%
  select(-count) %>%
  mutate(habitat = factor(Water,
                          levels=c("Fresh",
                                   "Low Salinity",
                                   "High Salinity")
  ))

##want to seperate two filters to plot seperately since they are distinct
pre_otu_rel_abund<-  otu_rel_abund %>%
  filter(str_detect(sample_id, "Pre"))

pre_taxon_pool <- pre_otu_rel_abund %>%
  group_by(habitat, taxon) %>%
  summarize(median=median(rel_abund), .groups="drop") %>%
  group_by(taxon) %>%
  summarize(pool = max(median) < 1.0,
            median = max(median),
            .groups="drop")

pre_otu_habitat_rel_abund <- inner_join(pre_otu_rel_abund, pre_taxon_pool, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", as.character(taxon))) %>%
  group_by(sample_id, habitat, taxon) %>%
  summarize(rel_abund = sum(rel_abund),
            median = max(median),
            .groups="drop") %>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, median, .desc=FALSE))

pre_otu_habitat_rel_abund_noother<-pre_otu_habitat_rel_abund %>%
  filter(taxon != "Other")
  
pre_otu_plot<- pre_otu_habitat_rel_abund_noother %>%
  mutate(rel_abund = if_else(rel_abund == 0,
                             2/3 * lod,
                             rel_abund)) %>%
  mutate(habitat = fct_rev(as_factor(habitat))) %>%
  ggplot(aes(y=taxon, x=rel_abund, color=habitat)) +
  #geom_point(position = position_jitterdodge(jitter.width = 0, dodge.width = 0.5),
  #           pch=21, stroke=0.1, size=1.5) +
  geom_vline(xintercept = lod, size=0.25) +
  stat_summary(fun.data=median_hilow, geom = "pointrange",
               fun.args=list(conf.int=0.5),
               position = position_dodge(width=0.5),
               size=0.25) +
  #coord_trans(x="log10") +
  #scale_x_continuous(limits=c(0.1, 100),
  #                   breaks=c(0.1, 1, 10, 100),
  #                   labels=c(0.1, 1, 10, 100)) +
  scale_color_manual(name=NULL,
                    values=c("steelblue", "orange2", "red", "purple"),
                    breaks=c("Fresh", "Low Salinity", "High Salinity"),
                    labels=c("Fresh", "Low Salinity", "High Salinity")) +
  labs(y=NULL,
       x="Relative Abundance (%)") +
  theme_classic() +
  theme(axis.text.y = element_markdown(size = 6),
        legend.text = element_markdown(),
        # legend.position = c(0.8, 0.6),
        legend.background = element_rect(color="black", fill = NA),
        legend.margin = margin(t=-5, r=3, b=3)
  )
pre_otu_plot

experiment_significance <- pre_otu_habitat_rel_abund_noother %>%
  nest(data = -taxon) %>%
  mutate(experiment_tests = map(.x=data, 
                                ~kruskal.test(rel_abund ~ habitat, data = .x) %>%
                                  tidy())) %>%
  unnest(experiment_tests) %>%
  mutate(p.experiment = p.adjust(p.value, method="BH")) %>%
  select(taxon, data, p.experiment) %>%
  filter(p.experiment < 0.05)

pairwise_test <- experiment_significance %>%
  mutate(max_abund = map(.x=data, ~get_max_abund(.x))) %>%
  mutate(pairwise_tests = map(.x=data, ~pairwise.wilcox.test(x=.x$rel_abund,
                                                             g=.x$habitat,
                                                             p.adjust.method = "BH") %>%
                                tidy())) %>%
  unnest(pairwise_tests) %>%
  filter(p.value < 0.05) %>%
  select(taxon, group1, group2, p.value, max_abund) %>%
  mutate(pos = as.numeric(taxon),
         y = if_else(group1 == "High Salininty", pos + 0.3, pos),
         yend = if_else(group2 == "Fresh", pos, pos - 0.3),
         x = case_when(group1 == "High Salinity" &
                         group2 == "Fresh" ~ as.numeric(max_abund) * 1,
                       group1 == "Low Salinity" &
                         group2 == "Fresh" ~ as.numeric(max_abund) * 2,
                       group1 == "High Salinity" &
                         group2 == "Low Salinity" ~ as.numeric(max_abund) * 3),
         xend = x,
         x_star = 1.1+x,
         y_star = case_when(group1 == "High Salinity" &
                              group2 == "Fresh" ~ pos + 0.05,
                            group1 == "Low Salinity" &
                              group2 == "Fresh" ~ pos - 0.05,
                            group1 == "High Salinity" &
                              group2 == "Low Salinity" ~ pos - 0.15)) 

pre_otu_plot+ 
  geom_segment(data = pairwise_test, 
               aes(x=x, xend=xend, y=y, yend=yend), inherit.aes = FALSE) +
  geom_text(data = pairwise_test, 
            aes(x=x_star, y=y_star), label="*", inherit.aes = FALSE) 


##Now repeat for Ster (0.2-2.7)
ster_otu_rel_abund<-  otu_rel_abund %>%
  filter(str_detect(sample_id, "Ster"))

ster_taxon_pool <- ster_otu_rel_abund %>%
  group_by(habitat, taxon) %>%
  summarize(median=median(rel_abund), .groups="drop") %>%
  group_by(taxon) %>%
  summarize(pool = max(median) < 1.0,
            median = max(median),
            .groups="drop")

ster_otu_habitat_rel_abund <- inner_join(ster_otu_rel_abund, ster_taxon_pool, by="taxon") %>%
  mutate(taxon = if_else(pool, "Other", as.character(taxon))) %>%
  group_by(sample_id, habitat, taxon) %>%
  summarize(rel_abund = sum(rel_abund),
            median = max(median),
            .groups="drop") %>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, median, .desc=FALSE))

ster_otu_habitat_rel_abund_noother<-ster_otu_habitat_rel_abund %>%
  filter(taxon != "Other")

ster_otu_plot <- ster_otu_habitat_rel_abund_noother %>%
  mutate(rel_abund = if_else(rel_abund == 0,
                             2/3 * lod,
                             rel_abund)) %>%
  mutate(habitat = fct_rev(as_factor(habitat))) %>%
  ggplot(aes(y=taxon, x=rel_abund, color=habitat)) +
  #geom_point(position = position_jitterdodge(jitter.width = 0, dodge.width = 0.5),
  #           pch=21, stroke=0.1, size=1.5) +
  geom_vline(xintercept = lod, size=0.2) +
  stat_summary(fun.data=median_hilow, geom = "pointrange",
               fun.args=list(conf.int=0.5),
               position = position_dodge(width=0.5),
               size=0.25) +
  #coord_trans(x="log10") +
  #scale_x_continuous(limits=c(0.1, 100),
  #                   breaks=c(0.1, 1, 10, 100),
  #                   labels=c(0.1, 1, 10, 100)) +
  scale_color_manual(name=NULL,
                    values=c("steelblue", "orange2", "red", "purple"),
                    breaks=c("Fresh", "Low Salinity", "High Salinity"),
                    labels=c("Fresh", "Low Salinity", "High Salinity")) +
  labs(y=NULL,
       x="Relative Abundance (%)") +
  theme_classic() +
  theme(axis.text.y = element_markdown(size = 6),
        legend.text = element_markdown(),
        # legend.position = c(0.8, 0.6),
        legend.background = element_rect(color="black", fill = NA),
        legend.margin = margin(t=-5, r=3, b=3)
  )
ster_otu_plot

experiment_significance <- ster_otu_habitat_rel_abund_noother %>%
  nest(data = -taxon) %>%
  mutate(experiment_tests = map(.x=data, 
                                ~kruskal.test(rel_abund ~ habitat, data = .x) %>%
                                  tidy())) %>%
  unnest(experiment_tests) %>%
  mutate(p.experiment = p.adjust(p.value, method="BH")) %>%
  select(taxon, data, p.experiment) %>%
  filter(p.experiment < 0.05)

pairwise_test <- experiment_significance %>%
  mutate(max_abund = map(.x=data, ~get_max_abund(.x))) %>%
  mutate(pairwise_tests = map(.x=data, ~pairwise.wilcox.test(x=.x$rel_abund,
                                                             g=.x$habitat,
                                                             p.adjust.method = "BH",exact=FALSE) %>%
                                tidy())) %>%
  unnest(pairwise_tests) %>%
  filter(p.value < 0.05) %>%
  select(taxon, group1, group2, p.value, max_abund) %>%
  mutate(pos = as.numeric(taxon),
         y = if_else(group1 == "High Salininty", pos + 0.3, pos),
         yend = if_else(group2 == "Fresh", pos, pos - 0.3),
         x = case_when(group1 == "High Salinity" &
                         group2 == "Fresh" ~ as.numeric(max_abund) * 1,
                       group1 == "Low Salinity" &
                         group2 == "Fresh" ~ as.numeric(max_abund) * 2,
                       group1 == "High Salinity" &
                         group2 == "Low Salinity" ~ as.numeric(max_abund) * 3),
         xend = x,
         x_star = 1.1+x,
         y_star = case_when(group1 == "High Salinity" &
                              group2 == "Fresh" ~ pos + 0.05,
                            group1 == "Low Salinity" &
                              group2 == "Fresh" ~ pos - 0.05,
                            group1 == "High Salinity" &
                              group2 == "Low Salinity" ~ pos - 0.15)) 

ster_otu_plot+ 
  geom_segment(data = pairwise_test, 
               aes(x=x, xend=xend, y=y, yend=yend), inherit.aes = FALSE) +
  geom_text(data = pairwise_test, 
            aes(x=x_star, y=y_star), label="*", inherit.aes = FALSE) 
