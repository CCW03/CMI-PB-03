rm(list=ls())

library(readr)
library(dplyr)
library(tidyverse)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
base_dir = "../"
dir_raw_training <- paste0(base_dir, "data/raw_training_dataset/")
dir_raw_prediction <- paste0(base_dir, "data/raw_prediction_dataset/")
select <- dplyr::select

## subject files
d2020_subject <- read_tsv(paste0(dir_raw_training, "2020LD_subject.tsv"))
d2021_subject <- read_tsv(paste0(dir_raw_training, "2021LD_subject.tsv"))
d2022_subject <- read_tsv(paste0(dir_raw_training, "2022LD_subject.tsv"))
d2023_subject <- read_tsv(paste0(dir_raw_prediction, "2023BD_subject.tsv"))

subject_training <- bind_rows(d2020_subject, d2021_subject, d2022_subject, d2023_subject)

## specimen files
d2020_specimen <- read_tsv(paste0(dir_raw_training, "2020LD_specimen.tsv"))
d2021_specimen <- read_tsv(paste0(dir_raw_training, "2021LD_specimen.tsv"))
d2022_specimen <- read_tsv(paste0(dir_raw_training, "2022LD_specimen.tsv"))
d2023_specimen <- read_tsv(paste0(dir_raw_prediction, "2023BD_specimen.tsv"))

specimen_training <-bind_rows(cbind(data = "2020", d2020_specimen), 
                              cbind(data = "2021", d2021_specimen),
                              cbind(data = "2022", d2022_specimen),
                              cbind(data = "2023", d2023_specimen))

## create new object subject_specimen
subject_specimen_training <- specimen_training %>%
  left_join(subject_training) %>%
  mutate(timepoint = planned_day_relative_to_boost)

subj_spec <- subject_specimen_training %>%select(data, subject_id,specimen_id)

## plasma_antibody_levels files
d2020_plasma_antibody_levels <- read_tsv(paste0(dir_raw_training, "2020LD_plasma_ab_titer.tsv")) %>%
  mutate(isotype_antigen = paste0(isotype,"_", antigen))  %>%
  dplyr::select(-antigen, -isotype)

d2021_plasma_antibody_levels <- read_tsv(paste0(dir_raw_training, "2021LD_plasma_ab_titer.tsv")) %>%
  mutate(isotype_antigen = paste0(isotype,"_", antigen))  %>%
  dplyr::select(-antigen, -isotype)

d2022_plasma_antibody_levels <- read_tsv(paste0(dir_raw_training, "2022LD_plasma_ab_titer.tsv")) %>%
  mutate(isotype_antigen = paste0(isotype,"_", antigen))  %>%
  dplyr::select(-antigen, -isotype)

d2023_plasma_antibody_levels <- read_tsv(paste0(dir_raw_prediction, "2023BD_plasma_ab_titer.tsv")) %>%
  mutate(isotype_antigen = paste0(isotype,"_", antigen))  %>%
  dplyr::select(-antigen, -isotype)

plasma_antibody_levels_common_features <- Reduce(intersect, list(unique(d2020_plasma_antibody_levels$isotype_antigen), 
                                                                 unique(d2021_plasma_antibody_levels$isotype_antigen),
                                                                 unique(d2022_plasma_antibody_levels$isotype_antigen),
                                                                 unique(d2023_plasma_antibody_levels$isotype_antigen))) 

plasma_antibody_levels_long <- bind_rows(cbind(data = "2020", d2020_plasma_antibody_levels), 
                                         cbind(data = "2021", d2021_plasma_antibody_levels),
                                         cbind(data = "2022", d2022_plasma_antibody_levels),
                                         cbind(data = "2023", d2023_plasma_antibody_levels)) %>%
  filter(isotype_antigen %in% plasma_antibody_levels_common_features)

plasma_antibody_levels_wide <- plasma_antibody_levels_long %>%
  dplyr::select(specimen_id, isotype_antigen, MFI_normalised) %>%
  pivot_wider(names_from = isotype_antigen, values_from = MFI_normalised, )

## plasma_cytokine_concentrations files
d2020_plasma_cytokine_concentrations <- read_tsv(paste0(dir_raw_training, "2020LD_plasma_cytokine_concentration_by_olink.tsv"))
d2021_plasma_cytokine_concentrations <- read_tsv(paste0(dir_raw_training, "2021LD_plasma_cytokine_concentration_by_olink.tsv"))
d2022_plasma_cytokine_concentrations <- read_tsv(paste0(dir_raw_training, "2022LD_plasma_cytokine_concentration_by_olink.tsv"))
d2023_plasma_cytokine_concentrations <- read_tsv(paste0(dir_raw_prediction, "2023BD_plasma_cytokine_concentration_by_olink.tsv"))

plasma_cytokine_concentrations_common_features <- Reduce(intersect, list(unique(d2020_plasma_cytokine_concentrations$protein_id), 
                                                                         unique(d2021_plasma_cytokine_concentrations$protein_id), 
                                                                         unique(d2022_plasma_cytokine_concentrations$protein_id),
                                                                         unique(d2023_plasma_cytokine_concentrations$protein_id))) 

plasma_cytokine_concentrations_long <-bind_rows(cbind(data = "2020", d2020_plasma_cytokine_concentrations), 
                                                cbind(data = "2021", d2021_plasma_cytokine_concentrations),
                                                cbind(data = "2022", d2022_plasma_cytokine_concentrations),
                                                cbind(data = "2023", d2023_plasma_cytokine_concentrations)) %>%
  filter(protein_id %in% plasma_cytokine_concentrations_common_features)

plasma_cytokine_concentrations_wide <- plasma_cytokine_concentrations_long %>%
  dplyr::select(specimen_id, protein_id, concentration) %>%
  pivot_wider(names_from = protein_id, values_from = concentration)

## pbmc_cell_frequency files
d2020_pbmc_cell_frequency <- read_tsv(paste0(dir_raw_training, "2020LD_pbmc_cell_frequency.tsv"))
d2021_pbmc_cell_frequency <- read_tsv(paste0(dir_raw_training, "2021LD_pbmc_cell_frequency.tsv"))
d2022_pbmc_cell_frequency <- read_tsv(paste0(dir_raw_training, "2022LD_pbmc_cell_frequency.tsv"))
d2023_pbmc_cell_frequency <- read_tsv(paste0(dir_raw_prediction, "2023BD_pbmc_cell_frequency.tsv"))

pbmc_cell_frequency_common_features <- Reduce(intersect, list(unique(d2020_pbmc_cell_frequency$cell_type_name), 
                                                              unique(d2021_pbmc_cell_frequency$cell_type_name), 
                                                              unique(d2022_pbmc_cell_frequency$cell_type_name), 
                                                              unique(d2023_pbmc_cell_frequency$cell_type_name))) 

pbmc_cell_frequency_long <-bind_rows(cbind(data = "2020", d2020_pbmc_cell_frequency), 
                                     cbind(data = "2021", d2021_pbmc_cell_frequency),
                                     cbind(data = "2022", d2022_pbmc_cell_frequency),
                                     cbind(data = "2023", d2023_pbmc_cell_frequency)) %>%
  filter(cell_type_name %in% pbmc_cell_frequency_common_features)

pbmc_cell_frequency_wide <- pbmc_cell_frequency_long %>%
  pivot_wider(names_from = cell_type_name, values_from = percent_live_cell)

d2020_pbmc_gene_expression <- read_tsv(paste0(dir_raw_training, "2020LD_pbmc_gene_expression.tsv"))
d2021_pbmc_gene_expression <- read_tsv(paste0(dir_raw_training, "2021LD_pbmc_gene_expression.tsv"))
d2022_pbmc_gene_expression <- read_tsv(paste0(dir_raw_training, "2022LD_pbmc_gene_expression.tsv"))
d2023_pbmc_gene_expression <- read_tsv(paste0(dir_raw_prediction, "2023BD_pbmc_gene_expression.tsv"))

colnames(d2020_pbmc_gene_expression) = c('versioned_ensembl_gene_id','specimen_id','tpm','raw_count')
colnames(d2021_pbmc_gene_expression) = c('versioned_ensembl_gene_id','specimen_id','tpm','raw_count')
colnames(d2022_pbmc_gene_expression) = c('versioned_ensembl_gene_id','specimen_id','tpm','raw_count')
colnames(d2023_pbmc_gene_expression) = c('versioned_ensembl_gene_id','specimen_id','tpm','raw_count')

pbmc_gene_expression_common_features <- Reduce(intersect, list(unique(d2020_pbmc_gene_expression$versioned_ensembl_gene_id), 
                                                               unique(d2021_pbmc_gene_expression$versioned_ensembl_gene_id), 
                                                               unique(d2022_pbmc_gene_expression$versioned_ensembl_gene_id), 
                                                               unique(d2023_pbmc_gene_expression$versioned_ensembl_gene_id))) 

pbmc_gene_expression_long <-bind_rows(cbind(data = "2020", d2020_pbmc_gene_expression), 
                                      cbind(data = "2021", d2021_pbmc_gene_expression),
                                      cbind(data = "2022", d2022_pbmc_gene_expression),
                                      cbind(data = "2023", d2023_pbmc_gene_expression)) %>%
  filter(versioned_ensembl_gene_id %in% pbmc_gene_expression_common_features)

pbmc_gene_expression_wide <- pbmc_gene_expression_long %>%
  dplyr::select(specimen_id, versioned_ensembl_gene_id, raw_count) %>%
  pivot_wider(names_from = versioned_ensembl_gene_id, values_from = raw_count)


master_database_data <- list(
  
  subject_specimen = subject_specimen_training,
  plasma_antibody_levels = list(
    
    wide = plasma_antibody_levels_wide,
    long = plasma_antibody_levels_long
  ),
  plasma_cytokine_concentrations = list(
    
    wide = plasma_cytokine_concentrations_wide,
    long = plasma_cytokine_concentrations_long
  ),
  pbmc_cell_frequency_wide = list(
    
    wide = pbmc_cell_frequency_wide,
    long = pbmc_cell_frequency_long
  ),
  pbmc_gene_expression_wide = list(
    
    wide = pbmc_gene_expression_wide,
    long = pbmc_gene_expression_long
  )
)

# saveRDS(master_database_data, file = "master_harmonized_training_data.RDS")