library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(multinichenetr)
library(Seurat)

# CCI analysis #
lung_all <- readRDS("/datastore/Olympia/lung_all_12923.rds")

lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% distinct(ligand, receptor) %>% mutate(ligand = make.names(ligand), receptor = make.names(receptor))
ligand_target_matrix = readRDS("/datastore/Olympia/ligand_target_matrix_nsga2r_final.rds")
colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% make.names()
rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% make.names()

lung_all_sce <- as.SingleCellExperiment(lung_all)
lung_all_sce = alias_to_symbol_SCE(lung_all_sce, "human") %>% makenames_SCE()

sample_id = "Sample"
group_id = "Group"
celltype_id = "celltype_annotation"
batches = NA
covariates = NA

celltypes <- c("Adventitial fibroblasts" , "Alveolar fibroblasts", "Alveolar macrophages", "AT1", "AT2",
               "CD16- Neutrophils","CD16+ Neutrophils", "CD8+ T cells", "Interstitial macrophages", "Monocyte-derived macrophages", "Naive CD4+ T cells","NK cells",
               "T reg", "Th1", "Venous endothelium")

lung_all_sce = lung_all_sce[, SummarizedExperiment::colData(lung_all_sce)[,celltype_id] %in% celltypes]
lung_all_sce$celltype_annotation <- lung_all_sce$celltype_annotation %>% make.names()
lung_all_sce$Sample <- lung_all_sce$Sample %>% make.names()
senders_oi = SummarizedExperiment::colData(lung_all_sce)[,celltype_id] %>% unique()
receivers_oi = SummarizedExperiment::colData(lung_all_sce)[,celltype_id] %>% unique()

min_cells = 5

lung_all_sce$Group[lung_all_sce$Group == "COVID-19"] <- "COVID19"

contrasts_oi = c("'COVID19-Pneumonia','Pneumonia-COVID19'")

contrast_tbl = tibble(
  contrast = c("COVID19-Pneumonia","Pneumonia-COVID19"), 
  group = c("COVID19","Pneumonia"))

lung_all_sce = lung_all_sce[, SummarizedExperiment::colData(lung_all_sce)[,group_id] %in% contrast_tbl$group]

logFC_threshold = 0.50
p_val_threshold = 0.05
fraction_cutoff = 0.10

p_val_adj = FALSE 
empirical_pval = FALSE

top_n_target = 250

prioritizing_weights_DE = c("de_ligand" = 1,
                            "de_receptor" = 1)
prioritizing_weights_activity = c("activity_scaled" = 2)

prioritizing_weights_expression_specificity = c("exprs_ligand" = 2,
                                                "exprs_receptor" = 2)

prioritizing_weights_expression_sufficiency = c("frac_exprs_ligand_receptor" = 1)

prioritizing_weights_relative_abundance = c( "abund_sender" = 0,
                                             "abund_receiver" = 0)
prioritizing_weights = c(prioritizing_weights_DE, 
                         prioritizing_weights_activity, 
                         prioritizing_weights_expression_specificity,
                         prioritizing_weights_expression_sufficiency, 
                         prioritizing_weights_relative_abundance)
n.cores = min(8, union(senders_oi, receivers_oi) %>% length())
multinichenet_output = multi_nichenet_analysis(sce = lung_all_sce, celltype_id = celltype_id, sample_id = sample_id, group_id = group_id, 
                                               lr_network = lr_network, ligand_target_matrix = ligand_target_matrix, contrasts_oi = contrasts_oi, contrast_tbl = contrast_tbl, batches = batches, covariates = covariates,
                                               prioritizing_weights = prioritizing_weights, min_cells = min_cells, logFC_threshold = logFC_threshold, p_val_threshold = p_val_threshold,  
                                               fraction_cutoff = fraction_cutoff, p_val_adj = p_val_adj, empirical_pval = empirical_pval, top_n_target = top_n_target, n.cores = n.cores, sender_receiver_separate = FALSE, verbose = TRUE)


### CIRCOS for groups ###
prioritized_tbl_oi_A_30 = get_top_n_lr_pairs(multinichenet_output$prioritization_tables, 50, groups_oi = "COVID19")
circos_COVID = make_circos_one_group_custom(prioritized_tbl_oi_A_30, colors, colors)
circos_COVID$COVID19       

# Interactions with macs being senders 
group_oi = "COVID19"


prioritized_tbl_oi_M_50 = get_top_n_lr_pairs(multinichenet_output$prioritization_tables, 25, groups_oi = group_oi, receivers_oi = c("AT1", "AT2", "Adventitial.fibroblasts",
                                                                                                                                    "Alveolar.fibroblasts", "Interstitial.macrophages"), senders_oi = "Alveolar.macrophages")

plot_oi = make_sample_lr_prod_activity_plots2plot(multinichenet_output$prioritization_tables, prioritized_tbl_oi_M_50)
plot_oi

# Interactions with macs being recievers 
prioritized_tbl_oi_M_50 = get_top_n_lr_pairs(multinichenet_output$prioritization_tables, 25, groups_oi = group_oi, receivers_oi = "Alveolar.macrophages", senders_oi = c("AT2"))
plot_oi = make_sample_lr_prod_activity_plots2plot(multinichenet_output$prioritization_tables, prioritized_tbl_oi_M_50)
plot_oi

combined_plot = make_ligand_activity_target_plot(group_oi, prioritized_tbl_oi_M_50, receiver_oi = c("Alveolar.macrophages"), multinichenet_output$prioritization_tables, multinichenet_output$ligand_activities_targets_DEgenes, contrast_tbl, multinichenet_output$grouping_tbl, multinichenet_output$celltype_info, ligand_target_matrix, plot_legend = FALSE)
combined_plot$combined_plot

# Interactions with neuts being recievers 
prioritized_tbl_oi_M_50 = get_top_n_lr_pairs(multinichenet_output$prioritization_tables, 25, groups_oi = group_oi, receivers_oi = "CD16..Neutrophils", senders_oi = c("Venous.endothelium"))
plot_oi = make_sample_lr_prod_activity_plots2plot(multinichenet_output$prioritization_tables, prioritized_tbl_oi_M_50)
plot_oi

combined_plot = make_ligand_activity_target_plot(group_oi, prioritized_tbl_oi_M_50, receiver_oi = c("CD16..Neutrophils"), multinichenet_output$prioritization_tables, multinichenet_output$ligand_activities_targets_DEgenes, contrast_tbl, multinichenet_output$grouping_tbl, multinichenet_output$celltype_info, ligand_target_matrix, plot_legend = FALSE)
combined_plot$combined_plot

# Interactions with neuts being senders
prioritized_tbl_oi_M_50 = get_top_n_lr_pairs(multinichenet_output$prioritization_tables, 25, groups_oi = group_oi, receivers_oi = c("Venous.endothelium"), senders_oi = "CD16..Neutrophils")
plot_oi = make_sample_lr_prod_activity_plots2plot(multinichenet_output$prioritization_tables, prioritized_tbl_oi_M_50)
plot_oi


