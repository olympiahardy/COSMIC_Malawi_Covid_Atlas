## Supplemental Data Figure 2 ##
DefaultAssay(lung.atlas) <- "SCT"
cosmic.lung.markers <- FindAllMarkers(lung.atlas, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
top.cosmic.lung.markers <- cosmic.lung.markers %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = avg_log2FC)

top.stromal.markers <- subset(top.cosmic.lung.markers, subset = cluster %in% c("AT1", "AT2", "Venous endothelium", "Lymphatic endothelium", "Adventitial fibroblasts","Alveolar fibroblasts", "Myofibroblast", "Lipofibroblast", 
                                                                               "Ciliated cells", "Basal", "Mesothelial", "Smooth muscle cells",
                                                                               "Secretory cells", "Soup", "Ribosomal high cells"))

top.immune.markers <- subset(top.cosmic.lung.markers, subset = cluster %in% c("Alveolar macrophages", "Interstitial macrophages", "Monocyte-derived macrophages",
                                                                              "CD16+ Neutrophils", "CD16- Neutrophils", "CD8+ T cells", "Naive CD4+ T cells", "NK cells", "gdT cells",
                                                                              "Th1", "T reg", "Plasma cells", "B cells", "Mast cells", "Ribosomal high cells", "Cycling cells",
                                                                              "Erythrocytes"))

p1 <- DotPlot(lung.atlas, features = top.immune.markers$gene, idents = c("Alveolar macrophages", "Interstitial macrophages", "Monocyte-derived macrophages",
                                                                         "CD16+ Neutrophils", "CD16- Neutrophils", "CD8+ T cells", "Naive CD4+ T cells", "NK cells",
                                                                         "Th1", "T reg", "Plasma cells", "B cells", "Mast cells", "Ribosomal high cells", "Cycling cells",
                                                                         "Erythrocytes"))

lung.atlas@active.ident <- as.ordered(factor(lung.atlas@active.ident,
                                             levels=c("Alveolar macrophages", "Interstitial macrophages", "Monocyte-derived macrophages",
                                                      "CD16+ Neutrophils", "CD16- Neutrophils", "CD8+ T cells", "NK cells", "gdT cells", "Naive CD4+ T cells",
                                                      "Th1", "T reg", "Plasma cells", "B cells", "Mast cells", "AT1", "AT2", "Venous endothelium", "Lymphatic endothelium", 
                                                      "Adventitial fibroblasts","Alveolar fibroblasts", "Myofibroblast", "Lipofibroblast", "Ciliated cells", "Basal", "Mesothelial",
                                                      "Smooth muscle cells","Secretory cells", "Soup", "Ribosomal high cells", "Cycling cells","Erythrocytes")))

top.immune.markers$cluster <- as.ordered(factor(top.immune.markers$cluster,
                                                levels=c("Alveolar macrophages", "Interstitial macrophages", "Monocyte-derived macrophages",
                                                         "CD16+ Neutrophils", "CD16- Neutrophils", "CD8+ T cells", "NK cells", "gdT cells", "Naive CD4+ T cells",
                                                         "Th1", "T reg", "Plasma cells", "B cells", "Mast cells", "Ribosomal high cells", "Cycling cells","Erythrocytes")))
top.immune.markers <- top.immune.markers %>% 
  arrange(factor(cluster))

DotPlot(lung.atlas, features = unique(top.immune.markers$gene), cols = "RdBu", idents = c("Alveolar macrophages", "Interstitial macrophages", "Monocyte-derived macrophages",
                                                                                          "CD16+ Neutrophils", "CD16- Neutrophils", "CD8+ T cells", "Naive CD4+ T cells", "NK cells", "gdT cells",
                                                                                          "Th1", "T reg", "Plasma cells", "B cells", "Mast cells", "Ribosomal high cells", "Cycling cells",
                                                                                          "Erythrocytes")) + RotatedAxis() + coord_flip()

top.stromal.markers$cluster <- as.ordered(factor(top.stromal.markers$cluster,
                                                 levels=c("AT1", "AT2", "Venous endothelium", "Lymphatic endothelium", 
                                                          "Adventitial fibroblasts","Alveolar fibroblasts", "Myofibroblast", "Lipofibroblast", "Ciliated cells", "Basal", "Mesothelial",
                                                          "Smooth muscle cells","Secretory cells", "Soup", "Ribosomal high cells")))
top.stromal.markers <- top.stromal.markers %>% 
  arrange(factor(cluster))

DotPlot(lung.atlas, features = unique(top.stromal.markers$gene), cols = "RdBu", idents = c("AT1", "AT2", "Venous endothelium", "Lymphatic endothelium", "Adventitial fibroblasts","Alveolar fibroblasts", "Myofibroblast", "Lipofibroblast", 
                                                                                           "Ciliated cells", "Basal", "Mesothelial", "Smooth muscle cells",
                                                                                           "Secretory cells", "Soup", "Ribosomal high cells")) + RotatedAxis() + coord_flip()


## Supplemental Data Figure 3 ##
DefaultAssay(nasal.atlas) <- "SCT"
nasal.atlas.markers <- FindAllMarkers(nasal.atlas, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")

top.nasal.markers <- nasal.atlas.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)
DotPlot(nasal.atlas, features = unique(top.nasal.markers$gene), cols = "RdBu") + RotatedAxis() + coord_flip()

DefaultAssay(blood.atlas) <- "SCT"
blood.atlas.markers <- FindAllMarkers(blood.atlas, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")

top.blood.markers <- blood.atlas.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)
DotPlot(blood.atlas, features = unique(top.blood.markers$gene), cols = "RdBu") + RotatedAxis() + coord_flip()

### Celltype Proportions ###
freq.table.group <- as.data.frame(prop.table(x = table(Idents(nasal.atlas), nasal.atlas$Group), margin = 2))
freq.table.group$Var1 <- as.ordered(factor(freq.table.group$Var1,
                                           levels=c("Macrophages", "Neutrophils", "CD4+ T cells", "CD8+ T cells", "Goblet cells", "Secretory cells", "Basal cells", "VEGFAhigh Squamous cells", 
                                                    "SPRR2Dhigh Squamous cells","Cilliated cells","Erthrocytes","Neurons")))

freq.table.group$Var2 <- as.ordered(factor(freq.table.group$Var2,
                                           levels=rev(c("Non-Pneumonia", "Pneumonia", "COVID-19"))))
ggplot(data=freq.table.group, aes(x=freq.table.group$Var2, y = freq.table.group$Freq, fill=freq.table.group$Var1)) + geom_bar(stat="identity", color="black") +
  labs(x="Group", y="Proportion of cells", fill="Cell Type") + scale_x_discrete(limits = rev(levels(freq.table.group$Var2))) + scale_fill_manual(values=nasal_UMAP_colours) + theme_minimal()


freq.table.group <- as.data.frame(prop.table(x = table(Idents(blood.atlas), blood.atlas$Group), margin = 2))
freq.table.group$Var1 <- as.ordered(factor(freq.table.group$Var1,
                                           levels=c("Monocytes", "Neutrophils", "CD4+ T cells", "CD8+ T cells", "NK cells","B cells","CMP/GMP", "Reticulocytes","Platelets")))

freq.table.group$Var2 <- as.ordered(factor(freq.table.group$Var2,
                                           levels=rev(c("Non-Pneumonia", "Pneumonia", "COVID-19"))))
ggplot(data=freq.table.group, aes(x=freq.table.group$Var2, y = freq.table.group$Freq, fill=freq.table.group$Var1)) + geom_bar(stat="identity", color="black") +
  labs(x="Group", y="Proportion of cells", fill="Cell Type") + scale_x_discrete(limits = rev(levels(freq.table.group$Var2))) + scale_fill_manual(values=blood_UMAP_colours) + theme_minimal()


## Supplemental Data Figure 4 ##

# Gene reg graph
prioritized_tbl_oi = get_top_n_lr_pairs(multinichenet_output$prioritization_tables, 150, rank_per_group = F, senders_oi = c("Adventitial.fibroblasts", "Alveolar.fibroblasts", "Venous.endothelium","AT1","AT2", "Alveolar.macrophages", "Interstitial.macrophages","CD16..Neutrophils"), receivers_oi = c("Adventitial.fibroblasts", "Alveolar.fibroblasts", "Venous.endothelium","AT1","AT2", "Alveolar.macrophages", "Interstitial.macrophages","CD16..Neutrophils"))

groups_oi <- c("COVID-19", "Pneumonia")
lr_target_prior_cor_filtered = multinichenet_output$prioritization_tables$group_prioritization_tbl$group %>% unique() %>% lapply(function(group_oi){
  lr_target_prior_cor_filtered = multinichenet_output$lr_target_prior_cor %>% inner_join(multinichenet_output$ligand_activities_targets_DEgenes$ligand_activities %>% distinct(ligand, target, direction_regulation, contrast)) %>% inner_join(contrast_tbl) %>% filter(group == group_oi)
  lr_target_prior_cor_filtered_up = lr_target_prior_cor_filtered %>% filter(direction_regulation == "up") %>% filter( (rank_of_target < top_n_target) & (pearson > 0.50 | spearman > 0.50))
  lr_target_prior_cor_filtered_down = lr_target_prior_cor_filtered %>% filter(direction_regulation == "down") %>% filter( (rank_of_target < top_n_target) & (pearson < -0.50 | spearman < -0.50))
  lr_target_prior_cor_filtered = bind_rows(lr_target_prior_cor_filtered_up, lr_target_prior_cor_filtered_down)
}) %>% bind_rows()

lr_target_prior_cor_filtered_covid <- lr_target_prior_cor_filtered

colors_sender["AT2"] = "pink" # the  original yellow with white font is not very readable
colors_dutch <- c("#FFC312", "#C4E538", "#12CBC4", "#FDA7DF", "#ED4C67", "#F79F1F", "#A3CB38" ,"#1289A7", "#D980FA", "#B53471", "#EE5A24", "#009432", "#0652DD",
                  "#9980FA", "#833471" ,"#EA2027" ,"#006266", "#1B1464", "#5758BB" ,"#6F1E51")
graph_plot = make_ggraph_ligand_target_links(lr_target_prior_cor_filtered = lr_target_prior_cor_filtered_covid, prioritized_tbl_oi = prioritized_tbl_oi, colors = colors_sender)
graph_plot$plot

## Supplemental Data Figure 5 ##

# Macs being senders
combined_plot = make_ligand_activity_target_plot(group_oi, prioritized_tbl_oi_M_50, receiver_oi = c("AT1", "AT2", "Adventitial.fibroblasts",
                                                                                                    "Alveolar.fibroblasts", "Interstitial.macrophages"), multinichenet_output$prioritization_tables, multinichenet_output$ligand_activities_targets_DEgenes, contrast_tbl, multinichenet_output$grouping_tbl, multinichenet_output$celltype_info, ligand_target_matrix, plot_legend = FALSE)
combined_plot$combined_plot

# Neuts being senders
combined_plot = make_ligand_activity_target_plot(group_oi, prioritized_tbl_oi_M_50, receiver_oi = c("Venous.endothelium"), multinichenet_output$prioritization_tables, multinichenet_output$ligand_activities_targets_DEgenes, contrast_tbl, multinichenet_output$grouping_tbl, multinichenet_output$celltype_info, ligand_target_matrix, plot_legend = FALSE)
combined_plot$combined_plot





