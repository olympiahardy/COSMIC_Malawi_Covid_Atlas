nasal_UMAP_colours <- c("#2ca02c", "#1f77b4", "#87422F", "#D4674A", "#9467bd", "#17becf","#FFC10D","#DBC00B", "#ff7f0e", "#d62728", "#11AD82", "#c7c7c7")
names(nasal_UMAP_colours) <- c("CD4+ T cells", "CD8+ T cells", "Goblet cells", "Secretory cells", "Neutrophils", "Basal cells", "VEGFAhigh Squamous cells", 
                               "SPRR2Dhigh Squamous cells","Macrophages","Erthrocytes","Cilliated cells","Neurons")


###Prepare data for plotting
circ_data <- prepare_circlize_data(nasal.atlas, scale = 0.75)
set.seed(1234)
cluster_colors<- nasal_UMAP_colours
group_colors<- c("#F1CDB1", "#CBDFB8", "#8FA9DF")
names(group_colors) <- c("COVID-19", "Non-Pneumonia", "Pneumonia")
infection_colors<-c("#74ADD1", "#D73027")
names(infection_colors) <- c("Positive", "Negative")

###plot and save figures
pdf('circlize_plot_nasal.pdf', width = 6, height = 6)
plot_circlize(circ_data,do.label = F, pt.size = 0.1, col.use = cluster_colors ,bg.color = 'white', kde2d.n = 200, repel = T, label.cex = 0, contour.levels = c(0,0),
              contour.nlevels = 0)
add_track(circ_data, group = "Group", colors = group_colors, track_num = 2, labs = "Group", srt = 180) ## can change it to one of the columns in the meta data of your seurat object
add_track(circ_data, group = "HIV_Status",colors = infection_colors, track_num = 3, labs = "HIV Status", srt = 180)
dev.off()

blood_UMAP_colours <- c("#9467bd", "#A33E3E", "#1f77b4", "#FFC10D", "#aec7e8", "#e377c2", "#2ca02c", "#573D70", "#A33E7B")
names(blood_UMAP_colours) <- c("Neutrophils", "Reticulocytes",  "CD8+ T cells", "Monocytes", "NK cells", "B cells","CD4+ T cells","CMP/GMP","Platelets")

###Prepare data for plotting - blood atlas
circ_data <- prepare_circlize_data(blood.atlas, scale = 0.75)
set.seed(1234)
cluster_colors<- blood_UMAP_colours
group_colors<- c("#F1CDB1", "#CBDFB8", "#8FA9DF")
names(group_colors) <- c("COVID-19", "Non-Pneumonia", "Pneumonia")
infection_colors<-c("#74ADD1", "#D73027")
names(infection_colors) <- c("Positive", "Negative")

###plot and save figures
pdf('circlize_plot_blood.pdf', width = 6, height = 6)
plot_circlize(circ_data,do.label = F, pt.size = 0.1, col.use = cluster_colors ,bg.color = 'white', kde2d.n = 200, repel = T, label.cex = 0, contour.levels = c(0,0),
              contour.nlevels = 0)
add_track(circ_data, group = "Group", colors = group_colors, track_num = 2, labs = "Group", srt = 180) ## can change it to one of the columns in the meta data of your seurat object
add_track(circ_data, group = "HIV_Status",colors = infection_colors, track_num = 3, labs = "HIV Status", srt = 180)
dev.off()

## Nasal stuff ##
cosmic.nasal.covid.list <- list()
cosmic.nasal.covid.list.up <- list()
cosmic.nasal.covid.list.down <- list()

nasal.de <- subset(nasal.atlas, subset = Group %in% c("COVID-19", "Pneumonia"))
nasal.de$celltype_source <- paste(nasal.de$celltype_final, nasal.de$Group, sep = "_")
# get the amount of clusters
num_cluster=length(unique(nasal.de$celltype_final))
Idents(nasal.de) <- nasal.de$celltype_final

cluster.ids <- unique(Idents(nasal.de))
Idents(nasal.de) <- "celltype_source"
DefaultAssay(nasal.de) <- "SCT"
for (x in 1:num_cluster) {
  n=cluster.ids[x]
  print(n)
  check_sanity = FALSE
  ## diff expression \
  COVID.response <- FindMarkers(nasal.de,ident.1 = paste(n,"_COVID-19", sep=""), assay = "RNA", slot = "data", ident.2 = paste(n,"_Pneumonia", sep=""), verbose = T, logfc.threshold = 0.25, min.pct = 0.1)
  COVID.response$gene <- rownames(COVID.response)
  
  COVID.response <- COVID.response[-c(grep("^RPS", COVID.response$gene),
                                      grep("^RPL", COVID.response$gene),
                                      grep("^RP", COVID.response$gene),
                                      grep("^HB", COVID.response$gene),
                                      grep("^MT-", COVID.response$gene),
                                      grep("MALAT1", COVID.response$gene),
                                      grep("^MTR", COVID.response$gene),
                                      grep("^RNA18S5", COVID.response$gene),
                                      grep("^RNA28S5", COVID.response$gene),
                                      grep("^RNA1\\dS5", COVID.response$gene)),]
  cosmic.nasal.covid.list[[n]] <- COVID.response
  COVID.response1 <- subset(COVID.response, COVID.response$avg_log2FC > 0.5 & COVID.response$p_val_adj < 0.05)
  COVID.response1 <- COVID.response1[order(COVID.response1$avg_log2FC, decreasing = T),]
  COVID.response2 <- subset(COVID.response, COVID.response$avg_log2FC < -0.5 & COVID.response$p_val_adj < 0.05)
  COVID.response2 <- COVID.response2[order(COVID.response2$avg_log2FC),]
  
  cosmic.nasal.covid.list.up[[n]] <- COVID.response1
  cosmic.nasal.covid.list.down[[n]] <- COVID.response2
}
names(cosmic.nasal.covid.list.up) <- cluster.ids
names(cosmic.nasal.covid.list.down) <- cluster.ids
names(cosmic.nasal.covid.list) <- cluster.ids

gene.labels <- c("CXCL5", "SPP1", "LGALS1", "HLA-DPB1", "HLA-DQA1", "IFI27", "TMSB10", "APOE", "C1QB", "IRF1", 
                 "CCL7","PITPNC1","SLAMF1", "TMEM158",   "LINC01619", "ADGRE2" ,   "ATP2B1" ,   "MET" ,      "ANPEP",     "FNIP2")

COVID.response1 <- subset(cosmic.nasal.covid.list$Macrophages, cosmic.nasal.covid.list$Macrophages$avg_log2FC > 0 &cosmic.nasal.covid.list$Macrophages$p_val_adj < 0.05)
COVID.response2 <- subset(cosmic.nasal.covid.list$Macrophages, cosmic.nasal.covid.list$Macrophages$avg_log2FC < 0 &cosmic.nasal.covid.list$Macrophages$p_val_adj < 0.05)
c <- rbind(COVID.response1, COVID.response2)
c <- rbind(cosmic.nasal.covid.list.up$Macrophages, cosmic.nasal.covid.list.down$Macrophages)


EnhancedVolcano(cosmic.nasal.covid.list$Macrophages,
                lab = rownames(cosmic.nasal.covid.list$Macrophages),
                selectLab = gene.labels,
                x = 'avg_log2FC',
                y = 'p_val',
                pCutoff = 0.01,
                parseLabels = TRUE,
                drawConnectors = T,
                max.overlaps = Inf,
                xlim = c(-3,3),
                colAlpha = 0.3)

## blood stuff ##
cosmic.blood.covid.list <- list()
cosmic.blood.covid.list.up <- list()
cosmic.blood.covid.list.down <- list()

blood.de <- subset(blood.atlas, subset = Group %in% c("COVID-19", "Pneumonia"))

blood.de$celltype_source <- paste(blood.de$celltype, blood.de$Group, sep = "_")
# get the amount of clusters
num_cluster=length(unique(blood.de$celltype))
Idents(blood.de) <- blood.de$celltype

cluster.ids <- unique(Idents(blood.de))
Idents(blood.de) <- "celltype_source"
DefaultAssay(blood.de) <- "SCT"
for (x in 1:num_cluster) {
  n=cluster.ids[x]
  print(n)
  check_sanity = FALSE
  ## diff expression \
  COVID.response <- FindMarkers(blood.de,ident.1 = paste(n,"_COVID-19", sep=""), assay = "RNA", slot = "data", ident.2 = paste(n,"_Pneumonia", sep=""), verbose = T, logfc.threshold = 0.25, min.pct = 0.1)
  COVID.response$gene <- rownames(COVID.response)
  
  COVID.response <- COVID.response[-c(grep("^RPS", COVID.response$gene),
                                      grep("^RPL", COVID.response$gene),
                                      grep("^RP", COVID.response$gene),
                                      grep("^HB", COVID.response$gene),
                                      grep("^MT-", COVID.response$gene),
                                      grep("MALAT1", COVID.response$gene),
                                      grep("^MTR", COVID.response$gene),
                                      grep("^RNA18S5", COVID.response$gene),
                                      grep("^RNA28S5", COVID.response$gene),
                                      grep("^RNA1\\dS5", COVID.response$gene)),]
  cosmic.blood.covid.list[[n]] <- COVID.response
  COVID.response1 <- subset(COVID.response, COVID.response$avg_log2FC > 0.5 & COVID.response$p_val_adj < 0.05)
  COVID.response1 <- COVID.response1[order(COVID.response1$avg_log2FC, decreasing = T),]
  COVID.response2 <- subset(COVID.response, COVID.response$avg_log2FC < -0.5 & COVID.response$p_val_adj < 0.05)
  COVID.response2 <- COVID.response2[order(COVID.response2$avg_log2FC),]
  
  cosmic.blood.covid.list.up[[n]] <- COVID.response1
  cosmic.blood.covid.list.down[[n]] <- COVID.response2
}
names(cosmic.blood.covid.list.up) <- cluster.ids
names(cosmic.blood.covid.list.down) <- cluster.ids
names(cosmic.blood.covid.list) <- cluster.ids

gene1 <- cosmic.blood.covid.list.up$Monocytes$gene[1:10]
gene2 <- cosmic.blood.covid.list.down$Monocytes$gene[1:10]
gene.labels <- c(gene1,gene2)

EnhancedVolcano(cosmic.blood.covid.list$Monocytes,
                lab = rownames(cosmic.blood.covid.list$Monocytes),
                selectLab = gene.labels,
                x = 'avg_log2FC',
                y = 'p_val',
                pCutoff = 0.01,
                parseLabels = TRUE,
                drawConnectors = T,
                max.overlaps = Inf,
                xlim = c(-1.5,1.5),
                colAlpha = 0.3)

gene1 <- cosmic.blood.covid.list.up$`CD4+ T cells`$gene[1:10]
gene2 <- cosmic.blood.covid.list.down$`CD4+ T cells`$gene[1:10]
gene.labels <- c(gene1,gene2)
EnhancedVolcano(cosmic.blood.covid.list$`CD4+ T cells`,
                lab = rownames(cosmic.blood.covid.list$`CD4+ T cells`),
                selectLab = gene.labels,
                x = 'avg_log2FC',
                y = 'p_val',
                pCutoff = 0.01,
                parseLabels = TRUE,
                drawConnectors = T,
                max.overlaps = Inf,
                xlim = c(-1.5,1.5),
                colAlpha = 0.3)






