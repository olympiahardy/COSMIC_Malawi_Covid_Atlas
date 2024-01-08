# Create broad cell types #
lung.atlas$overall_UMAP <- Idents(lung.atlas)

lung.atlas$overall_UMAP <- as.character(lung.atlas$overall_UMAP)

# Macrophages and neutrophils
lung.atlas$overall_UMAP[lung.atlas$overall_UMAP == "Interstitial macrophages"] <- "Macrophages"
lung.atlas$overall_UMAP[lung.atlas$overall_UMAP == "Monocyte-derived macrophages"] <- "Macrophages"
lung.atlas$overall_UMAP[lung.atlas$overall_UMAP == "Alveolar macrophages"] <- "Macrophages"
lung.atlas$overall_UMAP[lung.atlas$overall_UMAP == "CD16+ Neutrophils"] <- "Neutrophils"
lung.atlas$overall_UMAP[lung.atlas$overall_UMAP == "CD16- Neutrophils"] <- "Neutrophils"

# T-cell
lung.atlas$overall_UMAP[lung.atlas$overall_UMAP == "Naive CD4+ T cells"] <- "CD4+ T cells"
lung.atlas$overall_UMAP[lung.atlas$overall_UMAP == "Th1"] <- "CD4+ T cells"
lung.atlas$overall_UMAP[lung.atlas$overall_UMAP == "T reg"] <- "CD4+ T cells"
lung.atlas$overall_UMAP[lung.atlas$overall_UMAP == "CD4+ T cells"] <- "T cell"
lung.atlas$overall_UMAP[lung.atlas$overall_UMAP == "CD8+ T cells"] <- "T cell"
lung.atlas$overall_UMAP[lung.atlas$overall_UMAP == "NK cells"] <- "T cell"

# Fibroblasts
lung.atlas$overall_UMAP[lung.atlas$overall_UMAP == "Alveolar fibroblasts"] <- "Fibroblasts"
lung.atlas$overall_UMAP[lung.atlas$overall_UMAP == "Adventitial fibroblasts"] <- "Fibroblasts"
lung.atlas$overall_UMAP[lung.atlas$overall_UMAP == "Myofibroblast"] <- "Fibroblasts"
lung.atlas$overall_UMAP[lung.atlas$overall_UMAP == "Lipofibroblast"] <- "Fibroblasts"

# Endothelium
lung.atlas$overall_UMAP[lung.atlas$overall_UMAP == "Lymphatic endothelium"] <- "Endothelium"
lung.atlas$overall_UMAP[lung.atlas$overall_UMAP == "Venous endothelium"] <- "Endothelium"

DimPlot(lung.atlas, group.by = "overall_UMAP", cols = overall_UMAP_colours, order = c("T cell", "Neutrophils", "AT1", "AT2", "Fibroblasts", "Endothelium", "Plasma cells", "Basal",
                                                                                      "Smooth muscle cells", "Ciliated cells"))

###Prepare data for plotting
circ_data <- prepare_circlize_data(lung.atlas, scale = 0.75)
set.seed(1234)
cluster_colors<- overall_UMAP_colours
group_colors<- c("#F1CDB1", "#CBDFB8", "#8FA9DF")
names(group_colors) <- c("COVID-19", "Non-Pneumonia", "Pneumonia")
infection_colors<-c("#74ADD1", "#D73027")
names(infection_colors) <- c("Positive", "Negative")

###plot UMAP and save figures
pdf('COSMIC_lung_atlas.pdf', width = 6, height = 6)
plot_circlize(circ_data,do.label = F, pt.size = 0.1, col.use = cluster_colors ,bg.color = 'white', kde2d.n = 200, repel = T, label.cex = 0, contour.levels = c(0,0),
              contour.nlevels = 0)
add_track(circ_data, group = "Group", colors = group_colors, track_num = 2, labs = "Group", srt = 180) ## can change it to one of the columns in the meta data of your seurat object
add_track(circ_data, group = "HIV_Status",colors = infection_colors, track_num = 3, labs = "HIV Status", srt = 180)
dev.off()

pdf('COSMIC_lung_immune_atlas.pdf', width = 6, height = 6)
DimPlot(lung.immune, group.by = "celltype.final", cols = immune_UMAP_colours) + NoLegend()
dev.off()

pdf('COSMIC_lung_stromal_atlas.pdf', width = 6, height = 6)
DimPlot(lung.stromal, cols = stromal_UMAP_colours) + NoLegend()
dev.off()

stromal_UMAP_colours <- c("#90A37F", "#d6bcc0", "#e07b91","#F7867C", "#b86cb9","#864E87","#B959BA", "#B098D0", "#11AD82", "#17becf", "#f0b98d",
                          "#bcbd22", "#D4674A","#c49c94", "#c7c7c7")

names(stromal_UMAP_colours) <- c("AT1", "AT2", "Venous endothelium", "Lymphatic endothelium", "Adventitial fibroblasts","Alveolar fibroblasts", "Myofibroblast", "Lipofibroblast", 
                                 "Ciliated cells", "Basal", "Mesothelial", "Smooth muscle cells",
                                 "Secretory cells", "Soup", "Ribosomal high cells")


immune_UMAP_colours2 <- c("#ff7f0e", "#ffbb78", "#DBBA27", "#9467bd", "#c5b0d5", "#165582", "#1f77b4",
                         "#2ca02c", "#aec7e8", "#98df8a", "#75AB6A","#f7b6d2", "#e377c2", "#8c564b", "#7f7f7f", "#d62728")

names(immune_UMAP_colours2) <- c("Alveolar macrophages", "Interstitial macrophages", "Monocyte-derived macrophages",
                                "CD16+ Neutrophils", "CD16- Neutrophils", "gdT cells", "CD8+ T cells", "Naive CD4+ T cells", "NK cells",
                                "Th1", "T reg", "Plasma cells", "B cells", "Mast cells", "Cycling cells",
                                "Erythrocytes")

celltype.colours <- c(stromal_UMAP_colours, immune_UMAP_colours2)



## MOVE THIS ##
### Refining T cell ####
lung_tcells <- subset(lung.atlas, idents = c("CD8+ T cells", "Naive CD4+ T cells", "NK cells",
                                             "Th1", "T reg"))


lung_tcells <- SCTransform(lung_tcells, vars.to.regress = c("percent.mt", "percent.ribo", "percent.hb", "S.Score", "G2M.Score"),
                           variable.features.n = 2000)

lung_tcells <- lung_tcells %>% 
  RunPCA(npcs = 40, assay = "SCT")
lung_tcells <- RunHarmony(lung_tcells, 
                          group.by.vars = c("Sample", "Sequencing"), 
                          reduction = "pca", assay.use = "SCT", reduction.save = "harmony")
ElbowPlot(lung_tcells, ndims = 40)

lung_tcells <- RunUMAP(lung_tcells, reduction = "harmony", dims = 1:32)
lung_tcells <- FindNeighbors(lung_tcells, reduction = "harmony", dims = 1:32, verbose=FALSE)
lung_tcells <- FindClusters(lung_tcells, verbose = FALSE, res=0.5)

DimPlot(object = lung_tcells, reduction = "umap", label = T) + theme_powerpoint()
DimPlot(object = lung_tcells, reduction = "umap", group.by="Sample") + theme_powerpoint()
DimPlot(object = lung_tcells, reduction = "umap", group.by="Sequencing") + theme_powerpoint()
DimPlot(object = lung_tcells, reduction = "umap", group.by="Group") + theme_powerpoint()
DimPlot(object = lung_tcells, reduction = "umap", group.by="orig.ident") + theme_powerpoint()

VlnPlot(lung_tcells, features = c("TRDC", "TRGC1", "TRGC2", "TRGV9", "TRDV2", "KLRD1", "IL7R", "KLRC1", "DUSP2", "GNLY", "KLRG1", "CD7", "CD4", "CD8A", "CD8B", "CD3D", "CD3E"))
DefaultAssay(lung_tcells) <- "SCT"
lung_tcells.markers <- FindAllMarkers(lung_tcells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT")
top.lung_tcells.markers <- lung_tcells.markers %>%
  group_by(cluster) %>%
  slice_max(n = 40, order_by = avg_log2FC)

gdtcell <- subset(lung_tcells, ident = 8)
gdtcell.barcodes <- rownames(gdtcell@meta.data)
lung.atlas@meta.data[gdtcell.barcodes, "celltype_annotation"] <- "gdT cells"

lung.immune@meta.data[gdtcell.barcodes, "celltype.final"] <- "gdT cells"

DimPlot(lung.immune, group.by = "celltype.final")


# Volcano plot macs ##
cosmic.covid.lrtd <- subset(lung.atlas, subset = Group == "Non-Pneumonia", invert = T)
########### Lung COVID vs Lung LRTD ############
cosmic.covid.lrtd$celltype_group <- paste(cosmic.covid.lrtd$celltype_annotation, cosmic.covid.lrtd$Group, sep = "_")
# get the amount of clusters
num_cluster=length(unique(cosmic.covid.lrtd$celltype_annotation))
Idents(cosmic.covid.lrtd) <- cosmic.covid.lrtd$celltype_annotation

cluster.ids <- unique(Idents(cosmic.covid.lrtd))
Idents(cosmic.covid.lrtd) <- "celltype_group"
DefaultAssay(de.cosmic.covid.healthy.obj) <- "RNA"

cosmic.covid.lrtd.list <- list()
cosmic.covid.lrtd.list.up <- list()
cosmic.covid.lrtd.list.down <- list()
for (x in 1:num_cluster) {
  n=cluster.ids[x]
  print(n)
  check_sanity = FALSE
  ## diff expression \
  COVID.response <- FindMarkers(cosmic.covid.lrtd,ident.1 = paste(n,"_COVID-19", sep=""), assay = "RNA", slot = "data", ident.2 = paste(n,"_Pneumonia", sep=""), verbose = T, logfc.threshold = 0, min.pct = 0.1)
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
  cosmic.covid.lrtd.list[[n]] <- COVID.response
  COVID.response1 <- subset(COVID.response, COVID.response$avg_log2FC > 0.5 & COVID.response$p_val_adj < 0.05)
  COVID.response1 <- COVID.response1[order(COVID.response1$avg_log2FC, decreasing = T),]
  COVID.response2 <- subset(COVID.response, COVID.response$avg_log2FC < -0.5 & COVID.response$p_val_adj < 0.05)
  COVID.response2 <- COVID.response2[order(COVID.response2$avg_log2FC),]
  
  cosmic.covid.lrtd.list.up[[n]] <- COVID.response1
  cosmic.covid.lrtd.list.down[[n]] <- COVID.response2
}
names(cosmic.covid.lrtd.list) <- cluster.ids
names(cosmic.covid.lrtd.list.up) <- cluster.ids
names(cosmic.covid.lrtd.list.down) <- cluster.ids


EnhancedVolcano(cosmic.covid.lrtd.list$`Alveolar macrophages`, lab = cosmic.covid.lrtd.list$`Alveolar macrophages`$gene,
                x = 'avg_log2FC', y = 'p_val_adj', xlab = bquote(~Log[2]~ 'fold change'),
                vline = 0,vlineType = "solid",
                pCutoff = 0.05, pointSize = 2.5, labSize = 4.5, 
                colAlpha = 1,legendLabSize = 15,legendPosition = 'bottom', legendIconSize = 5.0,drawConnectors = TRUE,
                widthConnectors = 0.5,legendLabels=c('Not sig.','Log (base 2) FC','p-value', 'p-value & Log (base 2) FC'),
                colConnectors = 'grey50', gridlines.major = TRUE, gridlines.minor = FALSE, border = 'partial',
                borderWidth = 0.5, borderColour = 'black')

gene.labels <- c("SPP1", "APOE", "HLA-DRB1", "CCL18", "CD74", "B2M", "SRGN", "MRC1", "S100A6", "SERF2", "IFI30", "C1QB", "C1QC", "HLA-DRA", "S100A10",
                 "TSIX",
                 "MACROD2",
                 "SFTPB",
                 "PTPRG",
                 "LDB2",
                 "NTM",
                 "ZNF385B",
                 "DLG2",
                 "PRKG1",
                 "SERPINE1",
                 "MECOM",
                 "GALNT18",
                 "LAMA2",
                 "CALD1",
                 "LSAMP")

lab_italics <- paste0("italic('", rownames(cosmic.covid.lrtd.list$`Alveolar macrophages`), "')")
EnhancedVolcano(cosmic.covid.lrtd.list$`Alveolar macrophages`,
                lab = rownames(cosmic.covid.lrtd.list$`Alveolar macrophages`),
                selectLab = gene.labels,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                FCcutoff = 0.5,
                parseLabels = TRUE,
                drawConnectors = T,
                xlim = c(-1,1))


## Vlnplots ##
library(Seurat)
library(tidyverse)

gamma_response <- read.table(file = "gamma.txt", sep = "\t", row.names = NULL)
alpha_response <- read.table(file = "alpha.txt", sep = "\t", row.names = NULL)
beta_response <- read.table(file = "beta.txt", sep = "\t", row.names = NULL)
lambda_response <- read.table(file = "lambda.txt", sep = "\t", row.names = NULL)
tnf_response <- read.table(file = "tnf.txt", sep = "\t", row.names = NULL)
gene <- rownames(lung.atlas)

gene <- rownames(de.cosmic.covid.healthy.obj)
gamma_genes <- intersect(gamma_response$V2, gene)
alpha_genes <- intersect(alpha_response$V2, gene)
tnf_genes <- intersect(tnf_response$V2, gene)
beta_genes <- intersect(beta_response$V2, gene)
lambda_genes <- intersect(lambda_response$V2, gene)
marker.list <- list()

marker.list[["Alpha"]] <- alpha_genes
marker.list[["Beta"]] <- beta_genes
marker.list[["Gamma"]] <- gamma_genes
marker.list[["Lambda"]] <- lambda_genes
marker.list[["TNF"]] <- tnf_genes

de.cosmic.covid.healthy.obj <- subset(lung.atlas, subset = Group %in% c("COVID-19", "Pneumonia"))
de.cosmic.covid.healthy.obj$celltype_source <- paste(de.cosmic.covid.healthy.obj$celltype_annotation, de.cosmic.covid.healthy.obj$Group, sep = "_")
Idents(de.cosmic.covid.healthy.obj) <- "celltype_source"
de.cosmic.covid.healthy.obj <-  AddModuleScore(de.cosmic.covid.healthy.obj ,
                                               features = marker.list,
                                               name = "FunctionScore",
                                               ctrl = 200)

for(i in 1:length(marker.list)){
  colnames(de.cosmic.covid.healthy.obj@meta.data)[colnames(de.cosmic.covid.healthy.obj@meta.data) == paste0("FunctionScore", i)] <- names(marker.list)[i]
}
MarkerNameVector <- c("Alpha", "Beta", "Gamma", "Lambda", "TNF")

Idents(de.cosmic.covid.healthy.obj) <- "celltype_annotation"
p1 <- VlnPlot(de.cosmic.covid.healthy.obj, features = "Alpha", group.by = "Group", idents = "Alveolar macrophages", y.max = 1.25, cols = c("#F1CDB1", "#8FA9DF"), pt.size = 0) + scale_y_continuous(breaks=c(-0.25, 0, 0.25, 0.5, 0.75, 1.0, 1.25), limits=c(-0.25, 1.25)) + stat_compare_means(label = "p.signif", comparisons = list(c("COVID-19", "Pneumonia")), ref.group = "Pneumonia", method = "wilcox.test") + stat_summary(fun = "median",
                                                                                                                                                                                                                                                                                                                                                                                                                                    geom = "crossbar", 
                                                                                                                                                                                                                                                                                                                                                                                                                                    width = 0.5,
                                                                                                                                                                                                                                                                                                                                                                                                                                    colour = "black")
p2 <- VlnPlot(de.cosmic.covid.healthy.obj, features = "Beta", group.by = "Group", idents = "Alveolar macrophages", y.max = 1.25, cols = c("#F1CDB1", "#8FA9DF"), pt.size = 0)+ scale_y_continuous(breaks=c(-0.25, 0, 0.25, 0.5, 0.75, 1.0, 1.25), limits=c(-0.25, 1.25)) + stat_compare_means(label = "p.signif", comparisons = list(c("COVID-19", "Pneumonia")),ref.group = "Pneumonia", method = "wilcox.test")+ stat_summary(fun = "median",
                                                                                                                                                                                                                                                                                                                                                                                                                                geom = "crossbar", 
                                                                                                                                                                                                                                                                                                                                                                                                                                width = 0.5,
                                                                                                                                                                                                                                                                                                                                                                                                                                colour = "black")
p3 <- VlnPlot(de.cosmic.covid.healthy.obj, features = "Lambda", group.by = "Group", idents = "Alveolar macrophages", y.max = 1.25, cols = c("#F1CDB1", "#8FA9DF"), pt.size = 0)+ scale_y_continuous(breaks=c(-0.25, 0, 0.25, 0.5, 0.75, 1.0, 1.25), limits=c(-0.25, 1.25)) + stat_compare_means(label = "p.signif", comparisons = list(c("COVID-19", "Pneumonia")), ref.group = "Pneumonia",method = "wilcox.test")+ stat_summary(fun = "median",
                                                                                                                                                                                                                                                                                                                                                                                                                                  geom = "crossbar", 
                                                                                                                                                                                                                                                                                                                                                                                                                                  width = 0.5,
                                                                                                                                                                                                                                                                                                                                                                                                                                  colour = "black")
p4 <- VlnPlot(de.cosmic.covid.healthy.obj, features = "Gamma", group.by = "Group", idents = "Alveolar macrophages", y.max = 1.25, cols = c("#F1CDB1", "#8FA9DF"), pt.size = 0)+ scale_y_continuous(breaks=c(-0.25, 0, 0.25, 0.5, 0.75, 1.0, 1.25), limits=c(-0.25, 1.25)) + stat_compare_means(label = "p.signif", comparisons = list(c("COVID-19", "Pneumonia")),ref.group = "Pneumonia", method = "wilcox.test")+ stat_summary(fun = "median",
                                                                                                                                                                                                                                                                                                                                                                                                                                 geom = "crossbar", 
                                                                                                                                                                                                                                                                                                                                                                                                                                 width = 0.5,
                                                                                                                                                                                                                                                                                                                                                                                                                                 colour = "black")
p5 <- VlnPlot(de.cosmic.covid.healthy.obj, features = "TNF", group.by = "Group", idents = "Alveolar macrophages", y.max = 1.25,  cols = c("#F1CDB1", "#8FA9DF"), pt.size = 0)+scale_y_continuous(breaks=c(-0.25, 0, 0.25, 0.5, 0.75, 1.0, 1.25), limits=c(-0.25, 1.25))+ stat_compare_means(label = "p.signif", comparisons = list(c("COVID-19", "Pneumonia")), ref.group = "Pneumonia",method = "wilcox.test")+ stat_summary(fun = "median",
                                                                                                                                                                                                                                                                                                                                                                                                                              geom = "crossbar", 
                                                                                                                                                                                                                                                                                                                                                                                                                              width = 0.5,
                                                                                                                                                                                                                                                                                                                                                                                                                              colour = "black")
plot.list <- list(p1, p2, p3, p4, p5)
cowplot::plot_grid(plotlist = plot.list, nrow = 1, align = "h", axis = "bt")
