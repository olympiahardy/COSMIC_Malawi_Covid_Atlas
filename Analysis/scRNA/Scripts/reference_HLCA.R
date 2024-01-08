library(SeuratDisk)
library(Seurat)
library(tidyverse)
library(harmony)

Convert("/datastore3/Olympia/HLCA_ref/lung_final_reference_test.h5ad", dest = "h5seurat", overwrite = TRUE)
lung_atlas <- LoadH5Seurat("/datastore3/Olympia/HLCA_ref/lung_final_reference_test.h5seurat", meta.data = F, misc = FALSE)

lung_atlas_meta_subset <- read.csv("/datastore3/Olympia/HLCA_ref/lung_final_reference_meta.csv", row.names = 1)

write10xCounts("/datastore3/Olympia/HLCA_ref/raw/", lung_atlas@assays$RNA@counts, version = "3", overwrite = T)

lung_atlas <- AddMetaData(lung_atlas, lung_atlas_meta_subset)

cosmic_lung <- readRDS("/datastore3/Olympia/COSMIC_data/Final_lung_objects/lung_all_12923.rds")

Idents(cosmic_lung) <- cosmic_lung$celltype_annotation
cosmic_celltypes <- c("Adventitial fibroblasts" , "Alveolar fibroblasts", "Alveolar macrophages", "AT1", "AT2",
                      "CD16- Neutrophils","CD16+ Neutrophils", "CD8+ T cells", "Interstitial macrophages", "Monocyte-derived macrophages", "Naive CD4+ T cells","NK cells",
                      "T reg", "Th1", "Venous endothelium")
cosmic_lung_subset <- subset(cosmic_lung, idents = cosmic_celltypes)

DefaultAssay(cosmic_lung_subset) <- "RNA"
cosmic_lung_subset$source <- "COSMIC"
lung_atlas$source <- "HLCA"




library(org.Hs.eg.db)
library(biomaRt)

matrix <- Read10X("/datastore3/Olympia/HLCA_ref/raw/")
#Copy code to load in C40 expression matrix and meta data here

mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
x <- getBM(attributes =c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", row.names(matrix), mart = mart)

#There are some ensembl IDs that will not have an equivalent hgnc symbol, so we have to remove them
x <- subset(x, hgnc_symbol != "")

#Some ensembl IDs can map to multiple hgnc symbols (and vice versa), so we remove any that have two or more matches
#You can also just pick one of the genes at random
non_dup_genes <- names(table(x$hgnc_symbol))[ which(table(x$hgnc_symbol) == 1 ) ]
x <- subset(x, hgnc_symbol %in% non_dup_genes)

x <- x[order(match(x$ensembl_gene_id, row.names(matrix))),]

#Subset out the genes that aren't in the list of transformed gene IDs
matrix <- matrix[ x$ensembl_gene_id,]

#Add the new gene names to the expression matrix
row.names(matrix) <- x$hgnc_symbol

lung_atlas <- CreateSeuratObject(counts = matrix , meta.data = lung_atlas_meta_subset, project = "HLCA")

lung_atlas[["percent.ribo"]] <- PercentageFeatureSet(lung_atlas, pattern = "^RP[SL]")
lung_atlas[["percent.hb"]] <- PercentageFeatureSet(lung_atlas, pattern = "^HB[^(P)]")
lung_atlas[["percent.mt"]] <- PercentageFeatureSet(lung_atlas, pattern = "^MT-")
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
lung_atlas <- lung_atlas %>%
  NormalizeData() %>% 
  CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes)

lung.list <- c(lung_atlas, cosmic_lung_subset)

min_median_umi <- min(unlist(lapply(lung.list, FUN = function(x) median(x$nCount_RNA))))
lung.list <- lapply(X = lung.list, FUN = function(x) {
  x <- SCTransform(x, vars.to.regress = c("percent.mt", "percent.ribo", "percent.hb", "S.Score", "G2M.Score"), vst.flavor = "v2", scale_factor=min_median_umi)
})

features <- SelectIntegrationFeatures(object.list = lung.list, nfeatures = 3000, assay = c("SCT", "SCT"))

merged_lung <- merge(x = lung.list[[1]],
                       y = lung.list[2:length(lung.list)],
                       add.cell.ids = c("HLCA", "COSMIC"),
                       merge.data = T)

merged_lung <- merged_lung %>% 
  RunPCA(npcs = 40, assay = "SCT", features = features)

merged_lung$disease[merged_lung$Group == "COVID-19"] <- "COVID-19"
merged_lung$disease[merged_lung$Group == "Pneumonia"] <- "pneumonia"
merged_lung$disease[merged_lung$Group == "Non-Pneumonia"] <- "normal"

merged_lung <- RunHarmony(merged_lung, 
                             group.by.vars = c("source", "disease"), 
                             reduction = "pca", assay.use = "SCT", reduction.save = "harmony")
ElbowPlot(merged_lung, ndims = 40)

merged_lung <- RunUMAP(merged_lung, reduction = "harmony", dims = 1:38)
merged_lung <- FindNeighbors(merged_lung, reduction = "harmony", dims = 1:38, verbose=FALSE)
merged_lung <- FindClusters(merged_lung, verbose = FALSE, res=0.2)

p1 <- DimPlot(object = merged_lung, reduction = "umap", label = T) + theme_powerpoint()
p2 <- DimPlot(object = merged_lung, reduction = "umap", group.by = "source") + theme_powerpoint()
p3 <- DimPlot(object = merged_lung, reduction = "umap", group.by = "disease") + theme_powerpoint()
p1+p2+p3

DimPlot(object = merged_lung, reduction = "umap", split.by = "source") + theme_powerpoint()

DefaultAssay(merged_lung) <- "SCT"
merged_lung.markers <- FindAllMarkers(merged_lung, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT")
top.merged_lung.markers <- merged_lung.markers %>%
  group_by(cluster) %>%
  slice_max(n = 30, order_by = avg_log2FC)

features <- SelectIntegrationFeatures(object.list = lung.list, nfeatures = 3000, assay = c("SCT", "SCT"))



write.csv(merged_lung.markers, "/datastore3/Olympia/HLCA_ref/integrated_HLCA_markers.csv")
saveRDS(merged_lung, "/datastore3/Olympia/HLCA_ref/integrated_HLCA_COSMIC.rds")

merged_lung <- subset(merged_lung, idents = 21, invert = T)

merged_lung <- RenameIdents(merged_lung, "0" = "T cells",
             "1" = "Neutrophils",
             "2" = "Alveolar macrophages",
             "3" = "Endothelium",
             "4" = "AT1",
             "5" = "Endothelium",
             "6" = "AT2",
             "7" = "AT2",
             "8" = "Alveolar Fibroblasts",
             "9" = "Monocyte derived marcophages",
             "10" = "T reg",
             "11" = "Adventitial Fibroblasts",
             "12" = "Endothelium",
             "13" = "Interstitial macrophages",
             "14" = "AT1",
             "15" = "Endothelium",
             "16" = "Interstitial macrophages",
             "17" = "Endothelium",
             "18" = "Fibroblasts",
             "19" = "T reg",
             "20" = "Neutrophils")



merged_lung_tcells <- subset(merged_lung, idents = "T cells")


merged_lung_tcells <- SCTransform(merged_lung_tcells, vars.to.regress = c("percent.mt", "percent.ribo", "percent.hb", "S.Score", "G2M.Score"))


merged_lung_tcells <- merged_lung_tcells %>% 
  RunPCA(npcs = 40, assay = "SCT")

merged_lung_tcells <- RunHarmony(merged_lung_tcells, 
                          group.by.vars = c("source", "disease"), 
                          reduction = "pca", assay.use = "SCT", reduction.save = "harmony")

ElbowPlot(merged_lung_tcells, ndims = 40)

merged_lung_tcells <- RunUMAP(merged_lung_tcells, reduction = "harmony", dims = 1:26)
merged_lung_tcells <- FindNeighbors(merged_lung_tcells, reduction = "harmony", dims = 1:26, verbose=FALSE)
merged_lung_tcells <- FindClusters(merged_lung_tcells, verbose = FALSE, res=0.4)

p1 <- DimPlot(object = merged_lung_tcells, reduction = "umap", label = T) + theme_powerpoint()
p2 <- DimPlot(object = merged_lung_tcells, reduction = "umap", group.by = "source") + theme_powerpoint()
p3 <- DimPlot(object = merged_lung_tcells, reduction = "umap", group.by = "disease") + theme_powerpoint()
p1+p2+p3

VlnPlot(merged_lung_tcells, features = c("CD8A", "CD8B", "IL7R", "KLRD1", "CD3D", "CD3E", "CD4"))

DefaultAssay(merged_lung_tcells) <- "SCT"
merged_lung_tcells.markers <- FindAllMarkers(merged_lung_tcells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT")
top.merged_lung_tcells.markers <- merged_lung_tcells.markers %>%
  group_by(cluster) %>%
  slice_max(n = 30, order_by = avg_log2FC)

merged_lung_tcells <- RenameIdents(merged_lung_tcells, "0" = "NK cells",
                            "1" = "CD4 T cells",
                            "2" = "CD8 T cells",
                            "3" = "CD4 T cells",
                            "4" = "CD4 T cells",
                            "5" = "CD4 T cells",
                            "6" = "NK cells",
                            "7" = "CD8 T cells",
                            "8" = "CD8 T cells",
                            "9" = "CD4 T cells",
                            "10" = "CD4 T cells",
                            "11" = "CD4 T cells",
                            "12" = "CD4 T cells",
                            "13" = "NK cells",
                            "14" = "CD4 T cells",
                            "15" = "CD8 T cells",
                            "16" = "CD4 T cells")

merged_lung_tcells$cell_type <- Idents(merged_lung_tcells)
t_cells_annotation <- as.data.frame(merged_lung_tcells$cell_type)

merged_lung$cell.type <- Idents(merged_lung)
all.meta <- merged_lung@meta.data
all.meta$barcode <- rownames(all.meta)
t_cells_annotation <- rownames_to_column(t_cells_annotation, "barcode")
names(t_cells_annotation) <- c("barcode", "cell.type")
res.all <- all.meta %>% left_join(t_cells_annotation, by = c("barcode")) %>%
  mutate(celltype=coalesce(cell.type,cell.type))
lung_all$celltype_annotation <- res.all$cell.type.final

res <- all.meta %>% left_join(t_cells_annotation, by = c("barcode")) %>%
  mutate(celltype=coalesce(cell.type.y,cell.type.x))

merged_lung$final_celltype <- res$celltype

Idents(merged_lung) <- merged_lung$final_celltype




