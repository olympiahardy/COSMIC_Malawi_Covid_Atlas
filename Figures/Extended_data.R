library(Seurat)
library(tidyverse)
library(data.table)
library(dendextend)
library(ComplexHeatmap)
library(ggplot2)
library(ggrepel)

## Extended data Figure 3 ##

covid.infection <- subset(lung.atlas, subset = Covid > 0)
covid.infection.nasal <- subset(nasal.atlas, subset = Covid > 0)
covid.infection.blood <- subset(blood.atlas, subset = Covid > 0)

sample.colours <- c("#F8766D", "#DE8C00", "#00BA38", "#00B4F0", "#619CFF", "#C77CFF", "#F564E3")
names(sample.colours) <- c("Cos-1", "Cos-6", "Cos-11", "Cos-12", "Cos-14", "Cos-15", "Cos-16")
DimPlot(covid.infection, group.by = "Sample", cols = sample.colours)
DimPlot(covid.infection.nasal, group.by = "Sample", cols = sample.colours)
DimPlot(covid.infection.blood, group.by = "Sample", cols = sample.colours)

## Extended data Figure 4 ##

## Lung cell type proportions ##
# Plotting cell type proportions #
Idents(lung.atlas) <- lung.atlas$celltype_annotation
lung.immune <- subset(lung.atlas, idents = c("Alveolar macrophages", "Interstitial macrophages", "Monocyte-derived macrophages",
                                             "CD16+ Neutrophils", "CD16- Neutrophils", "CD8+ T cells", "NK cells", "gdT cells", "Naive CD4+ T cells",
                                             "Th1", "T reg", "Plasma cells", "B cells", "Mast cells", "Ribosomal high cells", "Cycling cells","Erythrocytes"))
freq.table.group <- as.data.frame(prop.table(x = table(Idents(lung.immune), lung.immune$Group), margin = 2))
freq.table.group$Var1 <- as.ordered(factor(freq.table.group$Var1,
                                           levels=c("Alveolar macrophages", "Interstitial macrophages", "Monocyte-derived macrophages",
                                                    "CD16+ Neutrophils", "CD16- Neutrophils", "CD8+ T cells", "NK cells", "gdT cells", "Naive CD4+ T cells",
                                                    "Th1", "T reg", "Plasma cells", "B cells", "Mast cells", "Ribosomal high cells", "Cycling cells","Erythrocytes")))

freq.table.group$Var2 <- as.ordered(factor(freq.table.group$Var2,
                                           levels=rev(c("Non-Pneumonia", "Pneumonia", "COVID-19"))))
ggplot(data=freq.table.group, aes(x=freq.table.group$Var2, y = freq.table.group$Freq, fill=freq.table.group$Var1)) + geom_bar(stat="identity", color="black") +
  labs(x="Group", y="Proportion of cells", fill="Cell Type") + scale_x_discrete(limits = rev(levels(freq.table.group$Var2))) + scale_fill_manual(values=celltype.colours) + theme_minimal()


lung.stromal <- subset(lung.atlas, idents = c("AT1", "AT2", "Venous endothelium", "Lymphatic endothelium", 
                                              "Adventitial fibroblasts","Alveolar fibroblasts", "Myofibroblast", "Lipofibroblast", "Ciliated cells", "Basal", "Mesothelial",
                                              "Smooth muscle cells","Secretory cells", "Soup", "Ribosomal high cells"))
freq.table.group <- as.data.frame(prop.table(x = table(Idents(lung.stromal), lung.stromal$Group), margin = 2))
freq.table.group$Var1 <- as.ordered(factor(freq.table.group$Var1,
                                           levels=c("AT1", "AT2", "Venous endothelium", "Lymphatic endothelium", 
                                                    "Adventitial fibroblasts","Alveolar fibroblasts", "Myofibroblast", "Lipofibroblast", "Ciliated cells", "Basal", "Mesothelial",
                                                    "Smooth muscle cells","Secretory cells", "Soup", "Ribosomal high cells")))

freq.table.group$Var2 <- as.ordered(factor(freq.table.group$Var2,
                                           levels=rev(c("Non-Pneumonia", "Pneumonia", "COVID-19"))))
ggplot(data=freq.table.group, aes(x=freq.table.group$Var2, y = freq.table.group$Freq, fill=freq.table.group$Var1)) + geom_bar(stat="identity", color="black") +
  labs(x="Group", y="Proportion of cells", fill="Cell Type") + scale_x_discrete(limits = rev(levels(freq.table.group$Var2))) + scale_fill_manual(values=celltype.colours) + theme_minimal()


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

de.cosmic.covid.healthy.obj <- subset(de.cosmic.covid.healthy.obj, subset = celltype_annotation %in% c("Alveolar macrophages", "Interstitial macrophages", "Monocyte-derived macrophages",
                                                                                                       "CD16+ Neutrophils", "CD16- Neutrophils", "CD8+ T cells", "Naive CD4+ T cells", "NK cells", "gdT cells",
                                                                                                       "Th1", "T reg", "Plasma cells", "B cells", "Mast cells", "AT1", "AT2"))
de.cosmic.covid.healthy.obj$celltype_annotation <- as.ordered(factor(de.cosmic.covid.healthy.obj$celltype_annotation, levels = c("Alveolar macrophages", "Interstitial macrophages", "Monocyte-derived macrophages",
                                                                                                                                 "CD16+ Neutrophils", "CD16- Neutrophils", "CD8+ T cells", "Naive CD4+ T cells", "NK cells", "gdT cells",
                                                                                                                                 "Th1", "T reg", "Plasma cells", "B cells", "Mast cells", "AT1", "AT2")))
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

FunctionScoreMatrix <- matrix(0,
                              ncol = length(unique(Idents(de.cosmic.covid.healthy.obj))),
                              nrow = length(marker.list))

colnames(FunctionScoreMatrix) <- unique(Idents(de.cosmic.covid.healthy.obj))
rownames(FunctionScoreMatrix) <- MarkerNameVector


for(ci in 1:ncol(FunctionScoreMatrix)){
  for(ri in 1:nrow(FunctionScoreMatrix)){
    FunctionVec <- as_tibble(de.cosmic.covid.healthy.obj@meta.data) %>% pull(MarkerNameVector[ri])
    fv <- mean(FunctionVec[de.cosmic.covid.healthy.obj$celltype_source == unique(Idents(de.cosmic.covid.healthy.obj))[ci]])
    FunctionScoreMatrix[ri, ci] <- fv
  }
}

cells <- gsub("\\_.*", "", colnames(FunctionScoreMatrix))    
annotation_col <- data.frame(Cluster = cells)
rownames(annotation_col) <- colnames(FunctionScoreMatrix)
annotation_col$Group <- gsub(".*\\_", "", rownames(annotation_col))    
annotation_col$Cluster <- as.ordered(factor(annotation_col$Cluster, levels = c("Alveolar macrophages", "Interstitial macrophages", "Monocyte-derived macrophages",
                                                                               "CD16+ Neutrophils", "CD16- Neutrophils", "CD8+ T cells", "Naive CD4+ T cells", "NK cells", "gdT cells",
                                                                               "Th1", "T reg", "Plasma cells", "B cells", "Mast cells", "AT1", "AT2")))
annotation_col$Group <- as.ordered(factor(annotation_col$Group,
                                          levels=c("COVID-19", "Pneumonia")))
annotation_col <- annotation_col[order(annotation_col$Cluster, annotation_col$Group),]


FunctionScoreMatrix <-
  FunctionScoreMatrix[, rownames(annotation_col)]

FunctionScoreMatrix <- t(apply(FunctionScoreMatrix, 1, rescale, to=c(-1, 1)))

my.breaks <- c(seq(-1, 0, by=0.1), seq(0.1, 1, by=0.1))
my.colors <- c(
  colorRampPalette(colors = c("#2166AC", "white"))(length(my.breaks)/2),
  colorRampPalette(colors = c("white", "#B2182B"))(length(my.breaks)/2))
## cellType_col <- data.frame(cell.type = CD8_Obj_CellType)
## rownames(cellType_col) <- colnames(FunctionScoreMatrix)

#Create colour scheme
cluster.cols <- c("#90A37F",
                  "#d6bcc0",
                  "#e07b91",
                  "#F7867C",
                  "#b86cb9",
                  "#864E87",
                  "#B959BA",
                  "#B098D0",
                  "#11AD82",
                  "#17becf",
                  "#f0b98d",
                  "#bcbd22",
                  "#D4674A",
                  "#c49c94",
                  "#c7c7c7",
                  "#ff7f0e")
names(cluster.cols) <- c("Alveolar macrophages", "Interstitial macrophages", "Monocyte-derived macrophages",
                         "CD16+ Neutrophils", "CD16- Neutrophils", "CD8+ T cells", "Naive CD4+ T cells", "NK cells", "gdT cells",
                         "Th1", "T reg", "Plasma cells", "B cells", "Mast cells", "AT1", "AT2")
group.cols <- c("#F1CDB1", "#8FA9DF")
names(group.cols) <-c("COVID-19", "Pneumonia")
signature.cols <- c("#8DD3C7", "#BEBADA", "#FCCDE5", "#6DCCFD", "blue")
names(signature.cols) <- c("Alpha Response", "Beta response", "Gamma Response", "Lambda Response", "TNF response")
metacols <- list(Cluster = cluster.cols,
                 Group = group.cols,
                 Signature = signature.cols)

signatureType_row <- data.frame(Signature.type = c(
  rep("Alpha Response", 1),
  rep("Beta response",1),
  rep("Gamma Response", 1),
  rep("Lambda response",1),
  rep("TNF response", 1)))

rownames(signatureType_row) <- MarkerNameVector
pheatmap(FunctionScoreMatrix,
         show_colnames = T,
         show_rownames = T,
         annotation_row = signatureType_row,
         annotation_col = annotation_col,
         annotation_colors = metacols,
         cluster_rows = F,
         cluster_cols = F,
         breaks = my.breaks,
         color = my.colors,
         border_color = "NA",
         fontsize = 8,
         width = 5,
         height = 3.8,
         gaps_col=head(as.numeric(cumsum(table(annotation_col$Cluster))), -1))

## Extended data Figure 5 ##
## Interferon 

interferon.list1 <- c("IFNL1",
                      "IFNL2",
                      "IFNL3",
                      "IL24",
                      "IFNLR1",
                      "IFNAR1",
                      "IFNAR2",
                      "CCL11",
                      "CCL2",
                      "IL6",
                      "IL12B",
                      "TNF",
                      "CSF2",
                      "IL2",
                      "IL1B",
                      "IL17A",
                      "CXCL8",
                      "IL1RN",
                      "CCL5",
                      "IL13",
                      "IL7",
                      "CXCL1",
                      "CXCL12",
                      "CCL4",
                      "IL15",
                      "IL1A",
                      "IL4",
                      "IL10",
                      "IFNA1",
                      "IL21",
                      "IL27",
                      "IFNG",
                      "IL18")
interferon.list2 <- c("CCRL2",
                      "CSF1",
                      "CXCL10",
                      "CXCL11",
                      "HLA-C",
                      "IFI27",
                      "IFI30",
                      "IFI35",
                      "IFI44",
                      "IFI44L",
                      "IFIH1",
                      "IFIT2",
                      "IFIT3",
                      "IFITM1",
                      "IFITM2",
                      "IFITM3",
                      "IL15",
                      "IL4R",
                      "IL7",
                      "IRF1",
                      "IRF2",
                      "IRF7",
                      "IRF9",
                      "ISG15",
                      "ISG20",
                      "LAMP3",
                      "LAP3",
                      "LGALS3BP",
                      "OAS1",
                      "OASL",
                      "OGFR",
                      "IL1B",
                      "IL17A",
                      "CCL11",
                      "CCL2",
                      "IL6",
                      "IL2",
                      "CXCL1",
                      "IFNG",
                      "IL10",
                      "IL13",
                      "IL1R1",
                      "TNF",
                      "TNFRSF1A",
                      "IL18",
                      "IL27",
                      "IL21",
                      "IL4",
                      "IL1A",
                      "CCL4",
                      "CXCL12")
interferon.list <- c(interferon.list1, interferon.list2)
interferon.list <- unique(interferon.list)


de.cosmic.covid.healthy.obj <- subset(merged_lung, subset = disease == "COVID-19")
de.cosmic.covid.healthy.obj$sourcegroup <- paste(de.cosmic.covid.healthy.obj$disease, de.cosmic.covid.healthy.obj$source, sep = "_")
Idents(de.cosmic.covid.healthy.obj) <- de.cosmic.covid.healthy.obj$final_celltype
########### COSMIC COVID vs HLCA COVID ############
de.cosmic.covid.healthy.obj$celltype_source <- paste(de.cosmic.covid.healthy.obj$final_celltype, de.cosmic.covid.healthy.obj$source, sep = "_")
# get the amount of clusters
num_cluster=length(unique(de.cosmic.covid.healthy.obj$final_celltype))
Idents(de.cosmic.covid.healthy.obj) <- de.cosmic.covid.healthy.obj$final_celltype
cluster.ids <- unique(Idents(de.cosmic.covid.healthy.obj))
Idents(de.cosmic.covid.healthy.obj) <- "celltype_source"
DefaultAssay(de.cosmic.covid.healthy.obj) <- "SCT"
for (x in 1:num_cluster) {
  n=cluster.ids[x]
  check_sanity = FALSE
  ## diff expression \
  COVID.response <- FindMarkers(de.cosmic.covid.healthy.obj,ident.1 = paste(n,"_COSMIC", sep=""), ident.2 = paste(n,"_HLCA", sep=""), features = interferon.list, verbose = FALSE, logfc.threshold = 0.25, min.pct = 0.1)#, latent.vars = c("Participant"))
  COVID.response$gene <- rownames(COVID.response)
  variable=paste(n)
  #COVID.response1 <- subset(COVID.response, COVID.response$avg_log2FC > 0 & COVID.response$p_val_adj < 0.01)
  #COVID.response1 <- COVID.response1[order(COVID.response1$avg_log2FC, decreasing = T),]
  #COVID.response2 <- subset(COVID.response, COVID.response$avg_log2FC < 0 & COVID.response$p_val_adj < 0.01)
  #COVID.response2 <- COVID.response2[order(COVID.response2$avg_log2FC),]
  all.interferon.gene.list[[variable]] <- COVID.response
  
}


de.cosmic.covid.healthy.obj <- subset(merged_lung, subset = source == "COSMIC", invert = T)
de.cosmic.covid.healthy.obj <- subset(de.cosmic.covid.healthy.obj, subset = disease == "pneumonia", invert = T)
Idents(de.cosmic.covid.healthy.obj) <- de.cosmic.covid.healthy.obj$final_celltype
########### HLCA COVID vs HLCA Healthy ############
de.cosmic.covid.healthy.obj$celltype_group <- paste(de.cosmic.covid.healthy.obj$final_celltype, de.cosmic.covid.healthy.obj$disease, sep = "_")
# get the amount of clusters
num_cluster=length(unique(de.cosmic.covid.healthy.obj$final_celltype))
Idents(de.cosmic.covid.healthy.obj) <- de.cosmic.covid.healthy.obj$final_celltype
cluster.ids <- unique(Idents(de.cosmic.covid.healthy.obj))
Idents(de.cosmic.covid.healthy.obj) <- "celltype_group"
DefaultAssay(de.cosmic.covid.healthy.obj) <- "SCT"

all.interferon.gene.list <- list()
for (x in 1:num_cluster) {
  n=cluster.ids[x]
  check_sanity = FALSE
  ## diff expression \
  COVID.response <- FindMarkers(de.cosmic.covid.healthy.obj,ident.1 = paste(n,"_COVID-19", sep=""), ident.2 = paste(n,"_normal", sep=""), features = interferon.list, verbose = FALSE, logfc.threshold = 0.25, min.pct = 0.1)#, latent.vars = c("Participant"))
  COVID.response$gene <- rownames(COVID.response)
  variable=paste(n)
  #COVID.response1 <- subset(COVID.response, COVID.response$avg_log2FC > 0 & COVID.response$p_val_adj < 0.01)
  #COVID.response1 <- COVID.response1[order(COVID.response1$avg_log2FC, decreasing = T),]
  #COVID.response2 <- subset(COVID.response, COVID.response$avg_log2FC < 0 & COVID.response$p_val_adj < 0.01)
  #COVID.response2 <- COVID.response2[order(COVID.response2$avg_log2FC),]
  all.interferon.gene.list[[variable]] <- COVID.response
  
}

de.cosmic.covid.healthy.obj <- subset(merged_lung, subset = disease == "pneumonia", invert = T)
de.cosmic.covid.healthy.obj$sourcegroup <- paste(de.cosmic.covid.healthy.obj$disease, de.cosmic.covid.healthy.obj$source, sep = "_")
de.cosmic.covid.healthy.obj <- subset(de.cosmic.covid.healthy.obj, subset = sourcegroup %in% c("COVID-19_HLCA", "normal_COSMIC"), invert = T)
Idents(de.cosmic.covid.healthy.obj) <- de.cosmic.covid.healthy.obj$final_celltype
########### Lung COVID vs HLCA Healthy ############
de.cosmic.covid.healthy.obj$celltype_source <- paste(de.cosmic.covid.healthy.obj$final_celltype, de.cosmic.covid.healthy.obj$source, sep = "_")
# get the amount of clusters
num_cluster=length(unique(de.cosmic.covid.healthy.obj$final_celltype))
Idents(de.cosmic.covid.healthy.obj) <- de.cosmic.covid.healthy.obj$final_celltype
cluster.ids <- unique(Idents(de.cosmic.covid.healthy.obj))
Idents(de.cosmic.covid.healthy.obj) <- "celltype_source"
DefaultAssay(de.cosmic.covid.healthy.obj) <- "SCT"

all.interferon.gene.list <- list()
for (x in 1:num_cluster) {
  n=cluster.ids[x]
  check_sanity = FALSE
  ## diff expression \
  COVID.response <- FindMarkers(de.cosmic.covid.healthy.obj,ident.1 = paste(n,"_COSMIC", sep=""), ident.2 = paste(n,"_HLCA", sep=""), features = interferon.list, verbose = FALSE, logfc.threshold = 0.25, min.pct = 0.1)#, latent.vars = c("Participant"))
  COVID.response$gene <- rownames(COVID.response)
  variable=paste(n)
  #COVID.response1 <- subset(COVID.response, COVID.response$avg_log2FC > 0 & COVID.response$p_val_adj < 0.01)
  #COVID.response1 <- COVID.response1[order(COVID.response1$avg_log2FC, decreasing = T),]
  #COVID.response2 <- subset(COVID.response, COVID.response$avg_log2FC < 0 & COVID.response$p_val_adj < 0.01)
  #COVID.response2 <- COVID.response2[order(COVID.response2$avg_log2FC),]
  all.interferon.gene.list[[variable]] <- COVID.response
  
}

## Extended data Figure 6 ##
luminex.panel <- c("IFNG",
                   "TNF",
                   "CXCL8",
                   "IL18",
                   "IL10",
                   "CXCL1",
                   "CCL4",
                   "CCL11",
                   "IL1RN",
                   "IL1B",
                   "IL21",
                   "IL12A",
                   "CSF2",
                   "IL2",
                   "IL15",
                   "IL6",
                   "CCL2",
                   "CCL5",
                   "IL27",
                   "IL4",
                   "IL13",
                   "IL1A",
                   "CXCL12",
                   "IL17A",
                   "IL7",
                   "IFNA1")

# Lung atlas
Idents(lung.atlas) <- "celltype_annotation"
lung.atlas.luminex <- subset(lung.atlas, idents = c("Alveolar macrophages", "Interstitial macrophages", "Monocyte-derived macrophages",
                                                    "CD16+ Neutrophils", "CD16- Neutrophils", "CD8+ T cells", "Naive CD4+ T cells", "NK cells", "gdT cells",
                                                    "Th1", "T reg", "Plasma cells", "B cells", "Mast cells"))

lung.atlas.luminex$celltype_annotation <- factor(lung.atlas.luminex$celltype_annotation)
cluster.ids <- unique(levels(lung.atlas.luminex$celltype_annotation))
lung.atlas.luminex$pseudo <- "Cluster"
Idents(lung.atlas.luminex) <- "pseudo"
clusterGroupBulk <- AverageExpression(lung.atlas.luminex, return.seurat = TRUE, add.ident = "Sample")

#Fetch expression of genes
DefaultAssay(clusterGroupBulk) <- "RNA"
matrix <- FetchData(clusterGroupBulk, vars = luminex.panel)
data <- t(matrix)

# Blood atlas
blood.atlas.luminex <- subset(blood.atlas, idents = "Reticulocytes", invert = T)
cluster.ids <- unique(levels(blood.atlas.luminex$celltype))
blood.atlas.luminex$pseudo <- "Cluster"
Idents(blood.atlas.luminex) <- "pseudo"
clusterGroupBulk <- AverageExpression(blood.atlas.luminex, return.seurat = TRUE, add.ident = "Sample")

#Fetch expression of genes
DefaultAssay(clusterGroupBulk) <- "RNA"
matrix <- FetchData(clusterGroupBulk, vars = luminex.panel)
data <- t(matrix)

# Nasal atlas
cluster.ids <- unique(levels(nasal.atlas$celltype_final))
nasal.atlas$pseudo <- "Cluster"
Idents(nasal.atlas) <- "pseudo"
clusterGroupBulk <- AverageExpression(nasal.atlas, return.seurat = TRUE, add.ident = "Sample")

#Fetch expression of genes
DefaultAssay(clusterGroupBulk) <- "RNA"
matrix <- FetchData(clusterGroupBulk, vars = luminex.panel)
data <- t(matrix)