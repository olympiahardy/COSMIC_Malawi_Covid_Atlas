###Prepare data for plotting - merged atlas
circ_data <- prepare_circlize_data(merged_lung, scale = 0.75)
set.seed(1234)
cluster_colors<- merged_UMAP_colours
disease_colors<- c("#F1CDB1", "#CBDFB8", "#8FA9DF")
names(disease_colors) <- c("COVID-19", "normal", "pneumonia")
source_colors<-c("#DB6C0B", "#4198AA")
names(source_colors) <- c("COSMIC", "HLCA")

###plot and save figures
pdf('circlize_plot_lung_cleaned.pdf', width = 6, height = 6)
plot_circlize(circ_data,do.label = F, pt.size = 0.1, col.use = merged_UMAP_colours ,bg.color = 'white', kde2d.n = 200, repel = T, label.cex = 0, contour.levels = c(0,0),
              contour.nlevels = 0)
add_track(circ_data, group = "disease", colors = disease_colors, track_num = 2, labs = "Disease", srt = 180) ## can change it to one of the columns in the meta data of your seurat object
add_track(circ_data, group = "source",colors = source_colors, track_num = 3, labs = "Source", srt = 180)
dev.off()

## GO ##

library(fgsea)

makeGeneList <- function(gl) {
  gl <- gl %>% dplyr::select(gene, avg_log2FC , p_val)
  gl$neglog10pval <- -log10(gl$p_val)
  rank <- unlist(gl$neglog10pval*sign(gl$avg_log2FC))
  rank[rank == Inf] = 300
  rank[rank == -Inf] = -300
  names(rank) <- gl$gene
  rank <- rev(sort(rank))
  gl <- rank
}


makeGeneList2 <- function(gl) {
  ranks <- gl$avg_log2FC
  names(ranks) <- gl$gene
  gl <- ranks
}

celltypes <- c("NK cells", 
               "CD4 T cells", 
               "CD8 T cells",
               "T reg",
               "Neutrophils",
               "Alveolar macrophages",
               "Interstitial macrophages",
               "Monocyte derived marcophages",
               "AT1",
               "AT2",
               "Endothelium",
               "Alveolar Fibroblasts",
               "Adventitial Fibroblasts",
               "Fibroblasts")

comparisons <- list()
comparisons$COSMIC <- list()
comparisons$HLCA <- list()

for (c in celltypes) {
  print(c)
  up.gene.list <- makeGeneList2(cosmic.hlca.covid.list.up[[c]])
  down.gene.list <- makeGeneList2(cosmic.hlca.covid.list.down[[c]])
  comparisons[["COSMIC"]][[c]]<- up.gene.list
  comparisons[["HLCA"]][[c]] <- down.gene.list
}

library(fgsea)

groups <- c("COSMIC", "HLCA")
res = list()
for (c in celltypes){
  for (g in groups){
    res[[c]][[g]] <- fgsea(pathways = pathwaysH, stats = comparisons[[g]][[c]], nperm = 10000, minSize = 0, maxSize =10000)
  }
}

for (c in celltypes){
  for (g in groups){
    res[[c]][[g]] <- as.data.frame(res[[c]][[g]])
  }
}

result = res

for(i in 1:length(comparisons)){
  result[[i]] <- lapply(result[[i]], function(x){
    x$ranking <- -log10(x$pval)*sign(x$NES)
    x <- x[order(x$pathway), ]
    return(x)
  })
}


result <- lapply(result ,function(x){
  x[['COSMIC']]$group = "COSMIC"
  x[['HLCA']]$group = "HLCA"
  return(x)
})

result2 <- lapply(result, function(x) {
  y <- do.call(rbind, x)
  return(y)
})


result3 <- lapply(result2, function(x) {
  x$group <- rownames(x)
  x$group <- gsub("cal.*","cal", x$group)
  return(x)
})

data<-rbindlist(lapply(1:length(result3)
                       , function(x){ setDT(result3[[x]])[
                         , id:=names(result3)[x]]})
                , use.names=TRUE, fill=TRUE)

datax <- as.data.frame(data)
datax <- apply(datax,2,as.character)

data3 <- data %>% select("NES", "pathway", "id", "group")
data3$group <- gsub("\\..*", " ", data3$group)
data3$pathway <- gsub("HALLMARK_|", "", data3$pathway)
data3$pathway <- gsub("_", " ", data3$pathway)
data4 <- data3 %>% pivot_wider(names_from = pathway, values_from = NES, values_fn = sum, values_fill = 0)
data4 <- as.data.frame(data4)
data4$sample_id <- paste(data4$id, data4$group, sep = "_")
row.names(data4) <- data4$sample_id
data4$sample_id <- NULL
data5 <- subset(data4, subset = group == "COSMIC ")
row.names(data5) <- data5$id
data5 <- data5 %>% select(-c("id", "group"))
data6 <- data4 %>% select("id", "group") 

color_palette <- function(colors) {
  # has_package("grDevices")
  grDevices::colorRampPalette(colors)(n = 1000)
}

sig_palette <- color_palette(c("blue", "white", "red"))


# Create dendrograms for rows and columns
row_dend <- as.dendrogram(hclust(dist(as.matrix(t(data5)))))
col_dend <- as.dendrogram(hclust(dist(t(as.matrix(t(data5))))))

# Create the heatmap with dendrograms
pheatmap(mat, add.dendrogram = "row", row.dendrogram = row_dend)
pheatmap(mat, add.dendrogram = "column", col.dendrogram = col_dend)

# Create the heatmap with dendrograms
pheatmap(as.matrix(t(data5), add.dendrogram = "row", row.dendrogram = row_dend))

p3 <- pheatmap(as.matrix(t(data5)), 
               treeheight_row = 25, treeheight_col = 10, fontsize=10, border_color="white", cellwidth=12)  ; p3

data$leadingEdge <- gsub('c\\(\\"',"",data$leadingEdge, perl = TRUE)
data$leadingEdge <- gsub('\\"\\)',"",data$leadingEdge, perl = TRUE)
data$leadingEdge <- gsub('"',"",data$leadingEdge, perl = TRUE)

data1 <- concat.split.multiple(
  data, split.cols = c("leadingEdge"),
  seps = ",", direction = "long")

data$pathway <- gsub("HALLMARK_|", "", data$pathway)
data4 <- data %>% pivot_wider(names_from = pathway, values_from = NES, values_fn = sum, values_fill = 0)
data4 <- as.data.frame(data4)
data4$group <- gsub("\\..*", " ", data4$group)
data5 <- data4 %>% select(c("id", "group", "padj", "leadingEdge",  "INFLAMMATORY_RESPONSE", "INTERFERON_ALPHA_RESPONSE", "INTERFERON_GAMMA_RESPONSE", "TNFA_SIGNALING_VIA_NFKB","COMPLEMENT","IL2_STAT5_SIGNALING","IL6_JAK_STAT3_SIGNALING", ))
data4 <- data5 %>% pivot_longer(cols = (c( "INTERFERON_GAMMA_RESPONSE") ))
data7 <-  data4 %>% filter(value != 0)  %>% filter(group == "COSMIC ")

data8 <- concat.split.multiple(
  data7, split.cols = c("leadingEdge"),
  seps = ",", direction = "long")

pathways <- data8 %>% select(leadingEdge) %>% unique()

set.seed(100)
p <- DotPlot(merged_lung, group.by = "source", idents = c("Alveolar macrophages"), assay = "RNA", cols = c("RdBu"), features = pathways$leadingEdge) + ggpubr::rotate_x_text() +
  theme(legend.position = 'right', plot.title = element_text(hjust = 0.5), text = element_text(size = 12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(face = "italic", angle = 45, hjust = 1, vjust = 1), 
        axis.text.y = element_text(face = "italic", hjust = 1, vjust = 1), 
        plot.margin = margin(t = 20, r = 0, b = 0, l = 0, unit = 'pt')) + coord_flip(); p