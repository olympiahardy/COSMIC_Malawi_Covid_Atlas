# Lung atlas broad
lung.atlas.cleaned$Group[lung.atlas.cleaned$Group == "Pneumonia"] <- "LRTD"
lung.atlas.cleaned$Group[lung.atlas.cleaned$Group == "Non-Pneumonia"] <- "No_LRTD"

SaveH5Seurat(lung.atlas.cleaned, filename = "COSMIC_Lung_Atlas.h5Seurat")
Convert("COSMIC_Lung_Atlas.h5Seurat", dest = "h5ad")

# Lung atlas immune

lung.immune.cleaned <- readRDS("/datastore3/Olympia/COSMIC_data/Final_lung_objects/lung_immune_object120923.rds")
Idents(lung.immune.cleaned) <- "celltype.final"

Idents(object = lung.immune.cleaned, cells = gdtcell.barcodes) <- "gdT cells"
lung.immune.cleaned <- RenameIdents(object = lung.immune.cleaned, `CD16+ Neutrophils` = "Mature neutrophils")
lung.immune.cleaned <- RenameIdents(object = lung.immune.cleaned, `CD16- Neutrophils` = "Immature neutrophils")
lung.immune.cleaned$celltype.annotation <- Idents(lung.immune.cleaned)
lung.immune.cleaned$celltype.annotation <- as.character(lung.immune.cleaned$celltype.annotation)
lung.immune.cleaned$Group[lung.immune.cleaned$Group == "Pneumonia"] <- "LRTD"
lung.immune.cleaned$Group[lung.immune.cleaned$Group == "Non-Pneumonia"] <- "No_LRTD"
m <- lung.immune.cleaned@meta.data
lung.immune.cleaned@meta.data[,26:49] <- NULL
DefaultAssay(lung.immune.cleaned) <- "RNA"
SaveH5Seurat(lung.immune.cleaned, filename = "COSMIC_Lung_Immune_Atlas.h5Seurat")
Convert("COSMIC_Lung_Immune_Atlas.h5Seurat", dest = "h5ad")

# Lung atlas stromal
lung.stromal.cleaned <- readRDS("/datastore3/Olympia/COSMIC_data/Final_lung_objects/lung_stromal_object120923.rds")
lung.stromal.cleaned <- subset(lung.stromal.cleaned, idents = "Alveolar macrophages", invert = T)
m <- lung.stromal.cleaned@meta.data
lung.stromal.cleaned@meta.data[,26:47] <- NULL
DefaultAssay(lung.stromal.cleaned) <- "RNA"
names(lung.stromal.cleaned@meta.data)[27] ="celltype.annotation"
lung.stromal.cleaned$celltype.annotation <- as.character(lung.stromal.cleaned$celltype.annotation)
lung.stromal.cleaned$Group[lung.stromal.cleaned$Group == "Pneumonia"] <- "LRTD"
lung.stromal.cleaned$Group[lung.stromal.cleaned$Group == "Non-Pneumonia"] <- "No_LRTD"
SaveH5Seurat(lung.stromal.cleaned, filename = "COSMIC_Lung_Stromal_Atlas.h5Seurat")
Convert("COSMIC_Lung_Stromal_Atlas.h5Seurat", dest = "h5ad")

# Nasal atlas
nasal.atlas.cleaned <- nasal.atlas
m <- nasal.atlas.cleaned@meta.data
nasal.atlas.cleaned@meta.data[,12:19] <- NULL
nasal.atlas.cleaned@meta.data[,16:17] <- NULL
nasal.atlas.cleaned@meta.data[,18:38] <- NULL
nasal.atlas.cleaned@meta.data[,19:21] <- NULL
nasal.atlas.cleaned@meta.data[,17] <- NULL
names(nasal.atlas.cleaned@meta.data)[17] ="celltype.annotation"
DefaultAssay(nasal.atlas.cleaned) <- "RNA"
nasal.atlas.cleaned$celltype.annotation <- as.character(nasal.atlas.cleaned$celltype.annotation)
nasal.atlas.cleaned$Group[nasal.atlas.cleaned$Group == "Pneumonia"] <- "LRTD"
nasal.atlas.cleaned$Group[nasal.atlas.cleaned$Group == "Non-Pneumonia"] <- "No_LRTD"

SaveH5Seurat(nasal.atlas.cleaned, filename = "COSMIC_Nasal_Atlas.h5Seurat")
Convert("COSMIC_Nasal_Atlas.h5Seurat", dest = "h5ad")

# Blood atlas
blood.atlas.cleaned <- blood.atlas
m <- blood.atlas.cleaned@meta.data
blood.atlas.cleaned@meta.data[,12:19] <- NULL
blood.atlas.cleaned@meta.data[,16] <- NULL
blood.atlas.cleaned@meta.data[,17] <- NULL
names(blood.atlas.cleaned@meta.data)[16] ="celltype.annotation"
DefaultAssay(blood.atlas.cleaned) <- "RNA"
blood.atlas.cleaned$celltype.annotation <- as.character(blood.atlas.cleaned$celltype.annotation)
blood.atlas.cleaned$Group[blood.atlas.cleaned$Group == "Pneumonia"] <- "LRTD"
blood.atlas.cleaned$Group[blood.atlas.cleaned$Group == "Non-Pneumonia"] <- "No_LRTD"

SaveH5Seurat(blood.atlas.cleaned, filename = "COSMIC_Blood_Atlas.h5Seurat")
Convert("COSMIC_Blood_Atlas.h5Seurat", dest = "h5ad")








