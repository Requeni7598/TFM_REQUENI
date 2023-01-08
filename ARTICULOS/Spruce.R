# install packages
# install.packages("Seurat","spruce")
devtools::install_github('satijalab/seurat-data')

# load packages
library(Seurat) # Esta es una de las tres técnicas previas propuestas.
                # Faltaría análisis bajo Giotto y stLearn.
library(SeuratData)
library(spruce)
# Mirar donde se instalan todos los repositorios de paquetes en nuestro ordenador.

# load data
# InstallData("stxBrain")
brain <- LoadData("stxBrain", type = "anterior1")
# Esta función es del paquete SeuratData

# normalize using sctransform
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
# Proponer otras técnicas de normalización en nuestro estudio.  

DefaultAssay(brain) <- "SCT"
brain <- RunPCA(brain)

# identify the most spatially variable features (SVGs) in the data
# using the markvariogram approach
brain <- FindSpatiallyVariableFeatures(brain, 
                                       assay = "SCT", 
                                       features = VariableFeatures(brain)[1:1000],
                                       selection.method = "markvariogram")
# Técnicas de variograma. Podemos tratar de aplicar conocimientos propios en la materia.
# Lleva más de 4 horas compilando

# identify the most spatially variable features (SVGs) in the data
# using Morans I statistic
brain <- FindSpatiallyVariableFeatures(brain, 
                                       assay = "SCT", 
                                       features = VariableFeatures(brain)[1:1000],
                                       selection.method = "moransi")

# Identify clusters using Seurat's graph clustering method
brain <- FindNeighbors(brain)
brain <- FindClusters(brain)

# plot results
SpatialDimPlot(brain, label = TRUE, label.size = 3)

# fit spruce
# HVGs
fit_K6_HVGs <- fit_spruce(brain, K = 6, emb = "HVGs")

brain$z_HVGs <- as.factor(fit_K6_HVGs$z)
Idents(brain) <- "z_HVGs"
SpatialDimPlot(brain, label = TRUE, label.size = 3)

#SVGs
fit_K6_SVGs <- fit_spruce(brain, K = 6, emb = "SVGs")

brain$z_SVGs <- as.factor(fit_K6_SVGs$z)
Idents(brain) <- "z_SVGs"
SpatialDimPlot(brain, label = TRUE, label.size = 3)

# PCs
fit_K6_PCs <- fit_spruce(brain, K = 6, emb = "PCs")

brain$z_PCs <- as.factor(fit_K6_PCs$z)
Idents(brain) <- "z_PCs"
SpatialDimPlot(brain, label = TRUE, label.size = 3)

