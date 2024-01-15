library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(stringr)

# visualization config
options(bitmapType="cairo")

PROJECT.DIR <- "/net/pupil1/home/people/s223162/sc_bioinfo/project/scRNAseq_project/"
DATA.DIR <- paste(PROJECT.DIR, "data/", sep="")

# load the Seurat objects
retina.scrna <- readRDS(paste(DATA.DIR, "retina_preprocessed.rds", sep=""))
retina.multimodal <- readRDS(paste(DATA.DIR, "retina_2.rds", sep=""))


# <<< 1. scRNA Clustering >>>

# normalization, scaling, and feature finding
retina.scrna <- SCTransform(
  retina.scrna, vars.to.regress=c("percent.mt"), ncells=3000,
  conserve.memory=TRUE, verbose=TRUE)

# PCA
retina.scrna <- RunPCA(retina.scrna)
ElbowPlot(retina.scrna)

# clustering
retina.scrna <- FindNeighbors(retina.scrna, dims=1:20)
retina.scrna <- FindClusters(retina.scrna, resolution=seq(0.1, 0.8, 0.1))
retina.scrna <- RunUMAP(retina.scrna, dims=1:20)

# save clustered object
saveRDS(retina.scrna, paste(DATA.DIR, "retina.sct.clust.rds", sep=""))
#retina.scrna <- readRDS(paste(DATA.DIR, "retina.sct.clust.rds", sep=""))