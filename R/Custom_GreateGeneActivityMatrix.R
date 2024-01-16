library(Seurat)

library(SeuratData)
library(patchwork)
library(dplyr)
library(tidyr)
library(multtest)
library(metap)
library(tibble)
library(purrr)
library(future)


CreateGeneActivityMatrix <- function(
    peak.matrix,
    annotation.file,
    seq.levels = c(1:22, "X", "Y"),
    include.body = TRUE,
    upstream = 2000,
    downstream = 0,
    verbose = TRUE
) {
  
  # convert peak matrix to GRanges object
  peak.df <- rownames(x = peak.matrix)
  peak.df <- do.call(what = rbind, args = strsplit(x = gsub(peak.df, pattern = ":", replacement = "-"), split = "-"))
  peak.df <- as.data.frame(x = peak.df)
  colnames(x = peak.df) <- c("chromosome", 'start', 'end')
  peaks.gr <- GenomicRanges::makeGRangesFromDataFrame(df = peak.df)
  
  # if any peaks start at 0, change to 1
  # otherwise GenomicRanges::distanceToNearest will not work
  BiocGenerics::start(peaks.gr[BiocGenerics::start(peaks.gr) == 0, ]) <- 1
  # get annotation file, select genes
  #   anno <- rtracklayer::import(con = annotation.file)
  # anno <- GenomeInfoDb::keepSeqlevels(x = anno, value = seq.levels, pruning.mode = 'coarse')
  
  gtf <- rtracklayer::import(con = annotation.file)
  gtf <- GenomeInfoDb::keepSeqlevels(x = gtf, value = seq.levels, pruning.mode = 'coarse')
  
  # change seqlevelsStyle if not the same
  if (!any(GenomeInfoDb::seqlevelsStyle(x = gtf) == GenomeInfoDb::seqlevelsStyle(x = peaks.gr))) {
    GenomeInfoDb::seqlevelsStyle(gtf) <- GenomeInfoDb::seqlevelsStyle(peaks.gr)
  }
  gtf.genes <- gtf[gtf$type == 'gene']
  
  # 
  # Extend definition up/downstream
  if (include.body) {
    gtf.body_prom <- Signac::Extend(x = gtf.genes, upstream = upstream, downstream = downstream)
  } else {
    gtf.body_prom <- SummarizedExperiment::promoters(x = gtf.genes, upstream = upstream, downstream = downstream)
  }
  gene.distances <- GenomicRanges::distanceToNearest(x = peaks.gr, subject = gtf.body_prom)
  keep.overlaps <- gene.distances[rtracklayer::mcols(x = gene.distances)$distance == 0]
  peak.ids <- peaks.gr[S4Vectors::queryHits(x = keep.overlaps)]
  gene.ids <- gtf.genes[S4Vectors::subjectHits(x = keep.overlaps)]
  
  # Some gtf rows will not have gene_name attribute
  # Replace it by gene_id attribute
  # 
  gene.ids$Name[is.na(gene.ids$Name)] <- gene.ids$gene_id[is.na(gene.ids$Name)]
  
  peak.ids$gene.name <- gene.ids$Name
  peak.ids <- as.data.frame(x = peak.ids)
  peak.ids$peak <- rownames(peak.matrix)[S4Vectors::queryHits(x = keep.overlaps)]
  annotations <- peak.ids[, c('peak', 'gene.name')]
  colnames(x = annotations) <- c('feature', 'new_feature')
  
  # collapse into expression matrix
  peak.matrix <- as(object = peak.matrix, Class = 'matrix')
  all.features <- unique(x = annotations$new_feature)
  
  if (future::nbrOfWorkers() > 1) {
    mysapply <- future.apply::future_sapply
  } else {
    mysapply <- ifelse(test = verbose, yes = pbapply::pbsapply, no = sapply)
  }
  newmat <- mysapply(X = 1:length(x = all.features), FUN = function(x){
    
    features.use <- annotations[annotations$new_feature == all.features[[x]], ]$feature
    
    submat <- peak.matrix[features.use, ]
    
    if (length(x = features.use) > 1) {
      return(Matrix::colSums(x = submat))
    } else {
      return(submat)
    }
  })
  
  newmat <- t(x = newmat)
  
  rownames(x = newmat) <- all.features
  
  colnames(x = newmat) <- colnames(x = peak.matrix)
  
  return(as(object = newmat, Class = 'dgCMatrix'))
}
