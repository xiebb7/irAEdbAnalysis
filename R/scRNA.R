#' seurat_pipeline method
#'
#' @param mat_data A gene-cell matrix of single cell expression data.
#' @param min.cells Include features detected in at least this many cells. Will subset the counts matrix as well. To reintroduce excluded features, create a new object with a lower cutoff. Default: 0.
#' @param min.features Include cells where at least this many features are detected. Default: 0.
#' @param normalization.method Method for normalization. Default: LogNormalize.
#' @param scale.factor Sets the scale factor for cell-level normalization. Default: 10000.
#' @param selection.method How to choose top variable features. Default: vst.
#' @param nfeatures Number of features to select as top variable features; only used when selection.method is set to 'dispersion' or 'vst'. Default: 2000.
#' @param mean.cutoff A two-length numeric vector with low- and high-cutoffs for feature means. Default: c(0.1, 8).
#' @param dispersion.cutoff A two-length numeric vector with low- and high-cutoffs for feature dispersions. Default: c(1, Inf).
#' @param features Vector of features names to scale/center. Default is variable features. Default: NULL.
#' @param vars.to.regress Variables to regress out (previously latent.vars in RegressOut). For example, nUMI, or percent.mito. Default: NULL.
#' @param npcs Total Number of PCs to compute and store (50 by default).
#' @param k.param Defines k for the k-nearest neighbor algorithm. Default: 20.
#' @param dims Dimensions of reduction to use as input. Default: 1:10.
#' @param resolution Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities. Default: 0.8.
#' @param n.neighbors This determines the number of neighboring points used in local approximations of manifold structure. Larger values will result in more global structure being preserved at the loss of detailed local structure. In general this parameter should often be in the range 5 to 50.
#' @param min.dist This controls how tightly the embedding is allowed compress points together. Larger values ensure embedded points are moreevenly distributed, while smaller values allow the algorithm to optimise more accurately with regard to local structure. Sensible values are in the range 0.001 to 0.5.
#'
#' @return The predict cell type.
#'
#' @importFrom Seurat CreateSeuratObject NormalizeData FindVariableFeatures ScaleData GetAssayData RunPCA FindNeighbors FindClusters RunUMAP
#'
#' @export
#'
#'
seurat_pipeline = function(mat_data,
                           min.cells = 0,
                           min.features = 0,
                           normalization.method = 'LogNormalize',
                           scale.factor = 10000,
                           selection.method = 'vst',
                           nfeatures = 2000,
                           mean.cutoff = c(0.1, 8),
                           dispersion.cutoff = c(1, Inf),
                           features = NULL,
                           vars.to.regress = NULL,
                           npcs = 50,
                           k.param = 20,
                           dims = 1:35,
                           resolution = 0.8,
                           n.neighbors = 30,
                           min.dist = 0.3
                           ){

  object = CreateSeuratObject(mat_data,
                              min.cells = min.cells,
                              min.features = min.features)

  object = NormalizeData(object,
                         normalization.method = normalization.method,
                         scale.factor = scale.factor,
                         verbose = F)

  object = FindVariableFeatures(object,
                                selection.method = selection.method,
                                nfeatures = nfeatures,
                                mean.cutoff = mean.cutoff,
                                dispersion.cutoff = dispersion.cutoff,
                                verbose = F)

  object = ScaleData(object,
                     features = features,
                     vars.to.regress = vars.to.regress,
                     verbose = F)

  object = RunPCA(object,
                  npcs = npcs,
                  verbose = F)

  object = FindNeighbors(object,
                         k.param = k.param,
                         dims = dims,
                         verbose = F)

  object = FindClusters(object,
                        resolution = resolution,
                        verbose = F)

  object = RunUMAP(object,
                   resolution = resolution,
                   n.neighbors = n.neighbors,
                   min.dist = min.dist,
                   dims = dims,
                   verbose = F)

  return(object)

}
