#' seurat_pipeline method
#'
#' @param input_data A gene-cell matrix or unprocessed Seurat object of single cell expression data.
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
#' @param k.param Defines k for the k-nearest neighbor algorithm. Default: 50.
#' @param dims Dimensions of reduction to use as input. Default: 1:10.
#' @param resolution Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities. Default: 0.8.
#' @param n.neighbors This determines the number of neighboring points used in local approximations of manifold structure. Larger values will result in more global structure being preserved at the loss of detailed local structure. In general this parameter should often be in the range 5 to 50.
#' @param min.dist This controls how tightly the embedding is allowed compress points together. Larger values ensure embedded points are moreevenly distributed, while smaller values allow the algorithm to optimise more accurately with regard to local structure. Sensible values are in the range 0.001 to 0.5.
#' @param qc Whether filter cells by gene number and mito ratio.
#' @param mad The fold multiple standard deviation of gene number to calculate max gene number of cell, Default: 2.
#' @param min.gene The minimal gene number of cell, Default:500.
#' @param mt.rate The threshold of mt gene expression, Default:10.
#'
#' @return The processed Seurat object.
#'
#' @importFrom Seurat CreateSeuratObject NormalizeData FindVariableFeatures ScaleData GetAssayData RunPCA FindNeighbors FindClusters RunUMAP
#'
#' @export
#'
#'
seurat_pipeline = function(input_data,
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
                           k.param = 50,
                           dims = 1:35,
                           resolution = 0.8,
                           n.neighbors = 30,
                           min.dist = 0.3,
                           qc = TRUE,
                           mad = 2,
                           min.gene = 500,
                           mt.rate = 10
                           ){

  if(class(input_data) == 'Seurat'){
    object = input_data
  }else{
    object = CreateSeuratObject(input_data,
                                min.cells = min.cells,
                                min.features = min.features)
  }

  if(qc){

    nFeature_highthresholds = mean(object$nFeature_RNA) + mad * sd(object$nFeature_RNA)
    object[["percent.mt"]] = PercentageFeatureSet(object, pattern = "^MT-|^mt-")
    object = subset(object, subset = nFeature_RNA > min.gene & nFeature_RNA < nFeature_highthresholds & percent.mt < mt.rate)

  }


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
                   n.neighbors = n.neighbors,
                   min.dist = min.dist,
                   dims = dims,
                   verbose = F)

  return(object)

}

#' Batch correct method
#'
#' @param object_list A list contains seurat object.
#' @param method Select a method to correct the batch effect among list, eg. seurat, harmony, and liger. Default: seurat.
#' @param normalization.method Method for normalization. Default: LogNormalize.
#' @param scale.factor Sets the scale factor for cell-level normalization. Default: 10000.
#' @param selection.method How to choose top variable features. Default: vst.
#' @param nfeatures Number of features to select as top variable features; only used when selection.method is set to 'dispersion' or 'vst'. Default: 2000.
#' @param mean.cutoff A two-length numeric vector with low- and high-cutoffs for feature means. Default: c(0.1, 8).
#' @param features Vector of features names to scale/center. Default is variable features. Default: NULL.
#' @param dispersion.cutoff A two-length numeric vector with low- and high-cutoffs for feature dispersions. Default: c(1, Inf).
#' @param nfeatures Number of features to select as top variable features among list, Default:2000.
#' @param reference A vector specifying the object/s to be used as a reference during integration. If NULL (default), all pairwise anchors are found (no reference/s). If not NULL, the corresponding objects in object.list will be used as references. When using a set of specified references, anchors are first found between each query and each reference. The references are then integrated through pairwise integration. Each query is then mapped to the integrated reference.
#' @param dims Which dimensions to use from the CCA to specify the neighbor search space, Default:30.
#' @param k.anchor How many neighbors (k) to use when picking anchors, Default:5.
#' @param k.filter How many neighbors (k) to use when filtering anchors, Default:200.
#' @param k.score How many neighbors (k) to use when scoring anchors, Default:30.
#' @param k.weight Number of neighbors to consider when weighting anchors, Default:100.
#' @param features_scale Vector of features names to scale/center. Default is variable features. Default: NULL.
#' @param vars.to.regress Variables to regress out (previously latent.vars in RegressOut). For example, nUMI, or percent.mito. Default: NULL.
#' @param npcs Total Number of PCs to compute and store (50 by default).
#' @param k.param Defines k for the k-nearest neighbor algorithm. Default: 50.
#' @param dims Dimensions of reduction to use as input. Default: 1:10.
#' @param resolution Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities. Default: 0.8.
#' @param n.neighbors This determines the number of neighboring points used in local approximations of manifold structure. Larger values will result in more global structure being preserved at the loss of detailed local structure. In general this parameter should often be in the range 5 to 50.
#' @param min.dist This controls how tightly the embedding is allowed compress points together. Larger values ensure embedded points are moreevenly distributed, while smaller values allow the algorithm to optimise more accurately with regard to local structure. Sensible values are in the range 0.001 to 0.5.
#'
#' @return The batch corrected Seurat object.
#'
#' @import Seurat harmony SeuratWrappers
#'
#' @export
#'
#'
batch_correct_pipeline = function(object_list,
                                  method,
                                  normalization.method = 'LogNormalize',
                                  scale.factor = 10000,
                                  selection.method = 'vst',
                                  nfeatures = 2000,
                                  mean.cutoff = c(0.1, 8),
                                  dispersion.cutoff = c(1, Inf),
                                  reference = NULL,
                                  dims = 1:30,
                                  k.anchor = 5,
                                  k.filter = 200,
                                  k.score = 30,
                                  k.weight = 100,
                                  features = NULL,
                                  features_scale = NULL,
                                  vars.to.regress = NULL,
                                  npcs = 50,
                                  k.param = 50,
                                  resolution = 0.8,
                                  n.neighbors = 30,
                                  min.dist = 0.3){

  if(method == 'seurat'){

    object_list = lapply(object_list, function(x) {

      x = NormalizeData(x,
                        normalization.method = normalization.method,
                        scale.factor = scale.factor,
                        verbose = F)

      x = FindVariableFeatures(x,
                               selection.method = selection.method,
                               nfeatures = nfeatures,
                               mean.cutoff = mean.cutoff,
                               dispersion.cutoff = dispersion.cutoff,
                               verbose = F)

    })

    features = SelectIntegrationFeatures(object_list, nfeatures = nfeatures)

    anchors = FindIntegrationAnchors(object.list = object_list,
                                     anchor.features = features,
                                     reference = reference,
                                     dims = dims,
                                     k.anchor = k.anchor,
                                     k.filter = k.filter,
                                     k.score = k.score
    )

    combined = IntegrateData(anchorset = anchors,
                             dims = dims,
                             k.weight = k.weight)

    combined = ScaleData(combined,
                         features = features_scale,
                         vars.to.regress = vars.to.regress,
                         verbose = F)

    combined = RunPCA(combined,
                      npcs = npcs,
                      verbose = F)

    combined = FindNeighbors(combined,
                             k.param = k.param,
                             dims = dims,
                             verbose = F)

    combined = FindClusters(combined,
                            resolution = resolution,
                            verbose = F)

    combined = RunUMAP(combined,
                       n.neighbors = n.neighbors,
                       min.dist = min.dist,
                       dims = dims,
                       verbose = F)

    return(combined)

  }

  if(method == 'harmony'){

    for(i in 1:length(object_list)){

      object_list[[i]]$harmony_batch = i

    }

    object_merge = merge(object_list[[1]], object_list[-1])

    object_merge = NormalizeData(object_merge,
                                 normalization.method = normalization.method,
                                 scale.factor = scale.factor,
                                 verbose = F)

    object_merge = FindVariableFeatures(object_merge,
                                        selection.method = selection.method,
                                        nfeatures = nfeatures,
                                        mean.cutoff = mean.cutoff,
                                        dispersion.cutoff = dispersion.cutoff,
                                        verbose = F)

    object_merge = ScaleData(object_merge,
                             features = features,
                             vars.to.regress = vars.to.regress,
                             verbose = F)

    object_merge = RunPCA(object_merge,
                          npcs = npcs,
                          verbose = F)

    object_merge$harmony_batch = as.factor(object_merge$harmony_batch)

    object_merge = RunHarmony(object_merge,
                              group.by.vars = 'harmony_batch')

    object_merge = FindNeighbors(object_merge,
                                 reduction = "harmony",
                                 k.param = k.param,
                                 dims = dims,
                                 verbose = F)

    object_merge = FindClusters(object_merge,
                                resolution = resolution,
                                verbose = F)

    object_merge = RunUMAP(object_merge,
                           reduction = "harmony",
                           n.neighbors = n.neighbors,
                           min.dist = min.dist,
                           dims = dims,
                           verbose = F)

    return(object_merge)

  }

  if(method == 'liger'){

    for(i in 1:length(object_list)){

      object_list[[i]]$liger_batch = i

    }

    object_merge = merge(object_list[[1]], object_list[-1])

    object_merge = NormalizeData(object_merge,
                                 normalization.method = normalization.method,
                                 scale.factor = scale.factor,
                                 verbose = F)

    object_merge = FindVariableFeatures(object_merge,
                                        selection.method = selection.method,
                                        nfeatures = nfeatures,
                                        mean.cutoff = mean.cutoff,
                                        dispersion.cutoff = dispersion.cutoff,
                                        verbose = F)

    object_merge = ScaleData(object_merge,
                             features = features,
                             split.by = 'liger_batch',
                             vars.to.regress = vars.to.regress,
                             do.center = FALSE,
                             verbose = F)

    object_merge = RunOptimizeALS(object_merge, k = 20, lambda = 5, split.by = "liger_batch")

    object_merge = RunQuantileNorm(object_merge, split.by = "liger_batch")

    object_merge = FindNeighbors(object_merge,
                                 reduction = "iNMF",
                                 k.param = k.param,
                                 dims = 1:20,
                                 verbose = F)

    object_merge = FindClusters(object_merge,
                                resolution = resolution,
                                verbose = F)

    object_merge = RunUMAP(object_merge,
                           reduction = "iNMF",
                           n.neighbors = n.neighbors,
                           min.dist = min.dist,
                           dims = 1:20,
                           verbose = F)

    return(object_merge)

  }

}

#' Automatic celltype annotation through reference mapping
#'
#' @param object A seurat object needs to be annotated.
#' @param reference A reference seurat object containing celltype info.
#'
#' @return The annotated Seurat object.
#'
#' @import SeuratDisk Seurat
#'
#' @export
#'
#'
reference_mapping = function(object,
                             reference){

  object = SCTransform(object, verbose = FALSE)

  anchors = FindTransferAnchors(
    reference = reference,
    query = object,
    normalization.method = "SCT",
    reference.reduction = "spca",
    dims = 1:50
  )

  object = MapQuery(
    anchorset = anchors,
    query = object,
    reference = reference,
    refdata = list(
      celltype.l1 = "celltype.l1",
      celltype.l2 = "celltype.l2",
      predicted_ADT = "ADT"
    ),
    reference.reduction = "spca",
    reduction.model = "wnn.umap"
  )

  return(object)

}



