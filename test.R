#########################################################################################
# R script to run FlowSOM_pre and FlowSOM
#
# Lukas Weber, August 2016
#########################################################################################


library(flowCore)
library(FlowSOM)
library(clue)



#################
### LOAD DATA ###
#################

# filenames

DATA_DIR <- "./FlowRepository_FR-FCM-ZZPH_files"

files <- list(
  Levine_32dim = file.path(DATA_DIR, "Levine_32dim.fcs"), 
  Levine_13dim = file.path(DATA_DIR, "Levine_13dim.fcs"), 
  Samusik_01   = file.path(DATA_DIR, "Samusik_01.fcs"), 
  Samusik_all  = file.path(DATA_DIR, "Samusik_all.fcs") 
)

# FlowCAP data sets are treated separately since they require clustering algorithms to be
# run individually for each sample

is_FlowCAP <- c(FALSE, FALSE, FALSE, FALSE)

# indices of protein marker columns

marker_cols <- list(
  Levine_32dim = 5:36, 
  Levine_13dim = 1:13, 
  Samusik_01   = 9:47, 
  Samusik_all  = 9:47
)
sapply(marker_cols, length)



# load data files: FlowSOM requires flowFrame objects

data <- vector("list", length(files))
names(data) <- names(files)

for (i in 1:length(data)) {
  f <- files[[i]]
  
  if (!is_FlowCAP[i]) {
    data[[i]] <- flowCore::exprs(flowCore::read.FCS(f, transformation = FALSE, truncate_max_range = FALSE))
    data[[i]] <- data[[i]][,marker_cols[[i]]]
    data[[i]] <- ifelse(data[[i]]<0, 0, data[[i]])
    
  } else {
    smp <- flowCore::exprs(flowCore::read.FCS(f, transformation = FALSE, truncate_max_range = FALSE))
    smp <- smp[, "sample"]
    d <- flowCore::read.FCS(f, transformation = FALSE, truncate_max_range = FALSE)
    data[[i]] <- flowCore::split(d, smp)
  }
}


#############################################################
### Run FlowSOM_pre: manually selected number of clusters ###
#############################################################

# run FlowSOM_pre with manually selected number of clusters

# grid size 20x20 (400 clusters) for Mosmann_rare (data set with very rare population);
# and grid size 10x10 (100 clusters, i.e. default) for all other data sets

# grid sizes
grid_size <- list(
  Levine_32dim = 10, 
  Levine_13dim = 10, 
  Samusik_01   = 10, 
  Samusik_all  = 10
)

seed <- 1000
out <- runtimes <- vector("list", length(data))
names(out) <- names(runtimes) <- names(data)

for (i in 1:length(data)) {
  
  if (!is_FlowCAP[i]) {
    set.seed(seed)
    runtimes[[i]] <- system.time({
      # fSOM <- FlowSOM::ReadInput(data[[i]], transform = TRUE, toTransform = 1:dim(data[[i]])[2],
      #                            transformFunction = flowCore::logicleTransform(w = 0.1, m=5, t = 5000, a=0),
      #                            scale = FALSE)
      fSOM <- FlowSOM::ReadInput(data[[i]], transform = FALSE, scale = FALSE)
      # fSOM$data <- asinh(fSOM$data/0.5)
      fSOM <- FlowSOM::BuildSOM(fSOM, colsToUse = 1:dim(data[[i]])[2],
                                xdim = grid_size[[i]], ydim = grid_size[[i]])
      fSOM <- FlowSOM::BuildMST(fSOM)
    })
    out[[i]] <- fSOM
    cat("data set", names(data[i]), ": run complete\n")
    
  } else {
    # FlowCAP data sets: run clustering algorithm separately for each sample
    out[[i]] <- runtimes[[i]] <- vector("list", length(data[[i]]))
    names(out[[i]]) <- names(runtimes[[i]]) <- names(data[[i]])
    
    for (j in 1:length(data[[i]])) {
      set.seed(seed)
      runtimes[[i]][[j]] <- system.time({
        fSOM <- FlowSOM::ReadInput(data[[i]][[j]], transform = FALSE, scale = FALSE)
        fSOM <- FlowSOM::BuildSOM(fSOM, colsToUse = marker_cols[[i]], 
                                  xdim = grid_size[[i]], ydim = grid_size[[i]])
        fSOM <- FlowSOM::BuildMST(fSOM)
      })
      out[[i]][[j]] <- fSOM
    }
    cat("data set", names(data[i]), ": run complete\n")
    
    # FlowCAP data sets: sum runtimes over samples
    runtimes_i <- do.call(rbind, runtimes[[i]])[, 1:3]
    runtimes_i <- colSums(runtimes_i)
    names(runtimes_i) <- c("user", "system", "elapsed")
    runtimes[[i]] <- runtimes_i
  }
}

# store output and runtimes for meta-clustering step below
out_pre_manual <- out
runtimes_pre_manual <- runtimes

# example of FlowSOM plots (one data set only)
FlowSOM::PlotStars(out[[1]])

# example showing how to extract cluster labels (one data set only)
str(out[[1]]$map)
head(out[[1]]$map$mapping)
dim(out[[1]]$map$mapping)

# extract cluster labels
clus <- vector("list", length(data))
names(clus) <- names(data)

for (i in 1:length(clus)) {
  if (!is_FlowCAP[i]) {
    clus[[i]] <- out[[i]]$map$mapping[, 1]
    
  } else {
    # FlowCAP data sets
    clus_list_i <- lapply(out[[i]], function(o) o$map$mapping[, 1])
    names(clus_list_i) <- names(out[[i]])
    
    # convert FlowCAP cluster labels into format "sample_number"_"cluster_number"
    # e.g. sample 1, cluster 3 -> cluster label 1.3
    names_i <- rep(names(clus_list_i), times = sapply(clus_list_i, length))
    clus_collapse_i <- unlist(clus_list_i, use.names = FALSE)
    clus[[i]] <- paste(names_i, clus_collapse_i, sep = "_")
  }
}

sapply(clus, length)


###########################################################################################
### Run FlowSOM (additional meta-clustering step): manually selected number of clusters ###
###########################################################################################

# run FlowSOM (additional meta-clustering step) with manually selected number of clusters

# using results from above (stored in object "out_pre_manual")

# number of clusters k
k0 <- 40
k <- list(
  Levine_32dim = k0, 
  Levine_13dim = k0, 
  Samusik_01   = k0, 
  Samusik_all  = k0
)

seed <- 1000
out <- runtimes <- vector("list", length(out_pre_manual))
names(out) <- names(runtimes) <- names(out_pre_manual)

for (i in 1:length(out_pre_manual)) {
  if (!is_FlowCAP[i]) {
    runtimes[[i]] <- system.time({
      # note: In the current version of FlowSOM, the recommended function 
      # FlowSOM::metaClustering_consensus() does not pass along the seed argument 
      # correctly, so results are not reproducible. To get around this, we use the 
      # dependency function ConsensusClusterPlus::ConsensusClusterPlus() instead. 
      # However, this will be fixed in the next update of FlowSOM (version 1.5); after 
      # the update the following (simpler) line of code can be used instead.
      #meta <- FlowSOM::metaClustering_consensus(out_pre_manual[[i]]$map$codes, k = k[[i]], seed = seed)
      
      meta <- suppressMessages(
        ConsensusClusterPlus::ConsensusClusterPlus(t(out_pre_manual[[i]]$map$codes), maxK = k[[i]], seed = seed)
      )
      meta <- meta[[k[[i]]]]$consensusClass
    })
    out[[i]] <- meta
    cat("data set", names(data[i]), ": run complete\n")
    
  } else {
    # FlowCAP data sets: run clustering algorithm separately for each sample
    out[[i]] <- runtimes[[i]] <- vector("list", length(data[[i]]))
    names(out[[i]]) <- names(runtimes[[i]]) <- names(data[[i]])
    
    for (j in 1:length(data[[i]])) {
      runtimes[[i]][[j]] <- system.time({
        # note: In the current version of FlowSOM, the recommended function 
        # FlowSOM::metaClustering_consensus() does not pass along the seed argument 
        # correctly, so results are not reproducible. To get around this, we use the 
        # dependency function ConsensusClusterPlus::ConsensusClusterPlus() instead. 
        # However, this will be fixed in the next update of FlowSOM (version 1.5); after 
        # the update the following (simpler) line of code can be used instead.
        #meta <- FlowSOM::metaClustering_consensus(out_pre_manual[[i]][[j]]$map$codes, k = k[[i]], seed = seed)
        
        meta <- suppressMessages(
          ConsensusClusterPlus::ConsensusClusterPlus(t(out_pre_manual[[i]][[j]]$map$codes), maxK = k[[i]], seed = seed)
        )
        meta <- meta[[k[[i]]]]$consensusClass
      })
      out[[i]][[j]] <- meta
    }
    cat("data set", names(data[i]), ": run complete\n")
    
    # FlowCAP data sets: sum runtimes over samples
    runtimes_i <- do.call(rbind, runtimes[[i]])[, 1:3]
    runtimes_i <- colSums(runtimes_i)
    names(runtimes_i) <- c("user", "system", "elapsed")
    runtimes[[i]] <- runtimes_i
  }
}

# combine runtimes
for (i in 1:length(runtimes)) {
  runtimes[[i]] <- runtimes_pre_manual[[i]] + runtimes[[i]]
}


# extract cluster labels
clus <- vector("list", length(data))
names(clus) <- names(data)

for (i in 1:length(clus)) {
  if (!is_FlowCAP[i]) {
    clus[[i]] <- out[[i]][out_pre_manual[[i]]$map$mapping[, 1]]
    
  } else {
    # FlowCAP data sets
    clus_list_i <- vector("list", length(out_pre_manual[[i]]))
    names(clus_list_i) <- names(out_pre_manual[[i]])
    for (j in 1:length(clus_list_i)) {
      clus_list_i[[j]] <- out[[i]][[j]][out_pre_manual[[i]][[j]]$map$mapping[, 1]]
    }
    
    # convert FlowCAP cluster labels into format "sample_number"_"cluster_number"
    # e.g. sample 1, cluster 3 -> cluster label 1.3
    names_i <- rep(names(clus_list_i), times = sapply(clus_list_i, length))
    clus_collapse_i <- unlist(clus_list_i, use.names = FALSE)
    clus[[i]] <- paste(names_i, clus_collapse_i, sep = "_")
  }
}

sapply(clus, length)


###########################################################################################
### compare with manual gating results ###
###########################################################################################

# helper functions to match clusters and evaluate
source("./cytometry-clustering-comparison-master/helpers/helper_match_evaluate_multiple.R")
#source("../helpers/helper_match_evaluate_single.R")
#source("../helpers/helper_match_evaluate_FlowCAP.R")
#source("../helpers/helper_match_evaluate_FlowCAP_alternate.R")

DATA_DIR <- "./FlowRepository_FR-FCM-ZZPH_files"

# which data sets required subsampling for this method (see parameters spreadsheet)
is_subsampled <- c(FALSE, FALSE, FALSE, FALSE)

# alternate FlowCAP results at the end
is_rare    <- c(FALSE, FALSE, FALSE, FALSE)
is_FlowCAP <- c(FALSE, FALSE, FALSE, FALSE)
n_FlowCAP <- 2

### load truth (manual gating population labels) ###

# files with true population labels (subsampled labels if subsampling was required for
# this method; see parameters spreadsheet)

files_truth <- list(
  Levine_32dim = file.path(DATA_DIR, "Levine_32dim_notransform.fcs"), 
  Levine_13dim = file.path(DATA_DIR, "Levine_13dim_notransform.fcs"), 
  Samusik_01   = file.path(DATA_DIR, "Samusik_01_notransform.fcs"), 
  Samusik_all  = file.path(DATA_DIR, "Samusik_all_notransform.fcs")
)

# extract true population labels

clus_truth <- vector("list", length(files_truth))
names(clus_truth) <- names(files_truth)

for (i in 1:length(clus_truth)) {
  if (!is_subsampled[i]) {
    data_truth_i <- flowCore::exprs(flowCore::read.FCS(files_truth[[i]], transformation = FALSE, truncate_max_range = FALSE))
  } else {
    data_truth_i <- read.table(files_truth[[i]], header = TRUE, stringsAsFactors = FALSE)
  }
  clus_truth[[i]] <- data_truth_i[, "label"]
}

for (i in 1:length(clus)) {
  #clus_truth[[i]] <- clus_truth[[i]][!clus_truth[[i]]=="NaN"]
  #clus[[i]] <- clus[[i]][which(clus_truth[[i]]!="NaN")]
  print(helper_match_evaluate_multiple(clus[[i]],clus_truth[[i]])$mean_F1)
}



# ####### tsne ########
# #xdata <- flowCore::exprs(flowCore::read.FCS(files[[1]], transformation = FALSE, truncate_max_range = FALSE))
# xdata <- flowCore::exprs(flowCore::read.FCS("/stornext/Bioinf/data/lab_speed/zonglun/test_cytofkit/FlowRepository_FR-FCM-ZZPH_files/Levine_32dim.fcs", transformation = FALSE, truncate_max_range = FALSE))
# xdata <- xdata[,marker_cols[[1]]]
# #xdata <- ifelse(xdata<0, 0, xdata)
# # fSOM <- FlowSOM::ReadInput(xdata, transform = TRUE, toTransform = 1:dim(xdata)[2],
# #                            transformFunction = flowCore::logicleTransform(w = 0.1, m=4.5, t = 5000, a=0),
# #                            scale = FALSE)
# # xdata <- fSOM$data
# #xdata <- asinh(xdata/0.5)
# set.seed(1000)
# data_sub_index <- sample(1:dim(xdata)[1],10000) 
# data_sub <- xdata[data_sub_index,]
# data_sub <- as.matrix(data_sub)
# dups <- duplicated(data_sub)
# data_sub <- data_sub[!dups, ]
# set.seed(1000)
# out_Rtsne <- Rtsne(data_sub, pca = FALSE, verbose = TRUE)
# 
# # prepare Rtsne output data for plot
# data_plot <- as.data.frame(out_Rtsne$Y)
# colnames(data_plot) <- c("tSNE_1", "tSNE_2")
# 
# cluster_sub <- clus[[1]][data_sub_index][!dups]
# data_plot[, "cluster"] <- as.factor(cluster_sub)
# 
# ggplot(data_plot, aes(x = tSNE_1, y = tSNE_2, color = cluster)) + 
#   geom_point(size = 0.5) + 
#   coord_fixed(ratio = 1) + 
#   ggtitle("t-SNE projection with FlowSOM clustering") + 
#   theme_bw()
