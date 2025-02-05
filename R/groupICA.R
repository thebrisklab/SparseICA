#' Create a List of fMRI Files for Group ICA Analysis
#'
#' This function scans a BIDS-formatted directory for subject-specific fMRI files that match a specified pattern 
#' and returns a list of these files for use in group ICA analysis.
#'
#' @param bids_path A character string specifying the path to the root directory of the BIDS-formatted dataset. This directory should contain subject folders (e.g., \code{sub-*}).
#' @param pattern A character string specifying the pattern to match fMRI files. The default is \code{"task-rest.*\\.dtseries\\.nii$"}.
#'
#' @return A named list where each element corresponds to a subject directory and contains a vector of matched fMRI file paths. The names of the list are the subject IDs.
#'
#' @examples
#' # Example usage:
#' # Assuming `bids_dir` is the path to a BIDS dataset:
#' # group_list <- create_group_list(bids_path = bids_dir, pattern = "task-rest.*\\.dtseries\\.nii$")
#' # Print the structure of the list:
#' # str(group_list)
#'
#' @export
create_group_list <- function(bids_path,pattern = "task-rest.*\\.dtseries\\.nii$") {
  
  # Validate the BIDS path
  if (!file.exists(bids_path) || !file.info(bids_path)$isdir) {
    stop("The specified BIDS path does not exist or is not a directory.")
  }
  
  # Get all subject directories in the BIDS folder
  subject_dirs <- list.dirs(bids_path, full.names = FALSE, recursive = FALSE)
  subject_dirs <- subject_dirs[grepl("sub-[a-zA-Z0-9]+$", basename(subject_dirs))]
  
  if (length(subject_dirs) == 0) {
    stop("No subject directories found in the specified BIDS path.")
  }
  
  # Initialize a list to store files for each subject
  result <- list()
  
  # Loop through each subject directory and search for files
  for (subject_dir in subject_dirs) {
    subject_id <- basename(subject_dir) # Extract subject ID
    
    # Find all files matching the pattern in the subject directory
    matched_files <- list.files(paste0(bids_path,"/",subject_dir), pattern = pattern, recursive = TRUE, full.names = FALSE)
    
    # Store the matched files in the result list under the subject ID
    result[[subject_id]] <- matched_files
  }
  
  return(result)
}


#' Generate Group-Level Principal Components (PCs) for fMRI Data
#'
#' This function computes subject-level principal components (PCs) from fMRI data and performs a group-level PCA for dimension reduction, designed for cortical surface fMRI data in BIDS format.
#'
#' @param bids_path A character string specifying the root directory of the BIDS-formatted dataset. 
#' @param subj_list A named list generated from \code{create_group_list} containing fMRI file paths for each subject. 
#' @param n.comp An integer specifying the number of components to retain during group-level PCA. Default is 30.
#' @param ncore An integer specifying the number of cores to use for parallel processing. Default is 1.
#' @param npc An integer specifying the number of components to retain during subject-level PCA. Default is 85.
#' @param iter_std An integer specifying the number of iterative standardization steps to apply to fMRI data. Default is 5.
#' @param brainstructures A character vector specifying the brain structures to include in the analysis. Options are \code{"left"} (left cortex), \code{"right"} (right cortex), and/or \code{"subcortical"} (subcortex and cerebellum). Can also be \code{"all"} (obtain all three brain structures). Default is \code{c("left", "right")}.
#' @param verbose A logical value indicating whether to print convergence information during execution. Default is \code{TRUE}.
#'
#' @return A numeric matrix containing the group-level principal components, with dimensions determined by the number of retained components (\code{n.comp}) and the concatenated data across all subjects.
#'
#' @details
#' NOTE: This function requires the \code{ciftiTools} package to be installed, and set up the path to the Connectome Workbench folder by \code{ciftiTools.setOption()}. See the package \code{ciftiTools} documentation for more information.
#'
#' @import ciftiTools
#' @import irlba
#' @import parallel
#' @export
gen_groupPC <- function(bids_path, subj_list, n.comp = 30, ncore=1,
                        npc = 85, iter_std = 5, brainstructures = c("left", "right"),verbose = TRUE) {
  
  if (verbose) message("# Start generating group PC of", length(subj_list), "subjects. #")
  
  process_subject <- function(subject_name, subject_files) {
    nfile <- length(subject_files)
    if (verbose) message("##", subject_name, "has", nfile, "cortical surface fMRI data. ##")
    
    dat <- c()
    for (j in 1:nfile) {
      xii <- read_cifti(cifti_fname = paste0(bids_path, "/", subject_name, "/", subject_files[j]),
                        brainstructures = brainstructures)
      xii_mat <- as.matrix(xii)
      
      # Iterative standardization
      for (k in 1:iter_std) {
        # Standardize across time
        xii_mat <- apply(xii_mat, 2, function(col) {
          (col - mean(col)) / sd(col)
        })
        # Standardize across voxel/vertex
        xii_mat <- t(apply(xii_mat, 1, function(row) {
          (row - mean(row)) / sd(row)
        }))
      }
      dat <- cbind(dat, xii_mat)
    }
    if (verbose) message("## Total number of TRs:", dim(dat)[2], ".##")
    
    # Perform PCA
    subj_PCA <- prcomp_irlba(dat, npc)
    PC_subj <- subj_PCA$x
    dimnames(PC_subj) <- NULL
    
    if (verbose) message("## Retained", npc, "PCs. ##")
    return(PC_subj)
  }
  
  # Parallel processing of subjects
  #cores <- detectCores() - 1  # Reserve one core for the system
  groupPC_list <- mclapply(names(subj_list), function(subject_name) {
    process_subject(subject_name, subj_list[[subject_name]])
  }, mc.cores = ncore)
  
  # Combine all subject PCs
  groupPC <- do.call(cbind, groupPC_list)
  if (verbose) message("# Finish subject-level PCA! The concatenated matrix has dimension", dim(groupPC), ". #")
  
  # Perform group PCA
  temp <- whitener(X = groupPC, n.comp = n.comp, use_irlba = TRUE)
  groupPC <- temp$Z
  
  if (verbose) message("# Finish group PCA. #")
  return(groupPC)
}


#' Perform Group Sparse Independent Component Analysis (Sparse ICA)
#'
#' This function performs Sparse ICA on group-level fMRI data. It processes BIDS-formatted fMRI datasets, performs PCA to reduce dimensionality, selects a tuning parameter \code{nu} (optionally using a BIC-like criterion), and executes Sparse ICA to estimate independent components.
#'
#' @param bids_path A character string specifying the root directory of the BIDS-formatted dataset.
#' @param subj_list A named list where each element corresponds to a subject and contains vectors of fMRI file names. If \code{NULL}, the subject list is generated automatically using \code{\link{create_group_list}}. Default is \code{NULL}.
#' @param nu A numeric value for the tuning parameter, or \code{"BIC"} to select \code{nu} using a BIC-like criterion. Default is \code{"BIC"}.
#' @param n.comp An integer specifying the number of components to estimate. Default is 30.
#' @param ncore An integer specifying the number of cores to use for parallel processing. Default is 1.
#' @param method A character string specifying the computation method for Sparse ICA. Options are \code{"C"} (default) for C-based computation or \code{"R"} for R-based computation.
#' @param npc An integer specifying the number of components to retain during subject-level PCA. Default is 85.
#' @param iter_std An integer specifying the number of iterative standardization steps to apply to fMRI data. Default is 5.
#' @param brainstructures A character vector specifying the brain structures to include in the analysis. Options are \code{"left"} (left cortex), \code{"right"} (right cortex), and/or \code{"subcortical"} (subcortex and cerebellum). Can also be \code{"all"} (obtain all three brain structures). Default is \code{c("left", "right")}.
#' @param restarts An integer specifying the number of random initializations for Sparse ICA. Default is 40.
#' @param positive_skewness A logical value indicating whether to enforce positive skewness on the estimated components. Default is \code{TRUE}.
#' @param use_irlba A logical value indicating whether to use the \code{irlba} method for fast truncated Singular Value Decomposition (SVD) during whitening. This can improve memory efficiency for intermediate datasets. Default is \code{TRUE}.
#' @param eps A numeric value specifying the convergence threshold. Default is \code{1e-6}.
#' @param maxit An integer specifying the maximum number of iterations for Sparse ICA. Default is 500.
#' @param BIC_plot A logical value indicating whether to generate a plot of BIC values for different \code{nu} candidates when selecting \code{nu}. Default is \code{TRUE}.
#' @param nu_list A numeric vector specifying candidate values for \code{nu} when selecting it using a BIC-like criterion. Default is \code{seq(0.1, 4, 0.05)}.
#' @param verbose A logical value indicating whether to print progress messages. Default is \code{TRUE}.
#' @param BIC_verbose A logical value indicating whether to print detailed messages during the BIC-based selection of \code{nu}. Default is \code{FALSE}.
#' @param converge_plot A logical value indicating whether to generate a plot showing the convergence trace during Sparse ICA. Default is \code{FALSE}.
#'
#' @return A list containing the results of the group Sparse ICA analysis, including:
#' \describe{
#'   \item{\code{loglik}}{The minimal log-likelihood value among the random initializations.}
#'   \item{\code{estS}}{A numeric matrix of estimated sparse independent components with dimensions P x Q.}
#'   \item{\code{estU}}{The estimated U matrix with dimensions Q x Q.}
#'   \item{\code{whitener}}{The whitener matrix used for data whitening.}
#'   \item{\code{converge}}{The trace of convergence for the U matrix.}
#'   \item{\code{best_nu}}{The selected \code{nu} value (if \code{nu = "BIC"}).}
#'   \item{\code{BIC}}{A numeric vector of BIC values for each \code{nu} candidate (if \code{nu = "BIC"}).}
#'   \item{\code{nu_list}}{The list of \code{nu} candidates used in the BIC-based selection (if \code{nu = "BIC"}).}
#' }
#'
#' @details
#' The function operates in four main steps:
#' \enumerate{
#'   \item If \code{subj_list} is not provided, it creates a list of subject-specific fMRI files using \code{\link{create_group_list}}.
#'   \item Performs subject-level PCA using \code{\link{gen_groupPC}} to reduce data dimensionality.
#'   \item Selects the tuning parameter \code{nu} using a BIC-like criterion (if \code{nu = "BIC"}) or uses the provided \code{nu}.
#'   \item Executes Sparse ICA on the group-level PCs to estimate independent components.
#' }
#'
#' @seealso \code{\link{create_group_list}}, \code{\link{gen_groupPC}}, \code{\link{BIC_sparseICA}}, \code{\link{sparseICA}}
#' @export
group_sparseICA <- function(bids_path, subj_list = NULL, nu = "BIC",
                            n.comp = 30, method = "C", ncore = 1,
                            npc = 85, iter_std = 5, brainstructures = c("left","right"),
                            restarts = 40, positive_skewness = TRUE,use_irlba = TRUE,
                            eps = 1e-6, maxit = 500, BIC_plot = TRUE, nu_list = seq(0.1,4,0.05),
                            verbose = TRUE, BIC_verbose = FALSE, converge_plot = FALSE){
  
  if(verbose){
    message("##################################")
    message("##### Start group Sparse ICA #####")
    message("##################################")
    message("+++++ Step 1: Create subject list +++++")
  }
  
  if(is.null(subj_list)){
    subj_list <- create_group_list(bids_path)
    if(verbose){
      message("# No input subject list, create using given path to BIDS. #")
      message("## Detect",length(subj_list),"subjects to be included. ##")
    }
  }else{
    message("# Use given subject list. #")
    message("## Detect",length(subj_list),"subjects to be included. ##")
  }
  
  if(verbose){
    message("+++++ Step 2: Perform subject level PCA +++++")
  }
  
  PC_group <- gen_groupPC(bids_path=bids_path,
                          subj_list = subj_list,
                          ncore = ncore,
                          n.comp = n.comp,
                          npc = npc,iter_std = iter_std,brainstructures = brainstructures,verbose = verbose)
  
  if(verbose){
    message("+++++ Step 3: Select the tuning parameter +++++")
  }
  
  if(nu == "BIC"){
    if (verbose) message("# Select nu by BIC. #")
    nu_selection <- BIC_sparseICA(xData = PC_group, n.comp = n.comp, whiten = "none",method = method, 
                                  use_irlba = use_irlba,eps = eps,maxit = maxit,
                                  BIC_plot = BIC_plot,nu_list = nu_list,verbose=BIC_verbose)
    my_nu <- nu_selection$best_nu
    if (verbose) message("## The selected nu is",my_nu,". ##")
  }else{
    my_nu <- nu
    if (verbose) message("# Use given nu = ",my_nu,". #")
  }
  
  if(verbose){
    message("+++++ Step 4: Perform Sparse ICA on group PC +++++")
  }
  
  my_group_sparseICA <- sparseICA(xData = PC_group,
                                  n.comp = n.comp, nu = my_nu,
                                  method = method,
                                  whiten = "none",restarts = restarts,use_irlba = use_irlba,
                                  positive_skewness = positive_skewness,eps = eps,maxit = maxit,
                                  verbose = verbose,converge_plot = converge_plot)
  
  my_group_sparseICA$estM <- NULL
  my_group_sparseICA$best_nu <- nu_selection$best_nu
  my_group_sparseICA$BIC <- nu_selection$BIC
  my_group_sparseICA$nu_list <- nu_selection$nu_list
    
  if(verbose){
    message("###################################")
    message("##### Finish group Sparse ICA #####")
    message("###################################")
  }
  
  return(my_group_sparseICA)
}