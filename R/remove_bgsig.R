#' Wrap for subtracting the background to estimate the experimentally generated signature
#' @param MutCatalogue 96-channel mutational catalog
#' @param bg_column the column name of the background catalog
#' @param ko_column the column names of the experiment catalog
#' @param sampling_number number of bootstrapping samples
#' @param start_num mutation burden to start
#' @param boundary Range of signature around centroid (default = 2)
#' @param outputname output file name
#' @export
Wrap_KOSig <- function(MutCatalogue, bg_column, ko_column, sampling_number, start_num, boundary, outputname){
  
  KOSig <- RemoveBackground_vector_single(MutCatalogue[,bg_column], MutCatalogue[,ko_column], sampling_number, start_num, boundary)
  KOSig$MutationType <- MutCatalogue[,"MutationType"]
  utils::write.table(KOSig, paste0(outputname, ".txt"), sep = "\t", col.names = T, row.names = F, quote = F)
  
}


#' Return single vector
#' @param background_profile background mutational catalog
#' @param sig_profile experimentally generated mutational catalog
#' @param sampling_number number of bootstrapping samples
#' @param start_num mutation burden to start
#' @param boundary Range of signature around centroid (default = 2)
#' @return data.frame including background signature and experiment signature
#' @export 
RemoveBackground_vector_single <- function(background_profile, sig_profile, sampling_number, start_num, boundary = 4){
  
  #Calculate number of mutation channels
  nTypes <- length(sig_profile)
  
  # Remove weak mutation types in sig_profile
  
  #set threshold value for weak mut types?
  removeWeakMutationTypes <- 0.01
  
  #calculate 1% of the total mutations in each channel ..?
  genomesOriginal <- as.data.frame(sig_profile) 
  Totalmutations <- sum(sig_profile)
  removeMutations_max <- removeWeakMutationTypes * Totalmutations 
  
  #Identify set of mutation types to remove - classes with cumulative sum (?) less than lower of 4 or 1% of ?? should be removed
  removeIdx <- which(cumsum(rowSums(genomesOriginal[,1:dim(genomesOriginal)[2], drop = FALSE])[order(rowSums(genomesOriginal[,1:dim(genomesOriginal)[2], drop = FALSE]))]) <= min(removeMutations_max, 4))
  
  #Remove any channels to be removed as IDd above
  if(length(removeIdx) > 0) {
    mutationTypesToRemoveSet <- order(rowSums(genomesOriginal[,1:dim(genomesOriginal)[2],drop = FALSE]))[removeIdx]
    genomesReducted <- as.data.frame(genomesOriginal[-mutationTypesToRemoveSet,])
    reducedMutationtypes <- dim(genomesReducted)[1] #Calculate new, reduced number of mutation classes
    } else { #else, return the original set
      mutationTypesToRemoveSet <- dim(genomesOriginal)[1]+1
      genomesReducted <- genomesOriginal
    }
  
  # Bootstrap signature profile ------
  
  #calculate av. mut profile across all samples?
  centroid_sig <- rowMeans(genomesReducted[,1:dim(genomesReducted)[2], drop = FALSE])
  
  #create matrix repeating centroid sig 'sampling_number' times
  RepSig <- matrix(rep(centroid_sig, sampling_number), ncol = sampling_number)
  
  #generate bootstrapped samples of a mutational catalogue, using the repeated centroid sigs as the catalogue, and the total mut sum of the centroid sig as the number of bootstraps?
  Sig_bootstraps <- bootstrapGenomesfun2(RepSig, sum(centroid_sig))
  
  
  i = start_num #input param? how is this decided? starting number of muts to include in sig profile?
  reachLimit <- FALSE
  diff_all_save <- NULL
  
  #While i is less than total mut sum of the centroid sig and limit not reached:
  while(i < sum(centroid_sig) & !reachLimit) {
    
    # Remove the same weak mutation types in background profile
    backgroundReducted <- background_profile[-mutationTypesToRemoveSet]
    RepControl <- matrix(rep(backgroundReducted, sampling_number), ncol = sampling_number)
    bg_bootstraps <- bootstrapGenomesfun2(RepControl, i)
    
    # Compute range of background profile
    centroid_background <- rowMeans(bg_bootstraps)
    sd_background <- apply(bg_bootstraps, 1, stats::sd)
    boundary_background <- centroid_background + boundary * sd_background
    
    # Compute range of signature
    centroid_sig <- rowMeans(Sig_bootstraps)
    sd_sig <- apply(Sig_bootstraps, 1, stats::sd)
    boundary_sig <- centroid_sig + boundary * sd_sig
    
    # Compute diff between boundary sig (range of input sig around its centroid) and the background centroid
    diff_all_boundary <- boundary_sig - centroid_background
    diff_all <- centroid_sig - centroid_background
    
    #i'm not so clear on what exactly is being done and why from this point on?
    if(length(which(diff_all_boundary < 0)) > 0) {
      reachLimit <- TRUE
    }
    
    if(length(which(diff_all_boundary < 0)) == 0) {
      diff_all_boundary_save <- diff_all_boundary
      diff_all_save <- diff_all
      diff_all_save[which(diff_all_save < 0)] <- 0
    }
    
    
    
    i = i + 1
    
  }
  
  if(length(diff_all_save) == 0) {
    stop("You need to reduce the start_number! ", 
         "Exiting...",call.=FALSE)
  }
  
  # Add Weak mutations for KO expoure
  exposure <- rep(0, nTypes)
  origArrayIndex <- 1
  for(i in 1:nTypes) {
    if(! i %in% mutationTypesToRemoveSet) {
      exposure[i] <- diff_all_save[origArrayIndex]
      origArrayIndex = origArrayIndex + 1
    }
  }
  
  # Add Weak mutations for background
  background_exposure<- rep(0, nTypes)
  origArrayIndex <- 1
  for(i in 1:nTypes) {
    if(! i %in% mutationTypesToRemoveSet) {
      background_exposure[i] <- centroid_background[origArrayIndex]
      origArrayIndex = origArrayIndex + 1
    }
  }
  
  
  return(data.frame("KO_exposure" = exposure, "background_exposure" = background_exposure))
  
  
}

#' Wrap for subtracting the background to estimate the experimentally generated signature
#' @param genomes mutational catalog for bootstrapping
#' @param n number of bootstrapping to draw
#' @return a matrix of bootstrapping sample catalogs
#' @export
bootstrapGenomesfun2 <- function(genomes, n){ #Generates bootstrapped samples of a mutational catalogue (genomes)
  
  #For each catalogue, generate multinomial sample of size n, with probabilities defined by mutational frequencies of input
  return(apply(genomes, 2, function(x) stats::rmultinom(1, n, x)))
}


