#' Calculate group and component PIPs for a fitted gsbart model

#' Computes posterior inclusion probabilities (PIPs) from a fitted \code{gsbart} regression object. The function can return PIPs at the group level and, when applicable, at the component or variable level.
#' 
#' @param fit A fitted \code{gsbart} model object. For continuous outcomes, this is typically the object returned by \code{gsbart_regression()}.
#' @param group Logical indicating whether group was specified in the fitted model. If \code{TRUE}, both group and component PIPs are returned.
#' @param vars A matrix or data frame of predictor variables used to fit the soft BART model.
#'

#' @return A list, vector, matrix, or data frame containing calculated group and/or component posterior inclusion probabilities, depending on the structure of \code{fit}, \code{group}, and \code{vars}.

#' @export
#' 
calPIPs_gsbart <- function(fit, group = NULL, vars = NULL){
  
  if(group){
    group_counts_mat <- as.matrix(do.call(rbind, lapply(1:length(fit$tree_counts), 
                                                        \(x) as.numeric(rowSums(fit$tree_counts[[x]]$Group_tree_counts)))))
    
    var_counts_mat <- fit$var_counts
    Grp_PIPs <- colMeans(group_counts_mat > 0)
    
    ##### conditional PIPs
    # create group-to-variable mapping
    num_vars = ncol(vars)
    gbart_hypers <- fit$hypers
    grps <- gbart_hypers$group + 1
    
    group_to_var_map <- lapply(1:max(grps), \(x) which(grps == x))
    
    Cond_PIPs_mat <- matrix(ncol = num_vars)
    for (grp_idx in 1:length(group_to_var_map)) {
      vars_idx <- group_to_var_map[[grp_idx]]
      
      reduced_mat <- var_counts_mat[which(group_counts_mat[,grp_idx] > 0), vars_idx]
      Cond_PIPs_mat[,vars_idx] <- colMeans(reduced_mat > 0)
    }
    
    if(is.null(colnames(vars))){
      colnames(Cond_PIPs_mat) <- paste0("Z", 1:num_vars)
    }else{colnames(Cond_PIPs_mat) <- colnames(vars)}
    
    ## extract the column names for easy readability.
    with_grp_PIPs <- as.data.frame(do.call(rbind, lapply(1:length(group_to_var_map), function(idx){
      names <- colnames(Cond_PIPs_mat)[group_to_var_map[[idx]]]
      values <- format(Cond_PIPs_mat[,names], nsmall = 8)
      
      grp_idx <- idx
      grp_values <- format(Grp_PIPs[idx], nsmall = 3)
      
      return(cbind(names, grp_idx, values, grp_values))
    })))
    rownames(with_grp_PIPs) <- NULL
    with_grp_PIPs_df <- data.frame("Variables"   = with_grp_PIPs[, "names"], 
                                   "Group"       = with_grp_PIPs[, "grp_idx"], 
                                   "GroupPIPs"   = with_grp_PIPs[, "grp_values"], 
                                   "CondPIPs"    = with_grp_PIPs[, "values"])
    return(with_grp_PIPs_df)
    
  }else{
    var_counts_mat <- fit$var_counts
    if(is.null(colnames(vars))){
      colnames(var_counts_mat) <- paste0("Z", 1:num_vars)
    }else{colnames(var_counts_mat) <- colnames(vars)}
    
    PIPs_df <- data.frame("Variables" = colnames(var_counts_mat),
                          "PIPs" = colMeans(var_counts_mat > 0))
    return(PIPs_df)
  }
}

