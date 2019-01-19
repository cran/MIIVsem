#' provides factor scores for miiv model
#'@keywords internal
factorScores <- function(model, data, lv_scores, lv_estimator){
  
  if(lv_estimator == "miiv"){
    
      fit        <- miive(model, data, var.cov = TRUE, overid.method = "stepwise.R2", overid.degree = 1)
      pe         <- estimatesTable(fit, sarg = TRUE)
      mod.syntax <- fit$v$model
    
      fs <- lavaan::lavPredict(
        lavaan::cfa(mod.syntax, data), 
        type = "lv", 
        method = lv_scores
      )
    
  } else if (lv_estimator == "pml"){
    
      fit <- lavaan::cfa(model, data)
      pe  <- lavaan::parameterEstimates(fit)
    
      fs <- lavaan::lavPredict(
        fit, 
        type = "lv", 
        method = lv_scores
      )
      
  } else if (lv_estimator == "svd"){
    
    
    pt <- lavaan::lavParTable(model)
    lv <- unique(pt[pt$op == "=~", "lhs"])
    
    fs <- matrix(0, length(lv), nrow(data))
    
    for (k in 1:length(lv)){
      Y  <- t(data[,pt[pt$op == "=~" & pt$lhs == lv[k], "rhs"]]) # selecting ROIs for that component
      vv <- svd(stats::cov(t(Y)))     # SVD is used, which on covariance matrix should be similar to PCA
      Lam <- vv$u[,1]                 # take the first component from decomposition
      fs[k,] <- t(Lam)%*%Y            # calculate scores across time
    }
    fs <- t(fs)
    colnames(fs) <- lv
    
    
  } else {
    
    stop(paste0("lvgimme error:", lv_estimator, " not supported"))
    
  }

  return(list(
    "fs" = as.data.frame(fs),
    "pe" = pe
    ))

}