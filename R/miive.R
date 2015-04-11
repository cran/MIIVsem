#' Model-implied instrumental variable (MIIV) estimation
#'
#' Estimate SEM models using MIIV-2SLS.
#'
#' @param model A model specified using lavaan model syntax. See the \code{model} argument within the \code{\link[lavaan]{lavaanify}} function for more information.
#' @param data A data frame, list or environment or an object coercible by as.data.frame to data frame.
#' @param instruments A user-supplied list of valid MIIVs for each equation. See Example 2 below. 
#' @param overid A user-specified degree of overidentification (\code{overid}). See Example 3 below. 
#'
#' @return model
#' @return dat
#' @return modeeqns
#' 
#' @details 
#' \itemize{
#'  \item{\code{overid}} {If the user-specified degree of overidentification (\code{overid}) exceeds the number of available MIIVs for a given equation, the maximum number of MIIVs will be used instead.  In this case, the \code{df} column for any equations in which the degrees of freedom for the \code{Sargan} test are less than the \code{overid} value will be appended with an \code{*}. A note willalso  be displayed to alert the user to the changes. In the example below, the \code{overid} parameter is set to 2, however the \code{y1} equation has only 1 additional MIIV avaialable.}
#' \item{\code{instruments}} {Using the \code{instruments} option you can specify the MIIVs directly for each equation in the model.  To utilize this option you must first define a list of instruments using the syntax displayed below. After the list is defined, set the \code{instruments} argument equal to the name of the list of MIIVs. Note, \code{instruments} are specified for an equation, and not for a specific endogenous variable.}
#' }
#' 
#' @references 
#' 
#' Bollen, K. A. 1996.	An	Alternative	2SLS Estimator	for	Latent	
#' Variable	Models.	\emph{Psychometrika}, 61, 109-121.
#' 
#' Bollen,	K. A. 2001.	Two-stage	Least	Squares	and	Latent	Variable	Models:	
#' Simultaneous	Estimation	and	Robustness	to	Misspecifications.
#' In	R.	Cudeck,	S.	Du	Toit,	and	D.	Sorbom	(Eds.),	Structural	
#' Equation	Modeling:	Present	and	Future,	A	Festschrift	in	Honor	of	Karl	
#' Joreskog	(pp. 119-138).	Lincoln,	IL: Scientific	Software.
#' 	
#'
#' @examples
#' 
#' # Example 1
#' 
#'  bollen1989a_model <- '
#'
#'    Eta1 =~ y1 + y2  + y3  + y4  
#'    Eta2 =~ y5 + y6  + y7  + y8    
#'    Xi1  =~ x1 + x2 + x3 
#'
#'    Eta1 ~ Xi1  
#'    Eta2 ~ Xi1 
#'    Eta2 ~ Eta1 
#'
#'    y1   ~~ y5
#'    y2   ~~ y4
#'    y2   ~~ y6
#'    y3   ~~ y7
#'    y4   ~~ y8
#'    y6   ~~ y8 
#'  '
#'  
#'   miive(model = bollen1989a_model, data = bollen1989a)
#'  
#'  
#' # Example 2
#' 
#'   my_instruments <- ' 
#'    y1 ~ x2 + x3                            
#'    y5 ~ y2 + y3 + y4 + x2                
#'    y2 ~ y3 + y7 + y8 + x2           
#'    y3 ~ y2 + y4 + y6 + y8        
#'    y4 ~ y3 + y6           
#'    y6 ~ y3 + y4 + y7 + x2            
#'    y7 ~ y2 + y4 + y6 + y8       
#'    y8 ~ y2 + y3 + y7 + x2          
#'    x2 ~ y1 + y5 + y2 + y3 + y4 + y6
#'    x3 ~ y1 + y5 + y2 + y3 + y4 + y6
#'  '
#'  
#' miive(model = bollen1989a_model, data = bollen1989a, 
#'       instruments = my_instruments)
#'  
#'  
#' # Example 3
#'  
#'  miive(model = bollen1989a_model, data = bollen1989a, overid = 2)
#'  
#'  
#' @export
miive<- function(model = model, data = data, instruments = NULL, overid = NULL){

  model  <- miivs(model)$eqns
  digits <- 3
  
 if ( !is.null(overid) && !is.null(instruments) )  {
        stop(paste("Cannot supply both instruments list and overid."))}
  
  if (!is.null(instruments)){
    table <- lavParTable(instruments)
    table <- table[table$op == "~", ]
    dv_list <- unique(table$lhs)
    iv_list <- list(DV_user = "", IV_user = "")
    iv_list <- replicate(length(dv_list), iv_list, simplify = FALSE)
  
    for (i in 1:length(dv_list)){
      iv_list[[i]]$DV_user <- dv_list[i]
      iv_list[[i]]$IV_user <- c(table[table$lhs == dv_list[i],]$rhs)
    }
  
    dv_user <- unlist(lapply(iv_list, "[", c("DV_user")), use.names=FALSE)
    iv_user <- unlist(lapply(iv_list, "[", c("IV_user")), use.names=FALSE)
  
    for (i in 1:length(model)){
    
      dv <- model[[i]]$DV
      index <- which(dv_user == dv)
    
      if (is.integer(index) && length(index) == 0L) {
        stop(paste("Cannot find instruments for ", dv))}
      
      iv_user  <- iv_list[[index]]$IV_user
      iv_miivs <- model[[i]]$IV
      n_pred   <- length(model[[i]]$P)
    
      if (length(iv_user) < n_pred) {
        stop(paste("Need at least ", n_pred," instruments for ", dv))}
    
    
      check <- iv_user[which(!(iv_user%in%iv_miivs))]
    
      if ( !(is.character(check) && length(check) == 0L) ) {
        stop(paste("Instruments for ", dv," are not valid."))}
  
      model[[i]]$IV3 <- iv_user
    
    }
    
  } #
  
  optim <- function(model, overid, data){
    for (i in 1:length(model)){
      y  <- as.matrix( cbind(data[,model[[i]]$DV] ) )
      P  <- model[[i]]$P
      k  <- overid + length(model[[i]]$P) # try absolute next
      
      if (k > length(model[[i]]$IV) ) {
        model[[i]]$IV2 <- model[[i]]$IV
        model[[i]]$NOTE <- 
          paste("* Maximum number of MIIVs is less than requested degree of 
                overidentification.\n  See df for degree of 
                overidentification.", sep="")
        model[[i]]$MSG <- "*"
        }

      else if (k <= length(model[[i]]$IV) ) {
      cd <- cor(data)
      cd <- cbind(cd[,P])
      cd[cd == 1] <- 0
      cd_m <- as.matrix(cd)
    
      cd_m <- cbind( apply(cd_m, 1, max) ) # get largest val from each row 
      ord_names <- rownames(cd_m)[order(cd_m, # return col vec
                                        decreasing=TRUE)][1:nrow(cd_m)]
    
      temp <- model[[i]]$IV
      temp <- temp[order(match(temp,ord_names))]
      model[[i]]$IV2 <- temp[1:k]
      }
      }
    return(model)
    }



  if (!is.null(overid)){ model <- optim(model, overid, data)}
 
  if (is.null(overid)){ 
    for (i in 1:length(model)){
      model[[i]]$IV2 <- model[[i]]$IV 
    }
  }
  
   if (!is.null(instruments)){ 
    for (i in 1:length(model)){
      model[[i]]$IV2 <- model[[i]]$IV3 
    }
  }
  
  dat    <- data.frame()
  
  for (i in 1:length(model)){
    y      <- as.matrix( cbind(data[,model[[i]]$DV] ) )
    H      <- as.matrix( cbind(1, data[,model[[i]]$P] ) )
    Z      <- as.matrix( cbind(1, data[,model[[i]]$IV2] ) )
    ZtZinv <- solve(crossprod(Z))
    ZtH    <- crossprod(Z, H)
    ZZZZH  <- Z %*% ZtZinv %*% ZtH
    ZZZinv <- solve(crossprod(ZZZZH))
    b      <- solve(crossprod(ZZZZH)) %*% t(ZZZZH) %*% y
    r      <- y - H  %*% b
    s2     <- sum(r^2)/(nrow(data) - ncol(H)) 
    cov    <- s2 * ZZZinv
    s      <- cbind(sqrt(diag(cov)))
    IVID   <- solve(crossprod(Z)) %*% crossprod(Z,r)
    SST    <- sum((r-mean(r))^2)
    SSE    <- sum((r - (Z %*% IVID))^2)
    Rsq    <- 1 - (SSE/SST)
    sarg   <- length(r)*Rsq
    length(sarg) <- nrow(s)
    sargdf <- ncol(Z)-ncol(H)
    length(sargdf) <- nrow(s)
    sargp  <- 1 - pchisq(sarg,df = sargdf) 
    length(sargp) <- nrow(s)
    temp <- as.data.frame(cbind(b,s,sarg,sargp))
    temp$t <- temp[,1]/temp[,2]
    temp$p <- 2* (pnorm(abs(temp$t),,lower.tail=FALSE))
    temp$p <- as.numeric(as.character(temp$p))
    temp$i <- cbind(c("Int", model[[i]]$P2))
    temp$d <- model[[i]]$DV
    
    
      if (model[[i]]$NOTE != "") {
        sargdf <- paste(sargdf, "*", sep="")
        temp$df <-  sargdf
        temp$df[-1] <- NA
        
      }
    
      if (model[[i]]$NOTE == "") {
        temp$df <- sargdf
      }
    
      if (model[[i]]$EQ == "y1") { 
        temp$i <- cbind(c("Int", model[[i]]$P2))
        temp$d <- model[[i]]$DV2
      }

    if (i == 1) {dat <- temp }
    if (i >  1) {dat <- rbind(dat,temp) }
  }
  colnames(dat) <- c("Estimate", "StdErr","Sargan",
                     "P(Chi)", "Z", "P(|Z|)", "EV", "DV", "df")
  dat <- dat[, c("DV", "EV", "Estimate", "StdErr", 
                 "Z", "P(|Z|)","Sargan", "df", "P(Chi)")]
  

   for (i in 1:length(model)){
    LHS <- paste(model[[i]]$DV, collapse = ", ")
    RHS <- paste(model[[i]]$P, collapse = ", ")
    Instruments <- paste(model[[i]]$IV2, collapse = ", ")
    Disturbance <- paste("e.",model[[i]]$C, collapse = ", ", sep="")
    Effects <- paste(model[[i]]$EF, collapse = ", ") 
    modtemp <- as.data.frame(cbind(LHS, RHS, Disturbance, Instruments))
    colnames(modtemp) <- c("LHS", "RHS", "Composite Disturbance", "MIIVs")
    if (i == 1) {modeqns <- modtemp }
    if (i >  1) {modeqns <- rbind(modeqns,modtemp) }
  }
  
  res <- list(model = model, dat = dat, modeqns = modeqns)
  class(res) <- "miive"
  res
}
