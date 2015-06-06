#' Model-implied instrumental variable (MIIV) estimation
#'
#' Estimate SEM models using MIIV-2SLS.
#'
#' @param model A model specified using lavaan model syntax. See the \code{model} argument within the \code{\link[lavaan]{lavaanify}} function for more information.
#' @param data A data frame, list or environment or an object coercible by as.data.frame to data frame.
#' @param instruments A user-supplied list of valid MIIVs for each equation. See Example 2 below. 
#' @param overid A user-specified degree of overidentification (\code{overid}). See Example 3 below. 
#' @param printmiivs A logical indicating whether or not to display MIIVs in output.
#' @param varcov Option for estimating conditional variance and coavariance paramaters. Default is (\code{NULL}).
#'
#'
#' @return model
#' @return dat
#' @return modeeqns
#' 
#' @details 
#' \itemize{
#' \item{\code{overid}} {If the user-specified degree of overidentification (\code{overid}) exceeds the number of available MIIVs for a given equation, the maximum number of MIIVs will be used instead.  In this case, the \code{df} column for any equations in which the degrees of freedom for the \code{Sargan} test are less than the \code{overid} value will be appended with an \code{*}. A note willalso  be displayed to alert the user to the changes. In the example below, the \code{overid} parameter is set to 2, however the \code{y1} equation has only 1 additional MIIV avaialable.}
#' \item{\code{instruments}} {Using the \code{instruments} option you can specify the MIIVs directly for each equation in the model.  To utilize this option you must first define a list of instruments using the syntax displayed below. After the list is defined, set the \code{instruments} argument equal to the name of the list of MIIVs. Note, \code{instruments} are specified for an equation, and not for a specific endogenous variable.}
#'  \item{\code{varcov}} {Currently, only \code{ML} and \code{ULS} options are supported through the \code{\link[lavaan]{lavaan}} package.}
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
miive <- function(model = model, data = data, instruments = NULL, overid = NULL,
                  varcov = NULL, printmiivs = FALSE){
  d  <- miivs(model)$eqns
  digits <- 3
  
 if ( !is.null(overid) && !is.null(instruments) )  {
        stop(paste("Cannot supply both instruments list and overid."))}
  
 if ( length(miivs(model)$constr) == 0 )  { 
      restrictions = FALSE }
 if ( length(miivs(model)$constr) > 0 )  { 
      restrictions = TRUE }
  
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
  
    for (i in 1:length(d)){
    
      dv <- d[[i]]$DV
      index <- which(dv_user == dv)
    
      if (is.integer(index) && length(index) == 0L) {
        stop(paste("Cannot find instruments for ", dv))}
      
      iv_user  <- iv_list[[index]]$IV_user
      iv_miivs <- d[[i]]$IV
      n_pred   <- length(d[[i]]$P)
    
      if (length(iv_user) < n_pred) {
        stop(paste("Need at least ", n_pred," instruments for ", dv))}
    
    
      check <- iv_user[which(!(iv_user%in%iv_miivs))]
    
      if ( !(is.character(check) && length(check) == 0L) ) {
        stop(paste("Instruments for ", dv," are not valid."))}
  
      d[[i]]$IV3 <- iv_user
    
    }
    
  } #
  
  optim <- function(d, overid, data){
    for (i in 1:length(d)){
      y  <- as.matrix( cbind(data[,d[[i]]$DV] ) )
      P  <- d[[i]]$P
      k  <- overid + length(d[[i]]$P) # try absolute next
      
      if (k > length(d[[i]]$IV) ) {
        d[[i]]$IV2 <- d[[i]]$IV
        d[[i]]$NOTE <- 
          paste("* Maximum number of MIIVs is less than requested degree of 
                overidentification.\n  See df for degree of 
                overidentification.", sep="")
        d[[i]]$MSG <- "*"
        }

      else if (k <= length(d[[i]]$IV) ) {
      cd <- cor(data)
      cd <- cbind(cd[,P])
      cd[cd == 1] <- 0
      cd_m <- as.matrix(cd)
    
      cd_m <- cbind( apply(cd_m, 1, max) ) # get largest val from each row 
      ord_names <- rownames(cd_m)[order(cd_m, # return col vec
                                        decreasing=TRUE)][1:nrow(cd_m)]
    
      temp <- d[[i]]$IV
      temp <- temp[order(match(temp,ord_names))]
      d[[i]]$IV2 <- temp[1:k]
      }
      }
    return(d)
    }

  if (!is.null(overid)){ d <- optim(d, overid, data)}
 
  if (is.null(overid)){ 
    for (i in 1:length(d)){
      d[[i]]$IV2 <- d[[i]]$IV 
    }
  }
  
   if (!is.null(instruments)){ 
    for (i in 1:length(d)){
      d[[i]]$IV2 <- d[[i]]$IV3 
    }
   }
  
  eqlist <- c()
  for(i in 1:length(d)){
    eqlist <- c( eqlist,
                 paste(d[[i]]$DV, "_Int", sep=""),
                 paste(d[[i]]$DV, "_",d[[i]]$P, sep="") )
  }

  resmat <- NULL
  if (restrictions == TRUE){
    
    con <- miivs(model)$constr
    Rnames <- list(unlist(lapply(con, "[[", c("NAME"))), eqlist)
    R <- matrix(0, nrow = length(con), ncol = length(eqlist),
                dimnames = Rnames)
    L <- matrix(0, nrow = length(con), ncol = 1,
                dimnames = list(unlist(lapply(con, "[[", c("NAME"))),
                                c(" ")))
    
    for (r in 1:length(con)){
     if (con[[r]]$FIX == 0 ){ # if the coefficients isn't fixed to a value
       R[r, paste(con[[r]]$SET[1],"_",con[[r]]$DV[1], sep="")] <- 1
       R[r, paste(con[[r]]$SET[2],"_",con[[r]]$DV[2], sep="")] <- -1
       L[r] <- 0
     }
     
     if (con[[r]]$FIX != 0 ){
       R[r, paste(con[[r]]$SET,"_",con[[r]]$DV, sep="")] <- 1
       L[r] <- con[[r]]$FIX
     }
    }
    
    resmat <- list(R=R, L=L)
    
    for (i in 1:length(d)){
      if (i == 1){
        y  <- cbind(data[,d[[i]]$DV])
        Z  <- as.matrix(cbind(1, data[,d[[i]]$IV]))
        H  <- as.matrix(cbind(1, data[,d[[i]]$P]))
      }
      if (i >= 2){
        y <- rbind(y, cbind(data[,d[[i]]$DV]))
        Z <- as.matrix(bdiag(Z, as.matrix(cbind(1, data[,d[[i]]$IV]))))
        H <- as.matrix(bdiag(H, as.matrix(cbind(1, data[,d[[i]]$P]))))
      }
    }

    Xhat  <- Z %*% solve(crossprod(Z)) %*% t(Z) %*% H
    top   <- cbind(t(Xhat) %*% Xhat, t(R))
    bot   <- cbind(R, matrix(0, ncol=nrow(R), nrow=nrow(R)))
    theta_full <- as.matrix(solve(rbind(top, bot)) %*% (rbind(t(Xhat) %*% y, L)))
    theta_res <- cbind(theta_full[1:length(eqlist),])
    rownames(theta_res) <- eqlist
    
    for(i in 1:length(d)){
      temp <- c( paste(d[[i]]$DV, "_Int", sep=""),
                 paste(d[[i]]$DV, "_",d[[i]]$P, sep="") )
      d[[i]]$ESTR <- theta_res[which(rownames(theta_res) %in% temp)]
    }
      
  }

  
  dat    <- data.frame()
  for (i in 1:length(d)){
    y <- as.matrix(cbind(data[,d[[i]]$DV] ) )
    H <- as.matrix(cbind(1, data[,d[[i]]$P] ) )
    Z <- as.matrix(cbind(1, data[,d[[i]]$IV] ) )
    P  <- Z %*% solve(crossprod(Z)) %*% t(Z) # projection matrix
    b  <- solve(t(H) %*% P %*% H) %*% t(H) %*% P %*% y
    d[[i]]$EST <- as.numeric(b)
    if (restrictions == TRUE) {b <- as.numeric(d[[i]]$ESTR)}
    RS <- y - H  %*% cbind(b)
    rownames(RS) <- rep(d[[i]]$DV, nrow(RS))
    L0 <- as.numeric(crossprod(RS) / (nrow(data)))
    ZH <- Z %*% solve(crossprod(Z)) %*% crossprod(Z,H)
    AC <- L0 * (solve(crossprod(ZH)) %*% crossprod(ZH) %*% solve(crossprod(ZH)))
    d[[i]]$SE <- sqrt(diag(AC))
    IVID   <- solve(crossprod(Z)) %*% crossprod(Z,RS)
    SST    <- sum((RS-mean(RS))^2)
    SSE    <- sum((RS - (Z %*% IVID))^2)
    Rsq    <- 1 - (SSE/SST)
    sarg   <- length(RS)*Rsq
    length(sarg) <- length(d[[i]]$SE)
    sargdf <- ncol(Z)-ncol(H)
    length(sargdf) <- length(sarg)
    sargp  <- 1 - pchisq(sarg,df = sargdf) 
    length(sargp) <- length(sarg)
    temp   <- as.data.frame(cbind(b,d[[i]]$SE,sarg,sargp))
    temp$t <- temp[,1]/temp[,2]
    temp$p <- 2* (pnorm(abs(temp$t),lower.tail=FALSE))
    temp$p <- as.numeric(as.character(temp$p))
    temp$i <- cbind(c("Int", d[[i]]$P2))
    temp$d <- d[[i]]$DV
    
    RSU   <- y - H  %*% cbind(d[[i]]$EST)
    rownames(RSU) <- rep(d[[i]]$DV, nrow(RSU))
    L02 <- as.numeric(crossprod(RSU) / (nrow(data)))
    AC2  <- L02 * (solve(crossprod(ZH)) %*% crossprod(ZH) %*% solve(crossprod(ZH)))
    
      if (i == 1){ 
        if (restrictions == TRUE){ 
          SIGR  <- RS 
          br.cov <- AC
        }
        SIG   <- RSU
        b.cov <- AC2
      }
      if (i > 1 ){ 
        if (restrictions == TRUE){ 
          SIGR  <- rbind(SIGR, RS)
          br.cov <- as.matrix(bdiag(br.cov, AC))
          }
        SIG   <- rbind(SIG , RSU)
        b.cov <- as.matrix(bdiag(b.cov, AC2))
      }
    
      if (d[[i]]$NOTE != "") {
        sargdf <- paste(sargdf, "*", sep="")
        temp$df <-  sargdf
        temp$df[-1] <- NA
      }
    
      if (d[[i]]$NOTE == "") {
        temp$df <- sargdf
      }
    
      if (d[[i]]$EQ == "y1") { 
        temp$i <- cbind(c("Int", d[[i]]$P2))
        temp$d <- d[[i]]$DV2
      }

      if (i == 1) {dat <- temp }
      if (i >  1) {dat <- rbind(dat,temp) }
  }
    colnames(dat) <- c("Estimate", "StdErr","Sargan",
                      "P(Chi)", "Z", "P(|Z|)", "EV", "DV", "df")
    dat <- dat[, c("DV", "EV", "Estimate", "StdErr", 
                   "Z", "P(|Z|)","Sargan", "df", "P(Chi)")]
    
    
    if (restrictions == TRUE){ 
      rescov.r  <- matrix(0, nrow = length(d), ncol = length(d))
      for(r in 1:length(d)){
          for(c in 1:length(d)){
            rescov.r[r,c] <- (t(SIGR[which(rownames(SIGR) == d[[r]]$DV)]) %*%
                                SIGR[which(rownames(SIGR) == d[[c]]$DV)]) /
                                nrow(data)
          }
      }
      dimnames(rescov.r) <- list(unlist(lapply(d, "[[", "DV")), 
                                 unlist(lapply(d, "[[", "DV")))
    }

    rescov  <- matrix(0, nrow = length(d), ncol = length(d))
      for(r in 1:length(d)){
        for(c in 1:length(d)){
          rescov[r,c] <- (t(SIG[which(rownames(SIG) == d[[r]]$DV)]) %*%
                              SIG[which(rownames(SIG) == d[[c]]$DV)]) /
                              nrow(data)
        }
    }
    dimnames(rescov) <- list(unlist(lapply(d, "[[", "DV")), 
                             unlist(lapply(d, "[[", "DV")))

    for (i in 1:length(d)){
      LHS <- paste(d[[i]]$DV, collapse = ", ")
      RHS <- paste(d[[i]]$P, collapse = ", ")
      Instruments <- paste(d[[i]]$IV2, collapse = ", ")
      Disturbance <- paste("e.",d[[i]]$C, collapse = ", ", sep="")
      Effects <- paste(d[[i]]$EF, collapse = ", ") 
      modtemp <- as.data.frame(cbind(LHS, RHS, Disturbance, Instruments))
      colnames(modtemp) <- c("LHS", "RHS", "Composite Disturbance", "MIIVs")
      if (i == 1) {modeqns <- modtemp }
      if (i >  1) {modeqns <- rbind(modeqns,modtemp) }
    }
  
  restests   <- NULL
  if (restrictions == TRUE){ 
    restlabels <- unlist(lapply(con, "[[", c("NAME")))
    
    lrtest.est <- nrow(data) * (log(det(rescov.r)) - log(det(rescov)))
    lrtest.df <- nrow(R)
    lrtest.p <- pchisq(lrtest.est, lrtest.df, lower.tail = FALSE)
    lrtest.lab <- "Likelihood Ratio Test: Asymptotic Chi-squared"
    lrtest <- list(lrtest.est, lrtest.p, lrtest.df, lrtest.lab)
  
    e <- cbind(unlist(lapply(d, "[[","EST")))
    waldtest.est <- t(R %*% e - L) %*% solve(R %*% b.cov %*% t(R)) %*% (R %*% e - L)
    waldtest.df <- nrow(R)
    waldtest.p <- pchisq(waldtest.est, waldtest.df, lower.tail = FALSE)
    waldtest.lab <- "Wald Test: Asymptotic Chi-squared"
    waldtest <- list(waldtest.est, waldtest.p, waldtest.df, waldtest.lab)
    
    restests <- list(lrtest,waldtest, restlabels)
  }
  
  if (!is.null(varcov)){
    
    estimator <- varcov
    
    fit <- lavaan(model, auto.fix.first = TRUE, auto.var=TRUE, 
                  auto.cov.lv.x = TRUE, estimator = estimator)
    index <- which(fit@ParTable$plabel != "")
    fit@ParTable <- lapply(fit@ParTable, function (x) x[index])
    for (i in 1:length(fit@ParTable$lhs)){
      for (j in 1:length(d)){
        for (k in 1:length(d[[j]]$P)){
          if (d[[j]]$EQ == "y1"){
            lavrhs <- d[[j]]$P2[k]
            lavlhs <- d[[j]]$DV2
              if (fit@ParTable$op[i]  == "~"    &&
                  fit@ParTable$lhs[i] == lavlhs &&
                  fit@ParTable$rhs[i] == lavrhs)
                  {
                    fit@ParTable$label[i] <-  ""
                    fit@ParTable$free[i] <-  0
                    fit@ParTable$ustart[i] <-  as.numeric(d[[j]]$EST[k+1])
                  }
          }
      
          if (d[[j]]$EQ != "y1"){
            lavlhs <- d[[j]]$P2[k]
            lavrhs <- d[[j]]$DV
              if (fit@ParTable$op[i]  == "=~"   && 
                  fit@ParTable$lhs[i] == lavlhs && 
                  fit@ParTable$rhs[i] == lavrhs)
                  {
                    fit@ParTable$label[i] <-  ""
                    fit@ParTable$free[i] <-  0
                    fit@ParTable$ustart[i] <-  as.numeric(d[[j]]$EST[k+1])
                  }
          }
        }
      }
    }
    
  lavsyntax <- lavExport(fit, export = FALSE)
  fit2 <- lavaan(lavsyntax, data = data, 
                 auto.fix.first = TRUE, 
                 auto.var=TRUE, auto.cov.lv.x = TRUE)
   pe  <-  parameterEstimates(fit2)
   cpe <-  pe[pe$op == "~~" & pe$lhs != pe$rhs,]
   cpe$x <- paste(cpe$lhs, cpe$op, cpe$rhs, sep=" ")
   cpe <- cpe[,c("x","est", "se", "z", "pvalue")]
   colnames(cpe) <- c( "","Estimate", 
                        "StdErr", "Z", "P(|Z|)")
   
   vpe <-  pe[pe$op == "~~" & pe$lhs == pe$rhs,]
   vpe <- vpe[,c("lhs","est", "se", "z", "pvalue")]
   colnames(vpe) <- c( "","Estimate", 
                       "StdErr", "Z", "P(|Z|)")
  vcov <- list(cov = cpe, var = vpe)
  }
  
  lavfit <- NULL
  if (!is.null(varcov)){lavfit <- fit2}

  ctrlopts <- list(printmiivs = printmiivs, restrictions = restrictions, 
                   varcov = varcov)
  
  res <- list(model = d, dat = dat, modeqns = modeqns, ctrlopts = ctrlopts, 
              restests = restests, vcov = vcov, lavfit = lavfit, resmat = resmat)
  class(res) <- "miive"
  res
}
  
 