#' @keywords internal
miivboot <- function(d, dat, restrictions, data, bootstrap.se, reps, R, L, B) {
  if (restrictions == TRUE & bootstrap.se == "residual"){
      
    stop(paste("MIIVsem does not currently support the residuals bootstrapping 
                  procedure for restricted models."))
      
  }
  
  if (restrictions == TRUE & bootstrap.se == "pairs"){  
    
    stop(paste("MIIVsem does not currently support the pairs bootstrapping 
                procedure for restricted models."))
    
#     b0   <- unlist(lapply(d, "[[", "EST"))
#     b0   <- matrix(b0, nrow = length(b0))
#     dvs0 <- unlist(lapply(d, "[[", "DVobs"))
#     
#     for (j in 1:reps){
#     
#     bdata <- data[sample(nrow(data),size=nrow(data),replace=TRUE),]
#     
#     for (i in 1:length(d)){
#       if (i == 1){
#         y1     <- cbind(bdata[,d[[i]]$DVobs])
#         X1     <- as.matrix(cbind(1,bdata[,d[[i]]$IV]))
#         Y1     <- as.matrix(cbind(1,bdata[,d[[i]]$IVobs]))
#       }
#       if (i >= 2){
#         y1    <- rbind(y1, cbind(bdata[,d[[i]]$DVobs]))
#         X1    <- as.matrix(bdiag(X1, as.matrix(cbind(1,bdata[,d[[i]]$IV]))))
#         Y1    <- as.matrix(bdiag(Y1, as.matrix(cbind(1,bdata[,d[[i]]$IVobs]))))
#       }
#       d[[i]]$IVobsInt <- c(paste(d[[i]]$DVobs,"_Int", sep= ""), d[[i]]$IVobs)
#     }
#     
#     y1 = as(Matrix(y1), "dgeMatrix")
#     X1 = as(Matrix(X1), "sparseMatrix")
#     Y1 = as(Matrix(Y1), "sparseMatrix")
#     
#     X1inv    <- solve(crossprod(X1)) 
#     Y1hat    <- X1%*% X1inv %*% t(X1)%*%Y1
#     b <- solve(rbind(cbind(t(Y1hat) %*% Y1hat, t(R)), cbind(R, matrix(0, 
#            ncol=nrow(R), nrow=nrow(R))))) %*% rbind(t(Y1hat) %*% y1, L)
#     b <- cbind(b[1:ncol(R)])
#     #b.unr <- solve(t(Y1hat)%*%Y1) %*% t(Y1hat)%*%y1
#     #b.unr <- cbind(b.unr[1:ncol(R)])
#     
#     Y         <- as(data[,all, drop=FALSE], "Matrix")
#     Y.c       <- as(apply(Y, 2, function(y) y - mean(y)), "Matrix")
#     E         <- Y.c %*% B
#     S.EE      <- (t(B) %*% crossprod(Y.c) %*% B) / N
#     I <- Diagonal(N)
#     omega <- solve(kronecker(diag(diag(S.EE[unlist(dvs0),unlist(dvs0)])), I) )
#     top   <- cbind(t(Y1hat) %*% omega %*% Y1hat, t(R))
#     bot   <- cbind(R, matrix(0, ncol=nrow(R), nrow=nrow(R)))
#     S.EEr <- as.matrix(solve(rbind(top, bot)) )
#     S.EEr <- (S.EEr[1:nrow(b),1:nrow(b)])
#     se <- cbind(sqrt(diag(S.EEr)))
#     t <- (b0 - b)/se
#     cbind(b,t)
    
  }
  
  if(restrictions == FALSE & bootstrap.se == "residual"){
  
  for (i in 1:length(d)){
    y  <- as.matrix(cbind(data[,d[[i]]$DVobs] ) )
    H  <- as.matrix(cbind(1, data[,d[[i]]$IVobs] ) )
    Z  <- as.matrix(cbind(1, data[,d[[i]]$IV] ) )
    P  <- Z %*% solve(crossprod(Z)) %*% t(Z) 
    b  <- solve(t(H) %*% P %*% H) %*% t(H) %*% P %*% y
    u1 <- y - H  %*% cbind(b)
    for ( k in 2:ncol(H)){
      if ( k ==2 ) { 
        u2 <- cbind(resid(lm(H[,k] ~ Z[,-1]))) 
        b2 <- cbind(coef(lm(H[,k] ~ Z[,-1])))
        }
      if ( k  >2 ) { 
        u2 <- cbind(u2, resid(lm(H[,k] ~ Z[,-1]))) 
        b2 <- cbind(b2, coef(lm(H[,k] ~ Z[,-1])))
        }
    }
    
    l  <- length(d[[i]]$IV)
    n  <- nrow(data)
    sc <- sqrt(n/(n-1))
    u2 <- u2 * sc
    uu <- cbind(u1,u2)
    
    for (j in 1:reps){
      #j2 <- j + 1
      u  <- uu[sample(nrow(uu),size=nrow(uu),replace=TRUE),]
      Z2 <- Z[sample(nrow(Z),size=nrow(Z),replace=TRUE),]
      for ( k in 2:ncol(H)){
        if ( k ==2 ) { 
          H2 <- Z2  %*% cbind(b2[,k-1]) + u[,k] 
        }
        if ( k  >2 ) { 
          H2 <- cbind(H2, Z2%*%cbind(b2[,k-1]) + u[,k]) 
        }
      }
      H2 <- cbind(1,H2)
      y2 <- H2  %*% cbind(b) + u[,1] 
      #y1 <- H  %*% cbind(b) + u[,1] # changed this
      # add this below
      ZZinv <- solve(crossprod(Z2))
      b3 <- solve(t(H2)%*%Z2%*%ZZinv%*%t(Z2)%*%H2)%*%t(H2)%*%Z2%*%ZZinv%*%t(Z2)%*%y2
      
      # add this above
      #y1 <- y2  %*% cbind(b) + u[,1] # changed this
      #y <- as.matrix(cbind(y1))
      #H <- as.matrix(cbind(1, y2))
      #P <- Z %*% solve(crossprod(Z)) %*% t(Z) 
      #b3 <- solve(t(H) %*% P %*% H) %*% t(H) %*% P %*% y
      RS <- y2 - H2  %*% cbind(b3)
      L0 <- as.numeric(crossprod(RS) / (nrow(data)))
      ZH <- Z2 %*% solve(crossprod(Z2)) %*% crossprod(Z2,H2)
      AC <- L0 * (solve(crossprod(ZH)) %*% crossprod(ZH) %*% solve(crossprod(ZH)))
      se <- sqrt(diag(AC))
      t <- (b - b3)/se
      dfz <- cbind(t(b3), t(t))
      varnames <- c(paste("Int_", d[[i]]$DVobs,sep=""), paste(d[[i]]$DVobs,"_",d[[i]]$IVobs,sep=""))
      tnames <- paste("test_", varnames,sep="")
      colnames(dfz) <- c(varnames,tnames)
      if (j == 1){bs <- dfz}
      if (j  > 1){bs <- rbind(bs,dfz)}
    }
    if (i == 1){bd <- bs}
    if (i  > 1){bd <- cbind(bd,bs)}
  }
  
  ttests <- bd[,grep("test_", colnames(bd))]
  ttests.sorted <- apply(ttests,2,sort,decreasing=F)
  crit.vals <- ttests.sorted[995,]

  betas <- bd[,-grep("test_", colnames(bd))]
  means <- as.matrix(colMeans(betas), nrow=1,ncol=ncol(bd))
  dimnames(means)[2] <- "btsrp.mean"
  means <- t(means)

  for (i in 1:ncol(betas)){
    if(i == 1) {ses <- sqrt(sum((betas[,i] - mean(betas[,i]))^2)/nrow(betas))}
    if(i  > 1) {ses <- c(ses, sqrt(sum((betas[,i] - mean(betas[,i]))^2)/nrow(betas)))}
  }

  ses <- as.matrix(ses, nrow=1, ncol=length(ses))
  colnames(ses)[1] <- list("btsrp.se")
  rownames(ses) <- colnames(means)
  ses <- t(ses)

  for (i in 1:ncol(ttests.sorted)){
    #p.one.tail <- length(which(ttests.sorted[,i] > crit.vals[i]))/nrow(ttests.sorted)
    #p.sym <- length(which(abs(ttests.sorted[,i]) > abs(crit.vals[i])))/nrow(ttests.sorted)
    p.non.sym <- 2* min(c(
      length(which(abs(ttests.sorted[,i]) < abs(crit.vals[i])))/nrow(ttests.sorted),
      length(which(abs(ttests.sorted[,i]) > abs(crit.vals[i])))/nrow(ttests.sorted)))
    #p <- rbind(p.one.tail, p.sym, p.non.sym)
    p <- p.non.sym
    if (i == 1){df.p <- p}
    if (i  > 1){df.p <- cbind(df.p,p)}
  }
  resid.df <- t(rbind(means,ses,df.p))
  colnames(resid.df) <- c("b", "se", "p")
  
  dat$StdErr   <- resid.df[,"se"]
  dat$Z <- NULL
  dat$`P(|Z|)` <- NULL
  dat$p <- resid.df[,"p"]
  dat <- dat[,c("DV", "EV", "Estimate", "StdErr", "p", "Sargan", "df", "P(Chi)")]
  colnames(dat) <- c("DV", "EV", "Estimate", "BtSD", "P(t)", "Sargan", "df", "P(Chi)")
  }
  
  if(restrictions == FALSE & bootstrap.se == "pairs"){

  for (j in 1:reps){
    bdata <- data[sample(nrow(data),size=nrow(data),replace=TRUE),]
    for (i in 1:length(d)){
      y <- as.matrix(cbind(bdata[,d[[i]]$DVobs] ) )
      H <- as.matrix(cbind(1, bdata[,d[[i]]$IVobs] ) )
      Z <- as.matrix(cbind(1, bdata[,d[[i]]$IV] ) )
      P <- Z %*% solve(crossprod(Z)) %*% t(Z) 
      b <- solve(t(H) %*% P %*% H) %*% t(H) %*% P %*% y
      RS <- y - H  %*% cbind(b)
      L0 <- as.numeric(crossprod(RS) / (nrow(bdata)))
      ZH <- Z %*% solve(crossprod(Z)) %*% crossprod(Z,H)
      AC <- L0 * (solve(crossprod(ZH)) %*% crossprod(ZH) %*% solve(crossprod(ZH)))
      se <- sqrt(diag(AC))
      b0 <- d[[i]]$EST
      t <- (b0 - b)/se
      dfz <- cbind(t(b), t(t))
      varnames <- c(paste("Int_", d[[i]]$DVobs,sep=""), paste(d[[i]]$DVobs,"_",d[[i]]$IVobs,sep=""))
      tnames <- paste("test_", varnames,sep="")
      colnames(dfz) <- c(varnames,tnames)
      if (i == 1){bs <- dfz}
      if (i  > 1){bs <- cbind(bs,dfz)}
    }
    if (j == 1){bd <- bs}
    if (j  > 1){bd <- rbind(bd,bs)}
  }
  
  ttests <- bd[,grep("test_", colnames(bd))]
  ttests.sorted <- apply(ttests,2,sort,decreasing=F)
  crit.vals <- ttests.sorted[995,]

  betas <- bd[,-grep("test_", colnames(bd))]
  means <- as.matrix(colMeans(betas), nrow=1,ncol=ncol(bd))
  dimnames(means)[2] <- "btsrp.mean"
  means <- t(means)
  
  for (i in 1:ncol(betas)){
    if(i==1) {ses <- sqrt(sum((betas[,i] - mean(betas[,i]))^2)/nrow(betas))}
    if(i >1) {ses <- c(ses, sqrt(sum((betas[,i] - mean(betas[,i]))^2)/nrow(betas)))}
  }

  ses <- as.matrix(ses, nrow=1, ncol=length(ses))
  colnames(ses)[1] <- list("btsrp.se")
  rownames(ses) <- colnames(means)
  ses <- t(ses)

  for (i in 1:ncol(ttests.sorted)){
    #p.one.tail <- length(which(ttests.sorted[,i] > crit.vals[i]))/nrow(ttests.sorted)
    #p.sym <- length(which(abs(ttests.sorted[,i]) > abs(crit.vals[i])))/nrow(ttests.sorted)
    p.non.sym <- 2* min(c(
      length(which(abs(ttests.sorted[,i]) < abs(crit.vals[i])))/nrow(ttests.sorted),
      length(which(abs(ttests.sorted[,i]) > abs(crit.vals[i])))/nrow(ttests.sorted)))
    #p <- rbind(p.one.tail, p.sym, p.non.sym)
    p <- rbind(p.non.sym)
    if (i ==1){df.p <- p}
    if (i  >1){df.p <- cbind(df.p,p)}
  }
  
  pairs.df <- t(rbind(means,ses,df.p))
  colnames(pairs.df) <- c("b", "se", "p")
  
  dat$StdErr   <- pairs.df[,"se"]
  dat$Z <- NULL
  dat$`P(|Z|)` <- NULL
  dat$p <- pairs.df[,"p"]
  dat <- dat[,c("DV", "EV", "Estimate", "StdErr", "p", "Sargan", "df", "P(Chi)")]
  colnames(dat) <- c("DV", "EV", "Estimate", "BtSD", "P(t)", "Sargan", "df", "P(Chi)")
  }

  return(dat)
}
