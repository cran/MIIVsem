#' Model-implied instrumental variable (MIIV) search 
#'
#' A key step in the MIIV-2SLS approach is to transform the SEM by replacing the latent variables with their scaling indicators minus their errors.  Upon substitution the SEM is transformed from a model with latent variables to one in observed variables with composite errors.  The miivs function automatically makes this transformation.
#'
#' @param model A model specified using lavaan model syntax. See the \code{model} argument within the \code{\link[lavaan]{lavaanify}} function for more information.
#'
#' @return eqns
#' @return modeqns
#' 
#' @details 
#' \itemize{
#'  \item \code{LHS} The "dependent" variable.
#'  \item \code{RHS} The right hand side variables of the transformed equation.
#'  \item \code{Composite Disturbance}  Elements of the composite errors in the transformed equation.
#'  \item \code{MIIVs} The model implied instrumental variables for each equation.
#' }
#' 
#' @references 
#' 
#'  Bollen,	K. A. and	D. J.	Bauer.	2004.	Automating	the	Selection	of 
#' 	Model-Implied	Instrumental	Variables.	\emph{Sociological	Methods	and	
#' 	Research}, 32, 425-52.
#' 	
#'
#' @examples
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
#'  miivs(bollen1989a_model)
#'  
#' @export
miivs <- function(model) {
    fit <- lavaan(model, auto.fix.first = TRUE, auto.var=TRUE)
    pt  <- parTable(fit)
    df  <- lavMatrixRepresentation(pt)
    LA  <- inspect(fit)$lambda # lambda matrix for y2, x2
    BA  <- inspect(fit)$beta
    TH  <- inspect(fit)$theta
    PS  <- inspect(fit)$psi
    
    # Get rid of zero'd beta values should they exist
    df  <- df[!(as.character(df$mat) == "beta" & (df$free == 0)), ]
    
    # Get original variable names
    latEnd <- unique(df[grep("beta", df$mat), ]$lhs)
    latExo <- unique(setdiff(df[grep("beta", df$mat), ]$rhs, latEnd))
    #manExo <- subset(df, lhs %in% latExo & mat %in% "lambda")$rhs #fix
    manExo <- df[df$lhs %in% latExo & df$mat %in% "lambda", ]$rhs
    #manEnd <- subset(df, lhs %in% latEnd & mat %in% "lambda")$rhs #fix
    manEnd <- df[df$lhs %in% latEnd & df$mat %in% "lambda", ]$rhs
    latAll <- c(latEnd, latExo)
    manAll <- c(latEnd, latExo)
    
    # Get counts for each variable type
    numLatEnd <- length(latEnd)
    numLatExo <- length(latExo)
    numManExo <- length(manExo)
    numManEnd <- length(manEnd)
    numLatAll <- length(latAll)
    numManAll <- length(manAll)
    
    # Get scaling and non-scaling indicators in x and y notation, dvs
    #y1 <- subset(df, lhs %in% latEnd & ustart == 1)$rhs #fix
    y1 <- df[df$lhs %in% latEnd & df$ustart == 1, "rhs"]
    y1 <- y1[!is.na(y1)]
    y2 <- setdiff(manEnd, y1)
    #x1 <- subset(df, lhs %in% latExo & ustart == 1)$rhs #fix
    x1 <- df[df$lhs %in% latExo & df$ustart == 1, "rhs"]
    x1 <- x1[!is.na(x1)]
    x2 <- setdiff(manExo, x1)
    xy <- c(y1, y2, x1, x2)
    dv <- c(y1, y2, x2)
    
    # Get lambda matrices for y2 and x2 only
    LA <- LA[c(y2, x2), ] 
    
    # Names functions
    #names <- subset(df[, c("rhs", "lhs")], df$ustart == 1) # e.g (Xi1, x1)
    names <- na.omit(df[df$ustart == 1, c("rhs", "lhs")])
    names <- rbind(names, data.frame("rhs" = xy, "lhs" = xy))
    
    # Replaces 
    repStr <- function(x) {
        as.character(names$rhs[match(x, names$lhs)])
    }
    
    # Replaces 
    repStrInv <- function(x) {
        as.character(names$lhs[match(x, names$rhs)])
    }
    
    repElemWithColName <- function(mat) {
      class(mat) <- "character"
        for (i in 1:nrow(mat)) {
            for (j in 1:ncol(mat)) {
                if (mat[i, j] != 0) {
                  mat[i, j] <- colnames(mat)[j]
                }
            }
        }
      t <- split(mat, row(mat))
      names(t) <- row.names(mat)
      t <- sapply(t,function(x) x[x != 0], simplify = FALSE, USE.NAMES = TRUE)
      return(t)
    }
    
    TH_temp <- repElemWithColName(TH)
    PS_t <- repElemWithColName(PS[latEnd, latEnd])
    names(PS_t) <- sapply(names(PS_t), function(x) repStr(x), simplify=FALSE)
    PS_temp <- sapply(PS_t, function(x) repStr(x), simplify = FALSE)
    BA_r <- BA[latEnd, latEnd]
    BA_r <- solve(diag(nrow(BA_r)) - BA_r)
    BA_x <- repElemWithColName(BA_r)
    BA_x <- BA_x[lapply(BA_x, length) > 0]
    LA_r <- LA[, latEnd]
    LA_r <- LA_r %*% BA_r
    LA_x <- repElemWithColName(LA_r)
    LA_x <- LA_x[lapply(LA_x, length) > 0]
    EF_x <- append(BA_x, LA_x)
    names(EF_x) <- sapply(names(EF_x), function(x) repStr(x))
    
    # Create error index
    #err_df <- subset(df[, c("lhs", "rhs")], df$mat == "theta")
    err_df <- na.omit(df[df$mat == "theta", c("lhs", "rhs")])
    err_df$lhs <- sapply(err_df$lhs, function(x) repStr(x))
    #err_df_add <- subset(df[, c("lhs", "rhs")], df$mat == "psi") # get psi rows
    err_df_add <- na.omit(df[df$mat == "psi", c("lhs", "rhs")]) # get psi rows
    
    # Get variances from psi matrix
    err_df_var <- err_df_add[err_df_add[, "lhs"] == err_df_add[, "rhs"], ]
    
    # Composite disturbance doesn't include covariances 
    err_df_var$rhs <- sapply(err_df_var$rhs, function(x) repStr(x))
    err_df <- rbind(err_df, err_df_var)

    # Get df of scaling var names and equivalent latent var names
    #parent_df <- subset(df[, c("lhs", "rhs")], df$mat == "lambda")
    #parent_df_add <- subset(df[, c("rhs", "lhs")], df$mat == "beta")
    parent_df <- na.omit(df[df$mat == "lambda", c("lhs", "rhs")])
    parent_df_add <- na.omit(df[df$mat == "beta", c("rhs", "lhs")])
    parent_df <- rbind(parent_df, setNames(parent_df_add , names(parent_df)))
    parent_df <- sapply(parent_df, function(x) repStr(x))
    
    # Get df of scaling var names and equivalent latent var names
    #lName_df <- subset(df[, c("rhs", "label")], df$label != "")
    lName_df <- na.omit(df[df$label != "", c("rhs", "label")])
    
    # Set up Master List for Equations
    eqns <- list(DV = "", DV2="", IV = "", L = "", EQ = "", PA = "", C = "", 
                 EF = "", PIV = "", W = "", P="", IV2="", P2="", IV3 = "",
                 MSG = "", NOTE = "")
    
    eqns <- replicate(length(dv), eqns, simplify = FALSE)

    for (i in 1:length(eqns)) {
        eqns[[i]]$DV <- dv[i]
        if (eqns[[i]]$DV %in% y1) {
            eqns[[i]]$EQ <- "y1"
        }
        if (eqns[[i]]$DV %in% y2) {
            eqns[[i]]$EQ <- "y2"
        }
        if (eqns[[i]]$DV %in% x2) {
            eqns[[i]]$EQ <- "x2"
        }
    }
    
    # Add labels to eqns
    if (nrow(lName_df) > 0){
    for (i in 1:length(eqns)) {
      for (j in 1:nrow(lName_df)) {
       if (eqns[[i]]$DV == lName_df[j,1]){
          eqns[[i]]$L <- lName_df[j,2]
       }
      }
    }
    }
    
    # Add parents to eqns
    for (i in 1:length(eqns)) {
        for (j in 1:nrow(parent_df)) {
            if (eqns[[i]]$DV == parent_df[j, 2]) {
                eqns[[i]]$PA <- c(eqns[[i]]$PA, as.character(parent_df[j, 1]))
                eqns[[i]]$PA <- eqns[[i]]$PA[eqns[[i]]$PA != ""]
                eqns[[i]]$P  <- eqns[[i]]$PA[eqns[[i]]$PA != eqns[[i]]$DV]
            }
        }
    }
    
    # Add latent var labels
    for (i in 1:length(eqns)) {
       if (eqns[[i]]$EQ == "y1"){
        eqns[[i]]$DV2 <- sapply(eqns[[i]]$DV, function(x) repStrInv(x))
        eqns[[i]]$P2 <- unname(sapply(eqns[[i]]$P, function(x) repStrInv(x)))
       }
       
       if (eqns[[i]]$EQ == "y2" | eqns[[i]]$EQ == "x2"){ #added
        eqns[[i]]$DV2 <- sapply(eqns[[i]]$P, function(x) repStrInv(x)) #added
        eqns[[i]]$P2 <- unname(sapply(eqns[[i]]$P, function(x) repStrInv(x)))
       }#added
    }
    
    # Add direct composite erros to equations
    for (i in 1:length(eqns)) {
        eqns[[i]]$C <- eqns[[i]]$PA  # add parent error
        eqns[[i]]$C <- c(eqns[[i]]$C, eqns[[i]]$DV)  # add direct error
        if (eqns[[i]]$EQ == "y1") {
            for (j in 1:nrow(err_df)) {
                if (eqns[[i]]$DV == err_df[j, 2]) {
                  eqns[[i]]$C <- c(eqns[[i]]$C, err_df[j, 1])
                  eqns[[i]]$C <- eqns[[i]]$C[eqns[[i]]$C != ""]
                  eqns[[i]]$C <- unique(eqns[[i]]$C)
                }
            }
        }
    }
    
    # Add Effects to Matrix
    for (i in 1:length(eqns)) {
        for (j in 1:length(EF_x)) {
            for (k in 1:length(eqns[[i]]$PA)) {
                if (eqns[[i]]$PA[k] == names(EF_x)[j]) {
                  eqns[[i]]$EF <- c(eqns[[i]]$EF, EF_x[j])
                }
            }
        }
        eqns[[i]]$EF <- c(eqns[[i]]$EF, eqns[[i]]$DV)
        eqns[[i]]$EF <- as.character(unlist(eqns[[i]]$EF))
        eqns[[i]]$EF <- eqns[[i]]$EF[eqns[[i]]$EF != ""]
        eqns[[i]]$EF <- unique(eqns[[i]]$EF)
    }
    
    # Add Effects
    effects <- list(DV = "", EF = "")
    effects <- replicate(numLatExo, effects, simplify = FALSE)
    for (i in 1:(length(effects))) {
        effects[[i]]$DV <- x1[i]
        effects[[i]]$EF <- x1[i]
    }
    effects <- append(lapply(eqns, "[", c("DV", "EF")), effects)
    
    # Potential Instruments
    for (i in 1:length(eqns)) {
        C <- unlist(eqns[[i]]$C)
        for (j in 1:length(effects)) {
            E <- unlist(effects[[j]]$EF)
            if (!any(C %in% E)) {
                eqns[[i]]$PIV <- c(eqns[[i]]$PIV, effects[[j]]$DV)
            }
        }
        eqns[[i]]$PIV <- eqns[[i]]$PIV[eqns[[i]]$PIV != ""]
    }
    
    # Final Instruments
    for (i in 1:length(eqns)) {
        C <- eqns[[i]]$C
        for (j in 1:length(C)) {
            W <- C[j]
            for (k in 1:length(TH_temp)) {
                if (W == names(TH_temp)[k]) {
                  W_temp <- as.character(unlist(TH_temp[k]))
                  eqns[[i]]$W <- c(eqns[[i]]$W, W_temp)
                  eqns[[i]]$W <- eqns[[i]]$W[eqns[[i]]$W != ""]
                }
            }
        }
        eqns[[i]]$IV <- setdiff(eqns[[i]]$PIV, eqns[[i]]$W)
    }
    
    # Eliminate Correlated zetas 
    for (i in 1:length(eqns)) {
      if (eqns[[i]]$EQ == "y1" ){
        D <- eqns[[i]]$DV
            for (k in 1:length(PS_temp)) {
                if (D == names(PS_temp)[k]) {
                  W_temp <- as.character(unlist(PS_temp[k]))
                  eqns[[i]]$W <- c(eqns[[i]]$W, W_temp)
                  eqns[[i]]$W <- eqns[[i]]$W[eqns[[i]]$W != ""]
                }
            }
        }
        eqns[[i]]$IV <- setdiff(eqns[[i]]$PIV, eqns[[i]]$W)
    }
    
    for (i in 1:length(eqns)){
      LHS <- paste(eqns[[i]]$DV, collapse = ", ")
      RHS <- paste(eqns[[i]]$P, collapse = ", ")
      Instruments <- paste(eqns[[i]]$IV, collapse = ", ")
      Disturbance <- paste("e.",eqns[[i]]$C, collapse = ", ", sep="")
      modtemp <- as.data.frame(cbind(LHS, RHS, Disturbance, Instruments))
      colnames(modtemp) <- c("LHS", "RHS", "Composite Disturbance", "MIIVs")
      if (i == 1) {modeqns <- modtemp }
      if (i >  1) {modeqns <- rbind(modeqns,modtemp) }
    } 
  
  search <- list(eqns = eqns, df = modeqns )
  class(search) <- "miivs"
  search
}
