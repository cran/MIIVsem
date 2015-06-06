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
#' 	Bauldry, S.	2014.	miivfind: A command for identifying model-implied instrumental 
#' 	variables for structural equation models in Stata.	\emph{Stata Journal}, 14:4.
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
  
    fit <- lavaan(model, auto.fix.first = TRUE, auto.var=TRUE, auto.cov.lv.x = TRUE)
    pt  <- parTable(fit)
    df  <- lavMatrixRepresentation(pt)
    LA  <- inspect(fit)$lambda # lambda matrix for y2, x2
    BA  <- inspect(fit)$beta
    TH  <- inspect(fit)$theta
    PS  <- inspect(fit)$psi
    
    # Get rid of zero'd beta values should they exist
    #df  <- df[!(as.character(df$mat) == "beta" & (df$free == 0)), ]
    
    # for manifest vars which are not indicators
    # indicators   <- unique(df[df$op == "=~", ]$rhs) 
    #origlatentvars <- unique(df[df$op == "=~", ]$lhs) 
    # manifestXvars <- setdiff(unique(df[df$op == "~", ]$rhs))
    # manifestYvars <- unique(df[df$op == "~", ]$lhs))
    #setdiff(latEnd <- unique(df[df$op == "=~", ]$lhs) )
    
    # Get original variable names
    if (is.element('beta', df$mat)){
      latEnd <- unique(df[df$op == "=~", ]$lhs) 
      latEnd <- intersect(unique(df[grep("beta", df$mat), ]$lhs), latEnd)
      latExo <- unique(df[df$op == "=~", ]$lhs) 
      latExo <- intersect(df[grep("beta", df$mat), ]$rhs, latExo)
      latExo <- setdiff(latExo, latEnd)
    }
    if (!is.element('beta', df$mat)){
      latEnd <- c()
      latExo <- unique(df[grep("lambda", df$mat), ]$lhs)
    }
    manExo <- df[df$lhs %in% latExo & df$mat %in% "lambda", ]$rhs
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
    # Take into account additional vars fixed to 1
    # For now the first indicator for each latvar is the scaling var.
    y1 <- df[df$lhs %in% latEnd & df$ustart == 1 & !duplicated(df$lhs), "rhs"]
    y1 <- y1[!is.na(y1)]
    y2 <- setdiff(manEnd, y1)
    x1 <- df[df$lhs %in% latExo & df$ustart == 1 & !duplicated(df$lhs), "rhs"]
    x1 <- x1[!is.na(x1)]
    x2 <- setdiff(manExo, x1)
    xy <- c(y1, y2, x1, x2)
    dv <- c(y1, y2, x2)
        
    # Get lambda matrices for y2 and x2 only
    LA <- LA[c(y2, x2), ] 
    
    # New names function    
    trans <- na.omit(df[df$ustart == 1, c("rhs", "lhs")])
    trans <- rbind(trans, data.frame("rhs" = c(y2,x2), "lhs" = c(y2,x2)))
    colnames(trans) <- c("obs", "lat")
    
    # Names functions
    names <- na.omit(df[df$ustart == 1, c("lhs", "rhs")])
    names <- rbind(names, data.frame("lhs" = xy, "rhs" = xy))
    
    # Replace latent with obs
    Lat2Obs <- function(x) {
        as.character(trans$obs[match(x, trans$lat)])
    }
    # Replace obs with lat
    Obs2Lat <- function(x) {
        as.character(trans$lat[match(x, trans$obs)])
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
    
    if (numLatEnd > 0){
      PS_t <- repElemWithColName(PS[latEnd, latEnd,drop = FALSE])
      names(PS_t) <- sapply(names(PS_t), function(x) Lat2Obs(x), simplify=FALSE)
      PS_temp <- sapply(PS_t, function(x) Lat2Obs(x), simplify = FALSE)
      BA_r <- BA[latEnd, latEnd, drop = FALSE]
      BA_r <- solve(diag(nrow(BA_r)) - BA_r)
      BA_x <- repElemWithColName(BA_r)
      BA_x <- BA_x[lapply(BA_x, length) > 0]
      LA_r <- LA[, latEnd,  drop = FALSE]
      LA_r <- LA_r %*% BA_r
      LA_x <- repElemWithColName(LA_r) 
      LA_x <- LA_x[lapply(LA_x, length) > 0]
      EF_x <- append(BA_x, LA_x)
      names(EF_x) <- sapply(names(EF_x), function(x) Lat2Obs(x))
    }
    
    # Create error index
    err_df <- na.omit(df[df$mat == "theta", c("lhs", "rhs")])
    err_df_add <- na.omit(df[df$mat == "psi", c("lhs", "rhs")]) # get psi rows
    
    # Get variances from psi matrix
    err_df_var <- err_df_add[err_df_add[, "lhs"] == err_df_add[, "rhs"], ]
    
    # Composite disturbance doesn't include covariances 
    err_df_var$rhs <- sapply(err_df_var$rhs, function(x) Lat2Obs(x))
    err_df <- rbind(err_df, err_df_var)

    # Get df of scaling var names and equivalent latent var names
    obsdf1 <- df[df$mat %in% c("beta"), 
                   c("lhs", "rhs", "label", "ustart")]
    obsdf2 <- df[df$mat %in% c("lambda") & df$rhs %in% c(y2, x2), 
                   c("rhs", "lhs", "label", "ustart")]
    obsdf2 <- setNames(obsdf2 , names(obsdf1))
    obsdf  <- rbind(obsdf1, obsdf2)
    obsdf[,1:2] <- sapply(obsdf[,1:2], function(x) Lat2Obs(x))
    obsdf$label[which(obsdf$label=="")] <- "NA"
    obsdf$ustart[which(obsdf$ustart=="")] <- "NA"
    
    # Get fixed coefficients
    rs <- rbind(df[!is.na(df$ustart) & df$rhs %in% c(y2, x2), ],
                df[df$op == "==", ])
     
    # Create contraints list
    constr <- list(DV = "", SET = "", FIX = "", NAME = "")
    constr <- replicate(nrow(rs), constr, simplify = FALSE)
    
    if (nrow(rs) > 0){
      for (i in 1:nrow(rs)){
        if (rs[i, "op"] == "=="){
          constr[[i]]$DV   <- c(Lat2Obs(df[df$plabel == rs$lhs[i], "lhs"]),
                                Lat2Obs(df[df$plabel == rs$rhs[i], "lhs"]))
          constr[[i]]$SET  <- c(df[df$plabel == rs$lhs[i], "rhs"],
                                df[df$plabel == rs$rhs[i], "rhs"])
          constr[[i]]$FIX  <- 0
          constr[[i]]$NAME <- paste(constr[[i]]$SET[1],
                                    " = ",
                                    constr[[i]]$SET[2],
                                    sep = "")
        }
        if (rs[i, "op"] != "=="){
         constr[[i]]$DV   <- Lat2Obs(rs[i, "lhs"])
         constr[[i]]$SET  <- rs[i, "rhs"]
         constr[[i]]$FIX  <- as.numeric(rs[i, "ustart"])
         constr[[i]]$NAME <- paste(constr[[i]]$SET,
                                    " = ",
                                    constr[[i]]$FIX,
                                    sep = "")
        }
      }
    }
    # Set up Master List for Equation
    eqns <- list(DV = "", DV2="", IV = "", L = "", EQ = "", PA = "", C = "", 
                 EF = "", PIV = "", W = "", P="", IV2="", P2="", EQNUM = "",
                 FIX = "", EST = "", ESTR = "", SE = "", MSG = "", NOTE = "")
    
    eqns <- replicate(length(dv), eqns, simplify = FALSE)

    for (i in 1:length(eqns)) {
        eqns[[i]]$DV    <- dv[i]
        eqns[[i]]$EQNUM <- i
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

    # Add predictos and equality constraints to eqns
    for (i in 1:length(eqns)) {
        for (j in 1:nrow(obsdf)) {
            if (eqns[[i]]$DV == obsdf[j, "lhs"]) {
              eqns[[i]]$PA   <- c(eqns[[i]]$PA, as.character(obsdf[j, "rhs"]))
              eqns[[i]]$PA   <- eqns[[i]]$PA[eqns[[i]]$PA != ""]
              eqns[[i]]$P    <- eqns[[i]]$PA[eqns[[i]]$PA != eqns[[i]]$DV]
              eqns[[i]]$L    <- c(eqns[[i]]$L, as.character(obsdf[j, "label"]))
              eqns[[i]]$L    <- eqns[[i]]$L[eqns[[i]]$L != ""]
              eqns[[i]]$FIX  <- c(eqns[[i]]$FIX, as.character(obsdf[j, "ustart"]))
              eqns[[i]]$FIX  <- eqns[[i]]$FIX[eqns[[i]]$FIX != ""]
            }
        }
    }
    
    # Add latent var labels
    for (i in 1:length(eqns)) {
       if (eqns[[i]]$EQ == "y1"){
        eqns[[i]]$DV2 <- sapply(eqns[[i]]$DV, function(x) Obs2Lat(x))
        eqns[[i]]$P2 <- unname(sapply(eqns[[i]]$P, function(x) Obs2Lat(x)))
       }
       
       if (eqns[[i]]$EQ == "y2" | eqns[[i]]$EQ == "x2"){
        eqns[[i]]$DV2 <- sapply(eqns[[i]]$P, function(x) Obs2Lat(x))
        eqns[[i]]$P2 <- unname(sapply(eqns[[i]]$P, function(x) Obs2Lat(x)))
       }
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
      if (numLatEnd > 0){ # just added
        for (j in 1:length(EF_x)) {
            for (k in 1:length(eqns[[i]]$PA)) {
                if (eqns[[i]]$PA[k] == names(EF_x)[j]) {
                  eqns[[i]]$EF <- c(eqns[[i]]$EF, EF_x[j])
                }
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
  
  search <- list(eqns = eqns, df = modeqns, constr = constr)
  class(search) <- "miivs"
  search
}
