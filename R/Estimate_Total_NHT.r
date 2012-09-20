Estimate.Total.NHT <- function(MatY.s, VecWk.s, VarEst= "SYG", MatPkl.s= NULL, PopSize= NULL, VecStrataId.s= NULL, VecStrataSizes.H= NULL, ShowStrata= FALSE)
{
  if(  is.vector(MatY.s)             ){tempname <- as.character(deparse(substitute(MatY.s))); MatY.s <- as.matrix(MatY.s, ncol =1); colnames(MatY.s) <- tempname            }
  if(! is.vector(VecWk.s)            ){stop("VecWk.s must be a vector.")                                                                                                         }
  if(any(is.na(VecWk.s))             ){stop("There are missing values in VecWk.s.")                                                                                              }
  if(any(VecWk.s<1)                  ){stop("There are invalid values in VecWk.s, some entries are < 1.")                                                                        }
  if(any(is.na(MatY.s))              ){stop("There are missing values in MatY.s.")                                                                                              }
  n                                   <- dim(MatY.s)[1]
  Q                                   <- dim(MatY.s)[2]
  Nhat                                <- sum(VecWk.s)
  fhat                                <- n/Nhat
  VecPk.s                             <- 1.0/VecWk.s
  if(n != as.integer(length(VecWk.s))){stop("Lengths of columns of MatY.s and the length of VecWk.s are different.")                                                            }
  if(is.null(MatPkl.s)               ){VarEst <- "Hajek"; cat("Note: MatPkl.s is not provided, VarEst is set to \"Hajek\".\n")                                                   }
  if(     VarEst == "Hajek"          ){cat("Note: when using Hajek's approximations, a high-entropy sampling design is assumed.\n")                                              }
  if(!is.null(MatPkl.s))if(any(is.na(MatPkl.s))){stop("There are missing values in MatPkl.s.")                                                                                   }
  if(!is.null(PopSize)               ){f <- n/PopSize                                                                                                                            }
  VecVarName                          <- colnames(MatY.s)
  OUTPUT                              <- data.frame(cbind(Statistic = rep("Total.NHT", times= Q), VariableName = colnames(MatY.s)))
  if(is.null(VecStrataId.s) & is.null(VecStrataSizes.H))
  {
    for(q in (1:Q))
    {
      OUTPUT$Estimate[q]              <- Est.Total.NHT(MatY.s[,q], VecPk.s)
      if(     VarEst == "HT"         ){OUTPUT$Variance[q] <- VE.HT.Total.NHT(   MatY.s[,q], VecPk.s, MatPkl.s)                                                                  }
      else if(VarEst == "SYG"        ){OUTPUT$Variance[q] <- VE.SYG.Total.NHT(  MatY.s[,q], VecPk.s, MatPkl.s)                                                                  }
      else if(VarEst == "Hajek"      ){OUTPUT$Variance[q] <- VE.Hajek.Total.NHT(MatY.s[,q], VecPk.s          )                                                                  }
      else                            {stop("The argument VarEst must be: \"HT\", \"SYG\" or \"Hajek\". If omitted, default is \"SYG\".")                                        }
      OUTPUT$StdErr[q]                <- sqrt(OUTPUT$Variance[q])
      OUTPUT$AbsErr[q]                <- OUTPUT$StdErr[q] * 1.959
      OUTPUT$LInfCI95[q]              <- OUTPUT$Estimate[q] - OUTPUT$AbsErr[q]; if(OUTPUT$LInfCI95[q]<0){OUTPUT$LInfCI95[q] <- 0; cat("Note: LInfCI95 < 0, it is set to zero.\n")}
      OUTPUT$LSupCI95[q]              <- OUTPUT$Estimate[q] + OUTPUT$AbsErr[q]
      OUTPUT$Range95[q]               <- OUTPUT$LSupCI95[q] - OUTPUT$LInfCI95[q]
      OUTPUT$PctCVE[q]                <- round(100.0*OUTPUT$StdErr[q]/OUTPUT$Estimate[q], digits= 3)
      if(!is.null(PopSize)           ){OUTPUT$DEff[q] <- round(OUTPUT$Variance[q]/VE.SYG.Total.NHT(MatY.s[,q],rep(f,   times=n),Pkl.Hajek.s(rep(f,   times=n))), digits = 5)    }
      else                            {OUTPUT$DEff[q] <- round(OUTPUT$Variance[q]/VE.SYG.Total.NHT(MatY.s[,q],rep(fhat,times=n),Pkl.Hajek.s(rep(fhat,times=n))), digits = 5)    }
    }
    OUTPUT$n                          <- n
    OUTPUT$Nhat                       <- round(Nhat, digits= 2)
    OUTPUT$fhat                       <- round(fhat, digits= 4)
    if(!is.null(PopSize)             ){OUTPUT$N <- PopSize; OUTPUT$f <- round(f, digits= 4)                                                                                      }
    cat("\n")
    return(OUTPUT)
  }
  else if(is.null(VecStrataId.s) & !is.null(VecStrataSizes.H)){stop("The argument VecStrataSizes.H is provided but the argument VecStrataId.s is missing.")                      }
  else if(!is.null(VecStrataId.s) & is.null(VecStrataSizes.H)){stop("The argument VecStrataId.s is provided but the argument VecStrataSizes.H is missing.")                      }
  else
  {
    if(is.unsorted(VecStrataId.s)    ){stop("VecStrataId.s must be a sorted vector (of course, all other arguments should follow that order).")                                  }
    if(any(is.na(VecStrataId.s))     ){stop("There are missing values in VecStrataId.s.")                                                                                        }
    if(!is.numeric(VecStrataId.s)    ){stop("The argument VecStrataId.s must be a numeric vector of integers.")                                                                  }
    VecStrataLbls.H                   <- unique(VecStrataId.s)
    H                                 <- length(VecStrataLbls.H)
    if(H != length(VecStrataSizes.H) ){stop("The length of VecStrataSizes.H and the number of different values in VecStrataId.s are different.")                                 }
    if(!is.null(PopSize))if(PopSize != sum(VecStrataSizes.H)){stop("The sum of strata sizes VecStrataSizes.H does not equal PopSize.")                                           }
    Vecnh.H                           <- table(VecStrataId.s)
    Vecfh.H                           <- Vecnh.H/VecStrataSizes.H
    VecWh.H                           <- VecStrataSizes.H/sum(VecStrataSizes.H)
    for(q in (1:Q))
    {
      OUTPUTSTRATA                    <- data.frame(cbind(Statistic = rep("Total.NHT", times= H), VariableName = VecVarName[q]))
      for(h in (1:H))
      {
        OUTPUTSTRATA$h[h]             <- h
        OUTPUTSTRATA$Stratum[h]       <- VecStrataLbls.H[h]
        OUTPUTSTRATA$Estimate[h]      <- Est.Total.NHT(MatY.s[VecStrataId.s == VecStrataLbls.H[h], q], VecPk.s[VecStrataId.s == VecStrataLbls.H[h]])
        if(     VarEst== "HT"        ){OUTPUTSTRATA$Variance[h] <- VE.HT.Total.NHT(   MatY.s[VecStrataId.s== VecStrataLbls.H[h],q], VecPk.s[VecStrataId.s== VecStrataLbls.H[h]],MatPkl.s[VecStrataId.s== VecStrataLbls.H[h],VecStrataId.s== VecStrataLbls.H[h]])}
        else if(VarEst== "SYG"       ){OUTPUTSTRATA$Variance[h] <- VE.SYG.Total.NHT(  MatY.s[VecStrataId.s== VecStrataLbls.H[h],q], VecPk.s[VecStrataId.s== VecStrataLbls.H[h]],MatPkl.s[VecStrataId.s== VecStrataLbls.H[h],VecStrataId.s== VecStrataLbls.H[h]])}
        else if(VarEst== "Hajek"     ){OUTPUTSTRATA$Variance[h] <- VE.Hajek.Total.NHT(MatY.s[VecStrataId.s== VecStrataLbls.H[h],q], VecPk.s[VecStrataId.s== VecStrataLbls.H[h]])}
        else                          {stop("The argument VarEst must be: \"HT\", \"SYG\" or \"Hajek\". If omitted, default is \"SYG\".")                                        }
        OUTPUTSTRATA$StdErr[h]        <- sqrt(OUTPUTSTRATA$Variance[h])
        OUTPUTSTRATA$AbsErr[h]        <- OUTPUTSTRATA$StdErr[h] * 1.959
        OUTPUTSTRATA$LInfCI95[h]      <- OUTPUTSTRATA$Estimate[h] - OUTPUTSTRATA$AbsErr[h]; if(OUTPUTSTRATA$LInfCI95[h]<0){OUTPUTSTRATA$LInfCI95[h] <- 0                         }
        OUTPUTSTRATA$LSupCI95[h]      <- OUTPUTSTRATA$Estimate[h] + OUTPUTSTRATA$AbsErr[h]
        OUTPUTSTRATA$Range95[h]       <- OUTPUTSTRATA$LSupCI95[h] - OUTPUTSTRATA$LInfCI95[h]
        OUTPUTSTRATA$PctCVE[h]        <- round(100.0*OUTPUTSTRATA$StdErr[h]/OUTPUTSTRATA$Estimate[h], digits= 3)
        OUTPUTSTRATA$DEff[h]          <- round(OUTPUTSTRATA$Variance[h]/VE.SYG.Total.NHT(MatY.s[VecStrataId.s== VecStrataLbls.H[h],q],rep(Vecfh.H[h],times=Vecnh.H[h]),Pkl.Hajek.s(rep(Vecfh.H[h],times=Vecnh.H[h]))), digits= 5)
        OUTPUTSTRATA$nh[h]            <- Vecnh.H[h]
        OUTPUTSTRATA$Nh[h]            <- VecStrataSizes.H[h]
        OUTPUTSTRATA$fh[h]            <- round(Vecfh.H[h], digits= 4)
        OUTPUTSTRATA$Wh[h]            <- VecWh.H[h]
      }
      if(ShowStrata                  ){cat("\n"); print(OUTPUTSTRATA, row.names = FALSE)                                                                                         }
      OUTPUT$Estimate[q]              <- sum(OUTPUTSTRATA$Estimate)
      OUTPUT$Variance[q]              <- sum(OUTPUTSTRATA$Variance)
      OUTPUT$StdErr[q]                <- sqrt(OUTPUT$Variance[q])
      OUTPUT$AbsErr[q]                <- OUTPUT$StdErr[q] * 1.959
      OUTPUT$LInfCI95[q]              <- OUTPUT$Estimate[q] - OUTPUT$AbsErr[q]; if(OUTPUT$LInfCI95[q]<0){OUTPUT$LInfCI95[q] <- 0; cat("Note: LInfCI95 < 0, it is set to zero.\n")}
      OUTPUT$LSupCI95[q]              <- OUTPUT$Estimate[q] + OUTPUT$AbsErr[q]
      OUTPUT$Range95[q]               <- OUTPUT$LSupCI95[q] - OUTPUT$LInfCI95[q]
      OUTPUT$PctCVE[q]                <- round(100.0*OUTPUT$StdErr[q]/OUTPUT$Estimate[q], digits = 3)
      if(!is.null(PopSize)           ){OUTPUT$DEff[q] <- round(OUTPUT$Variance[q]/VE.SYG.Total.NHT(MatY.s[,q],rep(f,   times=n),Pkl.Hajek.s(rep(f,   times=n))), digits = 5)    }
      else                            {OUTPUT$DEff[q] <- round(OUTPUT$Variance[q]/VE.SYG.Total.NHT(MatY.s[,q],rep(fhat,times=n),Pkl.Hajek.s(rep(fhat,times=n))), digits = 5)    }
    }
    OUTPUT$n                          <- n
    OUTPUT$Nhat                       <- round(Nhat, digits= 2)
    OUTPUT$fhat                       <- round(fhat, digits= 4)
    if(!is.null(PopSize)             ){OUTPUT$N <- PopSize; OUTPUT$f <- round(f, digits= 4)                                                                                      }
    cat("\n")
    return(OUTPUT)
  }
}
