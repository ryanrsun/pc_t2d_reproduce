library(csmGmm)
library(dplyr)
library(magrittr)
library(data.table)

# input and output
dataFname <- "/rsrch3/home/biostatistics/rsun3/pdac_pleio/data/final_pdac_ukb_diabetes_huerta.txt"
outputDir <- "/rsrch3/home/biostatistics/rsun3/pdac_pleio/output"
outnameRoot <- "pdac_ukb_t2d_huerta"

# controls convergence of EM algorithms
newEps <- 10^(-5)

dat <- fread(dataFname) %>%
  select(Zpdac, Zdiab) %>%
  as.matrix(.)

# adjust so test statistics are not too large for R
for (col_it in 1:2) {
  tooBig <- which(dat[, col_it] > 8.1)
  tooSmall <- which(dat[, col_it] < -8.1)
  if (length(tooBig) > 0) {
    dat[tooBig, col_it] <- 8.1
  }
  if (length(tooSmall) > 0) {
    dat[tooSmall, col_it] <- -8.1
  }
}

# initial parameters
initPiList <- list(c(0.82), c(0.02, 0.02), c(0.02, 0.02), c(0.1))
initMuList <- list(matrix(data=0, nrow=2, ncol=1), matrix(data=c(0, 3, 0, 6), nrow=2),
                   matrix(data=c(3, 0, 6, 0), nrow=2), matrix(data=c(8, 8), nrow=2))

# run csmGmm
newRes <- symm_fit_ind_EM(testStats = dat[, 1:2], sameDirAlt = FALSE, initMuList = initMuList, initPiList = initPiList, eps=newEps)

# save output
setwd(outputDir)
write.table(newRes$lfdrResults, paste0(outnameRoot, "_newlfdr.txt"), append=F, quote=F, row.names=F, col.names=T)
write.table(do.call(cbind, newRes$muInfo), paste0(outnameRoot, "_muInfo.txt"), append=F, quote=F, row.names=F, col.names=T)
write.table(do.call(cbind, newRes$piInfo), paste0(outnameRoot, "_piInfo.txt"), append=F, quote=F, row.names=F, col.names=T)




