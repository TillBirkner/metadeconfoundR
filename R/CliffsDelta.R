CliffsDelta <- function(x,
                        y,
                        fileBackedCliff
                        ) {

  # nx <- length(x)
  # ny <- length(y)
  #
  #
  #
  # # dom <- bigmemory::big.matrix(nrow = nx,
  # #                     ncol = ny,
  # #                     backingfile = ""
  # #                     )
  # #dominance matrix implemented as file-backed big.matrix,
  #   # that can be handled even if it is to big to be stored in memory
  #
  # # 20231204 TB
  # backingfileCompl <- tempfile()
  # backingfileCompl <- path.expand(backingfileCompl)
  # splitTemp <- strsplit(backingfileCompl, "[\\/\\\\]")
  # backingfile <- splitTemp[[1]][length(splitTemp[[1]])]
  # dom <- bigmemory::big.matrix(nrow = nx,
  #                              ncol = ny,
  #                              backingfile = backingfile,
  #                              descriptorfile = paste0(backingfile, ".desc", collapse = ""),
  #                              backingpath = tempdir()
  # )
  # # dom <- c()
  # # if (fileBackedCliff == TRUE) {
  # #   dom <- bigmemory::big.matrix(nrow = nx,
  # #                       ncol = ny,
  # #                       backingfile = ""
  # #                       )
  # # } else {
  # #   dom <- matrix(nrow = nx,
  # #          ncol = ny)
  # # }
  # #
  # #
  # # END 20231204 TB
  # for (i in 1:nx) {
  #   print(i)
  #   dom[i, ] <- -sign(y - x[i])
  # }
  #
  #
  # # computation of cliff's delta
  # nxny <- length(dom)
  # preSumPSc <- vector(length = ncol(dom))
  # preSumBiggers <- vector(length = ncol(dom))
  # for (i in 1:ncol(dom)) {
  #   preSumPSc[i] <- sum(dom[, i] < 0)
  #   preSumBiggers[i] <- sum(dom[, i] > 0)
  # }
  # rm(dom)
  # unlink(paste0(backingfileCompl, "*", collapse = "")) # 20231204 TB
  # PSc <- sum(preSumPSc)/nxny
  # PSbigger <- sum(preSumBiggers)/nxny
  # dc <- PSc - PSbigger
  # return(dc)

  nx <- length(x)
  ny <- length(y)
  nxny <- nx*ny
  preSumPScSum <- 0
  preSumBiggersSum <- 0
  for (i in 1:nx) {
    rowI <- -sign(y - x[i])
    preSumPScSum <- preSumPScSum + sum(rowI < 0)
    preSumBiggersSum <- preSumBiggersSum + sum(rowI > 0)
  }

  PSc <- preSumPScSum/nxny
  PSbigger <- preSumBiggersSum/nxny
  dc <- PSc - PSbigger
  return(dc)
}
