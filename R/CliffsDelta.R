#' @import bigmemory

CliffsDelta <- function(x,
                        y
                        ) {

  nx <- length(x)
  ny <- length(y)

  # dom <- big.matrix(nrow = nx,
  #                     ncol = ny,
  #                     backingfile = "big_dom_test_backing.bin",
  #                     backingpath = "/home/tbirkne/Desktop/",
  #                     descriptorfile = "big_dom_test_backing.desc"
  # )

  dom <- bigmemory::big.matrix(nrow = nx,
                      ncol = ny,
                      backingfile = ""
                      )
  #dominance matrix implemented as file-backed big.matrix,
    # that can be handled even if it is to big to be stored in memory
  for (i in 1:nx) {
    dom[i, ] <- -sign(y - x[i])
  }

  # computation of cliff's delta
  nxny <- length(dom)
  preSumPSc <- vector(length = ncol(dom))
  preSumBiggers <- vector(length = ncol(dom))
  for (i in 1:ncol(dom)) {
    preSumPSc[i] <- sum(dom[, i] < 0)
    preSumBiggers[i] <- sum(dom[, i] > 0)
  }
  PSc <- sum(preSumPSc)/nxny
  PSbigger <- sum(preSumBiggers)/nxny
  dc <- PSc - PSbigger
  #system("rm /home/tbirkne/Desktop/big_dom_test_backing*")
  gc()
  return(dc)

}
