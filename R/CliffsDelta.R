CliffsDelta <- function(x,
                        y
                        ) {

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
