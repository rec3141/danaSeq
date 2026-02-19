val2col <- function(y,
                    pal   = "viridis",     # name ("viridis","rainbow","Spectral", …) or a function(n)
                    zlim  = range(y, na.rm = TRUE),
                    alpha = 1,
                    nbase = 256,
                    na.color = "#00000000") {
  # resolve palette -> function
  pal_fun <- if (is.character(pal)) {
    base_pals <- c("rainbow","heat.colors","terrain.colors","topo.colors","cm.colors")
    if (pal %in% base_pals) get(pal, mode = "function") else
      function(n) grDevices::hcl.colors(n, palette = pal)
  } else if (is.function(pal)) {
    pal
  } else stop("`pal` must be a palette name or a function")
  
  nbase <- max(2L, as.integer(nbase))
  base_cols <- pal_fun(nbase)
  ramp <- grDevices::colorRamp(base_cols)
  
  # output (pre-fill with NA color)
  out <- rep(na.color, length(y))
  
  # handle degenerate zlim (constant vector → mid color)
  if (!all(is.finite(zlim)) || zlim[2] <= zlim[1]) {
    idx <- is.finite(y)
    mid <- ramp(0.5)
    out[idx] <- grDevices::rgb(mid[1], mid[2], mid[3],
                               alpha = round(pmin(1, pmax(0, alpha)) * 255),
                               maxColorValue = 255)
    return(out)
  }
  
  # finite y only
  idx <- is.finite(y)
  if (any(idx)) {
    p <- (y[idx] - zlim[1]) / (zlim[2] - zlim[1])
    p <- pmin(1, pmax(0, p))  # clamp to [0,1]
    rgbm <- ramp(p)           # numeric matrix, no NAs now
    out[idx] <- grDevices::rgb(rgbm[,1], rgbm[,2], rgbm[,3],
                               alpha = round(pmin(1, pmax(0, alpha)) * 255),
                               maxColorValue = 255)
  }
  out
}
