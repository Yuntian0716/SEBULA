#' Central-matching doublet detection with FDR control
#'
#' This function performs doublet detection on a dataset by applying Box-Cox transformation,
#' estimating null distributions, computing False Discovery Rate (FDR), and classifying
#' cells as "Doublet" or "Singlet" based on a specified threshold.
#'
#' @param dat A data frame containing at least the following columns:
#'   - `obs`: Observed values to be analyzed.
#'   - `barcode`: Unique cell identifiers.
#' Additional columns are allowed and will be carried through to the output.
#' @param truncation_point Numeric. The threshold below which cells are automatically labeled as "Singlet". Default is `0`.
#' @param pct0 Numeric vector of length 2. The lower and upper quantiles used for null distribution estimation. Default is `c(0.2, 0.6)`.
#' @param nulltype Integer. The type of null estimation model to use:
#'   - `2`: Symmetric model.
#'   - `3`: Asymmetric model allowing different variance estimates for left and right tails.
#' @param thres Numeric. FDR threshold for doublet classification. Default is `0.2`.
#' @param lfdr_method Character. Method used to compute local FDR (`cm.lfdr.truncated`) from the
#'   FDR ordering. Supported options are:
#'   - `"raw"`: Discrete-derivative recovery from cumulative FDR, then clipped to `[0, 1]`.
#'   - `"isoreg"`: Apply isotonic regression (`stats::isoreg`) to the clipped `"raw"` lfdr values
#'     to enforce monotonicity, then use the smoothed values as `cm.lfdr.truncated`.
#'   Default is `"raw"`.
#'
#' @return A list containing:
#'   - `result`: A data frame with the original data and additional columns:
#'       - `box.cox.obs.truncated`: Transformed observations.
#'       - `cm.FDR.truncated`: Computed FDR values.
#'       - `cm.label.truncated`: Predicted classification ("Doublet" or "Singlet").
#'       - `cm.lfdr.truncated`: Local FDR estimates.
#'   - Additional parameters used in the estimation.
#'
#' @import dplyr
#' @import readr
#' @import stats
#' @import ggplot2
#' @import utils
#' @importFrom MASS boxcox
#' @importFrom splines ns
#'
#' @export

doublet_cm <- function(dat, truncation_point = NULL, pct0 = c(0.2, 0.6), nulltype = 2, thres = 0.2, lfdr_method = c("raw", "isoreg")) {
  
  # Ensure 'obs' exists in dataset
  if (!"obs" %in% colnames(dat)) stop("Column 'obs' not found in input data.")
  lfdr_method <- match.arg(lfdr_method)
  ## -------------------------------------------------------------
  ## (0) automatic truncation-point selection via CE
  ## -------------------------------------------------------------
  if (is.null(truncation_point)) {
    q25     <- quantile(dat$obs, 0.25, type = 1)
    cand_tp <- sort(unique(c(0, dat$obs[dat$obs <= q25])))
    
    CE_vec <- vapply(
      cand_tp,
      .cross_entropy_for_tp,
      numeric(1),
      dat      = dat,
      pct0     = pct0,
      nulltype = nulltype
    )
    
    num_cells_retained <- vapply(
      cand_tp,
      function(tp) sum(dat$obs > tp),
      integer(1)
    )
    
    my_df <- data.frame(
      truncation_point = cand_tp,
      CE = CE_vec,
      num_cells_retained = num_cells_retained
    )
    
    ce_plot <- ggplot(my_df, aes(truncation_point, CE)) +
      geom_line(colour = "#3366CC", linewidth = 1) +
      geom_point(colour = "#3366CC", size = 2) +
      geom_text(aes(label = num_cells_retained), vjust = -1,
                size = 3.2, colour = "grey30") +
      labs(
        title = "Cross-entropy trace",
        x = "Candidate truncation point", y = "Cross-entropy"
      ) +
      scale_x_continuous(
        breaks = seq(min(my_df$truncation_point),
                     max(my_df$truncation_point), by = 1)
      ) +
      theme_classic(base_size = 14) +
      theme(plot.title = element_text(hjust = 0.5))
    
    print(ce_plot)
    
    ## ---- early exit so no downstream computation runs ----
    return(invisible(list(ce_table = my_df, ce_plot = ce_plot)))
  }
  
  ## -------------------------------------------------------------
  ## (1) apply truncation filter
  ## -------------------------------------------------------------
  dat$filter  <- dat$obs > truncation_point
  data.filtered <- dplyr::filter(dat, filter)
  
  if (nrow(data.filtered) == 0)
    stop("No barcodes exceed the chosen truncation point.")
  
  ## -------------------------------------------------------------
  ## (2) Box–Cox transform
  ## -------------------------------------------------------------
  x       <- data.filtered$obs
  bc      <- MASS::boxcox(lm(x ~ 1, y = TRUE))
  lambda  <- bc$x[ which.max(bc$y) ]
  newx    <- (x^lambda - 1) / lambda
  data.filtered$new.obs <- newx
  
  ## -------------------------------------------------------------
  ## (3) empirical-null fit (symmetric or asymmetric)
  ## -------------------------------------------------------------
  est <- cm_null_estimate(newx, bre = 120, df = 7,
                          pct0 = pct0, nulltype = nulltype)
  
  # Compute Empirical CDF
  # ecdf_function <- ecdf(newx)
  # ecdf_values <- ecdf_function(newx)
  # result <- data.frame(Observation = newx, ECDF = ecdf_values)
  
  # estimate the mixture density
  brks   <- seq(min(newx), max(newx), length = 120)
  bin_w <- diff(brks)[1]                                   # constant bin width
  cent   <- (brks[-1] + brks[-length(brks)]) / 2        # bin centres
  counts <- hist(newx, breaks = brks, plot = FALSE)$counts
  
  ## Poisson GLM with 7-df natural spline basis (identical to cm_null_estimate)
  mix_fit <- glm(counts ~ splines::ns(cent, df = 7), family = poisson)
  mix_hat <- fitted(mix_fit)                            # fitted mixture counts
  mix_hat_den <- mix_hat / (sum(counts) * bin_w)       # to density for plotting
  mix_mass <- mix_hat/sum(counts)                      # to probability mass per bin
  mix_cdf <- cumsum(mix_mass)
  cdf_fun  <- approxfun(x = cent, y = mix_cdf,
                        yleft = 0, yright = 1, rule = 2)
  
  obs_cdf  <- cdf_fun(newx)
  
  result <- data.frame(
    Observation = newx,
    CDF   = obs_cdf
  )
  ## -------------------- assess the goodness of fit using cross entropy -----------------------------------------
  
  x_flt <- dat$obs[dat$obs > truncation_point]
  if (length(unique(x_flt)) >= 30) {
    bc      <- MASS::boxcox(lm(x_flt ~ 1, y = TRUE))
    lambda  <- bc$x[ which.max(bc$y) ]
    z       <- (x_flt^lambda - 1) / lambda
    
    est     <- cm_null_estimate(z, bre = 120, df = 7,
                                pct0 = pct0,
                                nulltype = nulltype)
    
    mu0 <- est$delta.hat
    if (nulltype == 2) {
      sigma <- est$sigma.hat
      window_vals <- z[z > (mu0 - sigma) & z < (mu0 + sigma)]
      neg_avg_ll <- - mean(dnorm(window_vals, mean = mu0, sd = sigma, log = TRUE))
    } else {
      sigma.left <- est$sigma.left
      sigma.right <- est$sigma.right
      sigma.avg <- (sigma.left + sigma.right) / 2
      window_vals <- z[z > (mu0 - sigma.avg) & z < (mu0 + sigma.avg)]
      loglik_vals <- sapply(window_vals, function(val) {
        if (val > mu0) {
          dnorm(val, mean = mu0, sd = sigma.right, log = TRUE)
        } else {
          dnorm(val, mean = mu0, sd = sigma.left, log = TRUE)
        }
      })
      neg_avg_ll <- - mean(loglik_vals)
    }
    
    message(sprintf("Negative average log-likelihood at chosen truncation point: %.4f", neg_avg_ll))
  }
  
  # Compute FDR based on nulltype
  if (nulltype == 2) {
    mu0 <- est$delta.hat
    sigma <- est$sigma.hat
    pi0.hat <- est$p0
    
    p_x_greater <- 1 - pnorm(newx, mean = mu0, sd = sigma)
    # Visualization of histogram and null fit
    par(
      mar = c(5, 5, 4, 2),
      cex.lab = 1.4,
      cex.axis = 1.2,
      cex.main = 1.5,
      lwd = 2
    )
    
    extra_lines <- c(mu0 - sigma, mu0 + sigma, mu0)
    xrng <- range(c(newx, extra_lines))
    padding <- 0.05 * diff(xrng)
    xlim_range <- c(xrng[1] - padding, xrng[2] + padding)
    window_vals <- newx[newx > (mu0 - sigma) & newx < (mu0 + sigma)]
    
    x_vals <- seq(xlim_range[1], xlim_range[2], length.out = 1000)
    gauss_density <- dnorm(x_vals, mean = mu0, sd = sigma)
    y_null_max <- max(gauss_density, na.rm = TRUE)
    ylim_range <- c(0, y_null_max)
    
    hist(newx,
         breaks = 50,
         probability = TRUE,
         col = "gray93",
         border = NA,
         main = "Observed Distribution with Fitted Component",
         xlab = "Box-Cox Transformed Value",
         ylab = "Density",
         xlim = xlim_range,
         ylim = ylim_range)
    
    col_null   <- "blue"  # blue
    col_mix    <- "#EE553D"  # red–orange
    col_mu     <- "orange"  # purple (vertical μ0 line)
    col_window <- "olivedrab3"  # green  (central-window limits)
    col_window2 <- "purple"  # purple  (evaluation-window limits)
    
    x_vals <- seq(xlim_range[1], xlim_range[2], length.out = 1000)
    gauss_density <- dnorm(x_vals, mean = mu0, sd = ifelse(nulltype == 2, sigma, sigma.right))
    lines(x_vals, gauss_density, col = col_null, lwd = 2.5)
    ## ---- spline-based mixture density (convert counts → density) ----
    lines(cent, mix_hat_den, col = col_mix, lwd = 3, lty = 2)
    grid(nx = NULL, ny = NULL, col = "gray90", lty = "dotted")
    abline(v = mu0, col = col_mu, lty = 2, lwd = 3)
    abline(v = quantile(newx, pct0[1]), col = col_window, lty = 3, lwd = 3)
    abline(v = quantile(newx, pct0[2]), col = col_window, lty = 3, lwd = 3)
    abline(v = mu0 - sigma, col = col_window2, lty = 3, lwd = 3)
    abline(v = mu0 + sigma, col = col_window2, lty = 3, lwd = 3)
    legend("topright",
           legend = c("Null (singlet) density",
                      "Mixture density (spline fit)",
                      expression(hat(mu)[0]),
                      "Central window",
                      "Evaluation window"),
           col    = c(col_null, col_mix, col_mu, col_window, col_window2),
           lty    = c(1, 2, 2, 3, 3),
           lwd    = c(2.5, 3, 3, 3, 3),
           bty = "n", cex = 1.0)
    
  } else {
    mu0 <- est$delta.hat
    sigma.left <- est$sigma.left
    sigma.right <- est$sigma.right
    pi0.hat <- est$p0
    p_x_greater <- ifelse(newx > mu0,
                          1 - pnorm(newx, mean = mu0, sd = sigma.right),
                          1 - pnorm(newx, mean = mu0, sd = sigma.left))
    
    sigma.avg <- (sigma.left + sigma.right) / 2
    extra_lines <- c(mu0 - sigma.avg, mu0 + sigma.avg, mu0,
                     quantile(newx, pct0[1]), quantile(newx, pct0[2]))
    xrng <- range(c(newx, extra_lines))
    padding <- 0.05 * diff(xrng)
    xlim_range <- c(xrng[1] - padding, xrng[2] + padding)
    window_vals <- newx[newx > (mu0 - sigma.avg) & newx < (mu0 + sigma.avg)]
    
    hist(newx,
         breaks = 50,
         probability = TRUE,
         col = "gray93",
         border = NA,
         main = "Observed Distribution with Asymmetric Null Fit",
         xlab = "Box-Cox Transformed Value",
         ylab = "Density",
         xlim = xlim_range)
    
    col_null   <- "blue"  # blue
    col_mix    <- "#EE553D"  # red–orange
    col_mu     <- "orange"  # purple (vertical μ0 line)
    col_window <- "olivedrab3"  # green  (central-window limits)
    
    # Sequence for left and right sides
    x_left <- seq(xlim_range[1], mu0, length.out = 500)
    x_right <- seq(mu0, xlim_range[2], length.out = 500)
    
    # Asymmetric Gaussian curves
    y_left <- dnorm(x_left, mean = mu0, sd = sigma.left)
    y_right <- dnorm(x_right, mean = mu0, sd = sigma.right)
    
    # Add the two density curves
    lines(x_left, y_left, col = col_null, lwd = 2.5)
    lines(x_right, y_right, col = col_null, lwd = 2.5)
    lines(cent, mix_hat_den, col = col_mix, lwd = 3, lty = 2)
    abline(v = mu0, col = col_mu, lty = 2, lwd = 3)
    abline(v = quantile(newx, pct0[1]), col = col_window, lty = 3, lwd = 3)
    abline(v = quantile(newx, pct0[2]), col = col_window, lty = 3, lwd = 3)
    abline(v = mu0 - sigma.avg, col = col_window2, lty = 3, lwd = 3)
    abline(v = mu0 + sigma.avg, col = col_window2, lty = 3, lwd = 3)
    # Add grid and legend
    grid(nx = NULL, ny = NULL, col = "gray90", lty = "dotted")
    legend("topright",
           legend = c("Estimated Asymmetric Null (Singlet)",
                      "Mixture density (spline fit)",
                      expression(hat(mu)[0]),
                      "Central window",
                      "Evaluation window"),
           col    = c(col_null, col_mix, col_mu, col_window, col_window2),
           lty    = c(1, 2, 2, 3, 3),
           lwd    = c(2.5, 3, 3, 3, 3),
           bty = "n", cex = 1.0)
    
  }
  
  pi0.hat <- pmin(1,pi0.hat)
  
  
  ## -------------------------------------------------------------
  ## (4) posterior FDR + labels
  ## -------------------------------------------------------------
  result$singlet.cdf.complement <- p_x_greater
  result$FDR <- (result$singlet.cdf.complement * pi0.hat) / (1 - result$CDF)
  result$FDR[is.na(result$FDR) | is.infinite(result$FDR)] <- 0
  
  # Store results
  data.filtered$box.cox.obs.truncated <- result$Observation
  data.filtered$cm.FDR.truncated <- result$FDR
  data.filtered$cm.label.truncated <- ifelse(result$FDR < thres, "Doublet", "Singlet")
  
  # Merge results back to original dataset
  res.cm <- data.filtered %>%
    dplyr::select(barcode, obs, box.cox.obs.truncated, cm.FDR.truncated, cm.label.truncated)
  
  res <- merge(dat, res.cm, by = c("barcode", "obs"), all.x = TRUE)
  
  # Assign Singlet to truncated values
  res$cm.label.truncated[is.na(res$cm.label.truncated)] <- "Singlet"
  res$cm.FDR.truncated[is.na(res$cm.FDR.truncated)] <- max(res$cm.FDR.truncated, na.rm = TRUE)
  
  res <- res[, !(colnames(res) == "filter")]
  
  #lfdr recovery
  help <- res %>% arrange(cm.FDR.truncated)
  
  FDR <- help$cm.FDR.truncated
  cFDR <- seq_along(FDR) * FDR
  lfdr <- c(FDR[1], sapply(2:length(FDR), function(x) cFDR[x] - cFDR[x - 1]))
  
  # clip to [0,1] (keep your existing behavior)
  lfdr_clipped <- pmin(1, pmax(0, lfdr))
  
  # OPTIONAL isotonic regression smoothing (new option)
  if (lfdr_method == "isoreg") {
    # use ranks 1..n as x
    ir <- stats::isoreg(seq_along(lfdr_clipped), lfdr_clipped)
    lfdr_final <- ir$yf
    # ensure still within [0,1] after numerical issues
    lfdr_final <- pmin(1, pmax(0, lfdr_final))
  } else {
    lfdr_final <- lfdr_clipped
  }
  
  # Preserve your existing merge-back behavior (barcode key)
  lfdr_df <- data.frame(
    barcode = help$barcode,
    cm.lfdr.truncated = lfdr_final
  )
  res <- res %>% dplyr::left_join(lfdr_df, by = "barcode")
  
  # Corrected pi0 estimate
  pi0.correct <- (pi0.hat * length(newx) + sum(dat$obs <= truncation_point)) / length(dat$obs)
  
  # Count doublets and print the summary
  doublet_count <- sum(res$cm.label.truncated == "Doublet")
  total_cells <- nrow(res)
  doublet_proportion <- doublet_count / total_cells
  
  cat("Number of doublets called:", doublet_count, "\n")
  cat("Proportion of doublets:", round(doublet_proportion * 100, 2), "%\n")
  
  # Return final dataset
  if (nulltype == 3){
    
    return(list(result = res, pi0.hat.all = pi0.correct, pi0.hat = pi0.hat, mu0.hat = mu0, sigma.left = sigma.left,  sigma.right =  sigma.right))
    
  } else{
    
    return(list(result = res, pi0.hat.all = pi0.correct, pi0.hat = pi0.hat, mu0.hat = mu0, sigma.hat = sigma))
    
  }
  
}


#' Estimate the null distribution using central matching
#'
#' This function estimates the null distribution of a dataset by performing truncation,
#' Poisson regression, and either a symmetric or asymmetric quadratic approximation.
#'
#' @param zz Numeric vector. The dataset to be analyzed.
#' @param bre Integer. Number of breaks for histogram binning.
#' @param df Integer. Degrees of freedom for natural splines in Poisson regression.
#' @param pct0 Numeric vector of length 2. Defines the lower and upper quantile thresholds for selecting data used in null estimation.
#' @param nulltype Integer. The type of null estimation model:
#'   - `2`: Symmetric model using quadratic fitting.
#'   - `3`: Asymmetric model allowing different variance estimates for left and right tails.
#'
#' @return A list containing:
#'   - `delta.hat`: Estimated mode of the null distribution (for nulltype = 2).
#'   - `sigma.hat`: Estimated standard deviation of the null (for nulltype = 2).
#'   - `sigma.left`, `sigma.right`: Estimated left and right standard deviations (for nulltype = 3).
#'   - `p0`: Estimated null proportion for normalization.
#'
#' @export

cm_null_estimate <- function(zz, bre = 120, df = 7, pct0 = c(0.2, 0.6), nulltype) {

  # Set truncation limits
  lo <- min(zz)
  up <- max(zz)

  # Truncate values within the selected range
  zzz <- pmax(pmin(zz, up), lo)

  # Compute histogram
  breaks <- seq(lo, up, length = bre)
  zh <- hist(zzz, breaks = breaks, plot = FALSE)
  x <- (breaks[-1] + breaks[-length(breaks)]) / 2
  yall <- y <- zh$counts
  K <- length(y)
  N <- length(zz)

  # Poisson regression fit using natural splines
  X <- cbind(1, splines::ns(x, df = df))
  f <- glm(y ~ splines::ns(x, df = df), family = poisson)$fitted.values
  l <- log(f)

  # Compute central matching estimation
  imax <- which.max(l)
  xmax <- x[imax]

  # Define quantile thresholds
  pctlo <- pct0[1]
  pctup <- pct0[2]

  # Select data for null distribution estimation
  lo0 <- quantile(zz, pctlo)
  hi0 <- quantile(zz, pctup)
  nx <- length(x)
  i0 <- which(x > lo0 & x < hi0)
  x0 <- x[i0]
  y0 <- l[i0]

  ## If nulltype = 3 (Asymmetric model)
  if (nulltype == 3) {
    X00 <- cbind((x0 - xmax)^2, pmax(x0 - xmax, 0)^2)
    lr <- lm(y0 ~ X00)
    co <- lr$coef

    # Check if quadratic fit is valid
    cm_error <- is.na(co[3]) || co[2] >= 0 || (co[2] + co[3] >= 0)

    if (cm_error) {
      stop("CM estimation failed: Consider using nulltype = 2.")
    }

    X0 <- cbind(1, (x - xmax)^2, pmax(x - xmax, 0)^2)
    sigs <- 1 / sqrt(-2 * c(co[2], co[2] + co[3]))

    delta <- xmax
    sigma.left <- sigs[1]
    sigma.right <- sigs[2]

    l0 <- as.vector(X0 %*% co)
    f0 <- exp(l0)
    p0 <- sum(f0) / sum(f)
    f0 <- f0 / p0

    result <- list(delta.hat = delta, sigma.left = sigma.left, sigma.right = sigma.right, p0 = p0)

  } else if (nulltype == 2) {
    ## If nulltype = 2 (Symmetric model)
    X00 <- cbind(x0 - xmax, (x0 - xmax)^2)
    lr <- lm(y0 ~ X00)
    co <- lr$coef
    X0 <- cbind(1, x - xmax, (x - xmax)^2)

    # Compute peak and standard deviation
    xmaxx <- -co[2] / (2 * co[3]) + xmax
    sighat <- 1 / sqrt(-2 * co[3])

    # Compute null distribution function
    l0 <- as.vector(X0 %*% co)
    f0 <- exp(l0)
    p0 <- sum(f0) / sum(f)
    f0 <- f0 / p0

    result <- list(delta.hat = xmaxx, sigma.hat = sighat, p0 = p0)
  } else {
    stop("Invalid nulltype. Choose 2 or 3.")
  }


  return(result)
}

.cross_entropy_for_tp <- function(tp, dat, pct0, nulltype) {
  # Filter observations above truncation point
  x_flt <- dat$obs[dat$obs > tp]
  if (length(unique(x_flt)) < 30) return(Inf)

  # Box-Cox transform
  bc      <- MASS::boxcox(lm(x_flt ~ 1, y = TRUE))
  lambda  <- bc$x[which.max(bc$y)]
  z       <- (x_flt^lambda - 1) / lambda

  # Estimate null
  est <- cm_null_estimate(z, bre = 120, df = 7,
                          pct0 = pct0,
                          nulltype = nulltype)

  # Define evaluation window: [mu - sigma, mu + sigma]
  if (nulltype == 2) {
    mu0 <- est$delta.hat
    sigma <- est$sigma.hat
  } else {
    mu0 <- est$delta.hat
    # Approximate by averaging left and right sigma
    sigma <- (est$sigma.left + est$sigma.right) / 2
  }

  lower <- mu0 - sigma
  upper <- mu0 + sigma

  # Values in evaluation window
  z_eval <- z[z > lower & z < upper]

  # Avoid too few points
  if (length(z_eval) < 30) return(Inf)

  # Negative average log-likelihood (cross-entropy)
  if (nulltype == 2) {
    ll <- sum(dnorm(z_eval, mean = mu0, sd = sigma, log = TRUE))
  } else {
    ll <- sum(dnorm(z_eval, mean = mu0,
                    sd = ifelse(z_eval > mu0, est$sigma.right, est$sigma.left),
                    log = TRUE))
  }

  n <- length(z_eval)
  neg_avg_loglik <- - ll / n

  return(neg_avg_loglik)
}


