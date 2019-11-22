#' @export
simul.sample <- function(hessian, par, fp, sim = 3000, likdat = NULL, SIR = FALSE, nsir = 50000,
                        with_replacement = TRUE, progress = NULL) {

    log_sum_exp <- function(x) {
        xmax <- which.max(x)
        log1p(sum(exp(x[-xmax] - x[xmax]))) + x[xmax] }

    # We first test is system is singular (TEMPORARY FIX)
    if (any(Matrix::diag(hessian) == 0)) {
        index <- which(Matrix::diag(hessian) == 0)
        for (i in index) { hessian[i, i] <- hessian[i, i] - 1e-2 }
        print('Hessian is exactly singular, Temporary Fix Used')
    }
    # From the hessian, simulate the model
    vcova <- Matrix::solve(-hessian)
    # Test if positive semi-definite
    eS <- eigen(vcova, symmetric = TRUE)
    ev <- eS$values
    if (!all(ev >= -1e-06 * abs(ev[1L]))) {
        vcova <- Matrix::nearPD(vcova, corr = FALSE)$mat
        vcova <- matrix(vcova@x, 
                        nrow = vcova@Dim[1], 
                        ncol = vcova@Dim[2] )
        }

    # Laplace approximation only
    if (SIR == FALSE) {
        # Sample parameters from multivariate normal
            samp <- mvtnorm::rmvnorm(n = sim, par, vcova)
    }

    # Sampling Importance Resampling
    # http://www.sumsar.net/blog/2013/12/shaping_up_laplace_approximation/
    if (SIR == TRUE) {
        if (is.null(likdat)) { print('SIR needs to include <likdat>'); break }
        # Sample parameters from multivariate t-distribution to get thicker tails
        par_sir <- mvtnorm::rmvt(n = nsir, delta = par, sigma = as.matrix(vcova), df = 2)
        # We had the mode to the resampled values (to be sure that it is included in the CI)
        par_sir <- rbind(par_sir, par)
        prp_dens <- mvtnorm::dmvt(par_sir, par, as.matrix(vcova), df = 2, log = TRUE)

        # You might get some warnings, that's OK, the function will return NA
        # and we replace them with -Inf in the llk. We suppress them.
        llk <- suppressWarnings( apply(par_sir, MARGIN = 1, FUN = function(x) {
            if (!is.null(progress)){
                progress()
            }
            ll_hts(x, fp, likdat)
        }))

        llk[!is.finite(llk)] <- -Inf
        dens_ratio <- llk - prp_dens
        wgt <- exp(dens_ratio - log_sum_exp(dens_ratio))
        resampleid <- sample(nrow(par_sir), size = sim, replace = with_replacement, prob = wgt)
        samp <- par_sir[resampleid,, drop = FALSE]
        nunique <- length(unique(resampleid))
        max_wgt <- max(wgt, na.rm = TRUE)
        if (max_wgt > 0.1) { print(paste('Caution, maximum weight is ', round(max_wgt, 2))) }
        print(paste(nunique, 'unique parameter sets resampled'))
    }

    samp
}

#' @export
simul.run <- function(samp, fp, progress = NULL){
    # Define the proper categories over which CI are required
    end_date <- fp$ss$proj_start + fp$ss$PROJ_YEARS - 1L
    diagno <- expand.grid(year = 2000:end_date,
                            outcome = 'aware',
                            agegr = c("15-24", "25-34", '35-49', '15-49', '15+'),
                            sex = c("both", "female", "male"),
                            hivstatus = 'positive')
    out_art <- expand.grid(year = 2000:end_date,
                            outcome = "artcov",
                            agegr = c("15-24", "25-34", '35-49', '15-49', '15+'),
                            sex = c("both", "female", "male"),
                            hivstatus = "positive")
    out_evertest <- expand.grid(year = 2000:end_date,
                            outcome = "evertest",
                            agegr = c("15-24", "25-34", '35-49', '15-49', '15+'),
                            sex = c("both", "female", "male"),
                            hivstatus = c("all", "negative", "positive"))
    out_nbtest <- expand.grid(year = 2000:end_date,
                            outcome = "numbertests",
                            agegrp = "15+",
                            sex = c("both", "female", "male"),
                            hivstatus = 'all')
    out_nbtest_pos <- expand.grid(year = 2000:end_date,
                                outcome = "numbertests",
                                agegrp = "15+",
                                sex = c("both", "female", "male"),
                                hivstatus = 'positive')

    diagno_ss <- add_ss_indices(diagno, fp$ss)
    evertest_ss <- add_ss_indices(out_evertest, fp$ss)
    nbtest_ss <- add_ss_indices(out_nbtest, fp$ss)
    nbtestpos_ss <- add_ss_indices(out_nbtest_pos, fp$ss)
    nbtest_ss$hivstatus <- as.character(nbtest_ss$hivstatus)
    nbtestpos_ss$hivstatus <- as.character(nbtestpos_ss$hivstatus)
    
    # Create parameters (proper scale, etc.), and simulate model
    for (i in 1:nrow(samp)){
        fp <- create_hts_param(samp[i,], fp)
        mod <- simmod(fp)
        diagno[, ncol(diagno) + 1] <- diagnosed(mod, fp, diagno_ss)
        out_evertest[, (ncol(out_evertest) + 1)] <- evertest(mod, fp, evertest_ss)
        out_nbtest[, (ncol(out_nbtest) + 1)] <- total_tests(mod, nbtest_ss)
        out_nbtest_pos[, (ncol(out_nbtest_pos) + 1)] <- total_tests(mod, nbtestpos_ss)
        if (!is.null(progress)){
            progress()
        }
    }

    # Store results in a list of data frames
    list(diagnoses = diagno, ever.test = out_evertest,
                    number.test = out_nbtest, number.test.pos = out_nbtest_pos,
                    param = samp)
}


# Simulation of the dataset
#' @export
simul.test <- function(opt, fp, sim = 3000, likdat = NULL,
                       SIR = FALSE, nsir = 50000, with_replacement = TRUE, 
                       sir_progress = NULL, run_progress = NULL){

    samp <- simul.sample(opt$hessian, opt$par, fp, sim = sim, likdat = likdat,
                        SIR = SIR, nsir = nsir, with_replacement = with_replacement, progress = sir_progress)
    simul.run(samp, fp, progress = run_progress)
}

# Get the confidence interval for the data, first for parameters, second for 
# the data itself
getCI <- function(df) {
  n_col <- ncol(df)
  # Function for the confidence interval of the estimates
  CI <- apply(X = df[,(6:n_col)], MARGIN = 1,
              FUN = quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  lower <- CI[1,]
  upper <- CI[2,]
  df1 <- data.frame(df[c(1:5)], lower, upper)
  return(df1)
}

