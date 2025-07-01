
# hima_survival <- function(X, M, OT, status, COV = NULL,
#                           topN = NULL,
#                           scale = TRUE,
#                           FDRcut = 0.05,
#                           verbose = FALSE) {
#   X <- matrix(X, ncol = 1)
#   M <- as.matrix(M)
# 
#   M_ID_name <- colnames(M)
#   if (is.null(M_ID_name)) M_ID_name <- seq_len(ncol(M))
# 
#   n <- nrow(M)
#   p <- ncol(M)
# 
#   if (is.null(COV)) {
#     q <- 0
#     MZ <- cbind(M, X)
#   } else {
#     COV <- as.matrix(COV)
#     q <- dim(COV)[2]
#     MZ <- cbind(M, COV, X)
#   }
# 
#   MZ <- process_var(MZ, scale)
#   if (scale && verbose) message("Data scaling is completed.")
# 
#   #########################################################################
#   ################################ STEP 1 #################################
#   #########################################################################
#   message("Step 1: Sure Independent Screening ...", "     (", format(Sys.time(), "%X"), ")")
# 
#   if (is.null(topN)) d_0 <- ceiling(n / log(n)) else d_0 <- topN # the number of top mediators that associated with exposure (X)
#   d_0 <- min(p, d_0) # if d_0 > p select all mediators
# 
#   beta_SIS <- matrix(0, 1, p)
# 
#   for (i in 1:p) {
#     ID_S <- c(i, (p + 1):(p + q + 1))
#     MZ_SIS <- MZ[, ID_S]
#     fit <- survival::coxph(survival::Surv(OT, status) ~ MZ_SIS)
#     beta_SIS[i] <- fit$coefficients[1]
#   }
# 
#   alpha_SIS <- matrix(0, 1, p)
#   XZ <- cbind(X, COV)
#   for (i in 1:p) {
#     fit_a <- lsfit(XZ, M[, i], intercept = TRUE)
#     est_a <- matrix(coef(fit_a))[2]
#     alpha_SIS[i] <- est_a
#   }
# 
#   ab_SIS <- alpha_SIS * beta_SIS
#   ID_SIS <- which(-abs(ab_SIS) <= sort(-abs(ab_SIS))[min(p, d_0)])
# 
#   d <- length(ID_SIS)
# 
#   if (verbose) message("        Top ", d, " mediators are selected: ", paste0(M_ID_name[ID_SIS], collapse = ", "))
# 
#   #########################################################################
#   ################################ STEP 2 #################################
#   #########################################################################
#   message("Step 2: De-biased Lasso estimates ...", "     (", format(Sys.time(), "%X"), ")")
# 
#   if (verbose) {
#     if (is.null(COV)) {
#       message("        No covariate was adjusted.")
#     } else {
#       message("        Adjusting for covariate(s): ", paste0(colnames(COV), collapse = ", "))
#     }
#   }
# 
#   ## estimation of beta
#   P_beta_SIS <- matrix(0, 1, d)
#   beta_DLASSO_SIS_est <- matrix(0, 1, d)
#   beta_DLASSO_SIS_SE <- matrix(0, 1, d)
#   MZ_SIS <- MZ[, c(ID_SIS, (p + 1):(p + q + 1))]
#   MZ_SIS_1 <- t(t(MZ_SIS[, 1]))
# 
#   for (i in 1:d) {
#     V <- MZ_SIS
#     V[, 1] <- V[, i]
#     V[, i] <- MZ_SIS_1
#     LDPE_res <- LDPE_func(ID = 1, X = V, OT = OT, status = status)
#     beta_LDPE_est <- LDPE_res[1]
#     beta_LDPE_SE <- LDPE_res[2]
#     V1_P <- abs(beta_LDPE_est) / beta_LDPE_SE
#     P_beta_SIS[i] <- 2 * (1 - pnorm(V1_P, 0, 1))
#     beta_DLASSO_SIS_est[i] <- beta_LDPE_est
#     beta_DLASSO_SIS_SE[i] <- beta_LDPE_SE
#   }
# 
#   ## estimation of alpha
#   alpha_SIS_est <- matrix(0, 1, d)
#   alpha_SIS_SE <- matrix(0, 1, d)
#   P_alpha_SIS <- matrix(0, 1, d)
#   XZ <- cbind(X, COV)
# 
#   for (i in 1:d) {
#     fit_a <- lsfit(XZ, M[, ID_SIS[i]], intercept = TRUE)
#     est_a <- matrix(coef(fit_a))[2]
#     se_a <- ls.diag(fit_a)$std.err[2]
#     sd_1 <- abs(est_a) / se_a
#     P_alpha_SIS[i] <- 2 * (1 - pnorm(sd_1, 0, 1)) ## the SIS for alpha
#     alpha_SIS_est[i] <- est_a
#     alpha_SIS_SE[i] <- se_a
#   }
# 
#   #########################################################################
#   ################################ STEP 3 #################################
#   #########################################################################
#   message("Step 3: Multiple-testing procedure ...", "     (", format(Sys.time(), "%X"), ")")
# 
#   PA <- cbind(t(P_alpha_SIS), t(P_beta_SIS))
#   P_value <- apply(PA, 1, max) # the joint p-values for SIS variable
# 
#   ## the multiple-testing  procedure
#   N0 <- dim(PA)[1] * dim(PA)[2]
# 
#   input_pvalues <- PA + matrix(runif(N0, 0, 10^{
#     -10
#   }), dim(PA)[1], 2)
#   nullprop <- null_estimation(input_pvalues, lambda = 0.5)
#   fdrcut <- HDMT::fdr_est(nullprop$alpha00,
#                           nullprop$alpha01,
#                           nullprop$alpha10,
#                           nullprop$alpha1,
#                           nullprop$alpha2,
#                           input_pvalues,
#                           exact = 0
#   )
# 
#   ID_fdr <- which(fdrcut <= FDRcut)
# 
#   IDE <- alpha_SIS_est[ID_fdr] * beta_DLASSO_SIS_est[ID_fdr]
# 
#   if (length(ID_fdr) > 0) {
#     out_result <- data.frame(
#       Index = M_ID_name[ID_SIS][ID_fdr],
#       alpha_hat = alpha_SIS_est[ID_fdr],
#       alpha_se = alpha_SIS_SE[ID_fdr],
#       beta_hat = beta_DLASSO_SIS_est[ID_fdr],
#       beta_se = beta_DLASSO_SIS_SE[ID_fdr],
#       IDE = IDE,
#       rimp = abs(IDE) / sum(abs(IDE)) * 100,
#       pmax = P_value[ID_fdr]
#     )
#     if (verbose) message(paste0("        ", length(ID_fdr), " significant mediator(s) identified."))
#   } else {
#     if (verbose) message("        No significant mediator identified.")
#     out_result <- NULL
#   }
# 
#   message("Done!", "     (", format(Sys.time(), "%X"), ")")
# 
#   return(out_result)
# }
# 
# 
# # Internal function: parallel computing check
# 
# checkParallel <- function(program.name, parallel, ncore, verbose) {
#   if (parallel && (ncore > 1)) {
#     if (ncore > parallel::detectCores()) {
#       message(
#         "You requested ", ncore, " cores. There are only ",
#         parallel::detectCores(), " in your machine!"
#       )
#       ncore <- parallel::detectCores()
#     }
#     if (verbose) {
#       message(
#         "    Running ", program.name, " with ", ncore, " cores in parallel...   (",
#         format(Sys.time(), "%X"), ")"
#       )
#     }
#     doParallel::registerDoParallel(ncore)
#   } else {
#     if (verbose) {
#       message(
#         "    Running ", program.name, " with single core...   (",
#         format(Sys.time(), "%X"), ")"
#       )
#     }
#     registerDoSEQ()
#   }
# }
# 
# 
# 
# ## Internal function: doOne code generater
# 
# doOneGen <- function(model.text, colind.text) {
#   L <- length(eval(parse(text = colind.text)))
#   script <- paste0(
#     "doOne <- function(i, datarun, Ydat){datarun$Mone <- Ydat[,i]; model <- ",
#     model.text, ";if('try-error' %in% class(model)) b <- rep(NA, ",
#     L, ") else { res=summary(model)$coefficients; b <- res[2,", colind.text,
#     "]};invisible(b)}"
#   )
#   return(script)
# }
# 
# 
# 
# ## Internal function: create iterator for bulk matrix by column
# 
# iblkcol_lag <- function(M, ...) {
#   i <- 1
#   it <- iterators::idiv(ncol(M), ...)
# 
#   nextEl <- function() {
#     n <- iterators::nextElem(it)
#     r <- seq(i, length = n)
#     i <<- i + n
#     M[, r, drop = FALSE]
#   }
#   obj <- list(nextElem = nextEl)
#   class(obj) <- c("abstractiter", "iter")
#   obj
# }
# 
# 
# 
# # Internal function: Sure Independent Screening for hima
# 
# globalVariables("n")
# globalVariables("M_chunk")
# 
# himasis <- function(Y, M, X, COV, glm.family, modelstatement,
#                     parallel, ncore, verbose, tag) {
#   L.M <- ncol(M)
#   M.names <- colnames(M)
# 
#   X <- data.frame(X)
#   X <- data.frame(model.matrix(~., X))[, -1]
# 
#   if (is.null(COV)) {
#     if (verbose) message("    No covariate is adjusted")
#     datarun <- data.frame(Y = Y, Mone = NA, X = X)
#     modelstatement <- modelstatement
#   } else {
#     COV <- data.frame(COV)
#     COV <- data.frame(model.matrix(~., COV))[, -1, drop = FALSE]
#     conf.names <- colnames(COV)
#     if (verbose) message("    Adjusting for covariate(s): ", paste0(conf.names, collapse = ", "))
#     datarun <- data.frame(Y = Y, Mone = NA, X = X, COV)
#     modelstatement <- eval(parse(text = (paste0(
#       modelstatement, "+",
#       paste0(colnames(COV), collapse = "+")
#     ))))
#   }
# 
#   if (glm.family == "gaussian") {
#     doOne <- eval(parse(text = doOneGen(paste0(
#       "try(glm(modelstatement, family = ",
#       glm.family, ", data = datarun))"
#     ), "c(1,4)")))
#   } else if (glm.family == "negbin") {
#     doOne <- eval(parse(text = doOneGen(paste0("try(MASS::glm.nb(modelstatement, data = datarun))"), "c(1,4)")))
#   } else {
#     stop(paste0("Screening family ", glm.family, " is not supported."))
#   }
# 
# 
#   checkParallel(tag, parallel, ncore, verbose)
# 
#   results <- foreach(
#     n = iterators::idiv(L.M, chunks = ncore),
#     M_chunk = iblkcol_lag(M, chunks = ncore),
#     .combine = "cbind"
#   ) %dopar% {
#     sapply(seq_len(n), doOne, datarun, M_chunk)
#   }
# 
#   colnames(results) <- M.names
#   return(results)
# }
# 
# 
# 
# # Internal function: process_var
# # Helper function to process variables
# 
# process_var <- function(var, scale) {
#   if (!is.null(var)) {
#     if (scale) {
#       return(scale(var))
#     } else {
#       return(as.matrix(var))
#     }
#   } else {
#     return(NULL)
#   }
# }
# 
# 
# 
# # Internal function: LDPE (optimized)
# # the code of Liu han's JRSSB paper for high-dimensional Cox model
# # ID: the index of interested parameter
# # X: the covariates matrix with n by p
# # OT: the observed time = min(T,C)
# # status: the censoring indicator I(T <= C)
# 
# LDPE_func <- function(ID, X, OT, status) {
#   set.seed(1029)
# 
#   # Dimensions
#   coi <- ID
#   d <- ncol(X)
#   n <- nrow(X)
# 
#   # Initialize penalty factor
#   PF <- rep(1, d)
#   PF[ID] <- 1
# 
#   # Semi-penalized initial estimator
#   fit <- glmnet(X, survival::Surv(OT, status), family = "cox", alpha = 1, standardize = FALSE, penalty.factor = PF)
#   cv.fit <- cv.glmnet(X, survival::Surv(OT, status), family = "cox", alpha = 1, standardize = FALSE, penalty.factor = PF)
#   betas <- as.vector(coef(fit, s = cv.fit$lambda.min))
# 
#   # Precompute sorted times and order
#   stime <- sort(OT)
#   otime <- order(OT)
# 
#   # Preallocate storage
#   la <- numeric(n)
#   lb <- matrix(0, n, d - 1)
#   Hs <- matrix(0, d, d)
# 
#   # Iterate through observations
#   for (i in seq_len(n)) {
#     if (status[otime[i]] == 1) {
#       ind <- which(OT >= stime[i])
#       X_ind <- X[ind, , drop = FALSE]
#       tmp_exp <- as.vector(exp(X_ind %*% betas))
# 
#       S0 <- sum(tmp_exp)
#       S1 <- colSums(X_ind * tmp_exp)
#       S2 <- t(X_ind) %*% (X_ind * tmp_exp)
# 
#       la[i] <- -(X[otime[i], coi] - S1[coi] / S0)
# 
#       if (coi == 1) {
#         lb[i, ] <- -(X[otime[i], -1] - S1[-1] / S0)
#       } else if (coi == d) {
#         lb[i, ] <- -(X[otime[i], -d] - S1[-d] / S0)
#       } else {
#         lb[i, ] <- -(X[otime[i], -c(coi)] - S1[-c(coi)] / S0)
#       }
# 
#       Hs <- Hs + (S0 * S2 - tcrossprod(S1)) / (n * S0^2)
#     }
#   }
# 
#   # De-biased Lasso step
#   fit_res <- glmnet(lb, la, alpha = 1, standardize = FALSE, intercept = FALSE, lambda = sqrt(log(d) / n))
#   what <- as.vector(coef(fit_res)[-1])
# 
#   # Final estimate and variance
#   if (coi == 1) {
#     S <- betas[coi] - (mean(la) - crossprod(what, colMeans(lb))) / (Hs[coi, coi] - crossprod(what, Hs[-1, coi]))
#     var <- Hs[coi, coi] - crossprod(what, Hs[-1, coi])
#   } else if (coi == d) {
#     S <- betas[coi] - (mean(la) - crossprod(what, colMeans(lb))) / (Hs[coi, coi] - crossprod(what, Hs[-d, coi]))
#     var <- Hs[coi, coi] - crossprod(what, Hs[-d, coi])
#   } else {
#     S <- betas[coi] - (mean(la) - crossprod(what, colMeans(lb))) / (Hs[coi, coi] - crossprod(what, Hs[-c(coi), coi]))
#     var <- Hs[coi, coi] - crossprod(what, Hs[-c(coi), coi])
#   }
# 
#   beta_est <- S
#   beta_SE <- sqrt(1 / (n * var))
# 
#   result <- c(beta_est, beta_SE)
#   return(result)
# }
# 
# 
# # Internal function: null_estimation
# # A function to estimate the proportions of the three component nulls
# # This is from HDMT package version < 1.0.4 (optimized here)
# 
# null_estimation <- function(input_pvalues, lambda = 0.5) {
#   ## Validate input
#   if (is.null(ncol(input_pvalues)) || ncol(input_pvalues) != 2) {
#     stop("`input_pvalues` must be a matrix or data frame with exactly 2 columns.")
#   }
#   input_pvalues <- as.matrix(input_pvalues)
#   if (anyNA(input_pvalues)) {
#     warning("`input_pvalues` contains NAs, which will be removed.")
#     input_pvalues <- input_pvalues[complete.cases(input_pvalues), ]
#   }
#   if (nrow(input_pvalues) == 0) {
#     stop("`input_pvalues` does not contain valid rows.")
#   }
# 
#   ## Precompute threshold values
#   pcut <- seq(0.1, 0.8, 0.1)
#   one_minus_pcut <- 1 - pcut
#   one_minus_pcut_sq <- one_minus_pcut^2
# 
#   ## Calculate fractions using vectorized operations
#   frac1 <- colMeans(outer(input_pvalues[, 1], pcut, `>=`)) / one_minus_pcut
#   frac2 <- colMeans(outer(input_pvalues[, 2], pcut, `>=`)) / one_minus_pcut
#   frac12 <- colMeans(outer(input_pvalues[, 2], pcut, `>=`) & outer(input_pvalues[, 1], pcut, `>=`)) / one_minus_pcut_sq
# 
#   ## Estimate alpha00
#   alpha00 <- min(frac12[pcut == lambda], 1)
# 
#   ## Estimate alpha1 and alpha2
#   alpha1 <- if (stats::ks.test(input_pvalues[, 1], "punif", 0, 1, alternative = "greater")$p > 0.05) 1 else min(frac1[pcut == lambda], 1)
#   alpha2 <- if (stats::ks.test(input_pvalues[, 2], "punif", 0, 1, alternative = "greater")$p > 0.05) 1 else min(frac2[pcut == lambda], 1)
# 
#   ## Estimate other proportions
#   if (alpha00 == 1) {
#     alpha01 <- alpha10 <- alpha11 <- 0
#   } else {
#     alpha01 <- alpha1 - alpha00
#     alpha10 <- alpha2 - alpha00
#     alpha01 <- max(0, alpha01)
#     alpha10 <- max(0, alpha10)
# 
#     if (alpha1 == 1 && alpha2 == 1) {
#       alpha00 <- 1
#       alpha01 <- alpha10 <- alpha11 <- 0
#     } else if ((1 - alpha00 - alpha01 - alpha10) < 0) {
#       alpha11 <- 0
#       alpha10 <- 1 - alpha1
#       alpha01 <- 1 - alpha2
#       alpha00 <- 1 - alpha10 - alpha01
#     } else {
#       alpha11 <- 1 - alpha00 - alpha01 - alpha10
#     }
#   }
# 
#   ## Return results
#   list(alpha10 = alpha10, alpha01 = alpha01, alpha00 = alpha00, alpha1 = alpha1, alpha2 = alpha2)
# }
# 
# 
# # Internal function: DLASSO_fun
# # A function perform de-biased lasso estimator used by function "hima_microbiome"
# 
# DLASSO_fun <- function(X, Y) {
#   set.seed(1029)
# 
#   n <- dim(X)[1]
#   p <- dim(X)[2]
#   fit <- glmnet(X, Y, alpha = 1)
#   cv.fit <- cv.glmnet(X, Y, alpha = 1)
#   beta_0 <- coef(fit, s = cv.fit$lambda.min)[2:(p + 1)]
#   #
#   fit <- glmnet(X[, 2:p], X[, 1], alpha = 1)
#   cv.fit <- cv.glmnet(X[, 2:p], X[, 1], alpha = 1)
#   phi_hat <- coef(fit, s = cv.fit$lambda.min)[2:p]
#   ##
#   R <- X[, 1] - X[, 2:p] %*% t(t(phi_hat))
#   E <- Y - X %*% t(t(beta_0))
#   beta_1_hat <- beta_0[1] + sum(R * E) / sum(R * X[, 1]) #  The de-biased lasso estimator
#   ##
#   sigma_e2 <- sum(E^2) / (n - length(which(beta_0 != 0)))
# 
#   sigma_beta1_hat <- sqrt(sigma_e2) * sqrt(sum(R^2)) / abs(sum(R * X[, 1]))
# 
#   results <- c(beta_1_hat, sigma_beta1_hat)
#   return(results)
# }
# 
# # Internal function: rdirichlet
# # A function generate random number from Dirichlet distribution.
# 
# rdirichlet <- function(n = 1, alpha) {
#   Gam <- matrix(0, n, length(alpha))
#   for (i in seq_along(alpha)) Gam[, i] <- stats::rgamma(n, shape = alpha[i])
#   Gam / rowSums(Gam)
# }
# 
# 
# 
# # Internal function: DACT (optimized)
# # A function to perform Divide-Aggregate Composite-null Test (DACT) by Liu et al. (2020).
# # p value is corrected by JC method Jin and Cai (2007).
# # This function is used in hima_efficient
# 
# DACT <- function(p_a, p_b) {
#   Z_a <- stats::qnorm(p_a, lower.tail = FALSE)
#   Z_b <- stats::qnorm(p_b, lower.tail = FALSE)
#   pi0a <- 1 - .nonnullPropEst(Z_a, 0, 1)
#   pi0b <- 1 - .nonnullPropEst(Z_b, 0, 1)
#   pi0a <- min(pi0a, 1)
#   pi0b <- min(pi0b, 1)
# 
#   p3 <- (pmax(p_a, p_b))^2
#   wg1 <- pi0a * (1 - pi0b)
#   wg2 <- (1 - pi0a) * pi0b
#   wg3 <- pi0a * pi0b
#   wg.sum <- wg1 + wg2 + wg3
#   wg.std <- c(wg1, wg2, wg3) / wg.sum
# 
#   p_dact <- wg.std[1] * p_a + wg.std[2] * p_b + wg.std[3] * p3
#   p_dact <- .JCCorrect(p_dact)
#   return(p_dact)
# }
# 
# .JCCorrect <- function(pval) {
#   z <- stats::qnorm(pval, lower.tail = FALSE)
#   res <- .nullParaEst(z)
#   stats::pnorm(z, mean = res$mu, sd = res$s, lower.tail = FALSE)
# }
# 
# .nonnullPropEst <- function(x, u, sigma) {
#   z <- (x - u) / sigma
#   xi <- seq(0, 1, length.out = 101)
#   tmax <- sqrt(log(length(x)))
#   tt <- seq(0, tmax, 0.1)
# 
#   weights <- 1 - abs(xi)  # Weights based on xi
#   epsest <- numeric(length(tt))  # Preallocate results
# 
#   for (j in seq_along(tt)) {
#     t <- tt[j]
#     f <- exp((t * xi)^2 / 2)
#     co <- rowMeans(cos(outer(t * xi, z, `*`)))  # Calculate cosine terms
#     epsest[j] <- 1 - sum(weights * f * co) / sum(weights)
#   }
# 
#   max(epsest)
# }
# 
# .nullParaEst <- function(x, gamma = 0.1) {
#   n <- length(x)
#   t <- seq(0.005, 5, length.out = 1000)  # Define range for smoother spacing
# 
#   gan <- n^(-gamma)
# 
#   # Compute cos(s * x) and sin(s * x) using outer for vectorization
#   cos_vals <- outer(t, x, FUN = function(t, x) cos(t * x))  # Matrix of size (length(t), length(x))
#   sin_vals <- outer(t, x, FUN = function(t, x) sin(t * x))  # Matrix of size (length(t), length(x))
# 
#   # Compute phi and its derivatives
#   phiplus <- rowMeans(cos_vals)  # Mean of each row
#   phiminus <- rowMeans(sin_vals)  # Mean of each row
#   dphiplus <- -rowMeans(sweep(sin_vals, 2, x, `*`))  # Broadcasting x across columns
#   dphiminus <- rowMeans(sweep(cos_vals, 2, x, `*`))  # Broadcasting x across columns
#   phi <- sqrt(phiplus^2 + phiminus^2)  # Magnitude of phiplus and phiminus
# 
#   # Find the first index where phi - gan <= 0
#   ind <- which((phi - gan) <= 0)[1]
# 
#   if (is.na(ind)) {
#     stop("Unable to find a suitable index where phi - gan <= 0.")
#   }
# 
#   # Extract values at the identified index
#   tt <- t[ind]
#   a <- phiplus[ind]
#   b <- phiminus[ind]
#   da <- dphiplus[ind]
#   db <- dphiminus[ind]
#   c <- phi[ind]
# 
#   # Compute final estimates
#   shat <- sqrt(-(a * da + b * db) / (tt * c^2))
#   uhat <- -(da * b - db * a) / (c^2)
#   epshat <- 1 - c * exp((tt * shat)^2 / 2)
# 
#   list(mu = uhat, s = shat)
# }

########################simulation studies#######################
library(Matrix)
library(glmnet)
library(ncvreg)
library(R.matlab)
library(future.apply)
library(rhdf5)
library(HIMA)

setwd("D:/simulation_coding")

N <- 200
h5_file <- H5Fopen("simulated_data.mat")

file_path <- "simulated_data.mat"
all_nodes <- h5ls(file_path)

main_groups <- all_nodes[grepl("^/#refs#/(\\d+|[A-Za-z]{1,2})", all_nodes$group), ]
main_group_ids <- unique(gsub("^/#refs#/(\\d+|[A-Za-z]{1,2}).*", "\\1", main_groups$group))
main_groups <- all_nodes[grepl("^/#refs#/(\\d+[A-Za-z]|\\d+|[A-Za-z]{1,2})", all_nodes$group), ]

main_group_ids <- unique(
  gsub("^/#refs#/((\\d+[A-Za-z])|(\\d+)|([A-Za-z]{1,2})).*", "\\1", main_groups$group)
)

datasets <- c("T", "X", "Z","M", "Delta", "ID_M")

all_data <- lapply(main_group_ids, function(id) {
  group_path <- paste0("/#refs#/", id)
  data <- lapply(datasets, function(ds) {
    path <- paste0(group_path, "/", ds)
    h5read(file_path, path)
  })
  names(data) <- datasets
  data
})
names(all_data) <- paste0("Group_", main_group_ids)

H5close()


#mat_data <- readMat("simulated_data.mat")

Index_M <- c("M1","M2","M3","M4","M5","M8");

SurvivalData <- list()
PhenoData <- list()
result_matrix <- matrix(nrow = N,ncol = 6)
index_f <- list()


for (i in 1:N) {
i
#data <- mat_data$datasets[[i]][[1]]
data <- all_data[[i]]
  
ID_M <- data[[6]]
PhenoData$Time <- data[[1]]
PhenoData$Treatment <- data[[2]]
PhenoData$Status <- as.integer(data[[5]])
PhenoData$Cov <- data[[3]]
M <- data[[4]]
colnames(M) <- paste0("M", 1:ncol(M))

SurvivalData$PhenoData <- as.data.frame(PhenoData)
SurvivalData$Mediator <- scale(M[,ID_M])


hima_survival.fit <- hima_survival(
  X = PhenoData$Treatment,
  M = SurvivalData$Mediator,
  OT = PhenoData$Time,
  status = PhenoData$Status,
  COV = PhenoData$Cov,
  scale = FALSE,
  FDRcut = 0.05,
  topN = 150,
  verbose = FALSE
)

#hima_survival.fit

result_matrix[i, ] <- as.integer(Index_M %in% hima_survival.fit$Index)

index_f[[i]] <- hima_survival.fit$Index


}

cleaned_list <- Filter(Negate(is.null), index_f)
FDR <- matrix(nrow = length(cleaned_list),ncol = 1)
power_individual <- matrix(nrow = length(cleaned_list),ncol = 6)
power_overall <- matrix(nrow = length(cleaned_list),ncol = 1)

for (i in 1:length(cleaned_list)) {
  
  
  TP <- sum(cleaned_list[[i]] %in% Index_M) 
  FP <- sum(!cleaned_list[[i]] %in% Index_M)  

  FDR[i,] <- FP / (TP + FP)


  # power_individual[i,] <- sapply(Index_M, function(var) {
  #   sum(sapply(cleaned_list[[i]], function(x) var %in% x)) / length(Index_M)
  # })

  power_overall[i,] <- mean(sapply(hima_survival.fit$Index, function(x) {
    sum(x %in% Index_M) / N
  }))
  
  power_individual[i,] <- sapply(Index_M, function(var) {
    sum(sapply(cleaned_list[[i]], function(x) var %in% x)) / N  # 分母为模拟次数
  })
  

  
}


fdr_M <- sum(FDR)/N

power_indi <- apply(power_individual, 2,sum)

power_all <- sum(power_overall)


#print(paste("FDR =", round(FDR, 3)))

# print(paste("power_i:", round(power_individual, 3)))
# print(paste("power_all:", round(power_overall, 3)))


# for(i in 1:1) {
# datasets <- mat_data$datasets[[i]]
# i=1
# PhenoData <- data.frame(
#   Treatment = datasets[[i]][[2]],
#   Time = datasets[[i]][[1]],
#   Status = datasets[[i]][[5]],
#   M = datasets[[i]][[4]],
#   COV = datasets[[i]][[3]]
# )
# 



