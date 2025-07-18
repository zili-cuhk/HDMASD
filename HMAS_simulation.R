#High-dimensional Mediation Analysis in Survival Models
# load R packages
library(survival)
library(ncvreg)
library(ggm)

# SIS for alpha
sis_alpha <- function(X, M, COV, p){
  s_alpha <- function(j){
    if (is.null(COV)) {
      MX <- data.frame(M = M[, j], X = X)
    } else {
      MX <- data.frame(M = M[, j], X = X, COV = COV)
    }
    fit <- glm(M ~., data = MX)
    s1 <- summary(fit)$cov.scaled[2,2]   #var for alpha
    s2 <- summary(fit)$coef[2]           #coefficients for alpha
    s3 <- summary(fit)$coef[2,4]         #p-value for alpha
    return(data.frame(s_var=s1, s_coef=s2, s_p=s3))
  }
  dat=data.frame(do.call(rbind, lapply(1:ncol(M), s_alpha)))
  alpha_sis <- t(dat)
  colnames(alpha_sis) = colnames(M)
  return(s_alpha=alpha_sis)
}

# SIS for beta (the regression Y~X+M)
sis_beta1 <- function(X, M, Y, COV, p){
  s_beta <- function(j){
    if (is.null(COV)) {
      MX <- data.frame(Y=Y, M = M[, j], X = X)
    } else {
      MX <- data.frame(Y=Y, M = M[, j], X = X, COV = COV)
    }
    fit <- coxph(Y ~., data = MX)
    s1 <- fit$var[1,1]                   #var for alpha
    s2 <- summary(fit)$coef[1]           #coefficients for alpha
    s3 <- summary(fit)$coef[1,5]         #p-value for alpha
    return(data.frame(s_var=s1, s_coef=s2, s_p=s3))
  }
  dat=data.frame(do.call(rbind, lapply(1:ncol(M), s_beta)))
  beta_sis <- t(dat)
  colnames(beta_sis) = colnames(M)
  return(s_beta=beta_sis)
}

# SIS for beta (pcor for M and Y)
sis_beta2 <- function(X,M,Y,COV,p){
  p_cor = matrix(0,p,1)
  for(j in 1:p){
    if(is.null(COV)){
      d=data.frame(Y=Y[,1], M=M[,j], X=X)
      p_cor[j,] <- pcor(c(1,2,3),cov(d, method='spearman'))  #method of correlation
    }else{
      d=data.frame(Y=Y[,1], M=M[,j], X=X, COV=COV)            
      p_cor[j,] <- pcor(c(1,2,3,4,5), cov(d, method='spearman'))  #with 2 covariates
    }
  }
  rownames(p_cor) = colnames(M)
  return(p_cor = p_cor)
}

#main function
hmas <- function(X, Y, M, COV, k,
                 penalty = c("MCP", "SCAD", "lasso"),
                 path = c('MY', 'MX'),
                 topN = NULL,
                 verbose = TRUE, 
                 ...) {
  penalty <- match.arg(penalty)
  
  n <- nrow(M)
  p <- ncol(M)

  
  if (is.null(topN)) {
    if (path == 'MY'){
      d <- ceiling(k*n/log(n))  #the top d mediators that associated with outcome, k=2
    }else{
      d <- ceiling(k*n/log(n))  #the top d mediators that associated with exposure
    }
  } else {
    d <- topN  
  }
  
  d <- min(p, d)  # revised on 24/04/2025 if d_0 > p select all mediators
  
  if(verbose) message("Step 1: Prelimenary Screening...", "     (", Sys.time(), ")")
  
  if (path == 'MY'){
    # cor <- sis_beta2(X=X, M=M, Y=Y, COV=COV, p=ncol(M))
    # cor <- abs(cor)
    # rownames(cor)=paste0('M', 1:ncol(M))
    # cor_sort <- sort(cor, decreasing=T)
    beta_s <- sis_beta1(X=X, M=M, Y=Y, COV=COV, p=ncol(M))
    SIS_beta <- beta_s[3,]
    SIS_beta_sort <- sort(SIS_beta)
    ID_SIS <- which(SIS_beta <= SIS_beta_sort[d])  # the index of top d significant mediators (Y~X+M)
    # ID_SIS <- which(cor >= cor_sort[d])    
  }else{
    alpha_s <- sis_alpha(X, M, COV, p)
    SIS_alpha1 <- alpha_s[3,]
    SIS_alpha_sort <- sort(SIS_alpha1)
    ID_SIS <- which(SIS_alpha1 <= SIS_alpha_sort[d])  # the index of top d significant mediators (M~X)
  }
  
  M_SIS <- M[, ID_SIS]
  
  if(verbose) cat("Top", length(ID_SIS), "mediators selected (ranked from most to least significant): ", colnames(M_SIS), "\n")
  
  XM <- cbind(M_SIS, X)
  C_M <- colnames(XM)
  
  
  if(verbose) message("Step 2: Penalized Variable Selection (", penalty, ") ...", "   s  (", 
                      Sys.time(), ")")
  
  if (is.null(COV)) {
    fit <- ncvsurv(XM, Y,
                   penalty = penalty,
                   penalty.factor = c(rep(1, ncol(M_SIS)), 0), ...)
  } else {
    COV <- data.frame(COV)
    COV <- data.frame(model.matrix(~., COV))[, -1,drop=F]
    conf.names <- colnames(COV)
    XM_COV <- cbind(XM, COV)
    fit <- ncvsurv(XM_COV, Y,
                   penalty = penalty,
                   penalty.factor = c(rep(1, ncol(M_SIS)), rep(0, 1 + ncol(COV))), ...)
  }
  
  
  lam <- fit$lambda[which.min(BIC(fit))]
  if(verbose) cat("lambda selected: ", lam, "\n")
  Coefficients <- coef(fit, lambda = lam)
  est <- Coefficients[1:length(ID_SIS)]
  ID_p_non <- which(est != 0)
  
  if(verbose) cat("Non-zero", penalty, "beta estimate(s) of mediator(s) found: ", names(ID_p_non), "\n")
  
  beta_p <- est[ID_p_non]  # The non-zero MCP estimators of beta
  ID_p <- ID_SIS[ID_p_non]  # The index of the ID of non-zero beta in the Cox regression
  MCP_M <- names(ID_p_non)
  
  
  if(verbose) message("Step 3: The adjusted Sobel significance test ...", 
                      "     (", Sys.time(), ")")
  
  if (path == 'MY'){
    sub_M <- M[, MCP_M]
    sis_alpha1 <- sis_alpha(X, sub_M, COV, ncol(sub_M))
    alpha_s <- sis_alpha1[, MCP_M]
  }else{
    alpha_s <- alpha_s[, ID_p]
  }
  
  alpha_est <- alpha_s[2, ]   #  the estimator for alpha
  var_alpha <- alpha_s[1, ]
  
  # # true alpha and beta
  # beta_t <- beta[ID_p]
  # alpha_t <- alpha[ID_p]
  # ab_true <- alpha_t * beta_t
  
  if (is.null(COV)) {
    YMX <- data.frame(Y = Y, M[, ID_p, drop = FALSE], X = X)
  } else {
    YMX <- data.frame(Y = Y, M[, ID_p, drop = FALSE], X = X, COV = COV)
  }
  
  cox_model <- coxph(Y ~ ., data = YMX)
  
  beta_est <- summary(cox_model)$coefficients[1: length(ID_p)]     #the estimator of beta
  DE <- summary(cox_model)$coefficients[(length(ID_p)+1), 1]
  DE <- exp(DE)
  
  var_beta <- diag(cox_model$var[1:length(ID_p),1:length(ID_p)])
  
  ab_est <- alpha_est * beta_est   # the estimator of alpha*beta
  
  # var(alpha*beta)
  var_ab <- (alpha_est^2) * var_beta + var_alpha * (beta_est^2)
  
  # confidence interval
  conf_low <- ab_est - 1.96 * sqrt(var_ab)
  conf_up <- ab_est + 1.96 * sqrt(var_ab)
  
  # sobel test for alpha*beta
  s.test <- abs(ab_est)/sqrt(var_ab)   #z-score for sobel test
  P_sobel <- 2 * (1-pnorm(s.test))     #p-value of sobel test
  P_bon_sobel <- p.adjust(P_sobel, 'bonferroni', length(ID_p)) #Bonferroni adjusted p-value
  P_bon_sobel[P_bon_sobel > 1] <- 1
  P_fdr_sobel <- p.adjust(P_sobel, 'fdr', length(ID_p))
  P_fdr_sobel[P_fdr_sobel > 1] <- 1
  
  #joint test for alpha*beta
  P_bon_alpha <- length(ID_p) * sis_alpha1[3, MCP_M]  # The adjusted p-value for alpha (bonferroni)
  P_bon_alpha[P_bon_alpha > 1] <- 1
  P_fdr_alpha <- p.adjust(sis_alpha1[3, MCP_M], 'fdr', length(ID_p))  # The adjusted p-value for alpha (bonferroni)
  P_fdr_alpha[P_fdr_alpha > 1] <- 1
  P_bon_beta <- length(ID_p) * summary(cox_model)$coefficients[1:length(ID_p), 5]  # The adjused p-value for beta
  P_bon_beta[P_bon_beta > 1] <- 1
  P_fdr_beta <- p.adjust(summary(cox_model)$coefficients[1:length(ID_p), 5], 'fdr', length(ID_p))  # The adjusted p-value for alpha (bonferroni)
  P_fdr_beta[P_fdr_beta > 1] <- 1
  PB <- rbind(P_bon_beta, P_bon_alpha)
  P_bon_joint <- apply(PB, 2, max)
  PF <- rbind(P_fdr_beta, P_fdr_alpha)
  P_fdr_joint <- apply(PF, 2, max)
  
  # results <- data.frame(alpha = alpha_est, beta = beta_est,
  #                       `alpha_est*beta_est` = ab_est, `alpha_t*beta_t` = ab_true, 
  #                       conf_low=conf_low, conf_up=conf_up,
  #                       P_bon_sobel=P_bon_sobel, P_fdr_sobel=P_fdr_sobel,
  #                       P_bon_joint=P_bon_joint, P_fdr_joint=P_fdr_joint,
  #                       var_ab=var_ab, var_alpha=var_alpha, var_beta=var_beta,
  #                       DE, check.names = FALSE)
  
  
  results <- data.frame(alpha = alpha_est, beta = beta_est,
                        `alpha_est*beta_est` = ab_est,
                        conf_low=conf_low, conf_up=conf_up,
                        P_bon_sobel=P_bon_sobel, P_fdr_sobel=P_fdr_sobel,
                        P_bon_joint=P_bon_joint, P_fdr_joint=P_fdr_joint,
                        var_ab=var_ab, var_alpha=var_alpha, var_beta=var_beta,
                        DE, check.names = FALSE)
  
  
  if(verbose) message("Done!", "     (", Sys.time(), ")")
  
  return(list(C_M, MCP_M, results))
}



##########################Real-data-analysis###############

########################simulation studies#######################
library(Matrix)
library(glmnet)
library(ncvreg)
library(R.matlab)
library(future.apply)
library(rhdf5)

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



mat_data <- readMat("simulated_data7.mat")

Index_M <- c("M1","M2","M3","M4","M5","M8");

SurvivalData <- list()
PhenoData <- list()
result_matrix <- matrix(nrow = N,ncol = 6)
index_f <- list()


for (i in 1:N) {
  print(paste("当前循环次数：", i))
  #data <- mat_data$datasets[[i]][[1]]
  data <- all_data[[i]]
  
  ID_M <- mat_data[["datasets"]][[i]][[1]][[1]]
  PhenoData$Time <- data[[1]]
  PhenoData$Treatment <- data[[2]]
  PhenoData$Status <- as.integer(data[[5]])
  PhenoData$Cov <- data[[3]]
  M <- data[[4]]
  colnames(M) <- paste0("M", 1:ncol(M))
  
  SurvivalData$PhenoData <- as.data.frame(PhenoData)
  SurvivalData$Mediator <- scale(M[,ID_M])



hima_survival.fit <- hmas(
  X = SurvivalData$PhenoData$Treatment,
  #Y= survival::Surv(PhenoData$Time, PhenoData$Status),
  Y= survival::Surv(SurvivalData$PhenoData$Time, SurvivalData$PhenoData$Status),
  M = SurvivalData$Mediator,
  COV = SurvivalData$PhenoData[, c("Cov.1", "Cov.1")],
  k=1,
  penalty = "MCP", # OR "SCAD","lasso"
  path = 'MY', # OR 'MY'
  topN = 150,
  verbose = TRUE, 
) 

result_matrix[i, ] <- as.integer(Index_M %in% hima_survival.fit[[2]])
index_f[[i]] <- hima_survival.fit[[2]]

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
  # 
  power_overall[i,] <- mean(sapply(hima_survival.fit[[2]], function(x) {
    sum(x %in% Index_M) / N
  }))
  
  power_individual[i,] <- sapply(Index_M, function(var) {
    sum(sapply(cleaned_list[[i]], function(x) var %in% x)) / N  # 分母为模拟次数
  })
  
}

fdr_M <- sum(FDR)/N
power_indi <- apply(power_individual, 2,sum)
power_all <- sum(power_overall)

