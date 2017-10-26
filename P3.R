library(data.table)
library(stringr)
library(reshape2)
library(ggplot2)
#install.packages("kernlab")
library(kernlab)
library(RWeka)

formula #established from P1

#1. Linear
svm.linear = ksvm(as.matrix(all_dt[,3:1802]),all_dt$Classes,type="C-svc",kernel="vanilladot", cross=5)
#b(svm.linear) #negative intercept
# for scaled data
alpha.idxs <- alphaindex(svm.linear)[[1]]  # Indices of SVs in original data
alphas <- alpha(svm.linear)[[1]]
weight.vector <- (ymatrix(svm.linear)[alpha.idxs] * alphas) %*% xmatrix(svm.linear)[[1]]
b(svm.linear)
linear.weight = data.table(weight = t(weight.vector), gene = all_gene)
linear.weight[,abs_weight := abs(weight.V1)]
setorder(linear.weight, -abs_weight)
#check overlap with RF
merge(linear.weight[1:100,.(gene,weight.V1)],combine_rf_out[1:100,.(gene,avg_freq)], by = "gene")

#2. Homo Quadratic: offset = 0 
svm.poly0 = ksvm(as.matrix(all_dt[,3:1802]),all_dt$Classes,type="C-svc",kernel="polydot",cross=5,kpar = list(offset = 0, scale =1, degree = 2))
svb(svm.poly0) #negative intercept
alpha.idxs2 <- alphaindex(svm.poly0)[[1]]  # Indices of SVs in original data
alphas2 <- alpha(svm.poly0)[[1]]
weight.vector2 <- (ymatrix(svm.poly0)[alpha.idxs2] * alphas2) %*% xmatrix(svm.poly0)[[1]]
b(svm.poly0)
poly0.weight = data.table(weight = t(weight.vector2), gene = all_gene)
poly0.weight[,abs_weight := abs(weight.V1)]
setorder(poly0.weight, -abs_weight)
#check overlap with RF
merge(poly0.weight[1:100,.(gene, weight.V1)],combine_rf_out[1:100,.(gene,avg_freq)], by = "gene")

#Another way to select feature
#method below is credited to: 
#https://johanndejong.wordpress.com/2016/01/17/svm-with-recursive-feature-elimination/
svm_rfe <- function(X, y, elim_frac = 1 / ncol(X), ...) {
  # keep track of the iteration during which
  # a feature was eliminated
  ii <- rep(NA, ncol(X))
  i <- 0
  while ( any(is.na(ii)) ) {
    # indices of remaining features
    not_elim_yet <- which(is.na(ii))
    # number of features to eliminate
    n_to_elim <- ceiling ( elim_frac * length(not_elim_yet) )
    # train the classifier on the remaining features
    fit <- ksvm(X[,not_elim_yet], y, ...)
    # compute the primal problem coefficients from the dual
    # problem coefficients
    sv_i <- alphaindex(fit)[[1]]
    w <- t( coef(fit)[[1]] ) %*% X[ sv_i, not_elim_yet ]
    # eliminate the features with the smallest squared weights
    to_elim <- not_elim_yet[ head(order( w * w ), n_to_elim) ]
    ii[to_elim] <- i
    i <- i + 1
  }
  # convert iterations into ranks
  i - ii
}

poly0.test = svm_rfe(as.matrix(all_dt[,3:1802]),all_dt$Classes, elim_frac = 1700/1800, type="C-svc",kernel="polydot",cross=5,kpar = list(offset = 0, scale =1, degree = 2))
poly0.dt = data.table(gene=all_gene, rank = poly0.test)
merge(poly0.dt[rank<3, .(gene)],combine_rf_out[1:100,.(gene,avg_freq)], by = "gene")



#3. Inhomo Quadratic: offset = 1
svm.poly1 = ksvm(as.matrix(all_dt[,3:1802]),all_dt$Classes,type="C-svc",kernel="polydot",cross=5,kpar = list(offset = 1, scale =1, degree = 2))
svm.poly1
alpha(svm.poly1) #support vectors
alphaindex(svm.poly1)
nSV(svm.poly1)
b(svm.poly1) #negative intercept
alpha.idxs3 <- alphaindex(svm.poly1)[[1]]  # Indices of SVs in original data
alphas3 <- alpha(svm.poly1)[[1]]
y.sv3 <- as.numeric(all_dt$Classes[alpha.idxs3])
weight.vector3 <- (y.sv3 * alphas3) %*% xmatrix(svm.poly1)[[1]]
hist(weight.vector3)
poly1.weight = data.table(weight = t(weight.vector3), gene = all_gene)
poly1.weight[,abs_weight := abs(weight.V1)]
setorder(poly1.weight, -abs_weight)
#check overlap with RF
merge(poly1.weight[1:100,.(gene, weight.V1)],combine_rf_out[1:100,.(gene,avg_freq)], by = "gene")

poly1.test = svm_rfe(as.matrix(all_dt[,3:1802]),all_dt$Classes, elim_frac = 1700/1800, type="C-svc",kernel="polydot",cross=5,kpar = list(offset = 1, scale =1, degree = 2))
poly1.dt = data.table(gene=all_gene, rank = poly1.test)
merge(poly1.dt[rank<3, .(gene)],combine_rf_out[1:100,.(gene,avg_freq)], by = "gene")


linear.test = svm_rfe(as.matrix(all_dt[,3:1802]),all_dt$Classes, elim_frac = 1700/1800, type="C-svc",kernel="vanilladot",cross=5)
linear.dt = data.table(gene=all_gene, rank = linear.test)
merge(linear.dt[rank<3, .(gene)],combine_rf_out[1:100,.(gene,avg_freq)], by = "gene")
