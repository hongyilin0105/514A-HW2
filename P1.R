library(data.table)
library(randomForest)
library(stringr)
library(reshape2)
#install.packages("RWeka")
#install.packages("rJava")
library(RWeka)
library(rJava)
library(RWekajars)
library(partykit)


raw_ctrl = read.table("ctrl.gex", sep = "")
raw_ctrl_dt = data.table(t(raw_ctrl))
colnames(raw_ctrl_dt) = as.character(raw_ctrl_dt[1,])
ctrl_subset_dt = raw_ctrl_dt[-1, 1:1802]
colnames(ctrl_subset_dt) = str_replace(colnames(ctrl_subset_dt),"-","_")
ctrl_subset_dt$Classes = as.factor(ctrl_subset_dt$Classes)

raw_case = read.table("case.gex", sep = "")
raw_case_dt = data.table(t(raw_case))
colnames(raw_case_dt) = as.character(raw_case_dt[1,])
case_subset_dt = raw_case_dt[-1, 1:1802]
colnames(case_subset_dt) = str_replace(colnames(case_subset_dt),"-","_")
case_subset_dt$Classes = as.factor(case_subset_dt$Classes)

equation = paste("Classes",paste(names(ctrl_subset_dt)[3:1802],collapse = "+"), sep = "~")
formula = as.formula(equation)
ctrl_subset_dt[,3:1802] = lapply(ctrl_subset_dt[,3:1802], as.numeric)
case_subset_dt[,3:1802] = lapply(case_subset_dt[,3:1802], as.numeric)
summary(ctrl_subset_dt[,3:10])

#impute missing data with column mean
impute.mean <- function(x) replace(x, is.na(x), mean(x, na.rm = TRUE))

ctrl_impute_dt = ctrl_subset_dt[, lapply(.SD, impute.mean)]
ctrl_impute_dt$Classes = as.factor(ctrl_impute_dt$Classes)
summary(ctrl_impute_dt[,1:5])

case_impute_dt = case_subset_dt[, lapply(.SD, impute.mean)]
case_impute_dt$Classes = as.factor(case_impute_dt$Classes)
summary(case_impute_dt[,1:5])

# combine case and ctrl
all_dt = rbindlist(list(ctrl_impute_dt, case_impute_dt), use.names = T, fill = T, idcol = F)
all_dt$Classes = as.factor(all_dt$Classes)

#1. fix # of features to 70%, change number of samples from 50% to 90%
all_gene = colnames(all_dt)[3:1802]

rf_J48 <- function(data, ntree = 500, mtry.ratio = 0.7, samp.ratio = 0.5) {
  samp_size = nrow(data)
  out_tree = list()
  for (n in 1:ntree){
    subset= sample(samp_size, samp_size*samp.ratio)
    all_f = colnames(data)[3:ncol(data)]
    select_f = all_f[sample(length(all_f), length(all_f)*mtry.ratio)]
    equation = paste("Classes",paste(select_f,collapse = "+"), sep = "~")
    formula = as.formula(equation)
    t = J48(formula,data, subset)
    out_tree[[n]] = capture.output(t)
  }
  out_tree
}

extract_rf_gene_freq <- function(rf, all_gene) {
  l = length(rf)
  gene_list = data.table()
  gene_count = data.table(gene = all_gene)
  for (i in 1:l) {
    temp = data.table(gene = unique(str_extract_all(rf[i], "GI_[0-9]+_.")[[1]]), i = 1)
    colnames(temp)[2]=i
    gene_count = merge(gene_count,temp, by = "gene", all.x=T)
  }
  gene_count[, count:=colSums(gene_count[,2:l],na.rm = T)]
  gene_count[,.(gene, count)]
}

#an example call
rf_50sampratio = rf_J48(all_dt) #2:54-3:10
gene_freq = extract_rf_gene_freq(rf_50sampratio, all_gene)
gene_freq[order(count, decreasing = T),]

# iterative sample ratio
rf_J48_iter_samp <- function(data=all_dt, ntree=ntree, mtry.ratio=mtry.ratio, start = 0.5, end = 0.9, by = 0.05) {
  samp_vector = seq(start, end, by)
  len = length(samp_vector)
  rf_out = mapply(rf_J48, samp_vector, MoreArgs = list(data=data, ntree=ntree, mtry.ratio = mtry.ratio), SIMPLIFY = T)
  rf_subset = as.list(data.frame(matrix(seq(1, length(rf_out), by=1), ncol=length(samp_vector))))
  detail_count = data.table(gene = all_gene)
  for (i in 1:len){
    detail_count = merge(detail_count, extract_rf_gene_freq(rf_out[rf_subset[[i]]], all_gene), by="gene", all.x = T, suffixes = c(as.character(i-1),as.character(i)))
  }
  detail_count
}

outcome = rf_J48_iter_samp(ntree=500, mtry.ratio=0.7) #500tree:6:23

#2. fix # of samples to 70%, change # of features from 50% to 90%
rf_J48_iter_feat <- function(data=all_dt, ntree=ntree, samp.ratio=samp.ratio, start = 0.5, end = 0.9, by = 0.05) {
  feat_vector = seq(start, end, by)
  len = length(feat_vector)
  rf_out = mapply(rf_J48, feat_vector, MoreArgs = list(data=data, ntree=ntree, samp.ratio = samp.ratio), SIMPLIFY = T)
  rf_subset = as.list(data.frame(matrix(seq(1, length(rf_out), by=1), ncol=length(samp_vector))))
  detail_count = data.table(gene = all_gene)
  for (i in 1:len){
    detail_count = merge(detail_count, extract_rf_gene_freq(rf_out[rf_subset[[i]]], all_gene), by="gene", all.x = T, suffixes = c(as.character(i-1),as.character(i)))
  }
  detail_count
}

outcome2 = rf_J48_iter_feat(ntree=500, samp.ratio=0.7) #500tree:9:21

#rank gene freq to see the 100 most significant gene
outcome$avg_freg_rf1 = rowMeans(outcome[,2:10])
outcome2$avg_freg_rf2 = rowMeans(outcome2[,2:10])
setorder(outcome, -avg_freg_rf1)
setorder(outcome2, -avg_freg_rf2)
common_gene = merge(outcome[1:100,.(gene, avg_freg_rf1)], outcome2[1:100,.(gene, avg_freg_rf2)], by = "gene")

combine_rf_out = merge(outcome, outcome2, by ="gene")
combine_rf_out[, avg_freq:= 0.5*(avg_freg_rf1+avg_freg_rf2)]
setorder(combine_rf_out, -avg_freq)
