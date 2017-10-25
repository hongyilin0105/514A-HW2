library(data.table)
library(stringr)
library(reshape2)
#install.packages("RWeka")
#install.packages("rJava")
library(RWeka)
library(rJava)
library(RWekajars)
library(partykit)
library(ggplot2)

#Load both ctrl and case data
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

#Construct formula to be used in models
equation = paste("Classes",paste(names(ctrl_subset_dt)[3:1802],collapse = "+"), sep = "~")
formula = as.formula(equation)
ctrl_subset_dt[,3:1802] = lapply(ctrl_subset_dt[,3:1802], as.numeric)
case_subset_dt[,3:1802] = lapply(case_subset_dt[,3:1802], as.numeric)
summary(ctrl_subset_dt[,3:10])

#impute missing data with column mean for ctrl data
impute.mean <- function(x) replace(x, is.na(x), mean(x, na.rm = TRUE))
ctrl_impute_dt = ctrl_subset_dt[, lapply(.SD, impute.mean)]
ctrl_impute_dt$Classes = as.factor(ctrl_impute_dt$Classes)
summary(ctrl_impute_dt[,1:5])
#impute missing data with column mean for case data
case_impute_dt = case_subset_dt[, lapply(.SD, impute.mean)]
case_impute_dt$Classes = as.factor(case_impute_dt$Classes)
summary(case_impute_dt[,1:5])

# combine case and ctrl data into all_dt
all_dt = rbindlist(list(ctrl_impute_dt, case_impute_dt), use.names = T, fill = T, idcol = F)
all_dt$Classes = as.factor(all_dt$Classes)

#1. fix # of features to 70%, change number of samples from 50% to 90%
all_gene = colnames(all_dt)[3:1802] #total list of gene names 

#Generate an RF with J48 trees based on # of trees, percentage of features to try
#at each split, and precentage of instances to use in each tree; rf_J48() returns
#a list of ntree strings with detailed structure of each tree in the RF
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

#Extract node/gene names using regex "GI_[0-9]+_." in each tree and deduplicate the names;
#Count the frequency of distinct genes in an RF and return a table of 2 columns
extract_rf_gene_freq <- function(rf, all_gene) {
  l = length(rf)
  gene_list = data.table()
  gene_count = data.table(gene = all_gene)
  for (i in 1:l) {
    temp = data.table(gene = unique(str_extract_all(rf[i], "GI_[0-9]+_.")[[1]]), i = 1)
    colnames(temp)[2]=i
    gene_count = merge(gene_count,temp, by = "gene", all.x=T)
  }
  gene_count[, count:=rowSums(gene_count[,2:l],na.rm = T)]
  gene_count[,.(gene, count)]
}

#an example call to generate and count frequency for a single RF
rf_50sampratio = rf_J48(all_dt, ntree = 10) #2:54-3:10
gene_freq = extract_rf_gene_freq(rf_50sampratio, all_gene)
gene_freq[order(count, decreasing = T),]

#iteratively build RF based on varying instances sample ratio;
#return a table of detailed gene counts from each RF for each gene
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

#Call on all_dt
outcome = rf_J48_iter_samp(ntree=500, mtry.ratio=0.7)#500tree: 8:09
write.csv(outcome, "rf1_outcome_1017.csv", row.names = F)


#2. fix # of samples to 70%, change # of features from 50% to 90%
#iteratively build RF based on varying percentage of features to drAW;
#return a table of detailed gene counts from each RF for each gene
rf_J48_iter_feat <- function(data=all_dt, ntree=ntree, samp.ratio=samp.ratio, start = 0.5, end = 0.9, by = 0.05) {
  feat_vector = seq(start, end, by)
  len = length(feat_vector)
  rf_out = mapply(rf_J48, feat_vector, MoreArgs = list(data=data, ntree=ntree, samp.ratio = samp.ratio), SIMPLIFY = T)
  rf_subset = as.list(data.frame(matrix(seq(1, length(rf_out), by=1), ncol=length(feat_vector))))
  detail_count = data.table(gene = all_gene)
  for (i in 1:len){
    detail_count = merge(detail_count, extract_rf_gene_freq(rf_out[rf_subset[[i]]], all_gene), by="gene", all.x = T, suffixes = c(as.character(i-1),as.character(i)))
  }
  detail_count
}

#Call on all_dt
outcome2=rf_J48_iter_feat(ntree=500, samp.ratio=0.7)#2：30
write.csv(outcome2, "rf2_outcome_1017.csv", row.names = F)

colnames(outcome)[2:10] = paste(as.character(seq(50, 90, by = 5)),"Sample", sep = "%")
colnames(outcome2)[2:10] = paste(as.character(seq(50, 90, by = 5)),"Feature", sep = "%")

#rank gene freq to see the 100 most significant gene
outcome$avg_freg_rf1 = rowMeans(outcome[,2:10])
outcome2$avg_freg_rf2 = rowMeans(outcome2[,2:10])
setorder(outcome, -avg_freg_rf1)
setorder(outcome2, -avg_freg_rf2)
common_gene = merge(outcome[1:100,.(gene, avg_freg_rf1)], outcome2[1:100,.(gene, avg_freg_rf2)], by = "gene")

combine_rf_out = merge(outcome, outcome2, by ="gene")
combine_rf_out[, avg_freq:= 0.5*(avg_freg_rf1+avg_freg_rf2)]
setorder(combine_rf_out, -avg_freq)
#Check the top 100 genes
rf_top100 = combine_rf_out[1:100,.(gene,avg_freq)]


outcome[, rank:=seq_len(.N) ]
outcome2[, rank:=seq_len(.N) ]

#Iteratively Add rank of each RF
setorder(outcome, -"50%Sample")
outcome[,rank1 := seq_len(.N)]
setorder(outcome, -"55%Sample")
outcome[,rank2 := seq_len(.N)]
setorder(outcome, -"60%Sample")
outcome[,rank3 := seq_len(.N)]
setorder(outcome, -"65%Sample")
outcome[,rank4 := seq_len(.N)]
setorder(outcome, -"70%Sample")
outcome[,rank5 := seq_len(.N)]
setorder(outcome, -"75%Sample")
outcome[,rank6 := seq_len(.N)]
setorder(outcome, -"80%Sample")
outcome[,rank7 := seq_len(.N)]
setorder(outcome, -"85%Sample")
outcome[,rank8 := seq_len(.N)]
setorder(outcome, -"90%Sample")
outcome[,rank9 := seq_len(.N)]
setorder(outcome, rank)

setorder(outcome2, -"50%Feature")
outcome2[,rank1 := seq_len(.N)]
setorder(outcome2, -"55%Feature")
outcome2[,rank2 := seq_len(.N)]
setorder(outcome2, -"60%Feature")
outcome2[,rank3 := seq_len(.N)]
setorder(outcome2, -"65%Feature")
outcome2[,rank4 := seq_len(.N)]
setorder(outcome2, -"70%Feature")
outcome2[,rank5 := seq_len(.N)]
setorder(outcome2, -"75%Feature")
outcome2[,rank6 := seq_len(.N)]
setorder(outcome2, -"80%Feature")
outcome2[,rank7 := seq_len(.N)]
setorder(outcome2, -"85%Feature")
outcome2[,rank8 := seq_len(.N)]
setorder(outcome2, -"90%Feature")
outcome2[,rank9 := seq_len(.N)]
setorder(outcome2, rank)

outcome_rank = outcome[,.(rank1,rank2,rank3,rank4,rank5,rank6,rank7,rank8,rank9)]
outcome_rank_sd = transform(outcome_rank, SD=apply(outcome_rank,1, sd, na.rm = TRUE))
avg_rank_sd = mean(outcome_rank_sd$SD)

outcome2_rank = outcome2[,.(rank1,rank2,rank3,rank4,rank5,rank6,rank7,rank8,rank9)]
outcome2_rank_sd = transform(outcome2_rank, SD=apply(outcome2_rank,1, sd, na.rm = TRUE))
avg_rank_sd2 = mean(outcome2_rank_sd$SD)



#visualize the gene rank deviation below
x <- seq(0, 650, length.out=1300)
df <- with(outcome_rank_sd, data.frame(x = x, y = dnorm(x, mean(SD), sd(SD))))
d1=ggplot(outcome_rank_sd, aes(x=SD, y = ..density..)) + geom_histogram(bins = 50, col = "black", fill = "grey") + 
  geom_line(data = df, aes(x = x, y = y), color = "red")+
  theme(plot.title = element_text(size = 14), axis.text.x  = element_text(size = 8))+
  labs(x="Standard Deviation of Rank", y="Density", title = "Gene Rank Deviation from RF set 1: Varying Instance Size")

df2 <- with(outcome2_rank_sd, data.frame(x = x, y = dnorm(x, mean(SD), sd(SD))))
d2 = ggplot(outcome2_rank_sd, aes(x=SD, y = ..density..)) + geom_histogram(bins = 50, col = "black", fill = "grey") + 
  geom_line(data = df2, aes(x = x, y = y), color = "red")+
  theme(plot.title = element_text(size = 14), axis.text.x  = element_text(size = 8))+
  labs(x="Standard Deviation of Rank", y="", title = "Gene Rank Deviation from RF set 2: Varying Feature Size")


multiplot(d1,d2, cols = 2)
#par(mfrow = c(1,2))
#hist(outcome_rank_sd$SD, breaks = 50, xlim = c(0,700), ylim = c(0,80))
#hist(outcome2_rank_sd$SD, breaks = 50, xlim = c(0,700), ylim = c(0,80))

#visualize the distribution of frequency
outcome_long = reshape(outcome, direction = "long", varying=list(names(outcome)[2:11]), v.names = "Freq", idvar = "gene", timevar = "Forest", times = names(outcome)[2:11])
p1 = ggplot(outcome_long[Forest != "avg_freg_rf1"], aes(x = Freq, fill = Forest)) +
  geom_histogram(stat = "count") +
  facet_grid(Forest~.)+
  labs(y = "Rank of Gene", x = "Gene", title = "Gene Importance Rank from 9 RFs")+
  theme(strip.text.y = element_blank(), legend.title = element_blank(), legend.text = element_text(size = 8))

outcome2_long = reshape(outcome2, direction = "long", varying=list(names(outcome2)[2:11]), v.names = "Freq", idvar = "gene", timevar = "Forest", times = names(outcome2)[2:11])
p2= ggplot(outcome2_long[Forest != "avg_freg_rf2"], aes(x = Freq, fill = Forest)) +
  geom_histogram(stat = "count") +
  facet_grid(Forest~.)+
  labs(y = "Count of Genes", x = "Frequency in each Random Forest", title = "Gene Frequency Distribution in RF w/ Varying Features")+
  theme(strip.text.y = element_blank(), legend.title = element_blank(), legend.text = element_text(size = 8))

multiplot(p1,p2, cols = 2)
#to call multiplot, requires the function below
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}