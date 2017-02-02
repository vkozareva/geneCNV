library(nnet)
library(ggplot2)
library(tidyr)
library(dplyr)
setwd("/Users/nigel/git/dmd/exon_data")
d = read.csv("coverage_matrix.csv", skip = 1)
head(d)

exColumns = grep("Ex",colnames(d), value = TRUE)
exonName = sapply(exColumns, function(x) as.numeric(strsplit(x, "Ex")[[1]][2]))
dn =d 
dn[,exColumns] = apply(dn[,exColumns], 1, function(x) x / sum(x))

totReads = apply(d[,exColumns], 1, sum)
d$Total = totReads
ggplot(d, aes(x=gender, fill=gender, y=Total)) + geom_violin()

cnts = d %>% gather("Exon", "Count", Ex1:Ex79)
cntsnorm = dn %>% gather("Exon", "Count", Ex1:Ex79)
cntsnorm$Exon = factor(cntsnorm$Exon, levels = exColumns)
cntsnorm = cntsnorm[totReads > 8000 & dn$gender=="F",]
ggplot(cnts, aes(x=Exon, y=Count, group=sample)) + geom_line() + theme_bw()
ggplot(cnts, aes(x=Exon, y=Count, group=sample)) + geom_violin() + theme_bw()

ggplot(cntsnorm, aes(x=Exon, y=Count, group=sample)) + geom_line() + theme_bw()
ggplot(cntsnorm, aes(x=Exon, y=Count)) + geom_violin(fill="blue") + 
  theme_bw() + theme(axis.text.x=element_text(angle = -45, hjust = 0)) + 
  labs(title = "Exon Relative Frequencies for Females with > 8000 Reads Total")

fd = d[d$gender == "F",]
sl = factor(fd$sample:fd$lane)
simple = multinom(as.matrix(fd[,exColumns]) ~ 1, summ=1, MaxNWts = 10000)
sample = multinom(as.matrix(fd[,exColumns]) ~ fd$sample, summ=1, MaxNWts = 10000)
flowcell = multinom(as.matrix(fd[,exColumns]) ~ fd$flow_cell_id, summ=1, MaxNWts = 10000)
sample_lane = multinom(as.matrix(fd[,exColumns]) ~ sl, summ=1, MaxNWts = 10000)
anova(simple, flowcell, sample, sample_lane)


fd %>% group_by(flow_cell_id) %>% summarise(sum(Ex1:Ex79))

nd = d %>% gather("Exon", "Count", Ex1:Ex79) %>% group_by(flow_cell_id, Exon, gender, subject) %>% summarise(Total = mean(Count)) 
head(nd)
spr = nd %>% spread(flow_cell_id, Total)
head(spr)
spr[,3:5] = apply(spr[,3:5], 2, function(z) z / sum(z))
ggplot(spr, aes(x=male, y = female)) + geom_point() + geom_abline(slope = 1, intercept = 0) + facet_wrap(~ flow)

spr = nd %>% spread(gender, Total)
head(spr)
ggplot(spr, aes(x=M, y = F)) + geom_point() + geom_abline(slope = 2, intercept = 0) + facet_wrap(~ flow_cell_id)
ggplot(spr, aes(x=M, y = F)) + geom_point() + geom_abline(slope = 2, intercept = 0) + geom_abline(color="red", slope = 1, intercept = 0) + facet_wrap(~ flow_cell_id)


simplemf = multinom(as.matrix(d[,exColumns]) ~ 1, summ=1, MaxNWts = 10000)