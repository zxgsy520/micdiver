#!/usr/bin/env Rscript

library(pROC)
library(dplyr)
library(tidyverse)
library(randomForest)
library(ggforce)
library(ggplot2)

options(bitmapType='cairo') #关闭服务器与界面的互动响应

add_help_args <- function(args){

    if(length(args) <= 1) {
        cat("Version: v1.0.0\n")
        cat("Author:Xingguo Zhang\n")
        cat("Email:invicoun@foxmail.com\n")
        cat("Function:Random Forest Analysis.\n")
        cat("Example:Rscript random_forest.R abundance_species.xls prefix\n")
        cat("Input file format:\n")
        cat("Tax Id\tsample1\tsample2\tsample3\n")
        cat("OTU1\t10\t9\t9")
        quit()
    }
}

run_plot_pca <- function(data, prefix){

    collen <- length(colnames(data))
    group <- data$group
 
    data_nogroup.pca <- prcomp(data[,1:collen-1])
    data_nogroup.pca <- data.frame(data_nogroup.pca$x)
    data_nogroup.pca <- data.frame(data_nogroup.pca, group)
    data_nogroup.pca <- merge(data_nogroup.pca, group, all.x=TRUE, by='row.names')

    mx <- max(data_nogroup.pca$PC1)*1.5#定义坐标轴长
    nx <- min(data_nogroup.pca$PC1)*1.5
    my <- max(data_nogroup.pca$PC2)*1.5
    ny <- min(data_nogroup.pca$PC2)*1.5
    p <- ggplot(data_nogroup.pca, aes(PC1, PC2)) + geom_point(aes(colour=group, shape=group), alpha=0.6, size=3) +
                labs(x="PCA1", y="PCA2") +  geom_vline(xintercept=0, color='gray', size=0.1) +
                geom_hline(yintercept=0, color='gray', size=0.1) + theme_bw() +
                geom_mark_ellipse(aes(fill=group, colour=group), alpha=0.2, tol=0.6, expand=unit(0.6, "mm")) +
                coord_cartesian(clip="off") +
                scale_x_continuous(limits=c(nx, mx)) + scale_y_continuous(limits=c(ny,my)) +
                theme(legend.text=element_text(size=rel(0.6)), legend.key.size=unit(0.3,'lines'), legend.title=element_blank()) + 
                scale_shape_manual(values = c(10,12,13,14,15,16,17,18,19,20,21,22,23,1:25))

    ggsave(paste(prefix, ".PCA.pdf",sep=""), p, width=120, height=100, units="mm")
    ggsave(paste(prefix, ".PCA.png",sep=""), p, width=120, height=100, units="mm")
    
}


run_randomForest <- function(data, prefix="txt", top=10){

    otu <- t(data)
    group <- rownames(otu)
    otu <- data.frame(otu, group)

    otu$group = factor(otu$group)	#把分组变量设为因子
    set.seed(315)
    grouprf <- randomForest(group ~ ., data=otu, ntree=500, importance=TRUE, proximity=TRUE)
    print(grouprf)
    #datarf <- round(importance(grouprf), 2)
    #print(datarf)
    #mean decrease accuracy表示随机森林预测准确性的降低程度，该值越大表示该变量的重要性越大。
    #mean decrease gini计算每个变量对分类树每个节点上观测值的异质性的影响，从而比较变量的重要性,该值越大表示该变量的重要性越大。
    datarf <- as_tibble(round(importance(grouprf), 2), rownames="Geneid") %>% arrange(desc(MeanDecreaseAccuracy))
    write.table(datarf, file=paste(prefix, ".random_forest.tsv", sep=""), sep="\t", row.names=FALSE, col.names=TRUE)
    
    toprf <- datarf[1:top,]
    topgene <- toprf$Geneid
    newdata <- data[which(rownames(data)%in%topgene),]
    
    newotu <- t(newdata)
    newotu <- data.frame(newotu, group)
    run_plot_pca(newotu, prefix)
    newotu$group <- factor(newotu$group)
    #print(topgene)
    #print(newotu)
    tgrouprf <- randomForest(group ~ ., data=newotu, ntree=500, importance=TRUE, proximity=TRUE)
    print(tgrouprf)
    
}

args <- commandArgs(T)

add_help_args(args)
data <- read.table(args[1], sep="\t", row.names=1, head=TRUE, check.names=FALSE, quote="")

top = 30
if(length(args) >= 3){
   top = args[3]
}
run_randomForest(data, args[2], top)
