# ----0、条件设置----
# 系统报错改为英文
Sys.setenv(LANGUAGE = "en")
# 禁止转化为因子
options(stringsAsFactors = FALSE)
# 清空环境
rm(list=ls())
# 设置工作目录
setwd("/home/data/t180324/R/R_project/转录组学分析/GSE162955/")

library(dplyr)
library(Seurat)
library(patchwork)
library(stringr)
library(GEOquery)
library(data.table)
library(tidyverse)
library(oligo) 
library(pd.hg.u133.plus.2)
library(isotree)
library(missForest)

# ----1、获取生存信息与临床信息----
proj = "GSE162955"

if(!file.exists(paste0(proj,"_series_matrix.txt.gz"))){
  eSet <- getGEO(proj,destdir = ".",getGPL = T)
  pd <- eSet[[paste0(proj,"_series_matrix.txt.gz")]]@phenoData@data
}else{
  eSet <- getGEO(filename=paste0("./",proj,"_series_matrix.txt.gz"), getGPL = T,destdir = ".")
  pd <- pData(eSet)
}

exp <- as.data.frame(exprs(eSet)) #表达矩阵
gpl <- getGEO(filename = './GPL17586.soft.gz')
gpl=Table(gpl)

# ----2、数据预处理----
gpl$symbol <- gsub("---",NA,str_split(gpl$gene_assignment," // ",simplify = T)[,2])
exp$ID <- rownames(exp)
exp <- na.omit(merge(exp,gpl[,c("ID","symbol")],by="ID"))[,-1]

exp$symbol[duplicated(exp$symbol)] #查看依然重复的基因名
#计算行平均值，按降序排列
index=order(rowMeans(exp[,!grepl("symbol",colnames(exp))]),decreasing = T)
#调整表达谱的基因顺序
exp_ordered=exp[index,]
#对于有重复的基因，保留第一次出现的那个，即行平均值大的那个
keep=!duplicated(exp_ordered$symbol)
#得到最后处理之后的表达谱矩阵
exp_max=exp_ordered[keep,]
exp_max

rownames(exp_max) <- exp_max[,grepl("symbol",colnames(exp))] #将基因名作为行名
exp_1 <- exp_max[,!grepl("symbol",colnames(exp))]
expr <- t(as.data.frame(apply(exp_1,1,as.numeric)))
colnames(expr) <- unlist(colnames(exp_1))

# 3、----分组信息获取----
pd <- pd[-c(1:2),]
table(pd$`tissue:ch1`)
Ctr=pd$geo_accession[grepl('Healthy',as.character(pd$`tissue:ch1`))] #对照组
IS=pd$geo_accession[grepl('Infarcted',as.character(pd$`tissue:ch1`))] #IS组
group=factor(c(rep('Ctr',length(Ctr)) ,
               rep('IS',length(IS))),levels=c('Ctr','IS'))
table(group)

save(expr,IS,Ctr,group,file = paste0("./",proj,".Rdata")) # 保存数据
load("./GSE162955.Rdata")

#---------三大R包差异分析----------
#----deseq2----
exp <- 2^expr-0.01

library(DESeq2)
colData <- data.frame(row.names=c(Ctr,IS),
                      condition=group)
colData_1 <- as.data.frame(colData[order(rownames(colData)),])
rownames(colData_1) <- sort(rownames(colData))
colnames(colData_1) <- colnames(colData)
colData <- colData_1

dds <- DESeqDataSetFromMatrix(
  countData = round(exp),
  colData = colData,
  design = ~ condition)
dds <- DESeq(dds)
class(dds)
res <- results(dds, contrast = c("condition",rev(levels(group))))

#constrast
c("condition",rev(levels(group)))
class(res)
DEG1 <- as.data.frame(res)
DEG1 <- DEG1[order(DEG1$pvalue),] 
DEG1 = na.omit(DEG1)
head(DEG1)

#添加change列标记基因上调下调
logFC_t = 0
pvalue_t = 0.05

k1 = (DEG1$pvalue < pvalue_t)&(DEG1$log2FoldChange < -logFC_t);table(k1)
k2 = (DEG1$pvalue < pvalue_t)&(DEG1$log2FoldChange > logFC_t);table(k2)
DEG1$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
table(DEG1$change)
head(DEG1)

#----edgeR----
library(edgeR)
dge <- DGEList(counts=exp,group=group)
dge$samples$lib.size <- colSums(dge$counts)
dge <- calcNormFactors(dge) 

design <- model.matrix(~group)

dge <- estimateGLMCommonDisp(dge, design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)

fit <- glmFit(dge, design)
fit <- glmLRT(fit) 

DEG2=topTags(fit, n=Inf)
class(DEG2)
DEG2=as.data.frame(DEG2)
head(DEG2)

k1 = (DEG2$PValue < pvalue_t)&(DEG2$logFC < -logFC_t);table(k1)
k2 = (DEG2$PValue < pvalue_t)&(DEG2$logFC > logFC_t);table(k2)
DEG2$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))

head(DEG2)
table(DEG2$change)
#----limma----
library(limma)
dge <- edgeR::DGEList(counts=exp)
dge <- edgeR::calcNormFactors(dge)
design <- model.matrix(~group)
v <- voom(dge,design, normalize="quantile")

design <- model.matrix(~group)
fit <- lmFit(v, design)
fit= eBayes(fit)

DEG3 = topTable(fit, coef=2, n=Inf)
DEG3 = na.omit(DEG3)

k1 = (DEG3$P.Value < pvalue_t)&(DEG3$logFC < -logFC_t);table(k1)
k2 = (DEG3$P.Value < pvalue_t)&(DEG3$logFC > logFC_t);table(k2)
DEG3$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
table(DEG3$change)
head(DEG3)

tj = data.frame(deseq2 = as.integer(table(DEG1$change)),
                edgeR = as.integer(table(DEG2$change)),
                limma_voom = as.integer(table(DEG3$change)),
                row.names = c("down","not","up")
);tj

# 整合分组进矩阵
expr_batch <- as.data.frame(t(expr[,!duplicated(colnames(expr))]))
rownames(expr_batch) <- gsub("-","",rownames(expr_batch))
expr_batch$result <- c(NA)

for (i in 1:nrow(expr_batch)){
  for (j in 1:length(levels(group))){
    if(rownames(expr_batch)[i]%in%get(levels(group)[j])){expr_batch$result[i] <- levels(group)[j]}
  }
}
expr_batch <- expr_batch[!is.na(expr_batch$result),]
# 根据分组和数据集进行预处理
exp <- data.frame(t(2^expr_batch[,!grepl("result",colnames(expr_batch))]))
expr <- data.frame(expr_batch)

library(tidyverse)
library(rstatix)
library(ggpubr)
library(ggtext)
color <- c("#B8D4A4","#C485CD")
# 调整矩阵格式
key_DEG <- c("VCAM1")
plot.info_all <- expr_batch[,c(key_DEG,"result")] %>%
  as.data.frame() %>%
  rownames_to_column("expr_batch") %>%
  pivot_longer(cols = 1:length(key_DEG)+1,
               names_to = "VAR",
               values_to = "VAL")
# 计算p值
stat.test <- plot.info_all %>% 
  group_by(VAR) %>%
  t_test(VAL ~ result) %>%
  adjust_pvalue(method = "none") %>% 
  add_significance("p")
stat.test$p <- DEG2$PValue[rownames(DEG2)%in%key_DEG] 
stat.test <- add_significance(stat.test,"p")
stat.test <- stat.test %>% 
  add_xy_position(x='VAR',dodge = 1)

# 画图
ggboxplot(
  plot.info_all[order(plot.info_all$VAR),],
  x = "VAR",
  y = "VAL",
  color = "result", 
  xlab = "",
  ylab = "",
  add = "jitter"
) +
  stat_pvalue_manual(stat.test,  label = "{p.signif}", tip.length = 0, hide.ns = F)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

#----火山图展示----
#需要突出显示的基因列表
geneList0 <- c(key_DEG)
geneList <- DEG2[geneList0,]
library('ggplot2')
p <- ggplot(# 数据、映射、颜色
  DEG2, aes(x = logFC, y = -log10(PValue), colour=change)) +
  geom_point(alpha=0.5, size=3.5) +
  scale_color_manual(values=c("#546de5","#d2dae2","#ff4757"))+
  #突出表示差异基因
  geom_point(data=geneList,aes(x = logFC, y = -log10(PValue)),colour="yellow",size=3.5)+
  #辅助线
  geom_vline(xintercept=0,lty=3,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.8) +
  labs(x="logFC",y="-log10 (PValue)")+   # 坐标轴# 坐标轴和图标题title="Volcano plot",
  theme_bw()+    #去除背景色
  theme(panel.grid = element_blank())+  #去除网格线
  #xlim(-2, 2)+   #设置坐标轴范围
  #图例
  theme(plot.title = element_text(hjust = 0.5,size=24), 
        legend.position="bottom", 
        legend.title = element_blank(),
        legend.text=element_text(size=18),
        legend.key.size = unit(1, 'cm'),
        legend.background = element_rect(fill="gray90", linetype="solid",colour ="gray"),
        axis.title.x =element_text(size=18), 
        axis.title.y=element_text(size=18),
        axis.text=element_text(size=14,face = "bold"))
p
# 再作标注基因名称的图
DEG_vol <- DEG2
DEG_vol$gene <- rownames(DEG_vol)
rownames(DEG_vol) <- 1:nrow(DEG_vol)

#标记出5个基因的label
geneList <- as.data.frame(geneList0)
geneList[,2] <- geneList
colnames(geneList) <- c('gene','label')

c <-merge(DEG_vol,geneList,by='gene',all.x=TRUE)  #增加label列，以突出显示指定基因

library(ggrepel)
p + geom_label_repel(data = c, 
                     aes(x = logFC, y = -log10(PValue), label = label),
                     size = 4,color="black",
                     #box.padding = unit(0.5, "lines"),
                     #point.padding = unit(0.8, "lines"), 
                     #segment.color = "black",   #连线的颜色
                     #segment.size = 0.4,  #连线的粗细
                     #arrow = arrow(length=unit(0.01, "npc")), #标签、点之间连线的箭头
                     show.legend = FALSE,max.overlaps = 100)

