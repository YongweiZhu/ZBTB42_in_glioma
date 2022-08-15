#
##
####    ZBTB42 promotes tumor progression in glioma
#   Supporter: Yongwei_Zhu
#   Email: zhuyongwei@csu.edu.cn
#   Version1.0
#   Update time: 20220414

setwd('/export3/zhuyw/Atask/20220414_ZBTB42/runtime/')
save.image("working.RData") #保存工作空间

## 基本配置
# 设置颜色
bioCol <- c('skyblue3', "deeppink3","deepskyblue3",'orange3', 
            'tomato3', 'turquoise3', 'seagreen3')


#
######      数据的读取与保存      ######
##
###    数据读取




##
###   数据保存
#   TCGA.GBMLGG
write.table(data.frame(SID=rownames(tmp), tmp),
            "./output/tableS1.1_TCGA.GBMLGG_ZBTB42_cutpoint.tsv",
            row.names = F, sep = "\t", quote = F)
write.table(data.frame(SID=rownames(top.table), top.table[,2:ncol(top.table)]),
            "./output/tableS1.2_TCGA.GBMLGG_ZBTB42_DEGs.tsv",
            row.names = F, sep = "\t", quote = F)
write.table(data.frame(GOID=rownames(ego@result), ego@result),
            "./output/tableS1.3_TCGA.GBMLGG_ZBTB42_GO_UP.tsv",
            row.names = F, sep = "\t", quote = F)
write.table(data.frame(GOID=rownames(ego@result), ego@result),
            "./output/tableS1.4_TCGA.GBMLGG_ZBTB42_GO_DOWN.tsv",
            row.names = F, sep = "\t", quote = F)
write.table(data.frame(SID=rownames(KEGG@result), KEGG@result),
            "./output/tableS1.5_TCGA.GBMLGG_ZBTB42_GSEA.KEGG.tsv",
            row.names = F, sep = "\t", quote = F)
write.table(data.frame(SID=rownames(top.table), top.table),
            "./output/tableS1.6_TCGA.GBMLGG_ZBTB42_GSVA.reactome.tsv",
            row.names = F, sep = "\t", quote = F)
write.table(data.frame(SID=rownames(top.table), top.table),
            "./output/tableS1.7_TCGA.GBMLGG_ZBTB42_GSVA.hallmaker.tsv",
            row.names = F, sep = "\t", quote = F)


#   TCGA.GBM
write.table(data.frame(SID=rownames(tmp), tmp),
            "./output/tableS2.1_TCGA.GBM_ZBTB42_cutpoint.tsv",
            row.names = F, sep = "\t", quote = F)
write.table(data.frame(SID=rownames(top.table), top.table[,2:ncol(top.table)]),
            "./output/tableS2.2_TCGA.GBM_ZBTB42_DEGs.tsv",
            row.names = F, sep = "\t", quote = F)
write.table(data.frame(GOID=rownames(ego@result), ego@result),
            "./output/tableS2.3_TCGA.GBM_ZBTB42_GO_UP.tsv",
            row.names = F, sep = "\t", quote = F)
write.table(data.frame(GOID=rownames(ego@result), ego@result),
            "./output/tableS2.4_TCGA.GBM_ZBTB42_GO_DOWN.tsv",
            row.names = F, sep = "\t", quote = F)
write.table(data.frame(SID=rownames(KEGG@result), KEGG@result),
            "./output/tableS2.5_TCGA.GBM_ZBTB42_GSEA.KEGG.tsv",
            row.names = F, sep = "\t", quote = F)
write.table(data.frame(SID=rownames(top.table), top.table),
            "./output/tableS2.6_TCGA.GBM_ZBTB42_GSVA.reactome.tsv",
            row.names = F, sep = "\t", quote = F)
write.table(data.frame(SID=rownames(top.table), top.table),
            "./output/tableS2.7_TCGA.GBM_ZBTB42_GSVA.hallmaker.tsv",
            row.names = F, sep = "\t", quote = F)


#   TCGA.LGG
write.table(data.frame(SID=rownames(tmp), tmp),
            "./output/tableS3.1_TCGA.LGG_ZBTB42_cutpoint.tsv",
            row.names = F, sep = "\t", quote = F)
write.table(data.frame(SID=rownames(top.table), top.table[,2:ncol(top.table)]),
            "./output/tableS3.2_TCGA.LGG_ZBTB42_DEGs.tsv",
            row.names = F, sep = "\t", quote = F)
write.table(data.frame(GOID=rownames(ego@result), ego@result),
            "./output/tableS3.3_TCGA.LGG_ZBTB42_GO_UP.tsv",
            row.names = F, sep = "\t", quote = F)
write.table(data.frame(GOID=rownames(ego@result), ego@result),
            "./output/tableS3.4_TCGA.LGG_ZBTB42_GO_DOWN.tsv",
            row.names = F, sep = "\t", quote = F)
write.table(data.frame(SID=rownames(KEGG@result), KEGG@result),
            "./output/tableS3.5_TCGA.LGG_ZBTB42_GSEA.KEGG.tsv",
            row.names = F, sep = "\t", quote = F)
write.table(data.frame(SID=rownames(top.table), top.table),
            "./output/tableS3.6_TCGA.LGG_ZBTB42_GSVA.reactome.tsv",
            row.names = F, sep = "\t", quote = F)
write.table(data.frame(SID=rownames(top.table), top.table),
            "./output/tableS3.7_TCGA.LGG_ZBTB42_GSVA.hallmaker.tsv",
            row.names = F, sep = "\t", quote = F)

write.table(data.frame(SID=rownames(DEGs.Union), DEGs.Union[, -1]),
            "./output/tableS4.1_DEgs.Union.tsv",
            row.names = F, sep = "\t", quote = F)

#
######      1-ZBTB42高低分组间的差异表达分析      ######
##   TCGA.GBMLGG   (G2, G3, G4)
##   TCGA.LGG   (G2, G3)
##   TCGA.GBM   (G4)

#   1.1_样本筛选
tmp_traits <- TCGA_GBMLGG_traits[TCGA_GBMLGG_traits$Grade %in% c("G2", "G3", "G4") & 
                                   substr(rownames(TCGA_GBMLGG_traits),14,15)<10 &
                                   !is.na(TCGA_GBMLGG_traits$OS) &
                                   !is.na(TCGA_GBMLGG_traits$OS.time) &
                                   TCGA_GBMLGG_traits$OS.time>0, ]
sample_list <- intersect(rownames(tmp_traits), colnames(tmp_exp))


#   1.2_分组
tmp <- data.frame(tmp_traits[sample_list, c("OS", "OS.time")],
                  RS=t(tmp_exp["ZBTB42", sample_list]))
colnames(tmp) <- c("OS", "OS.time", "RS")
#colnames(tmp) <- c("OS", "OS.time", "ZBTB42", "ZBTB42.l")

res.cut <- surv_cutpoint(tmp, time = "OS.time", event = "OS",
                         variables = "RS")
summary(res.cut)
plot(res.cut)

tmp$RS.l <- ifelse(tmp$RS>(6.7), 'High', 'Low')

fit_km <- survfit(Surv(OS.time, OS) ~ RS.l , data = tmp)
ggsurvplot(fit_km, pval=TRUE,
           conf.int =TRUE,
           xlab ="Time of days",
           ggtheme =theme_light(),
           risk.table = T, 
           palette = bioCol[c(4:3)],
           #legend.labs=c('High', 'Low'),
           legend.title='Risk Score')


#   1.3_差异分析
counts <- floor(2^tmp_exp[, sample_list]-1)
group <- tmp[sample_list, 'RS.l']      # High_vs_Low

top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
top.table$change = as.factor(ifelse(abs(top.table$logFC) > 0.6, ifelse(top.table$logFC > 0.6, 'UP','DOWN'),'NOT'))
top.table[top.table$adj.P.Val>0.05, "change"] <- "NOT"
top.table$symbol <- gsub(top.table$Gene, pattern = '_', replacement = '-')

top.table$GSEA.factor <- top.table$logFC*(-log10(top.table$adj.P.Val))      #   GSEA过程中需要的参数

top.table.filter <- top.table[top.table$adj.P.Val < 0.01 &
                                abs(top.table$logFC) > 1.5, ]
table(top.table$change)
table(top.table.filter$change)

#   GO富集分析
gene_list


#   GSEA富集分析(KEGG)
genelist <- top.table$GSEA.factor
names(genelist) <- top.table$symbol
genelist <- sort(genelist, decreasing = T, )
head(genelist)

gmtFile <- read.gmt("/export2/zhuyw/database/GSEA/c2.cp.kegg.v7.2.symbols.gmt")

KEGG <- GSEA(genelist, TERM2GENE = gmtFile, pvalueCutoff = 0.05)  #  GSEA分析
gseaplot2(KEGG, 1:10,
          pvalue_table = T)

#
######      2-ZBTB42高低分组间的基因集富集分析      ######
##   2.1_数据读取
# reactome
tmpdata <- read.table("/export2/zhuyw/database/TCGA/GBMLGG/TCGA_GBMLGG_c2.cp.reactome.v7.5.1_ssGSEA.scores.tsv", 
                      header = T, row.names = 1, sep = "\t")

# hallmaker
tmpdata <- read.table("/export2/zhuyw/database/TCGA/GBMLGG/TCGA_GBMLGG_hallmaker.v7.5_ssGSEA.scores.tsv", 
                      header = T, row.names = 1, sep = "\t")

GSEA_matrix <- tmpdata     #备份数据
tmpdata <- tmpdata[, sample_list]
group <- factor(tmp_traits[sample_list, "RS.l"])

mm <- model.matrix(~0 + group)
fit <- lmFit(tmpdata, mm)
head(coef(fit))
contr <- makeContrasts(groupHigh - groupLow, levels = colnames(coef(fit)))
tmpdata <- contrasts.fit(fit, contr)
tmpdata <- eBayes(tmpdata)

top.table <- topTable(tmpdata, sort.by = "P", number=Inf)
head(top.table, 10)

top.table$Pathway <- rownames(top.table)
top.table <- top.table[,c("Pathway", names(top.table)[1:6])]
summary(top.table$logFC)
top.table$change = as.factor(ifelse(abs(top.table$logFC) > 0.12, ifelse(top.table$logFC > 0.35, 'UP','DOWN'),'NOT'))
top.table$change[top.table$P.Value>0.05] <- "NOT"
table(top.table$change)

gene_list <- top.table[top.table$change %in% c("UP", "DOWN"), "Pathway"]

#   热图可视化
tmpdata <- tmp_traits[sample_list,
                      c("OS", "Histology", "Grade", "Gender", "IDH.status", "X1p.19q.codel", "IDH.codel.subtype", "Age.l", "RS.l")]
tmpdata <- tmpdata[order(tmpdata$RS.l), ]

x <- list()
for(i in 1:length(colnames(tmpdata))){
  x[[i]] <- as.factor(tmpdata[, i])
}
anno_col <- as.data.frame(x)
colnames(anno_col) <- colnames(tmpdata)
rownames(anno_col) <- rownames(tmpdata)

x <- list()
for(i in 1:length(colnames(tmpdata))){
  n = length(unique(sort(tmpdata[, i])))+2
  x[[i]] <- cols[3:n]
  names(x[[i]]) <- unique(sort(tmpdata[, i]))
}
names(x) <- colnames(tmpdata)
anno_colors <- x

pheatmap(scale(GSEA_matrix[gene_list, rownames(tmpdata)], center = F),
         cluster_rows = T,
         cluster_cols = FALSE,
         color=colorRampPalette(colors = c("deepskyblue3", "white", "tomato3"))(10),
         #display_numbers=data_mark, 
         fontsize_number=3.6,
         #gaps_row = c(2,13),
         #legend_breaks = c(-4,4),
         annotation_col = anno_col,
         annotation_colors = anno_colors,
         show_colnames = F)



#
######      3-ZBTB42高低分组间的免疫微环境的差异分析      ######
##   3.1_数据读取
# 2016.Cell.Reports_28.ICI
tmpdata <- read.table("/export2/zhuyw/database/TCGA/GBMLGG/TCGA_GBMLGG_2016.Cell.Reports.Sigs_ssGSEA.scores.tsv", 
                      header = T, row.names = 1, sep = "\t")
tmpdata <- as.data.frame(t(tmpdata[, sample_list]))

# 2021.Cancer.cell_29.ICI
tmpdata <- read.table("/export2/zhuyw/database/TCGA/GBMLGG/TCGA_GBMLGG_2021.Cancer.cel.Sigs_ssGSEA.scores.tsv", 
                      header = T, row.names = 1, sep = "\t")
tmpdata <- as.data.frame(t(tmpdata[, sample_list]))

# cibersort_22.ICI
tmpdata <- read.table("/export2/zhuyw/database/TCGA/TCGA_ICI/CIBERSORT/TCGA.GBMLGG.cibersort.relative.tsv", 
                      header = T, row.names = 2, sep = "\t")
tmpdata <- avereps(tmpdata[, 3:24], ID = tmpdata$SampleID)
tmpdata <- as.data.frame(tmpdata[sample_list, ])

# ESTIMATE_3.ICI
tmpdata <- read.table("/export2/zhuyw/database/TCGA/TCGA_ICI/ESTIMATE/TCGA.all.estimate.RNASeqV2.scores.tsv", 
                      header = T, row.names = 1, sep = "\t")
tmpdata <- as.data.frame(tmpdata[intersect(rownames(tmpdata), sample_list), 3:5])



#   箱线图可视化
tmpdata$Type <- tmp_traits[rownames(tmpdata), "RS.l"]
tmpdata <- melt(tmpdata, id.vars=c("Type"))
colnames(tmpdata)=c("Type", "Immune", "Fraction")
tmpdata$Type <- as.factor(tmpdata$Type)
tmpdata$Fraction <- as.numeric(tmpdata$Fraction)
#tmpdata$Fraction <- scale(tmpdata[,3], center = F)
#data$Fraction[data$Fraction>12]=12
ggboxplot(tmpdata, x="Immune", y="Fraction", fill = "Type", palette = bioCol[4:3], #fill="ICIcluster",
          ylab="Scale of Fraction",
          xlab="",
          legend.title="RiskScore") + 
  #palette=bioCol)
  rotate_x_text(50) + #scale_fill_brewer(palette = bioCol[4:3]) +#"Set3") +
  stat_compare_means(aes(group=Type),
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", "ns")),label = "p.signif")

#
######      4-ZBTB42高低分组间的干性指数的差异分析      ######
##   4.1_数据读取
# mRNAsi
tmpdata <- read.table("/export2/zhuyw/database/TCGA/PanCan.Stem/StemnessScores_RNAexp.tsv", 
                      header = T, sep = "\t")
tmpdata <- avereps(tmpdata[, 5:6], ID = tmpdata$TCGAlong.id)
tmpdata <- as.data.frame(tmpdata[intersect(sample_list, rownames(tmpdata)), ])

# mDNAsi
tmpdata <- read.table("/export2/zhuyw/database/TCGA/PanCan.Stem/StemnessScores_DNAmeth.tsv", 
                      header = T, sep = "\t")
tmpdata <- avereps(tmpdata[, 5:8], ID = tmpdata$TCGAlong.id)
tmpdata <- as.data.frame(tmpdata[intersect(sample_list, rownames(tmpdata)), ])


##   4.2_Verhaak单细胞数据集查看(干性特征，缺氧特征)
FeaturePlot(tmp_scRNA[,tmp_scRNA$condition_revalue %in% c("normoxia-3d", "normoxia-9d")], 
            features = c("PGBD1"), max.cutoff = 0.2,
            split.by = 'cell_type')

VlnPlot(tmp_scRNA[,tmp_scRNA$condition_revalue %in% c("normoxia-3d", "normoxia-9d")],
        features = c("PGBD1"),  
        split.by = "condition_revalue", cols = bioCol[2:7])

DotPlot(tmp_scRNA[,tmp_scRNA$condition_revalue %in% c("normoxia-3d", "normoxia-9d")],
        features = c("ZBTB42"),  
        split.by = "condition_revalue", cols = bioCol[2:7])


##   4.3_GSC细胞中常见的干性基因标志物
# CD133, CD44, Musashi-1, CD15, L1CAM, Integrin α6, Nestin, CD36, A2B5, LGR5, B23, GPD1
gene_list <- c("SOX2", "OLIG2", "PROM1", "CD44", "MSI1", "FUT4", "L1CAM", 
               "ITGA6", "NES", "CD36", "GFAP", "LGR5", "NPM1", "GPD1", "POU5F1")

gene_list <- intersect(gene_list, rownames(tmp_exp))
gene_list <- c(gene_list, "ZBTB42")

tmpdata <- t(tmp_exp[gene_list, sample_list])
M=cor(tmpdata)    #相关性矩阵
res1=cor.mtest(tmpdata,   #相关性矩阵的检验
               conf.level = 0.95)
corrplot(M, p.mat = res1$p,
         order="original",
         method = "square",
         type = "upper",
         tl.cex=0.8, pch=F, tl.srt = 45,
         insig = "label_sig",
         sig.level=c(.001, .01, .05), pch.cex=0.8, pch.col="black",
         #cl.lim = c(-0.3, 0.3), is.corr = F,
         col=colorRampPalette(c("deepskyblue3", "white", "tomato3"))(10),
         tl.col="black")


#
######      5-ZBTB42相关的hub基因的筛选      ######
tmp1 <- read.table("./output/tableS1.2_TCGA.GBMLGG_ZBTB42_DEGs.tsv", 
                   header = T, row.names = 1, sep = "\t")
tmp1 <- tmp1[abs(tmp1$logFC)>1 & tmp1$adj.P.Val<0.05, ]
table(tmp1$change)

tmp2 <- read.table("./output/tableS2.2_TCGA.GBM_ZBTB42_DEGs.tsv", 
                  header = T, row.names = 1, sep = "\t")
tmp2 <- tmp2[abs(tmp2$logFC)>1 & tmp2$adj.P.Val<0.05, ]
table(tmp2$change)

tmp3 <- read.table("./output/tableS3.2_TCGA.LGG_ZBTB42_DEGs.tsv", 
                   header = T, row.names = 1, sep = "\t")
tmp3 <- tmp3[abs(tmp3$logFC)>1 & tmp3$adj.P.Val<0.05, ]
table(tmp3$change)

##   Venn图
x <- list(A=tmp1[tmp1$change %in% c("DOWN"), "symbol"],
          B=tmp2[tmp2$change %in% c("DOWN"), "symbol"],
          C=tmp3[tmp3$change %in% c("DOWN"), "symbol"])
VennDiagram::venn.diagram(x, "test01.tiff")

gene_list <- intersect(rownames(tmp1), rownames(tmp2))
gene_list <- intersect(gene_list, rownames(tmp3))

DEGs.Union <- data.frame(Genes=gene_list)
rownames(DEGs.Union) <- DEGs.Union$Genes
DEGs.Union$VS1 <- tmp1[rownames(DEGs.Union), 'change']
DEGs.Union$VS2 <- tmp2[rownames(DEGs.Union), 'change']
DEGs.Union$VS3 <- tmp3[rownames(DEGs.Union), 'change']

tmp <- DEGs.Union[,2:4]
for(i in 1:ncol(tmp))
{
  tmp[ ,i] <- ifelse(is.na(tmp[,i]), 0,
                     ifelse(tmp[,i] %in% "UP", 1, -1))
}
colnames(tmp) <- paste(colnames(tmp), 'l', sep = '.')
tmp$change.sum <- apply(tmp, 1, sum)

DEGs.Union <- data.frame(DEGs.Union, tmp)

gene_list <- gsub(DEGs.Union[DEGs.Union$change.sum %in% c(3), "Genes"], 
                  pattern = "_", replacement = "-")

#  UniCox + LASSO + MultiCox
##   5.1_TCGA.LGG 数据集
# Single term deletions
# Model:
#  Surv(OS.time, OS) ~ KCNIP3 + IGFBP2 + SLC16A3 + DUSP26 + NUAK2 + 
#  EMILIN1 + FAM115C + C5orf62 + HOXA1 + APOBEC3C + CRTAC1 + 
#  PABPC5 + C2orf85 + AIM1 + BEND4 + IL1RAPL1 + GDF15 + CYP2S1 + 
#  KCTD4 + IL15 + SAMD9L + CALN1 + KLRC2 + PIK3R6 + C6orf15 + 
#  SCNN1B + GRIK2 + IGFN1 + AKR1C2 + ST8SIA2
gene_list <- c("KCNIP3", "IGFBP2", "SLC16A3", "DUSP26", "NUAK2", 
  "EMILIN1", "FAM115C", "C5orf62", "HOXA1", "APOBEC3C", "CRTAC1", 
  "PABPC5", "C2orf85", "AIM1", "BEND4", "IL1RAPL1", "GDF15", "CYP2S1", 
  "KCTD4", "IL15", "SAMD9L", "CALN1", "KLRC2", "PIK3R6", "C6orf15", 
  "SCNN1B", "GRIK2", "IGFN1", "AKR1C2", "ST8SIA2")

##   5.2_TCGA.GBM数据集
# Single term deletions
# Model:
#  Surv(OS.time, OS) ~ KCNIP3 + IGFBP2 + CRTAC1 + IL15 + SAMD9L
gene_list <- c("KCNIP3", "IGFBP2", "CRTAC1", "IL15", "SAMD9L")








###
######      调试-张炜_cuprotosis      ######
tmp1 <- read.table("./cuprotosis_cluster_tmp.tsv", 
                   header = T, row.names = 1, sep = "\t")

tmp <- tmp_traits
tmp <- tmp[tmp$PID %in% tmp1$PID, ]
tmp <- tmp1 %>% left_join(tmp, by = "PID")
tmp <- tmp[!is.na(tmp$OS) &
           !is.na(tmp$OS.time) & 
           tmp$OS.time > 0, ]

fit_km <- survfit(Surv(OS.time, OS) ~ cluster.id , data = tmp)
ggsurvplot(fit_km, pval=TRUE,
           conf.int =TRUE,
           xlab ="Time of days",
           ggtheme =theme_light(),
           risk.table = T, 
           palette = bioCol[c(4:2)],
           #legend.labs=c('High', 'Low'),
           legend.title='Cluster')

#  设置2组DEGs分析
#   DEGs.VS1: cluster1_vs_cluster2
#   DEGs.VS2: cluster3_vs_cluster2

#sample_list <- tmp$SID
counts <- floor(2^tmp_exp[, tmp$SID]-1)
group <- tmp[, 'cluster.id']      # BadSur_vs_GoodSur

top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
top.table$change = as.factor(ifelse(abs(top.table$logFC) > 0.6, ifelse(top.table$logFC > 0.6, 'UP','DOWN'),'NOT'))
top.table[top.table$adj.P.Val>0.05, "change"] <- "NOT"
top.table$symbol <- gsub(top.table$Gene, pattern = '_', replacement = '-')

top.table$GSEA.factor <- top.table$logFC*(-log10(top.table$adj.P.Val))      #   GSEA过程中需要的参数

top.table.filter <- top.table[top.table$adj.P.Val < 0.05 &
                                abs(top.table$logFC) > 1, ]
table(top.table$change)
table(top.table.filter$change)

#   差异分析结果
DEGs.VS1 <- top.table.filter
DEGs.VS2 <- top.table.filter

length(intersect(DEGs.VS1[DEGs.VS1$change %in% c("UP"), "symbol"],
                 DEGs.VS2[DEGs.VS2$change %in% c("UP"), "symbol"]))

length(intersect(DEGs.VS1[DEGs.VS1$change %in% c("DOWN"), "symbol"],
                 DEGs.VS2[DEGs.VS2$change %in% c("DOWN"), "symbol"]))








                   
