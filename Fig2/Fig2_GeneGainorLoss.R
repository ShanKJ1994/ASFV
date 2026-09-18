library(data.table)
library(tidyverse)
###################################

setwd("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/miniprot//")

ASFV=fread("./uniprot.annotation.final.stat.txt",header = F)
#ASFV=fread("./uniprot.annotation.final.stat20260901.txt",header = F)

colnames(ASFV)=c("GnmID",	"Start",	"End",	"Frameshift",	"Inframe_stp","Start_codon","GeneName","Ident",	"ProteinID",	"Group")
ASFV=filter(ASFV, !GnmID %in% c("OR944085.1", "OR944086.1", "OR944087.1", "OR944088.1", "OR944089.1", "OR944090.1", "OR944091.1"))
head(ASFV)
#ASFV=filter(ASFV,Frameshift==0, Inframe_stp==0,Start_codon ==1)

ASFV[GeneName == "100_c1"]$GeneName="100_1L2L"
ASFV[GeneName == "100_c2"]$GeneName="100_1R"
ASFV[GeneName == "100_c3"]$GeneName="100_3L"
ASFV[GeneName == "110_c9"]$GeneName="110_8L"
ASFV[GeneName == "110_c5"]$GeneName="110_1L"
ASFV[GeneName == "300_c1"]$GeneName="300_1L"
ASFV[GeneName == "300_c2"]$GeneName="300_2R"
ASFV[GeneName == "300_c3"]$GeneName="300_3L4L"
ASFV[GeneName == "360_c8"]$GeneName="360_19Rb"
ASFV[GeneName == "360_c9"]$GeneName="360_14L"
ASFV[GeneName == "360_c10"]$GeneName="360_15R"
ASFV[GeneName == "360_c3"]$GeneName="360_16R"
ASFV[GeneName == "360_c5"]$GeneName="360_18R"
ASFV[GeneName == "505_c1"]$GeneName="505_1R"
ASFV[GeneName == "505_c6"]$GeneName="505_3R"
ASFV[GeneName == "505_c4"]$GeneName="505_4R"
ASFV[GeneName == "505_c3"]$GeneName="505_11L"
ASFV[GeneName == "505_c5"]$GeneName="505_6R7R8R"
ASFV[GeneName == "505_c2"]$GeneName="505_2R5R9R10R"

head(ASFV)

Geno=fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/miniprot/meta.tsv",header = T)
colnames(Geno)[c(1,3)]=c("GnmID","Gnm_len")
head(Geno)

ASFV=merge(ASFV,Geno,by="GnmID")
head(ASFV)

# ASFV[!ASFV$genotype %in% c("1",  
#                           "1/2",
#                           "2",  
#                           "9"  )]$genotype="Other"

test=ASFV %>% group_by(GnmID,Gnm_len,genotype) %>% dplyr::summarise(count=n())
table(test$genotype)

# # 1. 定义指定基因型配色
# color_def <- c(
#   "1"    = "#E31A1C",
#   "1/2"  = "#33A02C",
#   "2"    = "#1F78B4",
#   "9"    =  "#E18727"
# )
# 
# # 2. 提取数据里剩下的基因型，统一灰色
# other_gt <- setdiff(unique(test$genotype), names(color_def))
# color_final <- c(color_def, setNames(rep("#AAAAAA", length(other_gt)), other_gt))

color_final <- c(
  "1" = "#E41A1C", "2" = "#377EB8", "9" = "#FF7F00","8"="#bf5500", "20" ="#FFC0CB","22"="#EE4C97FF","1/2" = "#4DAF4A",
  "Outgroup" = "black", "unknown" = "grey70"
)

# 3. 绘图
ggplot(test,mapping=aes(x=log10(Gnm_len),y=count))+
  geom_point(data=test,mapping=aes(color=genotype))+
  geom_smooth(method = "glm", se = F, linewidth = 1,color="black") +
  theme_classic()+
  ggpubr::stat_cor(method = "spearman", label.y = 205, col="red",cor.coef.name = "rho")+
  labs(y="gene number")+
  scale_color_manual(values = color_final)
ggsave(
  filename = "FigS_length_geneNum.pdf",
  width = 4,
  height = 3,
  device = "pdf"
)
# 
# 
# head(test)
# ggplot(test,mapping=aes(x=genotype,y=count,color=genotype))+
#   geom_boxplot(outliers = F)+
#   geom_jitter()+
#   #geom_smooth(method = "glm", se = F, linewidth = 1,color="black") +
#   theme_classic()+
#   #ggpubr::stat_cor(method = "spearman", label.y = 205, col="red",cor.coef.name = "rho")+
#   labs(y="gene number")+
#   scale_color_manual(values = color_final)
# 
# ggsave(
#   filename = "FigS_length_geneNum.pdf",
#   width = 4,
#   height = 3,
#   device = "pdf"
# )
# 
# ggplot(test,mapping=aes(x=genotype,y=log10(Gnm_len),color=genotype))+
#   geom_boxplot(outliers = F)+
#   geom_jitter()+
#   #geom_smooth(method = "glm", se = F, linewidth = 1,color="black") +
#   theme_classic()+
#   #ggpubr::stat_cor(method = "spearman", label.y = 205, col="red",cor.coef.name = "rho")+
#   labs(y="gene number")+
#   scale_color_manual(values = color_final)

table(test$genotype)
# Geno=unique(ASFV[,c("GnmID","Gnm_len","genotype")])
# ggplot(filter(Geno,genotype %in% c("1/2","2","1","9")),aes(x=genotype,y=log10(Gnm_len),fill=genotype))+
#   geom_boxplot()+
#   scale_color_manual(values=color_final)+
#   theme_classic()+
#   ggpubr::stat_compare_means(method = 'wilcox.test', comparisons = list(c("1","2"),
#                                                                         c("1","1/2"),
#                                                                         c("1/2","2"),
#                                                                         c("9","2"),
#                                                                         c("9","1/2"),
#                                                                         c("1","9")),label = 'p.signif')#+ 
# #coord_flip()

#############################################

#gene number in one genome

#############################################

mid=ASFV %>% group_by(GnmID,genotype) %>% dplyr::summarise(count=n())
head(mid)


ggplot(mid,aes(x=count,fill=genotype))+
  geom_histogram(binwidth=1)+
  theme_classic()+
  labs(x="gene number in one genome",y="genome number")+
  scale_fill_manual(values = color_final)

ggsave(
  filename = "FigS_geneNum_distribution.pdf",
  width = 4,
  height = 3,
  device = "pdf"
)

quantile(mid$count,0.025)
quantile(mid$count,0.975)
median(mid$count)
#test<-mid %>% group_by(genotype) %>% dplyr::summarise(median=median(count),L=quantile(count,0.025),H=quantile(count,0.975))

#############################################

#gene loss in groups

#############################################

# 1. 获取完整基因组列表（基准）
all_gnm <- unique(ASFV$GnmID)

mid_core <- ASFV %>%
  group_by(GeneName) %>%  # 按 OrthoGroup 分组
  summarise(
    gene_count=n(),
    # 统计1：该 OrthoGroup 已有的唯一基因组
    unique_GnmID_count = n_distinct(GnmID, na.rm = TRUE),  # 已有数量
    #unique_GnmID_list = paste(sort(unique(GnmID[!is.na(GnmID)])), collapse = ", "),  # 已有列表（排序后更整齐）
    # 统计2：该 Group 缺失的基因组
    loss_GnmID_count = length(setdiff(all_gnm, unique(GnmID[!is.na(GnmID)]))),  # 缺失数量
    loss_Gnm_list = paste(
      sort(setdiff(all_gnm, unique(GnmID[!is.na(GnmID)]))),
      collapse = ", "  # 缺失列表
    ),
    # 可选统计：该 OrthoGroup 总记录数（原数据中该组有多少行）
    total_records = n(),
    # 可选统计：该 OrthoGroup 涉及的唯一基因数（如果需要）
    unique_gene_count = n_distinct(GeneName, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  # 优化显示：缺失列表为空时显示“无”（更直观）
  mutate(loss_Gnm_list = ifelse(loss_Gnm_list == "", "无", loss_Gnm_list)) %>%
  # 可选：按“缺失数量”降序排序（优先看缺失多的组）
  arrange(desc(loss_GnmID_count))

df_core <- filter(ASFV,GeneName %in% filter(mid_core,unique_GnmID_count== gene_count,loss_GnmID_count/357<=0.05 )$GeneName)
unique(df_core$GeneName)
unique(filter(ASFV,!GeneName %in% filter(mid_core,unique_GnmID_count== gene_count,loss_GnmID_count/357<=0.05 )$GeneName)$GeneName)
unique(ASFV$Group)
#############################################

#gene loss in gnm position

#############################################

#head(mid)
pos=filter(ASFV,!GeneName %in% df_core$GeneName)

pos$pos <- pmin(pos$Start, pos$End)

total_len <- 190000
left_var_end    <- 40000
right_var_start <- total_len - 25000
left_ratio_end    <- left_var_end / total_len
right_ratio_start <- right_var_start / total_len

# 先预计算直方图y峰值
p_pre <- ggplot(pos, aes(x = pos / Gnm_len)) + geom_histogram(binwidth = 0.005)
build_p <- ggplot_build(p_pre)
y_max <- max(build_p$data[[1]]$count, na.rm = TRUE)
bar_top <- y_max
genome_h <- 30

ggplot(pos, aes(x = pos / Gnm_len, fill = genotype)) +
  geom_histogram(binwidth = 0.005) +
  # 顶层基因组矩形
  annotate("rect", xmin = 0, xmax = left_ratio_end,
           ymin = bar_top, ymax = bar_top + genome_h, fill = "#d9d9d9") +
  annotate("rect", xmin = left_ratio_end, xmax = right_ratio_start,
           ymin = bar_top, ymax = bar_top + genome_h, fill = "#8c8c8c") +
  annotate("rect", xmin = right_ratio_start, xmax = 1,
           ymin = bar_top, ymax = bar_top + genome_h, fill = "#d9d9d9") +
  annotate("text", x = left_ratio_end/2, y = bar_top + genome_h/2, label = "left\nvariable", size = 2.8) +
  annotate("text", x = (left_ratio_end + right_ratio_start)/2, y = bar_top + genome_h/2, label = "core", size = 2.8, color = "white") +
  annotate("text", x = (right_ratio_start + 1)/2, y = bar_top + genome_h/2, label = "right\nvariable", size = 2.8) +
  theme_classic() +
  scale_fill_manual(values = color_final) +
  labs(x = "relative genome position", y = "gene number") +
  scale_x_continuous(limits = c(0, 1)) # +
  # expand_limits(y = c(0, bar_top + genome_h))

ggsave(
  filename = "FigA_noncoregene_distribution.pdf",
  width = 4,
  height = 3,
  device = "pdf"
)


#############################################

#noncore gene number heatmap

#############################################

ASFV$GnmID2=paste0(ASFV$genotype,":",ASFV$GnmID)
head(ASFV)
ASFV$GeneName=sub("_YP\\d+.\\d_","_",ASFV$GeneName)

#ASFV12=filter(ASFV,genotype %in% c("1","2","1/2"),GnmID2 %in% mid$Var)
ortho=filter(ASFV,!GeneName %in% df_core$GeneName) %>% group_by(GnmID,GnmID2,GeneName) %>% dplyr::summarise(count=n())
head(ortho)
ortho=reshape2::dcast(ortho,GnmID+GnmID2~GeneName,value.var = "count" ,fill = 0)
# filter(ortho,count>2)
# ortho=data.frame(fread("Orthogroups.GeneCount.csv")[,-c("Total")])

ortho[1:4,1:4]
rownames(ortho)=ortho$GnmID2
ortho=ortho[,-c(1,2)]
head(ortho)
#ortho2 <- ortho[, colSums(ortho, na.rm = TRUE) >= 10]
#ortho2[1:4,1:4]

library(pheatmap)
# 7 组 ColorBrewer 定性色，一眼区分
#my_colors <- c( "white", "#377EB8", "#E41A1C", "#984EA3", #"#FF7F00",
#                "#FFFF33", "black")
#max(ortho)
group_pos=filter(ASFV,!GeneName %in% df_core$GeneName)

#order=unique(group_pos$GeneName)

blues <- colorRampPalette(c("white","lightblue", "skyblue","blue", "purple"))(100)  
max(ortho[order(rownames(ortho)),])
pheatmap((ortho[order(rownames(ortho)),]),
         #gaps_row = gaps_split[1:(length(gaps)-1)]  ,# border_color = 'white', 
         #color = my_colors,
         color = blues,
         legend = TRUE,
         #filename = "Groups.GeneCount.uniprot.all.jpg",
         height = 35, width = 25,
         cluster_rows = FALSE,
         cluster_cols = F,
         angle_col = 90,
         fontsize = 5)
#dev.off()


ASFV$GnmID2=paste0(ASFV$genotype,":",ASFV$GnmID)
head(ASFV)
ASFV$GeneName=sub("_YP\\d+.\\d_","_",ASFV$GeneName)


######################################

#PCA

######################################

library(data.table)
library(tidyverse)
library(vegan)
library(pairwiseAdonis)
# library(devtools)


pca_base <- prcomp(t(ortho), center = F)   # 已提前标准化
summary(pca_base)
pca=pca_base$rotation
data.frame(pca)
pca<- data.table(GnmID= rownames(pca), PC1 = data.frame(pca)$PC1, PC2 = data.frame(pca)$PC2)

Geno=unique(ASFV[,c("GnmID2","genotype")])

# 左连接，新增 genotype 列（缺失值自动为 NA）
pca[Geno, genotype := i.genotype, on = .(GnmID == GnmID2)]
# 把 NA 改成 "unknown"
pca[is.na(genotype), genotype := "unknown"]

# pca2=merge(pca,unique(ASFV[,c("GnmID","Gnm_len")]),by=c("GnmID")) # 你这行本次PCA绘图暂时不用，保留注释

#========= 1. 定义配色（严格匹配你的规则）=========
color_final <- c(
  "1" = "#E41A1C", "2" = "#377EB8", "9" = "#FF7F00","8"="#bf5500", "20" ="#FFC0CB","22"="#EE4C97FF","1/2" = "#4DAF4A",
  "Outgroup" = "black", "unknown" = "grey70"
)
#========= 2. PERMANOVA 组间差异检验（1,1/2,2,9）=========
pca_sub <- pca[genotype %in% c("1","1/2","2","9")]
# 欧式距离，基于PC1+PC2
dist_eu <- vegdist(pca_sub[, .(PC1, PC2)], method = "euclidean")
set.seed(123)
perma_global <- adonis2(dist_eu ~ genotype, data = pca_sub, permutations = 999)
perma_pair <- pairwise.adonis2(dist_eu ~ genotype, data = pca_sub, permutations = 999, p.adjust.m = "fdr")
disp <- betadisper(dist_eu, group = pca_sub$genotype)
anova(disp)

# 提取标注文本
r2_val <- round(perma_global$R2[1], 3)
p_val  <- round(perma_global$`Pr(>F)`[1], 3)
stat_label <- paste0("PERMANOVA\nR² = ", r2_val, ", P = ", p_val)

#========= 3. PCA绘图代码（修改完善版）=========
ggplot(data = pca,aes(x=PC1,y=PC2,col=genotype))  +
  labs(x = "PC1 (91.39%)", y = "PC2 (0.39%)") +
  geom_point(size = 1.8) +
  stat_ellipse(
    data = pca[genotype %in% c("1", "1/2", "2", "9")], # data.table筛选，和filter等价
    level = 0.95,
    linewidth = 0.8
  ) +
  # 图内右上角添加PERMANOVA全局检验结果
  annotate("text", x = -0.058, y = 0.105, label = stat_label, size = 3.5) +
  scale_color_manual(values = color_final) +
  theme_classic()

# 导出图片示例
ggsave("FigB_PCA_noncore_gene_ellipse.pdf", width = 5, height = 3, device = "pdf")


############################################

#tree + noncore heatmap

############################################
library(ape)
library(ggtree)
library(treeio)
library(tidytree)

# ---------------------- 1. 读取并定根进化树 ----------------------
tree_file <- "C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/ASFV_364_2_tree.root.nwk" 
tree <- read.tree(tree_file)
# 定义要删掉的tip名称
remove_tips <- c("OR944085.1", "OR944086.1", "OR944087.1", 
                 "OR944088.1", "OR944089.1", "OR944090.1", "OR944091.1")

# 删除指定tip
tree_pruned <- drop.tip(tree, tip = remove_tips)

# 可选：输出修剪后的新树
#write.tree(tree_pruned, file = "C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/ASFV_357_2_tree.root.nwk")
root_tips <- c("LC506465.1", "LC637659.1")
if (!all(root_tips %in% tree$tip.label)) {
  missing_tips <- setdiff(root_tips, tree$tip.label)
  stop(paste("以下根序列不在进化树中：", paste(missing_tips, collapse = ", ")))
}
tree_rooted <- root(tree_pruned, outgroup = root_tips, resolve.root = TRUE)
tree_rooted  <- drop.tip(tree_rooted , tip = root_tips)


genotype_df <-read.table("./genotype.tsv",header = T)
#root_df=data.frame(SampleID=c("LC506465.1", "LC637659.1"),Genotype=c("Outgroup","Outgroup"))
colnames(genotype_df)=c("SampleID","Genotype")
#genotype_df=rbind(genotype_df,root_df)
# 按进化树顺序重排，填充NA为Unknown
genotype_df_ordered <- genotype_df[match(tree_rooted$tip.label, genotype_df$SampleID), ]
genotype_df_ordered$Genotype[is.na(genotype_df_ordered$Genotype)] <- "unknown"
# 构建样本ID→基因型的映射
tip_genotype <- setNames(genotype_df_ordered$Genotype, tree_rooted$tip.label)
# 构建基因型→颜色的映射（与之前一致）
genotype_levels <- unique(genotype_df_ordered$Genotype)
genotype_colors <- c(
  "1" = "#E41A1C", "2" = "#377EB8", "9" = "#FF7F00","8"="#bf5500", "20" ="#FFC0CB","22"="#EE4C97FF","1/2" = "#4DAF4A",
  "Outgroup" = "black", "unknown" = "grey70"
)


p_tree <- ggtree(tree_rooted,  
                 layout = "rectangular",  
                 branch.length = "none", 
                 #panel.width = 2,
                 lwd = 1) +
  #theme_tree2() +
  geom_tippoint(aes(color = tip_genotype[label]), size = 2, shape = 19) +
  geom_tree(aes(color = tip_genotype[label]), size = 1) +
  #geom_tiplab(size = 1.5, offset = 0.5) +  # 标签偏移，避免和条形图重叠
  scale_color_manual(values = genotype_colors, name = "Genotype")# +
  #theme_tree2(legend.position = "none")


ortho=filter(ASFV,!GeneName %in% df_core$GeneName) %>% group_by(GnmID,GeneName) %>% dplyr::summarise(count=n())
head(ortho)
table(ASFV$genotype)
ortho=reshape2::dcast(ortho,GnmID~GeneName,value.var = "count" ,fill = 0)
rownames(ortho)=ortho$GnmID
heatmap_df <- ortho[,-c(1)]
heatmap_df[is.na(heatmap_df)]

heatmap_df <-heatmap_df[rownames(heatmap_df) %in% tree_rooted$tip.label,] 
heatmap_matrix_ordered <- heatmap_df[tree_rooted$tip.label, ]
all(rownames(heatmap_matrix_ordered) == tree_rooted$tip.label)  

# 填充热图NA值
heatmap_matrix_ordered[is.na(heatmap_matrix_ordered)] <- 0


# 整合热图
# 精确匹配这 5 个前缀
#to_remove <- paste0("^(MGF)(", paste(c(110,100,300,360,505,148,285,86), collapse="|"), ")")
ASFV$pos <- pmax(ASFV$Start, ASFV$End)
mid=filter(ASFV,!GeneName %in% df_core$GeneName) %>% group_by(GeneName) %>% dplyr::summarise(max=max(pos/Gnm_len))
mid=mid[order(mid$max,decreasing = F),]
#colnames(heatmap_matrix_ordered)[!colnames(heatmap_matrix_ordered) %in% mid$GeneName]

p_tree
plot=(as.matrix(heatmap_matrix_ordered[,mid$GeneName]))
colnames(plot)=gsub("YP00.*_","",colnames(plot))
#colnames(plot) <- gsub("ASFV_G_", "", colnames(plot))
table(as.vector(plot))

final_plot <- p_tree %>%
  gheatmap(
    plot,
    colnames       = TRUE,
    #    rownames       = TRUE,
    font.size = 2,
    colnames_angle = 90,
    color          = "white"          # 仅控制瓷砖边框颜色
  ) +
  # 3. 再补一次 tip 标签（基因组名）
  # geom_tiplab(aes(label = label),
  #             size  = 2,
  #             angle = 0,
  #             hjust = 0,      # 向左对齐
  #             x     = max(node.depth(tree_rooted)) + 1) + # 热图右侧再往外一点
  scale_fill_gradientn(
    colours = colorRampPalette(c("white","lightblue", "skyblue","blue", "purple"))(7),
    #limits  = c(0, 12),
    name    = "Gene number"
  )

print(final_plot)      # 一定要显式打印

ggsave("FigS_noncore_gene_heatmap.pdf",final_plot,width =20, height =20, device = "pdf")


final_plot <- p_tree %>%
  gheatmap(
    plot,
    colnames       = TRUE,
    #    rownames       = TRUE,
    font.size = 2,
    colnames_angle = 90,
    color          = "white"          # 仅控制瓷砖边框颜色
  ) +
  # 3. 再补一次 tip 标签（基因组名）
  # geom_tiplab(aes(label = label),
  #             size  = 2,
  #             angle = 0,
  #             hjust = 0,      # 向左对齐
  #             x     = max(node.depth(tree_rooted)) + 1) + # 热图右侧再往外一点
  scale_fill_gradientn(
    colours = colorRampPalette(c("white","lightblue", "skyblue","blue", "purple"))(7),
    #limits  = c(0, 12),
    name    = "Gene number"
  )

print(final_plot)      # 一定要显式打印



ortho=filter(ASFV,Group %in% c(100,110,300,360,505)) %>% group_by(GnmID,Group) %>% dplyr::summarise(count=n())
head(ortho)
ortho=reshape2::dcast(ortho,GnmID~Group,value.var = "count" ,fill = 0)
rownames(ortho)=ortho$GnmID
heatmap_df <- ortho[,-c(1)]
heatmap_df[is.na(heatmap_df)]

heatmap_df <-heatmap_df[rownames(heatmap_df) %in% tree_rooted$tip.label,] 
heatmap_matrix_ordered <- heatmap_df[tree_rooted$tip.label, ]
all(rownames(heatmap_matrix_ordered) == tree_rooted$tip.label)  

# 填充热图NA值
heatmap_matrix_ordered[is.na(heatmap_matrix_ordered)] <- 0

p_tree
plot=(as.matrix(heatmap_matrix_ordered))
colnames(plot)=gsub("YP00.*_","",colnames(plot))
#colnames(plot) <- gsub("ASFV_G_", "", colnames(plot))
max(plot)
str(plot)
colnames(plot)=paste0("MGF",colnames(plot))
final_plot <- p_tree %>%
  gheatmap(
    data.frame(plot)[("MGF100")],
    offset=0,
    width=0.06,
    colnames=TRUE,
    font.size=3,
    color=NA
  ) +
  # 3. 再补一次 tip 标签（基因组名）
  # geom_tiplab(aes(label = label),
  #             size  = 2,
  #             angle = 0,
  #             hjust = 0,      # 向左对齐
  #             x     = max(node.depth(tree_rooted)) + 1) + # 热图右侧再往外一点
  scale_fill_gradientn(
    colours = colorRampPalette(c("#FFEBEE","#C62828"))(7),
    #limits  = c(0, 12),
    name    = "MGF-100 number"
  )

print(final_plot)      # 一定要显式打印


final_plot <- final_plot +
  ggnewscale::new_scale_fill()


final_plot  <- gheatmap(
  final_plot ,
  data.frame(plot)[("MGF110")],
  offset=2,
  width=0.06,
  colnames=TRUE,
  font.size=3,
  color=NA
)+
  
  scale_fill_gradientn(
    colours = colorRampPalette(c("#E8EAF6","#3949AB"))(7),
    #limits  = c(0, 12),
    name    = "MGF-110 number"
  )

print(final_plot)      # 一定要显式打印


final_plot <- final_plot +
  ggnewscale::new_scale_fill()

#min(data.frame(plot)[("MGF110")])
final_plot  <- gheatmap(
  final_plot ,
  data.frame(plot)[("MGF300")],
  offset=4,
  width=0.06,
  colnames=TRUE,
  font.size=3,
  color=NA
)+
  
  scale_fill_gradientn(
    
    colours = colorRampPalette(c("#FFFDE7","#F9A825"))(7),
    #limits  = c(0, 12),
    name    = "MGF-300 number"
  )

print(final_plot)      # 一定要显式打印



final_plot <- final_plot +
  ggnewscale::new_scale_fill()

#min(data.frame(plot)[("MGF110")])
final_plot  <- gheatmap(
  final_plot ,
  data.frame(plot)[("MGF360")],
  offset=6,
  width=0.06,
  colnames=TRUE,
  font.size=3,
  color=NA
)+
  
  scale_fill_gradientn(
    colours = colorRampPalette(c("#EFEBE9","#6D4C41"))(7),
    #limits  = c(0, 12),
    name    = "MGF-360 number"
  )

print(final_plot)      # 一定要显式打印


final_plot <- final_plot +
  ggnewscale::new_scale_fill()

#min(data.frame(plot)[("MGF110")])
final_plot  <- gheatmap(
  final_plot ,
  data.frame(plot)[("MGF505")],
  offset=8,
  width=0.06,
  colnames=TRUE,
  font.size=3,
  color=NA
)+
  
  scale_fill_gradientn(
    colours = colorRampPalette(c("#E0F2F1","#00897B"))(7),
    #limits  = c(0, 12),
    name    = "MGF-505 number"
  )

print(final_plot)      # 一定要显式打印
max(plot)
ggsave("FigC_MGF_heatmap.pdf",final_plot,width =20, height =20, device = "pdf")


# final_plot <- p_tree %>%
#   gheatmap(
#     plot,
#     colnames       = TRUE,
#     font.size      = 2,
#     colnames_angle = 90,
#     color          = "white"
#   ) +
#   scale_fill_gradientn(
#     # 按照您的要求，依次填入4个区间的颜色（每个区间从浅到深）
#     colours = c(
#       "white", "#1F78B4",  # 1-3
#       "#FDBF6F", "#FF7F00",  # 4-10
#       "#B2DF8A", "#33A02C",  # 11-15
#       "#FB9A99", "#E31A1C"   # 16-20
#     ),
#     # 定义区间的断点 (共5个断点，划分出4个区间)
#     breaks = c(0, 5, 10, 15, 20), 
#     # 锁定坐标轴范围，确保颜色准确映射
#     limits = c(0, 20),         
#     name   = "Gene number"
#   )
# final_plot
############################################

#MGF gene number

############################################


# 1. 数据筛选聚合（你的原有代码）
# test <- filter(ASFV,
#                genotype %in% c(1,2,9,"1/2"),
#                Group %in% c(100,110,300,360,505)) %>%
#   group_by(GnmID, Group, genotype) %>%
#   dplyr::summarise(count = n(), .groups = "drop")
# 
# ggplot(test, aes(x = genotype, y = count, fill = genotype)) +
#   #geom_violin() +
#   geom_boxplot(outliers = F) +
#   #geom_jitter(aes(col = genotype),alpha=0.1) +
#   theme_classic() +
#   labs(y = "gene number per genome") +
#   scale_fill_manual(values = color_final) +
#   facet_wrap(~Group, scale = "free") +
#   theme(
#     legend.position = "none",
#     strip.background = element_blank(),
#     strip.text.x = element_text(size = 11)
#   )
# 



# 1. 加载必要的包
library(dplyr)
library(ggplot2)
library(rstatix)
library(multcompView)

# 2. 你的原始数据（这里直接接上你的代码）
test <- filter(ASFV,
               genotype %in% c(1,2,9,"1/2"),
               Group %in% c(100,110,300,360,505)) %>%
  group_by(GnmID, Group, genotype) %>%
  dplyr::summarise(count = n(), .groups = "drop")

# 1. 使用原生 R 函数进行两两比较（极其稳健，不会报错）
stat_results_list <- test %>%
  group_split(Group) %>%  # 按 Group 拆分成列表
  map(~ {
    # 对每个子集进行原生 wilcox 检验，自动进行 BH 校正
    wt <- pairwise.wilcox.test(
      x = .x$count, 
      g = .x$genotype, 
      p.adjust.method = "BH",
      exact = FALSE  # 防止样本量相同时报错
    )
    
    # 将结果整理成数据框
    if (!is.null(wt$p.value)) {
      as.data.frame(wt$p.value, stringsAsFactors = FALSE) %>%
        rownames_to_column("group1") %>%
        pivot_longer(-group1, names_to = "group2", values_to = "p") %>%
        mutate(Group = unique(.x$Group)) %>%
        drop_na(p)  # 去掉因为组数不足产生的 NA
    } else {
      NULL
    }
  })

# 2. 将列表合并为一个完整的数据框
stat_results <- bind_rows(stat_results_list)




test <- filter(ASFV,
               genotype %in% c(1,2,9,"1/2"),
               Group %in% c(100,110,300,360,505)) %>%
  group_by(GnmID, Group, genotype) %>%
  dplyr::summarise(count = n(), .groups = "drop")  %>%
  group_by(Group, genotype) %>%
  dplyr::summarise(median = median(count),L=quantile(count,0.025),H=quantile(count,0.975), .groups = "drop")

test$genotype=factor(test$genotype,levels = c(1,"1/2",2,9))
ggplot(test, aes(x = genotype, y = median, fill = genotype)) +
  geom_col() +
  geom_errorbar(aes(ymin = L, ymax = H), width = 0.2, linewidth = 0.6) +
  theme_classic() +
  labs(y = "gene number per genome") +
  scale_fill_manual(values = color_final) +
  scale_x_discrete(limits = c("1","1/2","2","9")) +
  facet_wrap(~Group, scale = "free",nrow=1) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text.x = element_text(size = 11)
  )
ggsave("FigD_MGF_gene_number.pdf",width = 15, height = 3, device = "pdf")
#boxplot还是point？


test <- filter(ASFV,
               genotype %in% c(1,2,9,"1/2"),
               Group %in% c("VFE66","VF69R","PNG3","VF93L")) %>%
  group_by(GnmID, Group, genotype) %>%
  dplyr::summarise(count = n(), .groups = "drop")  %>%
  group_by(Group, genotype) %>%
  dplyr::summarise(median = median(count),L=quantile(count,0.025),H=quantile(count,0.975), .groups = "drop")

test[test$Group=="PNG3",]$Group="Uncharacterized membrane\nprotein PNG3"
test[test$Group=="VF69R",]$Group="Uncharacterized membrane\nprotein VF69R"
test[test$Group=="VF93L",]$Group="Uncharacterized membrane\nprotein VF93L"
test[test$Group=="VFE66",]$Group="Host translation \ninhibitor VFE66"

str(test)
test$genotype=as.character(test$genotype)
test$genotype=factor(test$genotype,levels = c("1","1/2","2","9"))
ggplot(test, aes(x = genotype, y = median, fill = genotype)) +
  geom_col() +
  geom_errorbar(aes(ymin = L, ymax = H), width = 0.2, linewidth = 0.6) +
  theme_classic() +
  labs(y = "gene number per genome") +
  scale_fill_manual(values = color_final) +
  scale_x_discrete(limits = c("1","1/2","2","9")) +
  facet_wrap(~Group, scale = "free",nrow=1) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text.x = element_text(size = 11)
  )
ggsave("FigS_function_gene_number.pdf",width = 12, height = 3, device = "pdf")

test=ASFV
test$Group=gsub("_YP\\w+.\\d+","",test$Group)
unique(test$Group[grepl("ASFV",test$Group)])

test <- filter(test,
               genotype %in% c(1,2,9,"1/2"),
               Group %in% c("KP86R","DP60R","PNG5","PNG7","IM014gp072","IM014gp185",
                            "ASFV_G_ACD_00190","ASFV_G_ACD_00270","ASFV_G_ACD_00330",
                            "ASFV_G_ACD_00350","ASFV_G_ACD_00360",
                            "ASFV_G_ACD_01980","ASFV_G_ACD_01990")) %>%
  group_by(GnmID, Group, genotype) %>%
  dplyr::summarise(count = n(), .groups = "drop")  %>%
  group_by(Group, genotype) %>%
  dplyr::summarise(median = median(count),L=quantile(count,0.025),H=quantile(count,0.975), .groups = "drop")

test$genotype=factor(test$genotype,levels = c(1,"1/2",2,9))
ggplot(test, aes(x = genotype, y = median, fill = genotype)) +
  geom_col() +
  geom_errorbar(aes(ymin = L, ymax = H), width = 0.2, linewidth = 0.6) +
  theme_classic() +
  labs(y = "gene number per genome") +
  scale_fill_manual(values = color_final) +
  scale_x_discrete(limits = c("1","1/2","2","9")) +
  facet_wrap(~Group, scale = "free") +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text.x = element_text(size = 11)
  )
ggsave("FigS_nonkownfunction_gene_number.pdf",width = 6, height = 6, device = "pdf")

# ============================================================================

# MGF gene example

# ============================================================================


library(data.table)
library(gggenes)
library(ggplot2)
library(dplyr)
library(tidyr)

# ============================================================================
# 1. 读取数据
# ============================================================================
setwd("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/miniprot/")
ASFV=fread("./ASFV_uniprot.annotation.final.stat.txt",header = T)
#colnames(ASFV)=c("GnmID",	"Start",	"End",	"Frameshift",	"Inframe_stp","Start_codon","GeneName","Ident",	"ProteinID",	"Group")
ASFV=filter(ASFV, !GnmID %in% c("OR944085.1", "OR944086.1", "OR944087.1", "OR944088.1", "OR944089.1", "OR944090.1", "OR944091.1"))

df <- ASFV
unique(Geno$genotype)
# ============================================================================
# 2. 方向计算（+/- 链正确）
# ============================================================================
df <- df %>%
  filter(!(genotype == "unknown" & GnmID != "NC_044951.1")) %>%
  mutate(
    strand = ifelse(Start < End, "+", "-"),
    # ✅ 关键：只根据 strand 交换坐标，保留原始位置
    xmin_fit = ifelse(strand == "+", Start, End),
    xmax_fit = ifelse(strand == "+", End, Start)
  )

df <- rbind(filter(df,!GnmID %in% unique(filter(ASFV,Group==110,genotype==1)$GnmID),!GnmID %in% unique(filter(ASFV,Group==110,genotype==2)$GnmID)),
            filter(df,GnmID == "NC_044941.1"),
            filter(df,GnmID == "C_AA040130.1"))

# ============================================================================
# 3. 核心/可变基因定义
# ============================================================================
total_per_genotype <- df %>%
  group_by(genotype) %>%
  summarise(total_genome = n_distinct(GnmID), .groups = "drop")

gene_freq <- df %>%
  group_by(genotype, GeneName) %>%
  summarise(n_genome = n_distinct(GnmID), .groups = "drop") %>%
  left_join(total_per_genotype, by = "genotype") %>%
  mutate(
    freq = n_genome / total_genome,
    core_type = ifelse(freq >= 0.95, "core", "variable")
  )

df <- left_join(df, gene_freq, by = c("genotype", "GeneName"))

# ============================================================================
# 4. 每个基因型只选1条
# ============================================================================
representative <- df %>%
  group_by(genotype) %>%
  summarise(
    has_NC = any(grepl("^NC_", GnmID)),
    selected = ifelse(has_NC, list(unique(GnmID[grepl("^NC_", GnmID)])[1]), list(sample(unique(GnmID), 1))),
    .groups = "drop"
  ) %>%
  unnest(selected) %>%
  pull(selected)

df_plot <- df %>% filter(GnmID %in% representative)

# ============================================================================
# 5. 配色：MGF=渐变 / core=黑色 / variable=灰色
# ============================================================================
target_groups <- c("100", "110", "300", "360", "505")
group_color <- c("100"="#C62828","110"="#3949AB","300"="#F9A825","360"="#6D4C41","505"="#00897B")
#group_light <- c("100"="#B2DFDB","110"="#C5CAE9","300"="#FFF9C4","360"="#D7CCC8","505"="#CFD8DC")
group_light <- c("100"="#FFEBEE","110"="#E8EAF6","300"="#FFFDE7","360"="#EFEBE9","505"="#E0F2F1")
# 
df_plot$cluster <- gsub("_.*", "", df_plot$GeneName)
df_plot$subcluster <- as.integer(gsub(".*_c", "", df_plot$GeneName))

df_plot <- df_plot %>%
  group_by(cluster) %>%
  mutate(
    n_sub = max(subcluster, na.rm = TRUE),
    frac = ifelse(n_sub > 1 & !is.na(n_sub), (subcluster - 1)/(n_sub - 1), 0.5)
  ) %>%
  rowwise() %>%
  mutate(
    fill_color = case_when(
      cluster %in% target_groups ~
        rgb(
          (1-frac)*col2rgb(group_light[cluster])[1] + frac*col2rgb(group_color[cluster])[1],
          (1-frac)*col2rgb(group_light[cluster])[2] + frac*col2rgb(group_color[cluster])[2],
          (1-frac)*col2rgb(group_light[cluster])[3] + frac*col2rgb(group_color[cluster])[3],
          maxColorValue = 255
        ),
      !cluster %in% target_groups & core_type == "core" ~ "gray70",
      !cluster %in% target_groups & core_type == "variable" ~ "black",
      TRUE ~ "gray70"
    )
  ) %>%
  ungroup()

# # ============================================================================
# # 6. ✅ 最终画图：原始位置 + 正确方向箭头
# # ============================================================================
# ggplot(df_plot,
#        aes(xmin = Start,
#            xmax = End,
#            y = GnmID,
#            fill = fill_color)) +
#   
#   geom_gene_arrow(
#     arrowhead_height = unit(7, "mm"),
#     arrowhead_width = unit(5, "mm"),
#     color = NA,
#     linewidth = 0
#   ) +
#   
#   scale_fill_identity() +
#   facet_wrap(~genotype, ncol = 1, scales = "free_y") +
#   labs(x="", y="Genome") +
#   theme_void() +
#   theme(
#     axis.text.y = element_text(size=11, face="bold", color="black"),
#     strip.text = element_text(size=12, face="bold"),
#     panel.background = element_rect(fill="transparent", color=NA),
#     plot.background = element_rect(fill="transparent", color=NA),
#     panel.spacing = unit(0.5, "lines")
#   )
# 
# 
# 
# 




## ============================================================================
# ✅ 只保留 MGF + 保留真实长度 + 紧凑排列 + 171000处空500
# ============================================================================
df_plot <- df_plot %>% filter(Group %in% c(100,110,300,360,505))

df_plot <- df_plot %>%
  group_by(genotype, GnmID) %>%
  arrange(Start) %>%
  mutate(
    gene_length = abs(Start - End),
    normal_gap = 30,
    after_171000 = Start >= 171000,
    gap = ifelse(lag(after_171000, default=F)==F & after_171000==T, 2000, normal_gap),
    gap = replace_na(gap, 0),
    offset = cumsum(lag(gene_length, default=0) + gap),
    x_start = offset,
    x_end = offset + gene_length,
    xmin_fit = ifelse(strand == "+", x_start, x_end),
    xmax_fit = ifelse(strand == "+", x_end, x_start),
    label_x = (xmin_fit + xmax_fit)/2,  # 基因名字居中位置
    
    # ==========================
    # ✅ 新增：符合条件加 (*)
    # ==========================
    is_abnormal = Frameshift == 1 | Inframe_stp == 1 | Start_codon == 0,
    gene2 = ifelse(is_abnormal, paste0(GeneName, "*"), GeneName)
  ) %>%
  ungroup()

# ============================================================================
# ✅ 最终画图：箭头 + 基因名字（带*）
# ============================================================================
df_plot=data.table(df_plot)
df_plot[GeneName == "100_c1"]$GeneName="100_1L2L"
df_plot[GeneName == "100_c2"]$GeneName="100_1R"
df_plot[GeneName == "100_c3"]$GeneName="100_3L"
df_plot[GeneName == "110_c9"]$GeneName="110_8L"
df_plot[GeneName == "110_c5"]$GeneName="110_1L"
df_plot[GeneName == "300_c1"]$GeneName="300_1L"
df_plot[GeneName == "300_c2"]$GeneName="300_2R"
df_plot[GeneName == "300_c3"]$GeneName="300_3L4L"
df_plot[GeneName == "360_c8"]$GeneName="360_19Rb"
df_plot[GeneName == "360_c9"]$GeneName="360_14L"
df_plot[GeneName == "360_c10"]$GeneName="360_15R"
df_plot[GeneName == "360_c3"]$GeneName="360_16R"
df_plot[GeneName == "360_c5"]$GeneName="360_18R"
df_plot[GeneName == "505_c1"]$GeneName="505_1R"
df_plot[GeneName == "505_c6"]$GeneName="505_3R"
df_plot[GeneName == "505_c4"]$GeneName="505_4R"
df_plot[GeneName == "505_c3"]$GeneName="505_11L"
df_plot[GeneName == "505_c5"]$GeneName="505_6R7R8R"
df_plot[GeneName == "505_c2"]$GeneName="505_2R5R9R10R"

df_plot$gene2=gsub("\\w+_\\w+","",df_plot$gene2)
ggplot(df_plot,
       aes(xmin = xmin_fit,
           xmax = xmax_fit,
           y = GnmID,
           fill = fill_color)) +
  
  # 基因箭头
  geom_gene_arrow(
    arrowhead_height = unit(7, "mm"),
    arrowhead_width = unit(5, "mm"),
    color = NA,
    linewidth = 0
  ) +
  
  # ✅ 显示带(*)的基因名
  geom_text(
    aes(x = label_x, label = gene2),  # 这里改成 gene2
    color = "black",
    size = 3.5,
    fontface = "bold",
    check_overlap = TRUE,
    vjust = 0.5,
    hjust = -0.5,
    angle=90
  ) +

  scale_fill_identity() +
  facet_wrap(~genotype, ncol = 1, scales = "free_y") +
  labs(x="", y="Genome") +
  theme_void() +
  theme(
    axis.text.y = element_text(size=11, face="bold", color="black"),
    strip.text = element_text(size=12, face="bold"),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.spacing = unit(0.8, "lines")
  )

ggsave(filename = "FigE_MGF_genome_arrow.pdf", 
       width = 16, height = 20, device = "pdf", bg = "white")

