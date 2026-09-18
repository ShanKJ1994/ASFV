# #############################################
# 
# #codon usage
# 
# #############################################
# 
# # ======================================
# # 1. 读入 featureCounts
# # ======================================
# fc <- read.delim("featureCounts_S0.txt", comment="#", sep="\t", stringsAsFactors=F)
# 
# # ======================================
# # 2. 只保留 病毒基因（去掉宿主 ENSSSCG 开头）
# # ======================================
# virus <- fc %>% filter(!grepl("^ENSSSCG", Geneid))
# 
# # 拆分：基因名、长度、counts
# gene_info <- virus %>% select(Geneid, Length)
# counts   <- virus %>% select(-c(Chr,Start,End,Strand,Length))
# # ======================================
# # 2. 只保留病毒基因（去掉宿主）
# # ======================================
# virus <- fc %>% filter(!grepl("^ENSSSCG", Geneid))
# 
# gene_names <- virus$Geneid
# gene_len   <- virus$Length
# count_mat  <- virus %>% select(-c(Chr, Start, End, Strand, Length, Geneid))
# 
# # ======================================
# # ✅ 绝对正确 tpm 函数（不会产生任何 NA/Inf）
# # ======================================
# # my_tpm <- function(counts, gene_length) {
# #   total <- colSums(counts)  # 每个样本总reads
# #   rpm <- t(t(counts) / total) * 1e6
# #   tpm <- rpm / gene_length
# #   return(as.matrix(tpm))
# # }
# 
# # TPM（顺序：先长度校正 → 再样本归一）
# my_tpm <- function(counts, gene_length) {
#   rpk <- counts / (gene_length / 1000)
#   total_rpk <- colSums(rpk)
#   tpm <- t(t(rpk) / total_rpk) * 1e6
#   return(as.matrix(tpm))
# }
# 
# # 计算 tpm
# tpm_mat <- my_tpm(count_mat, gene_len)
# rownames(tpm_mat) <- gene_names
# 
# # ======================================
# # 合并重复，取平均
# # ======================================
# groups <- sub("_\\d+$", "", colnames(tpm_mat))
# 
# tpm_avg <- do.call(cbind, lapply(unique(groups), function(g) {
#   rowMeans(tpm_mat[, groups == g, drop = FALSE])
# }))
# colnames(tpm_avg) <- unique(groups)
# tpm_avg<-data.frame(tpm_avg[,-c(1,5:7)])
# tpm_avg$gene=rownames(tpm_avg)
# write.csv(
#   tpm_avg,
#   "PAM_ASFV_gene_Class_RNAseq.csv",
#   row.names = F,
#   fileEncoding = "UTF-8"
# )

# ======================================
# 按表达量分类 Early / Late
# ======================================
setwd("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/RNAseq/")

tpm_avg=read.csv("PAM_ASFV_gene_Class_RNAseq.csv")
rownames(tpm_avg)=tpm_avg$gene
class_df <- tpm_avg[,-4]
head(tpm_avg)

# ======================================
# 1. 你的表达矩阵（基因 × 时间点）
# ======================================
#tpm_data <- tpm_avg
rownames(class_df)=class_df$gene
mat_exp <- class_df %>% select(ASFV_06H, ASFV_12H, ASFV_24H)

# 计算 Z-score (scale = "row" 和热图逻辑完全一样)
mat_z <- t(scale(t(mat_exp))) 


library(Mfuzz)

# 将数据框转换为矩阵格式
mat_exp_matrix <- as.matrix(mat_z )

# 使用转换后的矩阵创建 ExpressionSet 对象
eset <- new("ExpressionSet", exprs = mat_z)
# ================= 2. 数据标准化 =================
# Mfuzz 要求数据必须经过标准化处理（通常使用 z-score）
eset_std <- eset

# ================= 3. 寻找最佳聚类数 (c) 和模糊参数 (m) =================
# 这一步非常重要，用于确定最合适的簇数
# Dmin 默认设为 1，m 的搜索范围通常在 1.1 到 2.0 之间
c_opt <- Dmin(eset_std, c = 2:10, m = seq(1.1, 2.0, 0.1), 
              repeats = 100, visu = TRUE)
# 查看并选择最佳的 c 和 m 组合（通常选择曲线趋于平缓的拐点）
# 假设通过图表我们选择 c = 3, m = 1.25 (请根据实际生成的图调整)
c_opt
best_c <- 4
best_m <- 0.6951314

# ================= 4. 执行 Mfuzz 软聚类 =================
cl <- mfuzz(eset_std, c = best_c, m = best_m)

# 查看聚类结果
print(cl$cluster)       # 每个基因所属的主要簇
print(cl$membership)    # 每个基因属于各个簇的概率（软聚类核心）

# ================= 5. 结果可视化 =================
# 绘制聚类中心及包含的基因表达趋势
# col = rainbow(best_c) 为不同簇分配颜色
#pdf("Mfuzz.pdf", width = 8, height = 8)
mfuzz.plot(
  eset_std, 
  cl = cl, 
  mfrow = c(3, 2), 
  time.labels = c("06H", "12H", "24H"),
  colo = rainbow(best_c),  # ✅ 注意：是colo，不是col！
  new.window = FALSE
)
#dev.off()

# # 绘制软聚类成员关系图（查看哪些基因同时属于多个簇）
# pdf("Mfuzz_Membership_Heatmap.pdf", width = 8, height = 6)
# pheatmap(
#   cl$membership,
#   cluster_rows = FALSE,
#   cluster_cols = FALSE,
#   color = colorRampPalette(c("blue","white","#D82828"))(100),
#   show_rownames = FALSE,
#   main = "Membership matrix"
# )
# dev.off()

library(pheatmap)

#==================== 1. 提取mfuzz聚类信息 ====================
gene_cluster <- cl$cluster
df_mfuzz <- data.frame(
  gene_id = names(gene_cluster),
  mfuzz_cluster = factor(gene_cluster, levels = sort(unique(gene_cluster)))
)
rownames(df_mfuzz) <- df_mfuzz$gene_id
unique(cl$cluster)
#==================== 2. 合并原有注释（Early/Other/Late） ====================
# 不需要Early/Other/Late就注释merge两行，启用下面备选
# head(class_df)
# mat_z=data.frame(mat_z)
# mat_z$gene=rownames(tpm_avg)
# anno_row_combine <- merge(mat_z,mid,by=c("gene"))

# 2. 核心：class转为factor，强制顺序 Early > Other > Late
df_z <- anno_row_combine %>% 
  mutate(mfuzz_cluster = factor(mfuzz_cluster, levels = c(1,2,4,3))) %>%
  # ✅ 按你要求的优先级依次排序
  arrange(
    mfuzz_cluster,          # 第1优先级：分组顺序 Early > Other > Late
    desc(ASFV_06H),  
    desc(ASFV_12H),  
    desc(ASFV_24H)   
  )
head(df_z )


df_z$Row.names=gsub("MGF_","MGF-",df_z$Row.names)
write.csv(
  df_z,
  "PAM_ASFV_gene_Class_heatmap.csv",
  row.names = F,
  fileEncoding = "UTF-8"
)

df_z=fread(  "PAM_ASFV_gene_Class_heatmap.csv",header = T)
# 3. 还原成【用于pheatmap的z-score矩阵】（行名=基因，列=时间点）
mat_z_ordered <- df_z %>% 
  select(Row.names, ASFV_06H, ASFV_12H, ASFV_24H) %>% 
  tibble::column_to_rownames("Row.names") %>% 
  as.matrix()

# 4. 同步得到排序后匹配的anno_row（给pheatmap的annotation_row用！）
anno_row <- df_z %>% 
  select(Row.names, mfuzz_cluster) %>% 
  tibble::column_to_rownames("Row.names")


pheatmap::pheatmap(
  mat_z_ordered,
  annotation_row = anno_row,
  #annotation_colors = my_annotation_colors,  # ✅ 应用颜色
  #scale = "row",
  cluster_rows = F,   # 按你的顺序
  cluster_cols = FALSE,
  show_rownames = T,
  color = colorRampPalette(c("navy","white","firebrick3"))(100),angle_col = 90,
  #filename = 'C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/RNAseq/Muffz_PAM_heatmap.pdf',
  width = 6, height= 12
)
dev.off()

mid=merge(filter(final,class!="Host"),anno_row_combine[,c(1,5)],by.x="gene",by.y = "Row.names")
mid$mfuzz_cluster2="Late"
mid[mid$mfuzz_cluster==2,]$mfuzz_cluster2="Early"
mid[mid$mfuzz_cluster%in%c(1),]$mfuzz_cluster2="Other"
table(mid$mfuzz_cluster2)

my_annotation_colors <- list(
  class = c(
    Late = "#FFCCCC",    # 淡红
    Early  = "#CCCCFF",    # 淡蓝
    Other = "#EEEEEE"     # 淡灰
  )
)
mid=fread("PAM_CAI_RNAseq_mock_Class.csv")
ggplot(mid, aes(x = mfuzz_cluster2, y = CAI, fill = mfuzz_cluster2)) +
  geom_boxplot(outlier.shape = NA, width = 0.7) +  # 隐藏重复点
  geom_jitter(mid, mapping=aes(x = mfuzz_cluster2, y = CAI, fill = mfuzz_cluster2),alpha = 0.4, size = 1.5, width = 0.1) +
  scale_fill_manual(values = my_colors) +  # ✅ 应用颜色
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  ) +
  # 三组之间全部 Wilcoxon 比较
  ggpubr::stat_compare_means(
    method = "wilcox.test",
    comparisons = list(
      #c("Early", "Host"),
      c("Early", "Late"),
      c("Early", "Other"),
      c("Other", "Late")
    ),
    label = "p.value"  # 用显著性星号; 想改成p值用 label="p.value"
  ) + ylim(0.5,0.8)+
  #ggtitle("CAI by Early/Late (tpm-based classification)") +
  xlab("") + ylab("CAI")
ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/RNAseq/Muffz_PAM_CAI.pdf', width= 3 , height= 3, units='in')

##########################################################
final=fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/RNAseq/pig_pichong_CAI.csv")
table(final[,4:5])
head(final)
# --------------------
# 热图同款淡色（完全一致）
my_colors <- c(
  Host  = "#C5DEBA",
  Early = "#CCCCFF",    # 淡红
  Other = "#EEEEEE",    # 淡灰
  Late  =  "#FFCCCC"    # 淡蓝
)
ggplot(filter(final,class!="Host"), aes(x = class, y = CAI, fill = class)) +
  geom_boxplot(outlier.shape = NA, width = 0.7) +  # 隐藏重复点
  geom_jitter(filter(final,class!="Host"), mapping=aes(x = class, y = CAI, fill = class),alpha = 0.4, size = 1.5, width = 0.1) +
  scale_fill_manual(values = my_colors) +  # ✅ 应用颜色
  theme_classic() +
  theme(
    axis.text.x = element_text(),
    legend.position = "none",
    strip.background = element_blank() # ✅ 删掉分面顶部的方框
  )+facet_wrap(~type)+
  # 三组之间全部 Wilcoxon 比较
  ggpubr::stat_compare_means(
    method = "wilcox.test",
    comparisons = list(
      #c("Early", "Host"),
      c("Early", "Late"),
      c("Early", "Other"),
      c("Other", "Late")
    ),
    label = "p.value"  # 用显著性星号; 想改成p值用 label="p.value"
  ) + #ylim(0.5,0.9)+
  #ggtitle("CAI by Early/Late (tpm-based classification)") +
  xlab("") + ylab("CAI")

table(final$class)

ggplot(filter(final,mfuzz_cluster2!="Host"), aes(x = mfuzz_cluster2, y = CAI, fill = mfuzz_cluster2)) +
  geom_violin()+
  geom_boxplot(outlier.shape = NA, width = 0.3) +  # 隐藏重复点
  #geom_jitter(filter(final,mfuzz_cluster2!="Host"), mapping=aes(x = mfuzz_cluster2, y = CAI, fill = mfuzz_cluster2),alpha = 0.4, size = 1.5, width = 0.1) +
  scale_fill_manual(values = my_colors) +  # ✅ 应用颜色
  theme_classic() +
  theme(
    axis.text.x = element_text(),
    legend.position = "none",
    strip.background = element_blank() # ✅ 删掉分面顶部的方框
  )+facet_wrap(~type,scale="free")+
  # 三组之间全部 Wilcoxon 比较
  ggpubr::stat_compare_means(
    method = "wilcox.test",
    comparisons = list(
      #c("Early", "Host"),
      c("Early", "Late"),
      c("Early", "Other"),
      c("Other", "Late")
    ),
    label = "p.value"  # 用显著性星号; 想改成p值用 label="p.value"
  ) + #ylim(0.45,0.8)+
  #ggtitle("CAI by Early/Late (tpm-based classification)") +
  xlab("") + ylab("CAI")
ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/RNAseq/Fig_CAI_Class.pdf', width= 5 , height= 3, units='in')

final=fread("pig_pichong_CAI.csv")
head(final)
final$virus="Virus"
final[grepl("ost",final$class),]$virus="Host"
ggplot(final, aes(x = virus, y = CAI, fill = virus)) +
  geom_violin() + 
  geom_boxplot(outlier.shape = NA, width = 0.35) +  # 隐藏重复点
  scale_fill_manual(values = c("#1F78B4","#FF7F00")) +  # ✅ 应用颜色
  theme_classic() +
  theme(
    axis.text.x = element_text(),
    legend.position = "none",
    strip.background = element_blank() # ✅ 删掉分面顶部的方框
  )+facet_wrap(~type)+
  # 三组之间全部 Wilcoxon 比较
  # ggpubr::stat_compare_means(
  #   method = "wilcox.test",
  #   comparisons = list(
  #     c("Virus", "Host")
  #   ),
  #   label = "p.value"  # 用显著性星号; 想改成p值用 label="p.value"
  # ) + #ylim(0.5,1)+
  #ggtitle("CAI by Early/Late (tpm-based classification)") +
  xlab("") + ylab("CAI")
#final %>% filter(type == "Ornithodoros") %>% wilcox.test(CAI ~ virus, data = .)
#final %>% filter(type == "Pig") %>% wilcox.test(CAI ~ virus, data = .)
ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/RNAseq/Fig_total_CAI.pdf', width= 4 , height= 3, units='in')


filter(final,class=="Host") %>% group_by(type) %>% dplyr::summarise(Host_median=median(CAI))

pichong=filter(final,type=="Ornithodoros",class!="Host")[,c(1,2,4)]
colnames(pichong)[2]="Ornithodoros"
pig=filter(final,type=="Pig",class!="Host")[,c(1,2,4)]
colnames(pig)[2]="Pig"
final=merge(pig,pichong,by=c("gene","mfuzz_cluster2"))
final$diff=final$Pig-final$Ornithodoros
final=final[order(final$diff,decreasing = T),]
head(final)


final$CAI_pig_minus_virus=0.762-final$Pig
final$CAI_pichong_minus_virus=0.730-final$Ornithodoros


# 假设你的数据框叫 final
long_data <- final[,c(1,5,6)] %>%
  # 宽变长：保留 gene, class，把 CAI_pig... 和 CAI_pichong... 变成两列
  pivot_longer(
    cols = c(CAI_pig_minus_virus, CAI_pichong_minus_virus),
    names_to = "host_type",
    values_to = "CAI_diff"
  ) %>%
  # 把列名改成更干净的名字
  mutate(
    host_type = case_when(
      host_type == "CAI_pig_minus_virus" ~ "Pig_VS_Virus",
      host_type == "CAI_pichong_minus_virus" ~ "Ornithodoros_VS_Virus"
    )
  )

# 查看结果
head(long_data)
ggplot(long_data, aes(x = host_type, y = CAI_diff, fill = host_type)) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA, width = 0.4) +  # 隐藏重复点
  #geom_jitter(alpha = 0.4, size = 1.5, width = 0.1) +
  scale_fill_manual(values = c("#6A3D9A","#B15928")) +  # ✅ 应用颜色
  theme_classic() +
  #facet_wrap(~type,scales = "free")+
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  ) +
  # 三组之间全部 Wilcoxon 比较
  ggpubr::stat_compare_means(
    method = "wilcox.test",
    comparisons = list(
      c("Pig_VS_Virus", "Ornithodoros_VS_Virus")
    ),
    label = "p.value"  # 用显著性星号; 想改成p值用 label="p.value"
  ) + #ylim(0.5,0.9)+
  #ggtitle("CAI by Early/Late (tpm-based classification)") +
  xlab("") + ylab("Distance (host median CAI - virus CAI)")

ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/RNAseq/Fig_CAI_diff.pdf', width= 2 , height= 3, units='in')









df=fread( "C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/RNAseq/CAI_genome.tsv")
#genotype_df=fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/miniprot/meta.tsv")
#df <- merge(df,genotype_df,by.x = c("GenomeID"),by.y = c("name"))
head(df)
head(genotype_df)
table(df$genotype)

test=filter(df,genotype %in% c(1,2,9))[,c(1,2,5)] 
colnames(test)[2]="CAI"
test$Type="Ornithodoros"

test2=filter(df,genotype %in% c(1,2,9))[,c(1,3,5)] 
colnames(test2)[2]="CAI"
test2$Type="Pig"

ggplot(rbind(test,test2) ,aes(x=genotype,y=CAI,fill=genotype))+
  #geom_violin()+
  facet_wrap(~Type,scale="free")+
  geom_boxplot(outliers = F,width=0.3)+
  geom_jitter(mapping=aes(color=genotype),alpha=0.2,size = 0.5)+
  scale_color_manual(values=c("#1F78B4","#E31A1C","#33A02C","#B15928"))+
  scale_fill_manual(values=c("#1F78B4","#E31A1C","#33A02C","#B15928"))+
  theme_classic()+
  theme(
    strip.background = element_blank(),legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  ) +
#ylim(c(0.617,0.622))#+
ggpubr::stat_compare_means(method = 'wilcox.test', comparisons = list(c("2","9"),c("1","9"),c("1","2")),label = 'p.signif')
ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/RNAseq/Fig_CAI_genotype.pdf', width= 6 , height= 3, units='in')
# 
# ggplot(filter(df,genotype %in% c(1,2,9)) ,aes(x=genotype,y=CAI_PAMs,fill=genotype))+
#   geom_boxplot(outliers = F,width=0.3)+
#   geom_jitter(mapping=aes(color=genotype),alpha=0.2)+
#   scale_color_manual(values=c("#1F78B4","#E31A1C","#33A02C","#B15928"))+
#   scale_fill_manual(values=c("#1F78B4","#E31A1C","#33A02C","#B15928"))+
#   theme_classic()+
#   theme(
#     legend.position = "none",
#     plot.title = element_text(hjust = 0.5)
#   ) +
#   ylim(c(0.659,0.664))#+
# #ggpubr::stat_compare_means(method = 'wilcox.test', comparisons = list(c("1/2","9"),c("2","9"),c("1","9"),c("1","2"),c("1","1/2"),c("1/2","2")),label = 'p.signif')
# ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/RNAseq/Fig_CAI_PAM_class.pdf', width= 2 , height= 3, units='in')

ggplot(df ,aes(x=CAI_Ornithodoros,y=CAI_PAMs,color=genotype))+
  geom_point()+
  ggpubr::stat_cor(method = "spearman",col="red",size=5,cor.coef.name = "rho")+
  geom_smooth(method = "glm",se=F,col="black")+
  scale_color_manual(values=c("#1F78B4","#E31A1C","#33A02C","black","#6A3D9A","#FF7F00","#B15928","#AAAAAA"))+
  theme_classic()+
  # 下面这两行就是去掉 facet 顶部方框
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(size = 11),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  ) 
ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/RNAseq/Fig_CAI_PAM_Ornithodoros_cor.pdf', width= 3 , height= 3, units='in')

mid=rbind(filter(df,genotype==1,CodonNum>=55000,CodonNum<=56000),
          filter(df,genotype==2,CodonNum>=57000),
          filter(df,genotype==9,CodonNum>=54000,CodonNum<=56000))
library(MASS)   
library(robustbase)
ggplot(na.omit(filter(mid,genotype %in% c(1,2,9,"1/2"),years>0)) ,aes(x=as.numeric(years),y=CAI_PAMs,color=genotype))+
  geom_point()+
  facet_wrap(~genotype,scales = "free")+
  ggpubr::stat_cor(method = "spearman",col="red",size=5,cor.coef.name = "rho")+
  geom_smooth(method = "lmrob",se=F,col="black")+
  scale_color_manual(values=c("#1F78B4","#E31A1C","#33A02C","black","#6A3D9A","#FF7F00","#B15928","#AAAAAA"))+
  theme_classic()+
  # 下面这两行就是去掉 facet 顶部方框
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(size = 11)
  )
ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/RNAseq/Fig_CAI_PAM_year.pdf', width= 5 , height= 3, units='in')

ggplot(filter(mid,genotype %in% c(1,2,9,"1/2"),years>0) ,aes(x=as.numeric(years),y=CAI_Ornithodoros,color=genotype))+
  geom_point()+
  facet_wrap(~genotype,scales = "free")+
  ggpubr::stat_cor(method = "spearman",col="red",size=5,cor.coef.name = "rho")+
  geom_smooth(method = "lmrob",se=F,col="black")+
  scale_color_manual(values=c("#1F78B4","#E31A1C","#33A02C","black","#6A3D9A","#FF7F00","#B15928","#AAAAAA"))+
  theme_classic()+
  # 下面这两行就是去掉 facet 顶部方框
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(size = 11)
  )
ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/RNAseq/Fig_CAI_Ornithodoros_year.pdf', width= 5 , height= 3, units='in')

ggplot(mid ,aes(x=CodonNum,fill=genotype))+
  geom_histogram(binwidth=100)+
  facet_wrap(~genotype,scales = "free")+
  scale_color_manual(values=c("#1F78B4","#E31A1C","#33A02C","black","#6A3D9A","#FF7F00","#B15928","#AAAAAA"))+
  theme_classic()
table(mid$genotype)

ggplot(filter(mid,genotype %in% c(9),years>0) ,aes(x=as.numeric(years),y=CAI_Ornithodoros,color=genotype))+
  geom_point()+
  #facet_wrap(~genotype,scales = "free")+
  ggpubr::stat_cor(method = "spearman",col="red",size=5,cor.coef.name = "rho")+
  geom_smooth(method = "lmrob",se=F,col="black")+
  scale_color_manual(values=c("#1F78B4","#E31A1C","#33A02C","black","#6A3D9A","#FF7F00","#B15928","#AAAAAA"))+
  theme_classic()+
  # 下面这两行就是去掉 facet 顶部方框
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(size = 11)
  )+
  ylim(0.6199,0.62122)


df=fread( "C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/RNAseq/each_gnm_gene_CAI.txt")
genotype_df=read.table("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/miniprot/meta.tsv",header = T,sep="\t")
head(genotype_df)
df=merge(df,genotype_df,by.x="GnmID",by.y="name")
head(df)

#nrow(filter(ASFV,genotype %in% c(1,2,9)))
all_gnm <- unique(filter(df,genotype %in% c(1,2,9))$GnmID)

mid_core <- df %>%
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



table(genotype_df$genotype)
#test=filter(mid,unique_GnmID_count== gene_count,loss_GnmID_count==0 )
df_core <- filter(df,GeneName %in% filter(mid_core,unique_GnmID_count== gene_count,loss_GnmID_count/357<=0.05)$GeneName)
df_core_median <- df_core %>% group_by(GeneName,genotype) %>% dplyr::summarise(CAI_PAMs_median=median(CAI_PAMs),CAI_Ornithodoros_median=median(CAI_Ornithodoros))

head(df_core_median)
test=filter(df_core_median,genotype %in% c(1,2,9))
mid1=test[,c(1,2,3)]
mid1$Type="Pig"
colnames(mid1)[3]="CAI"
mid3=test[,c(1,2,4)]
mid3$Type="Ornithodoros"
colnames(mid3)[3]="CAI"

ggplot(rbind(mid1,mid3) ,aes(x=genotype,y=CAI,fill=genotype))+
  geom_violin()+
  facet_wrap(~Type,scale="free")+
  geom_boxplot(outliers = F,width=0.3)+
  #scale_y_continuous(breaks=seq(0.45,0.75,0.1),limits = c(0.45,0.75))+
  #  geom_jitter(mapping=aes(color=genotype),alpha=0.2)+
  ggpubr::stat_compare_means(method = 'wilcox.test', #paired = T,
                             comparisons = list(c(1,2),c(9,2),c(1,9))
  )+ 
  theme_classic()+
  # 下面这两行就是去掉 facet 顶部方框
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text.x = element_text(size = 11)
  ) 
ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/RNAseq/Fig_CAI_genotype_core.pdf', width= 5 , height= 3, units='in')


# 
# 
# df <- rbind(mid1, mid3)
# value_col <- "CAI"  # 确认这里是你的真实列名
# 
# comp_list <- list(c("1","2"), c("9","2"), c("1","9"))
# res_all <- list()
# 
# # 遍历每一个Type
# for(one_type in unique(df$Type)){
# 
#   subdf <- df[df$Type == one_type, ] # 筛选同一个Type的数据
# 
#   for(pair in comp_list){
#     g1 <- pair[1]
#     g2 <- pair[2]
# 
#     v1 <- subdf[[value_col]][as.character(subdf$genotype) == g1]
#     v2 <- subdf[[value_col]][as.character(subdf$genotype) == g2]
# 
#     # 样本量校验，少于2个样本直接返回NA，防止代码崩溃
#     if(length(v1) <2 || length(v2) <2){
#       p_out <- NA_real_
#       note_out <- "样本不足"
#     }else{
#       wt <- wilcox.test(v1, v2)
#       # 如果需要配对检验：wt <- wilcox.test(v1, v2, paired = TRUE)
#       p_out <- wt$p.value
#       note_out <- ""
#     }
# 
#     res_all[[length(res_all)+1]] <- data.frame(
#       Type = one_type,
#       group1 = g1,
#       group2 = g2,
#       p.value = p_out,
#       method = "wilcox.test",
#       note = note_out,
#       stringsAsFactors = FALSE
#     )
#   }
# }
# 
# # 合并结果
# stat_table <- do.call(rbind, res_all)
# stat_table
# #

df_noncore <- filter(df,!GeneName %in% df_core$GeneName)

mid=df_noncore %>% group_by(genotype,GnmID) %>% dplyr::summarise(count=n())
ggplot(filter(mid,genotype %in% c(1,2,9)) ,aes(x=genotype,y=count,fill=genotype))+
  #geom_violin()+
  #facet_wrap(~Type,scale="free")+
  geom_boxplot(outliers = F,width=0.3)+
  #scale_y_continuous(breaks=seq(0.45,0.75,0.1),limits = c(0.45,0.75))+
  geom_jitter(mapping=aes(color=genotype),alpha=0.2)+
  ggpubr::stat_compare_means(method = 'wilcox.test', #paired = T,
                             comparisons = list(c(1,2),c(9,2),c(1,9))
  )+ 
  theme_classic()+
  # 下面这两行就是去掉 facet 顶部方框
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text.x = element_text(size = 11)
  ) +
  #ylim(20,50)+
  xlab("") + ylab("noncore gene number")
ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/RNAseq/Fig_CAI_genotype_noncore_number.pdf', width= 2 , height= 3, units='in')

genotype_df=fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/miniprot/meta.tsv")
mid=merge(mid,filter(genotype_df,genotype %in% c(1,2,9))[,c(1,3)],by.x="GnmID",by.y="name")
ggplot(mid ,aes(x=log10(Length),y=count,col=genotype))+
  #geom_violin()+
  #facet_wrap(~Type,scale="free")+
  geom_point(mid ,mapping=aes(x=log10(Length),y=count))+
  ggpubr::stat_cor(method = "spearman",col="red",size=5,cor.coef.name = "rho")+
  geom_smooth(method = "lmrob",se=F,col="black")+
  theme_classic()+
  # 下面这两行就是去掉 facet 顶部方框
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text.x = element_text(size = 11)
  ) +
  #ylim(20,50)+
  xlab("genome length (log10)") + ylab("noncore gene number")
ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/RNAseq/Fig_CAI_genotype_noncore_number.pdf', width= 2 , height= 3, units='in')



df_core_median <- df_core %>% group_by(GeneName) %>% dplyr::summarise(CAI_PAMs_median=median(CAI_PAMs),CAI_Ornithodoros_median=median(CAI_Ornithodoros))
df_core_median$Core="core"

df_noncore <- filter(df,!GeneName %in% df_core $GeneName)
df_noncore_median <- df_noncore %>% group_by(GeneName) %>% dplyr::summarise(CAI_PAMs_median=median(CAI_PAMs),CAI_Ornithodoros_median=median(CAI_Ornithodoros))
df_noncore_median$Core="noncore"
head(df_noncore_median)



mid2=rbind(df_noncore_median,df_core_median)
mid1=mid2[,c(1,2,4)]
mid1$Type="Pig"
colnames(mid1)[2]="CAI"
mid3=mid2[,c(1,3,4)]
mid3$Type="Ornithodoros"
colnames(mid3)[2]="CAI"

ggplot(rbind(mid1,mid3) ,aes(x=Core,y=CAI,fill=Core))+
  geom_violin()+
  facet_wrap(~Type,scale="free")+
  geom_boxplot(outliers = F,width=0.3)+
  scale_y_continuous(breaks=seq(0.45,0.75,0.1),limits = c(0.45,0.75))+
  #  geom_jitter(mapping=aes(color=genotype),alpha=0.2)+
  ggpubr::stat_compare_means(method = 'wilcox.test', #paired = T,
                             comparisons = list(c("core","noncore"))
  )+ 
  theme_classic()+
  # 下面这两行就是去掉 facet 顶部方框
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text.x = element_text(size = 11)
  ) 

ggsave('C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/RNAseq/Fig_CAI_Core.pdf', width= 4 , height= 3, units='in')



