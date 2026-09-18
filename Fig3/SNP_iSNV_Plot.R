#############################################

#iSNV

#############################################
library(Biostrings)
setwd("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/SNP/")

#SNP<-fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/SNP/Common_snp_call_result.tsv")

SNP<-fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/SNP/snp_call_result.tsv")

head(SNP)
SNP$SNP=paste0(SNP$ancestor_base,">",SNP$sample_base)
count=SNP %>% filter(SNP %in% c("C>T","G>A","A>G","T>C","G>T","C>A","T>G","A>C","G>C","C>G","A>T","T>A")) %>% group_by(ancestor_id,aln_pos,genotype,SNP,motif) %>% dplyr::summarise(count=n()) 
snv_map <- tibble(
  SNP =  c("C>T","G>A",
           "A>G","T>C",
           "G>T","C>A",
           "T>G","A>C",
           "G>C","C>G",
           "A>T","T>A"),
  mut = c( "C:G>T:A","C:G>T:A",
           "A:T>G:C","A:T>G:C",
           "C:G>A:T","C:G>A:T",
           "A:T>C:G","A:T>C:G",
           "C:G>G:C","C:G>G:C",
           "A:T>T:A","A:T>T:A")
)

count=merge(count,snv_map,by=c("SNP"))

#=====================
# 2. 标准化函数：C>T不变；G>A反向互补；其他丢弃
#=====================
std_motif <- function(snp_type, seq){
  if(snp_type == "C>T"){
    return(seq)
  } else if(snp_type == "G>A"){
    return(as.character(reverseComplement(DNAString(seq))))
  } else {
    return(NA_character_)
  }
}


#=====================
# 3. 筛选C:G>T:A组 + 批量标准化motif
#=====================
df_cgta <- count %>%
  filter(mut == "C:G>T:A") %>%
  rowwise() %>%
  mutate(motif_std = std_motif(SNP, motif)) %>%
  ungroup() %>%
  filter(!is.na(motif_std)) # 剔除不是C>T/G>A的脏数据


df_tc <- filter(df_cgta,mut=="C:G>T:A") %>%
  mutate(
    pos5 = substr(motif_std, 5, 5),
    pos6 = substr(motif_std, 6, 6),
    dinuc = paste0(pos5, pos6),
    is_TC = ifelse(dinuc == "TC" , 1, 0)
  )

df_tc  %>%  group_by(genotype) %>% dplyr::summarise(count=n())
#==== 2. 按基因型计算TC比例 ====
df_tc_ratio <- df_tc %>%
  group_by(genotype) %>%
  summarise(
    n_total = n(),
    n_TC = sum(is_TC),
    TC_ratio = mean(is_TC), # 等价 n_TC / n_total
    .groups = "drop"
  )

# 查看结果
print(df_tc_ratio)


# # Genotype1
bt1 <- binom.test(x = 10, n = 33, p = 0.25)
# Genotype2
bt2 <- binom.test(x = 21, n = 62, p = 0.25)
# Genotype9
bt9 <- binom.test(x = 4, n = 22, p = 0.25)

count %>% group_by(genotype) %>% dplyr::summarise(sum=sum(count))
# # Genotype1
bt1 <- binom.test(x = 10, n = 33, p = 0.25)
binom.test(x = 33, n = 37, p = 0.5)
binom.test(x = 37, n = 46, p = 2/6)

# Genotype2
bt2 <- binom.test(x = 21, n = 62, p = 0.25)
binom.test(x = 62, n = 62+34, p = 0.5)
binom.test(x = 62+34, n = 157, p = 2/6)

# Genotype9
bt9 <- binom.test(x = 4, n = 22, p = 0.25)
binom.test(x = 22, n = 22+5, p = 0.5)
binom.test(x = 22+5, n = 33, p = 2/6)

count

#==== 3. 绘制bar图 ====
ggplot(df_tc_ratio, aes(x = factor(genotype), y = TC_ratio)) +
  geom_col(width = 0.6, fill = "#4088B4") + # 沿用你固定配色
  geom_text(aes(label = sprintf("%.2f\n(%d/%d)", TC_ratio, n_TC, n_total)),
            vjust = -0.3, size = 3.8) +
  scale_y_continuous(limits = c(0,0.4), expand = c(0,0)) +
  labs(
    x = "Genotype",
    y = "Proportion of TC dinucleotide (pos5–6)",
    title = "TC motif proportion at position 5–6"
  ) +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5))

#=====================
# 2. 标准化函数：C>T不变；G>A反向互补；其他丢弃
#=====================
std_motif <- function(snp_type, seq){
  if(snp_type == "C>T"){
    return(seq)
  } else if(snp_type == "G>A"){
    return(as.character(reverseComplement(DNAString(seq))))
  } else {
    return(NA_character_)
  }
}


#=====================
# 3. 筛选C:G>T:A组 + 批量标准化motif
#=====================
df_cgta <- count %>%
  filter(mut == "C:G>T:A",nchar(motif)==11) %>%
  rowwise() %>%
  mutate(motif_std = std_motif(SNP, motif)) %>%
  ungroup() %>%
  filter(!is.na(motif_std)) # 剔除不是C>T/G>A的脏数据

# 提取对齐好的标准化序列，直接喂给seqlogo
library(patchwork)
# 按基因型拆分
g1 <- filter(df_cgta, genotype == 1)$motif_std
g2 <- filter(df_cgta, genotype == 2)$motif_std
g9 <- filter(df_cgta, genotype == 9)$motif_std

p1 <- ggseqlogo(g1,  method="prob") + ggtitle("Genotype 1") + theme_bw()
p2 <- ggseqlogo(g2,  method="prob") + ggtitle("Genotype 2") + theme_bw()
p3 <- ggseqlogo(g9,  method="prob") + ggtitle("Genotype 9") + theme_bw()

p_combine <- p1 + p2 + p3 + plot_annotation(title="C:G>T:A motif")
ggsave("FigB_SNP_seqlogo_CGTA_genotype_patchwork.pdf", p_combine, width=9, height=3)



count=count %>% group_by(genotype,mut) %>% dplyr::summarise(count=n())
count
Self1<-c("#E31A1C","#1F78B4","#B2DF8A","#FDBF6F","#6A3D9A","#B15928" )
count$mut=factor(count$mut,levels = c("C:G>T:A",
                                      "A:T>G:C",
                                      "C:G>A:T",
                                      "A:T>C:G",
                                      "C:G>G:C",
                                      "A:T>T:A"))
ggplot(count,aes(x=mut,y=count,fill=mut))+
  geom_bar(position = 'dodge',stat = 'identity')+
  scale_fill_manual(values = Self1)+
  facet_wrap(~genotype,scales = "free")+
  theme_classic()+
  labs(y="SNP number")+
  theme(axis.text.x = element_text(),#angle = 90, hjust = 1, vjust = 0.5
        legend.position = "none" ,#  移除图例
        strip.background = element_blank()) # ✅ 删掉分面顶部的方框 
ggsave("FigA_SNP_count_CGTA_genotype_patchwork.pdf",  width=9, height=3)


library(Biostrings)
library(rtracklayer)
library(ggseqlogo)
library(ggplot2)
library(dplyr)
library(data.table)

iSNV=fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/SNP/iSNVs.norm.txt",header = T)
head(iSNV)

iSNV=filter(iSNV,nchar(Ref)==1,nchar(Alt)==1,Alt_count>=5,Ref_count/Depth>=0.9)#
iSNV$Alt_Frac=iSNV$Alt_count/iSNV$Depth
#min(iSNV$Alt_Frac)
ggplot(iSNV,aes(x=Alt_Frac))+
  geom_histogram(binwidth = 0.01)+
  theme_classic()+
  labs(x="mutation fraction",y="SNV number")
# iSNV$bin="[0.05,0.25)"
# iSNV[Alt_Frac>=0.25 & Alt_Frac<0.5 ]$bin="[0.25,0.5)"
# iSNV[Alt_Frac>=0.5 & Alt_Frac<0.75 ]$bin="[0.5,0.75)"
# iSNV[Alt_Frac>=0.75 & Alt_Frac<0.9 ]$bin="[0.75,0.9)"
iSNV$mut=paste0(iSNV$Ref,">",iSNV$Alt)

snv_map <- tibble(
  mut =  c("C>T","G>A",
           "A>G","T>C",
           "G>T","C>A",
           "T>G","A>C",
           "G>C","C>G",
           "A>T","T>A"),
  SNP = c( "C:G>T:A","C:G>T:A",
           "A:T>G:C","A:T>G:C",
           "C:G>A:T","C:G>A:T",
           "A:T>C:G","A:T>C:G",
           "C:G>G:C","C:G>G:C",
           "A:T>T:A","A:T>T:A")
)

iSNV=merge(iSNV,snv_map,by=c("mut"))

iSNV_count=iSNV %>% group_by(mut) %>% dplyr::summarise(count=n())
Self1<-c("#1F78B4","#A6CEE3","#33A02C","#B2DF8A","#E31A1C","#FB9A99","#FF7F00","#FDBF6F","#6A3D9A","#CAB2D6","#B15928","#FFFF99" )
iSNV_count$mut=factor(iSNV_count$mut,levels = c("C>T","G>A","A>G","T>C","G>T","C>A","T>G","A>C","G>C","C>G","A>T","T>A"))
ggplot(iSNV_count,aes(x=mut,y=count,fill=mut))+
  geom_bar(position = 'dodge',stat = 'identity')+
  scale_fill_manual(values = Self1)+
  #facet_wrap(~bin,scales = "free")+
  theme_classic()+
  labs(y="iSNV number")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none"   )# 移除图例


iSNV_count2=iSNV %>% group_by(SNP) %>% dplyr::summarise(count=n())
Self1<-c("#E31A1C","#1F78B4","#B2DF8A","#FDBF6F","#6A3D9A","#B15928" )
iSNV_count2$SNP=factor(iSNV_count2$SNP,levels = c("C:G>T:A",
                                                  "A:T>G:C",
                                                  "C:G>A:T",
                                                  "A:T>C:G",
                                                  "C:G>G:C",
                                                  "A:T>T:A"))
ggplot(iSNV_count2,aes(x=SNP,y=count,fill=SNP))+
  geom_bar(position = 'dodge',stat = 'identity')+
  scale_fill_manual(values = Self1)+
  #facet_wrap(~bin,scales = "free")+
  theme_classic()+
  labs(y="iSNV number")+
  theme(axis.text.x = element_text(),#angle = 90, hjust = 1, vjust = 0.5
        legend.position = "none"   )# 移除图例

binom.test(177,sum(iSNV_count2$count),0.5)
ggsave("FigC_iSNV_count.pdf",  width=3, height=3)

# iSNV_count=filter(iSNV,Alt_Frac<=0.5) %>% group_by(Sample,mut) %>% dplyr::summarise(count=n())
# Self1<-c("#1F78B4","#A6CEE3","#33A02C","#B2DF8A","#E31A1C","#FB9A99","#FF7F00","#FDBF6F","#6A3D9A","#CAB2D6","#B15928","#FFFF99" )
# iSNV_count$mut=factor(iSNV_count$mut,levels = c("C>T","G>A","A>G","T>C","G>T","C>A","T>G","A>C","G>C","C>G","A>T","T>A"))
# ggplot(iSNV_count,aes(x=mut,y=count,fill=mut))+
#   geom_bar(position = 'dodge',stat = 'identity')+
#   scale_fill_manual(values = Self1)+
#   facet_wrap(~Sample,scales = "free")+
#   theme_bw()+
#   labs(y="SNV number")

library(Biostrings)
library(data.table)
library(dplyr)
library(ggseqlogo)

# 1. 读取参考序列
ref_genome <- readDNAStringSet("sequence.fasta")
ref_seq <- ref_genome[["NC_044959.2"]]
ref_len <- length(ref_seq) # 参考序列总长度

# 2. 导入SNP表格
iSNV=fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/SNP/iSNVs.norm.txt",header = T)
head(iSNV)

iSNV=filter(iSNV,nchar(Ref)==1,nchar(Alt)==1,Alt_count>=5,Ref_count/Depth>=0.9)#

snp_dt <- iSNV

# 定义提取函数：坐标合法才提取，越界返回NA
safe_subseq <- function(s, e, ref, maxL){
  if(s >=1 && e <= maxL && s <= e){
    as.character(subseq(ref, start = s, end = e))
  }else{
    NA_character_
  }
}

# 计算上下游区间
snp_dt[, start := Pos - 5]
snp_dt[, end   := Pos + 5]

# 逐行安全提取侧翼序列
snp_dt[, flank_seq := mapply(safe_subseq, start, end, MoreArgs = list(ref=ref_seq, maxL=ref_len))]

# ✅ 校验：第6位碱基是否等于Ref
snp_dt[!is.na(flank_seq), check_ref := substr(flank_seq, 6,6)]

# 过滤：保留非NA、严格11bp侧翼（上下游完整5nt）
snp_filter <- snp_dt[!is.na(flank_seq) & nchar(flank_seq) == 11]
snp_filter$mut=paste0(snp_filter$Ref,">",snp_filter$Alt)

# 构造DNAString，批量做反向互补
snp_filter[, rc_seq := as.character(reverseComplement(DNAStringSet(flank_seq)))]

# 按mut选择序列：G>A使用反向互补序列，其余使用原始序列
snp_filter[, motif_seq := ifelse(mut == "G>A", rc_seq, flank_seq)]

# ✅ 校验转换是否正确（可选，核对G>A组第6位是否变为C）
snp_filter[mut == "G>A", check_motif_ref := substr(motif_seq,6,6)]


# 查看有多少位点因为坐标越界被剔除
cat("原始位点总数：", nrow(snp_dt), "\n")
cat("有效11bp侧翼位点：", nrow(snp_filter), "\n")
cat("越界/不完整位点数量：", nrow(snp_dt) - nrow(snp_filter), "\n")

seq_all <- filter(snp_filter,mut %in% c("C>T","G>A"))$motif_seq
table(substr(seq_all, 5, 6))
binom.test(42+68,sum(table(substr(seq_all, 5, 6))),0.5)
ggseqlogo(seq_all, seq_type = "dna",  method="prob") +
  theme_bw()
ggsave("FigD_iSNV_Seqlogo.pdf", width = 2, height = 3)


df <- read.table(text="
Gene	Time	FC_log2	lfcSE	Sample
APOBEC3B	6h	0.1595134	0.07620862	ASFV/mock
APOBEC3B	12h	1.0902016	0.07829925	ASFV/mock
APOBEC3B	24h	2.5025041	0.08949237	ASFV/mock
", header = TRUE, sep = "\t")


#APOBC3B expression
# 时间轴顺序
df$Time <- factor(df$Time, levels = c("6h","12h","24h"))


ggplot(df, aes(x = Time, y = FC_log2, fill = Time)) +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = FC_log2 - lfcSE, ymax = FC_log2 + lfcSE),
                position = position_dodge(width = 0.8), width = 0.25) +
  #scale_fill_manual(values = cols) +
  # 标注padj，x对应Time位置，向右偏移0.2对准ASFV柱子
  annotate("text", x = 1 , y = 0.1595134 + 0.07620862 + 0.12, label = "P=0.061", size = 3.2) +
  annotate("text", x = 2 , y = 1.0902016 + 0.07829925 + 0.12, label = "P=4.08e-43", size = 3.2) +
  annotate("text", x = 3 , y = 2.5025041 + 0.08949237 + 0.12, label = "P=1.26e-170", size = 3.2) +
  labs( y = "Fold change of APOBEC3B expression level\n(ASFV vs. Mock, log2)") +
  theme_classic() +
  theme(
    axis.text = element_text(size =11),
    axis.title = element_text(size =12),
    legend.position = "top"
  )
ggsave("FigE_APOBEC3B_bar.pdf", width = 2, height = 3)


# library(DESeq2)
library(tidyverse)
setwd("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/RNAseq/")
# 
# ## 读featureCounts结果
# fc <- read.delim("featureCounts_S0.txt", header = TRUE, sep = "\t", check.names = FALSE)
# colnames(fc)
# countData <- fc %>%
#   dplyr::select(-Chr,-Start,-End,-Strand,-Length) %>%
#   column_to_rownames("Geneid")
# 
# # =====================
# # 样本顺序：00h(R1,R2,R3), 06h(R1,R2,R3), 12h(R1,R2,R3), 24h(R1,R2,R3)
# # =====================
# sample_names <- colnames(countData)
# time_group <- c(
#   rep("00h",3),
#   rep("06h",3),
#   rep("12h",3),
#   rep("24h",3),
#   rep("Mock_06h",3),
#   rep("Mock_12h",3),
#   rep("Mock_24h",3)
# )
# 
# colData <- data.frame(
#   sample = sample_names,
#   time = factor(time_group, levels = c("00h","06h","12h","24h","Mock_06h","Mock_12h","Mock_24h")) # 00h作为参照水平
# )
# rownames(colData) <- sample_names
# 
# # 构建DESeq对象；featureCounts输出可能存在小数，round强制整数
# dds <- DESeqDataSetFromMatrix(
#   countData = round(countData),
#   colData   = colData,
#   design    = ~ time
# )
# 
# # 过滤极低表达基因：总count≥10
# dds <- dds[ rowSums(counts(dds)) >= 10, ]
# 
# dds <- DESeq(dds)
# 
# # 两两对比：全部时间点 VS 00h
# res_06v00  <- results(dds, contrast = c("time","06h","Mock_06h"), alpha = 0.05)
# res_12v00  <- results(dds, contrast = c("time","12h","Mock_12h"), alpha = 0.05)
# res_24v00  <- results(dds, contrast = c("time","24h","Mock_24h"), alpha = 0.05)
# 
# 
# # 转成data.frame，Geneid列
# res_06v00_df <- as.data.frame(res_06v00) %>% rownames_to_column("Geneid")
# res_12v00_df <- as.data.frame(res_12v00) %>% rownames_to_column("Geneid")
# res_24v00_df <- as.data.frame(res_24v00) %>% rownames_to_column("Geneid")
# 
# res_06v00_df$Sample="06vMock06"
# res_12v00_df$Sample="12vMock12"
# res_24v00_df$Sample="24vMock24"
# 
# 
# DEG=rbind(
#       res_06v00_df,
#       res_12v00_df,                  
#       res_24v00_df)     #APOBEC1              #APOBEC2                                                          #APOBEC4
# filter(DEG,Geneid %in% c("ENSSSCG00000000668","ENSSSCG00000042447","ENSSSCG00000054928","ENSSSCG00000000091","ENSSSCG00000015561"))
# # 输出单独csv
# write.csv(DEG, "DEG_PAM.csv", row.names = FALSE)
# 
# DEG=fread("DEG_PAM.csv")
# head(DEG)
# filter(DEG,Geneid %in% c("ENSSSCG00000000091"))
# 
# # Gene	Time	FC_log2	lfcSE	Sample
# # APOBEC3B	06h	0	0	Mock
# # APOBEC3B	06h	0.1595134	0.07620862	ASFV
# # APOBEC3B	12h	0	0	Mock
# # APOBEC3B	12h	1.0902016	0.07829925	ASFV
# # APOBEC3B	24h	0	0	Mock
# # APOBEC3B	24h	2.5025041	0.08949237	ASFV


