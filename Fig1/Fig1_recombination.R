library(ggplot2)
############################################################

#重组事件+重组基因组

############################################################
# ---- 数据（黑名单过滤后） ----
df <- data.frame(
  methods = factor(c(2, 3, 4, 5, 6, 7)),
  RE_number = c(183, 149, 126, 101, 71, 36),
  recomb_genomes = c(357, 357, 357, 357, 357, 345)
)

# ---- 双轴缩放参数 ----
# 左轴: 0 ~ 560 (RE number)
# 右轴: 345 ~ 370 (recombinant genomes)
y_left_max  <- 200
y_right_min <- 340
y_right_max <- 365
scale_coef  <- y_left_max / (y_right_max - y_right_min)  # 560 / 25 = 22.4

# 将右轴数据缩放到左轴坐标系
df$recomb_scaled <- (df$recomb_genomes - y_right_min) * scale_coef

# ---- 绘图 ----
p <- ggplot(df, aes(x = methods)) +
  
  # 蓝色柱: RE number (左轴)
  geom_col(aes(y = RE_number), fill = "#4088B4", width = 0.6, alpha = 0.85) +
  
  # 红色折线+点: recombinant genomes (右轴, 已缩放)
  geom_line(aes(y = recomb_scaled, group = 1), color = "#E08251", linewidth = 0.8) +
  geom_point(aes(y = recomb_scaled), color = "#E08251", size = 2.5) +
  
  # 柱顶标注 RE number
  geom_text(aes(y = RE_number, label = RE_number),
            vjust = -0.5, size = 3.2, color = "grey20") +
  
  # 红点旁标注 recombinant genomes
  geom_text(aes(y = recomb_scaled, label = recomb_genomes),
            vjust = -1, size = 3.2,
            color = "#E08251", fontface = "bold") +
  
  # 左轴 + 右轴（反向变换）
  scale_y_continuous(
    name = "Number of recombinant events",
    limits = c(0, y_left_max),
    expand = c(0, 0),
    sec.axis = sec_axis(
      trans = ~ . / scale_coef + y_right_min,
      name = "Number of recombinant genomes",
      breaks = seq(340, 365, 5)
    )
  ) +
  
  # 横轴
  labs(x = "RDP4 methods number supported recombinant events") +
  
  # 主题
  theme_classic(base_size = 12) +
  theme(
    axis.title.y.left  = element_text(color = "#4088B4", face = "bold", margin = margin(r = 8)),
    axis.title.y.right = element_text(color = "#E08251", face = "bold", margin = margin(l = 8)),
    axis.text.y.left   = element_text(color = "#4088B4"),
    axis.text.y.right  = element_text(color = "#E08251"),
    axis.line.y.left   = element_line(color = "#4088B4"),
    axis.line.y.right  = element_line(color = "#E08251"),
    axis.ticks.y.left  = element_line(color = "#4088B4"),
    axis.ticks.y.right = element_line(color = "#E08251"),
    axis.title.x       = element_text(face = "bold", margin = margin(t = 8)),
    plot.margin        = margin(10, 12, 10, 10)
  )

p
# 同时保存 PNG
ggsave("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/recombination/0RDP/ASFV_all_RDP/SCG/FigA_RDP_methods_dual_axis.pdf",
       p, width = 4.5, height = 3, dpi = 300)



library(ggplot2)
library(reshape2)

# ============================================================
# 1. 读取数据
#    格式: 第一列 CenterPos, 后续每列一个参考序列的相似性值
#    支持 tab 分隔 (.txt) 或 csv
# ============================================================
# 替换为你的实际文件路径
#query MN318203.3 基因型22
sim_file <- "C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/recombination/0RDP/SimPlot20260829_genotype20.csv"

# 自动判断分隔符
if (grepl("\\.csv$", sim_file, ignore.case = TRUE)) {
  sim_df <- read.csv(sim_file, header = TRUE, check.names = FALSE)
} else {
  sim_df <- read.delim(sim_file, header = TRUE, sep = "\t", check.names = FALSE)
}

# ============================================================
# 2. 宽表转长表 (ggplot2 需要长格式)
#    CenterPos 作为 x, 变量名作为分组, 值作为 y
# ============================================================
sim_long <- sim_df %>%
  tidyr::pivot_longer(
    cols      = -CenterPos,     # 除 CenterPos 外的所有列
    names_to  = "Sequence",
    values_to = "Similarity"
  )

# ============================================================
# 3. 绘图 — SimPlot 风格步进折线图
# ============================================================
p <- ggplot(sim_long, aes(x = CenterPos, y = Similarity, color = Sequence)) +
  
  # 步进线 (step plot), direction = "hv" 先水平后垂直, 匹配 SimPlot 风格
  #geom_step( direction = "hv") +
  geom_step( direction = "hv",linewidth = 1) +
  # 颜色: 蓝 + 红 (可按序列数扩展)
  scale_color_manual(values = c(
    "MN336500.3"  = "#0000FF",   # 蓝 基因型22
    "NC_044941.1" = "#FF0000"    # 红 基因型1
    # 如有更多序列, 在此追加: "序列名" = "颜色"
  )) +
  
  # X 轴: 基因组位置
  scale_x_continuous(
    name = "Genome alignment position (nt)",
    limits = c(0, max(sim_df$CenterPos) * 1.02),
    expand = c(0, 0),
    labels = scales::comma
  ) +
  
  # Y 轴: 相似性
  scale_y_continuous(
    name = "Similarity",
    limits = c(0.8, 1.0),
    expand = c(0, 0),
    breaks = seq(0.5, 1.0, 0.1)
  ) +
  
  # 主题: 白底 + 浅灰网格, 匹配 SimPlot 风格
  theme_classic() +
  theme(
    legend.title     = element_blank(),
    legend.position  = "top"
  )
p

out_pdf <- sub("\\.[^.]+$", "_simplot_FigB.pdf", sim_file)
ggsave(out_pdf, p, width = 10, height = 5, dpi = 300)



############################################################

#重组事件+重组基因组

############################################################
# 加载所需包
library(ggplot2)
library(data.table)
library(dplyr)
library(stringr)
library(gggenes)
library(patchwork)
library(tidyverse)
#=====================================================================
# 1 重组断点数据处理（完全沿用你的筛选逻辑，只取Break_Begin/Break_End原始断点）
#=====================================================================
RE <- fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/recombination/0RDP/ASFV.csv", header = T,skip = 2)
#RE <- fread("C:/Users/Ke-jia Shan/Desktop/ASFVtest.csv", header = T,skip = 2)

methods <- c("RDP", "GENECONV", "Bootscan", "Maxchi", "Chimaera", "SiSscan", "3Seq")

RE$`Recombinant Sequence(s)` 
# 7列全部非空筛选
RE_filtered <- RE[RDP != "" & GENECONV != "" & Bootscan != "" & Maxchi != "" & Chimaera != "" & SiSscan != "" & `3Seq` != ""]

# 断点列重命名+清洗为纯数字整数
setnames(RE_filtered, c(7,8), c("Break_Begin","Break_End"))
RE_filtered[, Break_Begin := as.integer(gsub("[^0-9]", "", Break_Begin))]
RE_filtered[, Break_End   := as.integer(gsub("[^0-9]", "", Break_End))]

# 统计每种方法P<0.05显著条数
RE_filtered[, support_count := {
  sig_matrix <- sapply(.SD, function(x) {
    x_num <- suppressWarnings(as.numeric(x))
    !is.na(x_num) & x_num < 0.05
  })
  rowSums(sig_matrix, na.rm = TRUE)
}, .SDcols = methods]

# 筛选≥3种方法支持的可靠重组事件
RE_sig <- RE_filtered[support_count >= 3] 

# 筛选≥3种方法支持的可靠重组事件
RE_sig <- RE_filtered[support_count >= 5] 
mid=Reduce(intersect, list(filter(data.frame(RE),grepl("NC_044941",`Minor.Parental.Sequence.s.`))$`Recombination.Event.Number`,
                                                     filter(data.frame(RE),grepl("MN336500",`Major.Parental.Sequence.s.`))$`Recombination.Event.Number`,
                                                    filter(data.frame(RE),grepl("MN318203",`Recombinant.Sequence.s.`))$`Recombination.Event.Number`))
filter(data.frame(RE_sig),`Recombination.Event.Number` %in% mid)[,3:4]


# ========== 仅使用原始两端断点，不计算中点 ==========
all_recomb_breakpoints <- c(RE_sig$Break_Begin, RE_sig$Break_End)

# Upset交集矩阵（保留变量，后续可直接画UpSetR）
upset_mat <- RE_sig[, ..methods] %>%
  mutate(across(all_of(methods), ~ifelse(suppressWarnings(as.numeric(.x))<0.05 & !is.na(suppressWarnings(as.numeric(.x))),1,0))) %>%
  as.data.frame()

#=====================================================================
# 2 滑动窗口计数绘图数据（100bp窗口，10bp步长）
#=====================================================================
chrom_length <- 190584
window_size <- 1000
step_size   <- 10

starts <- seq(1, chrom_length - window_size + 1, by = step_size)
ends   <- starts + window_size - 1

counts <- sapply(1:length(starts), function(i){
  sum(all_recomb_breakpoints >= starts[i] & all_recomb_breakpoints <= ends[i])
})

df_count <- data.frame(
  pos = (starts + ends) / 2,
  count = counts
)

#=====================================================================
# 3 上图：重组断点密度折线图
#=====================================================================
max(df_count$pos)
min(df_count$pos)
total_len <- 190000
left_var_end    <- 40000
right_var_start <- total_len - 25000
left_ratio_end    <- left_var_end 
right_ratio_start <- right_var_start 
genome_h <- 3
bar_top<-40

p1 <- ggplot(df_count, aes(x = pos, y = count)) +
  geom_line(color = "darkred", linewidth = 0.9) +
  scale_x_continuous(limits = c(0, chrom_length)) +
  labs(x = "", y = "Breakpoint number per 100bp window") +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(t = 10, r = 10, b = 0, l = 10)
  )+
  # 顶层基因组矩形
  annotate("rect", xmin = 1, xmax = left_ratio_end,
           ymin = bar_top, ymax = bar_top + genome_h, fill = "#d9d9d9") +
  annotate("rect", xmin = left_ratio_end, xmax = right_ratio_start,
           ymin = bar_top, ymax = bar_top + genome_h, fill = "#8c8c8c") +
  annotate("rect", xmin = right_ratio_start, xmax = total_len ,
           ymin = bar_top, ymax = bar_top + genome_h, fill = "#d9d9d9") +
  annotate("text", x = left_ratio_end/2, y = bar_top + genome_h/2, label = "left\nvariable", size = 2.8) +
  annotate("text", x = (left_ratio_end + right_ratio_start)/2, y = bar_top + genome_h/2, label = "core", size = 2.8, color = "white") +
  annotate("text", x = right_ratio_start+(total_len-right_ratio_start )/2, y = bar_top + genome_h/2, label = "right\nvariable", size = 2.8) 

p1
#=====================================================================
# 4 ASFV基因注释处理 + 原版渐变配色逻辑
#=====================================================================
ASFV <- fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/miniprot/ASFV_uniprot.annotation.final.stat.txt",header = T)
ASFV=filter(ASFV, !GnmID %in% c("OR944085.1", "OR944086.1", "OR944087.1", "OR944088.1", "OR944089.1", "OR944090.1", "OR944091.1"))

head(ASFV)
all_gnm <- unique(ASFV$GnmID)
mid <- ASFV %>%
  group_by(GeneName) %>%  # 按 OrthoGroup 分组
  summarise(
    gene_count=n(),
    # 统计1：该 OrthoGroup 已有的唯一基因组
    unique_GnmID_count = n_distinct(GnmID, na.rm = TRUE),  # 已有数量
    #unique_GnmID_list = paste(sort(unique(GnmID[!is.na(GnmID)])), collapse = ", "),  # 已有列表（排序后更整齐）
    # 统计2：该 OrthoGroup 缺失的基因组
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
core=filter(mid,unique_GnmID_count== gene_count,loss_GnmID_count/357<=0.05)

result <- filter(ASFV,GeneName %in% core$GeneName) %>%
  group_by(GeneName) %>%
  filter(Ident == max(Ident)) %>%  # 保留 Ident 最大的行
  slice(1) %>%                     # 如果有多个最大值，只保留第一行
  ungroup()

mid=data.frame(result[,7:8])


ASFV <- ASFV  %>% filter(GnmID == "C_AA040130.1")
# 基础列构建
df_plot <- ASFV %>%
  mutate(
    cluster = as.character(Group),
    subcluster = str_extract(GeneName, "(?<=_c)\\d+") %>% as.integer() %>% replace_na(1)
  )


# 分组配色字典（固定色值不变）
target_groups <- c("100", "110", "300", "360", "505","core")
group_color   <- c("100"="#c82020","110"="#2667a0","300"="#359035","360"="#803b8b","505"="#e67000","core"="#FFFF33")
group_light   <- c("100"="#c82020","110"="#2667a0","300"="#359035","360"="#803b8b","505"="#e67000","core"="#FFFF33")
df_plot[GeneName %in% core$GeneName & !cluster %in% target_groups ]$cluster="core"
# subcluster深浅渐变RGB混色代码
df_plot <- df_plot %>%
  group_by(cluster) %>%
  mutate(
    n_sub = max(subcluster, na.rm = TRUE),
    frac = ifelse(n_sub > 1 & !is.na(n_sub), (subcluster - 1)/(n_sub - 1), 0.5)
  )

df_plot <- df_plot  %>%
  rowwise() %>%
  mutate(
    fill_color = case_when(
      cluster %in% target_groups ~ rgb(
        (1-frac)*col2rgb(group_light[cluster])[1] + frac*col2rgb(group_color[cluster])[1],
        (1-frac)*col2rgb(group_light[cluster])[2] + frac*col2rgb(group_color[cluster])[2],
        (1-frac)*col2rgb(group_light[cluster])[3] + frac*col2rgb(group_color[cluster])[3],
        maxColorValue = 255),
      TRUE ~ "gray70"
    )
  ) %>%
  ungroup()

#=====================================================================
# 排布规则：每5个基因一行，从上到下依次换行
#=====================================================================
df_plot <- df_plot %>%
  mutate(gene_order_idx = row_number()) %>%
  mutate(track_row = (gene_order_idx - 1) %/% 5)

#=====================================================================
# 5 下图：gggenes基因结构图
#=====================================================================
p2 <- ggplot(df_plot, aes(xmin = Start, xmax = End, y = track_row, fill = fill_color)) +
  geom_gene_arrow(
    arrowhead_height = unit(3, "mm"),
    arrowhead_width = unit(2, "mm"),
    linewidth = 0
  ) +
  scale_fill_identity() +
  scale_x_continuous(limits = c(0, chrom_length)) +
  scale_y_reverse() +  
  geom_text(
    data = subset(df_plot, Group %in% c("RIR1","CD2H","EP153","148R","PP220")),
    aes(x = (Start + End)/2,y=-0.7+track_row, label = Group),
    size = 3.2,
    hjust = 0.5,
    vjust = 0.5
  ) +
  labs(x = "Position (C_AA040130.1)", y = "") +
  theme_void() +
  theme(
    axis.line.x = element_line(),
    axis.text.x = element_text(size = 10),
    axis.ticks.x = element_line(),
    plot.margin = margin(t = 0, r = 10, b = 10, l = 10),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )


p2

#=====================================================================
# 6 上下拼图 X轴对齐
#=====================================================================

p<-p1 / p2 + plot_layout(heights = c(3, 2))
ggsave("FigC_breakpoint.pdf",p , width = 10, height = 5, dpi = 300)

####################################################################

# ############################################################################################################################################################
# 
# #根据重组区域的数目着色
# 
# ############################################################################################################################################################
# # 加载所需包
# library(ggplot2)
# library(data.table)
# library(dplyr)
# library(stringr)
# library(gggenes)
# library(patchwork)
# 
#=====================================================================
# 1 重组断点数据处理（完全沿用你的筛选逻辑）
#=====================================================================
RE <- fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/recombination/0RDP/ASFV_RE_Upset.csv", header = T)
methods <- c("RDP", "GENECONV", "Bootscan", "Maxchi", "Chimaera", "SiSscan", "3Seq")

# 7列全部非空筛选
RE_filtered <- RE[RDP != "" & GENECONV != "" & Bootscan != "" & Maxchi != "" & Chimaera != "" & SiSscan != "" & `3Seq` != ""]

# 断点列重命名+清洗为纯数字整数
setnames(RE_filtered, c(7,8), c("Break_Begin","Break_End"))
RE_filtered[, Break_Begin := as.integer(gsub("[^0-9]", "", Break_Begin))]
RE_filtered[, Break_End   := as.integer(gsub("[^0-9]", "", Break_End))]

# 统计每种方法P<0.05显著条数
RE_filtered[, support_count := {
  sig_matrix <- sapply(.SD, function(x) {
    x_num <- suppressWarnings(as.numeric(x))
    !is.na(x_num) & x_num < 0.05
  })
  rowSums(sig_matrix, na.rm = TRUE)
}, .SDcols = methods]

# RE_filtered$Recom_ID=paste0(RE_filtered$Break_Begin,"_", RE_filtered$Break_End)
# head(RE_filtered)
# cd2h_region_str <- ASFV[grepl("CD2H", GeneName, ignore.case = TRUE)]$recomb_regions
# target_ids <- strsplit(cd2h_region_str, ",")[[1]]
#
# matched_recomb <- RE_filtered[Recom_ID %in% target_ids]
#
# print(matched_recomb)

# 筛选≥3种方法支持的可靠重组事件
RE_sig <- RE_filtered[support_count >= 3][, .(
  Break_Begin,
  Break_End
)]
RE_sig <- unique(RE_sig[Break_Begin!=Break_End])
chrom_length <- 190584

# ========== 滑动窗口绘图用断点向量不变
all_recomb_breakpoints <- c(RE_sig$Break_Begin, RE_sig$Break_End)

#=====================================================================
# 2 滑动窗口计数绘图数据（1000bp窗口，10bp步长）
#=====================================================================
window_size <- 1000
step_size   <- 10

starts <- seq(1, chrom_length - window_size + 1, by = step_size)
ends   <- starts + window_size - 1

counts <- sapply(1:length(starts), function(i){
  sum(all_recomb_breakpoints >= starts[i] & all_recomb_breakpoints <= ends[i])
})

df_count <- data.frame(
  pos = (starts + ends) / 2,
  count = counts
)

#=====================================================================
# 3 上图：重组断点密度折线图（完全不变）
#=====================================================================
p1 <- ggplot(df_count, aes(x = pos, y = count)) +
  geom_line(color = "darkred", linewidth = 0.9) +
  scale_x_continuous(limits = c(1, chrom_length), expand = c(0, 0)) +
  labs(x = "", y = "Breakpoint number per 1000bp window") +
  theme_bw() +
  #geom_vline(xintercept = 57969,color = "darkred", linewidth = 0.9) +
  #geom_vline(xintercept = 60305,color = "darkred", linewidth = 0.9) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(t = 10, r = 10, b = 0, l = 10)
  )

#=====================================================================
# 4 ASFV基因注释 + 【核心修复：正确的重组区间重叠判断函数】
#=====================================================================
ASFV <- fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/miniprot/ASFV_uniprot.annotation.final.stat.txt",header = T) %>% filter(GnmID == "C_AA040130.1")
setDT(ASFV)

# ====================== 重写正确判断函数 ======================

calc_gene_recomb_hit <- function(gL, gR) {
  if(is.na(gL) || is.na(gR)) {
    return(list(count = 0, regions = NA_character_))
  }

  hit_keys <- character(0)
  for (idx in 1:nrow(RE_sig)) {
    s <- RE_sig$Break_Begin[idx]
    e <- RE_sig$Break_End[idx]

    if (s < e) {
      overlap <- !(gR < s | gL > e)
    } else {
      seg1_overlap <- !(gR < s | gL > chrom_length)
      seg2_overlap <- !(gR < 1  | gL > e)
      overlap <- seg1_overlap || seg2_overlap
    }
    if (overlap) {
      key <- paste0(s, "_", e)
      hit_keys <- unique(c(hit_keys, key))
    }
  }
  count_val <- length(hit_keys)
  region_str <- paste(hit_keys, collapse = ",")
  return(list(count = count_val, regions = region_str))
}
#ASFV[, recomb_hit := mapply(calc_gene_recomb_hit, gL, gR)]
# ==============================================================

# 基因区间归一化（负链自动转gL-gR线性区间）
ASFV[, gL := pmin(Start, End)]
ASFV[, gR := pmax(Start, End)]

res <- lapply(1:nrow(ASFV), function(i){
  calc_gene_recomb_hit(ASFV$gL[i], ASFV$gR[i])
})

# 拆分结果写入两列
ASFV[, recomb_hit := sapply(res, `[[`, "count")]
ASFV[, recomb_regions := sapply(res, `[[`, "regions")]

# 查看计数分布
cat("===== 所有基因重组计数分布 =====\n")
table(ASFV$recomb_hit)

# 标注阈值
label_threshold <- 10
ASFV[, gene_label := ifelse(recomb_hit >= label_threshold, GeneName, NA_character_)]

# PP220核验
cat("========== PP220 重组计数结果 ==========\n")
print(ASFV[grepl("PP220", GeneName, ignore.case = TRUE), .(GeneName, gL, gR, recomb_hit,recomb_regions)])
print(ASFV[grepl("EP153", GeneName, ignore.case = TRUE), .(GeneName, gL, gR, recomb_hit,recomb_regions)])
print(ASFV[grepl("CD2H", GeneName, ignore.case = TRUE), .(GeneName, gL, gR, recomb_hit,recomb_regions)])
print(ASFV[grepl("148R", GeneName, ignore.case = TRUE), .(GeneName, gL, gR, recomb_hit,recomb_regions)])
print(ASFV[grepl("RIR1", GeneName, ignore.case = TRUE), .(GeneName, gL, gR, recomb_hit,recomb_regions)])

#=====================================================================
# 按重组次数渐变配色：白→黄→橙→红
#=====================================================================
df_plot <- ASFV %>%
  mutate(gene_order_idx = row_number(),
         track_row = (gene_order_idx - 1) %/% 5)

hit_min <- min(df_plot$recomb_hit)
hit_max <- max(df_plot$recomb_hit)

df_plot <- df_plot %>%
  mutate(
    norm_val = if (hit_max == hit_min) 0 else (recomb_hit - hit_min)/(hit_max - hit_min),
    fill_color = if (hit_max == hit_min) {
      "gray80"
    } else {
      rgb(1, 1 - norm_val * 0.9, 1 - norm_val, maxColorValue = 1)
    }
  )

#=====================================================================
# 5 下图：基因箭头 + 居中白色加粗标签（≥10自动显示）
#=====================================================================
p2 <- ggplot(df_plot, aes(xmin = Start, xmax = End, y = track_row, fill = fill_color)) +
  geom_gene_arrow(
    arrowhead_height = unit(3, "mm"),
    arrowhead_width = unit(2, "mm"),
    linewidth = 0
  ) +
  # 居中标签，不会被X轴截断删除
  geom_gene_label(
    aes(label = gene_label),
    size = 1.7,
    align = "middle",
    colour = "white",
    fontface = "bold",
    fill = NA
  ) +
  scale_fill_identity() +
  scale_x_continuous(limits = c(1, chrom_length), expand = c(0, 0)) +
  scale_y_reverse() +
  labs(x = "Position (C_AA040130.1)", y = "") +
  theme_void() +
  theme(
    axis.line.x = element_line(),
    axis.text.x = element_text(size = 10),
    axis.ticks.x = element_line(),
    plot.margin = margin(t = 0, r = 10, b = 10, l = 10),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

#=====================================================================
# 拼图输出
#=====================================================================
p1 / p2 + plot_layout(heights = c(3, 2))
print(final_plot)

# 保存（需要取消注释）
# ggsave("ASFV_Recomb_Fixed_GeneMap.png", final_plot, width = 18, height = 6, dpi = 300)
# ggsave("ASFV_Recomb_Fixed_GeneMap.pdf", final_plot, width = 18, height = 6)





####################################################################

#gene在多少个重组区间

####################################################################
write.table(df_plot,"C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/recombination/0RDP/ASFV_RE_gene.tsv",row.names = F,quote = F,sep="\t")
RE <- fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/recombination/0RDP/ASFV_RE_gene.tsv",sep="\t")

library(data.table)
library(ggplot2)

# 按重组事件数降序排序
RE_order <- RE[order(recomb_hit,decreasing = T)]
RE_order$Group <- sub("_YP.*", "", RE_order$Group)

RE_order$Gene<-paste0(RE_order$Group,":",RE_order$Start,"-",RE_order$End)

# # 绘图
# ggplot(filter(RE_order,recomb_hit>=5), aes(x = factor(Gene, levels = Gene), y = recomb_hit)) +
#   geom_col(fill = "#4088B4", width = 0.7) +
#   labs(
#     y = "Number of recombination events",
#     x = "Gene",
#     title = "Recombination events count per gene (descending order)"
#   ) +
#   theme_bw() +
#   theme(
#     plot.title = element_text(hjust = 0.5),
#     # X轴标签竖直90度摆放
#     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6.5),
#     # 清理多余网格线
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor = element_blank()
#   )
# 
# 
# unique(RE_order$Group )

head(filter(RE_order,recomb_hit>=5)[,c(1,8,10,15,22)])
Gene_Function=fread("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/recombination/0RDP/ASFV_all_RDP/SCG/gene_function.csv")
plot=merge(filter(RE_order)[,c(7,8,10,15,22)],Gene_Function[,c(1,3)],by.x = "GeneName",by.y = "GeneName")
plot <- plot[order(recomb_hit,decreasing = T)]
plot[grepl("^100",plot$Gene),]$Gene <- gsub("100","MGF-100",plot[grepl("^100",plot$Gene),]$Gene )
plot[grepl("^110",plot$Gene),]$Gene <- gsub("110","MGF-110",plot[grepl("^110",plot$Gene),]$Gene )
plot[grepl("^300",plot$Gene),]$Gene <- gsub("300","MGF-300",plot[grepl("^300",plot$Gene),]$Gene )
plot[grepl("^360",plot$Gene),]$Gene <- gsub("360","MGF-360",plot[grepl("^360",plot$Gene),]$Gene )
plot[grepl("^505",plot$Gene),]$Gene <- gsub("505","MGF-505",plot[grepl("^505",plot$Gene),]$Gene )
plt=filter(plot,recomb_hit>=10)
ggplot(filter(plt), aes(x = factor(Gene, levels = Gene), y = recomb_hit,fill = Function)) +
  geom_col(width = 0.7) +
  labs(
    y = "Number of recombination events",
    x = "Gene"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    # X轴标签竖直90度摆放
    axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1, size = 10),
    # 清理多余网格线
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position  = "top"
  ) +
  scale_fill_manual(
    values = c(
      "Viral replication and assembly" = "#58A8F6",
      "Membrane/envelope-associated"   = "#219180",
      "Anti-host immunity"             = "#E53C38",
      "Unknown function"               = "#A9A9A9",
      "Anti-host immunity; Viral replication and assembly" = "#F9A74D",
      "Viral replication and assembly; Membrane/envelope-associated" = "#66C755",
      "Anti-host immunity; Membrane/envelope-associated" = "#B545A2",
      "Anti-host immunity; Viral replication and assembly; Membrane/envelope-associated" = "#2B644A"
    )
  )
dev.off()
ggsave("FigE.pdf", width = 13, height = 5, dpi = 300)



table(plot$category)
