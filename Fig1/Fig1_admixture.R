library(ggtree)
library(ggplot2)
library(dplyr)
library(reshape2)
library(ape)
library(treeio)
library(data.table)
# ====================== 1. 基础数据加载（你的代码，仅微调） ======================
setwd("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/miniprot/")

# 读取并定根进化树
tree_file <- "C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/ASFV_357_2_tree.root.nwk" 
tree <- read.tree(tree_file)
root_tips <- c("LC506465.1", "LC637659.1")
if (!all(root_tips %in% tree$tip.label)) {
  missing_tips <- setdiff(root_tips, tree$tip.label)
  stop(paste("以下根序列不在进化树中：", paste(missing_tips, collapse = ", ")))
}
tree_rooted <- root(tree, outgroup = root_tips, resolve.root = TRUE)
# 删除指定tip
tree_rooted  <- drop.tip(tree_rooted , tip = root_tips)

# 处理基因型数据
genotype_df=read.table("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/miniprot/meta.tsv",header = T,sep="\t")
head(genotype_df)
colnames(genotype_df)[1:2] <- c("SampleID","Genotype")
table(genotype_df$Genotype)
root_df <- data.frame(SampleID = c("LC506465.1", "LC637659.1"), Genotype = c("Outgroup","Outgroup"))

# root_df <- data.frame(
#   SampleID = root_tips,
#   Genotype = "Outgroup",
#   Length = 190000,
#   Geo_Location = "Outgroup",
#   Host = "Outgroup",
#   date=NA,
#   years = NA,
#   stringsAsFactors = FALSE
# )

#genotype_df <- rbind(genotype_df, root_df)
genotype_df_ordered <- genotype_df[match(tree_rooted$tip.label, genotype_df$SampleID), ]
genotype_df_ordered$Genotype[is.na(genotype_df_ordered$Genotype)] <- "unknown"
tip_genotype <- setNames(genotype_df_ordered$Genotype, tree_rooted$tip.label)
head(genotype_df)
# 基因型配色
genotype_colors <- c(
  "1" = "#E41A1C", "2" = "#377EB8", "9" = "#FF7F00","8"="#bf5500", "20" ="#FFC0CB","22"="#EE4C97FF","1/2" = "#4DAF4A",
  "Outgroup" = "black", "unknown" = "grey70"
)

# p_tree <- ggtree(tree_rooted,  
#                  layout = "rectangular",  
#                  #branch.length = "none", 
#                  lwd = 1) +
#   #theme_tree2() +
#   geom_tippoint(aes(color = tip_genotype[label]), size = 2, shape = 19) +
#   geom_tree(aes(color = tip_genotype[label]), size = 1) +
#   #geom_tiplab(size = 1.5, offset = 0.5) +  # 标签偏移，避免和条形图重叠
#   scale_color_manual(values = genotype_colors, name = "Genotype") +
#   theme_tree2(legend.position = "none")
# p_tree


# ====================== 2. 处理Admixture数据（适配facet_plot） ======================
# 读取重组数据
recomb_path <- "C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/recombination"
setwd(recomb_path)
k_value <- 4
q_data <- read.table(paste0("asfv_admix_final.", k_value, ".Q"), header = FALSE)
samples <- read.table("asfv_admix_final.fam", header = FALSE)[, 2]

# 数据整理
colnames(q_data) <- paste0("Ancestry_", 1:k_value)
q_data$Sample <- samples
q_data$Sample[grepl("^\\d+",q_data$Sample)]=paste0("NC_",q_data$Sample[grepl("^\\d+",q_data$Sample)])
q_data$Sample[grepl("^AA\\d+",q_data$Sample)]=paste0("C_",q_data$Sample[grepl("^AA\\d+",q_data$Sample)])
q_data <-filter(q_data, !Sample %in% c("OR944085.1", "OR944086.1", "OR944087.1", "OR944088.1", "OR944089.1", "OR944090.1", "OR944091.1"))

#genotype_df$Sample <- sub("^[^_]*_", "", genotype_df$SampleID) # 匹配Sample格式
#genotype_df$Sample <-  genotype_df$SampleID

unique(q_data$Sample)
q_melt <- reshape2::melt(q_data, id.vars = "Sample")
unique(q_melt$Sample)
q_melt$Sample[!q_melt$Sample %in% genotype_df$Sample]
q_melt <- merge(q_melt, genotype_df, by.x = "Sample", by.y = "SampleID")
colnames(q_melt)[1]="SampleID"

# 关键：保证样本顺序与进化树完全一致（facet_plot的核心要求）
#tree_rooted$tip.label <- sub("^[^_]*_", "", tree_rooted$tip.label)
q_melt$SampleID <- factor(q_melt$SampleID, levels = tree_rooted$tip.label)
unique(q_melt$SampleID)
# ====================== 3. 核心：用facet_plot画树+条形图（100%稳定） ======================
# 第一步：画基础进化树

q_melt$SampleID[!q_melt$SampleID %in% tree_rooted$tip.label]
unique(q_melt$SampleID)
# # 第二步：用facet_plot添加堆叠条形图（ggtree官方推荐，无坐标冲突）
# p_final <- facet_plot(
#   p_tree,
#   panel = "Admixture",       # 注释面板名称
#   data = q_melt,             # 重组数据
#   geom = geom_bar,           # 用geom_bar
#   aes(x = value, fill = variable),
#   stat = "identity",
#   position = "stack",        # 堆叠
#   width = 1,                 # 条形宽度
#   color = "black",           # 边框
#   linewidth = 0.1,           # 边框粗细
#   orientation = "y"          # 水平条形图（关键！适配树的布局）
# ) +
#   # 配色+图例调整
#   scale_fill_manual(values = c("#E41A1C","purple", "#377EB8","#FF7F00","grey"), name = paste0("Ancestry (K=",k_value,")")) +
#   theme(
#     legend.position = "right",
#     legend.key.size = unit(0.5, "cm"),
#     legend.text = element_text(size = 8),
#     legend.title = element_text(size = 9),
#     strip.background = element_blank(),
#     # 调整面板间距（树和条形图的距离）
#     panel.spacing = unit(0.2, "cm")
#   )
# 
# # 显示图形
# print(p_final)
# 
# ============================================
# 1. 设置工作目录和文件路径
# ============================================

ani_file <- "C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/miniprot/fastani_results.3k.txt"

# 读取ANI结果（保留为data.table）
ANI <- fread(ani_file, header = FALSE)
setnames(ANI, c("Query", "Reference", "ANI_value", "Query_frags", "Reference_frags"))

# 过滤：只保留树中存在的样本
ANI_filtered <- ANI[Query %in% tree_rooted$tip.label & Reference %in% tree_rooted$tip.label, ]
cat("过滤后ANI数据行数:", nrow(ANI_filtered), "\n")

# 处理空数据
if (nrow(ANI_filtered) == 0) {
  warning("过滤后ANI数据为空！将创建全NA的矩阵")
  ANI_wide <- data.frame(Query =tree_rooted$tip.label)
  for (tip in tree_tips) {
    ANI_wide[[tip]] <- NA
  }
} else {
  # data.table的dcast
  ANI_wide <- dcast(ANI_filtered, Query ~ Reference, value.var = "ANI_value", fill = NA)
  ANI_wide <- as.data.frame(ANI_wide)
}

# 确保行顺序一致
ANI_wide <- ANI_wide[match(tree_rooted$tip.label, ANI_wide$Query), ]

# 补充缺失行
missing_rows <- setdiff(tree_rooted$tip.label, ANI_wide$Query)
if (length(missing_rows) > 0) {
  cat("补充缺失的行:", paste(missing_rows, collapse = ", "), "\n")
  missing_df <- data.frame(Query = missing_rows)
  for (col in colnames(ANI_wide)[-1]) {
    missing_df[[col]] <- NA
  }
  ANI_wide <- rbind(ANI_wide, missing_df)
  ANI_wide <- ANI_wide[match(tree_tips, ANI_wide$Query), ]
}

# 补充缺失列
existing_cols <- intersect(tree_rooted$tip.label, colnames(ANI_wide))
missing_cols <- setdiff(tree_rooted$tip.label, existing_cols)
if (length(missing_cols) > 0) {
  cat("补充缺失的列:", paste(missing_cols, collapse = ", "), "\n")
  for (col in missing_cols) {
    ANI_wide[[col]] <- NA
  }
}
ANI_wide <- ANI_wide[, c("Query", tree_rooted$tip.label), drop = FALSE]

# 转换为矩阵
ANI_matrix <- as.matrix(ANI_wide[, -which(colnames(ANI_wide) == "Query")])
rownames(ANI_matrix) <- ANI_wide$Query

# 处理对称性和NA值
ANI_sym <- pmax(ANI_matrix, t(ANI_matrix), na.rm = TRUE)
diag(ANI_sym) <- 100
#ANI_sym[is.na(ANI_sym)] <- 87

cat("ANI矩阵维度:", dim(ANI_sym), "\n")


############################################
## Packages
############################################

library(ggtree)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggnewscale)


############################################
## 1. Prepare annotation data
############################################

annotation_df <- genotype_df %>%
  select(SampleID,
         Geo_Location,
         years)


# match tree tip order
annotation_df <- annotation_df[
  match(tree_rooted$tip.label,
        annotation_df$SampleID),
]



#==== 1. 构建国家‑大洲映射表（覆盖你图中全部地点）====
country_continent_map <- tibble(
  Geo_Location = c(
    "Malawi",
    "South Africa",
    "Namibia",
    "Zambia",
    "Zaire",
    "Kenya",
    "Uganda",
    "Democratic Republic of the Congo",
    "Tanzania",
    "France",
    "Portugal",
    "Spain",
    "Cameroon",
    "Benin",
    "Ghana",
    "Dominican Republic",
    "Italy",
    "China",
    "Russia",
    "Zimbabwe",
    "Madagascar",
    "Mauritius",
    "Mozambique",
    "Georgia",
    "Lithuania",
    "Poland",
    "Viet Nam",
    "Serbia",
    "Belgium",
    "Czech Republic",
    "Moldova",
    "Ukraine",
    "Germany",
    "India",
    "Mongolia",
    "Philippines",
    "Estonia",
    "Singapore",
    "Hungary",
    "Timor-Leste",
    "South Korea",
    "Nigeria",
    "Latvia",
    "Indonesia",
    "Unknown"
  ),
  Continent = c(
    "Africa",
    "Africa",
    "Africa",
    "Africa",
    "Africa",
    "Africa",
    "Africa",
    "Africa",
    "Africa",
    "Europe",
    "Europe",
    "Europe",
    "Africa",
    "Africa",
    "Africa",
    "North America", #多米尼加，加勒比‑北美洲
    "Europe",
    "Asia",
    "Europe",
    "Africa",
    "Africa",
    "Africa",
    "Africa",
    "Asia", #格鲁吉亚 亚欧交界，病毒文献通常归亚洲
    "Europe",
    "Europe",
    "Asia",
    "Europe",
    "Europe",
    "Europe",
    "Europe",
    "Europe",
    "Europe",
    "Asia",
    "Asia",
    "Asia",
    "Europe",
    "Asia",
    "Europe",
    "Asia",
    "Asia",
    "Africa",
    "Europe",
    "Asia",
    "Unknown"
  )
)

#==== 2. 左连接，给你的注释数据框添加大洲列 ====
annotation_df <- annotation_df %>%
  left_join(country_continent_map, by = "Geo_Location")

rownames(annotation_df) <- annotation_df$SampleID

annotation_df <- annotation_df %>%
  select(Geo_Location,
         Continent,
         years)



############################################
## 2. Deal with missing Geo information
############################################

annotation_df$Geo_Location[
  is.na(annotation_df$Geo_Location) |
    annotation_df$Geo_Location == ""
] <- "Unknown"


# year must be numeric
annotation_df$years <- as.numeric(annotation_df$years)


# 大洲调色板
continent_colors <- c(
  "Africa"  = "#23452F",      # 深紫
  "Asia"    = "#009E73",      # 蓝‑绿色（水鸭绿，避开#4DAF4A）
  "Europe"  = "#D6E7A3",      # 玫紫豆沙
  "North America" = "#AA9A59",# 深海蓝
  "Unknown" = "white"
)


annotation_df$Continent[
  is.na(annotation_df$Continent) |
    annotation_df$Continent == ""
] <- "Unknown"

annotation_df$Continent <- factor(
  annotation_df$Continent,
  levels = names(continent_colors)
)


############################################
## 3. Generate Geo colors
############################################

#==== 1. 获取地理位置列表，移除Unknown ====
geo_levels <- unique(annotation_df$Geo_Location)
geo_names <- setdiff(geo_levels,"Unknown")
n_geo <- length(geo_names)

#==== 2. 高区分度基础调色板（28个高辨识度离散色，无相近色）====
# 注意：China 单独锁定为中国红
base_palette <- c(
  "Malawi"="#4eb899",
  "South Africa"="#73c990",
  "Namibia"="#a0d98b",
  "Zambia"="#d7e883",
  "Zaire"="#f8ee7c",
  "Kenya"="#e9e2a0",
  "Uganda"="#dcd8b4",
  "Democratic Republic of the Congo"="#c6c4cc",
  "Tanzania"="#a8a2cb",
  "France"="#b98eb8",
  "Portugal"="#d87c8c",
  "Spain"="#e97262",
  "Cameroon"="#ec6b50",
  "Benin"="#c96b62",
  "Ghana"="#a86872",
  "Dominican Republic"="#8882a0",
  "Italy"="#6285a8",
  "China"="#DE2910", # ==== 固定中国红 ====
  "Russia"="#8c9288",
  "Zimbabwe"="#bb9b6e",
  "Madagascar"="#f0a55b",
  "Mauritius"="#f2a84e",
  "Mozambique"="#d9b84c",
  "Georgia"="#b9c248",
  "Lithuania"="#9fcc52",
  "Poland"="#b1c578",
  "Viet Nam"="#d2c8b0",
  "Serbia"="#f0d2cd",
  "Belgium"="#f4c8dd",
  "Czech Republic"="#efd8e6",
  "Moldova"="#e7e2e7",
  "Ukraine"="#d9d7dc",
  "Germany"="#d1cfd4",
  "India"="#c3b9cf",
  "Mongolia"="#b4a1c7",
  "Philippines"="#a58ac0",
  "Estonia"="#9b78b8",
  "Singapore"="#b096b7",
  "Hungary"="#b9c1b9",
  "Timor-Leste"="#c9e2c7",
  "South Korea"="#d8e8bc",
  "Nigeria"="#e7ec9e",
  "Latvia"="#f3ed78",
  "Indonesia"="#f8e960"
)

#==== 3. 自动匹配你现有的国家名称；超出上面预设数量时再生成扩展色 ====
if(n_geo <= length(base_palette)){
  geo_palette <- base_palette[geo_names]
}else{
  #剩下多余样本生成高饱和度区分色，保证不和中国红冲突
  extra_colors <- colorRampPalette(brewer.pal(12, "Paired"))(n_geo - length(base_palette))
  geo_palette <- c(base_palette,extra_colors)[geo_names]
}

#==== 4. 追加Unknown为白色 ====
geo_colors <- c(geo_palette, Unknown = "white")


library(ggnewscale)

p_tree <- ggtree(tree_rooted,  
                 layout = "rectangular",  
                 #branch.length = "none", 
                 panel.width = 2,
                 lwd = 1) +
  #theme_tree2() +
  geom_tippoint(aes(color = tip_genotype[label]), size = 2, shape = 19) +
  geom_tree(aes(color = tip_genotype[label]), size = 1) +
  #geom_tiplab(size = 1.5, offset = 0.5) +  # 标签偏移，避免和条形图重叠
  scale_color_manual(values = genotype_colors, name = "Genotype") +
  theme_tree2(legend.position = "none")
p_tree

p_final <- facet_plot(
  p_tree,
  panel = "Admixture",
  data = q_melt,
  geom = geom_bar,
  aes(x = value,
      fill = variable),
  stat = "identity",
  position = "stack",
  width =1,
  #color = "black",
  linewidth = 0,panel.width = 0.2,
  orientation = "y"
)+
  
  scale_fill_manual(
    values=c(
      "#E41A1C",
      "purple",
      "#377EB8",
      "#FF7F00",
      "grey"
    ),
    name=paste0(
      "Ancestry (K=",
      k_value,
      ")"
    )
  )+
  
  theme(
    legend.position="right",
    legend.key.size=unit(0.5,"cm"),
    legend.text=element_text(size=8),
    legend.title=element_text(size=9),
    strip.background=element_blank(),
    panel.spacing=unit(0.2,"cm")
  )



#################################
# Geo annotation
#################################
p_final <- p_final +
  ggnewscale::new_scale_fill()

p_final <- gheatmap(
  p_final,
  annotation_df["Geo_Location"],
  offset=0.02,
  width=0.06,
  colnames=TRUE,
  font.size=3,
  color=NA
)+
  
  scale_fill_manual(
    values=geo_colors,
    breaks=setdiff(names(geo_colors),
                   "Unknown"),
    name="Location"
  )

p_final <- p_final +
  ggnewscale::new_scale_fill()


p_final <- gheatmap(
  p_final,
  annotation_df["years"],
  offset=0.03,
  width=0.06,
  colnames=TRUE,
  font.size=3,
  color=NA
)+
  
  scale_fill_gradientn(
    colours=c(
      "#FFF7BC",
      "#FEC44F",
      '#fd7e83',
      "#D95F0E",
      "#B30000",
      "#4B0082"
    ),
    name="Sampling year",
    na.value="white"
  )


p_final


# 检查有没有匹配失败NA条目
table(annotation_df$Continent, useNA = "ifany")

head(annotation_df)

p_final <- p_final +
  ggnewscale::new_scale_fill()

p_final <- gheatmap(
  p_final,
  annotation_df["Continent"],
  offset=0.01,
  width=0.06,
  colnames=TRUE,
  font.size=3,
  color=NA
)+
  scale_fill_manual(
    values=continent_colors,
    name="Continent",
    na.value="white"
  )

p_final


# 保存高分辨率图形（根据样本数调整尺寸）
ggsave(
  filename = "Fig1_ASFV_Tree+Admixture_Final.pdf",
  plot = p_final,
  width = 12,  # 宽度：树(8) + 条形图(4)
  height = 25, # 高度适配364个样本
  dpi = 300,
  device = "pdf",
  path = "C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/"
)

p_final <- p_final +
  ggnewscale::new_scale_fill()

str(ANI_sym)
p_hm <- gheatmap(
  p_tree,
  ANI_sym,
  offset = 0.01,
  #colnames_level = tree$tip.label,
  #width =2,
  colnames =F,
  #panel.width = 2,
  #color = "grey90",
  font.size = 3,
  color = NA   
) +
  scale_fill_gradientn(
    colors = c("#FFF7BC", "#FEC44F", "#fd7e83", "#D95F0E", "#B30000"),
    limits = c(88, 100),           # 你的数据范围 ~27–100%，压低限增强对比
    oob = scales::squish,
    na.value = "white",
    name = "3kb ANI (%)"
  ) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8)#,
    #panel.background = element_rect()
  )

print(p_hm)

ggsave(
  filename = "FigS1_ASFV_Tree+ANI.pdf",
  plot = p_hm,
  width = 14,  # 宽度：树(8) + 条形图(4)
  height = 7, # 高度适配364个样本
  dpi = 300,
  device = "pdf",
  path = "C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/"
)
# ============================================================
# 2. 读取并整理 Admixture 数据
#    最终结构：
#    Panel 1  = Phylogenetic tree
#    Panel 2  = Admixture K=1
#    Panel 3  = Admixture K=2
#    Panel 4  = Admixture K=3
#    Panel 5  = Admixture K=5
#    Panel 6  = Admixture K=6
#    Panel 7  = Admixture K=7
#    Panel 8  = Admixture K=8
#    Panel 9  = Admixture K=9
#    Panel 10 = Admixture K=10
# ============================================================
#BiocManager::install("ggtreeExtra")
library(ggplot2)
library(ggtree)
library(ggtreeExtra)
library(dplyr)
library(data.table)
library(grid)

# -----------------------------
# 工作目录
# -----------------------------
recomb_path <- "C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/recombination"
setwd(recomb_path)

# -----------------------------
# K值：1-10，但排除K=4
# -----------------------------
k_values <- setdiff(2:10, 4)

# -----------------------------
# 读取 fam 文件
# -----------------------------
fam <- read.table(
  "asfv_admix_final.fam",
  header = FALSE,
  stringsAsFactors = FALSE
)

samples <- fam[, 2]


# ============================================================
# 3. 准备样本名称
# ============================================================

samples <- samples

# 数字开头的样本名前面加 NC_
samples[grepl("^\\d+", samples)] <- paste0(
  "NC_",
  samples[grepl("^\\d+", samples)]
)

# AA开头的样本改成 C_
samples[grepl("^AA\\d+", samples)] <- paste0(
  "C_",
  samples[grepl("^AA\\d+", samples)]
)

# -----------------------------
# 需要排除的样本
# -----------------------------
remove_samples <- c(
  "OR944085.1",
  "OR944086.1",
  "OR944087.1",
  "OR944088.1",
  "OR944089.1",
  "OR944090.1",
  "OR944091.1"
)

samples <- samples[!samples %in% remove_samples]


# ============================================================
# 4. 设置 Admixture 配色
# ============================================================

# 最多需要10个祖源成分
ancestry_colors <- c(
  "#E41A1C",  # 1
  "#984EA3",  # 2
  "#377EB8",  # 3
  "#4DAF4A",  # 4
  "#FF7F00",  # 5
  "#A65628",  # 6
  "#F781BF",  # 7
  "#999999",  # 8
  "#66C2A5",  # 9
  "#FFD92F"   # 10
)

names(ancestry_colors) <- paste0(
  "Ancestry_",
  1:10
)


# ============================================================
# 5. 循环读取 K=1-10，跳过K=4
# ============================================================

q_melt_list <- list()

for (k in k_values) {
  
  cat("Processing Admixture K =", k, "\n")
  
  # -----------------------------
  # 检查Q文件
  # -----------------------------
  q_file <- paste0(
    "asfv_admix_final.",
    k,
    ".Q"
  )
  
  if (!file.exists(q_file)) {
    stop(
      paste0(
        "找不到文件：",
        q_file
      )
    )
  }
  
  # -----------------------------
  # 读取Q文件
  # -----------------------------
  q_data <- read.table(
    q_file,
    header = FALSE,
    stringsAsFactors = FALSE
  )
  
  # -----------------------------
  # 检查样本数量
  # -----------------------------
  if (nrow(q_data) != length(fam[, 2])) {
    stop(
      paste0(
        "K=", k,
        " 的Q文件行数与FAM文件样本数不一致！\n",
        "Q文件：", nrow(q_data),
        " 行\n",
        "FAM文件：", length(fam[, 2]),
        " 个样本"
      )
    )
  }
  
  # -----------------------------
  # 给Q文件列命名
  # -----------------------------
  colnames(q_data) <- paste0(
    "Ancestry_",
    1:k
  )
  
  # -----------------------------
  # 加入Sample
  # 注意：这里使用原始fam顺序
  # -----------------------------
  q_data$Sample <- fam[, 2]
  
  # -----------------------------
  # 修改Sample名称
  # -----------------------------
  q_data$Sample[grepl("^\\d+", q_data$Sample)] <-
    paste0(
      "NC_",
      q_data$Sample[grepl("^\\d+", q_data$Sample)]
    )
  
  q_data$Sample[grepl("^AA\\d+", q_data$Sample)] <-
    paste0(
      "C_",
      q_data$Sample[grepl("^AA\\d+", q_data$Sample)]
    )
  
  # -----------------------------
  # 删除指定样本
  # -----------------------------
  q_data <- q_data[
    !q_data$Sample %in% remove_samples,
    ,
    drop = FALSE
  ]
  
  # -----------------------------
  # 转成长格式
  # -----------------------------
  q_melt <- reshape2::melt(
    q_data,
    id.vars = "Sample",
    variable.name = "variable",
    value.name = "value"
  )
  
  # -----------------------------
  # 与genotype_df匹配
  # -----------------------------
  q_melt <- merge(
    q_melt,
    genotype_df,
    by.x = "Sample",
    by.y = "SampleID",
    all.x = FALSE,
    all.y = FALSE
  )
  
  # -----------------------------
  # 统一SampleID名称
  # -----------------------------
  colnames(q_melt)[1] <- "SampleID"
  
  # -----------------------------
  # 非常关键：
  # 按进化树tip顺序排列
  # -----------------------------
  q_melt$SampleID <- factor(
    q_melt$SampleID,
    levels = tree_rooted$tip.label
  )
  
  # -----------------------------
  # 保存
  # -----------------------------
  q_melt_list[[paste0("K", k)]] <- q_melt
}


# ============================================================
# 6. 检查每个K是否都成功匹配
# ============================================================

cat("\n==============================\n")
cat("Admixture data summary\n")
cat("==============================\n")

for (k in k_values) {
  
  dat <- q_melt_list[[paste0("K", k)]]
  
  cat(
    "K=", k,
    ": ",
    length(unique(dat$SampleID)),
    " samples\n",
    sep = ""
  )
}


# ============================================================
# 7. 先画基础进化树
# ============================================================

p_tree <- ggtree(tree_rooted,  
                 layout = "rectangular",  
                 #branch.length = "none", 
                 lwd = 1) +
  #theme_tree2() +
  geom_tippoint(aes(color = tip_genotype[label]), size = 2, shape = 19) +
  geom_tree(aes(color = tip_genotype[label]), size = 1) +
  #geom_tiplab(size = 1.5, offset = 0.5) +  # 标签偏移，避免和条形图重叠
  scale_color_manual(values = genotype_colors, name = "Genotype") +
  theme_tree2(legend.position = "none")
p_tree

# ============================================================
# 8. 依次添加9个Admixture panel
# ============================================================

p_final <- p_tree


for (k in k_values) {
  
  cat(
    "Adding Admixture panel K =",
    k,
    "\n"
  )
  
  dat <- q_melt_list[[paste0("K", k)]]
  
  p_final <- facet_plot(
    p_final,
    
    # 每个K使用不同panel名称
    panel = paste0(
      "K=",
      k
    ),
    
    # 当前K的数据
    data = dat,
    
    # 条形图
    geom = geom_bar,
    
    # Admixture比例
    aes(
      x = value,
      fill = variable
    ),
    
    # 柱子的高度
    stat = "identity",
    
    # 堆叠
    position = "stack",
    
    # 每个样本一个bar
    width = 1,
    
    # 黑色边框
    color = "black",
    
    # 边框非常细
    linewidth = 0.05,
    
    # 关键：横向bar
    orientation = "y"
  )
}


# ============================================================
# 9. 统一Admixture颜色
# ============================================================

p_final <- p_final +
  scale_fill_manual(
    values = ancestry_colors,
    name = "Ancestry"
  )


# ============================================================
# 10. 美化
# ============================================================

p_final <- p_final +
  theme(
    
    # -------------------------
    # 图例
    # -------------------------
    legend.position = "right",
    
    legend.key.size = unit(
      0.45,
      "cm"
    ),
    
    legend.text = element_text(
      size = 8
    ),
    
    legend.title = element_text(
      size = 9,
      face = "bold"
    ),
    
    # -------------------------
    # Panel标题
    # -------------------------
    strip.background = element_blank(),
    
    strip.text = element_text(
      size = 9,
      face = "bold"
    ),
    
    # -------------------------
    # panel之间距离
    # -------------------------
    panel.spacing = unit(
      0.08,
      "cm"
    ),
    
    # -------------------------
    # 去掉背景
    # -------------------------
    panel.background = element_blank(),
    
    # -------------------------
    # 去掉网格
    # -------------------------
    panel.grid = element_blank(),
    
    # -------------------------
    # 坐标轴文字
    # -------------------------
    axis.text = element_text(
      size = 7
    ),
    
    axis.title = element_text(
      size = 8
    )
  )



# 保存高分辨率图形（根据样本数调整尺寸）
ggsave(
  filename = "FigS2_ASFV_Tree+Admixture_Final.pdf",
  plot = p_final,
  width = 13,  # 宽度：树(8) + 条形图(4)
  height = 10, # 高度适配364个样本
  dpi = 300,
  device = "pdf",
  path = "C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/"
)


# ============================================================
# Admixture CV Error 折线图
# ============================================================
library(ggplot2)

# -----------------------------
# 数据
# -----------------------------
cv_data <- data.frame(
  K = 1:10,
  CV_error = c(
    0.45990,  # K=1
    0.30727,  # K=2
    0.08704,  # K=3
    0.04441,  # K=4
    0.04000,  # K=5
    0.03940,  # K=6
    0.03216,  # K=7
    0.03134,  # K=8 (最低)
    0.04119,  # K=9
    0.03969   # K=10
  )
)

# -----------------------------
# 找出最低CV error的K
# -----------------------------
best_k <- cv_data$K[which.min(cv_data$CV_error)]
best_cv <- min(cv_data$CV_error)

# -----------------------------
# 绘图
# -----------------------------
p_cv <- ggplot(cv_data, aes(x = K, y = CV_error)) +

# 0.05 参考横线
geom_hline(
  yintercept = 0.05,
  linetype = "dashed",
  color = "#7F8C8D",
  linewidth = 0.6
) +

  # 折线
  geom_line(
    color = "#2C3E50",
    linewidth = 0.8
  ) +
  
  # 数据点
  geom_point(
    shape = 21,
    fill = "white",
    color = "#2C3E50",
    size = 3,
    stroke = 1
  ) +
  
  # 高亮最低点 K=8 → 填充蓝色
  geom_point(
    data = subset(cv_data, K == best_k),
    aes(x = K, y = CV_error),
    shape = 21,
    fill = "#2980B9",
    color = "#2980B9",
    size = 4,
    stroke = 1
  ) +
  
  # 最低点标注 → 蓝色
  annotate(
    "text",
    x = best_k,
    y = best_cv + 0.015,
    label = paste0("K = ", best_k),
    #label = paste0("K = ", best_k, "\nCV = ", sprintf("%.5f", best_cv)),
    color = "#2980B9",
    size = 3.5,
    fontface = "bold",
    hjust = 0.5
  ) +
  
  # K=4 数据点高亮 → 填充红色
  geom_point(
    data = subset(cv_data, K == 4),
    aes(x = K, y = CV_error),
    shape = 21,
    fill = "#E74C3C",
    color = "#E74C3C",
    size = 4,
    stroke = 1.2
  ) +
  
  # K=4 标注 → 红色
  annotate(
    "text",
    x = 4,
    y = 0.07,
    label = "K = 4\n(genotype-based selection)",
    color = "#E74C3C",
    size = 3.5,
    fontface = "bold",
    hjust = 0.5
  ) +
  
  # 坐标轴
  scale_x_continuous(
    breaks = 1:10,
    limits = c(0.5, 10.5)
  ) +
  
  scale_y_continuous(
    breaks = c(seq(0,0.5,0.1),0.05),
    limits = c(0, 0.5),
    expand = c(0, 0)
  ) +
  
  # 标签
  labs(
    x = "Number of ancestral populations (K)",
    y = "Cross-validation error"
  ) +
  
  # 主题
  theme_classic() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10, color = "black"),
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6),
    plot.margin = margin(10, 15, 10, 10)
  )

print(p_cv)

# 保存
ggsave("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/FigS_CV_error_plot.pdf", p_cv, width = 13, height = 5, dpi = 300)
# ggsave("CV_error_plot.png", p_cv, width = 6, height = 4, dpi = 300)





library(vcfR)
library(ggplot2)

vcf <- read.vcfR("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/recombination/asfv_snps.vcf")
gt_matrix <- extract.gt(vcf, as.numeric = TRUE) #0=ref, 1=alt
#矩阵：行=SNP位点；列=样本；转置：行=样本
geno_t <- t(gt_matrix)

#计算协方差PCA
pca <- prcomp(geno_t, scale. = FALSE)
pc_df <- as.data.frame(pca$x[,c(1,2)])
colnames(pc_df) <- c("PC1","PC2")
pc_df$SampleID <- rownames(pc_df)

#方差解释
var <- round(pca$sdev^2/sum(pca$sdev^2)*100, 2)

#合并分组信息

pc_df <- merge(pc_df, genotype_df, by="SampleID")
annotation_df$SampleID=rownames(annotation_df)
pc_df <- merge(pc_df, annotation_df[,c("SampleID","Continent")], by="SampleID")

library(vegan)
library(pairwiseAdonis)
# 欧式距离，基于PC1+PC2

dist_eu <- vegdist(data.table(pc_df)[, .(PC1, PC2)], method = "euclidean")
set.seed(123)
perma_global <- adonis2(dist_eu ~ Genotype, data = pc_df, permutations = 999)
perma_pair <- pairwise.adonis2(dist_eu ~ Genotype, data = pc_df, permutations = 999, p.adjust.m = "fdr")
disp <- betadisper(dist_eu, group = pc_df$Genotype)
anova(disp)

# 提取标注文本
r2_val <- round(perma_global$R2[1], 3)
p_val  <- round(perma_global$`Pr(>F)`[1], 3)
stat_label <- paste0("PERMANOVA\nGenotype: R² = ", r2_val, ", P = ", p_val)


ggplot(pc_df[,1:4], aes(x=PC1,y=PC2,color=Genotype))+
  geom_point(size=3.5, alpha=0.8)+
  scale_color_manual(values=genotype_colors)+
  ggrepel::geom_text_repel(data =filter(pc_df,(Genotype==1 & PC2<25) | (Genotype==22 & PC2>0) | (Genotype=="1/2" & PC1>20)), aes(x=PC1,y=PC2, label = SampleID,color=Genotype))+
  # stat_ellipse(
  #   data = pc_df[pc_df$Genotype %in% c("1",  "9"),],
  #   level = 0.95,
  #   linewidth = 0.8
  # ) +
  annotate("text", x = -0.058, y = -45, label = stat_label, size = 3.5) +
  #annotate("text", x = -0.058, y = -50, label = stat_label2, size = 3.5) +
  labs(x=paste0("PC1 (",var[1],"%)"),
       y=paste0("PC2 (",var[2],"%)"))+
  theme_classic()+theme(panel.grid=element_blank())

ggsave("FigS1B_PCA_plot.pdf", p_cv, width = 6, height = 4, dpi = 300)


# perma_global <- adonis2(dist_eu ~ Continent, data = pc_df, permutations = 999)
# perma_pair <- pairwise.adonis2(dist_eu ~ Continent, data = pc_df, permutations = 999, p.adjust.m = "fdr")
# disp <- betadisper(dist_eu, group = pc_df$Continent)
# anova(disp)
# 
# # 提取标注文本
# r2_val <- round(perma_global$R2[1], 3)
# p_val  <- round(perma_global$`Pr(>F)`[1], 3)
# stat_label2 <- paste0("Continent: R² = ", r2_val, ", P = ", p_val)
# 

# ggplot(pc_df, aes(x=PC1,y=PC2,color=Continent))+
#   geom_point(size=3.5, alpha=0.8)+
#   scale_color_manual(values=continent_colors)+
#   # stat_ellipse(
#   #   data = pc_df[pc_df$Genotype %in% c("1",  "9"),],
#   #   level = 0.95,
#   #   linewidth = 0.8
#   # ) +
#   annotate("text", x = -0.058, y = -45, label = stat_label, size = 3.5) +
#   annotate("text", x = -0.058, y = -50, label = stat_label2, size = 3.5) +
#   labs(x=paste0("PC1 (",var[1],"%)"),
#        y=paste0("PC2 (",var[2],"%)"))+
#   theme_classic()+theme(panel.grid=element_blank())

library(dplyr)

# 1. PCA 筛选出的样本（你的条件）
highlight_samples <- filter(
  pc_df,
  (Genotype == 1    & PC2 < 25) |
    (Genotype == 22   & PC2 > 0)  |
    (Genotype == "1/2" & PC1 > 20) #|
   # (Genotype == "9"  & PC1 < 60)
)$SampleID

# 2. 只保留确实在树上的（避免名不对齐被静默忽略）
highlight_samples <- highlight_samples[highlight_samples %in% tree_rooted$tip.label]

# 3. 绘图 + 高亮
p_tree <- ggtree(tree_rooted,
                 layout = "rectangular",
                 panel.width = 2,
                 lwd = 1) +
  geom_tippoint(aes(color = tip_genotype[label]), size = 2, shape = 19) +
  geom_tree(aes(color = tip_genotype[label]), size = 1) +
  scale_color_manual(values = genotype_colors, name = "Genotype") +
  theme_tree2(legend.position = "none") +
  # 关键：给筛选样本套一个红色粗圈，不遮挡原有 genotype 点
  geom_point2(
    aes(subset = (label %in% highlight_samples)),
    shape = 21,
    fill = NA,
    color = "black",
    size = 4,
    stroke = 1.3
  )
p_tree
