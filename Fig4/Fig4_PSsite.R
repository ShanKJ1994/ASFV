
#############################################

#Pos  selection

#############################################
library(data.table)
library(UpSetR)


Site=fread("C:/Users/Ke-jia Shan/Desktop/ASFV_Figs/positive_selection_sites.csv",header = T)

#mid=filter(merge(Site,mid,by="gene"),site==min,significant=="TRUE",gene %in% filter(test,SeqNum>280)$Gene)  %>% group_by(gene_site) %>% dplyr::summarise(count=n())

# 长表转宽表：每个gene_site一行，每个algorithm一列，TRUE/FALSE
site_wide <- dcast(Site, gene_site ~ algorithm, value.var = "significant", fun.aggregate = any)

# UpSetR 需要data.frame，布尔转为整数1/0
upset_input <- site_wide[, !"gene_site"] %>% as.data.frame()
upset_input[] <- lapply(upset_input, as.integer)
rownames(upset_input) <- site_wide$gene_site

# ----------------------
# 绘制UpSet图
# ----------------------
upset(upset_input,
      sets = c("slac","fel","fubar","meme"),   # 你的算法集合
      #keep.order = TRUE,
      order.by = "freq",        # 交集柱子按频数从大到小排序
      sets.bar.color = c("#FB9A99","#4088B4","#E08251","#9D71A6"),
      mainbar.y.label = "Number of sites under positive selection",
      sets.x.label = "Sites supported by single algorithm",
      text.scale = c(1.2, 1.2, 1, 1, 1.2, 1)
      )

ggsave(
  filename = "C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/recombination/0RDP/ASFV_all_RDP/SCG/Upset.pdf",
  width = 7,
  height = 4,
  device = "pdf"
)

#data.table(upset_input)[fel+fubar+ meme+ slac>0]

library(ggplot2)

# 氨基酸配色映射
aa_color <- c(
  P = "#66C2A5",
  D = "#8DA0CB",
  E = "#E78AC3",
  R = "#A6D854",
  K = "#FFD92F",
  S = "#E5C494",
  G = "#A65628",
  N = "#A6CEE3",
  A = "#CAB2D6",
  V = "#FC8D62",
  X = "#B3B3B3",
  Y = "#4DAF4A"
)

df <- data.frame(
  AA = names(aa_color),
  col = aa_color
)

# 色块图例预览图
p <- ggplot(df, aes(x = AA, y = 1, fill = AA)) +
  geom_tile(colour = "black", linewidth = 0.3, height = 1) +
  geom_text(aes(label = AA), size = 4.5, colour = "black") +
  scale_fill_manual(values = aa_color) +
  labs(title = "AA color legend") +
  theme_bw() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  )

print(p)

# 可直接复制用于ggplot绘图的配色代码
cat("scale_color_manual(values = aa_color)\nscale_fill_manual(values = aa_color)\n\n")
dput(aa_color)


if(!require(tidyverse)) install.packages("tidyverse")
library(tidyverse)

# ========== 1. 原始数据 ==========
df_raw <- tibble(
  Gene = c("VF199", "P54", "RPB1", "VF145", "VF205", "VF602", "VF717"),
  Domain_range = c("150-170", "32-52", NA, NA, "10-75", NA, NA),
  Domain = c("Transmembrane", "Transmembrane", NA, NA, "CoiledCoil", NA, NA),
  Gene_length = c(199, 183, 1450, 145, 205, 530, 717),
  MutationSite = c("85,125,127", "139", "18", "145", "192", "41", "392,613")
)

# ========== 2. 拆分突变位点 ==========
mut_df <- df_raw %>%
  separate_rows(MutationSite, sep = ",") %>%
  mutate(MutationSite = as.integer(MutationSite)) %>%
  filter(!is.na(MutationSite))

# ========== 3. 拆分功能域起止坐标 ==========
domain_df <- df_raw %>%
  filter(!is.na(Domain_range)) %>%
  separate(Domain_range, into = c("start", "end"), sep = "-", convert = TRUE) %>%
  select(Gene, start, end, Domain)

# ========== 4. 基因全长背景框 ==========
gene_full_df <- df_raw %>%
  select(Gene, end = Gene_length) %>%
  mutate(start = 0)

# ========== 5. 循环：每个基因单独绘图，x轴范围自适应基因长度 ==========
gene_list <- unique(df_raw$Gene)

for (g in gene_list) {
  # 筛选当前基因数据
  g_full <- filter(gene_full_df, Gene == g)
  g_domain <- filter(domain_df, Gene == g)
  g_mut <- filter(mut_df, Gene == g)
  max_x <- g_full$end
  
  p <- ggplot() +
    # 基因全长浅米色背景
    geom_rect(
      data = g_full,
      aes(xmin = start, xmax = end, ymin = -0.25, ymax = 0.25),
      fill = "#D2C8B4", colour = NA
    ) +
    # 功能域色块
    geom_rect(
      data = g_domain,
      aes(xmin = start, xmax = end, ymin = -0.25, ymax = 0.25, fill = Domain),
      colour = NA
    ) +
    # 突变位点竖线（贯穿蛋白条）
    geom_segment(
      data = g_mut,
      aes(x = MutationSite, xend = MutationSite,
          y = -0.32, yend = 0.32),
      colour = "black",
      linewidth = 0.8
    ) +
    scale_fill_manual(values = c(
      "Transmembrane" = "#FFD700",
      "CoiledCoil" = "#66C2A5"
    )) +
    scale_x_continuous(limits = c(0, max_x), expand = c(0.02, 0)) +
    scale_y_continuous(limits = c(-0.6, 0.6), breaks = 0, labels = g) +
    labs(x = "Amino acid position", y = "") +
    theme_classic() +
    theme(
      panel.grid = element_blank(),
      legend.position = "bottom",
      plot.margin = margin(10,10,10,10)
    )
  
  # 保存PDF，宽度按基因长度动态缩放（可选固定宽度：width=8）
  ggsave(
    filename = paste0("C:/Users/Ke-jia Shan/Desktop/博后/博后课题/ASFV/recombination/0RDP/ASFV_all_RDP/SCG/",g, "_protein_domain.pdf"),
    plot = p,
    width = ifelse(max_x > 800, 12, 6),
    height = 3,
    device = "pdf"
  )
}

#正选择基因
# gene	branch	omega	pvalue
# IKBL	PP355086	1000	0
# IKBL	Genoype2_ancestor	1000	0
# P30	PP737715	1000	0.00102
# VF145	OL322099	1000	0.00278

