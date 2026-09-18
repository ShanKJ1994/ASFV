#==================== 配置（直接改路径，不要动其他逻辑）====================
library(ape)
library(seqinr)

# 文件路径
meta_file  <- "/gpfs2/shankj/ASFV/20251010/SNP/meta.tsv"
aln_file   <- "/gpfs2/shankj/ASFV/20251010/Recombination/0RDP/ASFV_all_RDP/trimmed20260803.fasta"
tree_file  <- "/gpfs2/shankj/ASFV/20251010/Recombination/0RDP/ASFV_all_RDP/PS/RDP4_methods_2_rmRE_includ-G1-G2-G9-G20.contree.root.txt"

# 输出目录（建议直接建一个子文件夹放拆分结果）
outdir     <- "/gpfs2/shankj/ASFV/20251010/SNP/split_by_geno"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
#===========================================================================

# 1.读meta（tab分隔）
meta <- read.delim(meta_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# 确保name列和fasta/tree ID完全一致！
head(meta)

# 2.读取比对fasta
aln <- read.fasta(aln_file, seqtype = "DNA", as.string = TRUE)
all_aln_ids <- names(aln)
# 3.读取进化树（ape读rooted树）
tr <- read.tree(tree_file)
all_tree_tips <- tr$tip.label

# 校验：meta里的name必须同时存在于fasta和树
meta <- meta[grepl("\\d+",meta$date) & meta$name %in% all_aln_ids & meta$name %in% all_tree_tips & !meta$name %in% c("OR944085.1", "OR944086.1", "OR944087.1", "OR944088.1", "OR944089.1", "OR944090.1", "OR944091.1"), ]


# 拿到所有unique genotype
geno_list <- unique(meta$genotype)
cat("发现基因型：", geno_list, "\n")

# 循环按genotype拆分
for(g in geno_list){
  # 当前基因型的样本ID
  keep_ids <- meta$name[meta$genotype == g]
  cat("\n===== Genotype", g, "｜共", length(keep_ids), "条 =====\n")
  
  # ---------- 裁剪fasta ----------
  sub_aln <- aln[keep_ids]
  fa_out <- file.path(outdir, paste0("G",g,"_aligned.fasta"))
  write.fasta(sequences = sub_aln, names = names(sub_aln), file.out = fa_out)
  
  # ---------- 裁剪进化树（drop其他tip） ----------
  # ape::drop.tip：保留keep_ids，删掉其余
  sub_tr <- drop.tip(tr, setdiff(all_tree_tips, keep_ids))
  tree_out <- file.path(outdir, paste0("G",g,"_tree.nwk"))
  write.tree(sub_tr, file = tree_out)
  
  cat("输出：", fa_out, "\n")
  cat("输出：", tree_out, "\n")
}

cat("\n✅ 全部拆分完成！结果目录：", outdir, "\n")
