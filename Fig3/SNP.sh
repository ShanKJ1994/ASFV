###################################################################################################################################################

#nextstrain for linegae

###################################################################################################################################################
source /lustre/user/lulab/shankj/.bashrc
cd /gpfs2/shankj/ASFV/20251010/SNP/
treetime ancestral --tree /gpfs2/shankj/ASFV/20251010/Recombination/0RDP/ASFV_all_RDP/PS/RDP4_methods_2_rmRE_includ-G1-G2-G9-G20.contree.root.txt --aln  /gpfs2/shankj/ASFV/20251010/Recombination/0RDP/ASFV_all_RDP/trimmed20260803.fasta  --outdir ./ --gtr infer 

/lustre/user/lulab/shankj/miniconda3/envs/GO/bin/Rscript /gpfs2/shankj/ASFV/20251010/SNP/Spilt_Genotype.R

cd /gpfs2/shankj/ASFV/20251010/SNP/split_by_geno/
mkdir ./G2
treetime clock  --dates /gpfs2/shankj/ASFV/20251010/SNP/meta.tsv --tree /gpfs2/shankj/ASFV/20251010/SNP/split_by_geno/G2_tree.nwk --aln  /gpfs2/shankj/ASFV/20251010/SNP/split_by_geno/G2_aligned.fasta  --outdir ./G2 --keep-root
mkdir ./G9
treetime clock  --dates /gpfs2/shankj/ASFV/20251010/SNP/meta.tsv --tree /gpfs2/shankj/ASFV/20251010/SNP/split_by_geno/G9_tree.nwk --aln  /gpfs2/shankj/ASFV/20251010/SNP/split_by_geno/G9_aligned.fasta  --outdir ./G9 --keep-root
mkdir ./G1
treetime clock  --dates  /gpfs2/shankj/ASFV/20251010/SNP/meta.tsv --tree /gpfs2/shankj/ASFV/20251010/SNP/split_by_geno/G1_tree.nwk --aln  /gpfs2/shankj/ASFV/20251010/SNP/split_by_geno/G1_aligned.fasta  --outdir ./G1 --keep-root

perl AncestorSNP.pl

perl CommonAncestorSNP.pl


