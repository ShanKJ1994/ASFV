#!/usr/bin/perl
use strict;
use warnings;
#use Text::CSV;

#===================== 【配置区，你修改这里】=====================
my $aln_fasta  = "/gpfs2/shankj/ASFV/20251010/SNP/ancestral_sequences.fasta";     # 比对好的fasta
my $meta_tsv   = "meta.tsv";      # metadata表格 tab分隔
my $out_tsv    = "Common_snp_call_result.tsv";

# genotype => 祖先序列ID（把这里替换成你实际的祖先fasta ID）
my %ancestor_map = (
    1 => "NODE_0000012",
    2 => "NODE_0000012",
    9 => "NODE_0000002",
);
#================================================================

# 1.读取metadata，建立 name -> genotype
my %name2geno;
open(my $mfh,"<",$meta_tsv) or die $!;
my $meta_header = <$mfh>;
while(my $l=<$mfh>){
    chomp $l;
    next if $l=~/^\s*$/;
    my @col = split /\t/, $l;
    my ($name,$genotype) = ($col[0],$col[1]);
    $name2geno{$name} = $genotype;
}
close $mfh;

# 2.读取比对fasta：ID => 比对序列（保留gap）
my %aln;
open(my $afh,"<",$aln_fasta) or die $!;
my $cur_id;
my $cur_seq;
while(my $l=<$afh>){
    chomp $l;
    if($l=~/^>(\S+)/){
        if(defined $cur_id){
            $aln{$cur_id} = uc($cur_seq);
        }
        $cur_id = $1;
        $cur_seq = "";
    }else{
        $l =~ s/\s+//g;
        $cur_seq .= $l;
    }
}
if(defined $cur_id){
    $aln{$cur_id}=uc($cur_seq);
}
close $afh;

# 3.预读祖先序列
my %ancestor_seq;
foreach my $g (keys %ancestor_map){
    my $aid = $ancestor_map{$g};
    die "Error: ancestor id $aid not found in aligned fasta!\n" unless exists $aln{$aid};
    $ancestor_seq{$g} = $aln{$aid};
}

# 获取比对总长度
my $aln_len;
{
    my $tmp_id = (keys %aln)[0];
    $aln_len = length($aln{$tmp_id});
}

# 4.输出结果表头
open(my $ofh,">",$out_tsv) or die $!;
print $ofh join("\t",
    qw/aln_pos genotype ancestor_id ancestor_base sample_id sample_base is_snp motif/
),"\n";

# 遍历每个比对位置 aln_pos 1‑based
for(my $pos=0; $pos < $aln_len; $pos++){
    my $aln_pos = $pos + 1;

    # 遍历所有样本（metadata中的name）
    foreach my $sample_id (keys %name2geno){
        next unless exists $aln{$sample_id};
        my $g = $name2geno{$sample_id};
        next unless exists $ancestor_map{$g}; # 只处理genotype=1/2/9

        my $anc_id  = $ancestor_map{$g};
        my $anc_base = substr($ancestor_seq{$g}, $pos, 1);
        my $sam_base = substr($aln{$sample_id}, $pos, 1);

        # 过滤：祖先或者样本为 N / - ，跳过这个位点不输出
        next if $anc_base =~ /[N\-]/i;
        next if $sam_base =~ /[N\-]/i;

        my $is_snp = ($anc_base ne $sam_base) ? 1 : 0;
        if ($is_snp ==1) {
        
        # 【新增】提取Motif：以当前pos为中心，往前取5，往后取5，总长度11
        my $motif = substr($ancestor_seq{$g}, $pos - 5, 11);

        print $ofh join("\t",
            $aln_pos,
            $g,
            $anc_id,
            $anc_base,
            $sample_id,
            $sam_base,
            $is_snp,
            $motif
        ),"\n";
		}
    }
}
close $ofh;

print "Done! output: $out_tsv\n";
print "Total alignment length: $aln_len\n";
