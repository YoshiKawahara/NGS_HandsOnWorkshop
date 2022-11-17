##########
# Bioconductorのパッケージ（ballgown, genefilter, DESeq2）のインストール
# 途中、Update all/some/none? [a/s/n]:などと聞かれたら"a"でupdateしておく
install.packages("BiocManager")
BiocManager::install(c("ballgown", "genefilter", "DESeq2"))

# Rのパッケージ（tidyverse, pheatmap）のインストール
install.packages("tidyverse")
install.packages("pheatmap")

# インストールしたパッケージの読み込み
library(ballgown)
library(genefilter)
library(DESeq2)
library(tidyverse)
library(pheatmap)

##########
##### 昼夜で遺伝子発現が有意に変化する遺伝子を抽出
##### ballgownによる遺伝子発現変動解析
# サンプル情報の読み込み
pheno_data <- read.csv("phenodata.csv")

# サンプル情報を表示
pheno_data

# ballgownの入力データの読み込み（StringTieの結果ディレクトリを指定）
bg <- ballgown(dataDir="ballgown", samplePattern="rice", pData=pheno_data)
bg

# サンプル間の分散が1以上のものに絞り込む、つまり発現変動が小さい遺伝子を除く
bg_filtered <- subset(bg, "rowVars(texpr(bg))>=1")
bg_filtered

# 転写産物レベルでサンプル（Day/Night）間での遺伝子発現変動を検定し、p-valueやq-valueを計算する
results_transcripts <- stattest(bg_filtered, feature="transcript", covariate="cond", getFC=TRUE, meas="FPKM")

# head関数で先頭の何行かを表示する（デフォルトは6行）
head(results_transcripts)

# 転写産物レベルの結果に遺伝子名、遺伝子ID、転写産物IDの列を追加する。
results_transcripts = data.frame(geneNames=ballgown::geneNames(bg_filtered),
                                 geneIDs=ballgown::geneIDs(bg_filtered),
                                 transcriptIDs=ballgown::transcriptNames(bg_filtered),
                                 results_transcripts)

head(results_transcripts)

# 転写産物レベルでの検定結果をp-valueの小さい順にソートする
results_transcripts = arrange(results_transcripts, pval)

# 有意（q-value<0.01）な発現変動を示す遺伝子を抽出する
subset(results_transcripts, results_transcripts$qval<0.01)


##########
##### ballgownによる発現データの可視化
# 全サンプルの発現量（FPKM）の平均が1以上の転写産物を取り出す
bg.expressed <- subset(bg, "rowMeans(texpr(bg))>=1")
bg.expressed

# 発現している転写産物のFPKM値（底が2の対数）を抽出
log2fpkm <- log2(texpr(bg.expressed, meas="FPKM") + 0.01)
head(log2fpkm)

# log2fpkmをggplotに適したフォーマットに整形する
log2fpkm.long <- as.data.frame(log2fpkm) %>%
  tidyr::gather(Dataset, FPKM)
head(log2fpkm.long)

# ggplotで各サンプルごとの遺伝子発現量の分布を示すboxplotを描画
ggplot(log2fpkm.long, aes(y=FPKM, x=Dataset)) + geom_boxplot() +
 theme(axis.text.x = element_text(angle=90, hjust=0, vjust=.5)) +
 xlab("") + ylab("log2(FPKM+0.01)")

# LHY遺伝子(ballgown id:19641)の発現量を抽出し、フォーマットを整形
log2fpkm["19641", ]
log2fpkm.LHY.long <- as.data.frame(log2fpkm["19641", ]) %>%
  tidyr::gather(Dataset, FPKM)

# サンプル情報(Day or Night)を付加
log2fpkm.LHY.long <- data.frame(log2fpkm.LHY.long, Condition=c(rep("Day", 4), rep("Night", 4)))
log2fpkm.LHY.long

# LHY遺伝子の昼と夜の発現量の分布を描画
ggplot(log2fpkm.LHY.long, aes(x=Condition, y=FPKM, fill=Condition)) +
  geom_boxplot() + geom_point() +
  scale_fill_manual(values=c("orange","blue")) +
  xlab("") + ylab("log2(FPKM+0.01)") + theme(legend.position = "bottom")

# ballgownのplotTranscripts関数によって、SIGA遺伝子(ballgown id: 31942)の発現量を転写産物構造と共に描画
plotTranscripts(ballgown::geneIDs(bg)[31942], bg, main=c('SIGA'), sample=c('rice_D_rep1','rice_D_rep2','rice_D_rep3','rice_D_rep4','rice_N_rep1','rice_N_rep2','rice_N_rep3','rice_N_rep4'))


##########
##### DESeq2による遺伝子発現変動解析
# カウントデータとサンプルデータを用意する
count_matrix <- as.matrix(read.csv("gene_count_matrix.csv", sep=",",row.names="gene_id"))
sample_info <- data.frame(CONDITION=as.factor(c("Day","Day","Day","Day","Night","Night","Night","Night")))
rownames(sample_info) <- as.factor(c("rice_D_rep1","rice_D_rep2","rice_D_rep3","rice_D_rep4","rice_N_rep1","rice_N_rep2","rice_N_rep3","rice_N_rep4"))

# DESeqの解析に用いるオブジェクトを作成
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                               colData = sample_info,
                               design =~ CONDITION)

# 可視化する前にカウントデータを変換する(variance stabilizing transformations)
vsd <- vst(dds, blind=FALSE)
head(assay(vsd),10)

# ヒートマップによってサンプル間の遺伝子発現の相関を可視化する
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(vsd)
colnames(sampleDistMatrix) <- NULL
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists)

# PCAによってサンプル間の遺伝子発現の相関を可視化する
plotPCA(vsd, intgroup=c("CONDITION"))

# DESeq2によるDEG解析
featureData <- data.frame(gene=rownames(count_matrix))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)

# 全サンプルの合計カウント数が10未満の低発現遺伝子を除く
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

dds$CONDITION <- factor(dds$CONDITION, levels=c("Day","Night"))

# DayとNightの条件間で発現比較を行う
dds <- DESeq(dds)
res <- results(dds, contrast=c("CONDITION","Day","Night"))

# 結果をP値でソートし、CSVファイルとして出力する。
resOrdered <- res[order(res$pvalue),]
resOrdered

write.table(as.data.frame(resOrdered), file="result.tsv", quote=F, sep="\t", row.names=T, col.names=NA)
