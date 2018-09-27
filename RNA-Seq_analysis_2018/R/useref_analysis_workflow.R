### ballgownやtidyverse等のパッケージのインストールと読み込み
# 途中、Update all/some/none? [a/s/n]:などと聞かれたら"a"でupdateしておく
source("https://bioconductor.org/biocLite.R")
biocLite("ballgown")
library(ballgown)

biocLite("genefilter")
library(genefilter)

install.packages("tidyverse")
library(tidyverse)

##########
##### ballgownを用いて昼夜で遺伝子発現が有意に変化する遺伝子を抽出
### サンプル情報の読み込み
pheno_data <- read.csv("phenodata.csv")
# サンプル情報を表示
pheno_data

### ballgownの入力データの読み込み（StringTieの結果ディレクトリを指定）
bg <- ballgown(dataDir="ballgown", samplePattern="rice", pData=pheno_data)
bg

### サンプル間の分散が1以上のものに絞り込む、つまり発現変動が小さい遺伝子を除く
bg_filtered <- subset(bg, "rowVars(texpr(bg))>=1")
bg_filtered

### 転写産物レベルでサンプル（Day/Night）間での遺伝子発現変動を検定し、p-valueやq-valueを計算する
results_transcripts <- stattest(bg_filtered, feature="transcript", covariate="cond", getFC=TRUE, meas="FPKM")

### head関数で先頭の何行かを表示する（デフォルトは6行）
head(results_transcripts)

### 転写産物レベルの結果に遺伝子名、遺伝子ID、転写産物IDの列を追加する。
results_transcripts = data.frame(geneNames=ballgown::geneNames(bg_filtered), geneIDs=ballgown::geneIDs(bg_filtered), transcriptIDs=ballgown::transcriptNames(bg_filtered), results_transcripts)

head(results_transcripts)

### 転写産物レベルでの検定結果をp-valueの小さい順にソートする
results_transcripts = arrange(results_transcripts, pval)

### 適当なq-valueの閾値で切った遺伝子を抽出する。
subset(results_transcripts, results_transcripts$qval<0.01)


##########
##### Visualization
### 全サンプルの発現量（FPKM）の平均が1以上の転写産物を取り出す
bg.expressed <- subset(bg, "rowMeans(texpr(bg))>=1")
bg.expressed

### 発現している転写産物のFPKM値（底が2の対数）を抽出
log2fpkm <- log2(texpr(bg.expressed, meas="FPKM") + 0.01)
head(log2fpkm)

### log2fpkmをggplotに適したフォーマットに整形する
log2fpkm.long <- as.data.frame(log2fpkm) %>%
  tidyr::gather(Dataset, FPKM)
head(log2fpkm.long)

### ggplotで各サンプルごとの遺伝子発現量の分布を示すboxplotを描画
ggplot(log2fpkm.long, aes(y=FPKM, x=Dataset)) + geom_boxplot() +
 theme(axis.text.x = element_text(angle=90, hjust=0, vjust=.5)) +
 xlab("") + ylab("log2(FPKM+0.01)")

### ggplotで各サンプルごとの遺伝子発現量の分布を示すviolinplotを描画
ggplot(log2fpkm.long, aes(x=Dataset, y=FPKM)) + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  theme(axis.text.x = element_text(angle=90, hjust=0, vjust=.5)) +
  xlab("") + ylab("log2(FPKM+0.01)")

### LHY遺伝子(ballgown id:19641)の発現量を抽出し、フォーマットを整形
log2fpkm["19641", ]
log2fpkm.LHY.long <- as.data.frame(log2fpkm["19641", ]) %>%
  tidyr::gather(Dataset, FPKM)

### サンプル情報(Day or Night)を付加
log2fpkm.LHY.long <- data.frame(log2fpkm.LHY.long, Condition=c(rep("Day", 4), rep("Night", 4)))
log2fpkm.LHY.long

### LHY遺伝子の昼と夜の発現量の分布を描画
ggplot(log2fpkm.LHY.long, aes(x=Condition, y=FPKM, fill=Condition)) +
  geom_boxplot() + geom_point() +
  xlab("") + ylab("log2(FPKM+0.01)") + theme(legend.position = "bottom")

### ballgownの描画機能を使い、SIGA遺伝子(ballgown id: 31942)の発現量を転写産物構造と共に描画
plotTranscripts(ballgown::geneIDs(bg)[31942], bg, main=c('SIGA'), sample=c('rice_D_rep1','rice_D_rep2','rice_D_rep3','rice_D_rep4','rice_N_rep1','rice_N_rep2','rice_N_rep3','rice_N_rep4'))
