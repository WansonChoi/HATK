########## <Loading Necessary Scripts(Modules)> ##########
source("assocplot_HATK.R")


########## <Argument Pasrsing> ##########
args = commandArgs(trailingOnly = TRUE)

# assocplot.R의 첫번째 함수의 파라미터 목록
#     : chr, title, locus, min.pos, max.pos, yrange, hitsnp, r.data, yax, pointcolor = "#000000", arrowsnp="", hitsnpB="", hitsnpC=""

# assocplot.R의 두 번째 함수의 파라미터 목록
#     : chr, min.pos, max.pos, pathToTheGeneBuild

p_assoc.logsitic_ = args[1]
out_ = args[2]

pointcol_ = args[3]
min.pos_ = as.numeric(args[4])
max.pos_ = as.numeric(args[5])
topsignal_ = args[6] # label of top signal marker
max.yaxis_ = as.numeric(args[7])

p_GeneBuild_ = args[8]

# for (i in 1:9) {
#   print(args[i])
# }


#### Main

pdf(paste0(out_, ".pdf"), width=10, height=7, pointsize=10) # implementation test (2018.4.25)

# (2018. 8. 9.)
layout(matrix(1:2,2,1,byrow=TRUE),heights=c(5,3)) # No more layout.

r.data <- read.table(p_r2totophit_)

# 논문 내용상 아무런 condition없이 logistic regression을 수행했을때는 DRB1_37
locus <- TRIM_LOGISTIC.ASSOC(p_assoc.logsitic_)
make.fancy.locus.plot.bare("6", "", locus, min.pos_, max.pos_, max.yaxis_, topsignal_, r.data, seq(0,max.yaxis_,by=5), pointcolor=pointcol_)

# # 그 다음, DRB1_37을 condition으로 하고 association test를 수행하면 DQB1_37이 떠오음
# locus <- read.table(file.path(Mfolder,"association_result/CD_ICHIP1_ICHIP2_GWAS/All_CD_DRB1_37.assoc.logistic.input"), header=T)
# make.fancy.locus.plot.bare("6","", locus, 29.60E6, max.pos, 20, hitsnp_57, r.data_57, seq(0,25,by=5), pointcolor="#FF0000")
# 
# # 얘네 둘다 condition으로 잡으면 DRB1_0403이 떠오름
# locus <- read.table(file.path(Mfolder,"association_result/CD_ICHIP1_ICHIP2_GWAS/All_CD_DRB1_37_DQB1_57.assoc.logistic.input"), header=T)
# make.fancy.locus.plot.bare("6","", locus, 29.60E6, max.pos, 10, hitsnp_0403, r.data_0403, seq(0,13,by=5), pointcolor="#FF00FF")
# 
# # 논문 내용상으로는 얘네 셋 까지 condition잡으면 더 이상 significant한 signal이 나타나지 않는다고 했던것 같음.

make.fancy.locus.plot.bottom("6", min.pos_, max.pos_, pathToTheGeneBuild=p_GeneBuild_)

dev.off()