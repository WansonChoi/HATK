
## v.8b: trying coloring HLA alleles.

require(gplots)
require(RColorBrewer)


########## <Argument Processing> ##########

args = commandArgs(trailingOnly = T)

args1.disease.map_ = args[1]
args2.disease.assoc_ = args[2]
args3.disease.alleleP_ = args[3]
args4.HLA_name_ = args[4]
args5.plot.outf_ = args[5]

# # static argument preparation for testing(set1 : "UC(Ulcerative Colitis)")
# args1.disease.map_ = "/Users/wansun/Git_Projects/HATK/hatk/src/plot/heatmap/7_UC_map.txt"
# args2.disease.assoc_ = "/Users/wansun/Git_Projects/HATK/hatk/src/plot/heatmap/7_UC_assoc.txt"
# args3.disease.alleleP_ = "/Users/wansun/Git_Projects/HATK/hatk/src/plot/heatmap/7_UC_alleleP.txt"
# args4.HLA_name_ = "DRB1"
# args5.plot.outf_ = "8_testws_UC"


# # static argument preparation for testing(set2 : "CD(Crohn's Disease)")
# args1.disease.map_ = "/Users/wansun/Git_Projects/HATK/hatk/src/plot/heatmap/7_testws_CD_map.txt"
# args2.disease.assoc_ = "/Users/wansun/Git_Projects/HATK/hatk/src/plot/heatmap/7_testws_CD_assoc.txt"
# args3.disease.alleleP_ = "/Users/wansun/Git_Projects/HATK/hatk/src/plot/heatmap/7_testws_CD_alleleP.txt"
# args4.HLA_name_ = "DRB1"
# args5.plot.outf_ = "8_testws_CD"


# # static argument preparation for testing(set3 : "UC(Ulcerative Colitis)")
# args1.disease.map_ = "./src/plot/heatmap/test_HEATMAP_NAME.map.txt"
# args2.disease.assoc_ = "./src/plot/heatmap/test_HEATMAP_NAME.assoc.txt"
# args3.disease.alleleP_ = "./src/plot/heatmap/test_HEATMAP_NAME.alleleP.txt"
# args4.HLA_name_ = "DRB1"
# args5.plot.outf_ = "TEST_heatmapheatmapheatmap"

# # static argument preparation for testing((2018. 8. 24.) Cancer Research)
# args1.disease.map_ = "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/HEATMAP_Cancer.map.txt"
# args2.disease.assoc_ = "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/HEATMAP_Cancer.assoc.txt"
# args3.disease.alleleP_ = "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/HEATMAP_Cancer.alleleP.txt"
# args4.HLA_name_ = "A"
# args5.plot.outf_ = "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/HEATMAP_Cancer"

# # (2018. 10. 29.) HATK Integration
# args1.disease.map_ = "/Users/wansun/Git_Projects/HATK/tests/_0_wholeProcess/HAPMAP_CEU/outputs/removelater.A.map.txt"
# args2.disease.assoc_ = "/Users/wansun/Git_Projects/HATK/tests/_0_wholeProcess/HAPMAP_CEU/outputs/removelater.A.assoc.txt"
# args3.disease.alleleP_ = "/Users/wansun/Git_Projects/HATK/tests/_0_wholeProcess/HAPMAP_CEU/outputs/removelater.A.alleleP.txt"
# args4.HLA_name_ = "A"
# args5.plot.outf_ = "/Users/wansun/Git_Projects/HATK/tests/_0_wholeProcess/HAPMAP_CEU/outputs/removelater.A.lklkl"



# ### Arguments checking
# print("Arguments checking")
# print(args1.disease.map_)
# print(args2.disease.assoc_)
# print(args3.disease.alleleP_)




########## <Loading Data> ##########

# maptable=as.matrix(read.table(paste0("7_",disease,"_map.txt"))) # argument[1]
maptable=as.matrix(read.table(args1.disease.map_, check.names = F)) # argument[1]
print(head(maptable))

# (old)
#               AA_DRB1_-25_ AA_DRB1_-24_ AA_DRB1_-17_ AA_DRB1_-16_ AA_DRB1_-1_ AA_DRB1_4_ AA_DRB1_9_ AA_DRB1_10_ AA_DRB1_11_ AA_DRB1_12_
# HLA_DRB1_0101 "K"          "L"          "T"          "A"          "A"         "R"        "W"        "Q"         "L"         "K"        
# HLA_DRB1_0301 "R"          "L"          "A"          "V"          "A"         "R"        "E"        "Y"         "S"         "T"        
# HLA_DRB1_0403 "K"          "F"          "A"          "A"          "A"         "R"        "E"        "Q"         "V"         "K"        
# HLA_DRB1_0404 "K"          "F"          "A"          "A"          "A"         "R"        "E"        "Q"         "V"         "K"        
# HLA_DRB1_0405 "K"          "F"          "A"          "A"          "A"         "R"        "E"        "Q"         "V"         "K"        
# HLA_DRB1_0406 "K"          "F"          "A"          "A"          "A"         "R"        "E"        "Q"         "V"         "K"        

# (new)


# P=as.matrix(read.table(paste0("7_",disease,"_assoc.txt"))) # argument[2]
P=as.matrix(read.table(args2.disease.assoc_, check.names = F)) # argument[2]
print(head(P))

#               AA_DRB1_-25_ AA_DRB1_-24_ AA_DRB1_-17_ AA_DRB1_-16_ AA_DRB1_-1_ AA_DRB1_4_ AA_DRB1_9_ AA_DRB1_10_ AA_DRB1_11_ AA_DRB1_12_
# HLA_DRB1_0101    0.5040396     4.641114     16.31211    0.4623068   -26.00802   3.916856  17.580375   0.5114493  -2.0530567  -0.1444203
# HLA_DRB1_0301   -0.4985299     4.641114    -16.27262   -0.5058454   -26.00802   3.916856  -5.199352   0.1467278   0.1302404   0.1444203
# HLA_DRB1_0403    0.5040396    -4.959793    -16.27262    0.4623068   -26.00802   3.916856  -5.199352   0.5114493  -8.7798919  -0.1444203
# HLA_DRB1_0404    0.5040396    -4.959793    -16.27262    0.4623068   -26.00802   3.916856  -5.199352   0.5114493  -8.7798919  -0.1444203
# HLA_DRB1_0405    0.5040396    -4.959793    -16.27262    0.4623068   -26.00802   3.916856  -5.199352   0.5114493  -8.7798919  -0.1444203
# HLA_DRB1_0406    0.5040396    -4.959793    -16.27262    0.4623068   -26.00802   3.916856  -5.199352   0.5114493  -8.7798919  -0.1444203


# alleleP=as.matrix(read.table(paste0("7_",disease,"_alleleP.txt"))) # argument[3]
alleleP=as.matrix(read.table(args3.disease.alleleP_, check.names = F)) # argument[3]
# (2018.5.21) "AA_DRB1_-25_" ????????? ???????????? label??? ???????????? "AA_DRB1_.25_" ????????? ???????????? ????????? ?????????. check.names = F ?????? ?????????.
print(head(alleleP))

#                        x
# HLA_DRB1_0101 -2.0555667
# HLA_DRB1_0301 -1.8645493
# HLA_DRB1_0403 -0.7473897




########## <Label Pre-processing> ##########

### Preprocessing Labels for `P` DataFrame.

# # (old)
# refine.hla.name <- function(x) { paste0(sprintf("HLA-%s*", args4.HLA_name_), substr(x,10,11), ":", substr(x,12,13)) } # argument[4]
# rownames(P)=sapply(rownames(P), refine.hla.name)
# #                AA_DRB1_-25_ AA_DRB1_-24_ AA_DRB1_-17_ AA_DRB1_-16_ AA_DRB1_-1_ AA_DRB1_4_ AA_DRB1_9_ AA_DRB1_10_ AA_DRB1_11_ AA_DRB1_12_
# # HLA-DRB1*01:01    0.5040396     4.641114     16.31211    0.4623068   -26.00802   3.916856  17.580375   0.5114493  -2.0530567  -0.1444203
# # HLA-DRB1*03:01   -0.4985299     4.641114    -16.27262   -0.5058454   -26.00802   3.916856  -5.199352   0.1467278   0.1302404   0.1444203
# # HLA-DRB1*04:03    0.5040396    -4.959793    -16.27262    0.4623068   -26.00802   3.916856  -5.199352   0.5114493  -8.7798919  -0.1444203
# refine.aa.name <- function(x) { gsub("\\.", "-", paste0(paste0(args4.HLA_name_, "#"), substr(x,9,nchar(x)-1))) }
# colnames(P)=sapply(colnames(P), refine.aa.name)
# #                  DRB1#-25  DRB1#-24  DRB1#-17   DRB1#-16   DRB1#-1   DRB1#4    DRB1#9   DRB1#10    DRB1#11    DRB1#12     DRB1#13   DRB1#14
# # HLA-DRB1*01:01  0.5040396  4.641114  16.31211  0.4623068 -26.00802 3.916856 17.580375 0.5114493 -2.0530567 -0.1444203 -17.9601894 -1.192668
# # HLA-DRB1*03:01 -0.4985299  4.641114 -16.27262 -0.5058454 -26.00802 3.916856 -5.199352 0.1467278  0.1302404  0.1444203  -0.9062282 -1.192668



# (new (by. wanson))

# White space procesing(To adjust the potision of "HLA-*-..." labels)
Labels_nchar = sapply(rownames(P), nchar)
N_WhiteSpace = sapply(Labels_nchar, function(x){max(Labels_nchar)-x})
Labels_prefix = paste0("HLA-", rownames(P), sapply(N_WhiteSpace, function(x){strrep('  ', x)}))

rownames(P) = Labels_prefix

# refine.hla.name <- function(x) { paste0("HLA-", x) } # argument[4]
# rownames(P)=sapply(rownames(P), refine.hla.name)

refine.aa.name <- function(x) { paste0(args4.HLA_name_, "#", x) }
colnames(P)=sapply(colnames(P), refine.aa.name)

cat("min=", min(P), " max=", max(P), "\n")
org.ncol=ncol(P)

# MERGE P AND ALLELE.P
for (i in 1:11) { 
P=cbind(P, alleleP)
maptable=cbind(maptable, rep("", nrow(maptable)))
colnames(P)[ncol(P)]=""
}

#if (disease == "UC")
#    P[P<0]=0
#    brew=rev(brewer.pal(11,"RdBu"))[6:11]
#}
#if (disease == "CD") {
#    P[P>0]=0
#    brew=rev(brewer.pal(11,"RdBu"))[1:6]
#}
#brew=rev(brewer.pal(11,"RdBu"))





########## <Main Plotting> ##########

brew=rev(brewer.pal(11,"Spectral"))

pdf(paste0(args5.plot.outf_, ".pdf"), width=7, height=5.2, pointsize=8) # argument[5]
par(mar=c(5,4,4,3))

mycol=colorRampPalette(brew)(100)
lmat=matrix(c(0,3,2,1,0,4), 3, 2, byrow=T)
lwid=c(1, 12)
lhei=c(1, 4, 0.6)

heatmap.2(P, Rowv=F, Colv=F, dendrogram="none", col=mycol,
          # density.info=c("histogram","density","none"), # ?????? ?????? ?????? ?????? ???????????? ????????????????????? ?????? ??????????????? ?????? ?????? ???????????? ???????????????????????? ?????????
          density.info="none",
          trace="none", 
          margin=c(6,6), # ??? margin??? plot area??? heatmap??? ???????????? margin????????????. ????????? AA??? Label??? ????????? ?????????.(?????? ???????????? ????????? ??????.)
          lmat=lmat,lwid=lwid,lhei=lhei,
          cellnote=maptable, notecol="#909090", cexRow=1.2, cexCol=1.2, adjCol=c(NA,0.4), adjRow=c(1.3,NA),
          
          colsep=c(org.ncol), sepcolor="white", sepwidth=.6,
          
          key.title="this is keys", key.par=list(mar=c(4, 4, 0.5 ,10)), ## Bottom, Left, Top, Right
          key.xlab=bquote(-log[10]~italic("P")), # "-log10P"?????? ????????? ??????
          key.xtickfun=function() {
            breaks=pretty(parent.frame()$breaks)
            list(at=parent.frame()$scale01(breaks), labels=abs(breaks))
          }
)
# (2018. 6. 20.)
#arrow(0,0,1,1,xpd=T)
dev.off()


# # (old verison) >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 
# brew=rev(brewer.pal(11,"Spectral"))
# 
# pdf(paste0(args5.plot.outf_, ".pdf"), width=7, height=5.2, pointsize=8) # argument[5]
# par(mar=c(5,3,4,3))
# mycol=colorRampPalette(brew)(100)
# lmat=matrix(c(0,3,2,1,0,4), 3, 2, byrow=T)
# lwid=c(1, 4)
# lhei=c(1, 4, 0.6)
# heatmap.2(P, Rowv=F, Colv=F, dendrogram="none", col=mycol, density.info="none",
#           trace="none", margin=c(6,6),
#           lmat=lmat,lwid=lwid,lhei=lhei,
#           cellnote=maptable, notecol="#909090", cexRow=1.2, cexCol=1.2, adjCol=c(NA,0.4), adjRow=c(1.3,NA),
#           colsep=c(org.ncol), sepcolor="white", sepwidth=.6,
#           key.title="this is keys", key.par=list(mar=c(4, 4, 0.5 ,10)), ## Bottom, Left, Top, Right
#           key.xlab=bquote(-log[10]~italic("P")),
#           key.xtickfun=function() {
#               breaks=pretty(parent.frame()$breaks)
#               list(at=parent.frame()$scale01(breaks), labels=abs(breaks))
#           }
#           )
# #arrow(0,0,1,1,xpd=T)
# dev.off()
# 
# # warnings() # (2018.5.21) ?????? ???????????? ??????????????? ???.
# 
# 
# # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<