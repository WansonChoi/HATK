## v.8b: trying coloring HLA alleles.


########## < Loading necessary packages > ##########

if("gplots" %in% rownames(installed.packages()) == F){
  install.packages("gplots")
}
library("gplots")

if("RColorBrewer" %in% rownames(installed.packages()) == F){
  install.packages("RColorBrewer")
}
library("RColorBrewer")

if("shape" %in% rownames(installed.packages()) == F){
  install.packages("shape")
}
library("shape")



########## < Argument Processing > ##########

args = commandArgs(trailingOnly = T)

args1.disease.map_ = args[1]
args2.disease.assoc_ = args[2]
args3.disease.alleleP_ = args[3]
args4.HLA_name_ = args[4]
args5.plot.outf_ = args[5]


### Example

# (2019. 04. 08.)
# args1.disease.map_ = "/Users/wansun/Git_Projects/HLA_Heatmap/tests/WTCCC_RA_DRB1/WTCCC_RA_DRB1.map.txt"
# args2.disease.assoc_ = "/Users/wansun/Git_Projects/HLA_Heatmap/tests/WTCCC_RA_DRB1/WTCCC_RA_DRB1.assoc.txt"
# args3.disease.alleleP_ = "/Users/wansun/Git_Projects/HLA_Heatmap/tests/WTCCC_RA_DRB1/WTCCC_RA_DRB1.alleleP.txt"
# args4.HLA_name_ = "DRB1"
# args5.plot.outf_ = "/Users/wansun/Git_Projects/HLA_Heatmap/tests/WTCCC_RA_DRB1/WTCCC_RA_DRB1_heatmapanalysis.pdf"


# ### Arguments checking
# print("Arguments checking")
# print(args1.disease.map_)
# print(args2.disease.assoc_)
# print(args3.disease.alleleP_)




########## < Loading Data > ##########

maptable=as.matrix(read.table(args1.disease.map_, check.names = F)) # argument[1]
# print(head(maptable))

P=as.matrix(read.table(args2.disease.assoc_, check.names = F)) # argument[2]
# print(head(P))

alleleP=as.matrix(read.table(args3.disease.alleleP_, check.names = F)) # argument[3]
# print(head(alleleP))



########## < Label Pre-processing > ##########

refine.hla.name <- function(x) { paste0("HLA-", x) } # argument[4]
rownames(P)=sapply(rownames(P), refine.hla.name)

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





########## < Main Plotting > ##########

brew=rev(brewer.pal(11,"Spectral"))
mycol=colorRampPalette(brew)(100)

pdf(paste0(args5.plot.outf_, ".pdf"), width=7, height=5.2, pointsize=8) # argument[5]
par(mar=c(5,4,4,3))

lmat=matrix(c(0,3,2,1,0,4), 3, 2, byrow=T)
lwid=c(1, 12)
lhei=c(1, 4, 0.6)

heatmap.2(P, Rowv=F, Colv=F, dendrogram="none", col=mycol,
          density.info="none",
          trace="none", 
          margin=c(6,6), 
          lmat=lmat,lwid=lwid,lhei=lhei,
          cellnote=maptable, notecol="#909090", cexRow=1.2, cexCol=1.2, adjCol=c(NA,0.4), adjRow=c(1.3,NA),
          
          colsep=c(org.ncol), sepcolor="white", sepwidth=.6,
          
          key.title="this is keys", key.par=list(mar=c(4, 4, 0.5 ,10)), ## Bottom, Left, Top, Right
          key.xlab=bquote(-log[10]~italic("P")),
          key.xtickfun=function() {
            breaks=pretty(parent.frame()$breaks)
            return(list(at=parent.frame()$scale01(as.numeric(breaks), parent.frame()$min.raw, parent.frame()$max.raw), labels=abs(breaks)))
          }
)
### Adding arrows (2019. 04. 08.)
# toward right
Arrows(0.53,-0.135, 0.84,-0.135, xpd=T, arr.type='triangle',
       lwd=0.5, arr.length=0.1, arr.width = 0.1)
text(x=0.885, y=-0.135, labels = "Risk", xpd=T)

# toward left
Arrows(0.44,-0.135, 0.16,-0.135, xpd=T, arr.type='triangle',
       lwd=0.5, arr.length=0.1, arr.width = 0.1)
text(x=0.105, y=-0.135, labels = "Protective", xpd=T)

dev.off()