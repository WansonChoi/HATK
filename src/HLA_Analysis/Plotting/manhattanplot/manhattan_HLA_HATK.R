########## <Loading Necessary Scripts(Modules)> ##########
source("assocplot_HATK.R")


########## <Argument Pasrsing> ##########
args = commandArgs(trailingOnly = TRUE)


p_assoc.logsitic_ = args[1]
out_ = args[2]

pointcol_ = args[3]
topcol_ = args[4]
min.pos_ = as.numeric(args[5])
max.pos_ = as.numeric(args[6])
topsignal_ = args[7] # label of top signal marker
max.yaxis_ = as.numeric(args[8])

p_GeneBuild_ = args[9]

# for (i in 1:9) {
#   print(args[i])
# }


#### Main

pdf(paste0(out_, ".pdf"), width=10, height=7, pointsize=10) # implementation test (2018.4.25)

# (2018. 8. 9.)
layout(matrix(1:2,2,1,byrow=TRUE),heights=c(5,3)) # No more layout.

# (2018. 8. 10.) No more `r.data` for colorRamp().

locus <- TRIM_LOGISTIC.ASSOC(p_assoc.logsitic_)
make.fancy.locus.plot.bare("6", "", locus, min.pos_, max.pos_, max.yaxis_, topsignal_, NULL, seq(0,max.yaxis_,by=5), pointcolor=pointcol_, topcolor=topcol_)
make.fancy.locus.plot.bottom("6", min.pos_, max.pos_, pathToTheGeneBuild=p_GeneBuild_)

dev.off()

print("Manhattan Plotting done.")