########## <Argument Pasrsing> ##########
args = commandArgs(trailingOnly = TRUE)


p_assoc.logsitic_ = strsplit(args[1], ',')[[1]]
plot_label_ = args[2]
out_ = args[3]

pointcol_ = args[4]
topcol_ = args[5]
min.pos_ = as.numeric(args[6])
max.pos_ = as.numeric(args[7])
topsignal_ = strsplit(args[8], ',')[[1]] # label of top signal marker
max.yaxis_ = as.numeric(strsplit(args[9], ',')[[1]])

p_GeneBuild_ = args[10]
p_src_ = args[11]

### Checking manually.

# for (i in 1:10) {
#   print(args[i])
# }

# print(p_assoc.logsitic_)
# print(length(p_assoc.logsitic_))
# print(topsignal_)
# print(length(topsignal_))
# print(max.yaxis_)
# print(length(max.yaxis_))


########## <Loading Necessary Scripts(Modules)> ##########
source(file.path(p_src_, "assocplot_HATK.R"))


#### Main

pdf(paste0(out_, ".pdf"), width=10, height=7, pointsize=10) # implementation test (2018.4.25)

NumberofLogistic = length(p_assoc.logsitic_)
layout(matrix(1:(NumberofLogistic+1), (NumberofLogistic+1),1, byrow=TRUE), heights=c(rep(5, NumberofLogistic),3))

# (2018. 8. 10.) No more `r.data` for colorRamp().

for (i in 1:NumberofLogistic) {
  
  print(p_assoc.logsitic_[i])
  print(topsignal_[i])
  print(max.yaxis_[i])

  locus <- TRIM_LOGISTIC.ASSOC(p_assoc.logsitic_[i])
  # manhattan plot
  make.fancy.locus.plot.bare("6", (if(i == 1) plot_label_ else ""), locus, min.pos_, max.pos_, max.yaxis_[i], topsignal_[i], NULL, seq(0,max.yaxis_[i],by=5), pointcolor=pointcol_, topcolor=topcol_)
  
}
# genomic position belt.
make.fancy.locus.plot.bottom("6", min.pos_, max.pos_, pathToTheGeneBuild=p_GeneBuild_)

dev.off()

print("Manhattan Plotting done.")