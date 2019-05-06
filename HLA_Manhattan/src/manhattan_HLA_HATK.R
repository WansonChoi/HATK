########## <Argument Pasrsing> ##########
args = commandArgs(trailingOnly = TRUE)


p_assoc.logsitic_ = strsplit(args[1], ',')[[1]]
# plot_label_ = args[]
out_ = args[2]

pointcol_ = args[3]
pointsize_ = as.numeric(args[4])
topcol_ = args[5]
min.pos_ = as.numeric(args[6])
max.pos_ = as.numeric(args[7])
topsignal_ = strsplit(args[8], ',')[[1]] # label of top signal marker
max.yaxis_ = as.numeric(strsplit(args[9], ',')[[1]])
yaxis_unit_ = as.numeric(args[10])

p_GeneBuild_ = args[11]
p_src_ = args[12]

### Checking manually.

# for (i in 1:13) {
#   print(args[i])
# }

# print(p_assoc.logsitic_)
# print(length(p_assoc.logsitic_))
# print(topsignal_)
# print(length(topsignal_))
# print(max.yaxis_)
# print(length(max.yaxis_))

print(args)

########## <Loading Necessary Scripts(Modules)> ##########
source(file.path(p_src_, "assocplot_HATK.R"))


#### Main

pdf(paste0(out_, ".pdf"), width=16, height=9, pointsize=pointsize_) # `pointsize` is now manipulated by the argument 'pointsize_'. (2019. 04. 09.)

NumberofLogistic = length(p_assoc.logsitic_)
layout(matrix(1:(NumberofLogistic+1), (NumberofLogistic+1),1, byrow=TRUE), heights=c(rep(7, NumberofLogistic), 2))
# (2019. 04. 01.) One manhattan graph : One fancy buttom = 7 : 2

# (2018. 8. 10.) No more `r.data` for colorRamp().

for (i in 1:NumberofLogistic) {
  
  print(p_assoc.logsitic_[i])
  print(topsignal_[i])
  print(max.yaxis_[i])

  locus <- TRIM_LOGISTIC.ASSOC(p_assoc.logsitic_[i])
  # manhattan plot
  make.fancy.locus.plot.bare("6", (if(i == 1) "" else ""), locus, min.pos_, max.pos_, max.yaxis_[i], topsignal_[i], NULL, seq(0,max.yaxis_[i],by=yaxis_unit_), topcolor=topcol_, p.color.AA = "#CF3721", p.color.HLA = "#F5BE41", p.color.SNP = "#626262") #
  
}
# genomic position belt.
make.fancy.locus.plot.bottom("6", min.pos_, max.pos_, pathToTheGeneBuild=p_GeneBuild_)

dev.off()

print("Manhattan Plotting done.")