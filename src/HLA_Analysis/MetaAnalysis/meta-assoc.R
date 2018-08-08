#!/usr/bin/env Rscript

ptm=proc.time() # for measuring running time

args <- commandArgs(TRUE)

kordata = args[1]
jpndata = args[2]
markerdata = args[3]
outf = paste(args[4], 'hmeta', sep='.') # modified by Wanson Choi (2018. 4. 29)

# print(kordata)
# print(jpndata)
# print(markerdata)
# print(outf)

### original codes
# kordata="../01-association/All_CD.assoc.logistic"
# jpndata="../../Data/JapanData/RIKEN_CD_HLA-imputation-binary.csv"
# markerdata="../../Data/Reference/K2panel/K2.markers"
# outf="meta.out.txt"

# ### for Testing(2018.4.17)
# kordata="./analysis/01-association/All_CD.assoc.logistic"
# jpndata="./Data/JapanData/RIKEN_CD_HLA-imputation-binary.csv"
# markerdata="./Data/Reference/K2panel/K2.markers"
# outf="meta.out.txt"


ko=as.matrix(read.table(kordata, header=T))
print(head(ko))
jp=as.matrix(read.csv(jpndata))
mk=as.matrix(read.table(markerdata))

out=NULL

complement=c("A", "T", "G", "C")
names(complement)=c("T", "A", "C", "G")

for (i in 1:nrow(ko)) {

    snp=ko[i,'SNP']
    
    ## Find alleles 
    k=which(mk[,1]==snp)
    Ak=mk[k,3]
    Bk=mk[k,4]
    j=which(jp[,1]==snp)
    
    if (length(j)==0) {
        next ## pass, since Japan data doesn't have this variable.....
    }
    
    Aj=jp[j,5]
    Bj=jp[j,4]
    
    ## Get summary stats.
    betak=log(as.numeric(ko[i,'OR']))
    sek=as.numeric(ko[i,'SE'])
    betaj=as.numeric(jp[j,'B_crude'])
    zj=qnorm(as.numeric(jp[j,'P_crude'])/2)
    sej=abs(betaj/zj)
    
    ## Flip, depending on alleles. 
    if (length(union(c(Ak,Bk), c(Aj,Bj)))==4) {
        Aj=complement[Aj]
        Bj=complement[Bj]
    }
    if (Ak==Aj && Bk==Bj) {
        ## Good
    } else if (Ak==Bj && Bk==Aj) {
        ## Perfect flip!!
        betaj=-betaj
    } else {
        stop(sprintf("Something's wrong: (%s,%s) vs (%s,%s) at %s\n",Ak,Bk,Aj,Bj,snp))
    }
    
    ## Now let's do meta-analysis.
    wk=1/sek**2
    wj=1/sej**2
    ivw=(betak*wk+betaj*wj)/(wk+wj)
    ivw.se=sqrt(1/(wk+wj))
    meta.z=ivw/ivw.se
    meta.p=2*pnorm(-abs(meta.z))
    ## Print output
    out=rbind(out, format(list(SNP=snp, BP=ko[i,'BP'], A1=Ak, A2=Bk, OR=exp(ivw), BETA=ivw, SE=ivw.se, MetaP=meta.p), digits=3))
}
write.table(out, outf, row.names=F, quote=F, sep="\t")

print(proc.time()-ptm)

