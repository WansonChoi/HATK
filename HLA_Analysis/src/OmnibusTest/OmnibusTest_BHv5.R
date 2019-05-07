###### Omnibus test code
###### Buhm Han, Xinli Hu
###### v1.00 : 2016/12/21
###### v2.00 : 2017/05/11 <- I fixed errors in conditioning.
###### v3.00 : 2017/08/01 <- *.covar file, i.e. collection (group) as covariates.
###### v4.00 : 2018/08/07 <- Now argument "data.prefix" is deprecated and divided into each new argument for those values.

# (Before)
# (1) file prefix
# (2) output prefix
# (3) phenotype name
# (4) rare threshold
# (5) conditional variables

# args <- commandArgs(TRUE)
# data.prefix=args[1] ## prefix of .fam, .aa, .phe, and .covar files
# outfile=args[2]
# phen.name=args[3] ## phenotype name defined in .phe file
# rare.thres=as.numeric(args[4]) ## frequency cut-off to ignore alleles in test
# condvar=args[5] ## comma separated, without spaces. 


# (After)
# (1) output prefix
# (2) .fam
# (3) .aa
# (4) .phe
# (5) pheno-name (2018. 11. 2.)
# (6) .covar
# (7) covar-name (2018. 11. 2.)
# (8) rare threshold
# (9) conditional variables.


args <- commandArgs(TRUE)
outfile=args[1]
data.fam = args[2]
data.aa = args[3]
data.phe = args[4]
data.phe.name=args[5] ## phenotype name defined in .phe file (must be single item.)
data.covar = args[6]
data.covar.name = args[7] # 
# rare.thres=as.numeric(args[8]) ## frequency cut-off to ignore alleles in test
condvar=args[8] ## comma separated, without spaces. 
# If no conditioning, provide "NA"
# Put "_All" for all search results: "AA_DRB1_All,AA_DQA1_All,AA_DQB1_All"  


########## <(1) Preprocessing given arguments.> ##########

# condvar trimming
condvars=unlist(strsplit(condvar,split=','))
if(condvars[1]=="NA") { condvars=c() }

# input data-file name trimming
famfile=paste0(data.fam)
phasedfile=paste0(data.aa) ## Amino acid called file
covfile=paste0(data.covar)
phenfile=paste0(data.phe)

# print(famfile)
# print(phasedfile)
# print(covfile)
# print(phenfile)



########## <(2) Loading Data> ##########

## 1. READ FAM, PHASED, COVAR, PHE
# (1) fam
fam <- as.matrix(read.table(famfile))
fam=rep(fam, each=2) ## make x2 for haploids

# (2) covar
covarexp=""
if (file.exists(covfile)) {
  
    covar <- as.matrix(read.table(covfile, header=T))[,-(1:2)]

    
    ### covar_name processing (2018. 11. 2.)
    if (data.covar.name !="NA"){
      
      covar_header = 1:(dim(covar)[2])
      names(covar_header) = colnames(covar)
      
      covar_target = strsplit(data.covar.name, ",(\\s+)?")[[1]]
      
      for (item in colnames(covar)) {
        covar_target = sub(item, covar_header[item], covar_target)
      }
      
      tryCatch({
        covar_target<-unlist(sapply(sub('-', ':', covar_target), function(x){eval(parse(text = x))}), use.names = F)
        }, error = function(e){stop(simpleError("[OmnibusTest::ERROR] The argument 'covar_name' is wrong. Some of covariate names can't be found in given *.covar file. Please check it again."))})
      
      # Subsetting columns
      covar = covar[,covar_target]
    }
    
    ### Duplicating for .bgl.phased format
    covar=apply(covar, 2, rep, each=2) ## make x2 for haploids
    
    ### Covar Expression
    # covarexp=" + covar.f[,1] + covar.f[,2]"
    covarexp = ""
    
    for (i in 1:dim(covar)[2]) {
      covarexp = paste0(covarexp, " + covar.f[,", i, "]")
    }

    print(dim(covar))

} else {
    # 없는 경우는 그냥 NA집어 넣어 놓음.
    covar=NA
}

# (3) pheno
pheno <- as.numeric(as.matrix(read.table(phenfile, header=T))[, data.phe.name])-1
#mode(pheno)="numeric"
#pheno=pheno-1 # 2,1,0 to 1,0,-1
pheno=rep(pheno, each=2) ## make x2 for haploids
is.live.phen=(pheno > -1) ## who is alive for analysis

# (4) phased
phased.in <- as.matrix(read.table(phasedfile)) # 얘가 phased된 데이터파일이라 로드하는데 시간이 꽤 김
variants=phased.in[-(1:5),2]
phased=t(phased.in[-(1:5),-(1:2)]) ## transpose: rows are individuals, cols are markers
colnames(phased)=variants



########## <(Optional) Condvar> ##########

## 2. DEFINE HAPLOTYPES TO BE CONDITIONED ON
# (2018.4.18) 일단 NA인 경우는 여기 안간다 치자.
is.live.cond=rep(TRUE, length(pheno)) # who is alive for analysis / (2017.04.18) pheno의 길이만큼 전부 true로 채워져 있음.
condexp=""
hap.cond=NULL
if (length(condvars)>0) {
	condvariants=c()
	for (thiscon in condvars) {
		if (grepl("AA",thiscon) || grepl("HLA",thiscon)) {
			if (grepl("_All",thiscon)) {
				thiscon <- strsplit(thiscon,"_All")[[1]][1] ## remove "_All" before search
            }
		    thiscon <- variants[grepl(thiscon,variants)] ## search variables
		}
		condvariants <- c(condvariants, thiscon)
	}
    #print(condvariants)
    condmatrix=phased[ ,condvariants, drop=FALSE]
    condexp=" + hap.cond.f"
    hap.cond=apply(condmatrix, 1, paste0, collapse='')
    is.live.cond=(apply(is.na(condmatrix), 1, sum)==0)
}


########## <(3) Regression> ##########


## 3. BEGIN REGRESSION
aavariants <- variants[which(grepl("AA_",variants))] ## test only amino acids
testvariants <- aavariants
#if (!is.null(aacondvariants)) {
#	testvariants <- aavariants[which(!aavariants %in% aacondvariants)]
#}
results <- matrix(ncol=6, nrow=length(testvariants))

for (i in 1:length(testvariants)) {
	  if (i%%10 == 0) {
		  cat(i,'\n')
  	}
	  #variant <- as.character(testvariants[i])
    variant=testvariants[i]
  	vcol <- which(variants==variant)
  	newaa <-phased[,vcol]
    is.live.aa=!is.na(newaa)

    is.live=(is.live.phen & is.live.cond & is.live.aa) ## who is finally alive
    n.is.live=sum(is.live)
    #print(length(is.live.phen))
    #print(length(is.live.cond))
    #print(length(is.live.aa))

    ## KEEP ALIVE PEOPLE
    pheno.f=pheno[is.live]
    if (!is.null(hap.cond)) {
        # NULL이 아닌 경우 <=> conditional argument가 주어져서 condvars처리하는 부분 거치고 온 경우.
        hap.cond.f=hap.cond[is.live]
    } else {
        # NULL인 경우 <=> conditional이 안 주어진 경우
        hap.cond.f=NULL
    }
    newaa.f=newaa[is.live]

    if(!is.na(covar)){
      covar.f=covar[is.live,]
    }

    ## DEFINE NEW HAPLOTYPES
  	residues <- unique(newaa.f)
  	#hap.new.f = apply(cbind(hap.cond.f,newaa.f),1,paste0,collapse='',sep='')
  	hap.new.f = apply(cbind(newaa.f),1,paste0,collapse='',sep='')

  	#ALTERNATIVE MODEL
  	if (length(unique(residues)) == 1) {
    		results[i,] <- c(variant,"NaN","NaN","NaN","NaN",paste0(residues,collapse=','))
  	}
  	if (length(unique(residues)) > 1) {

        nullexp <- paste0("glm(pheno.f ~ 1", covarexp, condexp, ", family=binomial(logit))")
        glm.null <- eval(parse(text=nullexp))
        nulldeviance <- summary(glm.null)$deviance
        nulldf <- summary(glm.null)$df[1]

        altexp <- paste0("glm(pheno.f ~ 1", covarexp, condexp, "+ hap.new.f, family=binomial(logit),maxit=100)")
    	glm.alt <- eval(parse(text=altexp))
    	altdeviance <- summary(glm.alt)$deviance
    	altdf <- summary(glm.alt)$df[1]

        #print(summary(glm.alt))
        #print(sum(pheno.f==1 & hap.new.f=="Y"))
        #print(sum(pheno.f==1 & hap.new.f=="V"))
        #print(sum(pheno.f==0 & hap.new.f=="Y"))
        #print(sum(pheno.f==0 & hap.new.f=="V"))

		#STATISTICS
    	deviancediff <- nulldeviance - altdeviance
    	dfdiff <- altdf - nulldf
    	log10pvalue <- pchisq(deviancediff, df=dfdiff, lower.tail=FALSE, log.p=TRUE)/log(10)


        ### (2018.4.19) 개인미팅 정리
        # 결국 deviance라는 통계 변량으로 Likelihood Raio test를 하고싶은거. deviance를 차를 구함으로써 LR test를 수행하는 방법이 있음. 저 qchisq라는 것 또한 결국 LR test를 하는 문맥.
        # log(10)으로 나눠주는건 가끔씩 너무 significant한 p-value가 나오면 너무 작은 값이라 처리를 못함. 이거를 log(10)으로 나눠서 캐치해내는거.

    	results[i,] <- c(variant, deviancediff, dfdiff, n.is.live, log10pvalue, paste0(residues,collapse=','))
  	}
}

outfilename = paste(outfile, data.phe.name, paste0(condvar,collapse='+'), "omnibus", sep='.', collapse='.')
write.table(results, outfilename, quote=F,sep='\t',row.names=F,col.names=c("Variant","deltaDeviance","deltaDF","N","P","Residues"))



