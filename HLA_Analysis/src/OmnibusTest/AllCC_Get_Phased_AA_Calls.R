########## < 1. Arguments > ##########

args = commandArgs(trailingOnly = T)

bgl.phased_ = args[1]
fam_ = args[2]
output_ = args[3]

# argument example

# # (1) original given by professor Han.
# bgl.phased_ = "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/UC-CD-HLA/data/Merged/merged.bgl.phased"
# fam_ = "/Users/wansun/Dropbox/_Sync_MyLaptop/Data/UC-CD-HLA/data/Merged/merged.fam"
# output_ = "/Users/wansun/Projects/20181004_Omnibus_test/merged.aa"
# 
# # (2) The result from HATK (2018. 10. 08.)
# bgl.phased_ = "/Users/wansun/Git_Projects/HATK/tests/_2_HLA2MARKER/Test_20181003/RemovethisPanel.MERGED.bgl.phased"
# fam_ = "/Users/wansun/Git_Projects/HATK/tests/_2_HLA2MARKER/Test_20181003/RemovethisPanel.MERGED.fam"
# output_ = "/Users/wansun/Git_Projects/HATK/tests/_4_HLA_Analysis/Omnibus/20181008_TEST"



# datadir = "/Users/buhmhan/Dropbox/Projects-open/UC-CD-HLA/Data/Merged"


########## < 2. Necessary functions > ##########

#DEFINE FXNs
getAAName <- function(aavar) {
	cut <- strsplit(aavar,"_")[[1]][1:3]
	AA <- paste0(cut,collapse="_")
	return(AA)
}

getAAPositions <- function(AA,allvar) {
	AAname <- paste0(AA,"_",sep='')
	thisaa <- variants[grepl(AAname,allvar)]
	positions <- c()
	if (length(thisaa)==1) {
		positions = thisaa
	}
	if (length(thisaa) > 1) {
		allelements <- lapply(thisaa,strsplit,split="_")
		aaidentity = rapply(allelements,function(x) tail(x,1))
	#	xdel <- which(aaidentity!="x")
	#	thisaa <- thisaa[xdel]
	#	aaidentity <- aaidentity[xdel]
		aalength <-nchar(aaidentity)
	#	keep <- thisaa[xdel]
		keep <- thisaa #
		positions <- thisaa[which(aalength==1)]
	}
	return(positions)
}

getResidueName <- function(vars) {
	residues <- c()
	if (length(vars) ==1 ) {
		residues = "check"
	}
	if (length(vars) > 1) {
		allelements <- lapply(vars,strsplit,split="_")
		residues = rapply(allelements,function(x) tail(x,1))
	}
	return(residues)
}

isP <- function(calls) {
	#check only one P
	index = "NA"
	hasP <- "P" %in% names(table(calls))
	if (hasP) {
		pCount <- table(calls)["P"]
		if (pCount == 1) {
			index = which(calls == "P")	
		}	
	}else{print("I found")}
	return(index)
}

isp <- function(calls) {
  # Modification of the function `isP` so that it can work with "HLA2MARKER".
	# check only one 'p'(new version) not 'P'(original version.)
	index = "NA"
	hasp <- "p" %in% names(table(calls))
	if (hasp) {
		pCount <- table(calls)["p"]
		if (pCount == 1) {
			index = which(calls == "p")	
		}	
	}else{print("I found")}
	return(index)
}

# setwd(datadir)


########## < 3. Main > ##########

#READ VARIANT AND FAM DATA (Loading Data)
phased <- read.table(bgl.phased_) # (*) args1
fam <- read.table(fam_) # (*) args2

#REMOVE INDIVIDUALS MISSING DATA
keepphased <- phased

#EXTRACT AA VARIANTS
variants <- as.character(phased[,2])
# allaavars <- variants[which(grepl("AA",variants)==TRUE)]
allaavars <- variants[which(grepl("^AA",variants)==TRUE)]
aanames <- unique(unlist(lapply(allaavars,getAAName)))

aavariants.1 <- unlist(lapply(aanames,getAAPositions,allva=variants))
aavariants <- variants[which(variants %in% aavariants.1)]

aaphasedindex <- which(variants %in% aavariants)
aaphased <- keepphased[aaphasedindex,3:ncol(keepphased)] # <- phased aa calling
#aanames <- aanames[1:5]

#GET AAs
testmatrix <- matrix(nrow = length(aanames),ncol=ncol(aaphased))
use.test <- 1:ncol(aaphased)

for (i in 1:length(aanames)) {
	cat(i,'\n')
	residue <- c()
	aa <- aanames[i]
	aa <- paste0(aa,"_",sep='')
	allvariants <- aavariants[grepl(aa,aavariants)]
	allresidues <- getResidueName(allvariants)

	phaseRow <- which(aavariants %in% allvariants)
	variant <- aa
	allphased <-t(aaphased[phaseRow,])

	if (allresidues[1] == "check") {
		residue = allphased
	}

	if (allresidues[1] != "check") {
		residueInd <- unlist(apply(allphased,1,isp)) # replaced by the function `isp` not `isP` (2018. 10. 08.).
		use <- which(residueInd != "NA")
		nouse <- which(residueInd == "NA")
		use.test <- use[which(use %in% use.test)]
		residue[use] <- allresidues[as.numeric(residueInd[use])]
		residue[nouse] = "NA"
	}
	testmatrix[i,] <- residue
}

AAmatrix <- cbind("M",aanames,testmatrix)
colnames(AAmatrix) <- colnames(phased)

# notAA <- which(grepl("AA",variants)==FALSE)
notAA <- which(grepl("^AA",variants)==FALSE)
notAAmatrix <- phased[notAA,]

newmatrix <- matrix(ncol=ncol(phased),nrow = nrow(AAmatrix)+nrow(notAAmatrix))
newmatrix[1:nrow(notAAmatrix),] <- as.matrix(notAAmatrix)
newmatrix[(nrow(notAAmatrix)+1):nrow(newmatrix),] <- as.matrix(AAmatrix)

write.table(newmatrix,paste(output_, "aa", sep = '.') ,col.names=F,row.names=F,quote=F) # (*) args3
	



