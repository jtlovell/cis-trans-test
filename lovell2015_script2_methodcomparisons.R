# Final Differential expression and cis-trans pipeline.
# Part 2: Comparison of methods for cis trans tests to define how our approach compares to other method
# Author: JT Lovell, S Schwartz
# Date: 17-April 2015
# Version: 5.2

# Notes: Here we compare the following methods - 
# 1. binomial tests for cis and trans
# 2. "full wald cis trans test" where we fit a deseq model with allele, generation and their interaction (base levels = F1,Fil)

####################
####################
####################
####################
# 2. Prepare datasets
rm(list=ls())
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/ph2015_eqtl")
pkg <- c("RCurl","plyr","mclust","qtl","DESeq2","GenomicRanges","car","MASS")
invisible(lapply(pkg, function(x) {cat(x,"..."); library(x, character.only=T, verbose=F, warn.conflicts=F,quietly=T)} ))
sessionInfo()

####################
# 2.1: Load DESeq objects
load("ct2015final_deseqinput.RData")

####################
# 2.2: Run Deseq for cis test (8) and trans test (9)
des.ase8 <- DESeq(dds.ase8, test="Wald")
des.ase9 <- DESeq(dds.ase9, test="Wald")
des.ase12 <- DESeq(dds.ase12, test="Wald")
####################
# 2.3 Extract results, plot histograms

####################
# function to do this:
DESeqResults<-function(DESeqobj, annotation, contrasts, prefix){
  out<-annotation
  for(i in 1:length(prefix)){
    res <- data.frame(results(DESeqobj, contrast=contrasts[i], pAdjustMethod = "BH"))
    colnames(res)<-paste(prefix[i],colnames(res),sep=".")
    dat<-res[,grep("padj", colnames(res))]
    res[,gsub("padj","sig",colnames(res)[grep("padj", colnames(res))])]<-ifelse(dat>0.05 | is.na(dat),"","sig")
    out<-cbind(out,res)
    hist(dat, main=prefix[i], breaks=100,xlim=c(0,1))
    cat(prefix[i], "n alpha<0.05 = ",length(dat[dat<0.05]),"\n")
  }
  return(out)
}
####################
res.ase8<-DESeqResults(DESeqobj=des.ase8,
                       annotation=ase.annot, 
                       contrasts=sapply("allelehal", list),
                       prefix=c("allele"))
res.ase9<-DESeqResults(DESeqobj=des.ase9,
                       annotation=ase.annot, 
                       contrasts=sapply("generationF0.allelehal", list),
                       prefix=c("allele.gen"))
res.ase12<-DESeqResults(DESeqobj=des.ase12,
                       annotation=ase.annot, 
                       contrasts=sapply(c("trans.filtrans.fil","trans.haltrans.hal"), list),
                       prefix=c("trans.fil","trans.hal"))

des.cis<-res.ase8[,c("id","allele.pvalue","allele.padj")]
des.cis$cis.des.rank<-rank(des.cis$allele.pvalue)
des.trans1<-res.ase9[,c("id","allele.gen.pvalue","allele.gen.padj")]
des.trans1$trans.des.rank<-rank(des.trans1$allele.gen.pvalue)
des.trans2<-res.ase12[,c("id","trans.fil.pvalue","trans.fil.padj","trans.hal.pvalue","trans.hal.padj")]
des.trans2$trans.fil.des.rank<-rank(des.trans2$trans.fil.pvalue)
des.trans2$trans.hal.des.rank<-rank(des.trans2$trans.hal.pvalue)
des.ps<-merge(des.cis,des.trans1, by="id")
des.ps<-merge(des.ps,des.trans2, by="id")
####################
# 2.4 prepare dataset for traditional binomial cis-trans test
quantnorm<-function(x) {
  n=sum(!is.na(x),na.rm=T)
  x=rank(x)/(n+1)
  x=qnorm(x)
  x[is.infinite(x)]=NA
  x
}
normalizedCounts <- t(counts(des.ase8)) / sizeFactors(des.ase8)
ct.info<-data.frame(ase.info) #info for counts
qv<-apply(normalizedCounts,1,quantnorm) #quantile norm counts
vsd.cistrt<-data.frame(ct.info,normalizedCounts) #combine info and qn counts
vsd.cistrt<-vsd.cistrt[,c("allele","generation",as.character(rownames(des.ase8)))] #drop unneeded columns

####################
# 2.5 Calculate cis vs. trans following Wittkopp 2004
witt2004<-function(g){
  cis<-wilcox.test(
    (j[j$generation=="F1" & j$allele=="fil",g]+1) / (j[j$generation=="F1" & j$allele=="hal",g]+1),
    mu=1)$p.value
  trans<-wilcox.test(
    (j[j$generation=="F1" & j$allele=="fil",g]+1) / (j[j$generation=="F1" & j$allele=="hal",g]+1),
    mu=mean(j[j$generation=="F0" & j$allele=="fil",g]+1) / mean(j[j$generation=="F0" & j$allele=="hal",g]+1))$p.value
  return(c(cis, trans))
}
j<-vsd.cistrt
out<-sapply(colnames(j)[-c(1,2)], witt2004)
wilcox.df<-data.frame(t(data.frame(out)))
wilcox.df$cis.wilcox.rank<-rank(wilcox.df[,1])
wilcox.df$trans.wilcox.rank<-rank(wilcox.df[,2])
wilcox.df$cis.wilcox.padj<-p.adjust(wilcox.df[,1], method="fdr")
wilcox.df$trans.wilcox.padj<-p.adjust(wilcox.df[,2], method="fdr")
wilcox.df$id<-rownames(wilcox.df)
all.ps<-merge(des.ps,wilcox.df,by="id")
####################
# 2.6 Calculate cis vs. trans following Wittkopp 2010

means.des<-vsd.cistrt %>% group_by(allele,generation) %>% summarise_each(funs(mean))
means.des<-data.frame(means.des)
means.long<-data.frame(t(means.des[,-c(1,2)]))
colnames(means.long)<-c("fil_f1","fil_f0","hal_f1","hal_f0")
means.long$lfc.f1<-log2(means.long$fil_f1)-log2(means.long$hal_f1)
means.long$lfc.f0<-log2(means.long$fil_f0)-log2(means.long$hal_f0)
means.long$id<-rownames(means.long)

means.long$cis.p<-apply(means.long[,c("fil_f1","hal_f1")],1, 
                        function(x) binom.test(round(as.numeric(x),0),p=.5)$p.value)
means.long$trans.p<-apply(means.long[,c("fil_f1","hal_f1","fil_f0","hal_f0")],1, 
                          function(x) binom.test(round(as.numeric(x[1:2]),0),p=as.numeric(x[3]/(x[3]+x[4])))$p.value)
means.long$cis.binom.rank<-rank(means.long$cis.p)
means.long$trans.binom.rank<-rank(means.long$trans.p)
means.long$cis.padj<-p.adjust(means.long$cis.p, method="fdr")
means.long$trans.padj<-p.adjust(means.long$trans.p, method="fdr")
all.ps<-merge(all.ps, means.long, by="id")


####################
# 2.7 plot various cis-trans test methods
pdf("ct2015_cistransmethodcomps.pdf", width=8, height=12)
par(mfrow=c(3,2))
with(all.ps,plot(cis.des.rank,cis.binom.rank, pch=".",
                 main=paste("des cis vs. binomial cis \n r2 =", round(summary(lm(cis.des.rank~cis.binom.rank, all.ps))$r.squared,2))))
with(all.ps,plot(cis.des.rank,cis.wilcox.rank, pch=".",
                 main=paste("des cis vs. wilcox cis \n r2 =", round(summary(lm(cis.des.rank~cis.wilcox.rank, all.ps))$r.squared,2))))
with(all.ps,plot(cis.binom.rank,cis.wilcox.rank, pch=".",
                 main=paste("binomial cis vs. wilcox cis \n r2 =", round(summary(lm(cis.binom.rank~cis.wilcox.rank, all.ps))$r.squared,2))))

with(all.ps,plot(trans.des.rank,trans.binom.rank, pch=".",
                 main=paste("des trans vs. binomial trans \n r2 =", round(summary(lm(trans.des.rank~trans.binom.rank, all.ps))$r.squared,2))))
with(all.ps,plot(trans.des.rank,trans.wilcox.rank, pch=".",
                 main=paste("des trans vs. wilcox trans \n r2 =", round(summary(lm(trans.des.rank~trans.wilcox.rank, all.ps))$r.squared,2))))
with(all.ps,plot(trans.binom.rank,trans.wilcox.rank, pch=".",
                 main=paste("binomial trans vs. wilcox trans \n r2 =", round(summary(lm(trans.binom.rank~trans.wilcox.rank, all.ps))$r.squared,2))))
dev.off()
