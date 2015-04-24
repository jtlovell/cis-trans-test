# Final Differential expression and cis-trans pipeline.
# Part 5: Test mode of inheritance for each gene. 
# Author: JT Lovell, S Schwartz
# Date: 17-April 2015
# Version: 5.2

####################
####################
####################
####################
# 5. Prepare datasets
rm(list=ls())
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/ph2015_eqtl")
pkg <- c("RCurl","plyr","mclust","qtl","DESeq2","GenomicRanges","car","MASS","SNPassoc","lsmeans","ggtern","dplyr","reshape2")
invisible(lapply(pkg, function(x) {cat(x,"..."); library(x, character.only=T, verbose=F, warn.conflicts=F,quietly=T)} ))
sessionInfo()

####################
# 5.1: Load DESeq objects and run DESeq
load("ct2015final_deseqinput.RData")

des.ase1 <- DESeq(dds.ase1, test="Wald")
des.ase.wet <- DESeq(dds.ase.wet, test="Wald")
des.ase.dry <- DESeq(dds.ase.dry, test="Wald")
####################

# Function to extract DESeq Results
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
resultsNames(des.ase.wet)
# 5.2 Get stats for within environment tests
res.ase.wet<-DESeqResults(DESeqobj=des.ase.wet,
                       annotation=ase.annot, 
                       contrasts=sapply(c("allele_hal_vs_fil","allelehal.generationF0"), list),
                       prefix=c("cis","trans"))
res.ase.dry<-DESeqResults(DESeqobj=des.ase.dry,
                       annotation=ase.annot, 
                       contrasts=sapply(c("allele_hal_vs_fil","allelehal.generationF0"), list),
                       prefix=c("cis","trans"))
asewet<-res.ase.wet
asedry<-res.ase.dry

asewet$main.sig<-with(asewet,ifelse(cis.padj>0.05 & trans.padj>0.05,"nosig",
                                    ifelse(cis.padj<=0.05 & trans.padj>0.05,"cis",
                                           ifelse(cis.padj>0.05 & trans.padj<=0.05,"trans",
                                                  ifelse(cis.padj<=0.05 & trans.padj<=0.05 ,"cis+trans",NA)))))
asedry$main.sig<-with(asedry,ifelse(cis.padj>0.05 & trans.padj>0.05,"nosig",
                                    ifelse(cis.padj<=0.05 & trans.padj>0.05,"cis",
                                           ifelse(cis.padj>0.05 & trans.padj<=0.05,"trans",
                                                  ifelse(cis.padj<=0.05 & trans.padj<=0.05 ,"cis+trans",NA)))))

####################
# 5.3 calculate LFC by hand within generations and treatments...
#wet
vsd.wet<-data.frame(t(counts(des.ase.wet,normalize=T)))
ct.info.wet<-data.frame(ase.info.wet) #info for counts
vsd.trans.wet<-data.frame(ct.info.wet,vsd.wet) #combine info and qn counts
vsd.trans.wet<-vsd.trans.wet[,c("allele","generation",as.character(colnames(vsd.wet)))] #drop unneeded columns
means.all.wet<-vsd.trans.wet %>% group_by(allele,generation) %>% summarise_each(funs(mean))
f1.lfc.wet<-log2(as.numeric(means.all.wet[1,-c(1:2)]))-log2(as.numeric(means.all.wet[3,-c(1:2)]))
f0.lfc.wet<-log2(as.numeric(means.all.wet[2,-c(1:2)]))-log2(as.numeric(means.all.wet[4,-c(1:2)]))
#dry
vsd.dry<-data.frame(t(counts(des.ase.dry,normalize=T))) #qn counts
ct.info.dry<-data.frame(ase.info.dry) #info for counts
vsd.trans.dry<-data.frame(ct.info.dry,vsd.dry) #combine info and qn counts
vsd.trans.dry<-vsd.trans.dry[,c("allele","generation",as.character(colnames(vsd.dry)))] #drop unneeded columns
means.all.dry<-vsd.trans.dry %>% group_by(allele,generation) %>% summarise_each(funs(mean))
f1.lfc.dry<-log2(as.numeric(means.all.dry[1,-c(1:2)]))-log2(as.numeric(means.all.dry[3,-c(1:2)]))
f0.lfc.dry<-log2(as.numeric(means.all.dry[2,-c(1:2)]))-log2(as.numeric(means.all.dry[4,-c(1:2)]))

lfcs<-data.frame(colnames(vsd.wet),f1.lfc.wet,f0.lfc.wet,f1.lfc.dry,f0.lfc.dry)
colnames(lfcs)[1]<-"id"
toplot.imb<-merge(lfcs,asedry[,c("id","main.sig")], by="id")
colnames(toplot.imb)[6]<-"main.sig.dry"
toplot.imb<-merge(toplot.imb,asewet[,c("id","main.sig")], by="id")
colnames(toplot.imb)[7]<-"main.sig.wet"

####################
#pdf("ct2015_allctplots.pdf")
# 5.4 make imbalance plots:
toplot.imb$f1.lfc.wet[toplot.imb$f1.lfc.wet>4]<-4
toplot.imb$f1.lfc.wet[toplot.imb$f1.lfc.wet<=(-4)]<-(-4)
toplot.imb$f0.lfc.wet[toplot.imb$f0.lfc.wet>4]<-4
toplot.imb$f0.lfc.wet[toplot.imb$f0.lfc.wet<=(-4)]<-(-4)
toplot.imb$f1.lfc.dry[toplot.imb$f1.lfc.dry>4]<-4
toplot.imb$f1.lfc.dry[toplot.imb$f1.lfc.dry<=(-4)]<-(-4)
toplot.imb$f0.lfc.dry[toplot.imb$f0.lfc.dry>4]<-4
toplot.imb$f0.lfc.dry[toplot.imb$f0.lfc.dry<=(-4)]<-(-4)
#postscript("ct2015_standardimbalanceplots_wet.pdf", height=8, width=8)
ggplot(toplot.imb, aes(x=f1.lfc.wet, y=f0.lfc.wet, col=main.sig.wet))+
  geom_vline(xintercept = 0)+
  geom_hline(xintercept = 0)+
  geom_point(alpha=.5, size=1)+theme_classic()+
  theme(legend.justification = c(0, 1), legend.position = c(0, 1),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())
ggplot(toplot.imb, aes(x=f1.lfc.dry, y=f0.lfc.dry, col=main.sig.dry))+
  geom_vline(xintercept = 0)+
  geom_hline(xintercept = 0)+
  geom_point(alpha=.5, size=1)+theme_classic()+
  theme(legend.justification = c(0, 1), legend.position = c(0, 1),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())
#dev.off()

####################
# 5.5 calculate significance groups, cis, trans, cis+trans, cis*trans, not significant
res.ase1<-DESeqResults(DESeqobj=des.ase1,
                       annotation=ase.annot, 
                       contrasts=sapply(resultsNames(des.ase1), list),
                       prefix=c("Intercept","time","trt","allele","gen","allele.gen","trt.allele","trt.gen","trt.allele.gen"))

####################
# 5.6  Rename columns, organize output dataframe
mstats<-res.ase1
colnames(mstats)<-gsub("allele.gen","trans",colnames(mstats))
colnames(mstats)<-gsub("allele","cis",colnames(mstats))
colnames(mstats)<-gsub("log2FoldChange","lfc", colnames(mstats))
colnames(mstats)<-gsub(".padj","_padj", colnames(mstats))
colnames(mstats)<-gsub(".sig","_sig", colnames(mstats))
colnames(mstats)<-gsub(".pvalue","_pvalue", colnames(mstats))
colnames(mstats)<-gsub(".lfc","_lfc", colnames(mstats))
colnames(mstats)<-gsub(".stat","_stat", colnames(mstats))
mstats<-mstats[-grep("_sig", colnames(mstats)),]
toplot<-mstats[complete.cases(mstats),]

####################
# 5.7 Make "main.sig" column with significance categories
#sig cis:
toplot$cistrt.transtrt_sig<-with(toplot,ifelse(trt.trans_padj<=0.05 & trt.cis_padj<=0.05,"cistrt.transtrt",0))
toplot$transtrt_sig<-with(toplot,ifelse(trt.trans_padj<=0.05 & trt.cis_padj>0.05,"transtrt",0))
toplot$cistrt_sig<-with(toplot,ifelse(cistrt.transtrt_sig==0 & transtrt_sig==0 & trt.cis_padj<=0.05,"cistrt",0))
toplot$cis.trans_sig<-with(toplot,ifelse(cistrt_sig==0 & transtrt_sig==0 & cistrt.transtrt_sig==0 & cistrt.transtrt_sig==0 & trans_padj<0.05 & cis_padj<0.05,"cisandtrans",0))
toplot$trans_sig<-with(toplot,ifelse(transtrt_sig==0 & cistrt.transtrt_sig==0 & cis.trans_sig==0 & trans_padj<0.05,"trans",0))
toplot$cis_sig<-with(toplot,ifelse(cistrt_sig==0 & cistrt.transtrt_sig==0 & cis.trans_sig==0 & cis_padj<0.05,"cis",0))

toplot$main.sig<-with(toplot, paste(cistrt.transtrt_sig,transtrt_sig,cistrt_sig,cis.trans_sig,trans_sig,cis_sig, sep=""))
toplot$main.sig<-gsub("0","",toplot$main.sig)
toplot$main.sig[toplot$main.sig==""]<-"nosig"
toplot$main.sig[toplot$main.sig=="nosig" & toplot$trt_padj<0.05]<-"trtonly"
table(toplot$main.sig)

####################
# 5.8 Make barplot of mainsig categories
tab<-table(toplot$main.sig)
sum(toplot$trt_sig=="trt")
tab<-tab[order(-tab)]
test<-as.numeric(tab)
x <- barplot(tab, xaxt="n", ylab="n. significant genes/term", 
             main="distribution of term significance across genes",
             sub="where trt interaction is sig. gene binned to int, not main", border = NA)
text(cex=.8, x=x+.25, y=-200, names(tab), xpd=TRUE, srt=45, pos=2)
text(cex=.8, x=x, y=tab+200, tab, xpd=TRUE)

####################
# 5.10: Prep data for analysis of Cis * treatment interactions categorization
quantnorm<-function(x) {
  n=sum(!is.na(x),na.rm=T)
  x=rank(x)/(n+1)
  x=qnorm(x)
  x[is.infinite(x)]=NA
  x
}
library(dplyr)
vsd<-data.frame(counts(des.ase1,normalize=T)) #qn counts
ct.info<-data.frame(ase.info) #info for counts
qv<-apply(vsd,1,quantnorm) #quantile norm counts
vsd.cistrt<-data.frame(ct.info,qv) #combine info and qn counts
vsd.cistrt<-vsd.cistrt[,c("Treatment","allele",as.character(rownames(vsd)))] #drop unneeded columns
means.all<-vsd.cistrt %>% group_by(Treatment,allele) %>% summarise_each(funs(mean))
means.all<-data.frame(means.all)

####################
# 5.11 groupings - cis*trt
means.all.sub<-means.all
means.all.sub<-means.all.sub[,as.character(toplot$id[toplot$trt.cis_padj<0.05])]
m.fd<-as.numeric(means.all.sub[1,])
m.hd<-as.numeric(means.all.sub[2,])
m.fw<-as.numeric(means.all.sub[3,])
m.hw<-as.numeric(means.all.sub[4,])
s.f<-m.fd-m.fw
s.h<-m.hd-m.hw
cistrt.cats<-list()
cistrt.cats[["g1.1"]]<-colnames(means.all.sub)[m.fd>=m.hd & m.fw>=m.hw & s.f>=0 & s.h>=0 & abs(s.f)>=abs(s.h)]
cistrt.cats[["g1.2"]]<-colnames(means.all.sub)[m.fd>=m.hd & m.fw>=m.hw & s.f<0 & s.h<0 & abs(s.f)>=abs(s.h)]
cistrt.cats[["g1.3"]]<-colnames(means.all.sub)[m.fd<m.hd & m.fw<m.hw & s.f>=0 & s.h>=0 & abs(s.f)>=abs(s.h)]
cistrt.cats[["g1.4"]]<-colnames(means.all.sub)[m.fd<m.hd & m.fw<m.hw & s.f<0 & s.h<0 & abs(s.f)>=abs(s.h)]
cistrt.cats[["g1.5"]]<-colnames(means.all.sub)[m.fd>=m.hd & m.fw>=m.hw & s.f>=0 & s.h>=0 & abs(s.f)<abs(s.h)]
cistrt.cats[["g1.6"]]<-colnames(means.all.sub)[m.fd>=m.hd & m.fw>=m.hw & s.f<0 & s.h<0 & abs(s.f)<abs(s.h)]
cistrt.cats[["g1.7"]]<-colnames(means.all.sub)[m.fd<m.hd & m.fw<m.hw & s.f>=0 & s.h>=0 & abs(s.f)<abs(s.h)]
cistrt.cats[["g1.8"]]<-colnames(means.all.sub)[m.fd<m.hd & m.fw<m.hw & s.f<0 & s.h<0 & abs(s.f)<abs(s.h)]

cistrt.cats[["g2.1"]]<-colnames(means.all.sub)[m.fd>=m.hd & m.fw<m.hw & s.f>=0 & s.h>=0]
cistrt.cats[["g2.2"]]<-colnames(means.all.sub)[m.fd<m.hd & m.fw>=m.hw & s.f>=0 & s.h>=0]
cistrt.cats[["g2.3"]]<-colnames(means.all.sub)[m.fd>=m.hd & m.fw<m.hw & s.f<0 & s.h<0]
cistrt.cats[["g2.4"]]<-colnames(means.all.sub)[m.fd<m.hd & m.fw>=m.hw & s.f<0 & s.h<0]

cistrt.cats[["g3.1"]]<-colnames(means.all.sub)[m.fd>=m.hd & m.fw>=m.hw & s.f>=0 & s.h<0]
cistrt.cats[["g3.2"]]<-colnames(means.all.sub)[m.fd>=m.hd & m.fw>=m.hw & s.f<0 & s.h>=0]
cistrt.cats[["g3.3"]]<-colnames(means.all.sub)[m.fd<m.hd & m.fw<m.hw & s.f>=0 & s.h<0]
cistrt.cats[["g3.4"]]<-colnames(means.all.sub)[m.fd<m.hd & m.fw<m.hw & s.f<0 & s.h>=0]

cistrt.cats[["g4.1"]]<-colnames(means.all.sub)[m.fd>=m.hd & m.fw<m.hw & s.f>=0 & s.h<0]
cistrt.cats[["g4.2"]]<-colnames(means.all.sub)[m.fd<m.hd & m.fw>=m.hw & s.f<0 & s.h>=0]
ctcat<-as.character(toplot$id)
for(i in 1:length(cistrt.cats)){
  x<-cistrt.cats[[i]]; ctcat[ctcat %in% x]<-names(cistrt.cats)[i]
}
ctcat[grep("Pahal",ctcat)]<-""
toplot$cistrt.group<-ctcat

####################
# 5.12: Prep data for analysis of trans interactions categorization
vsd.trans<-data.frame(ct.info,qv) #combine info and qn counts
vsd.trans<-vsd.trans[,c("allele","generation",as.character(rownames(vsd)))] #drop unneeded columns
means.all<-vsd.trans %>% group_by(allele,generation) %>% summarise_each(funs(mean))
means.all<-data.frame(means.all)

####################
# 5.13: Make trans interactions categorization
means.all.sub<-means.all
means.all.sub<-means.all.sub[,as.character(toplot$id[toplot$trans_padj<0.05])]
m.f0f<-as.numeric(means.all.sub[1,])
m.f0h<-as.numeric(means.all.sub[3,])
m.f1f<-as.numeric(means.all.sub[2,])
m.f1h<-as.numeric(means.all.sub[4,])
s.f<-m.f0f-m.f1f
s.h<-m.f0h-m.f1h

trans.cats<-list()
trans.cats[["g1.1"]]<-colnames(means.all.sub)[m.f0f>=m.f0h & m.f1f>=m.f1h & s.f>=0 & s.h>=0 & abs(s.f)>=abs(s.h)]
trans.cats[["g1.2"]]<-colnames(means.all.sub)[m.f0f>=m.f0h & m.f1f>=m.f1h & s.f<0 & s.h<0 & abs(s.f)>=abs(s.h)]
trans.cats[["g1.3"]]<-colnames(means.all.sub)[m.f0f<m.f0h & m.f1f<m.f1h & s.f>=0 & s.h>=0 & abs(s.f)>=abs(s.h)]
trans.cats[["g1.4"]]<-colnames(means.all.sub)[m.f0f<m.f0h & m.f1f<m.f1h & s.f<0 & s.h<0 & abs(s.f)>=abs(s.h)]
trans.cats[["g1.5"]]<-colnames(means.all.sub)[m.f0f>=m.f0h & m.f1f>=m.f1h & s.f>=0 & s.h>=0 & abs(s.f)<abs(s.h)]
trans.cats[["g1.6"]]<-colnames(means.all.sub)[m.f0f>=m.f0h & m.f1f>=m.f1h & s.f<0 & s.h<0 & abs(s.f)<abs(s.h)]
trans.cats[["g1.7"]]<-colnames(means.all.sub)[m.f0f<m.f0h & m.f1f<m.f1h & s.f>=0 & s.h>=0 & abs(s.f)<abs(s.h)]
trans.cats[["g1.8"]]<-colnames(means.all.sub)[m.f0f<m.f0h & m.f1f<m.f1h & s.f<0 & s.h<0 & abs(s.f)<abs(s.h)]

trans.cats[["g2.1"]]<-colnames(means.all.sub)[m.f0f>=m.f0h & m.f1f<m.f1h & s.f>=0 & s.h>=0]
trans.cats[["g2.2"]]<-colnames(means.all.sub)[m.f0f<m.f0h & m.f1f>=m.f1h & s.f>=0 & s.h>=0]
trans.cats[["g2.3"]]<-colnames(means.all.sub)[m.f0f>=m.f0h & m.f1f<m.f1h & s.f<0 & s.h<0]
trans.cats[["g2.4"]]<-colnames(means.all.sub)[m.f0f<m.f0h & m.f1f>=m.f1h & s.f<0 & s.h<0]

trans.cats[["g3.1"]]<-colnames(means.all.sub)[m.f0f>=m.f0h & m.f1f>=m.f1h & s.f>=0 & s.h<0]
trans.cats[["g3.2"]]<-colnames(means.all.sub)[m.f0f>=m.f0h & m.f1f>=m.f1h & s.f<0 & s.h>=0]
trans.cats[["g3.3"]]<-colnames(means.all.sub)[m.f0f<m.f0h & m.f1f<m.f1h & s.f>=0 & s.h<0]
trans.cats[["g3.4"]]<-colnames(means.all.sub)[m.f0f<m.f0h & m.f1f<m.f1h & s.f<0 & s.h>=0]

trans.cats[["g4.1"]]<-colnames(means.all.sub)[m.f0f>=m.f0h & m.f1f<m.f1h & s.f>=0 & s.h<0]
trans.cats[["g4.2"]]<-colnames(means.all.sub)[m.f0f<m.f0h & m.f1f>=m.f1h & s.f<0 & s.h>=0]
transcat<-as.character(toplot$id)
for(i in 1:length(trans.cats)){
  x<-trans.cats[[i]]; transcat[transcat %in% x]<-names(trans.cats)[i]
}
transcat[grep("Pahal",transcat)]<-""
toplot$trans.group<-transcat

####################
# 5.14: Remake means for each interaction and merge with stats
vsd<-data.frame(counts(des.ase1,normalize=T)) #qn counts
ct.info<-data.frame(ase.info) #info for counts
qv<-apply(vsd,1,quantnorm) #quantile norm counts
vsd.cistrt<-data.frame(ct.info,qv) #combine info and qn counts
vsd.cistrt<-vsd.cistrt[,c("Treatment","allele",as.character(rownames(vsd)))] #drop unneeded columns
means.all<-vsd.cistrt %>% group_by(Treatment,allele) %>% summarise_each(funs(mean))
ct.means<-data.frame(t(means.all[,-c(1:2)]))
colnames(ct.means)<-c("ctmeans_fil.dry","ctmeans_hal.dry","ctmeans_fil.wet","ctmeans_hal.wet")
ct.means$id<-rownames(ct.means)

vsd<-data.frame(counts(des.ase1,normalize=T)) #qn counts
ct.info<-data.frame(ase.info) #info for counts
qv<-apply(vsd,1,quantnorm) #quantile norm counts
vsd.cistrt<-data.frame(ct.info,qv) #combine info and qn counts
vsd.cistrt<-vsd.cistrt[,c("allele","generation",as.character(rownames(vsd)))] #drop unneeded columns
means.all<-vsd.cistrt %>% group_by(allele,generation) %>% summarise_each(funs(mean))
trans.means<-data.frame(t(means.all[,-c(1:2)]))
colnames(trans.means)<-c("transmeans_fil.F0","transmeans_fil.F1","transmeans_hal.F0","transmeans_hal.F1")
trans.means$id<-rownames(trans.means)

all.means<-merge(trans.means,ct.means, by="id")
toplot2<-merge(toplot,all.means, by="id")
toplot2$cistrt.main<-substr(toplot2$cistrt.group,2,2)
toplot2$cistrt.sub<-substr(toplot2$cistrt.group,4,4)
toplot2$trans.main<-substr(toplot2$trans.group,2,2)
toplot2$trans.sub<-substr(toplot2$trans.group,4,4)

####################
# 5.15 Make all sub plots for the final figure
#pdf("ct2015_allfig3subplots.pdf")
tab<-table(toplot2$main.sig)
par(mar=c(4,2,2,1)+.1)
# a) barplot
x <- barplot(tab,  
             col=rainbow(10), horiz=F, xaxt="n" )
text(cex=.8, x=x, y=tab, tab, xpd=TRUE, pos=3)
text(cex=.8, x=x, y=-100, names(tab), xpd=TRUE, pos=2, srt=45)
colnames(toplot2)<-gsub("log2FoldChange","lfc",colnames(toplot2))
toplot2$totint.var<-with(toplot2,abs(cis_lfc)+abs(trt_lfc)+abs(trans_lfc))
toplot2$cis_propvar<-with(toplot2, abs(cis_lfc)/totint.var)
toplot2$trt_propvar<-with(toplot2, abs(trt_lfc)/totint.var)
toplot2$trans_propvar<-with(toplot2, abs(trans_lfc)/totint.var)
tern2<-toplot2[union(grep("cis",toplot2$main.sig),grep("trans",toplot2$main.sig)),]
table(tern2$main.sig)
# b) ternary plot
ggtern(data=tern2,aes(x=trt_propvar,y=cis_propvar,z=trans_propvar)) + 
  geom_point(aes(col=main.sig),size=1.5, alpha=1) + 
  scale_color_manual(values=rainbow(10)[-c(6,10)])+
  geom_density2d(color="black",linetype=1, alpha=.8) +
  theme_rgbw()+  theme_notitles() + theme(legend.key = element_blank())+
  guides(colour = guide_legend(override.aes = list(size=3, alpha=1)))+
  ggtitle("comparison of trans, trt and interaction effects")

# c) cis*trt barplot
long2<-melt(toplot2[grep("g",toplot2$cistrt.group),
                    c("id","ctmeans_fil.dry","ctmeans_hal.dry","ctmeans_fil.wet","ctmeans_hal.wet","cistrt.main","cistrt.sub")],
            , idvars=c("id","cistrt.main","cistrt.sub"))
long2$variable<-gsub("ctmeans_","",long2$variable)
long2$treatment<-sapply(long2$variable,function(x) strsplit(x,"[.]")[[1]][2])
long2$allele<-sapply(long2$variable,function(x) strsplit(x,"[.]")[[1]][1])
long2$cistrt.sub[long2$cistrt.main==1 & long2$cistrt.sub==3]<-1
long2$cistrt.sub[long2$cistrt.main==1 & long2$cistrt.sub==4]<-2
long2$cistrt.sub[long2$cistrt.main==1 & long2$cistrt.sub==5 | long2$cistrt.sub==7]<-3
long2$cistrt.sub[long2$cistrt.main==1 & long2$cistrt.sub==6 | long2$cistrt.sub==8]<-4
cistrt.unique<-paste(long2$cistrt.main,long2$cistrt.sub, sep="_")
tab<-table(cistrt.unique)/4
par(mar=c(4,2,2,1)+.1)
x <- barplot(tab,  
             col=c(rep("red",4),rep("green",4),rep("blue",4), rep("cyan",2)), horiz=F, xaxt="n" )
text(cex=.8, x=x, y=tab, tab, xpd=TRUE, pos=3)
text(cex=.8, x=x, y=-1, names(tab), xpd=TRUE, pos=1, srt=45)
long.cistrt<-long2[!duplicated(long2$id),]
# d) cis*trt line
ggplot(long2, aes(x=treatment, y=value, col=allele, group=interaction(allele,id)))+
  geom_line(alpha=0.3)+facet_grid(cistrt.sub~cistrt.main)+theme_bw()

# e) trans interaction bar plots
long2<-melt(toplot2[grep("g",toplot2$trans.group),
                    c("id","transmeans_fil.F0","transmeans_fil.F1","transmeans_hal.F0","transmeans_hal.F1","trans.main","trans.sub")],
            , idvars=c("id"))
long2$variable<-gsub("transmeans_","",long2$variable)
long2$generation<-sapply(long2$variable,function(x) strsplit(x,"[.]")[[1]][2])
long2$allele<-sapply(long2$variable,function(x) strsplit(x,"[.]")[[1]][1])
long2$trans.sub[long2$trans.main==1 & long2$trans.sub==3]<-1
long2$trans.sub[long2$trans.main==1 & long2$trans.sub==4]<-2
long2$trans.sub[long2$trans.main==1 & long2$trans.sub==5 | long2$trans.sub==7]<-3
long2$trans.sub[long2$trans.main==1 & long2$trans.sub==6 | long2$trans.sub==8]<-4
trans.unique<-paste(long2$trans.main,long2$trans.sub, sep="_")
tab<-table(trans.unique)/4
par(mar=c(4,2,2,1)+.1)
x <- barplot(tab,  
             col=c(rep("red",4),rep("green",4),rep("blue",4), rep("cyan",2)), horiz=F, xaxt="n" )
text(cex=.8, x=x, y=tab, tab, xpd=TRUE, pos=3)
text(cex=.8, x=x, y=-10, names(tab), xpd=TRUE, pos=1, srt=45)
# f) trans line plots
ggplot(long2, aes(x=generation, y=value, col=allele, group=interaction(allele,id)))+
  geom_line(alpha=0.3)+facet_grid(trans.sub~trans.main)+theme_bw()

#dev.off()

#combine stats with Fig 2 and dominance stats
gxe.cats<-read.csv("ct2015_gxecats.csv")
fig2.cats<-read.csv("ct2015_genecategories_fig2.csv")
dominance.cats<-read.csv("ct2015_dominancecats.csv")
final.stats<-merge(fig2.cats,dominance.cats, by="id",all.x=T)
final.stats<-merge(final.stats,toplot2, by="id",all.x=T)
write.csv(final.stats, file="ct2015_finalstats.csv", row.names=F)

gxe.dom<-merge(gxe.cats,dominance.cats, by="id", all.x=T)