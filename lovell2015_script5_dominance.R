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
pkg <- c("RCurl","plyr","mclust","qtl","DESeq2","GenomicRanges","car","MASS","SNPassoc","lsmeans")
invisible(lapply(pkg, function(x) {cat(x,"..."); library(x, character.only=T, verbose=F, warn.conflicts=F,quietly=T)} ))
sessionInfo()

####################
# 5.1: Load DESeq objects and run DESeq
load("ct2015final_deseqinput.RData")
des.all3 <- DESeq(dds.all3, test="Wald")
normalized.counts<-counts(des.all3, normalize=T)
parentgenes<-rownames(des.all3)

# 5.1: Function to calculate dominance etc. 
testMode<-function(normalized.counts, design, form, allele.col="allele", gene){
  #quantnorm<-function(x) {n=sum(!is.na(x),na.rm=T);  x=rank(x)/(n+1);  x=qnorm(x); x[is.infinite(x)]=NA;  x }
  design[,allele.col][design[,allele.col]=="f1"]<-"H/F"
  design[,allele.col][design[,allele.col]=="hal"]<-"H/H"
  design[,allele.col][design[,allele.col]=="fil"]<-"F/F"
  df<-data.frame(normalized.counts[gene,], design)
  colnames(df)[1]<-"vsd"
  df<-setupSNP(df, colSNPs=which(colnames(df)==allele.col), sep = "/")
  ass<-association(formula=as.formula(form), data=df) 
  out<-ass[c(which(rownames(ass)=="Codominant")+1,which(rownames(ass)=="Dominant")+1,
             which(rownames(ass)=="Recessive")+1,which(rownames(ass)=="Overdominant")+1,        
             which(rownames(ass)=="log-Additive")+1),c("p-value","AIC")]
  rownames(out)<-c("Codominant","Dominant","Recessive","Overdominant","log-Additive")
  outnames<-c(paste(rownames(out),"pvalue", sep="."),paste(rownames(out),"AIC", sep="."))
  out<-data.frame(t(data.frame(as.numeric(out)))); colnames(out)<-c(outnames)
  return(out)
}

out<-sapply(parentgenes, function(x){
  testMode(normalized.counts=normalized.counts, 
           design=allF1.info,
           form="vsd ~ timecat + Treatment + allele",
           allele.col="allele", 
           gene=x)
} )

# 5.2: Process the output, generate a table of reaction types. 
dominance<-data.frame(t(out))
save(dominance, file="ct2015_dominancecatuntransformed.RData")
load("ct2015_dominancecat.RData")
head(dominance)
res.parents<-dominance
for(i in colnames(res.parents)[-11]) res.parents[,i]<-as.numeric(unlist(res.parents[,i]))
res.parents$id<-rownames(res.parents)
#res.parents<-res.parents[res.parents$id %in% lfcs2$id[lfcs2$cat2!="conserved"],]
modes<-c("codomominant","dominant","recessive","overdominant","log.additive")
res.parents$best.mode<-apply(res.parents[,6:10],1,function(x) modes[which(x==min(x))][1])
res.parents$min.p<-apply(res.parents[,1:5],1,function(x) x[which(x==min(x))][1])
res.parents$best.mode[res.parents$min.p>0.05]<-"ambiguous"
bm<-table(res.parents$best.mode)
ps<-unique(res.parents$best.mode)

# 5.3: Generate quantile normalized data for line plots (fig 3a)
normalized.counts=counts(des.all3, normalize=T) 
design=allF1.info
form="vsd ~ allele + Treatment + timecat"
allele.col="allele"
quantnorm<-function(x) {n=sum(!is.na(x),na.rm=T);  x=rank(x)/(n+1);  x=qnorm(x); x[is.infinite(x)]=NA;  x }
design$allele<-factor(design$allele)
res.parents$jgi.id<-rownames(res.parents)

out.list<-list()
for(j in c("log.additive", "dominant",     "codomominant", "recessive",    "overdominant")){
  x<-res.parents$jgi.id[res.parents$best.mode==j]
  print(j)
  out.all<-data.frame()
  for(i in 1:length(x)){
    gene=x[i]
    df<-data.frame(quantnorm(normalized.counts[gene,]), design)
    colnames(df)[1]<-"vsd"
    lm.out<-aov(vsd~Treatment + timecat + allele , data=df)
    ls.out<-lsmeans(lm.out,"allele")
    alleles<-rownames(ls.out)
    out<-data.frame(gene,j,data.frame(summary(ls.out)))
    out.all<-rbind(out.all,out)
    print(i)
  }
  out.list[[j]]<-out.all
}
out.all<-ldply(out.list, data.frame)
colnames(out.all)[3]<-"inheritance"
write.csv(out.all, file="ct2015_lsmhalfilf1.notqned.csv", row.names=F)
pdf("ct2015_inheritancemodebar.pdf")
res.parents2<-res.parents
res.parents2$best.mode<-gsub("log.additive","codominant",res.parents2$best.mode)
res.parents2$best.mode<-gsub("codomominant","codominant",res.parents2$best.mode)
bm2<-table(res.parents2$best.mode)
par(mar=c(5, 8, 4, 6) + 0.1)
x <- barplot(bm2,  
             col=rainbow(6), horiz=T,
             main="Inheritance modes",
             xlim=c(0,10000), axes=T, las=1)
text(cex=.8, x=bm2, y=x, bm2, xpd=TRUE, pos=4)
mode.list<-vector()
for(j in c("dominant",     "codominant", "recessive",    "overdominant")){
  sub<-res.parents[res.parents2$best.mode==j,]
  sub<-sub[with(sub, order(min.p)),"jgi.id"]
  sub<-sub[1:100]
  mode.list<-c(mode.list,sub)
}
pdf("ct2015_inheritancerxn.pdf")
inher.toplot<-out.all[out.all$gene %in% mode.list,]
ggplot(inher.toplot, aes(x=allele, y=lsmean, group=gene, col=inheritance)) + geom_line()+
  theme_classic() + facet_grid(inheritance~.) +
  scale_x_discrete(limits=c("fil","f1","hal"))
inher.toplot2<-data.frame(out.all$inheritance[out.all$allele=="fil"],
                          out.all$gene[out.all$allele=="fil"],
                          out.all$lsmean[out.all$allele=="fil"],
                          out.all$lsmean[out.all$allele=="f1"],
                          out.all$lsmean[out.all$allele=="hal"])
colnames(inher.toplot2)<-c("mode","jgi.id","fil.lsmean","f1.lsmean","hal.lsmean")
inher.toplot2$hvf<-inher.toplot2$hal.lsmean-inher.toplot2$fil.lsmean
inher.toplot2$hvf1<-inher.toplot2$hal.lsmean-inher.toplot2$f1.lsmean
inher.toplot2$f1vf<-inher.toplot2$f1.lsmean-inher.toplot2$fil.lsmean
inher.toplot2$mode<-gsub("log.additive","codominant",inher.toplot2$mode)
inher.toplot2$mode<-gsub("codomominant","codominant",inher.toplot2$mode)
resultsNames(des.all3)
fvf1<-results(des.all3,contrast=list("allelefil","allelef1"))
hvf1<-results(des.all3,contrast=list("allelehal","allelef1"))
hvf<-results(des.all3,contrast=list("allelehal","allelefil"))
fvf1<-data.frame(fvf1); fvf1$jgi.id<-row.names(fvf1)
hvf1<-data.frame(hvf1); hvf1$jgi.id<-row.names(hvf1)
hvf<-data.frame(hvf); hvf$jgi.id<-row.names(hvf)
colnames(fvf1)<-paste("fvf1",colnames(fvf1),sep="_")
colnames(hvf1)<-paste("hvf1",colnames(hvf1),sep="_")
colnames(hvf)<-paste("hvf",colnames(hvf),sep="_")
colnames(hvf1)[7]<-"jgi.id"
colnames(fvf1)[7]<-"jgi.id"
colnames(hvf)[7]<-"jgi.id"
all.f1res<-merge(fvf1,hvf1, by="jgi.id")
all.f1res<-merge(res.parents2,all.f1res, by="jgi.id")
#all.f1res$hvf1_log2FoldChange<-(all.f1res$hvf1_log2FoldChange*-1)
pdf("ct2015_modeinher_dotplot.pdf")
ggplot(all.f1res, aes(x=fvf1_log2FoldChange, y=hvf1_log2FoldChange, col=best.mode))+
  geom_vline(xintercept = 0)+
  geom_hline(xintercept = 0)+
  geom_point(alpha=1)+theme_classic()+
  theme(legend.justification = c(0, 1), legend.position = c(0, 1),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())
ggplot(inher.toplot2, aes(x=f1vf, y=hvf1, col=mode))+
  geom_vline(xintercept = 0)+
  geom_hline(xintercept = 0)+
  geom_point(alpha=1)+theme_classic()+
  theme(legend.justification = c(0, 1), legend.position = c(0, 1),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())
dev.off()
colnames(out.all)[1]<-"jgi.id"
out.all2<-merge(out.all,res.parents, by="jgi.id")

out.wide<-cbind(out.all$lsmean[out.all$allele=="F/F"],out.all$lsmean[out.all$allele=="H/F"],out.all$lsmean[out.all$allele=="H/H"])
out.wide = (out.wide-min(out.wide))/(max(out.wide)-min(out.wide))
out.wide<-data.frame(out.wide)
colnames(out.wide)<-c("lsm_fil","lsm_f1","lsm_hal")
ggtern(out.wide, aes(x=lsm_fil,y=lsm_f1,z=lsm_hal))+
  geom_point()
out<-ass[c(which(rownames(ass)=="Codominant")+1,which(rownames(ass)=="Dominant")+1,
           which(rownames(ass)=="Recessive")+1,which(rownames(ass)=="Overdominant")+1,        
           which(rownames(ass)=="log-Additive")+1),c("p-value","AIC")]
rownames(out)<-c("Codominant","Dominant","Recessive","Overdominant","log-Additive")
outnames<-c(paste(rownames(out),"pvalue", sep="."),paste(rownames(out),"AIC", sep="."))
out<-data.frame(t(data.frame(as.numeric(out)))); colnames(out)<-c(outnames)

genes.log.additive<-unique(as.character(res.parents$jgi.id[res.parents$best.mode=="log.additive"]))
genes.overdominant<-unique(as.character(res.parents$jgi.id[res.parents$best.mode=="overdominant"]))
genes.recessive<-unique(as.character(res.parents$jgi.id[res.parents$best.mode=="recessive"]))
genes.dominant<-unique(as.character(res.parents$jgi.id[res.parents$best.mode=="dominant"]))
genes.codominant<-unique(as.character(res.parents$jgi.id[res.parents$best.mode=="codomominant"]))
genes.ambiguous<-unique(as.character(res.parents$jgi.id[res.parents$best.mode=="ambiguous"]))
genes.allde<-unique(as.character(res.parents$jgi.id))
dominance.cats<-res.parents[,c("id","best.mode")]
write.csv(dominance.cats,file="ct2015_dominancecats.csv", row.names=F)
write(genes.log.additive, file="ct2015_genes_modeinher_logadditive.txt", sep="\n")
write(genes.overdominant, file="ct2015_genes_modeinher_overdominant.txt", sep="\n")
write(genes.recessive, file="ct2015_genes_modeinher_recessive.txt", sep="\n")
write(genes.dominant, file="ct2015_genes_modeinher_dominant.txt", sep="\n")
write(genes.codominant, file="ct2015_genes_modeinher_codominant.txt", sep="\n")
write(genes.ambiguous, file="ct2015_genes_modeinher_ambiguous.txt", sep="\n")
write(genes.allde, file="ct2015_genes_modeinher_allde.txt", sep="\n")