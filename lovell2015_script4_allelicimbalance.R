# Final Differential expression and cis-trans pipeline.
# Part 4: Compare allelic imbalance between total counts of hal and fil in wet and dry conditions
# Author: JT Lovell, S Schwartz
# Date: 17-April 2015
# Version: 5.2

####################
####################
####################
####################
# 4. Prepare datasets
rm(list=ls())
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/ph2015_eqtl")
pkg <- c("RCurl","plyr","mclust","qtl","DESeq2","GenomicRanges","car","MASS")
invisible(lapply(pkg, function(x) {cat(x,"..."); library(x, character.only=T, verbose=F, warn.conflicts=F,quietly=T)} ))
sessionInfo()

####################
# 4.1: Load DESeq objects and run DESeq
load("ct2015final_deseqinput.RData")
des.all1 <- DESeq(dds.all1, test="Wald")
des.all2 <- DESeq(dds.all2, test="Wald")
resultsNames(des.all2)

####################
# 4.2: Extract results from the DESeq model fit- g + e + gxe
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
resultsNames(des.all1)
res.all1<-DESeqResults(DESeqobj=des.all1,
                          annotation=all.annot, 
                          contrasts=sapply(resultsNames(des.all1), list),
                          prefix=c("Intercept","time","trt","geno","genotrt"))
####################
# 4.3 make main significance table
ms<-res.all1[,grep(".sig",colnames(res.all1))][,3:5]
ms$trt.sig[ms$trt.sig=="sig"]<-"E"
ms$geno.sig[ms$geno.sig=="sig"]<-"G"
ms$genotrt.sig[ms$genotrt.sig=="sig"]<-"GxE"
ms<-paste(ms[,1],ms[,2],ms[,3],sep="")
ms[grep("GxE",ms)]<-"GxE"
tab<-table(ms)
res.all1$main.sig<-ms

# 4.4: Create a log2FoldChange dataset
hvf_dry<-results(des.all2,contrast=list("trtgenoDry_FIL2","trtgenoDry_HAL2"))
hvf_wet<-results(des.all2,contrast=list("trtgenoWet_FIL2","trtgenoWet_HAL2"))
wvd_hal<-results(des.all2,contrast=list("trtgenoWet_HAL2","trtgenoDry_HAL2"))
wvd_fil<-results(des.all2,contrast=list("trtgenoWet_FIL2","trtgenoDry_FIL2"))
lfcs<-data.frame(rownames(hvf_dry),hvf_dry[,"log2FoldChange"],hvf_wet[,"log2FoldChange"])
lfcs2<-data.frame(rownames(hvf_dry),hvf_dry[,"log2FoldChange"],hvf_wet[,"log2FoldChange"],
                  wvd_hal[,"log2FoldChange"],wvd_fil[,"log2FoldChange"],
                  hvf_dry[,"padj"],hvf_wet[,"padj"],
                  wvd_hal[,"padj"],wvd_fil[,"padj"])
colnames(lfcs2)<-c("id","hvf.dry_lfc","hvf.wet_lfc","wvd_hal_lfc","wvd_fil_lfc","hvf.dry_padj","hvf.wet_padj","wvd_hal_padj","wvd_fil_padj")

# 4.4: Categorize by significance
lfcs2$cat<-"conserved"
lfcs2$cat[lfcs2$hvf.dry_padj<=0.05 & lfcs2$hvf.wet_padj<=0.05 & lfcs2$wvd_hal_padj<=0.05 & lfcs2$wvd_fil_padj<=0.05]<-"all.sig"
lfcs2$cat[lfcs2$hvf.dry_padj<=0.05 & lfcs2$hvf.wet_padj<=0.05 & lfcs2$wvd_hal_padj>0.05 & lfcs2$wvd_fil_padj>0.05]<-"gen in wet&dry"
lfcs2$cat[lfcs2$hvf.dry_padj<=0.05 & lfcs2$hvf.wet_padj>0.05 & lfcs2$wvd_hal_padj>0.05 & lfcs2$wvd_fil_padj>0.05]<-"gen in one"
lfcs2$cat[lfcs2$hvf.dry_padj>0.05 & lfcs2$hvf.wet_padj<=0.05 & lfcs2$wvd_hal_padj>0.05 & lfcs2$wvd_fil_padj>0.05]<-"gen in one"

lfcs2$cat[lfcs2$hvf.dry_padj>0.05 & lfcs2$hvf.wet_padj>0.05 & lfcs2$wvd_hal_padj<=0.05 & lfcs2$wvd_fil_padj<=0.05]<-"trt in hal&fil"
lfcs2$cat[lfcs2$hvf.dry_padj>0.05 & lfcs2$hvf.wet_padj>0.05 & lfcs2$wvd_hal_padj<=0.05 & lfcs2$wvd_fil_padj>0.05]<-"trt in hal"
lfcs2$cat[lfcs2$hvf.dry_padj>0.05 & lfcs2$hvf.wet_padj>0.05 & lfcs2$wvd_hal_padj>0.05 & lfcs2$wvd_fil_padj<=0.05]<-"trt in fil"

lfcs2$cat[lfcs2$hvf.dry_padj>0.05 & lfcs2$hvf.wet_padj<=0.05 & lfcs2$wvd_hal_padj<=0.05 & lfcs2$wvd_fil_padj<=0.05]<-"both trt one gen"
lfcs2$cat[lfcs2$hvf.dry_padj<=0.05 & lfcs2$hvf.wet_padj>0.05 & lfcs2$wvd_hal_padj<=0.05 & lfcs2$wvd_fil_padj<=0.05]<-"both trt one gen"
lfcs2$cat[lfcs2$hvf.dry_padj<=0.05 & lfcs2$hvf.wet_padj<=0.05 & lfcs2$wvd_hal_padj>0.05 & lfcs2$wvd_fil_padj<=0.05]<-"both gen one trt"
lfcs2$cat[lfcs2$hvf.dry_padj<=0.05 & lfcs2$hvf.wet_padj<=0.05 & lfcs2$wvd_hal_padj<=0.05 & lfcs2$wvd_fil_padj>0.05]<-"both gen one trt"

lfcs2$cat[lfcs2$hvf.dry_padj>0.05 & lfcs2$hvf.wet_padj<=0.05 & lfcs2$wvd_hal_padj>0.05 & lfcs2$wvd_fil_padj<=0.05]<-"opposites"
lfcs2$cat[lfcs2$hvf.dry_padj>0.05 & lfcs2$hvf.wet_padj<=0.05 & lfcs2$wvd_hal_padj<=0.05 & lfcs2$wvd_fil_padj>0.05]<-"opposites"
lfcs2$cat[lfcs2$hvf.dry_padj<=0.05 & lfcs2$hvf.wet_padj>0.05 & lfcs2$wvd_hal_padj>0.05 & lfcs2$wvd_fil_padj<=0.05]<-"opposites"
lfcs2$cat[lfcs2$hvf.dry_padj<=0.05 & lfcs2$hvf.wet_padj>0.05 & lfcs2$wvd_hal_padj<=0.05 & lfcs2$wvd_fil_padj>0.05]<-"opposites"
lfcs2$cat2<-"conserved"

# 4.5: Categorize by direction
lfcs2$cat2[lfcs2$hvf.dry_lfc<0 & lfcs2$hvf.wet_lfc<0]<-"sig in wet & dry, same direction"
lfcs2$cat2[lfcs2$hvf.dry_lfc>=0 & lfcs2$hvf.wet_lfc>=0]<-"sig in wet & dry, same direction"
lfcs2$cat2[lfcs2$hvf.dry_lfc<0 & lfcs2$hvf.wet_lfc>=0]<-"sig in wet & dry, opposite direction"
lfcs2$cat2[lfcs2$hvf.dry_lfc>=0 & lfcs2$hvf.wet_lfc<0]<-"sig in wet & dry, opposite direction"
lfcs2$cat2[lfcs2$hvf.dry_padj>0.05 & lfcs2$hvf.wet_padj>0.05]<-"conserved"
lfcs2$cat2[lfcs2$hvf.dry_padj>0.05 & lfcs2$hvf.wet_padj<=0.05]<-"only wet"
lfcs2$cat2[lfcs2$hvf.dry_padj<=0.05 & lfcs2$hvf.wet_padj>0.05]<-"only dry"
lfcs3<-lfcs2
lfcs3$hvf.wet_lfc[lfcs3$hvf.wet_lfc>5]<-5
lfcs3$hvf.dry_lfc[lfcs3$hvf.dry_lfc>5]<-5
lfcs3$hvf.wet_lfc[lfcs3$hvf.wet_lfc<(-5)]<-(-5)
lfcs3$hvf.dry_lfc[lfcs3$hvf.dry_lfc<(-5)]<-(-5)

#4.6: Plot allelic imbalance
pdf("ct2015_directionandsignificance.pdf")
ggplot(lfcs3, aes(x=hvf.wet_lfc,y=hvf.dry_lfc, col=cat2))+
  geom_vline(xintercept = 0)+
  geom_hline(xintercept = 0)+
  geom_point(alpha=.8, size=1.5)+
  scale_color_manual(name="significance and direction",
                     values=c("grey","red","blue","green","pink"),
                     labels=paste(names(table(lfcs2$cat2)),"\n n =",as.numeric(table(lfcs2$cat2))))+
  guides(colour = guide_legend(override.aes = list(size=3, alpha=1)))+
  theme_classic()+
  theme(legend.justification = c(0, 1), legend.position = c(0, 1),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())+
  scale_x_continuous("Log2 Fold Change between Hal and Fil (wet)")+
  scale_y_continuous("Log2 Fold Change between Hal and Fil (dry)")+
  ggtitle("Direction and significance of diff. expr.  Hal vs. Fil")
dev.off()

lfc_gxe<-merge(lfcs3,res.all1, by="id")
tab<-t(with(lfc_gxe,table(main.sig,cat2)))
write.csv(lfc_gxe, file="ct2015_gxecats.csv", row.names=F)
names(tab)[1]<-"conserved"
x <- barplot(tab, 
             col=rep(c("grey","red","blue","green","pink"),5),
             main="distribution of term significance across genes",
             sub="where trt interaction is sig. gene binned to int, not main", border = NA,
             horiz=T)
text(cex=.8, x=0, y=x, names(tab), xpd=TRUE, pos=2,srt=0)
text(cex=.8, x=tab+200, y=x, tab, xpd=TRUE)

barplot(tab2, )
legend()
tab<-table(lfcs2$cat2)
lfcs2$cat3<-lfcs2$cat2
lfcs2$cat3<-with(lfcs2, ifelse(cat2=="sig in wet & dry, same direction" & hvf.wet_lfc>0,"same_direction_up",
                               ifelse(cat2=="sig in wet & dry, same direction" & hvf.wet_lfc<0,"same_direction_down",
                                      ifelse(cat2=="only wet" & hvf.wet_lfc>0,"only_wet_up",
                                             ifelse(cat2=="only wet" & hvf.wet_lfc<0,"only_wet_down",
                                                    ifelse(cat2=="only dry" & hvf.dry_lfc>0,"only_dry_up",
                                                           ifelse(cat2=="only dry" & hvf.dry_lfc<0,"only_dry_down",
                                                                  ifelse(cat2=="sig in wet & dry, opposite direction" & hvf.dry_lfc>0,"opposite_direction_up",
                                                                         ifelse(cat2=="sig in wet & dry, opposite direction" & hvf.dry_lfc<0,"opposite_direction_down",NA)))))))))

# 4.7: Plot category membership
tab<-table(lfcs2$cat3)
par(mar=c(5, 10, 4, 2) + 0.1)
names(tab)<-gsub(",","\n",names(tab))
pdf("ct2015_trtcatbar.pdf")
x <- barplot(tab,  
             col=c("grey","red","blue","green","pink"), horiz=T,
             main="Direction and significance of diff. expr.  Hal vs. Fil",
             xlim=c(0,12000), xaxt="n", las=1)
text(cex=.8, x=tab, y=x, tab, xpd=TRUE, pos=4)
dev.off()

# 4.8: Write stats file with gene categorizations
fig2.cats<-lfcs2[,c("id","cat","cat2","cat3")]
colnames(fig2.cats)<-c("id","fig2_initialcat","fig2_colorcat","fig2_directioncat")
write.csv(fig2.cats, file="ct2015_genecategories_fig2.csv",row.names=F)
