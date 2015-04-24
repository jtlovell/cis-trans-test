# Analytical scripts that accompany Lovell et al. (in Review)
# Part 1: Preparation of DESeq2 datasets
# Authors: JT Lovell, S Schwartz
# Date: 17-April 2015
# Version: 5.2

####################
####################
####################
####################
# 1. Prepare datasets
rm(list=ls())
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/ph2015_eqtl")
pkg <- c("RCurl","plyr","mclust","qtl","DESeq2","GenomicRanges","car","MASS")
invisible(lapply(pkg, function(x) {cat(x,"..."); library(x, character.only=T, verbose=F, warn.conflicts=F,quietly=T)} ))
sessionInfo()

####################
# 1.1: Read in counts
counts.fil<-read.delim("Phallii-FIL.counts") #Fil ASE - find these datasets with the accompanying MS
counts.hal<-read.delim("Phallii-HAL.counts") #Hal ASE
counts<-read.delim("Phallii_eQTL.counts") #raw counts total

####################
# 1.2: Find lines for analysis
colnames(counts.fil)[grep("FIL2_H328_83",colnames(counts.fil))]<-"FIL2_H328_383" #incorrectly labelled individual
colnames(counts.hal)[grep("FIL2_H328_83",colnames(counts.hal))]<-"FIL2_H328_383"
# there is some obvious contamination - remove these outliers
bad.lines<-vector()
for(i in c("HAL","FIL","F1","FH")){
  dat.hal<-counts.hal[,grep(i,colnames(counts.hal))]
  dat.fil<-counts.fil[,grep(i,colnames(counts.fil))]
  mean.hal<-apply(dat.hal,2,mean) ;  mean.fil<-apply(dat.fil,2,mean)
  plot(mean.hal,mean.fil, main=i)
  mod<-lm(mean.fil~mean.hal)
  ol<-outlierTest(mod)
  if(ol$bonf.p<0.0001){
    out<-names(ol$bonf);    bad.lines<-c(bad.lines, out)
  }
}
lines<-colnames(counts.fil)[-which(colnames(counts.fil)%in%bad.lines | grepl("FH", colnames(counts.fil)))]
counts.fil<-counts.fil[,lines]
counts.hal<-counts.hal[,lines]
parent.lines<-lines[!grepl("F1",lines) & lines %in% colnames(counts)]
countsF1<-counts[,lines[lines %in% colnames(counts)]]
counts<-counts[,parent.lines]

####################
# 1.3: Find genes for analysis
filmean<-apply(counts.fil[,grep("FIL", colnames(counts.fil))],1, mean)
halmean<-apply(counts.hal[,grep("HAL", colnames(counts.hal))],1, mean)
halopmean<-apply(counts.fil[,grep("HAL", colnames(counts.fil))],1, mean)
filopmean<-apply(counts.hal[,grep("FIL", colnames(counts.hal))],1, mean)
parentmeans<-apply(counts,1, mean)
asegenes<-Reduce(intersect, list(names(filmean)[filmean>5],
                                 names(halmean)[halmean>5],
                                 names(halmean)[halopmean/halmean<=0.01],
                                 names(filmean)[filopmean/filmean<=0.01]))
parentgenes<-names(parentmeans)[parentmeans>5]
length(asegenes)
counts.fil<-counts.fil[asegenes,grep("HAL",colnames(counts.fil), invert=T)]
counts.hal<-counts.hal[asegenes,grep("FIL",colnames(counts.hal), invert=T)]
counts<-counts[parentgenes,]
countsF1<-countsF1[parentgenes,]

####################
# 1.4: add in information data
info<-read.csv("JGI_Samples_3_6_14.csv") #experimental design information from Lowry/TEJ, via Scott
info$jgi.id<-paste(info$id,info$N,info$Plant.Number, sep="_")
#calculate the numeric time variable
info$time.sampled<-sapply(strsplit(as.character(info$Tissue_Time_July),":"), function(x) { x <- as.numeric(x); (x[1]-10)*60+(x[2]-58) })
info$timecat<-cut(info$time.sampled, c(-10,99.9,210), labels=c("early","late"))     
info<-info[,c("id","Treatment","jgi.id","VWC_July","time.sampled","timecat")]
info$generation<-ifelse(grepl("F1",info$jgi.id),"F1", 
                         ifelse(grepl("FIL",info$jgi.id) | grepl("HAL",info$jgi.id) ,"F0",
                                ifelse(grepl("FH",info$jgi.id),"F2","NA")))

####################
# 1.5: combine counts and info, make sure every looks good, then split

####################
# 1.5.1: get hal-specific info, merge with counts
c.hal<-data.frame(t(counts.hal))
c.hal$jgi.id<-rownames(c.hal)
c.hal$unique<-paste(c.hal$jgi.id,"hal",sep="_")
all.hal<-merge(c.hal, info,by="jgi.id")

####################
# 1.5.2: split out hal counts, add allele column
ch<-all.hal[,grep("Pahal",colnames(all.hal))]
infoh<-all.hal[,-grep("Pahal",colnames(all.hal))]
infoh$allele<-"hal"

####################
# 1.5.3: get fil-specific info, merge with counts
c.fil<-data.frame(t(counts.fil))
c.fil$jgi.id<-rownames(c.fil)
c.fil$unique<-paste(c.fil$jgi.id,"fil",sep="_")
all.fil<-merge(c.fil, info,by="jgi.id")

####################
# 1.5.4: split out fil counts, add allele column
cf<-all.fil[,grep("Pahal",colnames(all.fil))]
infof<-all.fil[,-grep("Pahal",colnames(all.fil))]
infof$allele<-"fil"

####################
# 1.5.5: merge two ase datasets.
identical(colnames(cf), colnames(ch))
ase.counts<-rbind(cf,ch)
ase.info<-rbind(infof,infoh)
ase.info$trtgenal<-paste(ase.info$Treatment, ase.info$generation, ase.info$allele, sep="_")
ase.info$generation<-factor(ase.info$generation, levels=c("F1","F0"))
ase.info$allele<-factor(ase.info$allele, levels=c("fil","hal"))
ase.info$Treatment<-factor(ase.info$Treatment, levels=c("Dry","Wet"))
ase.info$genal<-paste( ase.info$generation, ase.info$allele, sep="_")
ase.info$trtal<-paste( ase.info$Treatment, ase.info$allele, sep="_")
ase.info$trtgen<-paste( ase.info$Treatment, ase.info$generation, sep="_")
ase.info$genal<-factor(ase.info$genal,levels=c("F1_fil","F1_hal","F0_fil","F0_hal"))
ase.info$trans.fil<-with(ase.info,ifelse(grepl("fil",allele) | grepl("F1",generation), "trans.fil", "not"))
ase.info$trans.hal<-with(ase.info,ifelse(grepl("hal",allele) | grepl("F1",generation), "trans.hal", "not"))

ase.info.wet<-ase.info[ase.info$Treatment=="Wet",]
ase.info.dry<-ase.info[ase.info$Treatment=="Dry",]
rownames(ase.counts)<-ase.info$unique
ase.counts<-data.frame(t(ase.counts))
ase.counts.wet<-ase.counts[,ase.info.wet$unique]
ase.counts.dry<-ase.counts[,ase.info.dry$unique]

####################
# 1.5.6: prepare all dataset
c.all<-data.frame(t(counts))
c.all$jgi.id<-rownames(c.all)
c.all$unique<-c.all$jgi.id
all.all<-merge(c.all, info,by="jgi.id")
all.counts<-all.all[,grep("Pahal",colnames(all.all))]
all.info<-all.all[,-grep("Pahal",colnames(all.all))]
all.info$allele<-ifelse(grepl("HAL", all.info$jgi.id),"hal","fil")
all.info$trtgeno<-paste(all.info$Treatment, all.info$id, sep="_")
rownames(all.counts)<-all.info$unique
all.counts<-data.frame(t(all.counts))

####################
#1.5.7: prepare allf1 dataset
c.allF1<-data.frame(t(countsF1))
c.allF1$jgi.id<-rownames(c.allF1)
c.allF1$unique<-c.allF1$jgi.id
all.c.allF1<-merge(c.allF1, info,by="jgi.id")
allF1.counts<-all.c.allF1[,grep("Pahal",colnames(all.c.allF1))]
allF1.info<-all.c.allF1[,-grep("Pahal",colnames(all.c.allF1))]
allF1.info$allele<-ifelse(grepl("HAL", allF1.info$jgi.id),"hal",
                          ifelse(grepl("FIL", allF1.info$jgi.id),"fil","f1"))
allF1.info$trtgeno<-paste(allF1.info$Treatment, allF1.info$id, sep="_")
rownames(allF1.counts)<-allF1.info$unique
allF1.counts<-data.frame(t(allF1.counts))

####################
# 1.6: prepare annotation data
annot<-read.csv("ph2015_v1.1annot.edited.csv")
annot<-annot[,c("chr","start","strand","id", "end")]

####################
# 1.7: prepare DESeq2 dataset
all.annot<-annot[annot$id %in% rownames(all.counts),]
ase.annot<-annot[annot$id %in% rownames(ase.counts),]
all.gr<-GRanges(seqnames = all.annot$id, 
            ranges = IRanges(start=all.annot$start,  end=all.annot$end, names=all.annot$id),
            strand = all.annot$strand)
ase.gr<-GRanges(seqnames = ase.annot$id, 
                ranges = IRanges(start=ase.annot$start,  end=ase.annot$end, names=ase.annot$id),
                strand = ase.annot$strand)
all.descounts<-data.matrix(all.counts)
ase.descounts<-data.matrix(ase.counts)
ase.descounts.wet<-data.matrix(ase.counts.wet)
ase.descounts.dry<-data.matrix(ase.counts.dry)
allF1.descounts<-data.matrix(allF1.counts)
se.all<-SummarizedExperiment(assays = all.descounts, rowData = all.gr, colData = DataFrame(all.info))
se.ase<-SummarizedExperiment(assays = ase.descounts, rowData = ase.gr, colData = DataFrame(ase.info))
se.ase.wet<-SummarizedExperiment(assays = ase.descounts.wet, rowData = ase.gr, colData = DataFrame(ase.info.wet))
se.ase.dry<-SummarizedExperiment(assays = ase.descounts.dry, rowData = ase.gr, colData = DataFrame(ase.info.dry))
se.allF1<-SummarizedExperiment(assays = allF1.descounts, rowData = all.gr, colData = DataFrame(allF1.info))
dds.all1 <- DESeqDataSet(se = se.all, 
                    design = ~ timecat + Treatment + allele + Treatment*allele)
dds.all2 <- DESeqDataSet(se = se.all, 
                     design = ~ timecat + trtgeno)
dds.all3 <-DESeqDataSet(se = se.allF1, 
                        design = ~ timecat + Treatment + allele + Treatment*allele)
dds.ase1 <- DESeqDataSet(se = se.ase, design = ~ timecat + Treatment + allele + generation +
                           allele*generation + allele*Treatment + generation*Treatment + allele*generation*Treatment)
dds.ase2 <- DESeqDataSet(se = se.ase, design = ~ timecat + trtgenal)
dds.ase3 <- DESeqDataSet(se = se.ase, design = ~ timecat + Treatment + genal + Treatment*genal)
dds.ase4 <- DESeqDataSet(se = se.ase, design = ~ timecat + trtal + generation + trtal*generation)
dds.ase5 <- DESeqDataSet(se = se.ase, design = ~ timecat + allele + Treatment + allele*Treatment)
dds.ase6 <- DESeqDataSet(se = se.ase, design = ~ timecat + generation*allele + Treatment + generation*allele*Treatment)
dds.ase7 <- DESeqDataSet(se = se.ase, design = ~ timecat + trans.hal + trans.fil + Treatment + trans.hal*Treatment + trans.fil*Treatment)
dds.ase8 <- DESeqDataSet(se = se.ase, design = ~ timecat + allele)
dds.ase9 <- DESeqDataSet(se = se.ase, design = ~ timecat + generation*allele)
dds.ase10 <- DESeqDataSet(se = se.ase, design = ~ timecat + trans.hal + trans.fil)
dds.ase11 <- DESeqDataSet(se = se.ase, design = ~ timecat + allele + generation:allele)
dds.ase12 <- DESeqDataSet(se = se.ase, design = ~ timecat + allele + trans.hal + trans.fil)
dds.ase.wet<-DESeqDataSet(se = se.ase.wet, design = ~ timecat + allele + generation:allele)
dds.ase.dry<-DESeqDataSet(se = se.ase.dry, design = ~ timecat + allele + generation:allele)
save(dds.all1,dds.all2,dds.all3,dds.all3, dds.ase1,dds.ase2,dds.ase3,dds.ase4, dds.ase5, dds.ase6, dds.ase7,dds.ase8,dds.ase9,dds.ase10,dds.ase11,dds.ase12, all.annot, ase.annot, all.info, ase.info,allF1.info,dds.ase.wet,dds.ase.dry,ase.info.wet,ase.info.dry, file="ct2015final_deseqinput.RData")
