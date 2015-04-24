#simulation of effects for cis trans test
#wts<-vector of weights indicating the effect size of the factors in each group
#nind<-number of individuals / group
#dis<-numeric indicator of dispersion around mean
#xbar<-intercept



sim.ct<-function(xbar=5, nind=80, sd.sim=.4, 
                 factors=list(allele=c("hal","fil"),generation=c("F0","F1"),treatment=c("wet","dry")),
                 model= "y~cis+trans.f+trans.h+trt+cis*trt + trans.f*trt + trans.h*trt",
                 wts=c(1,0,0,0,0,0,0), type="nbinom"){
  terms<-gsub(" ","",model); terms<-gsub("y~","",terms)
  terms<-unlist(strsplit(terms,"[+]"))
  model<-as.formula(model)                   
  grps<-expand.grid(factors)
  rownames(grps)<-paste(grps$allele, grps$generation, grps$treatment, sep=".")
  grps$effect<-0
  wts.fil.F0.wet<-c(1,1,0,1,1,1,0)
  wts.fil.F0.dry<-c(1,1,0,0,0,0,0)
  wts.hal.F0.wet<-c(0,0,1,1,0,0,1)
  wts.hal.F0.dry<-c(0,0,1,0,0,0,0)
  wts.fil.F1.wet<-c(1,1,1,1,1,1,1)
  wts.fil.F1.dry<-c(1,1,1,0,0,0,0)
  wts.hal.F1.wet<-c(0,1,1,1,0,1,1)
  wts.hal.F1.dry<-c(0,1,1,0,0,0,0)
  grps["fil.F0.wet","effect"]<-sum(wts*wts.fil.F0.wet)
  grps["hal.F0.wet","effect"]<-sum(wts*wts.hal.F0.wet)
  grps["fil.F1.wet","effect"]<-sum(wts*wts.fil.F1.wet)
  grps["hal.F1.wet","effect"]<-sum(wts*wts.hal.F1.wet)
  grps["fil.F0.dry","effect"]<-sum(wts*wts.fil.F0.dry)
  grps["hal.F0.dry","effect"]<-sum(wts*wts.hal.F0.dry)
  grps["fil.F1.dry","effect"]<-sum(wts*wts.fil.F1.dry)
  grps["hal.F1.dry","effect"]<-sum(wts*wts.hal.F1.dry)
  grps$means<-xbar+grps$effect
  ind.pergroup<-nind/dim(grps)[1]

  if(type=="normal"){
    for(i in 1:ind.pergroup){
      dat<-data.frame(rnorm(8,mean=grps$means,sd=sd.sim))
      colnames(dat)<-paste("ind",i,sep=".")
      grps<-cbind(grps,dat)
    }
  }else{
    for(i in 1:ind.pergroup){
      dat<-data.frame(rnegbin(8,mu=exp(grps$means),theta=exp(sd.sim)))
      colnames(dat)<-paste("ind",i,sep=".")
      grps<-cbind(grps,dat)
    }
  }
  library(reshape2)
  library(ggplot2)
  library(MASS)
  
  # Specify id.vars: the variables to keep but not split apart on
  grps.long<-melt(grps, id.vars=c("allele","generation","treatment","means","effect"))
  
  grps.long$cis<-grps.long$allele
  grps.long$trans.f<-ifelse(grps.long$generation == "F1" | grps.long$allele =="fil","trans.fil","not")
  grps.long$trans.h<-ifelse(grps.long$generation == "F1" | grps.long$allele =="hal","trans.hal","not")
  grps.long$treatment<-as.factor(grps.long$treatment)
  grps.long$trans.h<-as.factor(grps.long$trans.h)
  grps.long$trans.f<-as.factor(grps.long$trans.f)
  if(type=="normal"){
    glm.out<-lm(value~1+cis+trans.f+trans.h+cis*treatment+trans.f*treatment+trans.h*treatment, data=grps.long)
    print(summary(glm.out))
    print(ggplot(grps.long, aes(x=allele,y=value,color=generation))+
            geom_boxplot()+ facet_grid(treatment~.)+theme_bw()
    )
  }else{
    glm.out<-glm.nb(value~cis+trans.f+trans.h+cis*treatment+trans.f*treatment+trans.h*treatment, link="log", data=grps.long)
    print(summary(glm.out)$coefficients)
    print(ggplot(grps.long, aes(x=allele,y=value,color=generation))+
            geom_boxplot()+ facet_grid(treatment~.)+theme_bw()
    )
  }
  return(glm.out)
}
#cis effect only
test<-sim.ct(wts=c(2,0,0,0,0,0,0), nind=80, type="normal")
#cis and trans.f effect
test<-sim.ct(wts=c(2,2,0,0,0,0,0), nind=80, type="normal")
#cis and trans.h effect
test<-sim.ct(wts=c(2,0,2,0,0,0,0), nind=80, type="normal")
#cis and cis*trt int
test<-sim.ct(wts=c(1,0,0,0,2,0,0), nind=80, type="normal")
#trans.f only
test<-sim.ct(wts=c(0,2,0,0,0,0,0), nind=80, type="normal")
#trans.h only
test<-sim.ct(wts=c(0,0,2,0,0,0,0), nind=80, type="normal")
#both trans, no cis
test<-sim.ct(wts=c(0,2,2,0,0,0,0), nind=80, type="normal")
#trade off
test<-sim.ct(wts=c(-1,1,-1,0,0,0,0), nind=80, type="normal")
