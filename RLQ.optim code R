#........................................................................................................#
#................Quantifying succession of phytoplankton trait-environment associations..................#
#........................................................................................................#

#........................................................................................................#
#...........Author of the script: Charlie Loewen.........................................................#
#...........Date latest modifications: 2019-10-28........................................................#
#........................................................................................................#

#########################################################################
################# For finding traits and environmental variables optimizing RLQ  model
#########################################################################

##################
# Clear any variables from the R environment
##################

rm(list=ls())
while(dev.cur() != 1){ dev.off(dev.cur()) }

##################
# Load R packages
##################

library(ade4)

##################
# Load data in R environment
##################

species=read.table("Phyto_data-site_sp.csv",header=T,sep=",",row.names = 1) # matrix of sites (rows) by species (columns)
traitopts=read.table("Phyto_data-sp_trait.csv",header=T,sep=",",row.names = 1) # matrix of species (rows) by traits (columns)
envopts=read.table("Phyto_data-site_env.csv",header=T,sep=",",row.names = 1) # matrix of sites (rows) by environmental variables (columns)

##################
# RLQ optimization 
##################

## PREP ##

# Create data frame of all possible trait and environmental variable combinations
combolist<-list(a = colnames(traitopts), b = colnames(traitopts), c = colnames(envopts), d = colnames(envopts))
allcombo<-expand.grid(combolist)
allcombo<-allcombo[!allcombo$a == allcombo$b, ]
allcombo<-allcombo[!allcombo$c == allcombo$d, ]
indx<-!duplicated(t(apply(allcombo, 1, sort))) # finds non - duplicates in sorted rows
allcombo<-allcombo[indx, ] # selects only the non - duplicates according to that index

## SEEDING ##

output<-as.data.frame(matrix(ncol=7))

RLQ.optim = function(i,j,k,l){
  TRAIT.seed<-cbind(traitopts[i],traitopts[j])
  ENV.seed<-cbind(envopts[k],envopts[l])
  Lsp<-dudi.coa(species, scannf = FALSE)
  Renv<-dudi.hillsmith(ENV.seed, row.w = Lsp$lw, scannf = FALSE)
  Qtrait<-dudi.hillsmith(TRAIT.seed, row.w = Lsp$cw, scannf = FALSE)
  rlq<-rlq(Renv, Lsp, Qtrait, scannf = FALSE)
  testrlq<-randtest(rlq, modeltype = 6)
  output<-c(i,j,k,l,summary(rlq)$corr,sum(summary(rlq)$corr),use.names=TRUE)
}

Seeding<-mapply(RLQ.optim,i=as.character(allcombo$a), j=as.character(allcombo$b), k=as.character(allcombo$c), l=as.character(allcombo$d))

transposed.list.Seeding<-do.call(rbind.data.frame, Seeding)
colnames(transposed.list.Seeding)<-c("i","j","k","l","corr1","corr2","corrSUM")

## STEPWISE SELECTION ##

# Confirm significance of initial seeding round

TRAIT.start<-cbind(traitopts["TRAIT.1"], # TRAIT.1 = top ranking trait variable from seeding round
                   traitopts["TRAIT.2"]) # TRAIT.2 = second ranking trait variable from seeding round
ENV.start<-cbind(envopts["ENV.1"], # ENV.1 = top ranking environmental variable from seeding round
                 envopts["ENV.2"]) # ENV.2 = second ranking environmental variable from seeding round
Lsp<-dudi.coa(species, scannf = FALSE)
Renv<-dudi.hillsmith(ENV.start, row.w = Lsp$lw, scannf = FALSE)
Qtrait<-dudi.hillsmith(TRAIT.start, row.w = Lsp$cw, scannf = FALSE)
rlq<-rlq(Renv, Lsp, Qtrait, scannf = FALSE)
randtest(rlq, modeltype = 6)

# Step1

output<-as.data.frame(matrix(ncol=6))

step.functi = function(i){
  TRAIT.stepi<-cbind(traitopts["TRAIT.1"], # TRAIT.1 = top ranking trait variable from seeding round
                      traitopts["TRAIT.2"], # TRAIT.2 = second ranking trait variable from seeding round
                      traitopts[i])
  ENV.stepi<-cbind(envopts["ENV.1"], # ENV.1 = top ranking environmental variable from seeding round
                    envopts["ENV.2"]) # ENV.2 = second ranking environmental variable from seeding round
  Lsp<-dudi.coa(species, scannf = FALSE)
  Renv<-dudi.hillsmith(ENV.stepi, row.w = Lsp$lw, scannf = FALSE)
  Qtrait<-dudi.hillsmith(TRAIT.stepi, row.w = Lsp$cw, scannf = FALSE)
  rlq<-rlq(Renv, Lsp, Qtrait, scannf = FALSE)
  testrlq<-randtest(rlq, modeltype = 6)
  output<-c(i,summary(rlq)$corr,sum(summary(rlq)$corr),testrlq$pvalue, use.names=TRUE)
}

step1a<-sapply(as.character(colnames(traitopts)), FUN=step.functi)
t.step1a<-t(step1a)
colnames(t.step1a)<-c("variable","corr1","corr2","corrSum","pval1","pval2")

step.functj = function(j){
  TRAIT.stepj<-cbind(traitopts["TRAIT.1"], # TRAIT.1 = top ranking trait variable from seeding round
                     traitopts["TRAIT.2"]) # TRAIT.2 = second ranking trait variable from seeding round
  ENV.stepj<-cbind(envopts["ENV.1"], # ENV.1 = top ranking environmental variable from seeding round
                   envopts["ENV.2"], # ENV.2 = second ranking environmental variable from seeding round
                   envopts[j])
  Lsp<-dudi.coa(species, scannf = FALSE)
  Renv<-dudi.hillsmith(ENV.stepj, row.w = Lsp$lw, scannf = FALSE)
  Qtrait<-dudi.hillsmith(TRAIT.stepj, row.w = Lsp$cw, scannf = FALSE)
  rlq<-rlq(Renv, Lsp, Qtrait, scannf = FALSE)
  testrlq<-randtest(rlq, modeltype = 6)
  output<-c(j,summary(rlq)$corr,sum(summary(rlq)$corr),testrlq$pvalue, use.names=TRUE)
}

step1b<-sapply(as.character(colnames(envopts)), FUN=step.functj)
t.step1b<-t(step1b)
colnames(t.step1b)<-c("variable","corr1","corr2","corrSum","pval1","pval2")

step1merge <- rbind(t.step1a,t.step1b
                    [, colnames(t.step1a)])

### Revise step.funci and step.funcj to include trait or environmental variable with highest corrSum value, and repeat above procedure (i.e. TRAIT.3 or ENV.3). Continue step selection until Pval <0.05 or corrSum <0.005 to optimize RLQ and find dominant trait-environment relationships.

#........................................................................................................#
