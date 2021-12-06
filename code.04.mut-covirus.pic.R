library(seqinr)
library(Biostrings)
library(seqRFLP)
library(ggtree)
library(castor)
library(stringr)
library(ape)
library(tidyr)
fun.id <- function(x){
  if(substr(x,1,5)=="join("){
    id <- strsplit(strsplit(x,"[()]")[[1]][2],"[.]")[[1]][1]
  } else(
    id <- strsplit(x,"[.]")[[1]][1]
  )
  id
}

load("~/virus/earlygene/coronavirus/Dp.all.SARS2.Rdata")

va <- Dp %>% cbind(data.frame(virusname="SARS2"))
aa <- va %>% dplyr::filter(type=="all") %>% arrange(desc(len))
aa <- aa %>% dplyr::filter(len>18000)
aa$id <- unlist(mclapply(mc.cores = 10,aa$virusid,fun.id))
meta <- read.csv("~/virus/earlygene/coronavirus/nextstrain.sars2.clade.tsv",stringsAsFactors = F,sep = "\t",header = T)
aa$clade <- meta$Nextstrain_clade[match(aa$id,meta$genbank_accession)]
aa$date <- meta$date[match(aa$id,meta$genbank_accession)]
aa$strain <- meta$strain[match(aa$id,meta$genbank_accession)]
source("~/Rfunction/style.print.R")
aa <- aa %>% dplyr::filter(nchar(date) == nchar("2021-05-31"))
aa$time <- as.double(as.Date(aa$date))
aa <- aa %>% dplyr::filter(!is.na(time))
cor.test(aa$Dp,aa$time,method = "p")
summary(lm(aa$Dp~aa$time)) -> t
t$r.squared

aa %>%
  ggplot(aes(time,Dp))+
  geom_point(size=0.03,alpha=0.03)+
  geom_smooth(method = "lm",se=F)+
  labs(x="Collection time",y=expression(paste(italic(D)[P]," (virus/human)")))+
  style.print()

aa %>% group_by(clade) %>%
  dplyr::summarize(rho=cor.test(time,Dp,method = "p")$estimate[[1]],
                   p=cor.test(time,Dp,method = "p")$p.value,
                   n=length(clade)) %>% arrange(p) %>% as.data.frame()

##PIC
strfun <- function(i){
  return(strsplit(i,"[|]")[[1]][1])
}

mytree <- read.tree("~/virus/earlygene/coronavirus/tree.SARS2.2021-10-04.all.nwk")
mytree$tip.label <- unlist(lapply(mytree$tip.label,strfun))
aa %>% dplyr::filter(strain %in% mytree$tip.label) -> ful
nonlable <- unique(c(mytree$tip.label[which(!(mytree$tip.label %in% ful$strain))] %>% unique(),mytree$tip.label[which(duplicated(mytree$tip.label))]))
mytree <- drop.tip(mytree,mytree$tip.label[match(nonlable,mytree$tip.label)])
mytree$edge.length <- mytree$edge.length + 1e-07
data.frame(stringsAsFactors = F,a=mytree$tip.label,b=1:length(mytree$tip.label)) -> ranktmp
ful$rank <- ranktmp$b[match(ful$strain,ranktmp$a)]
ful %>% dplyr::filter(!is.na(rank)) %>% arrange(rank) -> dd.pic
dd.pic[,c(10,1,11)] -> ff
rownames(ff) <- ff$strain
ff[,c(-1)] -> ff
ff <- as.matrix(ff)
get_independent_contrasts(mytree,ff)$PICs %>% as.data.frame() -> ff
tt <- summary(lm(Dp~time-1,ff))
ff %>%
  ggplot(aes(time,Dp))+
  geom_point(size=0.3,alpha=0.3)+
  geom_smooth(method = "lm",se=F)+
  labs(x="PIC of time",y=expression(paste("PIC of ",italic(D)[P]," (virus/human)")))+
  geom_hline(yintercept = 0, linetype=2,color="red")+
  geom_vline(xintercept = 0, linetype=2,color="red")+
  style.print()

##test for other  other virus
Dpfile <- system("ls ~/virus/earlygene/coronavirus/Dp.*",intern = T)
treefile <- system("ls ~/virus/earlygene/coronavirus/tree.*",intern = T)
inffile <- system("ls ~/virus/earlygene/coronavirus/mut*csv",intern = T)
#1,3,4,5,6,7
n=7
Dpfile[n]
load(Dpfile[n])

ff <- Dp %>% dplyr::filter(type=="all") %>% arrange(desc(len))
##pic (ful length Dp ~ time)
library(timeDate)
ful <- ff
ful$virusid <- as.vector(ful$virusid)
ful$id <- sub('..$','',ful$virusid)
geneall <- read.csv(inffile[n],stringsAsFactors = F,sep = ",",header = T)
str_replace(string = ful$virusid,pattern = "NC_",replacement = "NC") -> ful$virusid
mytree <- read.tree(treefile[n])
str_replace(string = mytree$tip.label,pattern = "NC_",replacement = "NC") -> mytree$tip.label
mytree$tip.label <- substr(mytree$tip.label,1,10)
ful %>% dplyr::filter(virusid %in% mytree$tip.label) -> ful
nonlable <- mytree$tip.label[which(!(mytree$tip.label %in% ful$virusid))] %>% unique()
mytree <- drop.tip(mytree,mytree$tip.label[match(nonlable,mytree$tip.label)])
mytree$edge.length <- mytree$edge.length + 1e-07
data.frame(stringsAsFactors = F,a=mytree$tip.label,b=1:length(mytree$tip.label)) -> ranktmp
ful$rank <- ranktmp$b[match(ful$virusid,ranktmp$a)]
ful %>% dplyr::filter(!is.na(rank)) %>% arrange(rank) -> dd
dd[,c(2,1,8)] -> dd
rownames(dd) <- dd$virusid
dd[,c(-1)] -> dd
dd <- as.matrix(dd)
nrow(dd)
Dpfile[n]
get_independent_contrasts(mytree,dd)$PICs %>% as.data.frame() -> ff
tt <- lm(Dp~reltime-1,ff)
summary(tt)
