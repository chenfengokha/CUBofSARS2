library(seqinr)
library(Biostrings)
library(seqRFLP)


seqfile <- system("ls ~/virus/earlygene/coronavirus/*.fasta",intern = T) 
inffile <- system("ls ~/virus/earlygene/coronavirus/inf*",intern = T)
codon <- read.table("~/fcs/codon.txt",header = TRUE,stringsAsFactors = F)

n=6

myseq <- readDNAStringSet(seqfile[n]) 
cdsinf <- read.csv(inffile[n],stringsAsFactors = F,sep = "\t",header = F)

load("~/codonpaper/code and data/1_2.xijofhost.finall.Rdata")
hum <- xijofhost %>% dplyr::filter(special =="Homo_sapiens")
Dp <- mclapply(mc.cores = 20,1:length(myseq),function(y){
  a <- myseq[y]
  query <- strsplit(names(a)," ")[[1]][1]
  inftmp <- cdsinf %>% dplyr::filter(V1==query)
  if(nrow(inftmp)>0){
    yfre <- lapply(1:nrow(inftmp),function(z){
      ref1 <- inftmp[z,9]
      ref2 <- inftmp[z,10]
      quer1 <- inftmp[z,7]
      quer2 <- inftmp[z,8]
      gene <- strsplit(inftmp[z,2],",")[[1]]
      virusid <- inftmp[z,]$V1
      if(ref1<ref2){
        if(ref1%%3==1){
          quer1 <- quer1
        } 
        if(ref1%%3==2){
          quer1 <- quer1+2
        } 
        if(ref1%%3==0){
          quer1 <- quer1+1
        }
        if(ref2%%3==1){
          quer2 <- quer2-1
        } 
        if(ref2%%3==2){
          quer2 <- quer2-2
        } 
        if(ref2%%3==0){
          quer2 <- quer2
        }
        myseq <- as.character(substr(as.character(a),quer1,quer2))
      }
      
      if(ref1>ref2){
        if(ref1%%3==1){
          quer2 <- quer2
        } 
        if(ref1%%3==2){
          quer2 <- quer2-2
        } 
        if(ref1%%3==0){
          quer2 <- quer2-1
        }
        if(ref2%%3==1){
          quer1 <- quer1+1
        } 
        if(ref2%%3==2){
          quer1 <- quer1+2
        } 
        if(ref2%%3==0){
          quer1 <- quer1
        }
        myseq <- as.character(substr(as.character(a),quer1,quer2)) %>% revComp()
      }
      a2 <- substring(as.character(myseq),seq(1,(nchar(myseq)-2),by=3),seq(3,nchar(myseq),by=3)) %>% table() %>% as.data.frame()
      names(a2)[1] <- "V1"
      a3 <- codon %>% dplyr::filter(aa!="W")
      a3$fre <- a2$Freq[match(as.vector(a3$codon),as.vector(a2$V1))]
      a3$fre[is.na(a3$fre)] <- 0
      a3 %>% cbind(data.frame(gene,virusid))
    })%>% rbind.fill() 
    yijall <- yfre %>% group_by(aa,codon,virusid) %>%
      dplyr::summarize(fre1=sum(fre)) %>%
      group_by(aa,virusid) %>%
      dplyr::mutate(fre2=sum(fre1)) %>%
      group_by(aa,codon,virusid) %>%
      dplyr::summarize(yij=fre1/fre2)
    yijgene <- yfre %>% group_by(aa,codon,gene,virusid) %>%
      dplyr::summarize(fre1=sum(fre)) %>%
      group_by(aa,gene,virusid) %>%
      dplyr::mutate(fre2=sum(fre1)) %>%
      group_by(aa,codon,gene,virusid) %>%
      dplyr::summarize(yij=fre1/fre2)
    yij <- rbind(cbind(yijall %>% as.data.frame(),data.frame(gene=unique(yijall$virusid),type="all",len=3*sum(yfre$fre))),
                 cbind(yijgene %>% as.data.frame(),data.frame(type="eachgene",len=1)))
    yij$xij <- hum$xij[match(yij$codon,hum$codon)]
    dd <- yij %>% dplyr::filter(!is.na(yij)) %>%
      group_by(aa,virusid,gene,type,len) %>%
      dplyr::summarize(Di=(sum((xij-yij)^2))^0.5) %>%
      as.data.frame() %>%
      group_by(virusid,gene,type,len) %>%
      dplyr::summarize(Dp = prod(Di)^(1/length(Di)))
  } else {
    dd <- data.frame(virusid="999",gene="999",type="999",Dp=999,len=999)
  }
  dd
}) %>% rbind.fill() %>% dplyr::filter(Dp!=999)
seqfile[n]
save(Dp,file = "~/virus/earlygene/coronavirus/Dp.SARS.Rdata")

## all analysis
load("~/virus/earlygene/coronavirus/Dp.229E.Rdata")
va <- Dp %>% cbind(data.frame(virusname="229E"))
aa <- va %>% dplyr::filter(type=="all") %>% arrange(desc(len))
aa <- aa %>% dplyr::filter(len>25000)
load("/mnt/data/home/chenfeng/virus/earlygene/coronavirus/Dp.HKU1.Rdata")
vb <- Dp %>% cbind(data.frame(virusname="HKU1"))
bb <- vb %>% dplyr::filter(type=="all") %>% arrange(desc(len))
bb <- bb %>% dplyr::filter(len>25000)
load("/mnt/data/home/chenfeng/virus/earlygene/coronavirus/Dp.MERS.Rdata")
vc <- Dp %>% cbind(data.frame(virusname="MERS"))
cc <- vc %>% dplyr::filter(type=="all") %>% arrange(desc(len))
cc <- cc %>% dplyr::filter(len>25000)
load("/mnt/data/home/chenfeng/virus/earlygene/coronavirus/Dp.NL63.Rdata")
vd <- Dp %>% cbind(data.frame(virusname="NL63"))
dd <- vd %>% dplyr::filter(type=="all") %>% arrange(desc(len))
dd <- dd %>% dplyr::filter(len>25000)
load("/mnt/data/home/chenfeng/virus/earlygene/coronavirus/Dp.OC43.Rdata")
ve <- Dp %>% cbind(data.frame(virusname="OC43"))
ee <- ve %>% dplyr::filter(type=="all") %>% arrange(desc(len))
ee <- ee %>% dplyr::filter(len>25000)
load("/mnt/data/home/chenfeng/virus/earlygene/coronavirus/Dp.SARS2.Rdata")
vf <- Dp %>% cbind(data.frame(virusname="SARS2"))
ff <- vf %>% dplyr::filter(type=="all") %>% arrange(desc(len))
ff <- ff %>% dplyr::filter(len>25000)
load("/mnt/data/home/chenfeng/virus/earlygene/coronavirus/Dp.SARS.Rdata")
vg <- Dp %>% cbind(data.frame(virusname="SARS"))
gg <- vg %>% dplyr::filter(type=="all") %>% arrange(desc(len))
gg <- gg %>% dplyr::filter(len>25000)
##ful length Dp ~ time
library(timeDate)
ful <- rbind(aa,bb,cc,dd,ee,ff,gg)
ful$virusid <- as.vector(ful$virusid)
ful$id <- sub('..$','',ful$virusid)
geneA <- read.csv("~/virus/earlygene/coronavirus/mut-229E.csv",stringsAsFactors = F,sep = ",",header = T)
geneB <- read.csv("~/virus/earlygene/coronavirus/mut-hku1.csv",stringsAsFactors = F,sep = ",",header = T)
geneC <- read.csv("~/virus/earlygene/coronavirus/mut-MERS.csv",stringsAsFactors = F,sep = ",",header = T)
geneD <- read.csv("~/virus/earlygene/coronavirus/mut-nl63.csv",stringsAsFactors = F,sep = ",",header = T)
geneE <- read.csv("~/virus/earlygene/coronavirus/mut-oc43.csv",stringsAsFactors = F,sep = ",",header = T)
geneF <- read.csv("~/virus/earlygene/coronavirus/mut-SARS1.csv",stringsAsFactors = F,sep = ",",header = T)
geneG <- read.csv("~/virus/earlygene/coronavirus/mut-SARS2.csv",stringsAsFactors = F,sep = ",",header = T)
geneall <- rbind(geneA,geneB,geneC,geneD,geneE,geneF,geneG)
ful$Date <- geneall$Release_Date[match(ful$id,geneall$Accession)]
source("~/Rfunction/style.print.R")
ful$time <- as.double(as.Date(ful$Date))
ful$virusname <- factor(ful$virusname,levels = c("HKU1","NL63","OC43","229E","SARS2","MERS","SARS"))
ful %>% 
  ggplot(aes(x=virusname,y=Dp))+
  geom_violin(scale = "width")+
  geom_boxplot(width = 0.1, outlier.colour = NA)+
  stat_summary(fun = 'mean', geom = 'point', shape = 18, colour = 'black')+
  labs(x="",y=expression(paste(italic(D)[P]," (virus/human)")))+
  scale_y_continuous(limits = c(0.2,0.7),breaks = c(0.2,0.4,0.6))+
  style.print()
pairwise.wilcox.test(ful$Dp, ful$virusname)
##group by symptoms
ful$type2 <-"nonill"
ful$type2[which(ful$virusname %in% c("SARS2","MERS","SARS"))] <-"ill"
ful$type2 <- factor(ful$type2,levels = c("nonill","ill"))
ful %>% 
  ggplot(aes(x=type2,y=Dp))+
  geom_violin(scale = "width")+
  geom_boxplot(width = 0.1, outlier.colour = NA)+
  stat_summary(fun.y = 'mean', geom = 'point', shape = 18, colour = 'black')+
  labs(x="",y=expression(paste(italic(D)[P]," (virus/human)")))+
  scale_y_continuous(limits = c(0.2,0.7),breaks = c(0.2,0.4,0.6))+
  style.print()
pairwise.wilcox.test(ful$Dp, ful$type2)$p.value

## each gene
gene <- rbind(va,vb,vc,vd,ve,vf,vg) %>% dplyr::filter(type=="eachgene") 
gene <- gene %>% separate("gene",c("gene","type2"),"_")
source("~/Rfunction/style.print.R")
gene %>% 
  ggplot(aes(x=type2,y=Dp,color=type2))+
  geom_violin(scale = "width") +
  geom_boxplot(width = 0.1, outlier.colour = NA) +
  stat_summary(fun = 'mean', geom = 'point', shape = 18, colour = 'black')+
  facet_wrap(~virusname,ncol = 4)+
  labs(x=expression(paste(italic(D)[P]," (virus/human)")),y="Frequency")+
  style.print()+theme(legend.position = "none")
gene$virusname <- as.vector(gene$virusname)
lapply(1:length(unique(gene$virusname)),function(x){
  v <- unique(gene$virusname)[x]
  a <- gene %>% dplyr::filter(virusname == v)
  s <- a %>% dplyr::filter(type2=="Structure")
  n <- a %>% dplyr::filter(type2=="Non-structure")
  data.frame(v,p=wilcox.test(s$Dp,n$Dp)$p.value,ms=mean(s$Dp),mn=mean(n$Dp))
}) %>% rbind.fill()



#each gene ~ time
geneden <- rbind(va,vb,vc,vd,ve,vf,vg) %>% dplyr::filter(type=="eachgene")
geneden$virusid <- as.vector(geneden$virusid)
geneden$id <- sub('..$','',geneden$virusid)
geneA <- read.csv("~/virus/earlygene/coronavirus/mut-229E.csv",stringsAsFactors = F,sep = ",",header = T)
geneB <- read.csv("~/virus/earlygene/coronavirus/mut-hku1.csv",stringsAsFactors = F,sep = ",",header = T)
geneC <- read.csv("~/virus/earlygene/coronavirus/mut-MERS.csv",stringsAsFactors = F,sep = ",",header = T)
geneD <- read.csv("~/virus/earlygene/coronavirus/mut-nl63.csv",stringsAsFactors = F,sep = ",",header = T)
geneE <- read.csv("~/virus/earlygene/coronavirus/mut-oc43.csv",stringsAsFactors = F,sep = ",",header = T)
geneF <- read.csv("~/virus/earlygene/coronavirus/mut-SARS1.csv",stringsAsFactors = F,sep = ",",header = T)
geneG <- read.csv("~/virus/earlygene/coronavirus/mut-SARS2.csv",stringsAsFactors = F,sep = ",",header = T)
geneall <- rbind(geneA,geneB,geneC,geneD,geneE,geneF,geneG)
geneden$Date <- geneall$Release_Date[match(geneden$id,geneall$Accession)]




geneden <- geneden %>% separate("gene",c("gene","type2"),"_")

source("~/Rfunction/style.print.R")
geneden$time <- as.double(as.Date(geneden$Date)) 

geneden$type3 <-"Nonsymptom"
geneden$type3[which(geneden$virusname %in% c("SARS2","MERS","SARS"))] <-"Severe symptom"
geneden %>%
  group_by(type3,type2,time) %>%
  dplyr::summarize(me=mean(Dp)) %>%
  ggplot(aes(time,me))+
  geom_point(color="grey")+
  geom_smooth(method = "lm",se=F)+
  facet_grid(type3~type2,scale="free")+
  #scale_y_continuous(limits = c(0.3,0.5),breaks = c(0.3,0.4,0.5))+
  labs(x="Days",y=expression(paste(italic(D)[P]," (virus/human)")))+
  style.print()

geneden %>% 
  group_by(type2,type3,time) %>%
  dplyr::summarize(me=mean(Dp)) %>%
  group_by(type2,type3) %>%
  dplyr::summarize(p=cor.test(time,me,method = "s")$p.value,
                   rho=cor.test(time,me,method = "s")$estimate)

geneden %>% 
  dplyr::filter(virusname=="SARS2") %>%
  group_by(virusname,type2,time) %>%
  dplyr::summarize(me=mean(Dp)) %>%
  ggplot(aes(time,me))+
  geom_point(color="grey")+
  geom_smooth(method = "lm",se=F)+
  facet_grid(virusname~type2,scale="free")+
  scale_y_continuous(limits = c(0.3,0.5),breaks = c(0.3,0.4,0.5))+
  labs(x="Days",y=expression(paste(italic(D)[P]," (virus/human)")))+
  style.print()

geneden %>% 
  dplyr::filter(virusname=="SARS2") %>%
  group_by(type2,virusname,time) %>%
  dplyr::summarize(me=mean(Dp)) %>%
  group_by(type2,virusname) %>%
  dplyr::summarize(p=cor.test(time,me,method = "s")$p.value,
                   rho=cor.test(time,me,method = "s")$estimate)
