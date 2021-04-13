library(seqinr)
library(Biostrings)
library(seqRFLP)
a <- readDNAStringSet("~/virus/earlygene/dengue/ref-zika.fasta")[[1]] %>% as.character()
seqfile <- system("ls ~/virus/earlygene/dengue/mut*fasta",intern = T) 
inffile <- system("ls ~/virus/earlygene/dengue/inf*",intern = T)
codon <- read.table("~/fcs/codon.txt",header = TRUE,stringsAsFactors = F)

n=5

myseq <- readDNAStringSet(seqfile[n]) 
cdsinf <- read.csv(inffile[n],stringsAsFactors = F,sep = "\t",header = F)

load("~/codonpaper/code and data/1_2.xijofhost.finall.Rdata")
hum <- xijofhost %>% dplyr::filter(special =="Homo_sapiens")
nat <- xijofhost %>% dplyr::filter(special =="Aedes_aegypti")
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
    yij$nat <- nat$xij[match(yij$codon,nat$codon)]
    dd <- yij %>% dplyr::filter(!is.na(yij)) %>%
      group_by(aa,virusid,gene,type,len) %>%
      dplyr::summarize(Di=(sum((xij-yij)^2))^0.5,
                       Dinat=(sum((nat-yij)^2))^0.5) %>%
      as.data.frame() %>%
      group_by(virusid,gene,type,len) %>%
      dplyr::summarize(Dp = prod(Di)^(1/length(Di)),
                       Dpnat = prod(Dinat)^(1/length(Dinat)))
  } else {
    dd <- data.frame(virusid="999",gene="999",type="999",Dp=999,len=999,Dpnat=999)
  }
  dd
}) %>% rbind.fill() %>% dplyr::filter(Dp!=999)
seqfile[n]
save(Dp,file = "~/virus/earlygene/dengue/Dp.zika.Rdata")

## load data
load("~/virus/earlygene/dengue/Dp.dengue-1.Rdata")
va <- Dp %>% cbind(data.frame(virusname="Dengue 1"))
load("~/virus/earlygene/dengue/Dp.dengue-2.Rdata")
vb <- Dp %>% cbind(data.frame(virusname="Dengue 2"))
load("~/virus/earlygene/dengue/Dp.dengue-3.Rdata")
vc <- Dp %>% cbind(data.frame(virusname="Dengue 3"))
load("~/virus/earlygene/dengue/Dp.dengue-4.Rdata")
vd <- Dp %>% cbind(data.frame(virusname="Dengue 4"))
load("~/virus/earlygene/dengue/Dp.zika.Rdata")
ve <- Dp %>% cbind(data.frame(virusname="Zika"))

##picture and test
gene <- rbind(va,vb,vc,vd,ve) %>% dplyr::filter(type=="eachgene") 

type <- data.frame(gene=c("NS5","NS3","NS1","NS4B","NS2A","NS2B","NS4A","E","prM","ancC"),
                   type=c(rep("Non-structural",7),rep("Structural",3)))
gene$type2 <- type$type[match(gene$gene,as.vector(type$gene))]
source("~/Rfunction/style.print.R")
gene$virusname <- factor(gene$virusname,levels = c("Dengue 1","Dengue 2","Dengue 3","Dengue 4","Zika"))
gene %>% ggplot(aes(type2,Dp,color=type2)) + 
  geom_violin(scale = "width") +
  geom_boxplot(width = 0.1, outlier.colour = NA) +
  stat_summary(fun = 'mean', geom = 'point', shape = 18, colour = 'black')+
  scale_y_continuous(limits = c(0.1,0.85),breaks = c(0.25,0.50,0.75))+
  facet_wrap(~virusname,ncol = 3)+
  labs(x="Gene types of HAdV-D",y=expression(paste(italic(D)[P]," (virus/human)")))+
  style.print() +
  theme(legend.position = "None")

lapply(1:length(unique(gene$virusname)),function(x){
  v <- unique(gene$virusname)[x]
  a <- gene %>% dplyr::filter(virusname == v)
  s <- a %>% dplyr::filter(type2=="Structural")
  n <- a %>% dplyr::filter(type2=="Non-structural")
  data.frame(v,p=wilcox.test(s$Dp,n$Dp)$p.value,ms=mean(s$Dp),mn=mean(n$Dp))
}) %>% rbind.fill()




