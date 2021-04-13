library(seqinr)
library(Biostrings)
library(seqRFLP)
hadvfile <- system("ls ~/virus/earlygene/hadv/mut-hadv-*.fasta",intern = T) 
cdsinffile <- system("ls ~/virus/earlygene/hadv/infhadv*",intern = T)
load("~/Rfunction/wiforite.Rdata")
load("~/codonpaper/review/codon.Rdata")

n=6
hadv <- readDNAStringSet(hadvfile[n]) 
cdsinf <- read.csv(cdsinffile[n],stringsAsFactors = F,sep = "\t",header = F)
load("~/codonpaper/code and data/1_2.xijofhost.finall.Rdata")
hum <- xijofhost %>% dplyr::filter(special =="Homo_sapiens")

##calculate Dp and save
Dp <- mclapply(mc.cores = 20,1:length(hadv),function(y){
  a <- hadv[y]
  query <- strsplit(names(a)," ")[[1]][1]
  inftmp <- cdsinf %>% dplyr::filter(V1==query)
  if(nrow(inftmp)>0){
    yfre <- lapply(1:nrow(inftmp),function(z){
      ref1 <- inftmp[z,9]
      ref2 <- inftmp[z,10]
      quer1 <- inftmp[z,7]
      quer2 <- inftmp[z,8]
      gene <- strsplit(inftmp[z,2],",")[[1]][2]
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
    ##Dp
    Dp <- yij %>% dplyr::filter(!is.na(yij)) %>%
      group_by(aa,virusid,gene,type,len) %>%
      dplyr::summarize(Di=(sum((xij-yij)^2))^0.5) %>%
      as.data.frame() %>%
      group_by(virusid,gene,type,len) %>%
      dplyr::summarize(Dp = prod(Di)^(1/length(Di)))
    
   
    Dp$ite <- ite$ite[match(as.vector(Dp$gene),as.vector(ite$gene))]
    
  } else {
    Dp <- data.frame(virusid="999",gene="999",type="999",Dp=999,len=999)
  }
  Dp
}) %>% rbind.fill() %>% dplyr::filter(Dp!=999)
hadvfile[n]
save(Dp,file = "~/virus/earlygene/hadv/Dp.hadvF.Rdata")

## load data
load("~/virus/earlygene/hadv/Dp.hadvA.Rdata")
va <- Dp %>% cbind(data.frame(virusname="HAdV-A"))
load("~/virus/earlygene/hadv/Dp.hadvB.Rdata")
vb <- Dp %>% cbind(data.frame(virusname="HAdV-B"))
load("~/virus/earlygene/hadv/Dp.hadvC.Rdata")
vc <- Dp %>% cbind(data.frame(virusname="HAdV-C"))
load("~/virus/earlygene/hadv/Dp.hadvD.Rdata")
vd <- Dp %>% cbind(data.frame(virusname="HAdV-D"))
load("~/virus/earlygene/hadv/Dp.hadvE.Rdata")
ve <- Dp %>% cbind(data.frame(virusname="HAdV-E"))
load("~/virus/earlygene/hadv/Dp.hadvF.Rdata")
vf <- Dp %>% cbind(data.frame(virusname="HAdV-F"))

##merge data
genehadv <- rbind(va,vb,vc,vd,ve,vf) %>% dplyr::filter(type=="eachgene") 
load("~/virus/earlygene/hadv/hadvgene.Rdata")
genehadv$type2 <- hadvgene$type[match(as.vector(genehadv$gene),hadvgene$gene)]

##picture and test
##1)violin
source("~/Rfunction/style.print.R")
genehadv1 <- genehadv %>% dplyr::filter(virusname=="HAdV-B")
genehadv1 %>% ggplot(aes(type2,Dp)) + 
  geom_violin(scale = "width") +
  geom_boxplot(width = 0.1, outlier.colour = NA) +
  stat_summary(fun = 'mean', geom = 'point', shape = 18, colour = 'black')+
  #scale_y_continuous(limits = c(0.1,1),breaks = c(0.25,0.50,0.75))+
  labs(x="Gene types of HAdV-A",y=expression(paste(italic(I)[TE]," (virus/human)")))+
  style.print()
E <- genehadv1 %>% dplyr::filter(type2=="E")
I <- genehadv1 %>% dplyr::filter(type2=="I")
L <- genehadv1 %>% dplyr::filter(type2=="L")
wilcox.test(E$Dp,I$Dp)$p.value
wilcox.test(E$Dp,L$Dp)$p.value
wilcox.test(I$Dp,L$Dp)$p.value

##2)cumulative distribution
genehadv %>%
  dplyr::filter(type2!="I") %>%
  ggplot(aes(Dp,color=type2))+ 
  stat_ecdf()+
  facet_wrap(~virusname)+
  #scale_y_continuous(limits = c(0.51,1),breaks = c(0.6,0.8,1))+
  labs(x=expression(paste(italic(D)[P]," (virus/human)")),y= "Frequency")+
  style.print()

