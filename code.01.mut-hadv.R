library(seqinr)
library(Biostrings)
library(seqRFLP)
##nul
from=c("A","T","G","C","a","g","t","c","N","n")
to  =c("T","A","C","G","t","c","a","g","N","n")
names(to)=from
revcom <- function(x){
  sep_DNA <- unlist(strsplit(x,""))
  complementary_DNA <- to[sep_DNA]
  rev_complementary <- rev(complementary_DNA)
  rev_complementary_DNA <- paste(rev_complementary,collapse = "")
  rev_complementary_DNA
}


hadvfile <- system("ls ~/virus/earlygene/hadv/mut-hadv-*.fasta",intern = T) 
cdsinffile <- system("ls ~/virus/earlygene/hadv/infhadv*",intern = T)
codon <- read.table("/mnt/data4/disk/chenfeng/codonpaper/expriment/fcs/codon.txt",header = TRUE,stringsAsFactors = F)

n=7
hadv <- readDNAStringSet(hadvfile[n]) 
cdsinf <- read.csv(cdsinffile[n],stringsAsFactors = F,sep = "\t",header = F)
load("~/codonpaper/code and data/1_2.xijofhost.finall.Rdata")
hum <- xijofhost %>% dplyr::filter(special =="Homo_sapiens")

##calculate Dp and save
Dp <- mclapply(mc.cores = 60,1:length(hadv),function(y){
  a <- hadv[y]
  query <- strsplit(names(a)," ")[[1]][1]
  inftmp <- cdsinf %>% dplyr::filter(V1==query)
  if(nrow(inftmp)>0){
    yfre <- lapply(1:nrow(inftmp),function(z){
      
      #gene <- strsplit(inftmp[z,2],",")[[1]][2]
      gene <- inftmp[z,2]
      virusid <- inftmp[z,]$V1
      myseq <- as.character(substr(as.character(a),inftmp[z,7],inftmp[z,8]))
      #
      if(inftmp[z,9] < inftmp[z,10]){
        a1 <- which(s2c(as.character(translate(DNAString(myseq), no.init.codon=TRUE,if.fuzzy.codon="X"))) == "*")
        a2 <- which(s2c(as.character(translate(DNAString(substr(myseq,2,nchar(myseq))), no.init.codon=TRUE,if.fuzzy.codon="X"))) == "*")
        a3 <- which(s2c(as.character(translate(DNAString(substr(myseq,3,nchar(myseq))), no.init.codon=TRUE,if.fuzzy.codon="X"))) == "*")
        
        l1 <- length(which(a1 < (floor(nchar(myseq)/3)-2)))
        l2 <- length(which(a2 < (floor(nchar(myseq)/3)-2)))
        l3 <- length(which(a3 < (floor(nchar(myseq)/3)-2)))
        b <- data.frame(res=c(l1,l2,l3),n=1:3,stringsAsFactors = F) %>% dplyr::filter(res==0)
      } else{
        a4 <- which(s2c(as.character(translate(DNAString(revcom(myseq)), no.init.codon=TRUE,if.fuzzy.codon="X"))) == "*")
        a5 <- which(s2c(as.character(translate(DNAString(substr(revcom(myseq),2,nchar(myseq))), no.init.codon=TRUE,if.fuzzy.codon="X"))) == "*")
        a6 <- which(s2c(as.character(translate(DNAString(substr(revcom(myseq),3,nchar(myseq))), no.init.codon=TRUE,if.fuzzy.codon="X"))) == "*")
        
        l1 <- length(which(a4 < (floor(nchar(myseq)/3)-2)))
        l2 <- length(which(a5 < (floor(nchar(myseq)/3)-2)))
        l3 <- length(which(a6 < (floor(nchar(myseq)/3)-2)))
        b <- data.frame(res=c(l1,l2,l3),n=4:6,stringsAsFactors = F) %>% dplyr::filter(res==0)
      }
      
      if(nrow(b)==0 | abs(inftmp[z,9] - inftmp[z,10])<100) {
        c <- 0
      } else if (b$n==1) {
        c <- substr(myseq,1,nchar(myseq))
      } else if (b$n==2) {
        c <- substr(myseq,2,nchar(myseq))
      } else if (b$n==3) {
        c <- substr(myseq,3,nchar(myseq))
      } else if (b$n==4) {
        c <- substr(revcom(myseq),1,nchar(myseq))
      } else if (b$n==5) {
        c <- substr(revcom(myseq),2,nchar(myseq))
      } else if (b$n==6) {
        c <- substr(revcom(myseq),3,nchar(myseq))
      } 
      if(c!=0) {
        a2 <- substring(as.character(c),seq(1,(nchar(c)-2),by=3),seq(3,nchar(c),by=3)) %>% table() %>% as.data.frame()
        names(a2)[1] <- "V1"
        a3 <- codon %>% dplyr::filter(aa!="W")
        a3$fre <- a2$Freq[match(as.vector(a3$codon),as.vector(a2$V1))]
        a3$fre[is.na(a3$fre)] <- 0
        fin <- a3 %>% cbind(data.frame(gene,virusid,stringsAsFactors = F)) 
      } else {
        fin <- data.frame(stringsAsFactors = F,aa=99,codon=99, fre=99,gene,virusid)
      }
      fin
    }) %>% rbind.fill() %>% dplyr::filter(aa!=99)
    ###
    if(nrow(yfre)==0){
      dd <- data.frame(virusid="999",gene="999",type="999",Dp=999,len=999,stringsAsFactors = F)
    } else {
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
      yij <- rbind(cbind(yijall %>% as.data.frame(),data.frame(gene=unique(yijall$virusid),type="all",len=3*sum(yfre$fre),stringsAsFactors = F)),
                   cbind(yijgene %>% as.data.frame(),data.frame(type="eachgene",len=1,stringsAsFactors = F)))
      
      dd <- lapply(1:length(unique(yij$gene)), function(j){
        ge <- unique(yij$gene)[j]
        yijj <- yij %>% dplyr::filter(gene == ge)
        
        hum$yij <- yijj$yij[match(hum$codon,yijj$codon)]
        hum$yij[which(is.na(hum$yij))] <- 0
        hum %>% 
          group_by(aa) %>%
          dplyr::summarize(Di=(sum((xij-yij)^2))^0.5) %>%
          as.data.frame() %>%
          dplyr::summarize(Dp = prod(Di)^(1/length(Di))) %>% cbind(unique(yijj[,c(3,5,6,7)]))
        
      }) %>% rbind.fill()
    }
    
  } else {
    dd <- data.frame(virusid="999",gene="999",type="999",Dp=999,len=999,stringsAsFactors = F)
  }
  dd
}) %>% rbind.fill() %>% dplyr::filter(Dp!=999)
hadvfile[n]
save(Dp,file = "~/virus/earlygene/hadv/Dp.hadvG.Rdata")

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
load("~/virus/earlygene/hadv/Dp.hadvG.Rdata")
vg <- Dp %>% cbind(data.frame(virusname="HAdV-G"))

##merge data
genehadv <- rbind(va,vb,vc,vd,ve,vf,vg) %>% dplyr::filter(type=="eachgene") 
load("~/virus/earlygene/hadv/hadvgene.Rdata")
genehadv$type2 <- hadvgene$type[match(as.vector(genehadv$gene),hadvgene$gene)]

##picture and test
##1)violin E,I,L
genehadv$virusname <- as.vector(genehadv$virusname)

lapply(1:length(unique(genehadv$virusname)),function(x){
  v <- unique(genehadv$virusname)[x]
  a <- genehadv %>% dplyr::filter(virusname == v)
  
  E <- a %>% dplyr::filter(type2=="E")
  I <- a %>% dplyr::filter(type2=="I")
  L <- a %>% dplyr::filter(type2=="L")
  
  data.frame(v,EIp=wilcox.test(E$Dp,I$Dp)$p.value,ELp=wilcox.test(E$Dp,L$Dp)$p.value,ILp=wilcox.test(I$Dp,L$Dp)$p.value,
             mE=mean(E$Dp),mI=mean(I$Dp),mL=mean(L$Dp),medE=median(E$Dp),medI=median(I$Dp),medL=median(L$Dp))
}) %>% rbind.fill()

##2)violin, E blongs to nonstructural gene; I and L belong to structural; U means unkonwn 
genehadv <- genehadv %>% dplyr::filter(gene!="U")
genehadv$type3 <- "nonstructure"
genehadv$type3[c(which(substr(genehadv$gene,1,1)=="L"),which(genehadv$gene %in% c("IX","IVa2")))] <- "structure"

lapply(1:length(unique(genehadv$virusname)),function(x){
  v <- unique(genehadv$virusname)[x]
  a <- genehadv %>% dplyr::filter(virusname == v) 
  ns <- a %>% dplyr::filter(type3=="nonstructure")
  s <- a %>% dplyr::filter(type3=="structure")  
  data.frame(v,p=wilcox.test(ns$Dp,s$Dp)$p.value,mns=mean(ns$Dp),ms=mean(s$Dp),medns=median(ns$Dp),meds=median(s$Dp),stringsAsFactors = F)
}) %>% rbind.fill()
