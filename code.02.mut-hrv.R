library(seqinr)
library(Biostrings)
library(seqRFLP)
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
seqfile <- system("ls ~/virus/earlygene/hrv/mut*.fasta",intern = T) 
inffile <- system("ls ~/virus/earlygene/hrv/inf*",intern = T)
codon <- read.table("/mnt/data4/disk/chenfeng/codonpaper/expriment/fcs/codon.txt",header = TRUE,stringsAsFactors = F)
n=3

myseq <- readDNAStringSet(seqfile[n]) 
cdsinf <- read.csv(inffile[n],stringsAsFactors = F,sep = "\t",header = F)

load("~/codonpaper/code and data/1_2.xijofhost.finall.Rdata")
hum <- xijofhost %>% dplyr::filter(special =="Homo_sapiens")
Dp <- mclapply(mc.cores = 60,1:length(myseq),function(y){
  a <- myseq[y]
  query <- strsplit(names(a)," ")[[1]][1]
  inftmp <- cdsinf %>% dplyr::filter(V1==query)
  if(nrow(inftmp)>0){
    yfre <- lapply(1:nrow(inftmp),function(z){
      
      gene <- strsplit(inftmp[z,2],",")[[1]]
      virusid <- inftmp[z,]$V1
      myseq <- as.character(substr(as.character(a),inftmp[z,7],inftmp[z,8]))
      ##
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
    
    ##
    if(nrow(yfre)==0){
      dd <- data.frame(virusid="999",gene="999",type="999",Dp=999,len=999,stringsAsFactors = F)
    } else{
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
seqfile[n]
save(Dp,file = "~/virus/earlygene/hrv/Dp.HRV-C.Rdata")

##load data

#load("~/virus/earlygene/hrv/Dp.HRV.Rdata")
load("~/virus/earlygene/hrv/Dp.HRV-A.Rdata")
va <- Dp %>% cbind(data.frame(virusname="HRV-A"))
load("~/virus/earlygene/hrv/Dp.HRV-B.Rdata")
vb <- Dp %>% cbind(data.frame(virusname="HRV-B"))
load("~/virus/earlygene/hrv/Dp.HRV-C.Rdata")
vc <- Dp %>% cbind(data.frame(virusname="HRV-C"))

##merge data
#ful <- (rbind(va,vb,vc,vd,ve,vf,vg) %>% dplyr::filter(len>200))$virusid %>% unique()


gene <- rbind(va,vb,vc) %>% dplyr::filter(type=="eachgene") %>% separate(gene,c("gene","type2"),sep = "_")
gene$type2[which(gene$type2 == "Non-tructure")] <- "Non-structure"
source("~/Rfunction/style.print.R")

####cumulative distribution and test

gene %>% 
  ggplot(aes(x=type2,y=Dp,color=type2))+
  geom_violin(scale = "width") +
  geom_boxplot(width = 0.1, outlier.colour = NA) +
  stat_summary(fun = 'mean', geom = 'point', shape = 18, colour = 'black')+
  facet_wrap(~virusname,ncol = 4)+
  labs(y=expression(paste(italic(D)[P]," (virus/human)")),x="Frequency")+
  scale_y_continuous(limits = c(0.16,0.72),breaks = c(0.2,0.4,0.6))+
  style.print()+theme(legend.position = "none")
##test
gene$virusname <- as.vector(gene$virusname)
lapply(1:length(unique(gene$virusname)),function(x){
  v <- unique(gene$virusname)[x]
  a <- gene %>% dplyr::filter(virusname == v)
  s <- a %>% dplyr::filter(type2=="Structure")
  n <- a %>% dplyr::filter(type2=="Non-structure")
  data.frame(v,p=wilcox.test(s$Dp,n$Dp)$p.value,ms=mean(s$Dp),mn=mean(n$Dp))
}) %>% rbind.fill()

## test1
gene %>% group_by(virusname,gene,type2) %>% dplyr::summarize(n=length(Dp)) %>% as.data.frame() %>% arrange(type2) %>% arrange(desc(virusname))
##test2
ref <- system("ls ~/virus/earlygene/hrv/ref*",intern = T)
mclapply(mc.cores = 3,1:3,function(x){
  a <- readDNAStringSet(ref[x])
  
  
  mclapply(mc.cores = 2,1:length(a),function(y){
    
    
    gene <- names(a[y])
    leng <- width(a)[y]
    data.frame(gene,leng,stringsAsFactors = F)
    
  }) %>% rbind.fill() %>% group_by(gene) %>% dplyr::filter(leng==max(leng)) %>%
    cbind(virus=strsplit(ref[x],"/")[[1]][9])
  
  
})%>% rbind.fill()
