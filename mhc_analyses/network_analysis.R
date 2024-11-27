library(here)
library("ape")
library("pegas")

read.dna("2023_outputs/derepseqsbyindiv.fasta",format = "fasta")->chmy_seqs #read in fasta file as a DNAbin frame, 33973 sequences
chmy_seqs #all seqs 190bp, stored in matrix, gets base composition for each base


chmy_haps<-haplotype(chmy_seqs)
chmy_haps #23 haplotypes, gives haplotype labels and frequencies (except in this case all frequencies are 1 since i'm loading only the unique sequences rather than all copies of sequences)

chmynet<-haploNet(chmy_haps)
plot(chmynet,xlim =c(81,82))
#plot(chmynet,size=attr(chmynet,"freq"),  scale.ratio = 2, cex = 0.8, fast=FALSE,show.mutation=2) #can make it so that circle correspond to circle frequency, but I didn't do that since all freqs are the same
#other plotting options:
#plot(chmynet) #looks terrible
(sz<-summary(chmy_haps)) #extract frequencies of haplotypes
(chmynet.labs<-attr(chmynet, "labels")) #get object of haplotype names that we'll reorder by frequency below
sz<-sz[chmynet.labs]
plot(chmynet,size=sz,show.mutation = 3) #plots bubble size based on number of individuals who have this haplotype

#set more options:
#setHaploNetOptions(haplotype.inner.color = "light blue",haplotype.outer.color = "navy",show.mutation = 3, labels = TRUE, scale.ratio=0.8,link.type=1,link.width=1)
setHaploNetOptions(show.mutation = 3, labels = TRUE, scale.ratio=0.6,link.type=1,link.width=1,pie.colors.function=topo.colors)
setHaploNetOptions()
P<-as.matrix(read.table("../MHC_analysis/variantfreqbyloc.csv"))#read in variant frequencies by location
plot(chmynet,pie=P,size=sz, cex = 1.5, fast=FALSE, legend= c(-60,80))
o<-replot(col.identifier = "purple",pch=30) #interactive replotting to move nodes



tiff("MHCnetworkwithlocs.tiff", width=2000, height=1500)
plot(chmynet,pie=P,size=sz*3, cex = 1.5, fast=FALSE, legend= c(-55,70),xy=o)
dev.off()


#include supertype information as color
supertypes<-as.data.frame(read.csv("2023_outputs/MHC_alleles_by_supertype_final.csv"))
cols<-c("74","")
setHaploNetOptions(show.mutation = 0,labels = TRUE, scale.ratio=1,link.type=1,link.width=1, haplotype.inner.color=adjustcolor(supertypes$color, alpha.f = 0.7))
library(RColorBrewer)
cols<-brewer.pal(3,"Set2")
pal<-colorR

pdf("2023_outputs/MHCnetworkwithsupertypes1.19.24.pdf", width=10, height=10)
plot(chmynet,size=sz*2.5, fast=FALSE,legend=FALSE, cex= 1.5)
legend(600,65,legend = c("Supertype 1","Supertype 2", "Supertype 3"), box.col= "white", fill = c('#E0F3DB','#9ECAE1','#08306B'), cex=1.5)
o<-replot() #interactive replotting to move nodes- HELPFUL TO DO THIS WHEN THE CIRCLES AREN'T THERE AND SAVE TO SUE HERE

dev.off()


library(ggplot2)

#get sequences corresponding to roman numeral labels:
h<-haplotype(chmy_seqs)
dna_seqs<-write.dna(h,"haplotypes11.25.txt",colsep = "\t",colw = 190) #colw needs to be 190 to get all of dna seq on one line
hap_seqs<-read.table("haplotypes11.25.txt",header = T)
colnames(hap_seqs)<-c("Haplotype_ID","DNA-seq")
hap_seqs