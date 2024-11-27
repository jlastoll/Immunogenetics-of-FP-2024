## packages###
library(here)
library(RColorBrewer)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("edgeR")
library(edgeR)
library(tidyr)
library(reshape2)
library(dplyr)
library(statmod)
#BiocManager::install("org.Hs.eg.db",version="3.16")
library(org.Hs.eg.db)
#BiocManager::install("topGO")
library(topGO)
library(GOplot)
library(ggplot2)
#BiocManager::install("ReactomePA")
library("ReactomePA")
library(data.table)
library(gplots)
#BiocManager::install('EnhancedVolcano')
library('EnhancedVolcano')
library('ggbreak') #need this to make breaks in axes of plots




#Read in counts file
x <- read.table("raw_counts/refguided_all.counts.matrix", header=TRUE,sep="")

#Read in samples file
samples <- read.delim("sample_names.txt", sep="", header = FALSE) #metadata file with sample name, file name, tumor score

#Rename column names for counts file to sample names
#Confirmed they are in the same order
colnames(x) <- samples$V2

#Only keep relevant samples and columns
samples <- samples %>% dplyr::select(V2, V3) %>% 
  dplyr::filter(V3 != "not") %>% #the column that said "not using" was split into two columns when read in, so I just filtered based on the column that had "not"
  as_tibble()

#Remove extraneous samples from counts table
x <- x[,colnames(x) %in% samples$V2]
x$GENE<-rownames(x)
#Establish groups
group <- factor(samples$V3)

#Make design matrix
design <- model.matrix(~group)
rownames(design)<-samples$V2
#Make DGEList
y <- DGEList(counts=x, group=group,annotation.columns = "GENE")

#Save full count table so don't need to reload y if wanting to re-run things
y.full <- y

#Filter out lowly expressed genes that will mess with comparisons
keep <- filterByExpr(y, design = design) 
y <- y[keep,,keep.lib.sizes=FALSE] 

#Calculate normalization factors
y <- calcNormFactors(y)

# Visualize result of filters:

par(mfrow=c(1,2))
lcpm <- cpm(y.full, log=TRUE)
boxplot(lcpm, las=2, main="")
title(main="A. Unnormalised data", ylab="Log-cpm")

lcpm <- cpm(y, log=TRUE)
boxplot(lcpm, las=2, main="")
title(main="B. Normalised data", ylab="Log-cpm")


#Estimate dispersion for pairwise comparisons
y <- estimateDisp(y, design, robust=TRUE)
plotBCV(y)
y$common.dispersion #0.1114653
#Fit the model for treatments
fit1 <- glmQLFit(y, design, robust=TRUE) 

#Likelihood ratio test
lrt1 <- glmLRT(fit1)
tt_lrt1 <- topTags(lrt1,n=nrow(y))
head(tt_lrt1$table)
ALLGENES<-tt_lrt1$table
#write.csv(ALLGENES, "ALLGENES_LRT_10.17.23.csv", sep = "\t")
table(tt_lrt1$table$FDR< 0.05) #65 de genes
withannot<-tt_lrt1$table[tt_lrt1$table$FDR< 0.05,] #just keep de genes

#Inspect results
o <- order(lrt1$table$PValue)
toplcpm<-cpm(y,log = TRUE)[o[1:65],]
topcpm<-cpm(y)[o[1:65],]
# for the filtering it seems to be keeping genes that have cpm of at least 0.45 in at least 5 samples, and that at least some samples have a minimum count of 10. This is appropriate based on the published EdgeR paper

#Do some checks:
#top DE gene:LOC102939794 
#mean cpm in tumor group
#mean cpm in control group
tum_meancpm<-mean(topcpm[1,c(samples$V2[samples$V3=="tumor"])]) #-0.8384349
cont_meancpm<-mean(topcpm[1,c(samples$V2[samples$V3=="control"])]) #48.9577
log(mean(topcpm[1,]),2)#4.637, very close to the logCPM reported in table
#manually compute log2FC between groups:
log((tum_meancpm/cont_meancpm),2) #-5.867693, downregulated in tumor group

#the toptable gives the log2CPM using ALL SAMPLES which is kind of dumb?
#log fold change is : log of tumor group cpm - log of control group cpm for a gene, so this all checks out

#write.csv(withannot, "sigDEresults_10.15.23.csv", sep = "\t")


################################ GENE ANNOTATION ######################################################
#read in manual annotations (for the LOC gene symbols in our reference I identified the gene symbol by looking at the reference genome annotation file, just using grep to look for LOC id )
manannot<-read.csv("ALLGENES_LRT_10.17.23.csv")
degenes<-manannot$GENE.LOOKUP #get vector of gene names for goana function, must be length of nrow(lrt) object,  12329 rows

#need to get entrezids for GO analysis
annot<-AnnotationDbi::select(org.Hs.eg.db,keys = degenes, columns = c("ENTREZID","SYMBOL","GENENAME"),keytype = "SYMBOL") #generates table with gene symbol and entrezid, but its 12331 rows see if there are dups:
dups<-annot$SYMBOL[duplicated(annot$SYMBOL)] #find duplicated genes, now there are 9, some of these are true duplicated genes 
truedups<-manannot$GENE.LOOKUP[duplicated(manannot$GENE.LOOKUP)] #7 true duplicate gene symbols
diff<-setdiff(dups,truedups)#get difference between these two vectors to find the replicated (), there are 2
#based on our genome annotation we want MEMO1 that is mediator of cell motility entrezid 51072
#for TEC we wante tec proteintyrosine kinase, entrezid 7006
#remove these particular 2 rows:
annot2<-annot[!(annot$SYMBOL=="TEC" & !(annot$ENTREZID=="7006")),]
annot2<-annot2[!(annot2$SYMBOL=="MEMO1" & !(annot2$ENTREZID=="51072")),]
dups<-annot2$SYMBOL[duplicated(annot2$SYMBOL)] #Double check that we removed the dups, yes we did
#add these new entrez IDs to the manannot table that contains FDR values
manannot$ENTREZID<-annot2$ENTREZID[match(manannot$GENE.LOOKUP,annot2$SYMBOL)]
manannot$GENENAME<-annot2$GENENAME[match(manannot$GENE.LOOKUP,annot2$SYMBOL)]




################################    DGE plots   #######################################################


de_gene_names<-row.names(toplcpm) #get just DE gene names
toplcpm<-as.matrix(toplcpm)
#to determine range and breaks, look over the distributions:
summary(toplcpm[,1:10])
ggplot(stack(as.data.frame(toplcpm)),aes(x=ind,y=values))+
  geom_boxplot()
#this creates barplots by saple with lcpms of values
#min of all treatments is -3.448, max is 9.3579
#mostly even distribution,



## violin plot

top5<-head(manannot,n=5L) #subset out top 5 genes
top10<-head(manannot,n=10L)#subset top 10


#DE_lcpms<-read.csv("rawlcpmcounts_3.17.22.csv")
DE_lcpms3<-subset(toplcpm, rownames(toplcpm) %in% top10$reference.gene.name)
DE_lcpms3<-as.data.frame(t(DE_lcpms3))
colnames(DE_lcpms3)<-top10$GENE.LOOKUP
sampleids<-rownames(DE_lcpms3)
DE_lcpms3<-tibble::rownames_to_column(DE_lcpms3,"sample")
DE_lcpms3$group<-group

DE_lcpms4<-gather(DE_lcpms3,c(top10$GENE.LOOKUP),key="gene",value="lcpm")
DE_lcpms4$lcpm<-as.numeric(DE_lcpms4$lcpm)

ggplot(DE_lcpms4, aes(x=group,y=lcpm,,fill=group))+
  geom_boxplot()+
  facet_grid(.~gene)+
  theme_bw()

group.colors<-c("#C7E9B4","#006D2C")


pdf("violintop10_1.18.pdf", width=12, height=7)

ggplot(DE_lcpms4, aes(x=group,y=lcpm,fill=group))+
  geom_violin()+
  facet_grid(.~gene)+
  theme_bw()+
  scale_fill_manual(values=group.colors)+
  stat_summary(fun.data = "mean_cl_boot", geom = "pointrange",
               colour = "black")
dev.off()

#plot 10 immunge genes IDENTIFIED IN SIG DGE GENES THROUGH LIT SEARCH AND GO ANALYSIS
immunegenes<-c('IFI27L2','HLA-DRA','C3AR1','STOM','TRIM27','CHIT1','GZMB','HEBP2','SDC3')
immuneonly<-manannot[c(manannot$GENE.LOOKUP %in% immunegenes),]

DE_lcpms3<-subset(toplcpm, rownames(toplcpm) %in% immuneonly$reference.gene.name)
DE_lcpms3<-as.data.frame(t(DE_lcpms3))
colnames(DE_lcpms3)<-c('IFI27L2','HLA-DRA_1','C3AR1','STOM','TRIM27','HLA-DRA_2','CHIT1','GZMB','HEBP2','SDC3')
sampleids<-rownames(DE_lcpms3)
DE_lcpms3<-tibble::rownames_to_column(DE_lcpms3,"sample")
DE_lcpms3$group<-group

DE_lcpms4<-gather(DE_lcpms3,c('IFI27L2','HLA-DRA_1','C3AR1','STOM','TRIM27','HLA-DRA_2','CHIT1','GZMB','HEBP2','SDC3'),key="gene",value="lcpm")
DE_lcpms4$lcpm<-as.numeric(DE_lcpms4$lcpm)

pdf("violinimmune_1.18.pdf", width=12, height=7)

ggplot(DE_lcpms4, aes(x=group,y=lcpm,fill=group))+
  geom_violin()+
  facet_grid(.~gene)+
  theme_bw()+
  scale_fill_manual(values=group.colors)+
  stat_summary(fun.data = "mean_cl_boot", geom = "pointrange",
               colour = "black")
dev.off()



##### Volcano Plot of DE genes#####
library(ggbreak) 
library(patchwork)
library(ggplot2)

res<-manannot[,c(2,5,9)] #keep gene symbol, logfc, FDR
res$GENE.LOOKUP<- toupper(make.unique(as.character(res$GENE.LOOKUP), sep = "_"))#make all symbols unique

#extract gene symbols with sig FDR, and logFC greater than 2
gene_labs<-filter(res, FDR<0.05, abs(logFC) > 2, !grepl("LOC", GENE.LOOKUP))
lab_italics <- paste0("italic('", res$GENE.LOOKUP, "')")
selectLab_italics<- paste0(
  "italic('",gene_labs$GENE.LOOKUP, "')")

# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
# this can be achieved with nested ifelse statements
keyvals <- ifelse(
  res$logFC < -1 & res$FDR<0.05, "#A6BDDB",
  ifelse( res$logFC > 1 & res$FDR<0.05, "#016C59",
          'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == "#016C59"] <- 'Significantly upregulated'
names(keyvals)[keyvals == "black"] <- 'Not significant'
names(keyvals)[keyvals == "#A6BDDB"] <- 'Significantly downregulated'

pdf("volcano_all_1.19.23_breaks.pdf", width=15, height=10)
par(cex.main=3)


EnhancedVolcano(res,
                lab = lab_italics,
                subtitle = NULL,
                x = 'logFC',
                y = 'FDR',
                selectLab = selectLab_italics, 
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~ 'FDR'),
                pCutoff = 10e-3,
                FCcutoff = 1.0,
                pointSize = 3.0,
                labSize = 5.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                parseLabels = TRUE,
                min.segment.length = 0.5,
                arrowheads = FALSE,
                colCustom = keyvals,
                colAlpha = 4/5,
                legendPosition = 'bottom',
                legendLabSize = 13,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black') + coord_flip()+
  scale_y_break(c(10, 20),scales=1.5,space=0.5)



dev.off()

scale_y_break(c(10, 20),scales=0.5,space=0.5)
#make again but without outlier gene to see more of the plot then can paste together in illustrator:
res_nooutlier<-res[5:13399,]#removed the one gene that had a a FDR of e-90, to make a subsetted plot for viz purploses
gene_labs<-filter(res_nooutlier, FDR<0.05, abs(logFC) > 2, !grepl("LOC", GENE.LOOKUP))
lab_italics <- paste0("italic('", res_nooutlier$GENE.LOOKUP, "')")
selectLab_italics<- paste0(
  "italic('",gene_labs$GENE.LOOKUP, "')")

# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
# this can be achieved with nested ifelse statements
keyvals <- ifelse(
  res_nooutlier$logFC < -1 & res_nooutlier$FDR<0.05, "#A6BDDB",
  ifelse( res_nooutlier$logFC > 1 & res_nooutlier$FDR<0.05, "#016C59",
          'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == "#016C59"] <- 'Significantly upregulated'
names(keyvals)[keyvals == "black"] <- 'Not significant'
names(keyvals)[keyvals == "#A6BDDB"] <- 'Significantly downregulated'

pdf("volcano_subset_1.19.24.pdf", width=15, height=10)
par(cex.main=3)

EnhancedVolcano(res_nooutlier,
                lab = lab_italics,
                x = 'logFC',
                y = 'FDR',
                selectLab = selectLab_italics, 
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~Log[10]~ 'FDR'),
                xlim = c(-7,5),
                ylim = c(-5,10,70),
                pCutoff = 10e-3,
                FCcutoff = 1.0,
                pointSize = 3.0,
                labSize = 6.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                parseLabels = TRUE,
                colCustom = keyvals,
                colAlpha = 4/5,
                legendPosition = 'bottom',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black') + coord_flip()


#need to set range points for y axis 
dev.off()




############################################################# GO ANALYSIS###############################################################




#Run topGo package to get GO mappings and enrichment tests:
topDiffGenes<-function(allScore) {
  return(allScore <0.05)
} #function that pulls top DE genes based on FDR
genelist <- manannot$FDR # these 2 lines create genelist for topgo below
names(genelist) <- manannot$ENTREZID

entrez2go<-annFUN.org("BP",mapping= "org.Hs.eg.db", ID="entrez")
GOdata<-new("topGOdata",
            ontology="BP",
            allGenes=genelist,
            geneSelectionFun=topDiffGenes,
            annot=annFUN.org,
            mapping="org.Hs.eg.db",
            ID="entrez",
            nodeSize=5)

resultKS <- runTest(GOdata, algorithm = "weight01", statistic = "ks")
resultKS
resultKSelim<-runTest(GOdata, algorithm = "elim", statistic = "ks")
resultKSelim


tab <- GenTable(GOdata, KS = resultKSelim, topNodes = length(resultKSelim@score), numChar = 120)
tab <- tab[,c(1,2,3,6)]
names(tab)[4] <- "Raw.P.Value"
genes.in.term <- genesInTerm(GOdata)
genes.in.term2 <- unlist(lapply(genes.in.term, function(x)paste(x, collapse = ",")))
genes.in.term2 <- cbind(names(genes.in.term), genes.in.term2)#the three above lines retrieve and converting the genes that have each GO term to then combine with results below (I think)
tmp <- match(tab$GO.ID, genes.in.term2[,1])
tab$Genes.In.Term <- genes.in.term2[tmp,2]


tab2 <- GenTable(GOdata, KS = resultKS, topNodes = length(resultKS@score), numChar = 120)
tab2 <- tab2[,c(1,2,3,6)]
names(tab2)[4] <- "Raw.P.Value"
genes.in.term <- genesInTerm(GOdata)
genes.in.term2 <- unlist(lapply(genes.in.term, function(x)paste(x, collapse = ",")))
genes.in.term2 <- cbind(names(genes.in.term), genes.in.term2)#the three above lines retrieve and converting the genes that have each GO term to then combine with results below (I think)
tmp <- match(tab2$GO.ID, genes.in.term2[,1])
tab2$Genes.In.Term <- genes.in.term2[tmp,2]

#####Make circle plot of enriched GO terms: ####

tab$Category <- rep("BP", length(tab$GO.ID))
#re-order and rename columns to suit GOplot 
tab <- tab[,c(6,1,2,4,5)]
colnames(tab) <- c("Category", "ID", "Term", "adj_pval", "genes")
tab$genes<-as.character(tab$genes)
tab$adj_pval<-as.numeric(tab$adj_pval)



## create a circle object 

manannot$ID <- manannot$ENTREZID
row.names(manannot) <- c()
manannot2 <- manannot[,c(10, 5:9)]
manannot2<-na.omit(manannot2)


circ <-circle_dat(tab, manannot2) 
#COLUMN NAMES ALL NEED TO BE VERY SPECIFIC- SEE HERE:https://github.com/wencke/wencke.github.io/blob/91963feeaf4847f1ccff049511fb230eac153cc2/R/GOCore.R

## plot GoCircle 

cbPalette1 <- brewer.pal(n=9, name = 'PuBuGn')
cbPalette2<- brewer.pal(n=9,name='YlGnBu')

pdf("circleplottop12_1.18.23.pdf",width=15,height=10)



GOCircle(circ, n=12,zsc.col = c("#081D58", "#7FCDBB", "#FFFFD9" ),lfc.col =c("#016C59","#D0D1E6"))

dev.off()


#GO plot with terms of interest
IDs <- c('GO:0030198',
         'GO:0016192',
         'GO:0016049',
         'GO:1990748',
         'GO:0008152',
         'GO:0045087',
         'GO:0002250',
         'GO:0006954',
         'GO:0006950',
         'GO:0098609')

pdf("circleplotinterestgenes_1.18.23.pdf",width=15,height=10)

GOCircle(circ, nsub = IDs,zsc.col = c( "#7FCDBB", "#FFFFD9" ),lfc.col =c("#016C59","#D0D1E6"))

dev.off()

#generate table with zscores from circ object:

zscore<-
circ %>% dplyr::select(term,zscore) %>% distinct() %>% dplyr::left_join(tab, join_by("term"=="Term"))
write.csv(zscore,"GSEAzscore_1.23.24.csv") #need to go into excel and import csv that way otherwise the gene names don't get properly imported


### get gene info from the sig go term

de_withGO<- AnnotationDbi::select(org.Hs.eg.db,keys = manannot$GENE.LOOKUP[manannot$FDR<0.05], columns = c("ENTREZID","SYMBOL","GENENAME","GO"),keytype = "SYMBOL") %>% dplyr::filter(ONTOLOGY=="BP") %>% dplyr::select(SYMBOL,GO) %>% dplyr::group_by(SYMBOL) %>% summarize(GO=paste(GO, collapse = ";")) #generate a table with GO terms for biological processses for the DE genes, with all GO terms collapsed into one cell 
write.csv(de_withGO,"de_withGOterms.csv")



####################################### Reactome pathway enrichment analysis#####################################################
#run GSEA with reactome pathways, manually curated and reviewed pathways, instead of GO terms
#use the manannot file I read in for go enrichment above, then make sure to run code to look up entrez ids for the gene symbols which we'll also need here

d<-manannot[,c(3,9)] #subset to make gene list that is the  entrez id, FDR needs to be ordered and ranked
d<-na.omit(d) #remove NAs leave 9439 genes
d<-d[!duplicated(d$ENTREZID),]#remove duplicated values (not supposed to have duplicate gene ids, even though this will remove some DE genes) drops down to 9416 genes
geneList=d[,2] #numeric vector
names(geneList)=as.character(d[,1])#name vector with entrez IDs
geneList=sort(geneList, decreasing=TRUE) #make sure sorted in decreasing order
y<-gsePathway(geneList,
              organism="human",
              pvalueCutoff = 0.5,
              verbose = TRUE
)
head(y)

