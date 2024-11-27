# Risk ratio of top alleles in C. mydas with FP
#https://www.cdc.gov/csels/dsepd/ss1978/lesson3/section5.html
#see link above for explanation of risk ratio. a value of RR=` means risk identical among grop, greater than 1 indicates increased risk for the primary group`
library(epiR)
library(tidyverse)
library(ggplot2)
library(svglite)
library(here)


data<-read.csv("./ALLMHCsamp-allvars.csv")
data<-na.omit(data)#wierd import issue created nas
data<-data[!(data$FP.status.binary=="u"),]
#data<-data[!grepl("FFS",data$sample),] #exclude FFS site data

# Subset to alleles that occur in 7 or more C. mydas
Cm_data <- subset(data, select = c(sample, FP.status.binary,cmhi_1,cmhi_2,cmhi_3,cmhi_4,cmhi_5,cmhi_6,cmhi_7,cmhi_8,cmhi_9,cmhi_10,cmhi_11,cmhi_12,cmhi_13,cmhi_14,cmhi_15))


Cm_data[Cm_data$FP.status.binary == 1, "FP.status.binary"] <- "FP+" # replace 1s with yes and 0s with nos for FP
Cm_data[Cm_data$FP.status.binary == 0, "FP.status.binary"] <- "FP-"


#loop through each allele to test association with FP

#build empty matrix
Cm_risks <-matrix(nrow=15,ncol=4) #rows are total number of alleles we're testing
colnames(Cm_risks)<-c("est","lower","upper","pval")
alleles<-colnames(Cm_data[3:17])
rownames(Cm_risks)<-alleles

#loop through every allele in dataframe
for (i in  colnames(Cm_data[3:17])) {
  Cm_data[Cm_data[[i]] == 1, i]<- paste(i," +") 
  Cm_data[Cm_data[[i]] == 0, i]<- paste(i," -") 
  x<-table(Cm_data[[i]], Cm_data$FP.status.binary)
  x2<-cbind(x[,2], x[,1]) # rearrange columns
  colnames(x2) <- c("FP+", "FP-")
  x2<-rbind(x2[2,], x2[1,]) # rearrange rows
  rownames(x2) <- c(paste(i," +"),paste(i, " -"))
  risk<-epi.2by2(x2, method = "cohort.count", conf.level = 0.95)
  RR<- risk$massoc.summary[1,2:4]
  p_corr<-risk$massoc.detail$chi2.strata.yates[1,4]
  RR$pval<-p_corr
  Cm_risks[i,1]<-RR[[1]]
  Cm_risks[i,2]<-RR[[2]]
  Cm_risks[i,3]<-RR[[3]]
  Cm_risks[i,4]<-RR[[4]]
}

write.csv(Cm_risks,"riskratioanalysis.csv")

#test by supertype:

#recode number of supertype alleles per individual to binary 0/1
data2<-data
data2$supertype.1<-lapply(data2$supertype.1, function(x) ifelse (x > 1, 1, x))
data2$supertype.2<-lapply(data2$supertype.2, function(x) ifelse (x > 1, 1, x))
data2$supertype.3<-lapply(data2$supertype.3, function(x) ifelse (x > 1, 1, x))

Cm_data <- subset(data2, select = c(sample, FP.status.binary,supertype.1,supertype.2,supertype.3))
Cm_data[3:5]<-as.data.frame(sapply(Cm_data[3:5], as.numeric))

Cm_data[Cm_data$FP.status.binary == 1, "FP.status.binary"] <- "FP+" # replace 1s with yes and 0s with nos for FP
Cm_data[Cm_data$FP.status.binary == 0, "FP.status.binary"] <- "FP-"

#loop through by supertype 1-3

#build empty matrix
sup_risks <-matrix(nrow=3,ncol=4) #rows are total number of supertypes we're testing
colnames(sup_risks)<-c("est","lower","upper","pval")
sups<-colnames(data2[26:28])
rownames(sup_risks)<-sups

#loop through every allele in dataframe
for (i in  colnames(Cm_data[3:5])) {
  Cm_data[Cm_data[[i]] == 1, i]<- paste(i," +") 
  Cm_data[Cm_data[[i]] == 0, i]<- paste(i," -") 
  x<-table(Cm_data[[i]], Cm_data$FP.status.binary)
  x2<-cbind(x[,2], x[,1]) # rearrange columns
  colnames(x2) <- c("FP+", "FP-")
  x2<-rbind(x2[2,], x2[1,]) # rearrange rows
  rownames(x2) <- c(paste(i," +"),paste(i, " -"))
  risk<-epi.2by2(x2, method = "cohort.count", conf.level = 0.95)
  RR<- risk$massoc.summary[1,2:4]
  p_corr<-risk$massoc.detail$chi2.strata.yates[1,4]
  RR$pval<-p_corr
  sup_risks[i,1]<-RR[[1]]
  sup_risks[i,2]<-RR[[2]]
  sup_risks[i,3]<-RR[[3]]
  sup_risks[i,4]<-RR[[4]]
}

#no significance

# graph a forest plot of alleles

ggplot(data=as.data.frame(Cm_risks), aes(x=row.names(Cm_risks), y=est, ymin=lower, ymax=upper)) +
  geom_pointrange(shape = 18) +
  geom_hline(yintercept=1, lty=2, color = "red") +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("allele") + ylab("relative risk estimate of FP (95% CI)") +
  theme_bw() # use a white background


