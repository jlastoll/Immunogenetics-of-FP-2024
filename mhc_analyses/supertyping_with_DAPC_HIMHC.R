library(adegenet)
library(ggplot2)
library(reshape2)
library(scales)
library(unikn)
library(here)

# Data frame of amino acid residue values for 23 MHC alleles
sea_turtle_full <- read.csv("HI_MHC_amino_acid_matrix_values_for_DAPC_newdata.csv",header=TRUE)
sea_turtle_full <- sea_turtle_full[,-1] # remove first column (allele identifiers) so that it can be read into find.clusters (below)


### K MEANS CLUSTERING

# "K" is the number of clusters to group the data into and so choosing a k is very important. First, perform K-means clustering over a number of values of K and repeat 10 times for each value to explore each k.

maxK <- 11 # argument for max.n.clust, which is an integer indicating the maximum number of clusters to be tried; find.clusters will evaluate 1 through k clusters (here, 11)

myMat <- matrix(nrow=100, ncol=maxK) # create empty matrix to be filled by the for loop below
colnames(myMat) <- 1:ncol(myMat) # give column names to matrix


for(i in 1:nrow(myMat)){
  grp <- find.clusters(sea_turtle_full, n.pca = 50, choose.n.clust = FALSE,  max.n.clust = maxK) # retains 50 PCAs (n.pca), user doesn't choose cluster number (choose.n.clust = FALSE), and the maximum number of clusters to be tried (max.n.clust) is 11 (so it does 1-11)
  myMat[i,] <- grp$Kstat # fill matrix with Kstats from groups
}

# Visualizing k means clustering

my_df <- melt(myMat) # turn matrix into a df
colnames(my_df)[1:3] <- c("Group", "K", "BIC") # column headings for df
my_df$K <- as.factor(my_df$K) # make the K values a factor
head(my_df) # lists the groups per each K and the BIC

p1 <- ggplot(my_df, aes(x = K, y = BIC))
p1 <- p1 + geom_boxplot()
p1 <- p1 + theme_bw()
p1 <- p1 + xlab("Number of groups (K)")
p1 # shows how BIC changes with the number of groups (K) that's chosen

## There is no clear elbow point of minimizes BIC values and my curve is not really elbow shaped, but 2-5 seems to encompass the majority of the change, so evaluate k = 2:5

### CROSS-VALIDATION
# This is a two-step validation procedure, where the number of principle components to retain is evaluated at each k 3:11
# Prior to validation, run find.clusters to get the group membership (required for the validation steps)
# Then, run a broad validation step using xvalDAPC, where 1-50 principle components are retained for 30 iterations. This will yield a principle component number that minimizes the mean square error (and is therefore a good choice for that k value)
# Finally, run a narrower validation step using xvalDAPC, but this time instead of retaining 1-50 principle components, center the calculation on the result from above. For example, if PC = 10 was the value that minimized MSE, then run the second validation on PCs 1 to 20 (where 10 is the center of that distribution), and run for 1000 iterations.
# The result will be the number of PCs to retain at a particular k, which can then be used to run DAPC().


# k = 2
# Set up group membership
grp_2 <- find.clusters(sea_turtle_full, n.pca = 50, n.iter = 100, n.clust = 2)

set.seed(1)
xval_2 <- xvalDapc(sea_turtle_full, grp = grp_2$grp, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 100, xval.plot = TRUE) # leaving the n.pca.max default 

# xval gives several metrics on which number of PCs is best but these metrics aren't always congruent with one another. We used the number of PCs asociated iwth the lowest RMSE as the "optimum" number of PCAs in the DAPC analysis.- here that value is 10

xval_2[2:6] # number of PCs achieving lowest MSE is 10

# 3_a: evaluate 10 as center of optimal distribution:

set.seed(1)
xval_2_a <- xvalDapc(sea_turtle_full, grp = grp_2$grp, n.pca = 1:20, n.rep = 1000, parallel = "multicore", ncpus = 4L) # evaluate 1 to 20 PCAs, since 10 is the center.

xval_2_a[1:6] # 10 is optimal number of PCAs at 2 clusters for correctly predicting subsamples with the lowest error.

# Result: 2 clusters, 10 PCs

dapc_2 <- dapc(sea_turtle_full, grp_2$grp, n.pca = 10) # retain 2 DFs

# Visualize DAPC results

pal2 <- seecol(pal_unikn_pref, n = 2) # set up palette

k2_scatter <- scatter(dapc_2, # data
                                bg = "white", # white background
                                pch = 20, # size and shape of points
                                cstar = 0, # cstar= lines btwn points, 0 for null
                                solid = 0.6,
                                cex = 1.5, #
                                clab = 0.375,
                                leg = TRUE,
                                scree.da = TRUE,
                                scree.pca = TRUE,
                                col = pal2)



# Compoplots shows us the alleles that have less than a certain probability of membership in a single cluster; i.e., those alleles that might belong to both. For supertyping, we want to reduce this down to having NO alleles that are like this, since we're trying to break them down into the least divisible unit.

# Alleles with less than 100% probability of membership:
temp90 <- which(apply(dapc_2$posterior,1, function(e) all(e<0.90))) # all those that have less than 90% probability of membership in a single cluster.
temp90 # none

temp95 <- which(apply(dapc_2$posterior,1, function(e) all(e<0.95))) # all those that have less than x% probability of membership in a single cluster.
temp95 # none

temp99 <- which(apply(dapc_2$posterior,1, function(e) all(e<0.99))) # all those that have less than x% probability of membership in a single cluster.
temp99 # none

temp100 <- which(apply(dapc_2$posterior,1, function(e) all(e<0.9999)))
# Visualize:
compo2 <- compoplot(dapc_2, col = pal2)

##################################################
# k = 3
# Set up group membership
grp_3 <- find.clusters(sea_turtle_full, n.pca = 50, n.iter = 100, n.clust = 3) # retain 50 PCs. at k = 3 clusters

set.seed(1)
xval_3 <- xvalDapc(sea_turtle_full, grp = grp_3$grp, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 100, xval.plot = TRUE) # leaving the n.pca.max default so that it evaluates up to 124 PCs-- the number of alleles.

# xval gives several metrics on which number of PCs is best but these metrics aren't always congruent with one another. We used the number of PCs asociated iwth the lowest RMSE as the "optimum" number of PCAs in the DAPC analysis.

xval_3[2:6] # 6

# 3_a: evaluate 6 as center of optimal distribution:

set.seed(1)
xval_3_a <- xvalDapc(sea_turtle_full, grp = grp_3$grp, n.pca = 1:12, n.rep = 1000, parallel = "multicore", ncpus = 4L) # evaluate 1 to 20 PCAs, since 10 is the center.

xval_3_a[2:6] # 5 is optimal number of PCAs at 3 clusters for correctly predicting subsamples with the lowest error.

# Result: 3 clusters, 3 PCs

dapc_3 <- dapc(sea_turtle_full, grp_3$grp, n.pca = 6) #  retain 2 DFs

# Visualize DAPC results

pal3 <- seecol(pal_unikn_pref, n = 3) # set up palette

my_df <- as.data.frame(dapc_3[[ length(dapc_3) ]]$ind.coord)
my_df$Group <- dapc_l[[ length(dapc_l) ]]$grp
head(my_df)

pdf("DAPC_scatter.pdf",height=500, width=500)
k3_scatter <- scatter(dapc_3, # data
                      bg = "white", # white background
                      pch = 19, # size and shape of points
                      cstar = 1, # cstar= lines btwn points, 0 for null
                      solid = 1,
                      cex = 1.5, #
                      clab = 0.375,
                      leg = TRUE,
                      scree.da = FALSE,
                      scree.pca = TRUE,
                      posi.pca = "bottomright",
                      posi.leg = "topleft",
                      col = pal3)
dev.off()vet


pdf("DAPC3_scatter_9.23.23.pdf",height=10, width=12)

myCol<-c("#E0F3DB","#9ECAE1","#08306B" )
scatter(dapc_3, scree.da=FALSE, bg="white", pch=20, cell=2, cstar=0, col=myCol, solid=1,
        cex=4,clab=0, leg=TRUE, txt.leg=paste("Supertype",1:3))

dev.off()
# Compoplots shows us the alleles that have less than a certain probability of membership in a single cluster; i.e., those alleles that might belong to both. For supertyping, we want to reduce this down to having NO alleles that are like this, since we're trying to break them down into the least divisible unit.

# Alleles with less than 100% probability of membership:
temp90 <- which(apply(dapc_3$posterior,1, function(e) all(e<0.90))) # all those that have less than 90% probability of membership in a single cluster.
temp90 # none

temp95 <- which(apply(dapc_3$posterior,1, function(e) all(e<0.95))) # all those that have less than x% probability of membership in a single cluster.
temp95 # none

temp99 <- which(apply(dapc_3$posterior,1, function(e) all(e<0.99))) # all those that have less than x% probability of membership in a single cluster.
temp99 # none

# Visualize:

pdf("DAPC3_compoplot.pdf",height=10, width=10)
compo3 <- compoplot(dapc_3, col = myCol) #all alleles cluster to groups with probability of 1.0
dev.off()
###########################################################################################################################

# k = 4
# Set up group membership
grp_4 <- find.clusters(sea_turtle_full, n.pca = 50, n.iter = 100, n.clust = 4) 

set.seed(1)
xval_4 <- xvalDapc(sea_turtle_full, grp = grp_4$grp, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 100, xval.plot = TRUE) # leaving the n.pca.max default

xval_4[2:6] # 5 PCs

# 4_a: center xval on 5 PCs:
set.seed(1)
xval_4_a <- xvalDapc(sea_turtle_full, grp = grp_4$grp, n.pca = 1:10, n.rep = 1000, parallel = "multicore", ncpus = 4L)

xval_4_a[2:6] # 4 PCs

# Result: 4 clusters, 4 PCs

dapc_4 <- dapc(sea_turtle_full, grp_4$grp, n.pca = 4) # retain 1 PCs automatically; retain 3 DFs

pal4 <- seecol(pal_unikn_pref, n = 4)

k4_scatter <- scatter(dapc_4, # data
        bg = "white", # white background
        pch = 20, # size and shape of points
        cstar = 0, # cstar= lines btwn points, 0 for null
        solid = 0.6,
        cex = 1.5, #
        clab = 0.375,
        leg = TRUE,
        scree.da = TRUE,
        scree.pca = TRUE,
        col = pal4)

temp90 <- which(apply(dapc_4$posterior,1, function(e) all(e<0.90))) # all those that have less than 90% probability of membership in a single cluster.
temp90 # none

temp95 <- which(apply(dapc_4$posterior,1, function(e) all(e<0.95))) # all those that have less than x% probability of membership in a single cluster.
temp95 # none

temp99 <- which(apply(dapc_4$posterior,1, function(e) all(e<0.99))) # all those that have less than x% probability of membership in a single cluster.
temp99 # none

compo4 <- compoplot(dapc_4, col = pal4)
###################################

# k = 5
# Set up group membership
grp_5 <- find.clusters(sea_turtle_full, n.pca = 50, n.iter = 100, n.clust = 5) #  retain 50 PCs. at k = 5 clusters

set.seed(1)
xval_5 <- xvalDapc(sea_turtle_full, grp = grp_5$grp, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 100, xval.plot = TRUE) # leaving the n.pca.max default so that it evaluates up to 124 PCs-- the number of alleles.

xval_5[2:6] # 6 PCs

# 5_a: evaluate  as center of optimal distribution:
set.seed(1)
xval_5_a <- xvalDapc(sea_turtle_full, grp = grp_5$grp, n.pca = 1:12, n.rep = 1000, parallel = "multicore", ncpus = 4L) 

xval_5_a[2:6] # 1 is optimal number of PCAs at 5 clusters for correctly predicting subsamples with the lowest error.

# Result: 5 clusters, 1 PCs

dapc_5 <- dapc(sea_turtle_full, grp_5$grp, n.pca = 1) # retain 1 PCs automatically; retain 4 DFs

pal5 <- seecol(pal_unikn_pref, n = 5)

k5_scatter <- scatter(dapc_5, # data
        bg = "white", # white background
        pch = 20, # size and shape of points
        cstar = 0, # cstar= lines btwn points, 0 for null
        solid = 0.6,
        cex = 1.5, #
        clab = 0.375,
        leg = TRUE,
        scree.da = TRUE,
        scree.pca = TRUE,
        col = pal5)

temp90 <- which(apply(dapc_5$posterior,1, function(e) all(e<0.90))) # all those that have less than 90% probability of membership in a single cluster.
temp90 # 5 alleles

temp95 <- which(apply(dapc_5$posterior,1, function(e) all(e<0.95))) # all those that have less than x% probability of membership in a single cluster.
temp95 # 7 alleles

temp99 <- which(apply(dapc_5$posterior,1, function(e) all(e<0.99))) # all those that have less than x% probability of membership in a single cluster.
temp99 # 9 alleles

compoplot(dapc_5, col = pal5)


### K = 2 through 5 has been evaluated and 3 appears to be the optimal cluster number.
# Get a dataframe of each allele and its supertype
supertypes <- as.data.frame(dapc_3[["grp"]]) # gives supertype membership of each allele.
supertypes

allele_names <- read.csv("2023_outputs/HI_MHC_amino_acid_matrix_values_for_DAPC_newdata.csv")
allele_names <- allele_names[1]

MHC_alleles_by_supertype <- cbind(supertypes, allele_names)
MHC_alleles_by_supertype

write.csv(MHC_alleles_by_supertype, "2023_outputs/MHC_alleles_by_supertype_final.csv")



