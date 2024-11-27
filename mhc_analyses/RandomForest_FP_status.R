## Classification random forest for FP occurence in C. mydas

library(randomForest)
library(caret)
library(dplyr)
library(here)

# read in data

data_full <- read.csv("ALLMHCsamp-allvars.csv") #csv file with sample id, FP status, columns for each allele where the value is binary occurance, plus whatever other variables of interest
data_full<-na.omit(data_full)#wierd import issue created nas
data_full<-data_full[!(data_full$FP.status.binary=="u"),] #drop unknown FP individuals
#data_full<-data_full[!(data_full$location=="FFS"),] #can drop FFS location but don't recommend- very little predictive power without it

# variables: sample ID, location, spp, year, FP, paps_regressed, season, allele count, all alleles, supertypes A, B, C, and SCL ("Carap_L_sl" from database download)

# Convert location, FP, and all alleles, all supertypes to factors
data_full[,] <- lapply(data_full, factor) # location and species
data_full$sample<-lapply(data_full$sample, as.character)
#drop variables that only have one level (necessary to fit model)
lapply(data_full[,2:29], unique)
#drop alleles 23,22,20,16 since they are the same level across samples:
data_full<-data_full[,-which(names(data_full) %in% c("cmhi_16","cmhi_20","cmhi_22","cmhi_23"))]
#leaves 32 individuals, 27 variables

# Convert allele count to numeric
data_full$total.unique.vars<- as.numeric(data_full$total.unique.vars)

sapply(data_full, class)
sapply(data_full, levels)

data<-data_full[,c(1,23,4:22,24:27,2:3)] #order to have sample name, then fp status, then other variables


# Check for imbalance in response variable (FP)
length(which(data$FP.status.binary==0)) #39Cm without FP
length(which(data$FP==1)) #15 Cm with FP, for ~27.7% FP. This is unbalanced enough that we need to account for it by over-sampling the under-represented class (FP positive) and under-sampling the over-represented class (FP negative). Set the sampsize parameter to 2/3 of the under-represented class, and use that sample size in conjunction with the strata option when building the random trees and when training.

# there are 15 FP positive individuals; 15*(2/3) = 9.9 = 10 individuals, so we'll sample 10 FP positive individuals and 10 FP negative individuals for the training data set so that the resutling tree forest is not biased.
#when excluding the FFS individuals, this becomes 17 vs 15 so we can change the sample size and we'll just do 15 v 15
sample_size <- c(15,15) # will use subsequently with strata.

# Random Forest analysis

# Optimize mtry parameter. Default for RF classification is sqrt(p) where p is the number of variables in the dataframe.
# let's also run mtry values of sqrt(p), 2*sqrt(p), 0.1(p), 0.2(p), p/3, in addition to p, where p is the number of variables.
# Run each mtry at ntree= 100 to 1000 (by increments of 100), looking for plateau where the out of bag error rate (OOB-ER) is minimized while also maximizing the RF algorithm. The mtry value that minimizes the OOB-ER will be chosen for subsequent analyses.

# We are using 25 variables to explain FP occurence (this excludes the first column, which is sample ID, and FP occurence which is the response)

results_optimization <- matrix(data=NA , nrow = 0, ncol = 3) # create matrix that the following loop will dump info into
for (i in seq(from = 100, to = 1000 , by = 100)){  # values of ntree
  print(i)
  for (j in c(5, 10, 2.5, 5, 8.3, 25)){    #values of mtry based on 30 total predictors
    rf_ij <- randomForest(x = data[,3:27], y = as.factor(data$FP.status.binary), importance=TRUE ,proximity=TRUE, ntree=i, mtry=j, strata=as.factor(data$FP.status.binary), sampsize=sample_size)
    results_optimization <- rbind(results_optimization, c(i,j,tail(rf_ij$err.rate,1)[1]))
  }
}

# Clean up the file format
results_optimization<-as.data.frame(results_optimization)
colnames(results_optimization)<-c("ntree", "mtry","OOB_ER")

# Now plot results to see if there's a plateau
plot(results_optimization$ntree[results_optimization$mtry == 5],results_optimization$OOB_ER[results_optimization$mtry == 5], type="l", col="black", xlab="ntree",ylab="OOB-ER",ylim=c(0,1))
lines(results_optimization$ntree[results_optimization$mtry == 10],results_optimization$OOB_ER[results_optimization$mtry == 10], col="blue")
lines(results_optimization$ntree[results_optimization$mtry == 2.5],results_optimization$OOB_ER[results_optimization$mtry == 2.5], col="green")
lines(results_optimization$ntree[results_optimization$mtry == 5],results_optimization$OOB_ER[results_optimization$mtry == 5], col="purple")
lines(results_optimization$ntree[results_optimization$mtry == 8.3],results_optimization$OOB_ER[results_optimization$mtry == 8.3], col="orange")
lines(results_optimization$ntree[results_optimization$mtry == 25],results_optimization$OOB_ER[results_optimization$mtry == 25], col="red")

# all mtry parameters seem to behave similarly , but red line is best (value of 25) at minimizing error rate and is overall flattest
#when excluding FFS data, all the lines perform very similarly


# Create model with default paramters as a baseline with grid search, method = "oob"
control <- trainControl(method="oob", number=10, search = "grid")
seed <- 7
metric <- "Accuracy"
set.seed(seed)
mtry <- sqrt(25) # square root of total number of variables; this is the default mtry
tunegrid <- expand.grid(.mtry=mtry)
rf_default <- train(FP.status.binary ~ ., data=data[,2:27], method="rf", metric=metric, tuneGrid=tunegrid, trControl=control)
print(rf_default) # mtry = 5.099 at 83.3% accuracy.
#without FFS samples, mtry=5.099 and accuracy was 71.875% accuracy
#without sex as variable, mtry=5 and accuracy is 66.6%

# now let's fine tune
## Grid search, where method = "oob"
control <- trainControl(method="oob", number=10, search="grid")
set.seed(seed)
tunegrid <- expand.grid(.mtry=c(1:25)) # doing full 26 mtry
rf_gridsearch <- train(FP.status.binary~., data=data[,2:27], method="rf", metric=metric, tuneGrid=tunegrid, trControl=control)
print(rf_gridsearch)
plot(rf_gridsearch) # most accurate: mtry 3 at 87.0.0%
#when excluding the FFS samples, this was all over the place and mtry was 14 and only had 78% accuracy
#without using sex as variable, mtry=1 and accuracy is 72%


# Now begin the full Random Forest analyses: going to use mtry = 5 (default), and mtry = 4 which was given from grid search tuning.

# run a large number of trees with the above mtry values, in order to know what value of ntree achieves convergence of importance values between forests.
# As a starting point, grow 10,000 trees and increase if necessary.

# forests with default mtry:

rf_default_a <- randomForest(x = data[,3:27], y = as.factor(data$FP.status.binary), importance=TRUE ,proximity=TRUE, ntree=10000, strata=as.factor(data$FP.status.binary), sampsize=sample_size)
#save(rf_default_a,file="rf_default_a.Rdata")

rf_default_b <- randomForest(x = data[,3:27], y = as.factor(data$FP.status.binary), importance=TRUE ,proximity=TRUE, ntree=10000, strata=as.factor(data$FP.status.binary), sampsize=sample_size)
#save(rf_default_b,file="rf_default_b.Rdata")

#Check correlation of predictor importance values between forests 
importance_rf_default_a <- data.frame(importance(rf_default_a,type=1)) #type=1 is mean decrease in accuracy for classification, so a large, positive value means that permuting the variable led to a big decrease in prediction accuracy (which is indicative of an important predictor)
colnames(importance_rf_default_a)<-c("importance")

importance_rf_default_b <- data.frame(importance(rf_default_b,type=1))
colnames(importance_rf_default_b)<-c("importance")

cor(importance_rf_default_a,importance_rf_default_b) # A correlation of 0.9960438 for predictor importance values between forests when mtry = 5.099 and ntree = 10,000

# forest with mtry = 4, which was the best mtry value above

rf_4_a <- randomForest(x = data[,2:27], y = as.factor(data$FP.status.binary), importance=TRUE ,proximity=TRUE, mtry=4, ntree=10000, strata=as.factor(data$FP.status.binary), sampsize=sample_size)
#save(rf_17_a,file="rf_17_a.Rdata")

rf_4_b <- randomForest(x = data[,2:27], y = as.factor(data$FP.status.binary), importance=TRUE ,proximity=TRUE, mtry=4, ntree=10000, strata=as.factor(data$FP.status.binary), sampsize=sample_size)
#save(rf_17_b,file="rf_17_b.Rdata")

#Check correlation of locus importance values between forests 
importance_rf_4_a <- data.frame(importance(rf_4_a,type=1)) #type=1 is mean decrease in accuracy for classification, so a large, positive value means that permuting the variable led to a big decrease in prediction accuracy (which is indicative of an important locus)
colnames(importance_rf_4_a)<-c("importance")

importance_rf_4_b <- data.frame(importance(rf_4_b,type=1))
colnames(importance_rf_4_b)<-c("importance")

cor(importance_rf_4_a,importance_rf_4_b) # A correlation of 0.9952817 for predictor importance values between forests when mtry = 4 and ntree = 10,000

# Build final model with mtry = 5.099, ntree = 10,000, as mtry = 5.099 had very slightly higher importance value correlation

# Separate into training and test data (70% and 30%, respectively)
set.seed(3)
train = sample(1:nrow(data), 0.7*nrow(data))
data.test = data[train,]
data.train = data[-train,]

# Create random forest model using training data set
# Make sure to exclude the sample ID column (data.train[,-1])
set.seed(3)
data.model = randomForest(FP.status.binary ~ .,
                          data = data.train[,-1],
                          importance=TRUE,
                          proximity=TRUE,
                          mtry=1, # based on tuning done above
                          ntree=10000, # based on tuning for ntree above
                          strata=as.factor(data$FP.status.binary), sampsize=sample_size) # take into account imbalance in FP postiives vs FP negatives
data.model # OOB-ER = 47.06%
#when running wihtout FFS data the OOB-ER = 10%, mtry of 4, 5
#when running without sex data the OOB-ER =29.41%

# Use the model built above to make predictions on the test data
data.predict = predict(data.model,newdata = data.test[,-1])

# Build a confusion matrix and calculate accuracy of predicted model
table(Prediction = data.predict,Truth = data.test$FP.status.binary)
sum(diag(table(Prediction = data.predict,Truth = data.test$FP.status.binary))) / sum(table(Prediction = data.predict,Truth = data.test$FP.status.binary)) # Accuracy of model with all variables: #81.08%, this becomes 31.8% without FFS data

# Importance Sampling
importance(data.model)
varImpPlot(data.model, main = "Cm FP")

#mean decrease accuracy- higher values indicate importance of that variable in predicting outcome (FP status in this case)
#similar for mean decrease Gini
#based on mean decrease accuracy, location, sample year, sex, then cmhi_15 and supertype 1 were most important
#based on mean decreaseGini samp year, unique vars, supertype1, sex, location were most important

# Let's also build a model based on the mtry=3, since it also did decently in the tuning.

# Separate into training and test data (70% and 30%, respectively)
set.seed(3)
train = sample(1:nrow(data), 0.7*nrow(data))
data.test = data[train,]
data.train = data[-train,]

# Create random forest model using training data set
# Make sure to exclude the sample ID column (data.train[,-1])
set.seed(3)
data.model_3 <- randomForest(FP.status.binary ~ .,
                              data = data.train[,-1],
                              importance=TRUE,
                              proximity=TRUE,
                              ntree=10000, # based on tuning for ntree above
                              strata=as.factor(data$FP.status.binary), sampsize=sample_size) # take into account imbalance in FP postiives vs FP negatives
data.model_3 # OOB-ER = 41.18%

# Use the model built above to make predictions on the test data
data.predict_3 = predict(data.model_3, newdata = data.test[,-1])

# Build a confusion matrix and calculate accuracy of predicted model
table(Prediction = data.predict_3, Truth = data.test$FP.status.binary)
sum(diag(table(Prediction = data.predict_3,Truth = data.test$FP.status.binary))) / sum(table(Prediction = data.predict_3,Truth = data.test$FP.status.binary)) # Accuracy of model with all variables: #81.08%
# accuracy and OOB-ER are comparable

# Importance Sampling
importance(data.model_3)
varImpPlot(data.model_3, main = "Cm FP  mtry3")


#wondering if data is skewed because of FFS data being included? try without. Overall the result of this was lower OOB error at 10% but much lower accuracy at only ~30%. the allele count was the only important thing, but it only had a mean decreaseacc of 1.0

###############################################################################
#try iterative model as described in Han et al 2016 conference paper
data<-data_full %>% select(sample, FP.status.binary, location, cmhi_15, supertype.1, cmhi_11, cmhi_2, cmhi_8, cmhi_13, cmhi_3, total.unique.vars, samp_year)
data<-data_full %>% select(sample, FP.status.binary, location, cmhi_15, cmhi_2, cmhi_8,  samp_year)


# let's also run mtry values of sqrt(p), 2*sqrt(p), 0.1(p), 0.2(p), p/3, in addition to p, where p is the number of variables.

results_optimization <- matrix(data=NA , nrow = 0, ncol = 3) # create matrix that the following loop will dump info into
for (i in seq(from = 100, to = 1000 , by = 100)){  # values of ntree
  print(i)
  for (j in c(2.2, 4.4, 0.5, 1, 1.66, 5)){    #values of mtry based on 30 total predictors
    rf_ij <- randomForest(x = data[,3:7], y = as.factor(data$FP.status.binary), importance=TRUE ,proximity=TRUE, ntree=i, mtry=j, strata=as.factor(data$FP.status.binary), sampsize=sample_size)
    results_optimization <- rbind(results_optimization, c(i,j,tail(rf_ij$err.rate,1)[1]))
  }
}

# Clean up the file format
results_optimization<-as.data.frame(results_optimization)
colnames(results_optimization)<-c("ntree", "mtry","OOB_ER")

# Now plot results to see if there's a plateau
plot(results_optimization$ntree[results_optimization$mtry == 3.1],results_optimization$OOB_ER[results_optimization$mtry == 3.1], type="l", col="black", xlab="ntree",ylab="OOB-ER",ylim=c(0,1))
lines(results_optimization$ntree[results_optimization$mtry == 6.3],results_optimization$OOB_ER[results_optimization$mtry == 6.3], col="blue")
lines(results_optimization$ntree[results_optimization$mtry == 1],results_optimization$OOB_ER[results_optimization$mtry == 1], col="green")
lines(results_optimization$ntree[results_optimization$mtry == 2],results_optimization$OOB_ER[results_optimization$mtry == 2], col="purple")
lines(results_optimization$ntree[results_optimization$mtry == 3.3],results_optimization$OOB_ER[results_optimization$mtry == 3.3], col="orange")
lines(results_optimization$ntree[results_optimization$mtry == 10],results_optimization$OOB_ER[results_optimization$mtry == 10], col="red")

# all mtry parameters seem to behave similarly , but red line is best (value of 25) at minimizing error rate and is overall flattest
#when excluding FFS data, all the lines perform very similarly


# Create model with default paramters as a baseline with grid search, method = "oob"
control <- trainControl(method="oob", number=10, search = "grid")
seed <- 7
metric <- "Accuracy"
set.seed(seed)
mtry <- sqrt(5) # square root of total number of variables; this is the default mtry
tunegrid <- expand.grid(.mtry=mtry)
rf_default <- train(FP.status.binary ~ ., data=data[,2:7], method="rf", metric=metric, tuneGrid=tunegrid, trControl=control)
print(rf_default) # mtry = 5.099 at 83.3% accuracy.

# now let's fine tune
## Grid search, where method = "oob"
control <- trainControl(method="oob", number=10, search="grid")
set.seed(seed)
tunegrid <- expand.grid(.mtry=c(1:5)) # doing full  mtry
rf_gridsearch <- train(FP.status.binary~., data=data[,2:7], method="rf", metric=metric, tuneGrid=tunegrid, trControl=control)
print(rf_gridsearch)
plot(rf_gridsearch) # most accurate: mtry 3 at 87.0.0%

# Separate into training and test data (70% and 30%, respectively)
set.seed(3)
train = sample(1:nrow(data), 0.7*nrow(data))
data.test = data[train,]
data.train = data[-train,]

# Create random forest model using training data set
# Make sure to exclude the sample ID column (data.train[,-1])
set.seed(3)
data.model = randomForest(FP.status.binary ~ .,
                          data = data.train[,-1],
                          importance=TRUE,
                          proximity=TRUE,
                          mtry=1, # based on tuning done above
                          ntree=10000, # based on tuning for ntree above
                          strata=as.factor(data$FP.status.binary), sampsize=sample_size) # take into account imbalance in FP postiives vs FP negatives
data.model # OOB-ER = 47.06%
#when running wihtout FFS data the OOB-ER = 10%, mtry of 4, 5
#when running without sex data the OOB-ER =29.41%

# Use the model built above to make predictions on the test data
data.predict = predict(data.model,newdata = data.test[,-1])

# Build a confusion matrix and calculate accuracy of predicted model
table(Prediction = data.predict,Truth = data.test$FP.status.binary)
sum(diag(table(Prediction = data.predict,Truth = data.test$FP.status.binary))) / sum(table(Prediction = data.predict,Truth = data.test$FP.status.binary)) # Accuracy of model with all variables: #81.08%, this becomes 31.8% without FFS data

# Importance Sampling
importance(data.model)
varImpPlot(data.model, main = "Cm FP") #when using the variables above, the top 5 are cmhi_15, location, cmhi_8, cmhi_2, samp_year

