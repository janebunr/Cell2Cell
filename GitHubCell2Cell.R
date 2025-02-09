library(parallel)
library(iterators)
library(doParallel)
library(regclass)
library(mice)
library(discretization)
library(regclass)
library(caret)
library(caretEnsemble)
library(nnet)
library(neuralnet)
library(e1071)
library(glmnet)
library(gbm)
library(pROC)

# function to combine infrequent levels in a variable for better fit in a model
combine_infrequent_levels <- function(x,threshold=20,newname="Combined") { 
  x <- factor(x)
  rare.levels <- names( which( sort( table(x) ) <= threshold ) )
  if(length(rare.levels)==0) { return(list(values=x,combined=NULL)) }
  
  levels(x)[ which(levels(x) %in% rare.levels) ] <- newname
  ST <- sort(table(x))
  if(ST[newname]<=threshold) {  #Combined will be the least frequent level
    levels.to.combine <- which( levels(x) %in% c(newname,names(ST)[2]))
    levels(x)[levels.to.combine] <- newname
    rare.levels <- c(rare.levels,names(ST)[2])}
  return(list(values=x,combined=rare.levels))
}


discretize_x_for_categorical_y <- function(DATA,threshold=0,train.rows=NA,equal=FALSE) {
  require(discretization)
  require(regclass)
  if(class(DATA[,1]) %in% c("character","factor")) {
    if(threshold == 0 | is.na(threshold)) { } else { 
      DATA[,1] <- combine_infrequent_levels(DATA[,1],threshold)$values }
    if( length(train.rows)>1 ) { HALFDISC <- DATA[train.rows,] } else { HALFDISC <- DATA }
    A <- aggregate(HALFDISC[,2]~HALFDISC[,1],FUN=function(x)mean(x==levels(x)[1])) 
    A <- A[order(A[,2]),]
    SUB <- HALFDISC
    SUB[,1] <- match(HALFDISC[,1],A[,1]) 
    disc.scheme <- mdlp(SUB)
    cutoffs <- sort( unlist( disc.scheme$cutp ) )
    A$value <- 1:nrow(A)
    if(cutoffs[1] != "All") {
      A$code <- factor( rep(letters[1:(length(cutoffs)+1)],nrow(A)) )[1:nrow(A)] 
      
      for (i in 1:length(cutoffs)) {
        if(i==1) { A$code[1:floor(cutoffs[1])] <- letters[1] } else { 
          A$code[ which(A$value > cutoffs[i-1] & A$value <= cutoffs[i]) ] <- letters[i] }
      }
      A$code[ which(A$value>max(cutoffs)) ] <- letters[i+1]
      names(A) <- c( "OldValue","yAverage","Rank","NewValue")
      results <- list(Details=A,newlevels= factor( A$NewValue[ match(DATA[,1],A[,1]) ] ) )
      y <- DATA[,2]
      x <- results$newlevels
      TEMP <- data.frame(y,x)
      
      mosaic(y~x,data=TEMP,xlab="New Levels",ylab=names(DATA)[2],inside=TRUE,equal=equal) } else {
        names(A) <- c( "OldValue","yAverage","Rank")
        A$Rank <- order(A$yAverage)
        A$NewValue <- factor(rep("A",nrow(A)))
        results <- list(Details=A,newlevels= factor( rep("A",nrow(A)) ) )
      }
    
    return( results ) } else {
      
      
      
      if( length(train.rows)>1) { HALFDISC <- DATA[train.rows,] } else { HALFDISC <- DATA }
      SUB <- HALFDISC
      disc.scheme <- mdlp(SUB)
      cuts <- unlist( disc.scheme$cutp )
      A <- aggregate(disc.scheme$Disc.data[,2]~disc.scheme$Disc.data[,1],FUN=function(x)mean(x==levels(x)[1]))
      details <- data.frame(Lower=rep(0,nrow(A)),Upper=rep(0,nrow(A)),yAverage=rep(0,nrow(A)),NewValue=rep(0,nrow(A)))
      
      for (i in 1:nrow(details)) {
        details$Lower[i] <- min( SUB[ which(disc.scheme$Disc.data[,1]==i),1] )
        details$Upper[i] <- max( SUB[ which(disc.scheme$Disc.data[,1]==i),1] )
        details$NewValue[i] <- letters[i]
      }
      details$yAverage <- A[,2]
      newvalues <- rep(0,nrow(DATA)) 
      for (i in 1:nrow(DATA)) {
        temp <- which( DATA[i,1] >= details$Lower )
        if(length(temp)>0) { temp <- max( which( DATA[i,1] >= details$Lower ) ) } else { temp <- 1 }
        newvalues[i] <- details$NewValue[temp]
      }
      if(cuts[1]!="All") { 
        results <- list(Details=details,newlevels=factor(newvalues))
        y <- DATA[,2]
        x <- results$newlevels
        TEMP <- data.frame(y,x)
        mosaic(y~x,data=TEMP,xlab="New Levels",ylab=names(DATA)[2],inside=TRUE,equal=equal) } else {
          details$NewValue <- "A"
          results <- list(Details=details,newlevels=factor(rep("A",nrow(DATA))))
        }
      return(results) }
}



#################
# importing data
#################
TRAIN <- read.csv("cell2celltrain.csv")
HOLDOUT <- read.csv("cell2cellholdout.csv")
SS <- read.csv("cell2cellsamplesubmission.csv")

###################################################
#Check extent of missingness if there's any pattern
###################################################
MISSING <- data.frame( Column=names(TRAIN),
                       MissingInTRAIN=as.numeric(unlist(lapply(TRAIN,function(x)sum(is.na(x))))),
                       MissingInHOLDOUT=as.numeric(unlist(lapply(HOLDOUT,function(x)sum(is.na(x))))),
                       stringsAsFactors = FALSE)
subset(MISSING,MissingInTRAIN>0 | MissingInHOLDOUT>0 )  # variables that have missing values in both training and holdout
#look at the pattern of missingness
md.pattern(TRAIN[,c("MonthlyRevenue","MonthlyMinutes","TotalRecurringCharge","DirectorAssistedCalls","OverageMinutes",
                    "RoamingCalls","PercChangeMinutes","PercChangeRevenues","AgeHH1","AgeHH2")])  


#Looking if missingness itself a predictor of churn?
# MonthlyMinutes
fisher.test(TRAIN$Churn,is.na(TRAIN$MonthlyMinutes))  #significant
mosaic(TRAIN$Churn~is.na(TRAIN$MonthlyMinutes) ,equal=TRUE)

# PercChangeRevenues
fisher.test(TRAIN$Churn,is.na(TRAIN$PercChangeRevenues))  #highly significant
mosaic(TRAIN$Churn~is.na(TRAIN$PercChangeRevenues) ,equal=TRUE)

#TotalRecurringCharge
fisher.test(TRAIN$Churn,is.na(TRAIN$TotalRecurringCharge))  #highly significant
mosaic(TRAIN$Churn~is.na(TRAIN$TotalRecurringCharge) ,equal=TRUE)

#DirectorAssistedCalls
fisher.test(TRAIN$Churn,is.na(TRAIN$DirectorAssistedCalls))  #highly significant
mosaic(TRAIN$Churn~is.na(TRAIN$DirectorAssistedCalls) ,equal=TRUE)

#OverageMinutes
fisher.test(TRAIN$Churn,is.na(TRAIN$OverageMinutes))  #highly significant
mosaic(TRAIN$Churn~is.na(TRAIN$OverageMinutes) ,equal=TRUE)

#RoamingCalls
fisher.test(TRAIN$Churn,is.na(TRAIN$RoamingCalls))  #highly significant
mosaic(TRAIN$Churn~is.na(TRAIN$RoamingCalls) ,equal=TRUE)

#PercChangeMinutes
fisher.test(TRAIN$Churn,is.na(TRAIN$PercChangeMinutes))  #highly significant
mosaic(TRAIN$Churn~is.na(TRAIN$PercChangeMinutes) ,equal=TRUE)

#PercChangeRevenues
fisher.test(TRAIN$Churn,is.na(TRAIN$PercChangeRevenues))  #highly significant
mosaic(TRAIN$Churn~is.na(TRAIN$PercChangeRevenues) ,equal=TRUE)

#AgeHH2
fisher.test(TRAIN$Churn,is.na(TRAIN$AgeHH2))  # NOT significant
mosaic(TRAIN$Churn~is.na(TRAIN$AgeHH2) ,equal=TRUE)

#AgeHH1
fisher.test(TRAIN$Churn,is.na(TRAIN$AgeHH1))  #NOT significant
mosaic(TRAIN$Churn~is.na(TRAIN$AgeHH1) ,equal=TRUE)


# Because missingness is an important factor in predicting churn. 
# We'll create new columns for missingness and replace missing values with its median for each variable
# Process data together
CELL <- rbind(TRAIN,HOLDOUT)
CELL$MissingMonthlyRev <- factor( ifelse(is.na(CELL$MonthlyRevenue),"Yes","No") )
CELL$MissingPercentMin <- factor( ifelse(is.na(CELL$PercChangeMinutes),"Yes","No") )
CELL$MissingMonthlyMin <- factor( ifelse(is.na(CELL$MonthlyMinutes),"Yes","No") )
CELL$MissingTotalRecur <- factor( ifelse(is.na(CELL$TotalRecurringCharge),"Yes","No") )
CELL$MissingDirecAssis <- factor( ifelse(is.na(CELL$DirectorAssistedCalls),"Yes","No") )
CELL$MissingOverageMin <- factor( ifelse(is.na(CELL$OverageMinutes),"Yes","No") )
CELL$MissingRoamingCalls <- factor( ifelse(is.na(CELL$RoamingCalls),"Yes","No") )
CELL$MissingPercChangeRev <- factor( ifelse(is.na(CELL$PercChangeRevenues),"Yes","No") )

for( columns in MISSING$Column[which(MISSING$MissingInTRAIN>0)] ) {
  missing.positions <- which( is.na(CELL[,columns]) )
  replacement.value <- median(CELL[-missing.positions,columns])
  CELL[missing.positions,columns] <- replacement.value
}


#######################################
# Handling rare levels#################
#######################################
service <- as.character(CELL$ServiceArea)
head(service) # [1] "SEAPOR503" "PITHOM412" "MILMIL414" "PITHOM412" "OKCTUL918" "OKCTUL918"

# Extracting area code
area.code <- substr(service,start=7,stop=9)
area.code <- combine_infrequent_levels(area.code,threshold=20)$values
area.code <- factor( area.code )
y <- CELL$Churn
D <- data.frame(area.code,y)
DC <- discretize_x_for_categorical_y(D,train.rows=1:nrow(TRAIN) )
CELL$AreaCode <- factor( DC$newlevels )

# Extracting region code
region.code <- as.factor(substr(service,start=1,stop=6))
region.code <- combine_infrequent_levels(region.code,threshold=20)$values
region.code <- factor( region.code )
y <- CELL$Churn
DD <- data.frame(region.code,y)
DDCC <- discretize_x_for_categorical_y(DD,train.rows = 1:nrow(TRAIN))
CELL$Region <- factor( DDCC$newlevels )
summary(DD)
# region.code       y        
# Combined: 2759   No  :36336  
# NYCBRO  : 2354   Yes :14711  
# HOUHOU  : 2111   NA's:20000  
# DALDAL  : 2040               
# NYCMAN  : 1648               
# BOSBOS  : 1410               
# (Other) :58725  

summary(CELL$Region)
# a     b     c     d 
# 4494 23080 33075 10398 


###########################################################
#Replace TRAIN and HOLDOUT with imputed and recoded values.
##########################################################
TRAIN <- CELL[1:nrow(TRAIN),]
HOLDOUT <- CELL[-(1:nrow(TRAIN)),]
CELL[,3] <- log10( CELL[,3] - min(CELL[,3]) + 1 )
CELL[,4] <- log10( CELL[,4] - min(CELL[,4]) + 1 )
CELL[,5] <- log10( CELL[,5] - min(CELL[,5]) + 1 )
CELL[,7] <- log10( CELL[,7] - min(CELL[,7]) + 1 )
CELL[,8] <- log10(.1 + CELL[,8])
CELL[,11] <- log10(CELL[,11]+.5)
CELL[,12] <- log10(CELL[,12]+.5)
CELL[,13] <- log10(CELL[,13]+.1)
CELL[,14] <- log10(CELL[,14]+.1)
CELL[,15] <- log10(CELL[,15]+.1)
CELL[,16] <- log10(CELL[,16]+.01)
CELL[,17] <- log10(CELL[,17]+.5)
CELL[,18] <- log10(CELL[,18]+.5)
CELL[,19] <- log10(CELL[,19]+1)
CELL[,20] <- log10(CELL[,20]+1)
CELL[,21] <- log10(CELL[,21]+1)
CELL[,22] <- factor( ifelse(CELL[,22]>0,"Yes","No") )
CELL[,23] <- log10(CELL[,23]+1)
CELL[,25] <- log10(CELL[,25]+1)
CELL[,26] <- log10(CELL[,26]+1)
CELL[,28] <- log10(CELL[,28])
CELL[,29] <- log10(CELL[,29])
CELL[,30] <- log10(5+CELL[,30]-min(CELL[,30]))

#31 and 32 are ages, and a 0 is recorded when the age is not known it appears
#Good chance for discretization
CELL[,45] <- factor(CELL[,45])
levels(CELL[,45]) <- c("None","One",rep("TwoOrMore",3))
CELL[,46] <- factor(CELL[,46])
levels(CELL[,46]) <- c("None","One",rep("TwoOrMore",3))
CELL[,49] <- factor(CELL[,49])
levels(CELL[,49]) <- c("None","One","Two",rep("ThreeOrMore",10))
CELL[,50] <- factor(CELL[,50],ordered = TRUE)
CELL[,52] <- factor(CELL[,52])
levels(CELL[,52]) <- c("None","One","Two",rep("ThreeOrMore",nlevels(CELL[,52])-3))
levels(CELL[,53])
median( as.numeric(as.character(CELL[-which(CELL[,53]=="Unknown"),53] )) )
CELL[,53] <- factor(CELL[,53],ordered=TRUE,levels=c("10","30","40","60","Unknown","80","100","130","150","180","200","240","250","300","400","500"))
CELL[,55] <- factor(CELL[,55],ordered=TRUE)
CELL <- CELL[c(-27)]
FULLTRAIN <- CELL[1:nrow(TRAIN),]
HOLDOUT <- CELL[-(1:nrow(TRAIN)),]

#split full training data in half
set.seed(479); subtrain.rows <- sample(1:nrow(FULLTRAIN),0.5*nrow(FULLTRAIN))
SUBTRAIN <- FULLTRAIN[subtrain.rows,]
SUBHOLDOUT <- FULLTRAIN[-subtrain.rows,]


######################################### PREDICTIVE MODEL ######################################### 
######################################### ######################################### ######################################### 
######################################### ######################################### ######################################### 
BEST <- list()
seed <- 111
set.seed(seed)

####################################
# Partition
####################################
rpartGrid <- expand.grid(cp=10^seq(from=-4,to=-1,length=30))
#train function
cluster <- makeCluster(detectCores() - 1) # convention to leave 1 core for OS
registerDoParallel(cluster)
fitControl <- trainControl(method="cv",number=5, verboseIter = TRUE,
                           summaryFunction = twoClassSummary, classProbs = TRUE, allowParallel = TRUE)
set.seed(seed); RPARTfit <- train(Churn~.,data=TRAIN,method="rpart",metric="ROC",trControl=fitControl,tuneGrid=rpartGrid)
stopCluster(cluster)
registerDoSEQ()

RPARTfit$results[rownames(RPARTfit$bestTune),] #0.8749798  #not so good   #########
#              cp       ROC      Sens     Spec       ROCSD      SensSD      SpecSD
# 8 0.0005298317 0.6266175 0.9621035 0.109306 0.006490512 0.006103816 0.005945355
plot(ROC~cp,data=RPARTfit$results,log="x")  #If tuned on AUC
BEST$RPARTfit<-RPARTfit$bestTune

RPART.predictions <- predict(RPARTfit,newdata=HOLDOUT,type="prob")


#####################################
# Random Forest
#####################################
forestGrid <- expand.grid(mtry=c(1,2,3,7))
cluster <- makeCluster(detectCores() - 1) # convention to leave 1 core for OS
registerDoParallel(cluster)
fitControl <- trainControl(method="cv",number=5, verboseIter = TRUE,
                           summaryFunction = twoClassSummary, classProbs = TRUE, allowParallel = TRUE)
set.seed(seed); FORESTfit.new <- train(Churn~.,data=SUBTRAIN,method="rf",metric="ROC",trControl=fitControl,tuneGrid=forestGrid)

stopCluster(cluster)
registerDoSEQ()

FORESTfit.new$results[rownames(FORESTfit.new$bestTune),] 
#   mtry       ROC      Sens       Spec       ROCSD    SensSD      SpecSD
# 4    7 0.6640235 0.9888815 0.04867437 0.008280076 0.0018111 0.004939742



#########################
# BOOSTED TREE 
#########################
gbmGrid <- expand.grid(n.trees=c(500,1000,2000),
                       interaction.depth=2:4,
                       shrinkage=c(.005,.01),
                       n.minobsinnode=c(5,10))
cluster <- makeCluster(detectCores() - 1) # convention to leave 1 core for OS
registerDoParallel(cluster)
fitControl <- trainControl(method="cv",number=5, verboseIter = TRUE,
                           summaryFunction = twoClassSummary, classProbs = TRUE, allowParallel = TRUE) 
set.seed(474); GBMfit <- train(Churn~.,data=SUBTRAIN,method="gbm",metric="ROC",trControl=fitControl,tuneGrid=gbmGrid,verbose=FALSE)
stopCluster(cluster)
registerDoSEQ()

GBMfit$results[rownames(GBMfit$bestTune),]   #.68
#     shrinkage interaction.depth n.minobsinnode n.trees       ROC      Sens      Spec      ROCSD      SensSD     SpecSD
# 36      0.01                 4             10    2000 0.6842619 0.9680758 0.1170632 0.00715037 0.002630766 0.00130764


################# Finding the best model #################################
MASTERRESULTS <- rbind(RPARTfit$results[,2:7],
                       FORESTfit.new$results[,2:7],
                       GBMfit$results[,5:10])

MASTERRESULTS$method <- c(rep("rpart",nrow(RPARTfit$results)),
                          rep("randomforest",nrow(FORESTfit.new$results)),
                          rep("gbm",nrow(GBMfit$results)))

MASTERRESULTS$param <- c(RPARTfit$results$cp,
                         FORESTfit.new$results$mtry,
                         paste(GBMfit$results$n.trees,GBMfit$results$shrinkage,GBMfit$results$interaction.depth,GBMfit$results$n.minobsinnode,sep=";"))

head( MASTERRESULTS[order(MASTERRESULTS$ROC,decreasing=TRUE),], 5 )
#           ROC      Sens       Spec       ROCSD      SensSD      SpecSD method           param
# 36  0.6842619 0.9680758 0.11706322 0.007150370 0.002630766 0.001307640    gbm  2000;0.01;4;10   # best model
# 331 0.6837387 0.9686813 0.11733515 0.006207841 0.002683096 0.003742042    gbm   2000;0.01;4;5
# 271 0.6819408 0.9729746 0.10129164 0.006915591 0.002965540 0.002676413    gbm   2000;0.01;3;5
# 301 0.6819376 0.9731947 0.10237933 0.006795763 0.003684426 0.002986524    gbm  2000;0.01;3;10
# 181 0.6785689 0.9823316 0.07763426 0.007657496 0.002081257 0.002819367    gbm 2000;0.005;4;10

varImp(GBMfit)
# only 20 most important variables shown (out of 112)
# 
# Overall
# CurrentEquipmentDays       100.000
# MonthsInService             69.802
# MonthlyMinutes              61.336
# PercChangeMinutes           54.557
# PercChangeRevenues          24.755
# TotalRecurringCharge        22.134
# CustomerID                  20.389
# UniqueSubs                  19.660
# OverageMinutes              19.048
# AgeHH1                      15.182
# Regiond                     13.589
# MonthlyRevenue              12.495
# MadeCallToRetentionTeamYes  11.623
# MissingPercentMinYes        11.579
# DroppedCalls                10.424
# ReceivedCalls                9.346
# PeakCallsInOut               9.096
# CustomerCareCalls            9.037
# DroppedBlockedCalls          8.289
# OffPeakCallsInOut            8.083



########### Model Accuracy ########

GBMfit.predictions <- predict(GBMfit,newdata=SUBHOLDOUT,type="prob")

GBMfit.predictions$CustomerID <- SUBHOLDOUT$CustomerID
GBMfit.predictions$Churn <- GBMfit.predictions$Yes
GBMfit.predictions <- GBMfit.predictions[-c(1:2)]
results <- data.frame(pred=ifelse(GBMfit.predictions$Yes > .5, "Yes", "No"),actual=SUBHOLDOUT$Churn)
confusionMatrix(results$pred, results$actual)
# Confusion Matrix and Statistics
# 
# Reference
# Prediction    No   Yes
# No          17586  6521
# Yes          582   835
# 
# Accuracy : 0.7217          ##
# 95% CI : (0.7162, 0.7272)
# No Information Rate : 0.7118          
# P-Value [Acc > NIR] : 0.0002321       
# 
# Kappa : 0.1072          
# Mcnemar's Test P-Value : < 2.2e-16       
# 
# Sensitivity : 0.9680          
# Specificity : 0.1135          
# Pos Pred Value : 0.7295          
# Neg Pred Value : 0.5893          
# Prevalence : 0.7118          
# Detection Rate : 0.6890          
# Detection Prevalence : 0.9445          
# Balanced Accuracy : 0.5407          
# 
# 'Positive' Class : No   











