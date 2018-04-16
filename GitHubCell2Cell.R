setwd("~/Documents/BAS479/FinalCase")

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

# importing data
TRAIN <- read.csv("cell2celltrain.csv")
HOLDOUT <- read.csv("cell2cellholdout.csv")
SS <- read.csv("cell2cellsamplesubmission.csv")


#Check extent of missingness if there's any pattern
MISSING <- data.frame( Column=names(TRAIN),
                       MissingInTRAIN=as.numeric(unlist(lapply(TRAIN,function(x)sum(is.na(x))))),
                       MissingInHOLDOUT=as.numeric(unlist(lapply(HOLDOUT,function(x)sum(is.na(x))))),
                       stringsAsFactors = FALSE)
subset(MISSING,MissingInTRAIN>0 | MissingInHOLDOUT>0 )  # variables that have missing values in both training and holdout
#look at the pattern of missingness
md.pattern(TRAIN[,c("MonthlyRevenue","MonthlyMinutes","TotalRecurringCharge","DirectorAssistedCalls","OverageMinutes",
                    "RoamingCalls","PercChangeMinutes","PercChangeRevenues","AgeHH1","AgeHH2")])  
########################
#Only 3 missing patterns.  
#1)  Call details are either all present or all missing (same with Percentchanges)
#2)  Call details are always missing when percent is missing.
#3)  Ages are always missing or always present, but this is unrelated to other missingness
########################

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



#Missingness of account information is an important predictor of churn.  This is somewhat fishy.  There may be some
#response leakage here?

#Suggesetion:  replace missing values with median and create a new columns for missingness

#Missingness of age is not an important predictor of churn.

#Suggestion:  replace missing values with median

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

# Handling rare levels
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



#################################################################################
#Replace TRAIN and HOLDOUT with imputed and recoded values.
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

CELL <- CELL[c(-27)]   ########

FULLTRAIN <- CELL[1:nrow(TRAIN),]
HOLDOUT <- CELL[-(1:nrow(TRAIN)),]

#split full training data in half
set.seed(479); subtrain.rows <- sample(1:nrow(FULLTRAIN),0.5*nrow(FULLTRAIN))
SUBTRAIN <- FULLTRAIN[subtrain.rows,]
SUBHOLDOUT <- FULLTRAIN[-subtrain.rows,]


############################################################################################################
########################### lurking variable
############################################################################################################
L_TRAIN <- FULLTRAIN

L_TRAIN$Churn <- as.numeric(L_TRAIN$Churn) - 1   # yes=1, no=0
aggregate(Churn~IncomeGroup,data=L_TRAIN,FUN=mean)
AOV <- aov(Churn~IncomeGroup,data=L_TRAIN)
summary(AOV)
TUKEY <- TukeyHSD(AOV)
TUKEY
library(multcompView)
multcompLetters4(AOV,TUKEY)
# IncomeGroup     Churn
# 1            0 0.3022984
# 2            1 0.2913193
# 3            2 0.2963272
# 4            3 0.2631227
# 5            4 0.2765852
# 6            5 0.2827311
# 7            6 0.2835432
# 8            7 0.2979411
# 9            8 0.2765065
# 10           9 0.2820421
# $IncomeGroup
# 0    7    2    1    6    5    9    4    8    3 
# "a"  "a" "ab" "ab" "ab" "ab" "ab" "ab" "ab"  "b" 
cluster <- makeCluster(detectCores() - 1) # convention to leave 1 core for OS
registerDoParallel(cluster)
RF <- randomForest(Churn~.-IncomeGroup,data=L_TRAIN)
stopCluster(cluster)
registerDoSEQ()
r <- L_TRAIN$Churn - predict(RF)
aggregate(r~L_TRAIN$IncomeGroup,FUN=mean)
AOV <- aov(r~CUSTLOYALTY$Income)
summary(AOV)
TUKEY <- TukeyHSD(AOV)
multcompLetters4(AOV,TUKEY)  


aggregate(Churn~ChildrenInHH,data=FULLTRAIN,FUN=mean)
AOV <- aov(Churn~ChildrenInHH,data=FULLTRAIN)
summary(AOV)

aggregate(Churn~HandsetRefurbished,data=FULLTRAIN,FUN=mean)
AOV <- aov(Churn~HandsetRefurbished,data=FULLTRAIN)
summary(AOV)
RF <- randomForest(Churn~.-HandsetRefurbished,data=FULLTRAIN)

r <- CELL$Churn - predict(RF)
t.test(r~CELL$Churn)

aggregate(Churn~HandsetWebCapable,data=CELL,FUN=mean)
AOV <- aov(Churn~HandsetWebCapable,data=CELL)
summary(AOV)

aggregate(Churn~TruckOwner,data=CELL,FUN=mean)
AOV <- aov(Churn~TruckOwner,data=CELL)
summary(AOV)

RF <- randomForest(Churn~.-ChildrenInHH,data=CELL)
r <- CELL$Churn - predict(RF)
t.test(r~CELL$Churn)



DOG$OutcomeType <- as.numeric(DOG$OutcomeType) - 1 # Eu=1, adopt=0
aggregate(OutcomeType~Breed, data=DOG, FUN=mean)
AOV <- aov(OutcomeType~Breed,data=DOG)
summary(AOV)



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

RPARTfit  #Look at output to see which cp was best
RPARTfit$bestTune #Gives best parameters
RPARTfit$results #Look at output in more detail (lets you see SDs)
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
FORESTfit.new  #Look at output to see which cp was best
FORESTfit.new$bestTune #Gives best parameters
FORESTfit.new$results #Look at output in more detail (lets you see SDs)
FORESTfit.new$results[rownames(FORESTfit.new$bestTune),] 

BEST$FORESTfit.new<-FORESTfit.new$bestTune
varImp(FORESTfit.new)

FORESTfit.new
FORESTfit.new.predictions <- predict(FORESTfit.new,newdata=HOLDOUT,type="prob")
head(FORESTfit.new.predictions)

FORESTfit.new.predictions$CustomerID <- HOLDOUT$CustomerID
FORESTfit.new.predictions$Churn <- FORESTfit.new.predictions$Yes
FORESTfit.new.predictions <- FORESTfit.new.predictions[-c(1:2)]
head(FORESTfit.new.predictions)
dim(FORESTfit.new.predictions)


write.csv(FORESTfit.new.predictions,file="SineenartFinalPrediction.csv",row.names = FALSE)










##########################################################################
#                               BOOSTED TREE 
##########################################################################
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
GBMfit
GBMfit$bestTune
GBMfit$results
BEST$GBMfit<-GBMfit$bestTune
GBMfit$results[rownames(GBMfit$bestTune),]   #.68

#     shrinkage interaction.depth n.minobsinnode n.trees       ROC      Sens      Spec      ROCSD      SensSD     SpecSD
# 36      0.01                 4             10    2000 0.6842619 0.9680758 0.1170632 0.00715037 0.002630766 0.00130764

SUBTRAIN$CurrentEquipmentDays
varImp(GBMfit)

MadeCallToRetentionTeamYes  11.623
MissingPercentMinYes 
Regiond
SUBTRAIN$CurrentEquipmentDays

GBMfit
GBMfit.predictions <- predict(GBMfit,newdata=HOLDOUT,type="prob")
head(GBMfit.predictions)

GBMfit.predictions$CustomerID <- HOLDOUT$CustomerID
GBMfit.predictions$Churn <- GBMfit.predictions$Yes
GBMfit.predictions <- GBMfit.predictions[-c(1:2)]
head(GBMfit.predictions)
dim(GBMfit.predictions)


write.csv(GBMfit.predictions,file="SineenartFinalPredictionGMBfit.csv",row.names = FALSE)


#####################################
#K-nearest neighbors
#####################################
# 
# #convert_categorical will take a data frame and replace categorical variables with sets of indicator variables
# #Need to specify which variable is the y variable.  
# convert_categorical <- function(DATA,y=NA,denominator=1) {
#   if(is.na(y)) { stop("Must specify a y variable in quotes") }
#   to.convert <- setdiff( which(unlist(lapply(DATA,class))=="factor"), which(names(DATA)==y) )  #column numbers of categorical variables
#   if(length(to.convert)==0) { return(DATA) } #nothing to be done
#   for (i in to.convert) {  #loop over non-numerical columns
#     formula <- as.formula( paste("~",names(DATA)[i],"-1",sep="") )  #set up formula for model.matrix
#     INDICATORS <- data.frame(model.matrix(formula,data=DATA))/denominator
#     DATA <- cbind(DATA,INDICATORS) #add columns for the indicator variables
#   }
#   DATA <- DATA[,-to.convert[-1]] #Get rid of categorical variables
# }
# #Tune values for number of neighbors
# paramGrid <- expand.grid(k=3:5)   
# set.seed(seed)
# TRAIN.NN <- TRAIN
# TRAIN.NN <- convert_categorical(TRAIN.NN,y="Churn")
# TRAIN.NN$Category <- NULL
# cluster <- makeCluster(detectCores() - 1) # convention to leave 1 core for OS
# registerDoParallel(cluster)
# #AUC
# KNN <- train(Churn~.,data=TRAIN, metric="ROC", method='knn', trControl=fitControl,tuneGrid=paramGrid)
# stopCluster(cluster)
# registerDoSEQ()
# KNN
# BEST$knn<-KNN$bestTune






MASTERRESULTS <- rbind(RPARTfit$results[,2:7],
                       FORESTfit.new$results[,2:7],
                       GBMfit$results[,5:10])   ##
# KNN$results[,2:7])  ##

MASTERRESULTS$method <- c(rep("rpart",nrow(RPARTfit$results)),
                          rep("randomforest",nrow(FORESTfit.new$results)),
                          rep("gbm",nrow(GBMfit$results)))
# rep("knn",nrow(KNN$results)))

MASTERRESULTS$param <- c(RPARTfit$results$cp,
                         FORESTfit.new$results$mtry,
                         paste(GBMfit$results$n.trees,GBMfit$results$shrinkage,GBMfit$results$interaction.depth,GBMfit$results$n.minobsinnode,sep=";"))
# KNN$results$k)
# 
# ##Look at results        #####################
# head( MASTERRESULTS[order(MASTERRESULTS$Accuracy,decreasing=TRUE),], 20 )
# 
head( MASTERRESULTS[order(MASTERRESULTS$ROC,decreasing=TRUE),], 20 )
# 
# head( MASTERRESULTS[order(MASTERRESULTS$ROC,decreasing=TRUE),], 20 ) 

holdout.classifications <- predict(FORESTfit.new,newdata=HOLDOUT)  #classes
head(HOLDOUT)
mean(holdout.classifications==HOLDOUT$Churn)  #accuracy on holdout
head(HOLDOUT)
holdout.predictions <- predict(FORESTfit.new,newdata=HOLDOUT,type="prob")  #probabilities   !!!!!!!!!!!

library(pROC)
head(holdout.predictions)
head(holdout.predictions[,2])
roc(HOLDOUT$Churn,holdout.classifications[,2])
levels(holdout.classifications)
# roc(HOLDOUT$Junk,holdout.predictions[,2]) #will always be [,2] since that is class that's last alphabetically
# 
# 
# 
# class(TRAIN$CallForwardingCalls)
# class(HOLDOUT$CallForwardingCalls)
# 
# class(TRAIN$ReferralsMadeBySubscriber)
# class(HOLDOUT$ReferralsMadeBySubscriber)
# 
# # Error: variables ‘CallForwardingCalls’, 
# # ‘RetentionCalls’, ‘RetentionOffersAccepted’, ‘ReferralsMadeBySubscriber’, 
# # ‘IncomeGroup’, ‘AdjustmentsToCreditRating’ were specified with different types from the fit
# 
class(TRAIN$RetentionCalls)
as.factor(HOLDOUT$RetentionCalls)
# class(TRAIN$RetentionOffersAccepted)
# as.integer(HOLDOUT$RetentionOffersAccepted)
class(TRAIN$ReferralsMadeBySubscriber)
as.factor(HOLDOUT$ReferralsMadeBySubscriber)
# class(TRAIN$IncomeGroup)
# as.integer(HOLDOUT$IncomeGroup)
class(TRAIN$AdjustmentsToCreditRating)
as.factor(HOLDOUT$AdjustmentsToCreditRating)
# class(TRAIN$RetentionCalls)
# class(HOLDOUT$RetentionCalls)



# Error: variables ‘RetentionCalls’, 
# ‘ReferralsMadeBySubscriber’, ‘AdjustmentsToCreditRating’ 
# were specified with different types from the fit
class(TRAIN$RetentionCalls)
HOLDOUT$RetentionCalls <- as.factor(HOLDOUT$RetentionCalls)
class(HOLDOUT$RetentionCalls)

class(TRAIN$ReferralsMadeBySubscriber)
HOLDOUT$ReferralsMadeBySubscriber <- as.factor(HOLDOUT$ReferralsMadeBySubscriber)

class(TRAIN$AdjustmentsToCreditRating)
HOLDOUT$AdjustmentsToCreditRating <- as.factor(HOLDOUT$AdjustmentsToCreditRating)


GBMfit$results[rownames(GBMfit$bestTune),] 

holdout.predictions <- predict(FORESTfit.new,newdata=SS,type="prob")  


GLMnet.predictions <- predict(GLMnet,newdata=SS,type="prob")
dim(GLMnet.predictions)
holdout.predictions <- holdout.predictions[,2]
write.csv(holdout.predictions,file="SineenartFinalPrediction.csv",row.names = TRUE)
?write.csv


NNET$bestTune
NNET$results[rownames(NNET$bestTune),] #0.9344 for me
NNET.predictions <- predict(NNET,newdata=HOLDOUT.TO.SHARE,type="prob")
dim(NNET.predictions)  #verify 4986 rows
write.csv(NNET.predictions,file="SineenartFinalPredictions.csv",row.names = FALSE)


#### XG Boost
fitControl <- trainControl(method="cv",number=5, verboseIter = TRUE,
                           summaryFunction = twoClassSummary, classProbs = TRUE, allowParallel = TRUE, returnData = FALSE) 
?trainControl
#fitControl <- trainControl(method = "cv", number = 5, classProbs = TRUE, verboseIter = TRUE)
xgboostGrid <- expand.grid(eta=0.1,nrounds=c(100,200,500,1000,1500),
                           max_depth=5,min_child_weight=1,gamma=0,colsample_bytree=0.8,subsample=0.8)
TRAINXGB <- SUBTRAIN
TRAINXGB$Churn <- as.numeric(TRAINXGB$Churn)-1
TRAINXGB <- sparse.model.matrix(Churn~.-1,data=TRAINXGB)
cluster <- makeCluster(detectCores() - 1) # convention to leave 1 core for OS
registerDoParallel(cluster)
XTREME <- train(x=TRAINXGB,y=SUBTRAIN$Churn,
                method="xgbTree",
                trControl=fitControl,
                tuneGrid=xgboostGrid,
                verbose=FALSE)

stopCluster(cluster)
registerDoSEQ()
XTREME$results[,c(7,8,10)]


xgboostGrid <- expand.grid(eta=0.1,nrounds=500,
                           max_depth=c(3,6,9),min_child_weight=1,gamma=c(.1,.5,1,5,10),
                           colsample_bytree=c(0.6,0.8,1),subsample=c(0.6,0.8,1.0))
dim(xgboostGrid) #135 models yikes
cluster <- makeCluster(detectCores() - 1) # convention to leave 1 core for OS
registerDoParallel(cluster)
XTREME <- train(x=TRAINXGB,y=SUBTRAIN$Churn,
                method="xgbTree",
                trControl=fitControl,
                tuneGrid=xgboostGrid,
                verbose=TRUE)
stopCluster(cluster)
registerDoSEQ()

TUNE <- XTREME$results[,c(2,3,4,6,8,10)]
head( TUNE[order(TUNE$ROC,decreasing=TRUE),], 15 )



xgboostGrid <- expand.grid(eta=0.01,nrounds=c(250,500,1000,2500,5000),
                           max_depth=9,min_child_weight=1,gamma=1,colsample_bytree=1,subsample=.6)
FINE.XTREME <- train(x=TRAINXGB,y=SUBTRAIN$Churn,
                     method="xgbTree",
                     trControl=fitControl,
                     tuneGrid=xgboostGrid,
                     verbose=FALSE)





predict <- predict(XTREME,newdata=SUBHOLDOUT,type=)
xgboost.predictions <- predict(XTREME,newdata=HOLDOUT,type="prob")

GBMfit.predictions <- predict(GBMfit,newdata=HOLDOUT,type="prob")
dim(xgboost.predictions)
xgboost.predictions$CustomerID <- HOLDOUT$CustomerID
xgboost.predictions$Churn <- xgboost.predictions$Yes
xgboost.predictions <- xgboost.predictions[-c(1:2)]
head(xgboost.predictions)
dim(xgboost.predictions)


write.csv(xgboost.predictions,file="SineenartFinalPredictionXGBoost.csv",row.names = FALSE)

