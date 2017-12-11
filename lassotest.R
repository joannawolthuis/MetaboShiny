library(caret)
library(data.table)
library(xlsx)
# --------------

mat <- fread(input="/Users/jwolthuis/Google Drive/MetaboShiny/backend/appdata/brazil_chicken/BrazilFirstEight.csv",header = T)

facsTable <- read.xlsx("~/Documents/umc/DSM_BR poultry worksheets/DSM_BR1_8.xlsx",sheetName = "Individual Data",header = T)

pastytable <- as.data.table(facsTable[1:80,-c(1,3,4)])
setkey(pastytable, "Card.ID")
# reorder according
endorder <- mat$Sample
pastytable <- pastytable[endorder]
pastytable <- apply(pastytable, MARGIN = 2, FUN=function(vec){
 as.numeric(as.factor(vec))
})

# add to mat

mat <- cbind(mat, pastytable[,-1])

# ==================================

mat$Label <- as.factor(mat$Label)

inTrain <- createDataPartition(y = mat$Label,
                               ## the outcome data are needed
                                  p = .6,
                                ## The percentage of data in the training set
                                 list = FALSE)
training <- mat[ inTrain,]
testing <- mat[-inTrain,]

# prepare training scheme
control <- trainControl(method="cv", number=2, returnResamp="all",
                       classProbs=TRUE, summaryFunction=twoClassSummary)

colnames(training)[1:10]
# train the model
set.seed(849)

model <- train(x = training[,-c(1,2)], 
               y = training$Label,
               method = "rf",
               trControl = control,
               )

# model <- train(x = training[,-c(1,2)], 
#                y = training$Label,
#                method = "glmnet",
#                trControl = control,
#                metric = "ROC",
#                tuneGrid = expand.grid(alpha = 1,
#                                       lambda = seq(0.001,0.1,by = 0.001)))


model

testPred <- predict(model, testing)

View(data.table(known = training$Label,
           prediction = testPred))

# estimate variable importance
importance <- varImp(model, scale=FALSE)

# summarize importance
print(importance)

# plot importance
plot(importance)

# ------------- FROM METABOSHINY --------------

mat <- dataSet$norm
labels <- dataSet$cls

inTrain <- createDataPartition(y = labels,
                               ## the outcome data are needed
                               p = .60,
                               ## The percentage of data in the training set
                               list = FALSE)

training <- mat[ inTrain,]
training.labels <- labels[inTrain]

testing <- mat[-inTrain,]
testing.labels <- labels[-inTrain]

# prepare training scheme
control <- trainControl(method="LOOCV", 
                        number=10, 
                        returnResamp="all",
                        classProbs=TRUE, 
                        summaryFunction=twoClassSummary)

# train the model
set.seed(849)

model <- train(x = training, 
               y = training.labels,
               method = "glmnet",
               trControl = control,
               metric = "ROC",
               tuneGrid = expand.grid(alpha = 1,
                                      lambda = seq(0.001,0.1,by = 0.001)))

# Using caret to perform CV

testPred <- predict(model, testing)

data.table(known = testing.labels,
           prediction = testPred)

# estimate variable importance
importance <- varImp(model, scale=FALSE)

# summarize importance
print(importance)
  
View(get_matches("backend/db/chebi.full.db", cpd = 157.0815375867, searchid=NULL))