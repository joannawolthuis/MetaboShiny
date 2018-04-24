library(mixOmics)
data(liver.toxicity)

X <- as.matrix(mSet$dataSet$norm)
Y <- as.factor(mSet$dataSet$cls)

X <- as.matrix(liver.toxicity$gene)
Y <- as.factor(liver.toxicity$treatment[, 4])             

## PLS-DA function
plsda.res <- plsda(X, Y, ncomp = 5) # where ncomp is the number of components wanted
set.seed(2543) # for reproducibility here, only when the `cpus' argument is not used
perf.plsda <- perf(plsda.res, validation = "Mfold", folds = 5, 
                   progressBar = TRUE, auc = TRUE, nrepeat = 10) 
# perf.plsda.srbct$error.rate  # error rates
plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, legend.position = "horizontal")

samp <- sample(1:3, nrow(X), replace = TRUE) 

# 1/3 of the data will compose the test set
test <- which(samp == 1) 

# rest will compose the training set
train <- setdiff(1:nrow(X), test) 

## For PLS-DA, train the model
plsda.train <- plsda(X[train, ], Y[train], ncomp = 4)

# then predict
test.predict <- predict(plsda.train, X[test, ], dist = "max.dist")

plotIndiv(prediction)
# store prediction for the 4th component
prediction <- test.predict$class$max.dist[,4] 

# calculate the error rate of the model
confusion.mat = get.confusion_matrix(truth = Y[test], predicted = prediction)
get.BER(confusion.mat)
