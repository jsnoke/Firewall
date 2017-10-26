############################
###
### Author: Joshua Snoke
### Title: test_code.R
### Description: Code to be run on the sensitive data maintained by Michael Hunter.
### Results to be included in paper on secure maximum likelihood estimation.
### Date: 14.06.17
### 
############################

rm(list = ls())
set.seed(2385) ## to preserve randomness from imputation for replicability, can remove for other cases

source("algorithm_functions.R") ## run all functions from other file

###
###
### load packages and functions
###
###
library(OpenMx)
library(optimx)


###
###
### load FAKE datasets
### REAL datasets are behind firewall and code was run on them by Micahel Hunter
### need to change to appropriate paths and filenames
###
###
#dataFakeDhs = read.csv("~/Box Sync/Firewall Project/SecureAlgorithm/dataExample/dataFakeDhs.csv")
dataFakeOuhsc1a = read.csv("~/Box Sync/Firewall Project/SecureAlgorithm/dataExample/dataFakeOuhsc1a.csv")
dataFakeOuhsc1b = read.csv("~/Box Sync/Firewall Project/SecureAlgorithm/dataExample/dataFakeOuhsc1b.csv")
dataFakeOuhsc27 = read.csv("~/Box Sync/Firewall Project/SecureAlgorithm/dataExample/dataFakeOuhsc27.csv")


###
###
### make combined dataset 
### to ensure correct observation linkage
###
###
## combine 1a and 1b
dataFakeOuhsc1Comb = rbind(cbind("dfInd" = "1a", dataFakeOuhsc1a), cbind("dfInd" = "1b", dataFakeOuhsc1b))[!duplicated(rbind(dataFakeOuhsc1a, 
                                                                                             dataFakeOuhsc1b)[, 1:2]), ]
dataFakeOuhsc1Comb = subset(dataFakeOuhsc1Comb, !is.na(FamilyID) & !is.na(ChildID))

## 1 child ids
id1 = levels(interaction(dataFakeOuhsc1Comb$FamilyID, dataFakeOuhsc1Comb$ChildID, drop = T, sep = " "))

## 2-7 child ids
dataFakeOuhsc27Drop = dataFakeOuhsc27[!duplicated(dataFakeOuhsc27[, 1:3]), ]
dataFakeOuhsc27Drop = subset(dataFakeOuhsc27Drop, !is.na(FamilyID) & !is.na(ChildID))
id27 = levels(interaction(dataFakeOuhsc27Drop$FamilyID, dataFakeOuhsc27Drop$ChildID, drop = T, sep = " "))

## combine
idMatrix = matrix(NA, nrow = (length(id1) + length(id27)), ncol = 2)
for(a in 1:length(id1)){
  idMatrix[a, ] = strsplit(id1[[a]], " ")[[1]]
}
for(a in 1:length(id27)){
  idMatrix[(a + length(id1)), ] = strsplit(id27[[a]], " ")[[1]]
}

idMatrix = idMatrix[!duplicated(idMatrix), ]
colnames(idMatrix) = c("FamilyID", "ChildID")

## create wide DF
wideWaveData = data.frame(matrix(NA, ncol = 10, nrow = nrow(idMatrix)))
wideWaveData[, 1:2] = idMatrix
colnames(wideWaveData) = c("FamilyID", "ChildID", "dfInd", "Wave1", "Wave2", "Wave3", "Wave4", "Wave5", "Wave6", "Wave7")
for(a in 1:nrow(dataFakeOuhsc27Drop)){
  wideWaveData[(wideWaveData[, 1] == dataFakeOuhsc27Drop$FamilyID[a] & wideWaveData[, 2] == dataFakeOuhsc27Drop$ChildID[a]), 
               (dataFakeOuhsc27Drop$WAVE[a] + 3)] = dataFakeOuhsc27Drop$fnsScore[a]
}
for(a in 1:nrow(dataFakeOuhsc1Comb)){
  wideWaveData[(wideWaveData[, 1] == dataFakeOuhsc1Comb$FamilyID[a] & wideWaveData[, 2] == dataFakeOuhsc1Comb$ChildID[a]), 3:4] = 
    c(dataFakeOuhsc1Comb$dfInd[a], dataFakeOuhsc1Comb$fnsScore[a])
}
wideWaveData$dfInd[is.na(wideWaveData$dfInd)] = sample(1:2, sum(is.na(wideWaveData$dfInd)), replace = T)


###
###
### impute missingness
### once joint and once marginal
###
###
createImputedDF = function(origDF){
  imputeComb = origDF
  
  ## first var
  missX1 = which(is.na(imputeComb[, 1]))
  imputeComb[is.na(imputeComb[, 1]), 1] = mean(imputeComb[, 1], na.rm = T)
  
  for(a in 2:ncol(imputeComb)){
    tempX = cbind(1, imputeComb[, 1:(a - 1)])
    tempY = imputeComb[, a]
    tempCovInv = solve(t(tempX[!is.na(tempY), ]) %*% tempX[!is.na(tempY), ])
    tempCoef = t(tempY[!is.na(tempY)]) %*% tempX[!is.na(tempY), ] %*% tempCovInv
    tempResiduals = tempY[!is.na(tempY)] - tempX[!is.na(tempY), ] %*% t(tempCoef)
    tempSigma = sqrt((sum(tempResiduals ^ 2)) / (length(tempY) - ncol(tempX) - 1))
    newY = tempX[is.na(tempY), ] %*% t(tempCoef) + rnorm(nrow(tempX[is.na(tempY), ])) * tempSigma
    imputeComb[is.na(imputeComb[, a]), a] = newY
  }
  
  tempX = cbind(1, imputeComb[, 2:ncol(imputeComb)])
  tempY = imputeComb[, 1]
  tempCovInv = solve(t(tempX[-(missX1), ]) %*% tempX[-(missX1), ])
  tempCoef = t(tempY[-(missX1)]) %*% tempX[-(missX1), ] %*% tempCovInv
  tempResiduals = tempY[-(missX1)] - tempX[-(missX1), ] %*% t(tempCoef)
  tempSigma = sqrt((sum(tempResiduals ^ 2)) / (length(tempY) - ncol(tempX) - 1))
  newY = tempX[missX1, ] %*% t(tempCoef) + rnorm(nrow(tempX[missX1, ])) * tempSigma
  imputeComb[missX1, 1] = newY
  
  return(imputeComb)
}

fullImpute = cbind("id" = wideWaveData$dfInd, createImputedDF(as.matrix(wideWaveData[, -c(1:3)]))) ## imputed with all data together

DF1 = wideWaveData[wideWaveData$dfInd == 1, c(1:2, 4)]
DF2 = wideWaveData[wideWaveData$dfInd == 2, c(1:2, 4)]
DF3 = wideWaveData[, c(1:2, 5:10)]

imputedDF1 = cbind(as.numeric(rownames(DF1)), "id" = 1, DF1$Wave1)
imputedDF1[is.na(imputedDF1[, 3]), 3] = mean(imputedDF1[, 3], na.rm = T)
colnames(imputedDF1) = c("rowname", "id", "Wave1")

imputedDF2 = cbind(as.numeric(rownames(DF2)), "id" = 2, DF2$Wave1)
imputedDF2[is.na(imputedDF2[, 3]), 3] = mean(imputedDF2[, 3], na.rm = T)
colnames(imputedDF2) = c("rowname", "id", "Wave1")

imputedDF3 = cbind(as.numeric(rownames(DF3)), "id" = wideWaveData$dfInd, createImputedDF(as.matrix(DF3[, 3:8])))
imputedDF3 = imputedDF3[order(imputedDF3[, 2]), ]

## data imputed jointly
dataList1 = dataList2 = list("vector", 3)
dataList1[[1]] = fullImpute[wideWaveData$dfInd == 1, 1:2]
dataList1[[2]] = fullImpute[wideWaveData$dfInd == 2, 1:2]
dataList1[[3]] = fullImpute[order(fullImpute[, 1]), c(1, 3:8)]

## data imputed marginally
dataList2[[1]] = imputedDF1[, -1]
dataList2[[2]] = imputedDF2[, -1]
dataList2[[3]] = imputedDF3[, -1]
combMargImpute  = cbind(rbind(imputedDF1[, -1], imputedDF2[, -1]), imputedDF3[, -(1:2)])


###
###
### test!
###
###
matHold = genRAMSimulationMatrices(7)
horzInd = matrix(c(1, rep(3, 6), 2, rep(3, 6)), nrow = 2, ncol = 7, byrow = T)
partitions = c(1, rep(2, 6))

## test with jointly imputed data
test = secureWrapFunction(matHold, partitions, dataList1, numObs = 244, compareWithOpenMx = TRUE, complex = T, horzInd = horzInd,
                          noPartData = fullImpute[, 2:8], compareWithNoPart = T, verbose = 1, hessian = T)
#test$results

## test with marginally imputed data
test2 = secureWrapFunction(matHold, partitions, dataList2, numObs = 244, compareWithOpenMx = TRUE, complex = T, horzInd = horzInd,
                           noPartData = combMargImpute[, 2:8], compareWithNoPart = T, verbose = 1)
#test2$results

## test with FIML, no imputation
# Actually this uses casewise deletion
names(matHold) = c("A", "S", "F", "M", "I")

mxMod1 = mxModel(model = "testModel", type = "RAM", manifestVars = colnames(na.omit(wideWaveData[, 4:10])),
                 matHold$A, matHold$S, matHold$F, matHold$I, matHold$M, 
                mxData(wideWaveData[, 4:10], type = "raw"),
                mxAlgebra(F %*% solve(I - A) %*% (S %*% t(S)) %*% t(solve(I - A)) %*% t(F), 
                          name = "expCov"),
                mxAlgebra(t(F %*% solve(I - A) %*% t(M)), name = "expMean"),
                mxExpectationNormal(covariance = "expCov", means = "expMean",
                                    dimnames = colnames(na.omit(wideWaveData[, 4:10]))),
                mxFitFunctionML())
fitMod1 = mxRun(mxMod1)
#summary(fitMod1)
estimatesOpenMx = list(fitMod1$output$estimate, fitMod1$output$hessian)
#estimatesOpenMx

###
###
### save results
###
###
savePut = list("full_impute" = test, "marg_impute" = test2, "fiml" = estimatesOpenMx)

#save(savePut, file = "test_run_output.RData")
save(savePut, file = "test_run_output_real.RData")



