############################
###
### Author: Joshua Snoke
### Title: Matrix input for partitioned estimation
### Description: Changing previous code to allow for OpenMx matrix objects as input
### Date: 12.2.15
############################
library(OpenMx)
library(optimx)
mxVersion()

#####
## OpenMx matrix objects - input
#####

matrA = mxMatrix( type="Full", nrow=6, ncol=6,
                  free =  c(F,F,F,F,F,F,
                            F,F,F,F,F,F,
                            F,F,F,F,F,F,
                            T,F,F,F,F,F,
                            F,F,F,F,F,F,
                            F,T,F,F,F,F),
                  values = c(0,0,0,0,0,0,
                             0,0,0,0,0,0,
                             1,0,0,0,0,0,
                             0.6,0,0,0,0,0,
                             0,1,0,0,0,0,
                             0,0.6,0,0,0,0),
                  labels = c(NA,NA,NA,NA,NA,NA,
                             NA,NA,NA,NA,NA,NA,
                             NA,NA,NA,NA,NA,NA,
                             "regx2",NA,NA,NA,NA,NA,
                             NA,NA,NA,NA,NA,NA,
                             NA,"regx4",NA,NA,NA,NA),
                  byrow = TRUE, name = "A" )
matrS = mxMatrix( type="Symm", nrow=6, ncol=6,
                  free =  c(T,T,F,F,F,F,
                            T,T,F,F,F,F,
                            F,F,T,F,F,F,
                            F,F,F,T,F,F,
                            F,F,F,F,T,F,
                            F,F,F,F,F,T),
                  values = c(1, 0.5, 0, 0, 0, 0,
                             0.5, 1, 0, 0, 0, 0,
                             0, 0, 1, 0,0,0,
                             0,0,0, 1, 0,0,
                             0,0,0,0, 1, 0,
                             0,0,0,0,0, 1),
                  labels = c("varFac1","covFac12",NA,NA,NA,NA,
                             "covFac12","varFac2",NA,NA,NA,NA,
                             NA,NA,"varx1",NA,NA,NA,
                             NA,NA,NA,"varx2",NA,NA,
                             NA,NA,NA,NA,"varx3",NA,
                             NA,NA,NA,NA,NA,"varx4"),
                  lbound = c(rep(0.1, 2), rep(NA, 4),
                             rep(0.1, 2), rep(NA, 4),
                             rep(NA, 2), rep(0.1, 1), rep(NA, 3),
                             rep(NA, 3), rep(0.1, 1), rep(NA, 2),
                             rep(NA, 4), rep(0.1, 1), rep(NA, 1),
                             rep(NA, 5), rep(0.1, 1)),
                  byrow=TRUE, name="S" )
matrF = mxMatrix( type = "Full", nrow = 4, ncol = 6,
                  free=FALSE,
                  values=c(0,0,1,0,0,0,
                           0,0,0,1,0,0,
                           0,0,0,0,1,0,
                           0,0,0,0,0,1),
                  byrow = TRUE, name = "F" )
matrM = mxMatrix( type="Full", nrow = 1, ncol = 6,
                  free=c(F, F, T, T, T, T),
                  values=c(0,0,2,2,2,2),
                  labels=c(NA, NA, "meanx1", "meanx2", "meanx3", "meanx4"),
                  name="M" )
matrI = mxMatrix(type = "Iden", nrow = 6, ncol = 6, name = "I")

matHold = list(matrA, matrS, matrF, matrM, matrI)


#####
## generating data
#####

newGenData = function(ramMatr, numObs, colNames = NULL){
  for(i in 1:length(ramMatr)){
    names(ramMatr)[i] = ramMatr[[i]]$name
    if(ramMatr[[i]]$name == "F")
      names(ramMatr)[i] = "Fm"
  }
  
  for(j in c("A", "S", "Fm", "M", "I")){
    if(is.na(match(j, names(ramMatr))) == TRUE)
      stop("Input A, S, Fm, M, or I matrices missing or misnamed")
  }
  
  fullFacCov = ramMatr$Fm$values %*% solve(ramMatr$I$values -ramMatr$A$values) %*% 
    ramMatr$S$values %*% t(solve(ramMatr$I$values - ramMatr$A$values)) %*% t(ramMatr$Fm$values)
  
  fullFacMean = t(ramMatr$Fm$values %*% solve(ramMatr$I$values - ramMatr$A$values) %*% 
                    as.numeric(t(ramMatr$M$values)))
  
  facModData = mvrnorm(n = numObs, mu = fullFacMean, Sigma = fullFacCov)
  if(is.null(colNames) == FALSE)
    colnames(facModData) = colNames
  else
    colnames(facModData) = paste('x', c(1:ncol(facModData)), sep = '')
  
  output = list(data = facModData, cov = fullFacCov, mean = fullFacMean)
  return(output)
}

facData = newGenData(matHold, 500)

dataList = list(data1 = facData$data[, c(1,3), drop = FALSE],
                data2 = facData$data[, c(2,4), drop = FALSE])

dataLocationMatrix = data.frame(cbind(rep(1, 500),
                                      rep(1, 500),
                                      c(rep(2, 200), rep(3, 300)),
                                      rep(4, 500),
                                      rep(4, 500),
                                      c(rep(5, 200), rep(6, 300)),
                                      c(rep(7, 100), rep(8, 200), rep(9, 200))) )

partition = split(unique(dataLocationMatrix), row.names(unique(dataLocationMatrix)))

whichRowPart = function(i, dataLocationMatrix, part){
  sum(dataLocationMatrix[i, ] != part)
}

whichRowAllPart = function(part, dataLocationMatrix){
  which(unlist(lapply(1:nrow(dataLocationMatrix), whichRowPart, 
                      dataLocationMatrix = dataLocationMatrix, part = part)) == 0)

  }

partObs = lapply(partition, whichRowAllPart, dataLocationMatrix = dataLocationMatrix)

#partition = c(rep(1, 2), rep(2, 2))
#partition = c(1, 2, 1, 2)


#####
## Estimation functions
#####

nodeCentralEst = function(inputSecret,inputShared){
  ##Full Cov and Mean
  matrices = inputSecret$model
  matReOrder = inputSecret$reOrder
  #inv = chol2inv(chol((matrices$I-matrices$A)))
  inv = solve(matrices$I - matrices$A)
  fullCov = matrices$Fm %*% inv %*% matrices$S %*% t(inv) %*% t(matrices$Fm)
  fullMean = t(matrices$Fm %*% inv %*% as.numeric(t(matrices$M)))
  ## reorder matrices
  fullCov = fullCov[matReOrder, matReOrder]
  fullMean = fullMean[, matReOrder, drop = F]
  inputShared$matrices = list(fullCov = fullCov, fullMean = fullMean)
  output = list(secret = inputSecret, shared = inputShared)
  return(output)
}

firstNodeExternalEst = function(inputSecret,inputShared){
  data = inputSecret
  meanSelf = inputShared$matrices$fullMean[,1:ncol(data)]
  meanOther = inputShared$matrices$fullMean[,-(1:ncol(data))]
  covSelf = inputShared$matrices$fullCov[1:ncol(data),1:ncol(data),drop=FALSE]
  covOther = inputShared$matrices$fullCov[-(1:ncol(data)),-(1:ncol(data)),drop=FALSE]
  covBoth = inputShared$matrices$fullCov[-(1:ncol(data)),1:ncol(data),drop=FALSE]
  #invSelf = chol2inv(chol(covSelf))
  invSelf = solve(covSelf)
  
  ##Data minus expected mean
  diffMean = t(t(data)-meanSelf)
  
  ##Log Likelihoods
  LL = matrix(NA,nrow=nrow(data),ncol=1)
  for(i in 1:nrow(data)){
    tempInv = determinant(covSelf, log=TRUE) 
    if(tempInv$sign <= 0) {
      tempInv = -Inf 
    }
      else {
        tempInv = tempInv$modulus
        attributes(tempInv) = NULL
      }
    LL[i,] = ncol(data)*log(2*pi)+ tempInv +
      t(as.matrix(diffMean[i,]))%*%invSelf%*%as.matrix(diffMean[i,])
  }
  
  ##Conditional Model parameters
  condCov = covOther - covBoth%*%invSelf%*%t(covBoth)
  
  condMean = matrix(NA,nrow=nrow(data),ncol=ncol(covOther))
  for(i in 1:nrow(data)){
    condMean[i,] = t(as.matrix(meanOther) + covBoth%*%invSelf%*%as.matrix(diffMean[i,]))
  }
  
  newLL = sum(LL) + inputShared$LL  
  output = list(secret = data,
                shared = list(matrices = list(fullCov = condCov, fullMean = condMean),
                              LL = newLL))
  return(output)
}

nodeExternalEst = function(inputSecret,inputShared){
  data = inputSecret
  meanSelf = inputShared$matrices$fullMean[,1:ncol(data)]
  meanOther = inputShared$matrices$fullMean[,-(1:ncol(data)),drop=FALSE]
  covSelf = inputShared$matrices$fullCov[1:ncol(data),1:ncol(data),drop=FALSE]
  covOther = inputShared$matrices$fullCov[-(1:ncol(data)),-(1:ncol(data)),drop=FALSE]
  covBoth = inputShared$matrices$fullCov[-(1:ncol(data)),1:ncol(data),drop=FALSE]
  #invSelf = chol2inv(chol(covSelf))
  invSelf = solve(covSelf)
  
  ##Data minus expected mean
  diffMean = data-meanSelf
  
  ##Log Likelihoods
  LL = matrix(NA,nrow=nrow(data),ncol=1)
  for(i in 1:nrow(data)){
    tempInv = determinant(covSelf, log=TRUE) 
    if(tempInv$sign <= 0) {
      tempInv = -Inf 
    }
    else {
      tempInv = tempInv$modulus
      attributes(tempInv) = NULL
    }
    LL[i,] = ncol(data)*log(2*pi) + tempInv +
      t(as.matrix(diffMean[i,]))%*%invSelf%*%as.matrix(diffMean[i,])
  }
  
  ##Conditional Model parameters
  condCov = covOther - covBoth%*%invSelf%*%t(covBoth)
  
  condMean = matrix(NA,nrow=nrow(data),ncol=ncol(covOther))
  for(i in 1:nrow(data)){
    condMean[i,] = meanOther[i,] + t(covBoth%*%invSelf%*%diffMean[i,])
  }
  
  newLL = sum(LL) + inputShared$LL  
  output = list(secret = data,
                shared = list(matrices = list(fullCov = condCov, fullMean = condMean),
                              LL = newLL))
  return(output)
}

lastNodeExternalEst = function(inputSecret,inputShared){
  data = inputSecret
  meanSelf = inputShared$matrices$fullMean
  covSelf = inputShared$matrices$fullCov
  #invSelf = chol2inv(chol(covSelf))
  invSelf = solve(covSelf)
  
  ##Data minus expected mean
  diffMean = data-meanSelf
  
  ##Log Likelihoods
  LL = matrix(NA,nrow=nrow(data),ncol=1)
  for(i in 1:nrow(data)){
    tempInv = determinant(covSelf, log=TRUE) 
    if(tempInv$sign <= 0) {
      tempInv = -Inf 
    }
    else {
      tempInv = tempInv$modulus
      attributes(tempInv) = NULL
    }
    LL[i,] = ncol(data)*log(2*pi)+ tempInv +
      t(as.matrix(diffMean[i,]))%*%invSelf%*%as.matrix(diffMean[i,])
  }
  
  newLL = sum(LL) + inputShared$LL
  output = list(secret = data,
                shared = newLL)
  return(output)
}

soloNodeExternalEst = function(inputSecret,inputShared){
  data = inputSecret
  meanSelf = inputShared$matrices$fullMean[,1:ncol(data)]
  covSelf = inputShared$matrices$fullCov
  #invSelf = chol2inv(chol(covSelf))
  invSelf = solve(covSelf)
  
  ##Data minus expected mean
  diffMean = t(t(data)-meanSelf)
  
  ##Log Likelihoods
  LL = matrix(NA,nrow=nrow(data),ncol=1)
  for(i in 1:nrow(data)){
    tempInv = determinant(covSelf, log=TRUE) 
    if(tempInv$sign <= 0) {
      tempInv = -Inf 
    }
    else {
      tempInv = tempInv$modulus
      attributes(tempInv) = NULL
    }
    LL[i,] = ncol(data)*log(2*pi)+ tempInv +
      t(as.matrix(diffMean[i,]))%*%invSelf%*%as.matrix(diffMean[i,])
  }
  
  newLL = sum(LL) + inputShared$LL
  output = list(secret = data,
                shared = newLL)
  return(output)
}

setNewParam = function(ramMatr, param){
  ## insert new parameters to RAM matrices
  tempMat = ramMatr
  for(i in 1:length(tempMat)){
    ## insert new params
    tempVec = tempMat[[i]]$labels[!is.na(tempMat[[i]]$labels)]
    
    tempVec2 = tempVec[!duplicated(tempVec)]
    
    tempMat[[i]]$values[match(tempVec2, tempMat[[i]]$labels)] = 
      param[is.element(names(param), tempMat[[i]]$labels)]
    
    ## insert duplicated (e.g. covariances)
    for(j in 1:length(tempMat[[i]]$free[duplicated(tempMat[[i]]$labels, MARGIN = 0)] == TRUE)){
      
      ind = which(duplicated(tempMat[[i]]$labels, MARGIN = 0))[j]
      if(tempMat[[i]]$free[ind] == TRUE){
        
        tempMat[[i]]$values[ind] = param[intersect(tempMat[[i]]$labels[ind], tempVec2)]
      }
    }
  }
  
  inputMat = list("A" = tempMat$A$values, "Fm" = tempMat$Fm$values, "M" = tempMat$M$values,
                  "S" = tempMat$S$values, "I" = tempMat$I$values)
  
  return(inputMat)
  
}

newFitFunc = function(param, dataList, parts, horzPart, ramMatr, reOrder){
  
  numPart = length(unique(parts))
  
  inputMat = setNewParam(ramMatr, param)
  
  ##Set Information
  nodeCentralInfo = list(secret = list(model = inputMat, reOrder = reOrder),
                         shared = list(matrices = vector("list",2), LL = c(0)))
  
  nodeExternalInfo = vector("list",numPart)
  for(i in 1:numPart){
    nodeExternalInfo[[i]] = list(secret = dataList[[parts[i]]],
                                 shared = list(matrices = vector("list",2), LL = c(0)))
  }
  
  ##Node Central Start Process
  nodeCentralInfo = nodeCentralEst(nodeCentralInfo$secret,nodeCentralInfo$shared)
  
  ##First node receives from central node
  if(numPart > 1){
    nodeExternalInfo[[1]] = firstNodeExternalEst(nodeExternalInfo[[1]]$secret,nodeCentralInfo$shared)
    nodeExternalInfo[[2]]$shared = nodeExternalInfo[[1]]$shared
    ##Iteration of Node External Process
    if(numPart > 2){
      for(i in 2:(numPart-1)){
        nodeExternalInfo[[i]] = nodeExternalEst(nodeExternalInfo[[i]]$secret, 
                                                nodeExternalInfo[[i]]$shared) ##Within node compute
        nodeExternalInfo[[i+1]]$shared = nodeExternalInfo[[i]]$shared ##Transfer
      } 
    }
    ##Last Node, no conditional comp needed
    nodeExternalInfo[[numPart]] = lastNodeExternalEst(nodeExternalInfo[[numPart]]$secret,
                                                      nodeExternalInfo[[numPart]]$shared)
  }
  ##For case of only one node
  else if(numPart == 1){
    nodeExternalInfo[[numPart]] = soloNodeExternalEst(nodeExternalInfo[[numPart]]$secret,
                                                      nodeCentralInfo$shared) 
  }
  
  output = nodeExternalInfo[[numPart]]$shared
  return(output)
}
#newFitFunc(freeValues, dataList, partition, ramMatr, matReOrder)

complexWrap = function(param, parts, obs, dataList, ramMatr){

  horzPartLL = vector("list", length(parts))
  
  for(a in 1:length(parts)){
    ##reorder for observed variables according to partitions
    matReOrder = rep(NA, length(unique(parts)))
    tempLen = c(1, 0)
    for(i in 1:length(unique(parts))){
      tempPart = which(partitions == i)
      tempLen[2] = tempLen[2] + length(tempPart)
      matReOrder[tempLen[1]:tempLen[2]] = tempPart
      tempLen[1] = 1 + tempLen[2]
    }
    
    horzPartLL[[a]] = newFitFunc(param, dataList, parts[[a]], horzPart[[a]], ramMatr, reOrder)
  }
  
  return(sum(horzPartLL))
  
}

newWrapFunc = function(ramMatr, partitions, observations, data){
  
  ## set names
  for(i in 1:length(ramMatr)){
    names(ramMatr)[i] = ramMatr[[i]]$name
    if(ramMatr[[i]]$name == "F")
      names(ramMatr)[i] = "Fm"
  }
  
  ## check and standardize names
  for(j in c("A", "S", "Fm", "M", "I")){
    if(is.na(match(j, names(ramMatr))) == TRUE)
      stop("Input A, S, Fm, M, or I matrices missing or misnamed")
  }
  
  ## set initial values and bounds based on user input matrices
  freeValues = c(ramMatr$A$values[ramMatr$A$free == T], ramMatr$S$values[ramMatr$S$free == T], 
                 ramMatr$Fm$values[ramMatr$Fm$free == T], ramMatr$M$values[ramMatr$M$free == T])
  upperB = c(ramMatr$A$ubound[ramMatr$A$free == T], ramMatr$S$ubound[ramMatr$S$free == T], 
             ramMatr$Fm$ubound[ramMatr$Fm$free == T], ramMatr$M$ubound[ramMatr$M$free == T])
  lowerB = c(ramMatr$A$lbound[ramMatr$A$free == T], ramMatr$S$lbound[ramMatr$S$free == T], 
             ramMatr$Fm$lbound[ramMatr$Fm$free == T], ramMatr$M$lbound[ramMatr$M$free == T])
  upperB[is.na(upperB) == TRUE] = Inf
  lowerB[is.na(lowerB) == TRUE] = -Inf
  
  names(freeValues) = c(ramMatr$A$labels[ramMatr$A$free == T], ramMatr$S$labels[ramMatr$S$free == T], 
                        ramMatr$Fm$labels[ramMatr$Fm$free == T], ramMatr$M$labels[ramMatr$M$free == T])
  upperB = upperB[!duplicated(names(freeValues))]
  lowerB = lowerB[!duplicated(names(freeValues))]
  freeValues = freeValues[!duplicated(names(freeValues))]
  
  ## optimize
  estimates = optimx(par = freeValues, fn = complexWrap,
                     dataList = data,
                     parts = partitions,
                     obs = observations,
                     ramMatr = ramMatr,
                     reOrder = matReOrder,
                     method = "L-BFGS-B",
                     lower = lowerB,
                     upper = upperB,
                     itnmax = 1000,
                     hessian = TRUE,
                     #control = list(trace = 5, maximize = TRUE,
                     gr = NULL)
  return(estimates)
}


#####
## test
#####

test2 = newWrapFunc(matHold, partition, partObs, dataList)
test2

test = newWrapFunc(matHold, rep(1, 4), facData)
test

mxMod1 = mxModel(model = "testModel", type = "RAM", matrA, matrM, matrS, matrI, matrF, 
                 mxData(facData$data, type = "raw"), 
                 mxAlgebra(F%*%solve(I-A)%*%S%*%t(solve(I-A))%*%t(F), name="expCov"),
                 mxAlgebra(t(F%*%solve(I-A)%*%t(M)), name="expMean"),
                 mxExpectationNormal(covariance="expCov",means="expMean",
                                     dimnames = colnames(facData$data)))
fitMod1 = mxRun(mxMod1)
fitMod1$output$estimate
test
test2
## add real values
## add openmx results
## table form


