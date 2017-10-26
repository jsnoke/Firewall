############################
###
### Author: Joshua Snoke
### Title: algorithm_functions.R
### Description: Internal functions and wrapper functions for secure algorithm
###               producing MLE across partitioned databases
### Date: 14.06.17
### 
############################

###
###
### function to estimate LL in non-partitioned setting
### for comparison
###
###
nonPartitionedLL = function(param, ramMatr, dataF){
  inputMat = setNewParam(ramMatr, param)
  
  inv = solve(inputMat$I - inputMat$A)
  fullCov = inputMat$Fm %*% inv %*% (inputMat$S %*% t(inputMat$S)) %*% t(inv) %*% t(inputMat$Fm)
  fullMean = t(inputMat$Fm %*% inv %*% as.numeric(t(inputMat$M)))
  
  numObs = nrow(dataF)
  numVar = ncol(dataF)
  selfCovInv = solve(fullCov)
  
  for(a in 1:numObs){
    LL = rep(NA, numObs)
    if(class(fullCov) == "matrix"){
      tempDet = determinant(fullCov, log = TRUE) 
      if(tempDet$sign <= 0) {
        tempDet = -Inf 
      } else {
        tempDet = tempDet$modulus
        attributes(tempDet) = NULL
      }
    } else{
      if(fullCov < 0){
        tempDet = -Inf
      } else{
        tempDet = log(fullCov)
      }
    }
    for(a in 1:numObs){
      LL[a] = numVar * log(2 * pi) + tempDet + (dataF[a, ] - fullMean) %*% selfCovInv %*% t((dataF[a, ] - fullMean))
    }
  }
  return(sum(LL))
}


###
###
### function generate the RAM matrices based on desired model
###
###
genRAMSimulationMatrices = function(p){
    matrA = mxMatrix("Full", nrow = (p + 2), ncol = (p + 2),
                     values = rbind(c(rep(0.0, p + 2)), 
                                    c(rep(0.0, p + 2)),
                                    cbind(rep(1.0, p), seq(0, (p - 1)), diag(0, p))), 
                     free = F,
                     labels = NA,
                     byrow = TRUE, name = "A")
    
    matrS = mxMatrix("Full", nrow = (p + 2), ncol = (p + 2),
                     values = rbind(c(0.5, 0, rep(0.0, p)),
                                    c(0.2, 0.5, rep(0.0, p)),
                                    cbind(rep(0.0, p), rep(0.0, p), diag(1.0, p))), 
                     free = rbind(c(rep(T, 1), rep(F, p + 1)),
                                  c(rep(T, 2), rep(F, p)),
                                  cbind(rep(F, p), rep(F, p), 
                                        matrix(as.logical(diag(1.0, p)), nrow = p))),
                     labels = rbind(c("sigma_i", rep(NA, p + 1)),
                                    c("COV_icept_slope", "sigma_s", rep(NA, p)),
                                    cbind(rep(NA, p), rep(NA, p), 
                                          matrix(c(rep(c("e", rep(NA, p)), p)), 
                                                 nrow = p, ncol = p, byrow = T)
                                          #matrix(NA, nrow = p, ncol = p, byrow = T)
                                          )),
                     lbound = rbind(c(rep(0, 1), rep(NA, p + 1)),
                                    c(rep(0, 2), rep(NA, p)),
                                    cbind(rep(NA, p), rep(NA, p), 
                                          matrix(c(rep(c(1e-5, rep(NA, p)), (p))), 
                                                 nrow = p, ncol = p, byrow = T))), 
                     byrow = TRUE, name = "S")
    #diag(matrS$labels[-(1:2), -(1:2)]) = paste("e", 1:p, sep = "")
    
    matrF = mxMatrix("Full", nrow = p, ncol = (p + 2), 
                     values = cbind(rep(0, p), rep(0, p), diag(1, p)), 
                     free = F, labels = NA, byrow = TRUE, name = "F")
    
    matrM = mxMatrix("Full", nrow = 1, ncol = (p + 2), values = rep(0.0, (p + 2)), 
                     free = c(rep(T, 2), rep(F, (p))), labels = c("mean_i", "mean_s", rep(NA, (p))), 
                     byrow = TRUE, name = "M")
    
    matrI = mxMatrix(type = "Iden", nrow = (p + 2), ncol = (p + 2), name = "I")
    
    matHold = list(matrA, matrS, matrF, matrM, matrI)
    return(matHold)
    
}


###
###
### the following are the internal functions of the secure likelihood algorithm
###
###
centralNodeOpenComputation = function(inputSelf){
    ##Full Cov and Mean
    matrices = inputSelf$model
    matReOrder = inputSelf$reOrder
    numObs = inputSelf$numObs
    partitions = inputSelf$partitions
    numNodes = nlevels(as.factor(partitions))
    
    inv = solve(matrices$I - matrices$A)
    tempCov = matrices$Fm %*% inv %*% (matrices$S %*% t(matrices$S)) %*% t(inv) %*% t(matrices$Fm)
    fullMean = t(matrices$Fm %*% inv %*% as.numeric(t(matrices$M)))
    ## reorder matrices
    #if(complex == F){ ## TEMP QUICK FIX
    #    tempCov = tempCov[matReOrder, matReOrder]
    #    fullMean = fullMean[, matReOrder, drop = F]  
    #}
    
    nodeCov = vector("list", numNodes)
    fullCond = vector("list", numNodes)
    tempPart = partitions
    #tempCov = inputSelf$cov
    nodeNoisyMean = vector("list", 2)
    pNoise = vector("list", numNodes)
    nodeNoisyMean[[1]] = matrix(NA, nrow = numObs, ncol = length(partitions[partitions == 1]))
    nodeNoisyMean[[2]] = matrix(NA, nrow = numObs, ncol = length(partitions[partitions != 1]))
    centralPartialMean2 = vector("list", numNodes)
    
    for(a in 1:numNodes){
        nodeCov[[a]] = tempCov[which(tempPart == a), which(tempPart == a)]
        centralPartialMean2[[a]] = tempCov[which(tempPart != a), which(tempPart == a)] %*% solve(nodeCov[[a]])
        fullCond[[a]] = tempCov = tempCov[which(tempPart != a), which(tempPart != a)] - 
            tempCov[which(tempPart != a), which(tempPart == a)] %*% solve(nodeCov[[a]]) %*% 
            tempCov[which(tempPart == a), which(tempPart != a)]
        
        tempPart = tempPart[-(which(tempPart == a))]
        
        pNoise[[a]] = matrix(runif((numObs * length(partitions[partitions == a])), 500, 2000), 
                             nrow = numObs, ncol = length(partitions[partitions == a]))
        if(a == 1){
            for(b in 1:nrow(pNoise[[a]])){
                nodeNoisyMean[[1]][b, ] = fullMean[, which(partitions == 1)] + pNoise[[a]][b, ]
            }
            nonFirstPart = partitions[-(which(partitions == 1))]
        } else{
            for(b in 1:nrow(pNoise[[a]])){
                nodeNoisyMean[[2]][b, which(nonFirstPart == a)] = fullMean[, which(nonFirstPart == a)] + 
                    pNoise[[a]][b, ]
            }      
        }
    }
    
    output = list("nodeCov" = nodeCov, "nodeNoisyMean" = nodeNoisyMean, 
                  "pNoise" = pNoise, "centralPartialMean2" = centralPartialMean2)
    return(output)
}

externalNodeComputation = function(inputSelf, inputReceived, inputProduced){
    numObs = nrow(inputSelf$data)
    numVar = ncol(inputSelf$data)
    selfCovInv = solve(inputReceived$nodeCov)
    rNoise = matrix(runif(numObs * numVar, 500, 2000), nrow = numObs, ncol = numVar)
    qNoise = matrix(runif(numObs * numVar, 500, 2000), nrow = numVar, ncol = numObs)
    
    if(is.null(inputProduced$denoisedCondMeanSelf)){
        selfNoisyMeanDiff1 = inputSelf$data - inputReceived$nodeNoisyMean + rNoise
        selfNoisyMeanDiff2 = inputSelf$data - inputReceived$nodeNoisyMean - rNoise
    } else {
        selfNoisyMeanDiff1 = inputSelf$data - inputProduced$denoisedCondMeanSelf + rNoise
        selfNoisyMeanDiff2 = inputSelf$data - inputProduced$denoisedCondMeanSelf - rNoise
    }
    
    externalPartialMean1 = selfCovInv %*% t(selfNoisyMeanDiff1)
    externalPartialMean2 = (selfCovInv %*% t(selfNoisyMeanDiff2)) + qNoise
    #externalPartialMean1 = selfCovInv %*% t(selfNoisyMeanDiff1)
    #externalPartialMean2 = selfCovInv %*% t(selfNoisyMeanDiff2) + t(qNoise)
    
    LL = rep(NA, numObs)
    if(class(inputReceived$nodeCov) == "matrix"){
        tempDet = determinant(inputReceived$nodeCov, log = TRUE) 
        if(tempDet$sign <= 0) {
            tempDet = -Inf 
        } else {
            tempDet = tempDet$modulus
            attributes(tempDet) = NULL
        }
    } else{
        if(inputReceived$nodeCov < 0){
            tempDet = -Inf
        } else{
            tempDet = log(inputReceived$nodeCov)
        }
    }
    for(a in 1:numObs){
        LL[a] = numVar * log(2 * pi) + tempDet + 
            selfNoisyMeanDiff1[a, , drop = F] %*% selfCovInv %*% t(selfNoisyMeanDiff2[a, , drop = F]) +
            (rNoise[a, , drop = F] %*% selfCovInv %*% t(rNoise[a, , drop = F]))
    }
    
    if(inputProduced$denoisyLL != 0){
        newLL = inputProduced$denoisyLL + sum(LL)
    } else{
        newLL = sum(LL)
    }
    
    output = list("rNoise" = rNoise, "qNoise" = qNoise, 
                  "externalPartialMean1" = externalPartialMean1, "externalPartialMean2" = externalPartialMean2, 
                  "noisyLL" = newLL)
    return(output)
}

centralAdjustmentComputation = function(inputProduced, inputReceived, node){
    if(node > 1){
        centralPartialMean1 = inputReceived$nodeNoisyMean[[node]] + 
            t(inputProduced$centralPartialMean2[[node]] %*% inputProduced$nodeCov[[node]] %*% 
                  inputReceived$externalPartialMean1[[node]])
    } else if(node == 1){
        centralPartialMean1 = inputProduced$nodeNoisyMean[[2]] + 
            t(inputProduced$centralPartialMean2[[node]] %*% inputProduced$nodeCov[[node]] %*% 
                  inputReceived$externalPartialMean1[[node]])
    }
    
    return(centralPartialMean1)
}

externalAdjustmentComputation = function(inputSelf, inputReceived){
    colSelf = c(1:ncol(inputSelf$data))
    denoisyLL = inputReceived$noisyLL - sum(inputReceived$pNoise * t(inputReceived$qNoise))
    
    if(is.null(inputReceived$mNoise)){
        denoisedCondMean = inputReceived$centralPartialMean1 - 
            t(inputReceived$centralPartialMean2 %*% t(inputReceived$rNoise - inputReceived$pNoise))
    } else{
        denoisedCondMean = inputReceived$centralPartialMean1 - inputReceived$mNoise -
            t(inputReceived$centralPartialMean2 %*% t(inputReceived$rNoise - inputReceived$pNoise))
    }
    denoisedCondMean1 = denoisedCondMean[, colSelf]
    temp = denoisedCondMean[, -(colSelf), drop = F]
    mNoise = matrix(runif(nrow(inputSelf$data) * ncol(temp), 500, 2000), nrow = nrow(inputSelf$data), ncol = ncol(temp))
    denoisedCondMean2 = temp + mNoise
    
    output = list("denoisyLL" = denoisyLL, "denoisedCondMeanSelf" = denoisedCondMean1, 
                  "denoisedCondMeanOther" = denoisedCondMean2, "mNoise" = mNoise)
    return(output)
}

firstNodeFinalComputation = function(inputReceived){
    denoisyLL = inputReceived$noisyLL - sum(inputReceived$pNoise * t(inputReceived$qNoise))
    
    return(denoisyLL)
}

centralFinalComputation = function(inputSelf, inputProduced, inputReceived){
    numVar = nlevels(as.factor(inputSelf$partitions))
    numObs = inputSelf$numObs
    llNoise = vector("list", numVar)
    for(b in 1:numVar){
        llNoise[[b]] = rep(NA, inputSelf$numObs)
        for(a in 1:numObs){
            llNoise[[b]][a] = inputProduced$pNoise[[b]][a, ] %*% inputReceived$externalPartialMean1[[b]][, a, drop = F] + 
                inputProduced$pNoise[[b]][a, ] %*% inputReceived$externalPartialMean2[[b]][, a, drop = F] + 
                inputProduced$pNoise[[b]][a, ] %*% solve(inputProduced$nodeCov[[b]]) %*% 
                t(inputProduced$pNoise[[b]][a, , drop = F])
        }
    }
    
    trueLL = inputReceived$noisyLL + sum(unlist(llNoise))
    return(trueLL)
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


###
###
### complete secure fit function given a set of parameters
###
###
secureFitFunction = function(param, dataList, partitions, ramMatr, reOrder, numObs){
    
    numPart = length(unique(partitions))
    #reOrder = c(1, 2, 3, 4, 5, 6)
    #inputMat = list("A" = matHold$A$values, "S" = matHold$S$values, "Fm" = matHold$Fm$values, "M" = matHold$M$values,
    #                "I" = matHold$I$values)
    inputMat = setNewParam(ramMatr, param)
    
    ##Set Information
    nodeCentralInfo = list("self" = list("model" = inputMat, "reOrder" = reOrder, "partitions" = partitions, "numObs" = numObs),
                           "produced" = list("nodeCov" = vector("list", numPart), "nodeNoisyMean" = vector("list", numPart), 
                                             "pNoise" = vector("list", numPart), 
                                             "centralPartialMean1" = vector("list", numPart), 
                                             "centralPartialMean2" = vector("list", numPart), "LL" = c(0)),
                           "received" = list("nodeNoisyMean" = vector("list", numPart), 
                                             "externalPartialMean1" = vector("list", numPart), 
                                             "externalPartialMean2" = vector("list", numPart), "noisyLL" = c(0)))
    
    nodeExternalInfo = vector("list", numPart)
    for(i in 1:numPart){
        nodeExternalInfo[[i]] = list("self" = list("data" = dataList[[i]]),
                                     "received" = list("nodeCov" = NULL, "noisyLL" = c(0), "nodeNoisyMean" = NULL, 
                                                       "pNoise" = NULL, "rNoise" = NULL, "qNoise" = NULL, "mNoise" = NULL,
                                                       "centralPartialMean1" = NULL, "centralPartialMean2" = NULL),
                                     "produced" = list("denoisyLL" = c(0), "noisyLL" = c(0), "denoisedCondMeanSelf" = NULL, 
                                                       "mNoise" = NULL, "denoisedCondMeanOther" = NULL, "rNoise" = NULL,
                                                       "qNoise" = NULL, "externalPartialMean1" = NULL, 
                                                       "externalPartialMean2" = NULL))
    }
    
    ## central setup computation
    nodeCentralInfo$produced[c("nodeCov", "nodeNoisyMean", 
                               "pNoise", "centralPartialMean2")] = centralNodeOpenComputation(nodeCentralInfo$self)
    
    ## pass to first data node
    nodeExternalInfo[[1]]$received[c("nodeCov", "nodeNoisyMean", 
                                     "pNoise")] = list(nodeCentralInfo$produced$nodeCov[[1]],
                                                       nodeCentralInfo$produced$nodeNoisyMean[[1]],
                                                       nodeCentralInfo$produced$pNoise[[numPart]])
    
    ## first data node computation
    nodeExternalInfo[[1]]$produced[c("rNoise", "qNoise", "externalPartialMean1", "externalPartialMean2", 
                                     "noisyLL")] = externalNodeComputation(nodeExternalInfo[[1]]$self,
                                                                           nodeExternalInfo[[1]]$received,
                                                                           nodeExternalInfo[[1]]$produced)
    
    ## pass
    nodeCentralInfo$received$externalPartialMean1[[1]] = nodeExternalInfo[[1]]$produced$externalPartialMean1
    nodeCentralInfo$received$externalPartialMean2[[1]] = nodeExternalInfo[[1]]$produced$externalPartialMean2
    
    nodeExternalInfo[[2]]$received[c("rNoise", "qNoise", 
                                     "noisyLL")] = nodeExternalInfo[[1]]$produced[c("rNoise", "qNoise", 
                                                                                    "noisyLL")]
    
    ## compute
    ## loop here
    for(k in 2:numPart){
        nodeCentralInfo$produced$centralPartialMean1[[(k - 1)]] = centralAdjustmentComputation(nodeCentralInfo$produced, 
                                                                                               nodeCentralInfo$received, (k - 1))
        
        ## pass
        nodeExternalInfo[[k]]$received[c("nodeCov",  "centralPartialMean1", "centralPartialMean2",
                                         "pNoise")] = list(nodeCentralInfo$produced$nodeCov[[k]],
                                                           nodeCentralInfo$produced$centralPartialMean1[[(k - 1)]],
                                                           nodeCentralInfo$produced$centralPartialMean2[[(k - 1)]],
                                                           nodeCentralInfo$produced$pNoise[[(k - 1)]])
        
        ## compute 
        nodeExternalInfo[[k]]$produced[c("denoisyLL", "denoisedCondMeanSelf", "denoisedCondMeanOther", 
                                         "mNoise")] = externalAdjustmentComputation(nodeExternalInfo[[k]]$self, 
                                                                                    nodeExternalInfo[[k]]$received)
        
        nodeExternalInfo[[k]]$produced[c("rNoise", "qNoise", "externalPartialMean1", "externalPartialMean2", 
                                         "noisyLL")] = externalNodeComputation(nodeExternalInfo[[k]]$self,
                                                                               nodeExternalInfo[[k]]$received,
                                                                               nodeExternalInfo[[k]]$produced) 
        
        ## pass
        nodeCentralInfo$received$externalPartialMean1[[k]] = nodeExternalInfo[[k]]$produced$externalPartialMean1
        nodeCentralInfo$received$externalPartialMean2[[k]] = nodeExternalInfo[[k]]$produced$externalPartialMean2
        if(k != numPart){
            nodeCentralInfo$received$nodeNoisyMean[[k]] = nodeExternalInfo[[k]]$produced$denoisedCondMeanOther
            
            nodeExternalInfo[[(k + 1)]]$received[c("rNoise", "qNoise", "noisyLL", 
                                                   "mNoise")] = nodeExternalInfo[[k]]$produced[c("rNoise", "qNoise", 
                                                                                                 "noisyLL", "mNoise")]
        }
    }
    
    ## stop loop here
    nodeExternalInfo[[1]]$received[c("qNoise", "noisyLL")] = nodeExternalInfo[[numPart]]$produced[c("qNoise", "noisyLL")]
    
    ## compute
    nodeExternalInfo[[1]]$produced[c("denoisyLL")] = firstNodeFinalComputation(nodeExternalInfo[[1]]$received)
    
    ## pass
    nodeCentralInfo$received[c("noisyLL")] = nodeExternalInfo[[1]]$produced[c("denoisyLL")]
    
    ## compute
    nodeCentralInfo$produced[c("LL")] = centralFinalComputation(nodeCentralInfo$self, nodeCentralInfo$produced, 
                                                                nodeCentralInfo$received)
    return(nodeCentralInfo$produced$LL)
    
}


###
###
### wrapper for complex partitions
###
###
secureComplexWrapper = function(param, complex, horzInd, dataList, partitions, ramMatr, reOrder, numObs){
  if(complex == F){
    llRun = secureFitFunction(param = param, dataList = dataList, partitions = partitions, ramMatr = ramMatr, reOrder = reOrder, 
                              numObs = numObs)
  } else{
    subDataList = vector("list", nrow(horzInd))
    llRun = rep(NA, nrow(horzInd))
    for(a in 1:nrow(horzInd)){
      subDataList[[a]] = list(dataList[[horzInd[a, 1]]][, -which(colnames(dataList[[horzInd[a, 1]]]) %in% c("id")), drop = F], 
                              #dataList[[horzInd[a, 2]]][id = dataList[[horzInd[a, 1]]][, "id"], 
                              dataList[[horzInd[a, 2]]][dataList[[horzInd[a, 2]]][, which(colnames(dataList[[horzInd[a, 2]]]) 
                                                                                          %in% c("id"))] == a,
                                                        -which(colnames(dataList[[horzInd[a, 2]]]) %in% c("id")), drop = F])
      
      tempPart = as.numeric(as.factor(horzInd[a, ]))
      llRun[a] = secureFitFunction(param = param, dataList = subDataList[[a]], partitions = tempPart, ramMatr = ramMatr, 
                                   reOrder = reOrder, numObs = nrow(subDataList[[a]][[1]]))
    }
  }
  return(sum(llRun))
}


###
###
### the wrapper function for the entire secure model estimation
### including comparison to non-partitioned fitting through optimx and openmx
###
###
secureWrapFunction = function(ramMatr, partitions, data, numObs, complex = FALSE, horzInd = NULL,
                              compareWithNoPart = FALSE, compareWithOpenMx = FALSE,
                              noPartData = NULL, hessian = FALSE, verbose = 0){
    
    tempTime = proc.time()
    ## set names
    for(i in 1:length(ramMatr)){
        names(ramMatr)[i] = ramMatr[[i]]$name
        if(ramMatr[[i]]$name == "F")
            names(ramMatr)[i] = ramMatr[[i]]$name = "Fm"
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
    
    ##reorder for observed variables according to partitions
    if(complex == T){
        matReOrder = NULL ## TEMP QUICK FIX
    } else{
        matReOrder = rep(NA, length(unique(partitions)))
        tempLen = c(1, 0)
        for(i in 1:length(unique(partitions))){
            tempPart = which(partitions == i)
            tempLen[2] = tempLen[2] + length(tempPart)
            matReOrder[tempLen[1]:tempLen[2]] = tempPart
            tempLen[1] = 1 + tempLen[2]
        }
    }
    
    ## optimize
    estimatesUser = optimx(par = freeValues, fn = secureComplexWrapper,
                           complex = complex,
                           horzInd = horzInd,
                           dataList = data,
                           partitions = partitions,
                           ramMatr = ramMatr,
                           reOrder = matReOrder,
                           numObs = numObs,
                           method = "L-BFGS-B",
                           lower = lowerB,
                           upper = upperB,
                           itnmax = 1000,
                           hessian = hessian,
                           control = list(trace = verbose, factr = 1e10),
                           #control = list(trace = 5, maximize = TRUE,
                           gr = NULL)
    
    outTime = proc.time() - tempTime
    
    if(compareWithNoPart == TRUE){
        if(is.null(noPartData)){
            stop("please supply non-partitioned data to compare fit or set compareWithNoPart to FALSE")
        }
        if(class(noPartData) != "list"){
            noPartData = list(noPartData)
        }
        
        matReOrder = 1:ncol(noPartData[[1]])
        #partitions = rep(1, ncol(noPartData[[1]]))
        
        estimatesNoPart = optimx(par = freeValues, fn = nonPartitionedLL,
                                 dataF = noPartData[[1]],
                                 ramMatr = ramMatr,
                                 method = "L-BFGS-B",
                                 lower = lowerB,
                                 upper = upperB,
                                 itnmax = 1000,
                                 hessian = hessian,
                                 control = list(trace = verbose, factr = 1e10),
                                 #control = list(trace = 5, maximize = TRUE,
                                 gr = NULL)
        
        #estimatesNoPart = unlist(estimatesNoPart[1:attr(estimatesNoPart, "npar")])
    } else{
        estimatesNoPart = NULL
    }
    
    if(compareWithOpenMx == TRUE){
        if(is.null(noPartData)){
            stop("please supply non-partitioned data to compare fit or set compareWithOpenMx to FALSE")
        }
        if(class(noPartData) != "list"){
            noPartData = list(noPartData)
        }
        
        tempTime2 = proc.time()
        
        mxMod1 = mxModel(model = "testModel", type = "RAM", manifestVars = colnames(noPartData[[1]]),
                         ramMatr$A, ramMatr$S, ramMatr$Fm, ramMatr$I, ramMatr$M, 
                         mxData(noPartData[[1]], type = "raw"), 
                         mxAlgebra(Fm %*% solve(I - A) %*% (S %*% t(S)) %*% t(solve(I - A)) %*% t(Fm), 
                                   name="expCov"),
                         mxAlgebra(t(Fm %*% solve(I - A) %*% t(M)), name="expMean"),
                         mxExpectationNormal(covariance="expCov",means="expMean",
                                             dimnames = colnames(noPartData[[1]])))
        fitMod1 = mxRun(mxMod1)
        estimatesOpenMx = fitMod1$output$estimate
        
        mxTime = proc.time() - tempTime2
        
    } else{
        estimatesOpenMx = NULL
        mxTime = NULL
    }
    
    outputMatrix = cbind("Partitioned" = unlist(estimatesUser[1:attr(estimatesUser, "npar")]),
                         "Non-Partitioned" = unlist(estimatesNoPart[1:attr(estimatesNoPart, "npar")]),
                         "OpenMX-Non-Partitioned" = estimatesOpenMx)
    
    return(list("results" = outputMatrix, "Time" = c("Partitioned" = outTime, "OpenMx" = mxTime), 
                "Hessian" = list(attr(estimatesUser, "details")["L-BFGS-B" ,"nhatend"], 
                                 attr(estimatesNoPart, "details")["L-BFGS-B" ,"nhatend"], fitMod1$output$hessian)))
    
}



