#####
## run
#####

source('~/Box Sync/Firewall Project/Code/matrixInputs.R')

genRAMSimulationMatrices = function(p){
    matrA = mxMatrix("Full", nrow = (p + 3), ncol = (p + 3),
                     values= rbind(c(rep(0.0, 2), 0.5, rep(0.0, p)), 
                                   c(rep(0.0, 2), 0.5, rep(0.0, p)),rep(0.0, (p+3)),
                                   cbind(rep(1.0, p), seq(0, (p-1)), rep(0, p),diag(0, p))), 
                     free=c(rep(F, 2), T, rep(F, (p+2)), T, rep(F, (p+(p+2)*(p+3)))),
                     labels=c(rep(NA, 2), "beta_I", rep(NA, (p+2)), "beta_S", rep(NA, (p+(p+2)*(p+3)))),
                     byrow=TRUE, name="A")
    
    matrS = mxMatrix("Full", nrow = (p + 3), ncol = (p + 3),
                     values = rbind(c(0.5, 0.2, rep(0.0, (p+1))),
                                    c(0.2, 0.5, rep(0.0, (p+1))),
                                    cbind(rep(0.0, (p + 1)), rep(0.0, (p + 1)), diag(1.0, (p+1)))), 
                     free = rbind(c(rep(T, 2), rep(F, (p + 1))),
                                  c(rep(T, 2), rep(F, (p + 1))),
                                  cbind(rep(F, (p + 1)), rep(F, (p + 1)), 
                                        matrix(as.logical(diag(1.0, (p+1))), nrow = (p+1)))),
                     labels = rbind(c("sigma_i","COV_icept_slope", rep(NA, (p + 1))),
                                    c("COV_icept_slope","sigma_s", rep(NA, (p + 1))),
                                    cbind(rep(NA, (p + 1)), rep(NA, (p + 1)), 
                                          matrix(c("VAR_x0", rep(c(rep(NA, (p + 1)), "e"), (p))), 
                                                 nrow = (p + 1), ncol = (p + 1), byrow = T))),
                     lbound = rbind(c(rep(0.0, 2), rep(NA, (p + 1))),
                                    c(rep(0.0, 2), rep(NA, (p + 1))),
                                    cbind(rep(NA, (p + 1)), rep(NA, (p + 1)), 
                                          matrix(c(rep(c(0.0, rep(NA, (p + 1))), (p)), 0.0), 
                                                 nrow = (p + 1), ncol = (p + 1), byrow = T))), 
                     byrow = TRUE, name = "S")
    
    matrF = mxMatrix("Full", nrow=(p+1),ncol=(p+3), 
                     values=cbind(rep(0, (p+1)), rep(0, (p+1)), diag(1, (p+1))), 
                     free=F,labels=NA, byrow=TRUE, name="F")
    
    matrM = mxMatrix("Full", nrow=1,ncol=(p+3), values=rep(0.0, (p+3)), 
                     free=rep(F, (p+3)),labels=rep(NA, (p+3)), 
                     byrow=TRUE, name="M")
    
    matrI = mxMatrix(type = "Iden", nrow = (p+3), ncol = (p+3), name = "I")
    
    matHold = list(matrA, matrS, matrF, matrM, matrI)
    return(matHold)

}

#test = genRAMSimulationMatrices(10)
#test

#p = 10

#simData = newGenData(test, 50)
#summary(simData$data)
#simData$cov
#simData$mean

#partition = c(rep(1, 5), rep(2, 6))

#dataList = list(data1 = simData$data[, c(1:5), drop = FALSE],
#                data2 = simData$data[, c(6:11), drop = FALSE])

#fitTest = newWrapFunc(test, partition, dataList, compareWithOpenMx = T, noPartData = simData)
#fitTest

## run 1 - vary n and k
p = 30 
k = c(1, 2, 10, 30)
n = c(100, 1000, 5000, 10000)
runTime1 = vector("list", length(k)*length(n))
runResult1 = vector("list", length(k)*length(n))

for(b in 1:length(n)){
    simMat = genRAMSimulationMatrices(p)
    simData = newGenData(simMat, n[b])
    for(a in 1:length(k)){
        
        partition = rep(1:k[a], each = (p/k[a]))
        dataList = vector("list", k[a])
        for(c in 1:length(dataList)){
            dataList[[c]] = simData$data[, c( (1 + (c-1)*(p/k[a])) : (c*(p/k[a])) ), drop = FALSE]
        }
        
        fitTest = newWrapFunc(simMat, partition, dataList, compareWithOpenMx = T, noPartData = simData)
        
        runTime1[[a + (b-1)*length(k)]] = fitTest$Time
        runResult1[[a + (b-1)*length(k)]] = fitTest$results
        
        cat(a)
    }
    cat(b, "\n")
}


## run 2 - vary p and k
p = c(10, 30, 60, 150, 300, 600)
k = c(1, 2, 10, 30)
n = 1000
runTime2 = vector("list", length(k)*length(p))
runResult2 = vector("list", length(k)*length(p))

for(b in 1:length(p)){
    simMat = genRAMSimulationMatrices(p[b])
    simData = newGenData(simMat, n)
    for(a in 1:length(k)){
        if((p[b] < k[a]) == TRUE){
            runTime2[[a + (b-1)*length(k)]] = NULL
            runResult2[[a + (b-1)*length(k)]] = NULL
        } else{
            partition = rep(1:k[a], each = (p[b]/k[a]))
            dataList = vector("list", k[a])
            for(c in 1:length(dataList)){
                dataList[[c]] = simData$data[, c( (1 + (c-1)*(p[b]/k[a])) : (c*(p[b]/k[a])) ), 
                                             drop = FALSE]
            }
            
            fitTest = newWrapFunc(simMat, partition, dataList, compareWithOpenMx = T, 
                                  noPartData = simData)
            
            runTime2[[a + (b-1)*length(k)]] = fitTest$Time
            runResult2[[a + (b-1)*length(k)]] = fitTest$results
        }
        cat(a)
    }
    cat(b, "\n")
}


## run 3 - vary k
p = 150
k = c(1, 2, 5, 10, 30, 75, 150)
n = 1000
runTime3 = vector("list", length(k))
runResult3 = vector("list", length(k))

simMat = genRAMSimulationMatrices(p)
simData = newGenData(simMat, n)
for(a in 1:length(k)){
    partition = rep(1:k[a], each = (p/k[a]))
    dataList = vector("list", k[a])
    for(c in 1:length(dataList)){
        dataList[[c]] = simData$data[, c( (1 + (c-1)*(p/k[a])) : (c*(p/k[a])) ), 
                                     drop = FALSE]
    }
    
    fitTest = newWrapFunc(simMat, partition, dataList, compareWithOpenMx = T, 
                          noPartData = simData)
    
    runTime3[[a]] = fitTest$Time
    runResult3[[a]] = fitTest$results
    
    cat(a, "\n")
}


#####
## makes sense
#####
load("~/Box Sync/Firewall Project/Code/PSD_simulations/varyK.RData")
load("~/Box Sync/Firewall Project/Code/PSD_simulations/varyKN.RData")
load("~/Box Sync/Firewall Project/Code/PSD_simulations/varyKP.RData")

kSimOut1
kSimOut2
kSimOut3

#kSimOut1[1, 1]$result.1$var
#kSimOut1[1, 2]$result.1$var
#kSimOut1[2, 1]$result.2$var
#kSimOut1[1, 17]$result.1$var

runTime1 = data.frame(matrix(NA, nrow = 160, ncol = 5))
colnames(runTime1) = c("time", "MAE", "p", "k", "n")

for(a in 1:4){
    for(b in 1:40){
        if(!is.na(kSimOut1[a, b][[1]][[2]][1])){
            runTime1[(b + (a - 1)*40), 1] = kSimOut1[a, b][[1]][[1]][3]
            runTime1[(b + (a - 1)*40), 3:5] = kSimOut1[a, b][[1]][[3]]
            
            runTime1[(b + (a - 1)*40), 2] = mean(abs(kSimOut1[a, b][[1]][[2]][,1] - 
                                                          kSimOut1[a, b][[1]][[2]][,2]))
        }
    }
}

runTime1$p = as.factor(runTime1$p)
runTime1$k = as.factor(runTime1$k)
#View(runTime1)
#summary(runTime1)
#table(runTime1$k, runTime1$n)
runTime1 = na.omit(runTime1)


runTime2 = data.frame(matrix(NA, nrow = 160, ncol = 5))
colnames(runTime2) = c("time", "MAE", "p", "k", "n")

for(a in 1:4){
    for(b in 1:40){
        if(!is.null(kSimOut2[a, b][[1]][[2]][1])){
            if(!is.na(kSimOut2[a, b][[1]][[2]][1])){
                runTime2[(b + (a - 1)*40), 1] = kSimOut2[a, b][[1]][[1]][3]
                runTime2[(b + (a - 1)*40), 3:5] = kSimOut2[a, b][[1]][[3]]
                
                runTime2[(b + (a - 1)*40), 2] = mean(abs(kSimOut2[a, b][[1]][[2]][,1] - 
                                                             kSimOut2[a, b][[1]][[2]][,2]))
            }
        }
    }
}

#runTime2$p = as.factor(runTime2$p)
runTime2$k = as.factor(runTime2$k)
#View(runTime2)
#summary(runTime2)
#table(runTime2$k, runTime2$p)
runTime2 = na.omit(runTime2)


runTime3 = data.frame(matrix(NA, nrow = 60, ncol = 5))
colnames(runTime3) = c("time", "MAE", "p", "k", "n")

for(a in 1:6){
    for(b in 1:10){
        if(!is.na(kSimOut3[a, b][[1]][[2]][1])){
            runTime3[(b + (a - 1)*10), 1] = kSimOut3[a, b][[1]][[1]][3]
            runTime3[(b + (a - 1)*10), 3:5] = kSimOut3[a, b][[1]][[3]]  
            
            runTime3[(b + (a - 1)*10), 2] = mean(abs(kSimOut3[a, b][[1]][[2]][,1] - 
                                                          kSimOut3[a, b][[1]][[2]][,2]))
        }
    }
}

runTime3$p = as.factor(runTime3$p)
#runTime3$k = as.factor(runTime3$k)
#View(runTime3)
#summary(runTime3)
runTime3 = na.omit(runTime3)



library(ggplot2)
library(grid)

varyKN = ggplot(data = runTime1, aes(x = log(n), y = log(time), color = k)) + 
    stat_summary(fun.y = "mean", geom = c("line")) + geom_point() + 
    ylab("Log Time") + xlab("Log Observations") + ggtitle("p = 30, k, n Varying")
#varyKN
varyKP = ggplot(data = runTime2, aes(x = log(p), y = log(time), color = k)) + 
    stat_summary(fun.y = "mean", geom = c("line")) + geom_point() + 
    ylab("Log Time") + xlab("Log Variables") + ggtitle("n = 500, k, p Varying")
#varyKP
varyK = ggplot(data = runTime3, aes(x = log(k), y = log(time), color = k)) + 
    stat_summary(fun.y = "mean", geom = c("line")) + geom_point() + 
    ylab("Log Time") + xlab("Log Partitions") + ggtitle("p = 150, n = 500, k Varying")
#varyK

errKN = ggplot(data = runTime1, aes(x = n, y = MAE, color = k)) + 
    stat_summary(fun.y = "mean", geom = c("line")) + geom_point() + 
    ylab("Paremeter MAE to OpenMx") + xlab("Observations") + ggtitle("p = 30, k, n Varying")
#errKN
errKP = ggplot(data = runTime2, aes(x = p, y = MAE, color = k)) + 
    stat_summary(fun.y = "mean", geom = c("line")) + geom_point() + 
    ylab("Paremeter MAE to OpenMx") + xlab("Variables") + ggtitle("n = 500, k, p Varying")
#errKP
errK = ggplot(data = runTime3, aes(x = k, y = MAE, color = k)) + 
    stat_summary(fun.y = "mean", geom = c("line")) + geom_point() + 
    ylab("Paremeter MAE to OpenMx") + xlab("Partitions") + ggtitle("p = 30, k, n Varying")
#errK


grid.newpage()
pushViewport(viewport(layout = grid.layout(3, 2)))
print(varyKN, vp = viewport(layout.pos.row = 1, layout.pos.col = 1) )
print(errKN, vp = viewport(layout.pos.row = 1, layout.pos.col = 2) )
print(varyKP, vp = viewport(layout.pos.row = 2, layout.pos.col = 1) )
print(errKP, vp = viewport(layout.pos.row = 2, layout.pos.col = 2) )
print(varyK, vp = viewport(layout.pos.row = 3, layout.pos.col = 1) )
print(errK, vp = viewport(layout.pos.row = 3, layout.pos.col = 2) )


### get O notation

summary(lm(log(time) ~ log(n), data = runTime1))

summary(lm(log(time) ~ log(p), data = runTime2))








