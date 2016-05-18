library(foreign)

cannabis = read.spss("NU321_cannabis_psychosis.sav", to.data.frame = T)
cannabisNoMissing = na.omit(cannabis)
colnames(cannabis)
summary(cannabis$stem_ca)
summary(cannabis$dsm5_CUD)
summary(cannabis$sx_psyc)

cannabisSimple = cannabis[, c("stem_ca", "dsm5_CUD", "sx_psyc")]
cannabisSimpleNoMissing = na.omit(cannabisSimple)

cannabis$familyid <- as.numeric(as.character(cannabis$familyid))
cannabis$twinid <- as.numeric(as.character(cannabis$twinid))

# ----------------------------------------------------------------
# MxFactor ordinal data
# ----------------------------------------------------------------
library(OpenMx)

# Lifetime cannabis use frequency
match(c("stem_ca"),names(cannabis)) # 14
cannabis[,c(14)]<-mxFactor(cannabis[,c(14)], levels=c(0:2))

# DSM-V cannabis use disorder
match(c("dsm5_CUD"),names(cannabis)) # 15
cannabis[,c(15)]<-mxFactor(cannabis[,c(15)], levels=c(0:3))

# WHO CIDI psychosis screening items sum score
match(c("sx_psyc"),names(cannabis)) # 37
cannabis[,c(37)]<-mxFactor(cannabis[,c(37)], levels=c(0:2))

# WHO CIDI mania screening items sum score
match(c("sx_mania"),names(cannabis)) # 36
cannabis[,c(36)]<-mxFactor(cannabis[,c(36)], levels=c(0:3))

cannabis 	<- cannabis[ which(cannabis$zygosity<=6 & cannabis$twinid<=2), ]
data1		<- cannabis[,c("familyid", "twinid", "sex", "zygosity", "A47C_AGE","stem_ca",
                        "dsm5_CUD","sx_psyc")]
summary(data1)

# A47C_AGE = age of cannabis initiation

twindat <- function(dat, famid, twinid, zygosity) {
    datA <- dat[dat[,twinid]==min(dat[,twinid]),]    #twin1
    datB <- dat[dat[,twinid]==max(dat[,twinid]),]    #twin2 
    DAT <- merge(datA, datB, by=famid, all.x=TRUE, all.y=TRUE, suffixes=c("_T1","_T2"))
    DAT[,paste(twinid,"_T1",sep="")] <- NULL
    DAT[,paste(twinid,"_T2",sep="")] <- NULL
    DAT[,zygosity] <- ifelse(is.na(DAT[,paste(zygosity,"_T1",sep="")]),DAT[,paste(zygosity,"_T2",sep="")],DAT[,paste(zygosity,"_T1",sep="")])
    DAT[,paste(zygosity,"_T1",sep="")] <- NULL  
    DAT[,paste(zygosity,"_T2",sep="")] <- NULL  
    return(DAT)
}
twin_data <- twindat(dat=data1, famid="familyid", twinid= "twinid", zygosity= "zygosity")
#describe(twin_data)
summary(twin_data)

#####
## model 1
#####
A = mxMatrix("Full", nrow=6,ncol=6,values=c(
    0.0,0.0,0.0,0.0,0.0,0.0,
    0.0,0.0,0.0,0.0,0.0,0.0,
    0.0,0.0,0.0,0.0,0.0,0.0,
    1.0,0.0,0.0,0.0,0.0,0.0,
    0.0,1.0,0.0,1.0,0.0,0.0,
    1.0,1.0,1.0,0.0,0.0,0.0), 
    free=c(
    F,F,F,F,F,F,
    F,F,F,F,F,F,
    F,F,F,F,F,F,
    T,F,F,F,F,F,
    F,T,F,T,F,F,
    T,T,T,F,F,F),
    labels=c(
    NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA,
    "a11",NA,NA,NA,NA,NA,
    NA,"a22",NA,"b21",NA,NA,
    "a31","a32","a33",NA,NA,NA), byrow=TRUE, name="A")

S = mxMatrix("Full", nrow=6,ncol=6,values=c(
    1.0,0.0,0.0,0.0,0.0,0.0,
    0.0,1.0,0.0,0.0,0.0,0.0,
    0.0,0.0,1.0,0.0,0.0,0.0,
    0.0,0.0,0.0,0.0,0.0,0.0,
    0.0,0.0,0.0,0.0,0.0,0.0,
    0.0,0.0,0.0,0.0,0.0,0.0), 
    free=c(
    F,F,F,F,F,F,
    F,F,F,F,F,F,
    F,F,F,F,F,F,
    F,F,F,F,F,F,
    F,F,F,F,F,F,
    F,F,F,F,F,F),
    labels=c(
    NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA),
    lbound=c(
        0.0,NA,NA,NA,NA,NA,
        NA,0.0,NA,NA,NA,NA,
        NA,NA,0.0,NA,NA,NA,
        NA,NA,NA,NA,NA,NA,
        NA,NA,NA,NA,NA,NA,
        NA,NA,NA,NA,NA,NA),
    byrow=TRUE, name="S")

Fm = mxMatrix("Full", nrow=3,ncol=6, values=c(
    0.0,0.0,0.0,1.0,0.0,0.0,
    0.0,0.0,0.0,0.0,1.0,0.0,
    0.0,0.0,0.0,0.0,0.0,1.0), 
    free=c(
    F,F,F,F,F,F,
    F,F,F,F,F,F,
    F,F,F,F,F,F),
    labels=c(
    NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA), byrow=TRUE, name="F")

M = mxMatrix("Full", nrow=1,ncol=6, values=c(
    0.0,0.0,0.0,0.0,0.0,0.0), 
    free=c(
    F,F,F,F,F,F),
    labels=c(
    NA,NA,NA,NA,NA,NA), byrow=TRUE, name="M")

I = mxMatrix("Iden", nrow = 6, ncol = 6, name = "I")

simMat = list(A, S, Fm, M, I)

partition = c(1, 2, 3)
dataList = list(cannabisSimple[, c(1), drop = F], cannabisSimple[, c(2), drop = F],
                cannabisSimple[, c(3), drop = F])
dataList = list(cannabisSimpleNoMissing[, c(1), drop = F], cannabisSimpleNoMissing[, c(2), drop = F],
                cannabisSimpleNoMissing[, c(3), drop = F])

model1 = newWrapFunc(simMat, partition, dataList, compareWithNoPart = T,
            compareWithOpenMx = T, noPartData = list(cannabisSimpleNoMissing))


#####
## model 2
#####
A = mxMatrix("Full", nrow=6,ncol=6,values=c(
    0.0,0.0,0.0,0.0,0.0,0.0,
    0.0,0.0,0.0,0.0,0.0,0.0,
    0.0,0.0,0.0,0.0,0.0,0.0,
    1.0,0.0,0.0,0.0,0.0,0.0,
    0.0,1.0,0.0,1.0,0.0,0.0,
    0.0,0.0,1.0,1.0,1.0,0.0), free=c(
    F,F,F,F,F,F,
    F,F,F,F,F,F,
    F,F,F,F,F,F,
    T,F,F,F,F,F,
    F,T,F,T,F,F,
    F,F,T,T,T,F),
    labels=c(
    NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA,
    "a11",NA,NA,NA,NA,NA,
    NA,"a22",NA,"b21",NA,NA,
    NA,NA,"a33","b31","b32",NA), byrow=TRUE, name="A")

S = mxMatrix("Full", nrow=6,ncol=6,values=c(
    1.0,0.0,0.0,0.0,0.0,0.0,
    0.0,1.0,0.0,0.0,0.0,0.0,
    0.0,0.0,1.0,0.0,0.0,0.0,
    0.0,0.0,0.0,0.0,0.0,0.0,
    0.0,0.0,0.0,0.0,0.0,0.0,
    0.0,0.0,0.0,0.0,0.0,0.0), 
    free=c(
    F,F,F,F,F,F,
    F,F,F,F,F,F,
    F,F,F,F,F,F,
    F,F,F,F,F,F,
    F,F,F,F,F,F,
    F,F,F,F,F,F),
    labels=c(
    NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA), byrow=TRUE, name="S")

Fm = mxMatrix("Full", nrow=3,ncol=6, values=c(
    0.0,0.0,0.0,1.0,0.0,0.0,
    0.0,0.0,0.0,0.0,1.0,0.0,
    0.0,0.0,0.0,0.0,0.0,1.0), 
    free=c(
    F,F,F,F,F,F,
    F,F,F,F,F,F,
    F,F,F,F,F,F),
    labels=c(
    NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA,
    NA,NA,NA,NA,NA,NA), byrow=TRUE, name="F")

M = mxMatrix("Full", nrow=1,ncol=6, values=c(
    0.0,0.0,0.0,0.0,0.0,0.0), 
    free=c(
    F,F,F,F,F,F),
    labels=c(
    NA,NA,NA,NA,NA,NA), byrow=TRUE, name="M")

I = mxMatrix("Iden", nrow = 6, ncol = 6, name = "I")

simMat = list(A, S, Fm, M, I)

partition = c(1, 2, 3)
dataList = list(cannabisSimpleNoMissing[, c(1), drop = F], cannabisSimpleNoMissing[, c(2), drop = F],
                cannabisSimpleNoMissing[, c(3), drop = F])

model2 = newWrapFunc(simMat, partition, dataList, compareWithNoPart = T,
            compareWithOpenMx = T, noPartData = list(cannabisSimpleNoMissing))

model2


#####
## print results
#####
model1$results
model2$results

library(xtable)
psdTable = xtable(rbind(model1$results, model2$results), digits = rep(6, 4))
psdTable



#####
## fix bug
#####
ramMatr = simMat
for(i in 1:length(ramMatr)){
    names(ramMatr)[i] = ramMatr[[i]]$name
    if(ramMatr[[i]]$name == "F")
        names(ramMatr)[i] = ramMatr[[i]]$name = "Fm"
}
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
matReOrder = rep(NA, length(unique(partition)))
tempLen = c(1, 0)
for(i in 1:length(unique(partition))){
    tempPart = which(partition == i)
    tempLen[2] = tempLen[2] + length(tempPart)
    matReOrder[tempLen[1]:tempLen[2]] = tempPart
    tempLen[1] = 1 + tempLen[2]
}


newFitFunc(freeValues, dataList, partition, ramMatr, matReOrder)








