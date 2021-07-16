# Helper functions to calculate intermediate steps for projection-based clustering
globalVariables(c("%>%", "ID", "idx", "y"))

Q.func <- function(idxA, idxB, dat,Gamma, G, id){
  nSub= length(unique(id))

  GA <- G[idxA,idxA]
  GA1 <- solve(GA)
  if(length(idxB)>0){
    GB <- G[idxB,idxB]
    GAB <- G[idxA,idxB]
    ZB <- as.matrix(dat[,paste0('Z',idxB)])
  }
  ZA <- as.matrix(dat[,paste0('Z',idxA)])

  Q <- array(NA, dim=c(nSub, length(idxA), length(idxA)))
  for(i in seq(nSub)){
    idx <- which(id==i)
    Gammai <- Gamma*diag(1,nrow=length(idx),ncol=length(idx))
    ZiA <- ZA[idx,]
    if(length(idxB)>0){
      ZiB <- ZB[idx,]
      term1 <- ZiA+ZiB%*%t(GAB)%*%GA1
      term2 <- ZiB%*%(GB-t(GAB)%*%solve(GA)%*%GAB)%*%t(ZiB) + Gammai
      Q[i,,] <- t(term1) %*% solve(term2) %*% term1

    }
    else{
      Q[i,,] <- Gamma*diag(1,nrow=length(idxA),ncol=length(idxA))

    }

  }
  Q
}

# k means algorithm -------------------------------------------------------

dK.func0 <- function(nCluster, clusterVec, Q, bAMatrix,regQ){
  dK <- matrix(NA, nrow=nCluster, ncol=ncol(bAMatrix))
  for(j in seq(nCluster)){
    ind <- which(clusterVec==j)
    if(length(ind)>1){
      term1 <- apply(Q[ind,,], 2:3, sum)
      term2 <- apply(sapply(ind, function(i)Q[i,,]%*%bAMatrix[i,]),1,sum)
    } else {
      term1 <- Q[ind,,]
      term2 <- Q[ind,,]%*%bAMatrix[ind,]
    }
    term11 <-
      tryCatch(expr = solve(term1),
               error = function(c) {
                 # matlib::Ginv(term1, tol=regQ)
                 solve(term1+diag(regQ,nrow(term1),ncol(term1)))
               })

    dK[j,] <- term11%*%term2
  }
  dK
}
dK.func1 <- function(nCluster, clusterVec, ls_Q, ls_bAMatrix, regQ){
  # used for summing over posterior samples
  nDraw <- dim(ls_bAMatrix)[3]
  sum_Q <- apply(ls_Q, 1:3, sum)
  sum_Qb <- array(0, dim=dim(ls_Q)[1:2])
  for(k in seq(dim(ls_Q)[1])){
    for(i in seq(nDraw))
      sum_Qb[k,] <- sum_Qb[k,] + c(ls_Q[k,,,i] %*% ls_bAMatrix[k,,i])
  }

  dK <- matrix(NA, nrow=nCluster, ncol=dim(ls_Q)[2])
  for(j in seq(nCluster)){
    ind <- which(clusterVec==j)
    if(length(ind)>1){
      term1 <- apply(sum_Q[ind,,], 2:3, sum)
      term2 <- apply(sapply(ind, function(i)sum_Qb[i,]),1,sum)
    } else {
      term1 <- sum_Q[ind,,]
      term2 <- sum_Qb[ind,]
    }
    term11 <-
      tryCatch(expr = solve(term1),
               error = function(c) {
                 # matlib::Ginv(term1, tol=regQ)
                 solve(term1+diag(regQ,nrow(term1),ncol(term1)))
               })

    dK[j,] <- term11%*%term2
  }
  dK
}

KL.func <- function(diA, biA, Qi){
  1/2*t(diA-biA)%*%Qi%*%(diA-biA)
}

C.func0 <- function(nCluster, dK, Q, bAMatrix){
  clusterVec1 <- numeric(nrow(bAMatrix))
  termMatrix <- matrix(NA,nrow=nrow(bAMatrix),ncol=nCluster)
  for(i in seq_along(clusterVec1)){
    biA <- bAMatrix[i,]
    term <- numeric(nCluster)
    for(j in seq(nCluster)){
      term[j] <- KL.func(dK[j,],biA, Q[i,,])
      termMatrix[i,j] <- term[j]
    }
    clusterVec1[i] <- which(term==min(term))
  }
  clusterVec1
}
C.func1 <- function(nCluster, dK, ls_Q, ls_bAMatrix){
  nSub <- dim(ls_Q)[1]
  nDraw <- dim(ls_bAMatrix)[3]
  clusterVec1 <- numeric(nSub)
  termMatrix <- matrix(NA,nrow=nSub,ncol=nCluster)
  for(i in seq(nSub)){
    term <- numeric(nCluster)
    for(j in seq(nCluster)){
      for(k in seq(nDraw))
        term[j] <-  term[j] + KL.func(dK[j,],ls_bAMatrix[i,,k], ls_Q[i,,,k])

      # term[j] <- KL.func(dK[j,],biA, Q[i,,])
      termMatrix[i,j] <- term[j]
    }
    clusterVec1[i] <- which(term==min(term))
  }
  clusterVec1
}

cluster.initialize <- function(clusterVec0,nCluster,Q,bAMatrix,regQ, seed){
  clusterVec <- clusterVec0
  temp <- which(table(clusterVec0)>1)
  if(length(temp)==1)  {
    iCluster <- temp
    # split the first cluster of previous iteration (who has more than 1 points)
    idx <- which(clusterVec0==iCluster)
    split <- cluster.initialize0(2,Q[idx,,], bAMatrix[idx,],seed)
    clusterVec[idx[split==2]] <- nCluster
  }
  else {
    clusterMat <- matrix(NA, nrow=length(clusterVec),ncol=length(temp))
    kl <- numeric(length(temp))
    for(j in seq_along(kl)){
      iCluster <- temp[j]

      # split the first cluster of previous iteration (who has more than 1 points)
      idx <- which(clusterVec0==iCluster)
      split <- cluster.initialize0(2,Q[idx,,], bAMatrix[idx,],seed)
      clusterVec[idx[split==2]] <- nCluster
      clusterMat[,j] <- clusterVec
      dK <- dK.func0(nCluster, clusterVec, Q, bAMatrix,regQ)
      kl[j] <- KLK.func(dK,clusterVec,bAMatrix,Q)
    }
    ii <- which(kl==min(kl))
    if(length(ii)==1) clusterVec <- clusterMat[,ii]
    else clusterVec <- clusterMat[,sample(ii,1)]

  }
  clusterVec
}
cluster.initialize0 <- function(nCluster,Q,bAMatrix, seed){
  set.seed(seed)
  dK <- bAMatrix[sort(sample(seq(nrow(bAMatrix)), nCluster, replace = F)),]
  clusterVec0 <- C.func0(nCluster,dK, Q, bAMatrix)
  clusterVec0
}

cluster.iter <- function(clusterVec0,nIter,nCluster,Q,bAMatrix,regQ, seed, do.plot){
  ls_Q <- Q
  ls_bAMatrix <- bAMatrix
  if(length(dim(Q))==4){
    Q <- ls_Q[,,,1]
    bAMatrix <- ls_bAMatrix[,,1]
    dK.func <- dK.func1
    C.func <- C.func1
  }
  else{
    dK.func <- dK.func0
    C.func <- C.func0
  }

  clusterMatrix <- matrix(NA, nrow= nrow(bAMatrix), ncol= nIter)
  dKArray <- array(NA, dim=c(nCluster, ncol(bAMatrix), nIter))
  if(is.null(clusterVec0))
    clusterVec <- cluster.initialize0(nCluster,Q,bAMatrix, seed)
  else clusterVec <- cluster.initialize(clusterVec0,nCluster, Q,bAMatrix,regQ, seed)

  dK <- dK.func(nCluster, clusterVec, ls_Q, ls_bAMatrix,regQ)
  clusterMatrix[,1] <- clusterVec
  dKArray[,,1] <- dK

  if(do.plot){
    plot(bAMatrix[,1],bAMatrix[,2], col=clusterVec, main=sprintf('Iter 1 with %d clusters',nCluster))
    text(dK[,1],dK[,2],seq(nCluster), col=seq(nCluster), cex=2)
  }

  for(i in seq(2,nIter)){
    #update
    clusterVec1 <- C.func(nCluster,dK, ls_Q, ls_bAMatrix)
    # if(sum(clusterVec1-clusterMatrix[,i-1])==0) break
    if(sum(clusterVec1-clusterMatrix[,i-1])==0 | length(unique(clusterVec1))<nCluster) break

    clusterVec <- clusterVec1
    dK <- dK.func(nCluster, clusterVec, ls_Q, ls_bAMatrix, regQ)
    clusterMatrix[,i] <- clusterVec
    dKArray[,,i] <- dK

    if(do.plot){
      plot(bAMatrix[,1],bAMatrix[,2], col=clusterVec, main=sprintf('Iter %d with %d clusters',i,nCluster))
      text(dK[,1],dK[,2],seq(nCluster), col=seq(nCluster), cex=2)
    }
  }
  clusterMatrix[,nIter] <- clusterVec
  dKArray[,,nIter] <- dK
  list(clusterMatrix,dKArray)
}

KLK.func <- function(dK,clusterVec,bAMatrix,Q){
  nSub <- nrow(bAMatrix)
  KLKvec <- numeric(nSub)
  for(i in seq(nSub)){
    KLKvec[i] <- KL.func(dK[clusterVec[i],], bAMatrix[i,], Q[i,,])
  }
  KLK <- sum(KLKvec)
  KLK
}

