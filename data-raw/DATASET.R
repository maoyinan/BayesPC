## code to prepare `DATASET` dataset goes here

nClusters <- 4
nSub <- 30*nClusters # subjects
nT <- 40 # time points per subject
dat <- data.frame(
  ID=rep(seq(nSub),each=nT),
  Group=rep(seq(nClusters),each=nSub/nClusters*nT),
  Time=rep(seq(nT),nSub)
) %>%
  mutate(t= (Time-min(Time))/diff(range((Time))) )
dat$Z1 <- 1
for(i in c(1:9)){
  dat[,paste0('Z',i+1)] <- cos(pi*i*dat$t)
}
nBasis <- 10

# simulate
set.seed(12345)
epsilon <- rnorm(nSub*nT,0,1/10)
i1wi <- sample(c(1:3),nSub,replace = T)+1
i2wi <- sample(c(7:9),nSub,replace = T)+1

dat$Record <- NA
for(j in 1:nSub){# filter ID==j
  idx <- which(dat$ID==j)
  gr <- dat$Group[idx[1]] # which cluster determines beta
  beta1 <- ifelse(gr%in%c(1,2),1,.1)
  beta2 <- ifelse(gr%in%c(1,3),1,.1)
  z1 <- dat[idx,paste0('Z',i1wi[j])]
  z2 <- dat[idx,paste0('Z',i2wi[j])]
  yi <- beta1*z1+beta2*z2 + epsilon[idx]
  dat$Record[idx] <- yi
}
dat$Record <- as.numeric(scale(dat$Record))


DATASET <- dat
use_data(DATASET, overwrite = TRUE)
