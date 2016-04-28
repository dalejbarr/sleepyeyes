library("parallel")

## start the computing cluster
## If you want to make a fancier cluster to speed up the analysis, change the call to =makeCluster()= below.
## for example, this was the setting for my lab to use multiple cores on multiple machines:
## cl <- makeCluster(rep(c("localhost", "gossip", "yap", "chatter"),
##                      c(6, 6, 8, 8)))
cl <- makeCluster(detectCores())
