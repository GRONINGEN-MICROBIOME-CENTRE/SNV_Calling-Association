getClusNumStability <- function(subsampleProportions, nIterClusStability = 10, distObj,psCut) {
  
  subsampIters <- sort(unlist(rep(subsampleProportions,nIterClusStability)))
  subsampNClusters <- sapply(subsampIters, nClusWithSubsample, distObj=distObj,psCut=psCut) # this is slow
  
  names(subsampNClusters) <- subsampIters
  nClusStability <- cbind.data.frame(propSamples = subsampIters,numClusters = subsampNClusters)
  
  return(nClusStability)
}

getClusNumStability2 <- function(subsampleProportions, nIterClusStability = 10, distObj,psCut, threads=1) {
  
  subsampIters <- sort(unlist(rep(subsampleProportions,nIterClusStability)))
  doParallel::registerDoParallel(threads)
  subsampNClusters <- foreach( x = seq(subsampIters)) %dopar% {
      nClusWithSubsample( subsampIters[x] ,  distObj=distObj,psCut=psCut )
  }
  subsampNClusters = unlist(subsampNClusters)
  #sapply(subsampIters, nClusWithSubsample, distObj=distObj,psCut=psCut) # this is slow
  
  names(subsampNClusters) <- subsampIters
  nClusStability <- tibble(propSamples = subsampIters,numClusters = subsampNClusters)
  
  return(nClusStability)
}

nClusWithSubsample <- function(distObj,subsampleProp,psCut){
  nSamples <- length(labels(distObj))
  indx <- sample(x = 1:nSamples,size = floor(nSamples*subsampleProp),replace = F)
  distDistinctSub <-as.dist(as.matrix(distObj)[indx,indx])
  res <- getClusPredStrengthResult(distDistinctSub,psCut=psCut, warn = FALSE)
  numClusters <- res[["optimalk"]]
  return(numClusters)
}

getClusNumStabilityPlots <- function(nClusStability,defaultGMax=10){
  
  nClusStabilityMean <- nClusStability %>%
    group_by(propSamples) %>%
    summarise(meanNumClusters=mean(as.numeric(as.character(numClusters))),
              medianNumClusters=median(as.numeric(as.character(numClusters))),
              nIter = n()) %>%
    ungroup()
  nClusStability$propSamples <- factor(nClusStability$propSamples)
  
  defaultGMax <- max(defaultGMax, max(nClusStability$numClusters))
  nClusStabilityNumIter <- max(nClusStabilityMean$nIter)
  
  # p1 <- ggplot(nClusStability,aes(x=propSamples))+
  #   geom_boxplot(aes(y=numClusters,group=propSamples),fill=NA)+
  #   xlab("Proportion of samples used \nto predict number of clusters")+
  #   ylab(paste0("Number of clusters predicted (",
  #               nClusStabilityNumIter," iterations)"))+
  #   scale_y_continuous(breaks = seq(1,defaultGMax+0.5,1),
  #                      limits = c(0.5,defaultGMax+0.5))
  
  p1b <- ggplot(nClusStability,aes(x=propSamples))+
    geom_jitter(aes(y=numClusters,group=propSamples),
                width = 0.2, height = 0.2,color="grey50",alpha=0.5)+
    geom_boxplot(aes(y=numClusters,group=propSamples),fill=NA,outlier.colour = NA)+
    geom_point(data=nClusStabilityMean,
               aes(x=factor(propSamples),y=meanNumClusters),
               color="darkred",shape=7,inherit.aes = F)+
    geom_line(data=nClusStabilityMean, group=1,
              aes(x=factor(propSamples),y=meanNumClusters),
              color="darkred",inherit.aes = F)+
    xlab("Proportion of samples used")+
    ylab(paste0("Number of clusters predicted",
                #" (",nClusStabilityNumIter," iterations)",
                "\nRed is mean number of clusters predicted"
    ))+
    scale_y_continuous(breaks = seq(1,defaultGMax+0.5,1),
                       limits = c(0.5,defaultGMax+0.5))
  
  
  # p2 <- nClusStability %>% group_by(propSamples,numClusters) %>% summarise(count=n()) %>%
  #   ggplot(aes(x=propSamples,y=numClusters,group=propSamples,fill=count,label=count))+
  #   geom_tile()+geom_text(color="grey50")+
  #   xlab("Proportion of samples used \nto predict number of clusters")+
  #   ylab("Number of clusters predicted")+
  #   scale_fill_continuous(paste0("Number of times \nprediction was \nmade out of ",
  #                                nClusStabilityNumIter,"\niterations"))+
  #   scale_y_continuous(breaks = seq(1,defaultGMax,1),
  #                      limits = c(0.05,defaultGMax+0.5)) # need 0.5 for block to show
  
  nClusStabilityCounts <- nClusStability %>%
    group_by(propSamples,numClusters) %>% summarise(count=n()) %>% ungroup() %>%
    mutate(numClusters = factor(numClusters,levels = seq(from = 1,to = defaultGMax,by = 1)),
           propSamples = factor(propSamples))
  
  p3 <- nClusStability %>%
    mutate(numClusters = factor(numClusters,levels = seq(from = 1,to = defaultGMax,by = 1)),
           propSamples = factor(propSamples)) %>%
    ggplot(aes(x=propSamples))+
    scale_y_discrete(drop=F)+ # don't drop unused levels
    geom_tile(data = nClusStabilityCounts,
              aes(y=numClusters,fill=count,group=propSamples),width=0.5,height=0.5)+
    #geom_boxplot(aes(y=numClusters,group=propSamples),fill=NA,outlier.colour = NA,color="grey50")+
    # geom_point(data=nClusStabilityMean,
    #            aes(x=factor(propSamples),y=meanNumClusters),
    #            color="darkred",shape=7,inherit.aes = F)+
    # geom_line(data=nClusStabilityMean, group=1,
    #           aes(x=factor(propSamples),y=meanNumClusters),
    #           color="darkred",inherit.aes = F)+
    geom_text(data = nClusStabilityCounts,
              aes(y=numClusters,group=propSamples,label=count),color="grey20")+
    xlab("Proportion of samples used")+
    ylab(paste0("Number of clusters predicted"#,
                #" (",nClusStabilityNumIter," iterations)",
                #"\nRed is mean number of clusters predicted"
    ))+
    scale_fill_gradient(paste0("Number of times \nprediction was \nmade out of ",
                               nClusStabilityNumIter,"\niterations"),
                        limits=c(0,nClusStabilityNumIter),
                        breaks=seq(0,nClusStabilityNumIter,
                                   by = floor(nClusStabilityNumIter)/5),
                        low = "mistyrose",high = "lightcoral")+
    #scale_y_continuous(breaks = seq(1,defaultGMax+0.5,1),
    #                   limits = c(0.5,defaultGMax+0.5))+
    theme_minimal()
  
  # nClusStability %>% group_by(propSamples,numClusters) %>% summarise(count=n()) %>%
  #   ggplot(aes(x=factor(propSamples),y=numClusters,
  #              group=propSamples,label=count))+
  #   scale_size_area(max_size = 10)+
  #   geom_count(aes(size=count))+geom_text(color="grey50")+
  #   xlab("Proportion of samples used to predict number of clusters")+
  #   ylab("Number of clusters predicted")+
  #   scale_color_continuous(paste0("Number of times \nprediction was \nmade out of ",
  #                                nClusStabilityNumIter,"\niterations"))+
  #   scale_y_continuous(breaks = seq(1,max(nClusStability$numClusters),1),
  #                      limits = c(1,max(nClusStability$numClusters)+0.5))
  
  return(list("boxJitter"=p1b,"boxHeat"=p3))
}

# Assess cluster membership stability -------------

clusterAssignSubsample <- function(subsampleProp,distObj,numClusters){
  nSamples <- length(labels(distObj))
  clusteringBoot <- fpc::clusterboot(data = distObj,distances = T,
                                     k = numClusters,diss = TRUE,usepam = TRUE,
                                     bootmethod = "subset",
                                     subtuning = floor(nSamples*subsampleProp),
                                     count = F,showplots = F,
                                     clustermethod = claraCBI)
  
  #clustering <- clusteringBoot$result$result
  
  clusStab_jaccSubsetMean <- clusteringBoot$subsetmean
  clusStab_recover <- clusteringBoot$subsetrecover/clusteringBoot$subsetbrd
  
  return(data.frame(clusterID = 1:length(clusStab_jaccSubsetMean),
                    nSamplesInCluster =  as.vector(table(clusteringBoot$partition)),
                    subsampleProp = rep(subsampleProp,length(clusStab_jaccSubsetMean)),
                    clusterStabilityJaccardMean = clusStab_jaccSubsetMean,
                    clusterStabilityPropRecover = clusStab_recover))
}

getClusMembStabPlots <- function(rare){
  p1 <- ggplot(rare,aes(x=subsampleProp,y=clusterStabilityPropRecover,
                        color=paste0(clusterID," (n=",nSamplesInCluster,")"),
                        group=clusterID))+
    geom_point()+geom_line()+ylim(c(0,1))+
    scale_color_discrete("clusterID \n(n= # samples \nin cluster)")+
    xlab("Proportion of samples used in cluster assignment")+
    ylab("Cluster Stability:\nMean proportion of times cluster was recovered")
  
  p2 <- ggplot(rare,aes(x=subsampleProp,y=clusterStabilityJaccardMean,
                        color=paste0(clusterID," (n=",nSamplesInCluster,")"),
                        group=clusterID))+
    geom_point()+geom_line()+ylim(c(0,1))+
    scale_color_discrete("clusterID \n(n= # samples \nin cluster)")+
    xlab("Proportion of samples used in cluster assignment")+
    ylab("Cluster Stability:\nJaccard subsetting mean")
  
  return(list(p1,p2))
}

getClusMembStability <- function(subsampleProportions, numClusters,distObj) {
  rare <- lapply(subsampleProportions,clusterAssignSubsample,
                 distObj=distObj,numClusters=numClusters)
  rare <- do.call(rbind,rare)
  rare$clusterID <- factor(rare$clusterID)
  return(rare)
}
getClusMembStability2 <- function(subsampleProportions, numClusters,distObj, threads=1) {
  doParallel::registerDoParallel(threads)
  rare <- foreach( x = seq(subsampleProportions), .combine = "rbind" ) %dopar% {
    clusterAssignSubsample( subsampleProportions[[x]] ,  distObj=distObj,numClusters=numClusters )
  }
  #rare <- lapply(subsampleProportions,clusterAssignSubsample,
  #               distObj=distObj,numClusters=numClusters)
  #rare <- do.call(rbind,rare)
  rare$clusterID <- factor(rare$clusterID)
  return(rare)
}

# Summarise clustering stability assessment --------
summariseClusteringStability <- function(nClusStability, clusMembStability, numClusters){
  nClusStabScoreNum <- getNClusStabScore(nClusStability)
  
  clusMembStabScoresNum <- sapply(1:numClusters,getClusMembStabScore,
                                  clusMembStability=clusMembStability)
  clusMembStabScores <- sapply(clusMembStabScoresNum,getStabScoreString)
  names(clusMembStabScores) <- paste0("clust",1:numClusters)
  
  list(numClusStabScore=getStabScoreString(nClusStabScoreNum),
       clusMembStabScores=clusMembStabScores)
}

# can give clustering stability a rating of: high, medium, or low
getStabScoreString <- function(numScore){
  STABILITY.SCORE.LEVELS <- c("Low","Medium","High")
  return(STABILITY.SCORE.LEVELS[numScore])
}



getNClusStabScore <- function(nClusStability){
  nClusStability$propSamples <- round(nClusStability$propSamples,1)
  nClusStabScore <- 1 # 1 = low
  # Number of clusters
  # same number of clusters predicted in every iteration when using 100% of data?
  noVar100 <- var(nClusStability[nClusStability$propSamples == 1,"numClusters"]) == 0
  nClusStabScore <- nClusStabScore + noVar100
  
  # same number of clusters predicted in every iteration when using 80% & 90% of data?
  # and this == to 100%
  if(nClusStabScore>1){
    noVar80 <- var(nClusStability[nClusStability$propSamples == 0.8,"numClusters"]) == 0
    noVar90 <- var(nClusStability[nClusStability$propSamples == 0.9,"numClusters"]) == 0
    nC100 <- nClusStability[nClusStability$propSamples == 1,"numClusters"][1]
    nC90 <- nClusStability[nClusStability$propSamples == 0.9,"numClusters"][1]
    nC80 <- nClusStability[nClusStability$propSamples == 0.8,"numClusters"][1]
    c2 <- all(noVar80, noVar90, nC90 == nC100, nC80 == nC100 )
    nClusStabScore <- nClusStabScore + c2
  }
  return(nClusStabScore)
}

# Cluster membership stability for each cluster
getClusMembStabScore <- function(clustID,clusMembStability){
  clusMembStabScore <- 1 # 1 = low
  clusMembStability$subsampleProp <- round(clusMembStability$subsampleProp,1)
  # both stability measures > 0.8 when using 90% of data? (medium achieved)
  clusMembStabScore <- clusMembStabScore +
    (clusMembStability[clusMembStability$subsampleProp==0.9,"clusterStabilityPropRecover"][clustID] > 0.8 &
       clusMembStability[clusMembStability$subsampleProp==0.9,"clusterStabilityJaccardMean"][clustID] > 0.8 )
  
  # both stability measures > 0.9 when using 70% of data? (high achieved)
  clusMembStabScore <- clusMembStabScore +
    (clusMembStability[clusMembStability$subsampleProp==0.7,"clusterStabilityPropRecover"][clustID] > 0.9 &
       clusMembStability[clusMembStability$subsampleProp==0.7,"clusterStabilityJaccardMean"][clustID] > 0.9 )
  return(clusMembStabScore)
}

saveClustStabilityPlots <- function(outDir, filePrefix, clusNumStabilityPlots, clusMembStabilityPlots) {
  
  # png(filename = paste(outDir,"/",filePrefix,'_clusNumStability-box.png',sep=''),
  #     res = 200,width = 12,height = 10,units = "cm")
  # print(clusNumStabilityPlots[[1]])
  # dev.off()
  
  png(filename = paste(outDir,"/",filePrefix,'_clusNumStability-heatmap.png',sep=''),
      res = 200,width = 12,height = 10,units = "cm")
  print(clusNumStabilityPlots[[2]])
  dev.off()
  
  png(filename = paste(outDir,"/",filePrefix,'_clusMembStability-recover.png',sep=''),
      res = 200,width = 12,height = 10,units = "cm")
  print(clusMembStabilityPlots[[1]])
  dev.off()
  
  png(filename = paste(outDir,"/",filePrefix,'_clusMembStability-jaccard.png',sep=''),
      res = 200,width = 12,height = 10,units = "cm")
  print(clusMembStabilityPlots[[2]])
  dev.off()
}


getMaxNumClustersToTry <- function(distObj,warn=T,filePrefix = "",
                                   defaultMaxNumClusters=10,
                                   minClusterSize = 5) {
  Gmax <- defaultMaxNumClusters
  n <- nrow(as.matrix(distObj))
  nf <- c(floor(n*0.5), n - floor(n*0.5))
  maxNclus <- floor(n/minClusterSize)
  
  if(min(min(nf)-1,maxNclus) < Gmax){
    Gmax <- min( min(nf)-1 , maxNclus)
    if(warn){
      warning(paste0("After filtering, species ",filePrefix," only has ",n,
                     " samples to identify clusters from.",
                     " Only testing prediction strength up to ",Gmax,
                     " clusters instead of up to ",defaultMaxNumClusters))
    }
  }
  return(Gmax)
}

predStrengthCustom <- function (distance, Gmin = 2, Gmax = 10, M = 50,
                                classification = "centroid", cutoff = 0.8, nnk = 1,
                                clustermethod = cluster::pam, ...)
{
  dist <- as.matrix(distance)
  n <- nrow(dist)
  nf <- c(floor(n*0.5), n - floor(n*0.5))
  indvec <- clcenters <- clusterings <- jclusterings <- classifications <- list()
  prederr <- list()
  
  for (k in Gmin:Gmax) {
    prederr[[k]] <- numeric(0)
    for (l in 1:M) {
      nperm <- sample(n, n)
      indvec[[l]] <- list()
      indvec[[l]][[1]] <- nperm[1:nf[1]]
      indvec[[l]][[2]] <- nperm[(nf[1] + 1):n]
      for (i in 1:2) {
        if(identical(clustermethod, cluster::pam) ){
          clusterings[[i]] <- as.vector(cluster::pam(as.dist(dist[indvec[[l]][[i]],indvec[[l]][[i]]]), k, diss=TRUE))
        }else{
          clusterings[[i]] <- as.vector(clustermethod(as.dist(dist[indvec[[l]][[i]],indvec[[l]][[i]]]), k, diss=TRUE, ...))$result
        }
        
        jclusterings[[i]] <- rep(-1, n)
        jclusterings[[i]][indvec[[l]][[i]]] <- clusterings[[i]]$clustering
        centroids <- clusterings[[i]]$medoids
        j <- 3 - i
        classifications[[j]] <- fpc::classifdist(distance, jclusterings[[i]],
                                                 method = classification, centroids = centroids,
                                                 nnk = nnk)[indvec[[l]][[j]]]
      }
      
      ps_f <- matrix(0, nrow = 2, ncol = k)
      
      for (i in 1:2) {
        for (kk in 1:k) {
          nik <- sum(clusterings[[i]]$clustering == kk)
          if (nik > 1) {
            a <- which(clusterings[[i]]$clustering[1:(nf[i] - 1)] == kk)
            ps_f[i,kk] <- sum(outer(classifications[[i]][a],classifications[[i]][a],'=='))-length(a)
            ps_f[i,kk] <- ps_f[i, kk]/(nik * (nik - 1))
          }
        }
      }
      
      ps <- ps_f
      
      prederr[[k]][l] <- mean(c(min(ps[1, ]), min(ps[2,
      ])))
    }
  }
  mean.pred <- numeric(0)
  if (Gmin > 1)
    mean.pred <- c(1)
  if (Gmin > 2)
    mean.pred <- c(mean.pred, rep(NA, Gmin - 2))
  for (k in Gmin:Gmax) mean.pred <- c(mean.pred, mean(prederr[[k]]))
  optimalk <- max(which(mean.pred > cutoff))
  out <- list(predcorr = prederr, mean.pred = mean.pred, optimalk = optimalk,
              cutoff = cutoff, method = clusterings[[1]]$clustermethod,
              Gmax = Gmax, M = M)
  class(out) <- "predstr"
  out
}

getClusPredStrengthResult <- function(distObj,warn=T,defaultMaxNumClusters=10,
                                      psCut=0.8,minClusterSize = 5,filePrefix = "",
                                      usePackagePredStrength = FALSE){
  Gmax <- getMaxNumClustersToTry(defaultMaxNumClusters = defaultMaxNumClusters,
                                 distObj = distObj,filePrefix = filePrefix,
                                 minClusterSize = minClusterSize, warn = warn)
  if(Gmax <= 1){
    res <- NULL
  }else{
    #resC <- predStrengthCustom(distDistinct,species=filePrefix,cutoff = psCut)
    # custom version ^ was used originally because the method didn't support
    # a dist matrix as input, now it does. The new method gives different results, especially
    # when the dataset is small. This isn't due to the difference in clustering method implementation
    # (cluster::pam -> fpc::claraCBI).
    # the results look very wrong from the package version when there is a small number of clusters
    # e.g. 10 clusters for 22 samples
    # this is because a cluster of size 1 has a good score in the new implementation and a bad score in the local
    # implementation
    # the implementation in the original paper is better represented by the new package version
    # but the package implementation is also unstable, even when using 100% of the input
    # but perhaps this reflects real uncertainty
    
    if(usePackagePredStrength){
      res <- fpc::prediction.strength(xdata = distObj,distances = T,cutoff = psCut,
                                      Gmax = Gmax,clustermethod = claraCBI, #pamLike=T, error
                                      Gmin = 2, M = 50,classification = "centroid")  # defaults
    }else{
      res <- predStrengthCustom(distObj,Gmax = Gmax,Gmin = 2,cutoff = psCut,
                                clustermethod = cluster::pam)
    }
    
  }
  return(res)
}

