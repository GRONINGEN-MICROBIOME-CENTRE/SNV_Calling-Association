Find_pheasable_samples = function(AF, Samples, Matrix, homozygous_threshold = 0.9, percentage_homozygous = 0.8){
  #Filter 1. Take quasy-pheasable samples
  #Quasy-pheasable: At least {percentage_homozygous} (e.g 80%) of SNP positions are not heterozygous (MAF<(1-homozygous_threshold), (eg. 0.9) ) IF homozygous_threshold == 1, equivalent than counting number of alleles 1/2
  #We allow up to (1-homozygous_threshold) of  alternative alleles per site.
  Up_threshold = homozygous_threshold
  low_threshold = 1 - homozygous_threshold
  #We are going to count allele frequencies of reference allele, so this can be the major AF or the minor allele frequency. If it is minor it should be lower than low_threshold, if it si major it should be higher than up_threshold
  
  #AF is a datataframe of unname columns. Column 1-2 are species ID and position. Columns 2-last are percentage of reads supporting reference allele in each position.
  #Get out columns of only NAs or with non-numeric variables
  AF[1:dim(Samples)[1],] -> AF_clean
  apply( AF_clean, 1, function(x){  x = x[!is.na(x)] ; length(x[x > Up_threshold | x < low_threshold ])/length(x) } ) -> Heterogeneity_90 #Percentage of good locations
  Samples %>% mutate(Homozygous_perc = Heterogeneity_90) -> Samples
  colnames(Samples) = c("ID","Cohort","Homozygous_perc" )
  
  cbind(Matrix, Samples) %>% filter(Homozygous_perc > homozygous_threshold)  %>% as_tibble() -> Matrix_clustering
  cbind(AF_clean, Samples) %>% filter(Homozygous_perc > homozygous_threshold)  %>% as_tibble() -> AF_clustering
  #Samples$Homozygous_perc > homozygous_threshold -> KEEP
  #KEEP[is.na(KEEP)] = F
  #sum(KEEP == T) %>% print()
  #AF_clean[KEEP, ]  -> AF_clustering
  
  Matrix_clustering %>% select(c("ID","Cohort","Homozygous_perc" )) -> Samples_clustering
  Matrix_clustering %>% select(-c("ID","Cohort","Homozygous_perc" )) %>% as_tibble() -> Matrix_clustering
  AF_clustering %>% select(-c("ID","Cohort","Homozygous_perc" )) %>% as_tibble() -> AF_clustering
  
  return(list(Matrix_clustering, Samples_clustering, AF_clustering))
}

Remove_SNPLocations_with_high_missingRate = function(Matrix, SNP, Porcentage_missing=0.6){
  ###Filter 2. Remove snps with over {percentage_missing} (eg. 0.6) of NA (samples with over 50% NA might not have any common SNP)
  apply(Matrix, 2, function(x){ length(x[is.na(x)])/length(x) } ) -> SNP_missingRate_cluster
  Matrix[,SNP_missingRate_cluster<Porcentage_missing] -> Matrix_clustering
  SNP[SNP_missingRate_cluster<Porcentage_missing,] -> SNP_filterd
  return( list(Matrix_clustering, SNP_filterd ))
}  

Remove_Samples_lowCoverage = function(Matrix, Samples, AF, Porcentage_missing = 0.6){
  ###Filter 3. Remove samples with high missing rate
  ##Might be that the species is not really covered in that sample due to low abundance or absence of the species
  apply(Matrix, 1, function(x){ length(x[is.na(x)])/length(x) } ) -> sample_missingRate_cluster
  Matrix[!sample_missingRate_cluster>Porcentage_missing,] -> Matrix_clustering
  Samples[!sample_missingRate_cluster>Porcentage_missing,] -> Samples_clustering
  AF[!sample_missingRate_cluster>Porcentage_missing,] -> AF_clustering
  return(list(Matrix_clustering, Samples_clustering, AF_clustering))
}

Remove_rare_SNP = function(Matrix, SNP, MAF_filter = 0.1){
  ###Filter 4. Remove SNPs with really low variability (MAF < MAF_filiter), (e.g 0.1)
  Matrix %>% apply(2, function(x){ x = x[!is.na(x)]
                  AFR = ( 2*sum(x==0) +  sum(x==1) ) /(2*length(x))
                  if(AFR > 0.5){ return(1-AFR) }else{ return(AFR) } } ) -> MAF_n
  Matrix[,MAF_n>MAF_filter] -> Matrix_clustering
  SNP[MAF_n>MAF_filter,] -> SNP_clustering
  return(list(Matrix_clustering, SNP_clustering))
}  

Perform_clustering = function(Distance, PC,  minClusterSize = 5, psCut=0.8, threads=1){
  #Functions used here are taken from the metaSNV2 package: https://github.com/metasnv-tool/metaSNV/blob/599ce46260e203100e7439999383b621de81d170/src/subpopr/R/clustering.R
  #Requires fpc package
  Distance = as.dist(Distance)
  print("Getting cluster strength")
  Time_clusterNumber <- Sys.time()
  getClusPredStrengthResult(Distance,minClusterSize = minClusterSize, psCut = psCut) -> Pred_results
  Time_clusterNumberEnd <- Sys.time()  
  print( paste0("Ideal cluster number completed: ", abs(Time_clusterNumber-Time_clusterNumberEnd)))

  numClusters <- Pred_results[["optimalk"]]
  
  #Clustering using the best cluster number
  print("Starting PAM")
  Time_clusterPAM <- Sys.time()
  clustering <- cluster::pam(Distance, numClusters, diss=TRUE)
  Time_clusterPAM_end <- Sys.time()
  print( paste0("PAM completed: ", abs(Time_clusterPAM-Time_clusterPAM_end)))
  K = clustering$clustering
  clustering$numberOfSamplesUsedForClusterDetection <- nrow(as.matrix(Distance))
  clustering$dist <- Distance
  clustering$psVals <- Pred_results[["mean.pred"]]
  clustering$psValsAll <- Pred_results[["predcorr"]]
  clustering$notes <- c()
  clustering$numClusters <- numClusters

  PC  %>% mutate(C = as.factor(K)) %>% ggplot(aes(x=V1, y=V2, col= C)) + geom_point() -> pcoa_plot
  nSamples = length(labels(Distance))
  print("Starting cluster consistency")
  Time_clusterConsistency <- Sys.time()
  #Get cluster stability
  if(nSamples >= 10){ # need at least 10 samples to do this
    #New version of scripts v2, I have used the foreach library to multithread the operations
    minSamplesToUse <- 10
    lowProp <- max(0.3, ceiling(10/nSamples*10)/10)
    # assess clustering stability
    subsampleProportions<-as.list(seq(from=lowProp,to=1,by=0.1))
    clusNumStabilityIter <- 10 # increases time a lot so 10 is enough
    # assess clustering stability for the number of clusters
    clusNumStability <- getClusNumStability2(subsampleProportions = subsampleProportions, nIterClusStability = clusNumStabilityIter, distObj=Distance, psCut = psCut, threads=threads)
    clusNumStabilityPlots <- getClusNumStabilityPlots(clusNumStability)
    # and assess clustering stability for the cluster membership
    clusMembStability <- getClusMembStability2(subsampleProportions = subsampleProportions, numClusters = numClusters, distObj=Distance, threads=threads)
    clusMembStabilityPlots <- getClusMembStabPlots(clusMembStability)
    clusteringStabilityAssessment <- summariseClusteringStability(nClusStability = clusNumStability, clusMembStability = clusMembStability, numClusters = numClusters)
  } else { clusNumStability = NA ; clusNumStabilityPlots = NA ; clusMembStability = NA ; clusteringStabilityAssessment = NA ; clusMembStabilityPlots =NA  }
  Time_clusterConsistencyEnd <- Sys.time()
  print( paste0("Consistency completed: ", abs(Time_clusterConsistency-Time_clusterConsistencyEnd)))

  return(list(K, clusNumStability, clusNumStabilityPlots, clusMembStability, clusMembStabilityPlots, clusteringStabilityAssessment, pcoa_plot  ) )  
  
}

Prevalence = function(X){
  1 - ( sum(is.na(X)) / length(X) )
}

Define_genotyping_snp_MAF = function(AF, Samples_clustering){
  #Gets a series of SNPs that are significantly different (P_bonferroni < 0.05) between clusters
  #Computes the mean AF of alternative alleles in each cluster
  #1. adjusting dimensions of AF and adding IDs and cluster information
  AF[1:(dim(Samples)[1]),] -> AF_t
  AF_t %>% mutate(ID = Samples$ID) -> AF_t
  left_join(Samples_clustering, AF_t, by="ID") -> AF_clusters
  #Linear model associating AF and cluster
  apply(select(AF_clusters, -colnames(Samples_clustering)) , 2 , function(x){
    y = x[!is.na(x)]
    if (length(y) == 0 ){ return(NA) }
    P = (sum(y==0)/length(y))
    if(P > 0.9 | P <0.1){ return(NA) }
    C = AF_clusters$Cluster[!is.na(x)]
    if (length(unique(C)) < 2){ return(NA) }
    lm(y ~ 1) -> model_0
    lm(y ~ as.factor(C)) -> model_1 
    anova(model_0, model_1) -> M
    return(M$`Pr(>F)`[2])
  }) -> Abundance_association
  tibble(SNP = names(Abundance_association), P = Abundance_association) -> Abundance_association_df
  #Removing P_bonferroni<0.05 SNPs
  Abundance_association_df %>% mutate(P_bonf = p.adjust(P,"bonferroni") ) %>% arrange(P) -> Abundance_association_df
  Abundance_association_df %>% filter(P_bonf < 0.05) -> SNPs_investigate
  #Getting Mean of each SNP per cluster and adding to the SNPs_investigate data frame
  AF_clusters %>% select(c("Cluster", as.character(SNPs_investigate$SNP))) %>% group_by(Cluster) %>% summarise_each( funs(mean(., na.rm = TRUE) )) -> Mean_per_cluster
  AF_clusters %>% select(c("Cluster", as.character(SNPs_investigate$SNP))) %>% group_by(Cluster) %>% summarise_each( funs(Prevalence(.) )) -> Prevalence_per_cluster
  t(Mean_per_cluster) %>% as.data.frame() %>% rownames_to_column("SNP") %>% as_tibble() -> Mean_per_cluster
  t(Prevalence_per_cluster) %>% as.data.frame() %>% rownames_to_column("SNP") %>% as_tibble() -> Prevalence_per_cluster
  colnames(Mean_per_cluster) = Mean_per_cluster[1,] ; colnames(Prevalence_per_cluster) = Prevalence_per_cluster[1,]
  left_join(Prevalence_per_cluster,Mean_per_cluster, by="Cluster", suffix=c("_prevalence","")) -> Mean_n_prev_per_cluster
  Mean_n_prev_per_cluster %>% filter(!Cluster == "Cluster") %>% mutate(SNP = Cluster) %>% select(-Cluster) -> Mean_n_prev_per_cluster 
  left_join(SNPs_investigate, Mean_n_prev_per_cluster) -> SNPs_investigate
  #Get the difference between the first and second most abundant clusters. E.g if 3 clusters and abundances are 0.2 0.3 and 0.5, the difference will be 0.5-0.3= 0.2, it also gives which is the most abundant cluster
  All_dif = tibble()
  for (SNP_position in  seq(dim(SNPs_investigate)[1]) ){
    SNP_position = SNPs_investigate[SNP_position,]
    SNP_position %>% select(-SNP) %>% select(colnames(select(Prevalence_per_cluster, -Cluster))) %>% as_vector() -> x
    x= arrange(tibble(Value = x, Position = seq(length(x))),desc(Value))
    D=head(x,n=1)$Value - head(x,n=2)$Value[2] 
    tibble(SNP = SNP_position$SNP, Minimum_abundance_dif = D, Biggest_abundance = head(x,n=1)$Position) -> Diff
    rbind(All_dif, Diff) -> All_dif
  }
  left_join(SNPs_investigate, All_dif) -> Genotype_SNP
  return(Genotype_SNP)
}

Genotype_samples = function(g_SNP, M, Samples, Clusters,Min_abundace_dif=0.8, Prevalence_threshold=0 ){
  #Filter SNPs where the abundance differences is not of at least Min_abundance_dif
  g_SNP %>% filter(Minimum_abundance_dif > Min_abundace_dif)  -> g_SNP
  if (dim(g_SNP)[1] < 10){
	print("Warning: Few genotyping SNPs were identified for the given minimal abundance difference. Samples were not genotyped")
	return(NA)
  }
  #Get prevalence of dominant cluster, filter
  Clusters = seq( 1:Clusters )
  Prev_columns = paste( Clusters, "_prevalence" , sep ="" )
  Prevalence = c()
  #Prev_columns =  paste( levels(as.factor(g_SNP$Biggest_abundance)) , "_prevalence" , sep ="" ) ; Prevalence = c()
  for (i in  rownames(g_SNP) ){ Position = g_SNP[i,] ; Column_s = Prev_columns[Position$Biggest_abundance] ; Prevalence = c(Prevalence, as_vector(select(Position, Column_s)) ) }
  g_SNP %>% mutate(Prevalence_in_major_cluster = Prevalence) %>% filter(Prevalence_in_major_cluster > Prevalence_threshold) -> g_SNP
  #Clusters = unique(g_SNP$Biggest_abundance)
  Genotyped = tibble(ID = Samples$ID)
  M[1:dim(Samples)[1],] -> M
  #Compute the mean number of Genotype alleles covered per sample and per cluster
  for (Cluster_n in sort(Clusters)){
    g_SNP %>% filter(Biggest_abundance == Cluster_n) -> g_SNP_c
    if (dim(g_SNP_c)[1] == 0){ Percentage_SNP = NA 
    } else {
    M %>% select(as.character(g_SNP_c$SNP)) %>% apply(1, function(x){ x=x[!is.na(x)] ; median(x)  } ) -> Percentage_SNP
    }
    Genotyped %>% mutate( P = Percentage_SNP ) -> Genotyped
    colnames(Genotyped)[Cluster_n+1] = as.character(Cluster_n)
  }
  return(Genotyped)
  
}

