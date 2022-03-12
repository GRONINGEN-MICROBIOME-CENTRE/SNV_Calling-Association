library(tidyverse)
library(vegan)
library(fpc)
library(cluster)
#For multithreads
library(foreach)
library(doParallel)

source(file = "/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dwang/SNP/subspecies/bin/Functions_metasnv.R")
Time_elapsed = tibble() # Measures the time taking by each operation

####Load files
#############################
#Prior to this script bcftools should be used to
#1- extract chromosomes belonging to a species
#2- Create an 012 genotype encoding
#3- Extract allele frequencies of alternative allele (so that it matches with 012)
############################
#Files should be located in Prefix directory
#Species indicates the name of the file

args = commandArgs(trailingOnly=TRUE)

Species    = args[1] #Species = "102454" #Akkermansia Municipilla
Output_dir = args[2]
Multi_th = args[6] #Number of threats to parallelize


setwd(Output_dir)

print( paste( c("Running analysis in ", Species), collapse="" ) )
if(!dir.exists(Output_dir)){dir.create(Output_dir)}

Read_data  = function(){
  Samples    = paste(args[3], ".012.indv", sep="")
  Matrix     = paste(args[3], ".012",      sep="")
  SNP        = paste(args[3], ".012.pos",  sep="")
  AF         = args[4]
  CohortInfo = args[5]

  ##Reading with readR. No column names on bcftools-generated files. NAs are encoded as either "-1" or "."
  read_tsv(Samples,    col_names = F)            -> Samples
  read_tsv(Matrix,     col_names = F)            -> Matrix
  read_tsv(SNP,        col_names = F, na = "-1") -> SNP
  read_tsv(AF,         col_names=F, na = ".")    -> AF
  read_tsv(CohortInfo, col_names=T)              -> CohortInfo

  snp_location = paste(SNP$X1, SNP$X2, sep=":")
  SNP$X2 = snp_location

  #Some formatting
  Matrix %>% select(-X1) -> Matrix
  colnames(Matrix) = SNP$X2

  AF %>% select(-c(X1, X2)) -> AF
  t(AF) %>% as_tibble() -> AF
  colnames(AF) = SNP$X2
  
  colnames(Samples)[1] = "ID"
  #AF and genotype matrices should be of the same dimensions samples x SNP

  #Add cohort information to Samples metadata
  Samples$Cohort<-CohortInfo$Cohort[match(Samples$ID,CohortInfo$ID)]
#  Samples %>% 
#	mutate(Cohort = ifelse( grepl("LLD", ID ), "LLD",
#	 ifelse(grepl("HV", ID ), "500FG",
#	 ifelse(grepl("IBD", ID ), "IBD", 
#	 ifelse(grepl("C", ID ), "300OB", NA) )) ) ) -> Samples
 
 return(list(Samples, AF, Matrix, SNP))
}

Sys.time() -> Time_read1
Read_data() -> AllData
Sys.time() -> Time_read2
Samples = AllData[[1]]
AF = AllData[[2]]
Matrix = AllData[[3]]
SNP = AllData[[4]]

Time_elapsed  = rbind(Time_elapsed, tibble(Operation = "Read", Samples = dim(Samples)[1], SNPs = dim(SNP)[1], Time = Time_read2-Time_read1, C=1  ) )
###############################################################
##########################FUNCTIONS############################
###############################################################

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

Perform_clustering = function(Distance, minClusterSize = 5, psCut=0.8, threads=1){
  #Functions used here are taken from the metaSNV2 package: https://github.com/metasnv-tool/metaSNV/blob/599ce46260e203100e7439999383b621de81d170/src/subpopr/R/clustering.R
  #Requires fpc package
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
  as_tibble(ape::pcoa(Distance)$vectors) %>% mutate(C = as.factor(K)) %>% ggplot(aes(x=Axis.1, y=Axis.2, col= C)) + geom_point() -> pcoa_plot
  
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
  }
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

Genotype_samples = function(g_SNP, M, Samples, Min_abundace_dif=0.8, Prevalence_threshold=0 ){
  #Filter SNPs where the abundance differences is not of at least Min_abundance_dif
  g_SNP %>% filter(Minimum_abundance_dif > Min_abundace_dif)  -> g_SNP
  #Get prevalence of dominant cluster, filter
  Prev_columns =  paste( levels(as.factor(g_SNP$Biggest_abundance)) , "_prevalence" , sep ="" ) ; Prevalence = c()
  for (i in  rownames(g_SNP) ){ Position = g_SNP[i,] ; Column_s = Prev_columns[Position$Biggest_abundance] ; Prevalence = c(Prevalence, as_vector(select(Position, Column_s)) ) }
  g_SNP %>% mutate(Prevalence_in_major_cluster = Prevalence) %>% filter(Prevalence_in_major_cluster > Prevalence_threshold) -> g_SNP
  
  Clusters = unique(g_SNP$Biggest_abundance)
  Genotyped = tibble(ID = Samples$ID)
  M[1:dim(Samples)[1],] -> M
  #Compute the mean number of Genotype alleles covered per sample and per cluster
  for (Cluster_n in sort(Clusters)){
    g_SNP %>% filter(Biggest_abundance == Cluster_n) -> g_SNP_c
    M %>% select(as.character(g_SNP_c$SNP)) %>% apply(1, function(x){ x=x[!is.na(x)] ; median(x)  } ) -> Percentage_SNP
    Genotyped %>% mutate( P = Percentage_SNP ) -> Genotyped
    colnames(Genotyped)[Cluster_n+1] = as.character(Cluster_n)
  }
  return(Genotyped)
  
}

print("Selection of samples and SNPs to be used in clustering analysis")
Time_SNP_sample <- Sys.time()

Find_pheasable_samples(AF, Samples, Matrix, homozygous_threshold = 0.9, percentage_homozygous = 0.8) -> R1
Matrix_clustering = R1[[1]] ; Samples_clustering = R1[[2]] ; AF_clustering = R1[[3]]
#I added this step before the third since if many samples have low coverage then the proportion of SNPs to remove because large NA is larger
Remove_Samples_lowCoverage(Matrix_clustering, Samples_clustering, AF_clustering, Porcentage_missing =  0.6) -> R3
Matrix_clustering = R3[[1]] ; Samples_clustering = R3[[2]] ; AF_clustering = R3[[3]]

Remove_SNPLocations_with_high_missingRate(Matrix_clustering,SNP, Porcentage_missing=0.6) -> R2
Matrix_clustering = R2[[1]] ; SNP_filtered = R2[[2]]

Remove_rare_SNP(Matrix_clustering, SNP_filtered, MAF_filter = 0.1) -> R4
Matrix_clustering = R4[[1]] ; SNP_filtered = R4[[2]]

AF_clustering %>% select(as.character(SNP_filtered$X2)) -> AF_clustering
Time_SNP_sample_end <- Sys.time()
print( paste0("Selection completed, time elapsed:", abs(Time_SNP_sample-Time_SNP_sample_end) ))
Time_elapsed  = rbind(Time_elapsed, tibble(Operation = "Selection", Samples = dim(Matrix_clustering)[1], SNPs = dim(Matrix_clustering)[2], Time = Time_SNP_sample_end-Time_SNP_sample, C=1  ) )



print("Check that all dimensions match")
dim(Matrix_clustering) %>% print()
dim(Samples_clustering) %>% print()
dim(SNP_filtered) %>% print()
dim(AF_clustering) %>% print()

print("Compute Distance")
Time_SNP_dist <- Sys.time()
#Distance based on genotype matrix
as.data.frame(Matrix_clustering) -> Matrix_clustering
#Somehow without this line the distance matrix has no labels
rownames(Matrix_clustering) = rownames(Matrix_clustering) 
# na.rm = T makes sure that only complete pairwise cases are used to compute distances
vegdist(Matrix_clustering, method = "manhattan", na.rm=T) -> Distance
PCOA = ape::pcoa(Distance)
as_tibble(PCOA$vectors) %>% mutate(Cohort = Samples_clustering$Cohort)  %>% ggplot(aes(x=Axis.1, y=Axis.2, color=Cohort)) + geom_point() -> PCA_color_by_cohort
Time_SNP_dist_end <- Sys.time()
print( paste0("Distance completed, time elapsed:", abs(Time_SNP_dist-Time_SNP_dist_end) ))
Time_elapsed  = rbind(Time_elapsed, tibble(Operation = "Distance", Samples = dim(Matrix_clustering)[1], SNPs = dim(Matrix_clustering)[2], Time = Time_SNP_dist_end-Time_SNP_dist, C = 1  ) )


print("Clustering - Slow")
Time_SNP_clustering <- Sys.time()
#This uses scripts from metaSNV2 for picking optimal number of clusters and assessing cluster and sample-in-cluster stabilities. I added multithreadening option
Clustering_results = Perform_clustering(Distance, threads = Multi_th)
Samples_clustering %>% mutate(Cluster = Clustering_results[[1]]) -> Samples_clustering
unlist(Clustering_results[[6]]) %>% as.data.frame() %>% rownames_to_column("Info") -> Report
Time_SNP_clustering_end <- Sys.time()
print( paste0("Clustering completed, time elapsed:", abs(Time_SNP_clustering-Time_SNP_clustering_end) ))
Time_elapsed  = rbind(Time_elapsed, tibble(Operation = "Clustering", Samples = dim(Matrix_clustering)[1], SNPs = dim(Matrix_clustering)[2], Time = Time_SNP_clustering_end-Time_SNP_clustering,C = Multi_th  ) )

#I have also tested using AF_clustering as input for distance and clustering. Results are similar.
print("Clustering done, number of found clusters:")
# print(unique(Clustering_results[[6]]))
print(length(unique(as.factor(Samples_clustering$Cluster))))

print("Writing Clustering report")
paste(Species, ".Clusters_discovery.tsv", sep = "")                        -> Output_ClusterDiscovery
paste(Species, ".Clusters_discovery_consistency.tsv", sep = "")            -> Output_ClusterDiscoveryConsistency
paste(Species, ".Clusters_discovery_AsssignmentConsistency.tsv", sep = "") -> Output_AssignmentDiscoveryConsistency
paste(Species, ".Stability_report.tsv", sep = "")                          -> Output_StabilityReport

paste(Species, ".Clusters_discovery_consistency.pdf", sep = "")            -> Output_ClusterDiscoveryConsistencyPlot
paste(Species, ".Clusters_discovery_AsssignmentConsistency.pdf", sep = "") -> Output_AssignmentDiscoveryConsistencyPlot
paste(Species, ".Clusters_discovery_PCoA.pdf", sep = "")                   -> Output_ClusterDiscoveryPlot


write_tsv(Samples_clustering, Output_ClusterDiscovery)
write_tsv(Clustering_results[[2]], Output_ClusterDiscoveryConsistency)
write_tsv(Clustering_results[[4]], Output_AssignmentDiscoveryConsistency)
write_tsv(Report, Output_StabilityReport)

ggsave(plot = Clustering_results[[7]], filename = Output_ClusterDiscoveryPlot)
ggsave(plot = Clustering_results[[3]]$boxHeat, filename = Output_ClusterDiscoveryConsistencyPlot)
ggsave(plot = Clustering_results[[5]][[2]], filename = Output_AssignmentDiscoveryConsistencyPlot)


if (length(unique(as.factor(Samples_clustering$Cluster))) < 2){
	print("Finishing script since only 1 cluster is supported")
	q()
}


#########################################################
#Define Genotype SNP and count how many in each sample###
#########################################################

#Similar strategy than in metaSNV2. We define SNPs that are enriched in specific clusters. Then we count how many alleles of that SNP we have present in each sample
print("Identifying genotyping SNPs")
Time_SNP_genotyping <- Sys.time()

write_rds(AF, paste(Species,".AFs.rds", sep = ""))
write_rds(Samples_clustering, paste(Species, ".Samples_clustering.rds",sep = ""))
Genotype_SNP = Define_genotyping_snp_MAF(AF, Samples_clustering)
print("Quantifying genotyping SNPs in all samples")
#do median of Allele frequency
Genotyped = Genotype_samples(Genotype_SNP, AF, Samples, Min_abundace_dif = 0.8, Prevalence_threshold = 0.8 )

#Check how do the percentages relate with the known clusters. And how many samples that we did not cluster can know be genotyped
left_join(Genotyped, Samples_clustering) %>% select(-c(Homozygous_perc, Cohort)) -> Check_genotype
Check_genotype %>% filter(is.na(Cluster) & ( ! is.na(`1`) | ! is.na(`2`)   )) %>% print(n=200)
Time_SNP_genotyping_end <- Sys.time()
print(paste0 ("Genotyping SNP completed, time: ", Time_SNP_genotyping - Time_SNP_genotyping_end))
Time_elapsed  = rbind(Time_elapsed, tibble(Operation = "Genotyping", Samples = dim(Matrix_clustering)[1], SNPs = dim(Matrix_clustering)[2], Time = Time_SNP_genotyping_end-Time_SNP_genotyping,C = 1  ) )


#########################################################################
##Check how many PCs needed to explain most Cluster variability##########
#########################################################################
print("Checking which PCs are associated with subspecies clustering")
#We will check which PCs represent subspecies differences. We will use those for correcting in GWAS
PCOA$vectors -> PCs
PCOA$values %>% as_tibble() %>% filter(Cum_corr_eig<0.2) %>% dim() -> N_com#Number of componetns to retrieve 20% of SNP variability

Inv_rank = function(x){ qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}
Association_PC = tibble()
for (i in seq(N_com[1])){
  PCs[,i] -> PC
  Samples_clustering %>% mutate(PC = Inv_rank(PC)) -> Samples_model
  lm(PC ~ 1, Samples_model) -> model_0
  lm(PC ~ as.factor(Cluster), Samples_model) -> model_r
  R2 = summary(model_r)$adj.r.squared
  P = anova(model_0, model_r)$`Pr(>F)`[2]
  rbind(Association_PC, tibble(PC = i, P_value = P, R2_of_PC = R2) ) -> Association_PC
}
Association_PC %>% arrange(P_value) -> Association_PC

#########################################
#############Final Report################
#########################################
print("Writing reports")

##2. Genotyping
paste(Species, ".Genotpe_SNP.tsv",       sep = "")       -> Output_gSNV #Genotype_SNP
paste(Species, ".Samples_genotyped.tsv", sep = "") -> Output_Genotyped #Check_genotype

write_tsv(Genotype_SNP , Output_gSNV)
write_tsv(Check_genotype , Output_Genotyped)

##3. Number of components
paste(Species, ".PC_association_Cluster.tsv", sep = "") -> Output_PCAssociation #Association_PC
write_tsv(Association_PC , Output_PCAssociation)
paste(Species, ".Time.tsv", sep = "") -> Output_time #Association_PC
write_tsv(Time_elapsed, Output_time)
