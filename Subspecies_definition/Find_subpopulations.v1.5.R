library(tidyverse)
library(vegan)
library(fpc)
library(cluster)
#For multithreads
library(foreach)
library(doParallel)
#New
library(data.table)

source(file = "/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dwang/SNP/subspecies/bin/Functions_metasnv.R")
source(file ="/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dwang/SNP/subspecies/bin/Functions.v1.5.R")
source(file = "/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dwang/SNP/subspecies/bin/scaledDist.R")

set.seed(46733)

#ml RPlus/4.0.3-foss-2018b-v21.12.10

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
Multi_th = as.numeric(args[6]) #Number of threats to parallelize
Number_SNP = as.numeric(args[7]) #Number of SNPs to use
Number_samples = as.numeric(args[8]) #Number of samples to use
setwd(Output_dir)


Welcome_message =  paste0("Arguments Provided:
Species code:", Species, "
Output Directory:", Output_dir,"
Prefix genetics : ", args[3] ,"
Allele Frequency table: ", args[4],"
Cohort information: ", args[5],"
Threads to use: ", args[6],"
SNPs to use: ", args[7]," 
Samples to use: ", args[8])
writeLines(Welcome_message)


print( paste( c("Running analysis in ", Species), collapse="" ) )
if(!dir.exists(Output_dir)){dir.create(Output_dir)}

Read_data  = function(){
  Samples    = paste(args[3], ".012.indv", sep="")
  Matrix     = paste(args[3], ".012",      sep="")
  SNP        = paste(args[3], ".012.pos",  sep="")
  AF         = args[4]
  CohortInfo = args[5]
  
  ##Reading with readR. No column names on bcftools-generated files. NAs are encoded as either "-1" or "."
  fread(Samples, header = F)                               -> Samples
  fread(Matrix, header = F, na.strings = "-1", sep = "\t", drop="V1") -> Matrix
  fread(SNP,  header = F, sep = "\t")                      -> SNP
  fread(AF, header = F, na.strings = ".", sep ="\t", drop=c("V1","V2"))       -> AF
  fread(CohortInfo, header = T)                            -> CohortInfo
  
  snp_location = paste(SNP$V1, SNP$V2, sep=":")
  SNP[,X2 := snp_location]
  
  #Some formatting
  transpose(AF) -> AF #Get Participants in rows
  setnames(Samples, old = "V1", new = "ID")
  colnames(Samples)[1] = "ID"
  #Get SNPs as column names
  colnames(Matrix) = SNP$X2
  colnames(AF) = SNP$X2
  #AF and genotype matrices should be of the same dimensions samples x SNP
  
  #Add cohort information to Samples metadata
  Samples[, Cohort := CohortInfo$Cohort[match(Samples$ID,CohortInfo$ID)] ]
  
  return(list(Samples, AF, Matrix, SNP))
}

Sys.time() -> Time_read1
Read_data() -> AllData
Sys.time() -> Time_read2
Samples = AllData[[1]]
AF = AllData[[2]]
Matrix = AllData[[3]]
SNP = AllData[[4]]

print( paste0( "SNPs: ", dim(SNP) ) )
print( paste0( "Matrix: ", dim(Matrix) ) )
print( paste0( "AF matrix: ", dim(AF) ) )
print( paste0( "Communities: ", dim(AF) ) )


Time_elapsed  = rbind(Time_elapsed, tibble(Operation = "Read", Samples = dim(Samples)[1], SNPs = dim(SNP)[1], Time = Time_read2-Time_read1, C=1  ) )


Run_samplePickUp_and_cluster =  function(AF, Samples, Matrix, Number_samples = 0, Number_SNP = 0, Iterations_consistency=10, Save=T){
  print("Selection of samples and SNPs to be used in clustering analysis")
  Time_SNP_sample <- Sys.time()
  print("Sample selection")
  #Samples pick up
  Find_pheasable_samples(AF, Samples, Matrix, homozygous_threshold = 0.9, percentage_homozygous = 0.8) -> R1
  Matrix_clustering = R1[[1]] ; Samples_clustering = R1[[2]] ; AF_clustering = R1[[3]]
  #I added this step before the third since if many samples have low coverage then the proportion of SNPs to remove because large NA is larger
  Remove_Samples_lowCoverage(Matrix_clustering, Samples_clustering, AF_clustering, Porcentage_missing =  0.6) -> R3
  Matrix_clustering = R3[[1]] ; Samples_clustering = R3[[2]] ; AF_clustering = R3[[3]]
  
  print("SNP selection")
  #SNPs pick up
  Remove_SNPLocations_with_high_missingRate(Matrix_clustering,SNP, Porcentage_missing=0.6) -> R2
  Matrix_clustering = R2[[1]] ; SNP_filtered = R2[[2]]
  
  Remove_rare_SNP(Matrix_clustering, SNP_filtered, MAF_filter = 0.1) -> R4
  Matrix_clustering = R4[[1]] ; SNP_filtered = R4[[2]]
  
  SNP_names = as.character(SNP_filtered$X2)
  AF_clustering[,..SNP_names] -> AF_clustering
  Time_SNP_sample_end <- Sys.time()
  print( paste0("Selection completed, time elapsed:", abs(Time_SNP_sample-Time_SNP_sample_end) ))
  Time_elapsed  = rbind(Time_elapsed, tibble(Operation = "Selection", Samples = dim(Matrix_clustering)[1], SNPs = dim(Matrix_clustering)[2], Time = Time_SNP_sample_end-Time_SNP_sample, C=1  ) )
  print("Sample and SNP subsampling")
  #Pick ups in number of samples and SNPs to use for clustering
  if (! Number_samples == 0){	
    #Keep just {Number_samples} samples
    if (dim(Samples_clustering)[1] > Number_samples ){
      print(paste0("Subsampling number of Samples to ", as.character(Number_samples)))
      sample(Samples_clustering$ID, Number_samples ,replace = F) -> Samples_choice
      Samples_row = Samples_clustering$ID %in% Samples_choice

      #Samples_clustering -> Samples_clustering_nosubs
      
      Samples_clustering[ID %chin% Samples_choice] -> Samples_clustering
      Matrix_clustering[Samples_row, ] -> Matrix_clustering
      AF_clustering[Samples_row, ] -> AF_clustering
    }
  }
  if (! Number_SNP == 0){
    #Keep a maximum of {Number_SNP} SNPs
    if (dim(Matrix_clustering)[2] > Number_SNP ){
      print(paste0("Subsampling number of SNPs to ", as.character(Number_SNP)))
      sample(colnames(Matrix_clustering) , Number_SNP ,replace = F) -> SNP_choice
      #Matrix_clustering -> Matrix_clustering_subset
      
      Keep_SNP = colnames(Matrix_clustering)[colnames(Matrix_clustering) %in% SNP_choice]
      Matrix_clustering[ , ..Keep_SNP ] -> Matrix_clustering
      AF_clustering[ , ..Keep_SNP  ] -> AF_clustering
      SNP_filtered[X2 %chin% SNP_choice] -> SNP_filtered
    }
  }
  
  
  print("Check that all dimensions match")
  dim(Matrix_clustering) %>% print()
  dim(Samples_clustering) %>% print()
  dim(SNP_filtered) %>% print()
  dim(AF_clustering) %>% print()
  
  print("Compute Distance")
  Time_SNP_dist <- Sys.time()
  #scaledManhattan(Matrix_clustering) -> Distance
  vegdist(Matrix_clustering, na.rm = T, method = "manhattan") -> Distance
  as.matrix(Distance) -> Distance
  rownames(Distance) = seq(1, dim(Matrix_clustering)[1] )
  K =  min(dim(Matrix_clustering)[1] - 1, dim(SNP_filtered)[1]-1 )
  #PCOA = ape::pcoa(Distance)
  PCOA = stats::cmdscale(Distance, eig = T, k = K)
  #as_tibble(PCOA$vectors) %>% mutate(Cohort = Samples_clustering$Cohort)  %>% ggplot(aes(x=Axis.1, y=Axis.2, color=Cohort)) + geom_point() -> PCA_color_by_cohort
  as_tibble(PCOA$points) %>% mutate(Cohort = Samples_clustering$Cohort)  %>% ggplot(aes(x=V1, y=V2, color=Cohort)) + geom_point() -> PCA_color_by_cohort
  as_tibble(PCOA$points)[,1:2] -> PCs
  colnames(PCs) = c("V1", "V2")
  PCs = tibble(PCs)
  Time_SNP_dist_end <- Sys.time()
  print( paste0("Distance completed, time elapsed:", abs(Time_SNP_dist-Time_SNP_dist_end) ))
  Time_elapsed  = rbind(Time_elapsed, tibble(Operation = "Distance", Samples = dim(Matrix_clustering)[1], SNPs = dim(Matrix_clustering)[2], Time = Time_SNP_dist_end-Time_SNP_dist, C = 1  ) )
  
  
  print("Clustering - Slow")
  Time_SNP_clustering <- Sys.time()
  #This uses scripts from metaSNV2 for picking optimal number of clusters and assessing cluster and sample-in-cluster stabilities. I added multithreadening option
  Clustering_results = Perform_clustering(Distance, PCs, threads = Multi_th, Iterations=Iterations_consistency)
  Samples_clustering[, Cluster := Clustering_results[[1]]  ] -> Samples_clustering
  unlist(Clustering_results[[6]]) %>% as.data.frame() %>% rownames_to_column("Info") -> Report
  Time_SNP_clustering_end <- Sys.time()
  print( paste0("Clustering completed, time elapsed:", abs(Time_SNP_clustering-Time_SNP_clustering_end) ))
  Time_elapsed  = rbind(Time_elapsed, tibble(Operation = "Clustering", Samples = dim(Matrix_clustering)[1], SNPs = dim(Matrix_clustering)[2], Time = Time_SNP_clustering_end-Time_SNP_clustering,C = Multi_th  ) )
  
  
  
  #I have also tested using AF_clustering as input for distance and clustering. Results are similar.
  print("Clustering done, number of found clusters:")
  # print(unique(Clustering_results[[6]]))
  print(length(unique(as.factor(Samples_clustering$Cluster))))
  if (Save == T){
    
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
  }
  
  return( list(PCOA, Clustering_results, Matrix_clustering,Samples_clustering, AF_clustering,SNP_filtered, Time_elapsed))
}

Run_samplePickUp_and_cluster(AF, Samples, Matrix, Number_samples, Number_SNP, Save=T, Iterations_consistency = 5) -> Output
Output[[1]] -> PCOA
Output[[2]] -> Clustering_results
Output[[3]] -> Matrix_clustering ; Output[[4]] -> Samples_clustering; Output[[5]] -> AF_clustering; Output[[6]] -> SNP_filtered
Output[[7]] -> Time_elapsed
K = length(unique(as.factor(Samples_clustering$Cluster)))

if ( K < 2){
  print("Finishing script since only 1 cluster is supported")
  paste(Species, ".Time.tsv", sep = "") -> Output_time #Association_PC
  write_tsv(Time_elapsed, Output_time)
  q()
}


#########################################################
#Define Genotype SNP and count how many in each sample###
#########################################################

#Similar strategy than in metaSNV2. We define SNPs that are enriched in specific clusters. Then we count how many alleles of that SNP we have present in each sample
print("Identifying genotyping SNPs")
Time_SNP_genotyping <- Sys.time()

write_rds(AF, paste(Species,".AFs.rds", sep = ""))
#write_rds(Samples_clustering, paste(Species, ".Samples_clustering.rds",sep = ""))
Genotype_SNP = Define_genotyping_snp_MAF(AF, Samples_clustering)
print("Quantifying genotyping SNPs in all samples")
#do median of Allele frequency
Genotyped = Genotype_samples(Genotype_SNP, AF, Samples, Min_abundace_dif = 0.8, Prevalence_threshold = 0.3, Clusters=K ) #Change prevalence threshold accordingly.

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
#PCOA$vectors -> PCs
PCOA$points %>% as_tibble() -> PCs
#PCOA$values %>% as_tibble() %>% filter(Cum_corr_eig<0.2) %>% dim() -> N_com#Number of componetns to retrieve 20% of SNP variability
PCOA$eig %>% as_tibble() %>% filter(! value < 0 ) %>% mutate(Variability = value/sum(value), CumVariability = cumsum(Variability) ) %>% print() %>% filter(CumVariability<0.5) %>% dim() -> N_com
Inv_rank = function(x){ qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}
Association_PC = tibble()
for (i in seq(N_com[1])){
  PCs[,i] -> PC
  Samples_clustering %>% mutate(PC = Inv_rank(as_vector(PC))) -> Samples_model
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
