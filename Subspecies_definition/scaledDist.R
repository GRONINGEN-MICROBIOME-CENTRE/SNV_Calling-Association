require(tidyverse)
require(parallelDist)
require(RcppArmadillo)
require(RcppXPtrUtils)

# Enocode missing values as -9.0 (any value except for 0,1,2)
scaledDistFuncPtr <- cppXPtr( 
  "double customDist(const arma::mat &A, const arma::mat &B) {
    auto C = (A!=-9) + (B!=-9); // Sum of shared non-missing element will be 2
    arma::uvec index = arma::find(C > 1.0); // Get indices of non-missing element shared by A and B
    if(index.n_elem>0) { // 0 can be changed to other cutoffs
      // return sqrt(arma::accu(arma::square(A.elem(index) - B.elem(index)))/(4*index.n_elem)); // scaled euclidean distance
      return arma::accu(arma::abs(A.elem(index) - B.elem(index)))/(2*index.n_elem); // scaled manhattan distance 
    } else {
      return -9.0;
    }
  }",
  depends = c("RcppArmadillo"))

sharedNumFuncPtr <- cppXPtr(
 "double customDist(const arma::mat &A, const arma::mat &B) {
  auto C = (A!=-9) + (B!=-9);
  arma::uvec index = arma::find(C > 1.0);
  return index.n_elem;
 }",
 depends = c("RcppArmadillo"))

## Remove samples with NA in distance matrix
dist_rmna<-function(inDist){
  while(sort(colSums(is.na(inDist)),decreasing = T)[1] > 0){
    rmid<-names(sort(colSums(is.na(inDist)),decreasing = T)[1])
    inDist<-inDist[-match(rmid, rownames(inDist)),-match(rmid, colnames(inDist))]
  }
  return(inDist)
}


scaledManhattan<-function(inMat){
	inMat = as.matrix(inMat)
	inMat[is.na(inMat)] <- -9
	distMat<- parallelDist(inMat, method="custom", func = scaledDistFuncPtr) %>% as.matrix
	rownames(distMat)<-rownames(inMat);colnames(distMat)<-rownames(inMat)
	distMat[distMat==-9.0]<-NA
	return(distMat)
}

sharedFeatureN<-function(inMat){
	inMat[is.na(inMat)] <- -9
	nMat<- parallelDist(inMat, method="custom", func = sharedNumFuncPtr) %>% as.matrix
	rownames(nMat)<-rownames(inMat);colnames(nMat)<-rownames(inMat)
	return(nMat)
}

# how to use:
# (1) Calculate scaled manhattan distance
# distMatrix <- scaledManhattan(yourMatrix) # Get scaled manhattan distance matrix, please note that the missing values in input matrix should be NA
# distMatrix.rmna <- dist_rmna(distMatrix) # Iteratively remove columns and rows with NA, and also try to keep more samples.

# (2) Calculate shared feature numbers (Ns) between all sample pairs, get the N matrix
# nMatrix <- sharedFeatureN(yourMatrix) 
