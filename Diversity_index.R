#Calculate the diversity index based on Leinster and Cobbold Ecology (2012)
#Assume a similarity matrix (Z; n x n) for n species (which can include all species)
#Assume a community with species a subset of these n species

#The user needs to construct a similarity matrix Z based on taxonomy or phylogeny or functional traits
#Only one Z is needed for all samples!!
#If the species is absent from the community, its relative abundance = 0

#Exponential Shannon index
ExpS <- function(p){
 #p is the vector of  abundances
  b <- p/sum(p)
  return(exp(-sum(b*log(b))))
}

#Rao's entropy (See also Stirling 2007)
Rao2 <- function(U, D = diag(0, nr = length(U), nc = length(U))){
	#U: an abundance vector  of n species. 
  #D: distance matrix
 	stopifnot(nrow(D) ==   ncol(D))
	stopifnot(nrow(D) == length(U))
	
 	#Each element of U must be positive
	stopifnot(any(U >= 0)) 
	
	nsp <- nrow(D) #Number of species
	
	#Convert to relative abundance
	b <- U/sum(U)
	b1<- rep(b, nsp)
	
	#Construct two matricies
	M1 <- matrix(b1, nr = nsp, nc = nsp, byrow = T)
	M2 <- matrix(b1, nr = nsp, nc = nsp, byrow = F)
	
	return(sum(M1*M2*D))
}

#First argument only takes one vector, suitable for many communities with huge amount of species
qDz_v <- function(U, Z = diag(1, nr = length(U), nc = length(U)), q = 0){
	#U: an abundance vector  of n species. 
  #Z: a n x n matrix of species similarity
  #q: value controlling preferences on rare vs. dominant species. 
  #if Z is an identity matrix and q = 0, the output should be the richness
  #if Z is an identity matrix and q = 1, the output should be exponential Shannon index
  
	stopifnot(nrow(Z) == ncol(Z))
	stopifnot(nrow(Z) == length(U))
    
	#Each element of U must be positive
	stopifnot(any(U >= 0))
	nsp <- nrow(Z) #Number of species
	
	#Convert to relative abundance
	b <- U/sum(U)
  Zp<- Z %*% b  #n x 1 matrix
    
  if (q != 1){
    Y <- ( sum(b*Zp**(q-1)) )**( 1/(1-q) )
  }else{
    cff <- prod(Zp**b)
    Y <- 1/cff
  }
	return(Y)
}

qDz <- function(U, Z, q = 0){
	#U: a n x nc matrix of nc communities, each community a subset of n species. Each column is a community 
  #   the [i,j] element of U represents the relative abundance of the ith species
  #   in the jth community (which can be zero)
  #Z: a n x n matrix of species similarity
  #q: value controlling preferences on rare vs. dominant species. 
  #if Z is an identity matrix and q = 0, the output should be the richness
  #if Z is an identity matrix and q = 1, the output should be exponential Shannon index
  
	stopifnot(nrow(Z) == ncol(Z))
	stopifnot(nrow(Z) == nrow(U))
    
	#Each element of U must be positive
	stopifnot(any(U >= 0))
	nsp <- nrow(Z) #Number of species
	nc  <- ncol(U) #Number of communities
	
	#Convert to relative abundance
	for (j in 1:nc) U[,j] <- U[,j]/sum(U[,j])
	
	Y   <- numeric(nc)  #Diversity index for nc communities
  
  for (j in 1:nc){ #This for loop can be optimized
    
    #Remove 0 relative abundance
    b   <- U[,j]
    
    #Get index of positive abundance
    w   <- which(b > 0)
    w   <- as.integer(w)
    b   <- b[w]
    
    #Change Z
    Zc  <- Z[w,w]
    Zp  <- Zc %*% b  #n x 1 matrix
    
    if (q != 1){
      Y[j] <- ( sum(b*Zp**(q-1)) )**( 1/(1-q) )
    }else{
      cff  <- prod(Zp**b)
      Y[j] <- 1/cff
    }
  }
	return(Y)
}


#An example on hypothetical communities
Nsp   <- 5
Edom  <- c(6, 1, 1, 1, 1)

#Construct similarity matrix
Dsim <- diag(1, nr=Nsp, nc=Nsp)

#Construct identity matrix
Di <- Dsim

Dsim[2,1] <- Dsim[1,2] <- .4
Dsim[3,1] <- Dsim[1,3] <- .2
Dsim[4,1] <- Dsim[1,4] <- .2
Dsim[5,1] <- Dsim[1,5] <- .2
Dsim[3,2] <- Dsim[2,3] <- .2
Dsim[4,2] <- Dsim[2,4] <- .2
Dsim[5,2] <- Dsim[2,5] <- .2
Dsim[4,3] <- Dsim[3,4] <- .8
Dsim[5,3] <- Dsim[3,5] <- .4
Dsim[5,4] <- Dsim[4,5] <- .4

#Converting species abundances to matrices
Edom <- matrix(Edom,  nc=1)

#Compare richness estimates based on identity matrix and similarity matrix (r0 vs. r1)
r0   <- qDz(Edom,  Di,   q=0)
r1   <- qDz(Edom,  Dsim, q=0)

q0   <- qDz(Edom,  Di,   q=1)
q1   <- qDz(Edom,  Dsim, q=1)
