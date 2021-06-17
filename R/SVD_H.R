library(tidyverse)
library(QZ)
library(onion)

######################## METHODS FOR QUATERNIONIC SVD ########################
q2cdq<-function(q){ #This function transforms the quaternion using the Cayley Dickson notation
  #input :
  #q = a quaternion  as a tibble with columns w, x, y, z
  #output :
  #a list of 2 tibbles containing a complex number
  A1<- tibble(cplx = complex(real = q %>% pull(w),imaginary = q %>% pull(x)))
  A2<- tibble(cplx = complex(real = q %>% pull(y),imaginary = q %>% pull(z)))
  return(list(A1,A2))
}

cdq2q<-function(q1,q2){ #This function transforms a quaternion written with Cayley Dickson notation 
                        #in a "classic" quaternion represented by a tibble of shape 1x4
  #Input : 
  #q1 and q2 : two tibbles containing one complex number
  q<-tibble(w=Re(q1[[1]]), x=Im(q1[[1]]),y=Re(q2[[1]]),z=Im(q2[[1]]))
  return(q)
}

MatQtoMatcdQ<-function(A){
  #This function transforms a quaternion matrix (tibble) in 2 complex matrixes (tibbles) 
  #using the Cayley Dickson notation 
  #input :
  #A = a quaternion matrix (tibble) with quaternions = tibbles with 1 line, 4 columns (w,x,y,z)
  #output :
  #a list of 2 tibbles containing A1 and A2 : complex 'matrixes' corresponding to A = A1 + A2.j
  N<-dim(A)[1]
  M<-dim(A)[2]
  A1<-tibble()
  A2<-tibble()
  for(i in 1:N){
    for(j in 1:M){
      A1[i,j]<-q2cdq(A[[i,j]][[1]])[1]
      A2[i,j]<-q2cdq(A[[i,j]][[1]])[2]
    }
  }
  return(list(A1,A2))
}

MatcdQtoMatQ<-function(A1,A2){
  #A1 and A2 are such that A1+A2.j=A with A being a quaternion matrix
  #This function gives A as a tibble (same shape as A1 and A2) of tibbles representing quaternions
  #Input : 
  #A1 and A2 : complex matrixes of same size
  if(!identical(dim(A1),dim(A2))) stop("A1 and A2 should have the same shape")
  N<-dim(A1)[1]
  M<-dim(A1)[2]
  #A<-matrix(,nrow=N,ncol=M)
  A<-tibble(.rows = N)
  for(i in 1:N){
    for(j in 1:M){
      A[[i,j]]<-cdq2q(q1 = A1[i,j],q2=A2[i,j])
    }
  }
 
  return(A)
}


MatHtoC<-function(A1,A2){
  #Given a quaternion matrix A with expression A = A1 + A2.j
  #with size NxM, and A1 and A2 complex tibble
  #build the associated complex representation tibble ChiA
  
  #input : A1 and A2, complex tibbles with same shape NxM
  #output : Chia, complex tibble of size 2Nx2M
  
  if(!identical(dim(A1),dim(A2))) stop("A1 and A2 should have the same shape")
  
  N<-dim(A1)[1]
  M<-dim(A1)[2]
  
  ChiA<-A1
  ChiA<-bind_cols(ChiA, A2)
  temp<-cbind(-t(H(as.matrix(A2))), t(H(as.matrix(A1))))
  ChiA<-as_tibble(rbind(as.matrix(ChiA), as.matrix(temp)))
  #ChiA<-bind_rows(as.tibble(-H(as.matrix(A2)))

  
  return(ChiA)
}

Hmatprod<-function(A1,A2,B1,B2){
  
    #Function that computes the matrix multiplication
    #for 2 quaternion matrices A=A1+A2j and B=B1+B2j
    
    #Input:
    #-----
    #A1 and A2, complex matrices such that A=A1+A2j
    #B1 and B2, complex matrices such that B=B1+B2j
    
    #Output:
    #------
    #AB1 and AB2, complex matrices such that A.B = AB1 + AB2j
  
  AB1<-as.matrix(A1)%*%as.matrix(B1)-as.matrix(A2)%*%t(H(as.matrix(B2)))
  AB2<-as.matrix(A1)%*%as.matrix(B2)+as.matrix(A2)%*%t(H(as.matrix(B1)))
  return(list(AB1,AB2))

}

########################## METHODS TESTING ##########################
quat<-tibble(w=1, x=6, y=2*w,z=3*y)
q2cdq(quat) #quat écrit en notation de Cayley-Dickson

qts<-tibble(w=1:5, x=6:10, y=2*w,z=3*y) #quaternion time series i.e. quaternion vector

q11<-qts %>% slice(1) #getting each quaternion from qts as 1 tibble
q12<-qts %>% slice(2)
q21<-qts %>% slice(3)
q22<-qts %>% slice(4)
q1<-list(q11,q21) #creating 2 lists q1,q2 of 2 quaternions in order to format them as a tibble of quaternions
q2<-list(q12,q22)
MatQ<-tibble(q1,q2) #MatQ is a tibble of quaternions (representing a matrix 2x2 of quaternions)


q3<-list(tibble(w=12, x=65, y=2*w,z=3*y), tibble(w=1, x=6, y=2*w,z=3*y))
MatQ2<-MatQ %>% add_column(q3=q3)
MatQ2 #tibble MatQ with one more column

MatQtoMatcdQ(MatQ) #Transorming quaternion matrix MatQ into 2 complex matrixes considerinf the Cayley-Dickson notation

MatHtoC(MatQtoMatcdQ(MatQ)[[1]],MatQtoMatcdQ(MatQ)[[2]]) #Creation of the adjoint complex matrix of MatQ (as a tibble)

Hmatprod(MatQtoMatcdQ(MatQ)[[1]],MatQtoMatcdQ(MatQ)[[2]],MatQtoMatcdQ(MatQ)[[1]],MatQtoMatcdQ(MatQ)[[2]]) #Square MatQ as its Cayley-Dickson notation

#Checking the Hmatprod method's result
a<-quaternion(Re = 1, i = 6, j = 2, k = 6)
b<-quaternion(Re = 2, i = 7, j = 4, k = 12)
c<-quaternion(Re = 3, i = 8, j = 6, k = 18)
d<-quaternion(Re = 4, i = 9, j = 8, k = 24)

a*a+b*c
a*b+b*d
c*a+d*c
c*b+d*d

################### SVDH METHOD #########################

SVDH<-function(A1,A2){
  
    #Function that computes the SVD of a quaternion 
    #matrix A = A_1 + A_2 j
    
    #input: A_1 and A_2, complex matrices of size NxM
    
    #output: 
    #U1,U2: two complex matrices of left singular vectors forming the quaternion matrix 
    #U = U1 + U2 j
    #Vh1,Vh2: two complex matrices of right singular vectors forming the quaternion matrix 
    #Vh = Vh1 + Vh2 j
    #D: real singular values (dimension = min(N,M))
    
    #Testing shape of A1, A2
  
  if(!identical(dim(A1),dim(A2))) stop("A1 and A2 should have the same shape")
  
  #Forming the complex representation
  ChiA<-MatHtoC(A1,A2)
  
  #Computing the SVD of ChiA
  svd<-svd(as.matrix(ChiA))
  Uc<-svd$u
  Dc<-svd$d
  Vc<-t(H(svd$v)) 
  
  #allocating singular values of ChiA to those of A
  D<-Dc[seq(1,length(Dc),2)]
  
  #allocating singular vectors of ChiA to those of A
  UC<-Uc[,seq(1,dim(Uc)[2],2)]
  VC<-Vc[,seq(1,dim(Vc)[2],2)]
  U1<-UC[1:N,]
  U2<- -t(H(UC[(N+1):(2*N),]))
  
  V1<-VC[1:M,]
  V2<- -t(H(VC[(M+1):(2*M),]))
  
  Vh1<-t(H(V1))
  Vh2<- -V2
  
  return(list('U1'=U1,'U2'=U2,'D'=D,'V1'=Vh1,'V2'=Vh2))
}

############ TESTING SVDH METHOD ON EXAMPLE 1 FROM PEI ET AL #################

Q_tib<-tribble(~w,~x,~y,~z,1,1,1,1,2,1,0,-1,1,0,-1,2,3,2,-2,1)
Q<-tribble(~c1,~c2, Q_tib %>% slice(1),Q_tib %>% slice(2),Q_tib %>% slice(3),Q_tib %>% slice(4))
Q #Quaternion matrix of shape 2x2  (exemple Pei et al.)
Qcd<-MatQtoMatcdQ(Q) #Q with Cayley Dickson notation

N<-2 #number of row in Q
M<-2 #number of columns in Q
svd<-SVDH(as.matrix(Qcd[[1]]),as.matrix(Qcd[[2]])) #computing the SVDH on Q (using Cayley Dickson notation)

k<-min(N,M)
DD<-matrix(0i,N,M)
diag(DD)<-complex(real=svd$D)
h<-Hmatprod(svd$U1,svd$U2,DD,matrix(0i,N,M))
UD1<-h[[1]] #multiplication of matrix U (left singular vectors) with D (singular values)
UD2<-h[[2]]
h2<-Hmatprod(UD1,UD2,H(svd$V1),H(svd$V2)) #multiplication of UD with transpose-conjugate V (right singular vectors)
Res1<-h2[[1]]
Res2<-h2[[2]] 
#Res1 and Res2 so that Res1 + Res2.j corresponds to the original matrix Q using the SVD

Q1<-as.matrix(Qcd[[1]])
Q2<-as.matrix(Qcd[[2]])

#Calculating the error between original matrix Q and recomposed matrix using the SVDH

n1<-Matrix::norm(Re(A1-Res1))**2+Matrix::norm(Im(A1-Res1))**2
n2<-Matrix::norm(Re(A2-Res2))**2+Matrix::norm(Im(A2-Res2))**2
sqrt(n1+n2)

################# TESTING THE SVDH METHOD AS DONE IN LE BIHAN AND MARS #####################

#Creation of the matrixes and computation of the SVD
N<-10
M<-10
A1<-matrix(complex(real = rnorm(N*M, mean=0,sd=1),imaginary = rnorm(N*M, mean=0,sd=1)),N,M)
A2<-matrix(complex(real = rnorm(N*M, mean=0,sd=1),imaginary = rnorm(N*M, mean=0,sd=1)),N,M)
svd<-SVDH(A1,A2)

#Recomposition of the matrix Q using multiplication of matrixes resulting from the SVDH
k<-min(N,M)
DD<-matrix(0i,N,M)
diag(DD)<-complex(real=svd$D)
h<-Hmatprod(svd$U1,svd$U2,DD,matrix(0i,N,M))
UD1<-h[[1]]
UD2<-h[[2]]
h2<-Hmatprod(UD1,UD2,H(svd$V1),H(svd$V2))
Res1<-h2[[1]]
Res2<-h2[[2]]

#Calculating the error between original matrix and its recomposition using SVDH
n1<-Matrix::norm(Re(A1-Res1))**2+Matrix::norm(Im(A1-Res1))**2
n2<-Matrix::norm(Re(A2-Res2))**2+Matrix::norm(Im(A2-Res2))**2
sqrt(n1+n2)

################## PCA ON QUATERNIONS METHOD #######################@

qPCA<-function(X, nb.axis){
  #Input :
  #X : a quaternion tibble
  #nb.axis : the number of axis to keep in PCA according to the user
  #Output : a list of several elements
  #total_eigenval : a numeric containing all the eigenvalues of X
  #eigenval : a numeric containing the nb.axis first eigenvalues
  #total_left.eigenvect : a tibble of quaternions with the left eigenvectors as columns
  #left.eigenvect : a tibble of quaternions with the nb.axis first left eigenvectors as columns
  #total_right.eigenvect : a tibble of quaternions with the right eigenvectors as columns
  #right.eigenvect : a tibble of quaternions with the nb.axis first right eigenvectors as columns
  #intertia : the percentage of inertia represented by the nb.axis first axis
  
  N<-dim(X)[1]
  M<-dim(X)[2]
  k<-min(N,M)
  if(nb.axis>k) stop("nb.axis can't be greater than the minimum dimension of X")
  cdX<-MatQtoMatcdQ(X)
  svd<-SVDH(cdX[[1]],cdX[[2]])
  
  eigenval<-svd$D
  left.eigenvect<-MatcdQtoMatQ(svd$U1,svd$U2)
  right.eigenvect<-MatcdQtoMatQ(svd$V1,svd$V2)
  
  inertia<-sum(eigenval[1:nb.axis])/sum(eigenval)
  
  return(list('total_eigenval'=eigenval,'eigenval'=eigenval[1:nb.axis],'total_left.eigenvect'=left.eigenvect,
  'left.eigenvect'=left.eigenvect[1:nb.axis],'total_right.eigenvect'=right.eigenvect,
  'right.eigenvect'=right.eigenvect[1:nb.axis],'inertia'=inertia*100))
}
qpca<-qPCA(X = Q,nb.axis=1)
