MCS<-function(data1=data1, nc=nc, method1="method1", method2="method2", index="index", print.stats=FALSE, st.data=FALSE,
     plot.hc=FALSE, circ=FALSE, convert=TRUE, plot.data=FALSE){
# Error messages
if(is.matrix(data1)==FALSE && is.numeric(data1)==FALSE) return(cat('Error: Data must be numeric string or matrix.\n'))
if(is.matrix(data1)) sample.size = length(data1[1,])*length(data1[,1])
if(is.matrix(data1)==FALSE) sample.size = length(data1)
if (nc<2 || nc>(sample.size-1)) return(cat('Error: Number of clusters must be between 2 and n-1.\n'))
if (method1!="ward.D2"&&method1!="single"&&method1!="complete"&&method1!="average"&&method1!="mcquitty"&&method1!="median"&&method1!=
"centroid"&&method1!="kmeans") return(cat('Error:',method1,'is not a plausible clustering method.\n'))
if (method2!="ward.D2"&&method2!="single"&&method2!="complete"&&method2!="average"&&method2!="mcquitty"&&method2!="median"&&method2!=
"centroid"&&method2!="kmeans") return(cat('Error:',method2,'is not a plausible clustering method.\n'))
if (index!="rand"&&index!="jaccard") return(cat('Error:',index,'is not a plausible similarity index.\n'))
################################################################################
nc2<- nc-1
y<- seq(2,nc)
membvecA<- vector(length=nc)
membvecB<- vector(length=nc)
MCS<- vector(length=nc2)
MCS_output<- matrix(data=NA, nrow=nc2, ncol=2)
k1<- 1
################### CIRCULAR DATA ##############################################
if(circ==TRUE){
  length.data1 <-length(data1)
  if (convert==TRUE) data1 <-(data1*pi)/180
  data1mat<- matrix(data1,nrow=length.data1,ncol=1,byrow=TRUE)
  if (plot.data==TRUE) circ.plot(data1mat, main="Plot of Circular Data")
    df1<- array(dim=c(length.data1,length.data1))
    for (ii in 1:length.data1){
       for (jj in 1:length.data1){
          df1[ii,jj]<- 0.5*(1 - cos(data1mat[ii] - data1mat[jj]))
       }
    }
  distmat<- as.dist(df1)
}
################### NON-CIRCULAR DATA ##########################################
if(circ==FALSE){
  if(st.data==FALSE) distmat<- dist(data1,method="euclidean")
  ################## STANDARDIZATION OF THE DATA ###############################
  if (st.data==TRUE) {
    numcols <-length(data1[1,])
    numrows <-length(data1[,1])
    standat<- as.matrix(data1,nrow=numrows,ncol=numcols,byrow=TRUE)
    colsums<- apply(standat,2,sum)
    colavg<- (colsums/numrows)
    colvar<- apply(standat,2,var)
    colstdev<- sqrt(colvar)
    standat1<- array(dim=c(numrows,numcols))
    for(ii in 1:numrows){
      for(jj in 1:numcols) standat1[ii,jj]<- (standat[ii,jj] - colavg[jj])/colstdev[jj]
    }
    distmat<- dist(standat1,method="euclidean")
  }
}
################################################################################
for(i in 1:nc2){
k1<- i + 1
################################## Method 1 ####################################
{if (method1=="kmeans") {
   cl<- kmeans(distmat,k1,iter.max = 10, nstart = 5)
   methodA<- as.vector(cl$cluster) }
 else {
   hc<-hclust(distmat,method=method1)
   if(i==1) {if(plot.hc==TRUE) plot(hc)}
   methodA<- cutree(hc, k1)
   methodA<-as.vector(methodA)} }
################################## Method 2 ####################################
{if (method2=="kmeans") {
   cl<- kmeans(distmat,k1,iter.max = 10, nstart = 5)
   methodB<- as.vector(cl$cluster) }
 else {
   hc<-hclust(distmat,method=method2)
   if(circ==FALSE){
      if(i==1) {if(method1!=method2){if(plot.hc==TRUE) plot(hc)}}}
   methodB<- cutree(hc, k1)
   methodB<-as.vector(methodB) }  }
########### Calculating the Matching counts matrix and the marginal totals #####
nij<- table(methodA,methodB)
nr<- max(methodA)
nc<- max(methodB)
n<- length(methodA)
nij<- matrix(nij,nrow=nr,ncol=nc,byrow=TRUE)
choose2 <- function(v){
            out <- numeric(0)
            for (i in 1:length(v)) out[i] <- ifelse(v[i] >= 2,
                choose(v[i], 2), 0)
            out
        }
M<- choose2(n) ### M=a+b+c+d = (n choose 2)
Totsumsqr<- sum((nij*nij))
rows<- apply(nij,1,sum)
cols<- apply(nij,2,sum)
sumpsqr<- (sum(rows**2))
P<- sumpsqr - n
sumqsqr<- (sum(cols**2))
Q<- sumqsqr - n
########## Calculating corrected Rand index Using Hubert and Arabie expectation #########
if(index=="rand") {
  cluster.stats<- function (d, clustering, alt.clustering = NULL){
  R<- 1 - ( P + Q + 2*n)/(2*M) + (Totsumsqr/M)
  ExpRand<- 1 -( P + Q + 2*n)/(2*M) + (1/M)*((sumpsqr*sumqsqr)/(n*(n-2)) +(n**2 - (sumpsqr + sumqsqr))/(n-1))
  CR<- ( R - ExpRand )/( 1.0 - ExpRand )
  out<- CR
  out
  }
  output1<- cluster.stats(distmat,membvecA,alt.clustering=membvecB)
}
#### Calculating corrected Jaccard index using Taylor series expansion with Hubert and Arabie expectation ###
if(index=="jaccard") {
  cluster.stats<- function (d, clustering, alt.clustering = NULL){
  Jaccard<- (Totsumsqr - n)/(sumpsqr + sumqsqr - n - Totsumsqr)
  EJN<- ((1/(2*M))*(P+n)*(Q+n) + (n^2/(n-1))- (1/(n-1))*(P+Q+2*n) - n)
  EJD<- (P+Q+n) - ((1/(2*M))*(P+n)*(Q+n) + (n^2/(n-1))- (1/(n-1))*(P+Q+2*n))
  EJ<- EJN/EJD
  CJ<- (Jaccard - EJ)/(1 - EJ)
  out<- CJ
  }
  output1<- cluster.stats(distmat,membvecA,alt.clustering=membvecB)
}
################################################################################
MCS[i]<- output1
MCS_output[i,1]<- k1
MCS_output[i,2]<- output1
if (print.stats==TRUE)cat('No. of Clusters = ',k1,'\t','MCS = ',MCS[i],'\n')}
if(index=="rand") plot(y,MCS, xlab="No. of Clusters", ylab="MCS using corrected Rand",type="l",pch=1,col='black', lwd=2,lty=2)
if(index=="jaccard") plot(y,MCS, xlab="Number of Clusters", ylab="MCS using corrected Jaccard",type="l",pch=1,col='black', lwd=2,lty=2)
}
