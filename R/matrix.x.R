matrix.x <-
function(comm,traits,scale=TRUE)
	{
	matrix.w<-sweep(comm, 1, rowSums(comm), "/")
	dist.traits<-as.matrix(vegdist(traits,method="euclidean",diag=TRUE,upper=TRUE))
	similar.traits<-1-(dist.traits/max(dist.traits))
	matrix.traits<-1/colSums(similar.traits)
	matrix.u<-sweep(similar.traits,1,matrix.traits,"*")
	if (scale=="TRUE"){
		matrix.trait.s<- apply(as.matrix(traits^2), 2, sum)
		traits.s<- sweep(as.matrix(traits), 2, as.matrix(sqrt(matrix.trait.s)), "/")
		dist.traits<-vegdist(traits.s,method="gower",diag=TRUE,upper=TRUE)
		similar.traits<-1 - as.matrix(dist.traits/max(dist.traits))
		matrix.traits<-1/colSums(similar.traits)
		matrix.u<-sweep(similar.traits,1,matrix.traits,"*")
		}
	matrix.X<-matrix.w%*%matrix.u	
	return(list(matrix.w=matrix.w,matrix.u=matrix.u,matrix.X=matrix.X))
		}

