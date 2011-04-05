matrix.t <-
function(comm,traits,scale=TRUE)
	{
	matrix.w<-sweep(comm, 1, rowSums(comm), "/")
	matrix.b<-traits
	matrix.T<-matrix.w%*%matrix.b
	if (scale=="TRUE"){
		matrix.traits <- apply(matrix.T^2, 2, sum)
		matrix.T<- sweep(matrix.T, 2, sqrt(matrix.traits), "/")
	}
	return(list(matrix.w=matrix.w,matrix.b=matrix.b,matrix.T=matrix.T))
		}

