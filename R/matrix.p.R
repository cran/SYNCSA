matrix.p <-
function(comm,dist.spp)
	{
	matrix.w<-sweep(comm, 1, rowSums(comm), "/")
	similar.phy<-1-(dist.spp/max(dist.spp))
	matrix.phy<-1/colSums(similar.phy)
	matrix.q<-sweep(similar.phy,1,matrix.phy,"*")
	matrix.P<-matrix.w%*%matrix.q
	return(list(matrix.w=matrix.w,matrix.q=matrix.q,matrix.P=matrix.P))
		}

