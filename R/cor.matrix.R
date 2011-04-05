cor.matrix <-
function(m1,m2,x,y,method="pearson",dist="euclidean",permutations=999,norm=FALSE)
	{
	permut<-permutations
	dist.y<-vegdist(y,method=dist)
	dist.x<-vegdist(x,method=dist)
	correlation<-cor(dist.x,dist.y,method = method)
	value<-matrix(NA,nrow=1,ncol=permut)
	for (i in 1:permut){
		m2.permut<-permut.row.matrix(m2)
		x.permut<-m1%*%m2.permut
		if (norm == "TRUE") {
     	   	matrix.permut <- apply(x.permut^2, 2, sum)
        	x.permut <- sweep(x.permut, 2, sqrt(matrix.permut), "/")
    		}
		dist.x.permut<-vegdist(x.permut,method=dist)
		cor.x.permut <- cor(dist.x.permut, dist.y,method = method)
		value[i]<-cor.x.permut
		}
	signific<-(sum(value>=correlation)+1)/(permut+1)
	return(list(Obs=correlation,p=signific))	
	}
