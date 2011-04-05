cent.norm <-
function(x){
	x.cent<-sweep(x, 2, colMeans(x),"-")
	x.norm<- apply(x.cent^2, 2, sum)
	x.cent.norm <- sweep(x.cent, 2, sqrt(x.norm), "/")
	return(x.cent.norm)
	}

