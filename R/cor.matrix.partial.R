cor.matrix.partial<-function (m1, m2, x, y, z, method = "pearson", dist = "euclidean", permutations = 999, norm = FALSE) {
    permut <- permutations
    part.cor <- function(rxy, rxz, ryz) {
        (rxy - rxz * ryz)/sqrt(1 - rxz * rxz)/sqrt(1 - ryz * 
            ryz)
    }
    dist.x <- vegdist(x, method = dist)
    dist.y <- vegdist(y, method = dist)
    dist.z <- vegdist(z, method = dist)
    rxy <- cor(dist.x, dist.y, method = method)
    rxz <- cor(dist.x, dist.z, method = method)
    ryz <- cor(dist.y, dist.z, method = method)
    statistic <- part.cor(rxy, rxz, ryz)
    	if((rxz==1|ryz==1)==TRUE){
    		statistic<-0	
    }
    value <- matrix(NA, nrow = 1, ncol = permut)
    for (i in 1:permut) {
        m2.permut <- permut.row.matrix(m2)
        x.permut <- m1 %*% m2.permut
        if (norm == "TRUE") {
            matrix.permut <- apply(x.permut^2, 2, sum)
            x.permut <- sweep(x.permut, 2, sqrt(matrix.permut), 
                "/")
        }
        dist.x.permut <- vegdist(x.permut, method = dist)
        rxy <- cor(dist.x.permut, dist.y, method = method)
        rxz <- cor(dist.x.permut, dist.z, method = method)
        value[i] <- part.cor(rxy, rxz, ryz)
     		if((rxz==1|ryz==1)==TRUE){
    		value[i]<-0	
    	}
    }
    signif <- (sum(value >= statistic) + 1)/(permut + 1)
    return(list(Obs = statistic, p = signif))
}