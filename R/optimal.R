optimal <-
function (comm, envir, traits, subset = 3, pattern="tcap",dist = "euclidean",method = "pearson", scale = TRUE, scale.envir = TRUE) {
    part.cor <- function(rxy, rxz, ryz) {
        (rxy - rxz * ryz)/sqrt(1 - rxz * rxz)/sqrt(1 - ryz * ryz)
    }
    colnames(traits) <- colnames(traits, do.NULL = FALSE, prefix = "T")
    if (scale.envir == "TRUE") {
        envir <- cent.norm(envir)
    }
    dist.y <- vegdist(envir, method = dist)
    m <- dim(traits)[2]
    if (subset > m) {
        stop("\n Subset must be lower than the number of traits\n")
    }
    PATTERNS <- c("tcap","tdap","tcap.tdap")
    pattern <- pmatch(pattern, PATTERNS)
    if (length(pattern)>1){
    		stop("\n Only one argument is accepted in pattern \n")
    	}
    if (is.na(pattern)){ 
        stop("\n Invalid pattern \n")
    	}	
    p <- 1:subset
    bin <- factorial(m)/(factorial(p) * factorial(m - p))
    comb <- matrix(NA, nrow = sum(bin), ncol = 1)
    for (i in 1:subset) {
        combinations <- combn(colnames(traits), i, simplify = TRUE)
        for (j in 1:bin[i]) {
            comb[(j + sum(bin[1:i - 1])), 1] <- paste(combinations[,j], collapse = " ")
        }
    }
    correlation <- matrix(NA, nrow = sum(bin), ncol = 1)
    for (i in 1:subset) {
        combinations1 <- combn(colnames(traits), i, simplify = TRUE)
        for (j in 1:bin[i]) {
            if (pattern==1){
            	T <- matrix.t(comm, as.matrix(traits[, combinations1[,j]]), scale = scale)
            	correlation[(j + sum(bin[1:i - 1])), 1] <- cor(vegdist(as.matrix(T$matrix.T),method = dist), dist.y, method = method)
            }
            if (pattern==2){
            	T <- matrix.t(comm, as.matrix(traits[, combinations1[,j]]), scale = scale)
            	X <- matrix.x(comm, as.matrix(traits[, combinations1[,j]]), scale = scale)
            	dist.x <- vegdist(X$matrix.X, method = dist)
            	dist.z <- vegdist(T$matrix.T, method = dist)
            	rxy <- cor(dist.x, dist.y, method = method)
            	rxz <- cor(dist.x, dist.z, method = method)
            	ryz <- cor(dist.y, dist.z, method = method)
            	correlation[(j + sum(bin[1:i - 1])), 1] <- part.cor(rxy,rxz, ryz)
        	}
        	if (pattern==3){
        		X <- matrix.x(comm, as.matrix(traits[, combinations1[,j]]), scale = scale)
            	correlation[(j + sum(bin[1:i - 1])), 1] <- cor(vegdist(as.matrix(X$matrix.X), method = dist), dist.y, method = method)
        	}
        }
    }
    result <- data.frame(Subset = comb, ro = correlation, stringsAsFactors = FALSE)
    return(result[order(result[, 2], decreasing = TRUE), ])
}