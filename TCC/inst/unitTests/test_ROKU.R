## Kadota original code
kadota_2006_bmc_bioinformatics <- function(x){ 
    x <- x[(!is.na(x))]
    x_length <- length(x)
    x <- x[(x != 0)]
    if(sum(x) <= 0) {
        return(log(x_length, base = 2)) 
    } else if(sd(x) == 0) {
        return(log(x_length, base = 2))
    } else {
        x_prime <- abs(x - affy::tukey.biweight(x))
        p <- x_prime / sum(x_prime) 
        e <- sum(p * log(p, base = 2)) 
        return(-e)
    }
}
## Kadota original code
kadota_2003_physiol_genomics_0.25 <- function(x){ 
    if(length(x) == sum(is.na(x))){
        x <- c(rep(0, length(x)))
    } else if(length(x) == sum(is.nan(x))){
        x <- c(rep(0, length(x)))
    }
    x_org <- x 
    x <- x[(!is.na(x))] 
    x <- x[(!is.nan(x))] 
    n_plus_s <- length(x)
    x.sort <- sort(x)
    x.order <- order(x)
    maice_Ut <- 0
    maice_i <- 0
    maice_j <- 0
    flag <- c(rep(0, length = n_plus_s))
    if (sd(x) != 0) {
        for (i in 1:(n_plus_s * 0.25 + 1)) {
            for (j in 1:(n_plus_s - i)) {
                if ((i + j - 2) <= n_plus_s * 0.25) {
                    n <- (n_plus_s + 1 - j) - i + 1
                    s <- n_plus_s - n
                    set_sd <- sd(x.sort[i:(n_plus_s + 1 - j)]) * sqrt((n - 1) / n)
                    Ut <- n * log(set_sd) + sqrt(2) * s * lfactorial(n) / n
                    if(maice_Ut > Ut){
                        maice_Ut <- Ut
                        maice_i <- i
                        maice_j <- j
                    }
                }
            }
        }
        if (maice_i > 1) {
            flag[x.order[1:(maice_i - 1)]] <- -1
        }
        if (maice_j > 1) {
            flag[x.order[(n_plus_s + 1 - maice_j + 1):n_plus_s]] <- 1
        }
        tmp <- replace(x_org, ((!is.nan(x_org)) & (!is.na(x_org))), flag) 
        return(tmp)
    } else {
        tmp <- replace(x_org, ((!is.nan(x_org)) & (!is.na(x_org))), flag)
        return(tmp)
    }
}


test_ROKU <- function() {
    x <- matrix(rnorm(100), ncol = 10)
    y <- som::normalize(x)

    roku.tcc <- ROKU(x)
    roku.kdt <- t(apply(y, 1, kadota_2003_physiol_genomics_0.25))
    checkEqualsNumeric(roku.kdt, roku.tcc[, 1:ncol(x)])
    
    outl.kdt <- apply(x, 1, kadota_2006_bmc_bioinformatics)
    checkEqualsNumeric(outl.kdt, roku.tcc[, ncol(x) + 1])

    colnames(x) <- paste("t", 1:ncol(x))
    rownames(x) <- paste("g", 1:nrow(x))
    roku.tccnm <- ROKU(x)
    checkEqualsNumeric(roku.tcc, roku.tccnm)
}
