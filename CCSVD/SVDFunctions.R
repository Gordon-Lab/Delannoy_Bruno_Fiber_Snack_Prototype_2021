############################ 
#
#
#
# 1. Compute the Fat SVD
#
#
#
############################
svdFatOg <- function(m){
  # Given a matrix m, return the fat form svd
  # We do this because R only returns the thin form svd
  # This makes it difficult to do spectral cleaning
  # 1. For the case of rows < columns
  nrows = nrow(m)
  ncols = ncol(m)
  l <- list()
  if(nrows < ncols){
    s <- svd(m)
    u <- s$u
    sigma <- diag(s$d)
    v <- s$v
    sigma2 <- cbind(sigma, matrix(0, nrow = nrow(m), ncol = ncol(m) - nrow(m)))
    v2t <- t(v)
    v2t <- rbind(v2t, matrix(0, nrow = ncol(v2t) - nrow(v2t), ncol = ncol(v2t)))
    v2 <- t(v2t)
    l[['u']] <- u
    l[['s']] <- sigma2
    l[['v']] <- v2
    return(l)
  } else if(nrows > ncols){
    s <- svd(m)
    u <- s$u
    sigma <- diag(s$d)
    v <- s$v
    sigma2 <- rbind(sigma, matrix(0, nrow = nrow(m)-ncol(m), ncol = ncol(m)))
    u2 <- cbind(u, matrix(0, nrow = nrow(m), ncol = nrow(m)-ncol(m)))
    l[['u']] <- u2
    l[['s']] <- sigma2
    l[['v']] <- v
    return(l)
  } else{
    s <- svd(m)
    l[['u']] <- s$u
    l[['s']] <- diag(s$d)
    l[['v']] <- s$v
    return(l)
  }
}

svdFat <- function(m){
  # Given a matrix m, return the fat form svd
  # We do this because R only returns the thin form svd
  # This makes it difficult to do spectral cleaning
  # 1. For the case of rows < columns
  nrows = nrow(m)
  ncols = ncol(m)
  l <- list()
  s <- svd(m, nu = nrows, nv = ncols)
  u <- s$u
  v <- s$v
  sigma <- diag(s$d)
  if(nrows < ncols){
    sigma2 <- cbind(sigma, matrix(0, nrow = nrows, ncol = (ncols - nrows)))
    l[['u']] <- u
    l[['s']] <- sigma2
    l[['v']] <- v
    return(l)
  } else if(nrows > ncols){
    sigma2 <- rbind(sigma, matrix(0, nrow = nrow(m)-ncol(m), ncol = ncol(m)))
    l[['u']] <- u
    l[['s']] <- sigma2
    l[['v']] <- v
    return(l)
  } else{
    s <- svd(m)
    l[['u']] <- s$u
    l[['s']] <- diag(s$d)
    l[['v']] <- s$v
    return(l)
  }
}
############################ 
#
#
#
# 2. Spectrally clean some matrix
#
#
#
############################
svdSpectralClean <- function(m, indices){
  # Given a matrix m, clean out various columns from the eigenspectrum, and return the new svd
  # 1. Compute the singular value decomposition
  svdecomp <- svdFat(m)
  # 2. Assign the various decomposed matrices
  U <- svdecomp$u
  S <- svdecomp$s
  V <- svdecomp$v
  # 3. Remove relevant indices
  U1 <- as.matrix(U[,-indices])
  S1 <- as.matrix(S[-indices, -indices])
  V1 <- as.matrix(V[, -indices])
  # 4. Recompute m
  M <- U1 %*% S1 %*% t(V1)
  # 5. Recompute SVD
  svdecomp1 <- svdFat(M)
  return(svdecomp1)
}

spectralClean <- function(m, indices){
  # Given a matrix m, clean out various columns from the eigenspectrum, and return the new matrix
  # 1. Compute the singular value decomposition
  svdecomp <- svdFat(m)
  # 2. Assign the various decomposed matrices
  U <- svdecomp$u
  S <- svdecomp$s
  V <- svdecomp$v
  # 3. Remove relevant indices
  U1 <- as.matrix(U[,-indices])
  S1 <- as.matrix(S[-indices, -indices])
  V1 <- as.matrix(V[, -indices])
  # 4. Recompute m
  M <- U1 %*% S1 %*% t(V1)
  return(M)
}

############################ 
#
#
#
# 3. Choose top percentile
#
#
#
############################
choosePercentile <- function(l, percentile, tail = 'both'){
  # Given a list of numbers, a percentile, and a choice for the top/bottom/both tail, return the values of that list
  l <- l[order(l, decreasing = TRUE)]
  if(tail == 'top'){
    return(head(l, floor(percentile*length(l))))
  } else if(tail == 'bot'){
    return(rev(tail(l, floor(percentile*length(l)))))
  } else{
    return(c(head(l, floor(percentile*length(l))),
             tail(l, floor(percentile*length(l)))))
  }
}

############################ 
#
#
#
# 4. Plot Eigenspectrum
#
#
#
############################
plotEigenSpectrum <- function(s){
  # Takes as input a svd object and returns a ggplot2 object with the eigenspectrum
  s2 <- diag(s$s)^2
  pve <- s2/sum(s2)
  pve <- pve[pve > 0.00001]
  df <- as.data.frame(pve)
  colnames(df) <- 'Eig'
  g <- ggplot(df, aes(x = Eig)) + 
    geom_histogram(aes(y = ..count..), bins = nrow(df), fill = 'lightsteelblue4', color = 'lightsteelblue2', alpha = 0.8) + 
    geom_density(aes(y = ((df[1,1] - df[nrow(df),1])/(nrow(df)))*..count..), alpha = 0.5, color = 'cornflowerblue', fill = 'lightsteelblue2') +
    ggtitle("Eigenspectrum") + 
    labs(x = 'Variance Explained') +
    theme(text = element_text(family = "Fira Mono"),
          axis.title.x = element_text(colour = "DarkBlue", size = 10),
          axis.title.y = element_text(colour = "DarkBlue", size = 10),
          axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 8),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8),
          plot.title = element_text(colour = "DarkBlue", size = 12))
  return(g)
}

############################ 
#
#
#
# 5. Iterative Spectral Cleaning
#
#
#
############################
iterativeSpectralCleaning <- function(m, thresh){
  # Given a matrix, iteratively decompose it into lists of component parts
  # 1. Establish the random matrix r
  r <- apply(m, 2, function(x){sample(x, length(x))})
  svd.r <- svdFat(r)
  eval.r <- diag(svd.r$s)^2
  pcv.r <- eval.r/sum(eval.r)
  noiseThreshold <- pcv.r[1]
  
  # 2. Establish an end criteria, e.g. when the PC explains less than random matrix PC1
  s <- svdFat(m)
  eval <- diag(s$s)^2
  pcv <- eval/sum(eval)
  pcv <- pcv[pcv > noiseThreshold]
  
  # 3. Loop through each eigenvector
  eigenvectorholder <- 1:length(pcv)
  l <- list()
  if(length(pcv) == 1){
    i=1
    svc <- svdFat(m)
    ev1p <- svc$u[,1]
    names(ev1p) <- rownames(m)
    ev1d <- svc$v[,1]
    names(ev1d) <- colnames(m)
    # Extract top and bottom % row projections
    topp <- choosePercentile(ev1p, thresh, 'top')
    botp <- choosePercentile(ev1p, thresh, 'bot')
    topd <- choosePercentile(ev1d, thresh, 'top')
    botd <- choosePercentile(ev1d, thresh, 'bot')
    l[[paste('PC', i, sep = '')]] <- list('LeftFeatures' = c(topp, rev(botp)),
                                          'RightFeatures' = c(topd, rev(botd)))
  } else{
    for(i in 1:length(pcv)){
      svc <- svdSpectralClean(m, eigenvectorholder[-i])
      ev1p <- svc$u[,1]
      names(ev1p) <- rownames(m)
      ev1d <- svc$v[,1]
      names(ev1d) <- colnames(m)
      # Extract top and bottom % row projections
      topp <- choosePercentile(ev1p, thresh, 'top')
      botp <- choosePercentile(ev1p, thresh, 'bot')
      topd <- choosePercentile(ev1d, thresh, 'top')
      botd <- choosePercentile(ev1d, thresh, 'bot')
      l[[paste('PC', i, sep = '')]] <- list('LeftFeatures' = c(topp, rev(botp)),
                                            'RightFeatures' = c(topd, rev(botd)))
    }
  }
  return(l)
}


############################ 
#
#
#
# 6. Module Eigengene Cross-Correlation Analysis
#
#
#
############################
############################ a. Function for choosing significant projections
corTestListAgainstMatrix <- function(l, m, cormethod = 'pearson'){
  # Given a list and a matrix, correlate the list against the columns of the matrix
  cortests <- lapply(1:ncol(m), function(x){cor.test(l, m[, x], method = cormethod)})
  plist <- unlist(lapply(cortests, function(x){x[[3]]}))
  corlist <- unlist(lapply(cortests, function(x){x[[4]]}))
  names(plist) <- colnames(m)
  names(corlist) <- colnames(m)
  output <- list(); output[['Rho']] <- corlist; output[['pVals']] <- plist
  return(output)
}

############################ b. Implementation of MECCA
MECCA <- function(m, cormethod){
  # Given a matrix, iteratively decompose it into lists of component parts
  # 1. Establish the random matrix r, significant loadings onto u, and significant loadings onto v
  randomStats <- matrix(0, nrow = 100, ncol = 5)
  for(i in 1:100){
    r <- apply(m, 2, function(x){sample(x, length(x))})
    rownames(r) <- rownames(m)
    svd.r <- svdFat(r)
    eval.r <- diag(svd.r$s)^2
    pcv.r <- eval.r/sum(eval.r)
    noiseThreshold <- pcv.r[1]
    # Establish what a significant loading would look like onto u, using the 5% criteria choosing the 2.5% tail on either side
    r.u.loadings <- choosePercentile(svd.r$u[,1], 0.05, 'top')
    r.u.05top <- min(r.u.loadings)
    r.u.loadings <- choosePercentile(svd.r$u[,1], 0.05, 'bot')
    r.u.05bot <- max(r.u.loadings)
    # Establish what a significant loading would look like onto u
    r.v.loadings <- choosePercentile(svd.r$v[,1], 0.1, 'top')
    r.v.05top <- min(r.v.loadings)
    r.v.loadings <- choosePercentile(svd.r$v[,1], 0.1, 'bot')
    r.v.05bot <- max(r.v.loadings)
    # Add these to a dataframe to average later
    randomStats[i,] <- c(noiseThreshold, r.u.05top, r.u.05bot, r.v.05top, r.v.05bot)
  }
  # Set noise thresholds for the leading eigenvector % variance explained, the loadings onto U, and the loadings onto V
  pveThreshold <- median(randomStats[,1])
  UloadingsThresholdtop <- median(randomStats[,2])
  UloadingsThresholdbot <- median(randomStats[,3])
  VloadingsThresholdtop <- median(randomStats[,4])
  VloadingsThresholdbot <- median(randomStats[,5])
  
  # 3. Establish an end criteria, e.g. when the PC explains less than random matrix PC1
  s <- svdFat(m)
  eval <- diag(s$s)^2
  pve <- eval/sum(eval)
  pve <- pve[pve > pveThreshold]
  
  # 4. Loop through each eigenvector
  eigenvectorholder <- 1:length(pve)
  l <- list()
  for(i in 1:length(pve)){
    svc <- svdSpectralClean(m, eigenvectorholder[-i])
    ev1p <- svc$u[,1]
    names(ev1p) <- rownames(m)
    ev1d <- svc$v[,1]
    names(ev1d) <- colnames(m)
    # Extract top and bottom % row projections
    toprow <- ev1p[ev1p > UloadingsThresholdtop]
    toprow <- toprow[order(toprow, decreasing = TRUE)]
    botrow <- ev1p[ev1p < UloadingsThresholdbot]
    botrow <- botrow[order(botrow, decreasing = FALSE)]
    # Extract top and bottom % col projections
    topcol <- ev1d[ev1d > VloadingsThresholdtop]
    topcol <- topcol[order(topcol, decreasing = TRUE)]
    botcol <- ev1d[ev1d < VloadingsThresholdbot]
    botcol <- botcol[order(botcol, decreasing = FALSE)]
    l[[paste('PC', i, sep = '')]] <- list('Plasma' = c(toprow, rev(botrow)),
                                          'Duo' = c(topcol, rev(botcol)))
  }
  return(l)
}

MECCA2 <- function(m, rho){
  # Given a matrix, iteratively decompose it into lists of component parts
  # 1. Establish the random matrix r, significant loadings onto u, and significant loadings onto v
  randomStats <- list()
  for(i in 1:20){
    r <- apply(m, 2, function(x){sample(x, length(x))})
    svd.r <- svdFat(r)
    eval.r <- diag(svd.r$s)^2
    pcv.r <- eval.r/sum(eval.r)
    randomStats[i] <- pcv.r[1]
  }
  # Set noise thresholds for the leading eigenvector % variance explained, the loadings onto U, and the loadings onto V
  pveThreshold <- median(unlist(randomStats))

  # 3. Establish an end criteria, e.g. when the PC explains less than random matrix PC1
  s <- svdFat(m)
  eval <- diag(s$s)^2
  pve <- eval/sum(eval)
  pve <- pve[pve > pveThreshold]
  
  # 4. Loop through each eigenvector
  eigenvectorholder <- 1:length(pve)
  l <- list()
  for(i in 1:length(pve)){
    svc <- svdSpectralClean(m, eigenvectorholder[-i])
    ev1p <- svc$u[,1]
    names(ev1p) <- rownames(m)
    ev1d <- svc$v[,1]
    names(ev1d) <- colnames(m)
    # Extract top and bottom % row projections
    corrow <- cor(ev1d, t(m))
    toprow <- corrow[,corrow > rho]
    toprow <- toprow[order(toprow, decreasing = TRUE)]
    botrow <- corrow[,corrow < -rho]
    botrow <- botrow[order(botrow, decreasing = FALSE)]
    # Extract top and bottom % col projections
    corcol <- cor(ev1p, m)
    topcol <- corcol[,corcol > rho]
    topcol <- topcol[order(topcol, decreasing = TRUE)]
    botcol <- corcol[,corcol < -rho]
    botcol <- botcol[order(botcol, decreasing = FALSE)]
    l[[paste('PC', i, sep = '')]] <- list('DuoProteins' = c(toprow, rev(botrow)),
                                          'ASVs' = c(topcol, rev(botcol)))
  }
  return(l)
}

################################### Finding derivatives
findLocalPeaks2 <- function (x){
  # Given a list of points and a local environment of size m, find the local maxima
  # The algorithm will be as follows:
  # 1. Find 2nd derivative intervals
  # 2. Find the maximum point in each 2nd derivative interval
  # 3. Return the index for the maximum of each 2nd derivative interval
  secondDeriv <- sign(diff(diff(x))) # Two differences, make sure to add 2 the index for finding inflection
  inflections <- diff(secondDeriv)
  negToPos <- which(diff(secondDeriv) == -2)
  posToNeg <- which(diff(secondDeriv) == 2)
  # Find which interval to start with and create a running interval list padded by 0 and length(x)
  minimum <- min(negToPos, posToNeg)
  if(minimum %in% negToPos){
    suppressWarnings(intervalList <- c(rbind(negToPos, posToNeg)))
  } else{
    suppressWarnings(intervalList <- c(0, rbind(posToNeg, negToPos)))
  }
  # Create a table with intervals
  intervalList <- intervalList + 2 # have to add 2 because the 2nd derivative loses 2 indices
  maxVals <- list()
  for(i in seq(from = 1, to = length(intervalList), by = 2)){
    if(i + 1 > length(intervalList)){
      break
    } else{
      maxVals[[i]] <- max(x[intervalList[i]:intervalList[i+1]])
    }
  }
  maxVals <- unlist(maxVals)
  indices <- match(maxVals, x)
  return(indices)
}

############################ 
calculateGamma <- function(PC, Eigenvalues){
  # Given a list of PCs and Eigenvalues, calculate the scale-freeness
  #x1 <- PC[1:30]
  #y1 <- Eigenvalues[1:30]
  x1 <- PC
  y1 <- Eigenvalues
  return(lm(y1 ~ x1)$coefficients[2])
}

calcGammaVelocity <- function(m){
  # Given a matrix, iteratively spectrally clean it and calculate the change in gamma
  # 1. Establish the random matrix r
  r <- apply(m, 2, function(x){sample(x, length(x))})
  svd.r <- svdFat(r)
  eval.r <- diag(svd.r$s)^2
  pcv.r <- eval.r/sum(eval.r)
  noiseThreshold <- pcv.r[1]
  
  # 2. Establish an end criteria, e.g. when the PC explains less than random matrix PC1
  s <- svdFat(m)
  eval <- diag(s$s)^2
  pcv <- eval/sum(eval)
  pcv <- pcv[pcv > noiseThreshold]
  
  # 3. Loop through each eigenvector
  eigenvectorholder <- 1:length(pcv)
  l <- list()
  for(i in 1:length(pcv)){
    newM <- spectralClean(m, eigenvectorholder[-i])
    ev1p <- svc$u[,1]
    names(ev1p) <- rownames(m)
    ev1d <- svc$v[,1]
    names(ev1d) <- colnames(m)
    # Extract top and bottom % row projections
    topp <- choosePercentile(ev1p, thresh, 'top')
    botp <- choosePercentile(ev1p, thresh, 'bot')
    topd <- choosePercentile(ev1d, thresh, 'top')
    botd <- choosePercentile(ev1d, thresh, 'bot')
    l[[paste('PC', i, sep = '')]] <- list('LeftFeatures' = c(topp, rev(botp)),
                                          'RightFeatures' = c(topd, rev(botd)))
  }
  return(l)
}


