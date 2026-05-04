


###################   meta-analysis functions   #####################################

effect.sizes <- function(study){
  check.metaGEM.study(study)
  summ <- t( apply( study$expr, 1, getES, g=study$class ) )
  summ <- data.frame( summ, keys=study$keys )
  return(summ)
}


getES <- function(v, g){
  stopifnot( identical( length(v), length(g) ) )
  
  x <- cleanNA( v[ which(g==1) ] )
  y <- cleanNA( v[ which(g==0) ] )
  
  # Checking sample size
  n1 <- length(x); n2 <- length(y)
  if( n1 < 2 | n2 < 2 )
    return( c(n1=NA, m1=NA, sd1=NA,
              n2=NA, m2=NA, sd2=NA,
              diff=NA, pooled.sd=NA,
              g=NA, se.g=NA) )
  
  # Means and difference
  m1   <- mean(x); m2 <- mean(y)
  diff <- m1 - m2
  
  # Standard deviations
  sd1  <- sd(x);  sd2 <- sd(y)
  sp   <- sqrt( ( (n1-1)*sd1^2 + (n2-1)*sd2^2 )/( n1 + n2 - 2 ) )
  
  # Hedges' g correction
  cf   <- 1 - 3/( 4*(n1 + n2) - 9 )
  g    <- cf * diff/sp
  se.g <- sqrt( (n1+n2)/(n1*n2) + 0.5*g^2 /(n1+n2-3.94) )
  
  return( c(n1=n1, m1=m1, sd1=sd1,
            n2=n2, m2=m2, sd2=sd2,
            diff=diff, pooled.sd=sp,
            g=g, se.g=se.g) )
  
}


summ.eff.within <- function(effects, option="fixed.iv"){
  stopifnot( all( c("g", "se.g", "keys") %in% colnames(effects) ) )
  effects <- effects[ , c("g", "se.g", "keys")]
  effects$keys  <- as.character(effects$keys)
  
  if( length( grep(",", effects$keys) ) > 0 ) stop("Multiple keys detected. Please expand geneID first")
  
  if(nrow(effects)==1){
    rownames(effects) <- effects$keys
    effects <- effects[ , c("g", "se.g")]
    return(effects)
  }
  
  ## Deal with singletons
  singles.keys <- names(which(table(effects$keys) == 1))
  singles.ind <- which(effects$keys %in% singles.keys)
  out <- effects[ singles.ind, ]
  rownames(out) <- out$keys;  out$keys <- NULL
  
  ## Next, deal with multiple keys within a study
  multis <- effects[ -singles.ind, ]
  multis$abs.z <- abs( multis$g/multis$se.g )
  
  if(nrow(multis) > 0){  
    
    tmp <- split(multis, multis$keys)
    
    if (option == "fixed.iv") {
      out2 <- sapply(tmp, function(m) {
        unlist(meta.summaries(m$g, m$se.g, method = "fixed")[c("summary", "se.summary")])
      })
      out2 <- t(out2)
      colnames(out2) <- c("g", "se.g")
    }
    
    if (option == "extreme") {
      out2 <- lapply(tmp, function(mat) mat[which.max(mat$abs.z), ])
      out2 <- do.call(rbind, out2)
      out2 <- out2[, c("g", "se.g")]
    }
    out <- rbind(out, out2)
  }
  
  out <- out[sort(rownames(out)), ]
  return(out)
}


## Modified to include various measures of heterogeneity in effect sizes
## 1. Q.het - Cochrane's Q
## 2. df.het - degrees of freedom (number of studies a gene is measured in)
## 3. pval.het - is heterogeneity significant?
pool.inverseVar <- function( g, se.g, method ){
  stopifnot( identical( rownames(g), rownames(se.g) ) )
  out <- matrix( nr=nrow(g), nc=8,
                 dimnames=list( rownames(g), c("n.studies", "summary", "se.summary", "tau2", "p.value", "Q", "df", "pval.het") ) )
  
  for(j in 1:nrow(g)){
    
    e  <- cleanNA(    g[j, ] )
    se <- cleanNA( se.g[j, ] )
    n  <- length(e)
    
    if(n==1){
      summ <- e;   se.summ <- se;   tau2 <- NA
      Q.het = NA
      df.het = NA
      pval.het = NA
    } else {
      fit <- meta.summaries(e, se, method = method)
      summ <- fit$summary
      se.summ <- fit$se.summary
      tau2 <- ifelse( method=="fixed", NA, fit$tau2 )
      Q.het = fit$het[1]
      df.het = fit$het[2]
      pval.het = fit$het[3]
      rm(fit)
    }
    
    pval     <- 2*pnorm( abs(summ/se.summ), lower.tail=FALSE )
    
    out[j, ] <- c(n, summ, se.summ, tau2, pval, Q.het, df.het, pval.het)
    rm(e, se, n, summ, se.summ, tau2, pval, Q.het, df.het, pval.het)
  }
  return(out)
}


combine.effect.sizes <- function (list.of.effects, between.method="random", within.method="fixed.iv", everything=TRUE){
  
  if( is.null(names(list.of.effects)) )
    names(list.of.effects) <- paste("data", 1:length(list.of.effects), sep="")
  
  study.effects <- lapply(list.of.effects, function(effects) {
    
    effects <- data.frame(effects)
    effects$keys <- as.character(effects$keys)
    
    ## remove probes that cannot be mapped or have insufficient observations to calculate effect size
    bad <- which( is.na(effects$g) | is.na(effects$keys) | effects$keys=="NA" )
    effects <- effects[ setdiff(1:nrow(effects), bad), ]
    
    ## expand probes that maps to multiple keys
    effects <- expand.df( effects )
    
    ## summarize multiple probes within a study
    effects <- summ.eff.within(effects, option = within.method)
  })
  
  
  tmp <- multimerge(study.effects)
  g    <- tmp[, paste(names(study.effects), "_g", sep = ""), drop=FALSE]
  se.g <- tmp[, paste(names(study.effects), "_se.g", sep = ""), drop=FALSE]
  
  pooled.estimates <- data.frame( pool.inverseVar(g, se.g, method=between.method ) )
  
  if (everything) {
    return(list(g=g, se.g=se.g, pooled.estimates=pooled.estimates))
  } else {
    return(pooled.estimates)
  }
}


ttest.Pvalues <- function(study){
  check.metaGEM.study(study)
  summ <- get.ttest.P( study$expr, study$class )[ , c("P.up", "P.down")]
  summ <- data.frame( summ, keys=study$keys )
  return(summ)
}

get.ttest.P <- function(mat, g){
  ## test statistic and DF calculated using equal variance assumption
  
  tstat <- mt.teststat( mat, g, test="t.equalvar" )
  
  tstat <- apply(tmp$X, 1, function(x) {
    group1 <- x[tmp$classlabel == 1]
    group2 <- x[tmp$classlabel == 2]
    wilcox.test(group1, group2, exact = FALSE)$p.value
  })
  
  df <- length(g) - 2
  
  P.both <- 2*pt( abs(tstat), df=df, lower=FALSE )
  P.down <- pt( tstat, df=df, lower=TRUE )
  P.up   <- pt( tstat, df=df, lower=FALSE )
  
  out <- cbind(P.both, P.down, P.up)
  rownames(out) <- rownames(mat)
  return(out)
}

combine.significances <- function(list.of.sigs){
  study.sigs <- lapply(list.of.sigs, function(sigs){
    sigs <- data.frame(sigs)
    sigs$keys <- as.character(sigs$keys)
    
    ## remove probes that cannot be mapped or have insufficient observations to calculate the p-value
    bad <- which( (is.na(sigs[ ,1]) & is.na(sigs[ ,2]) ) | is.na(sigs$keys) | sigs$keys=="NA" )
    sigs <- sigs[ setdiff(1:nrow(sigs), bad), ]
    
    ## expand probes that maps to multiple keys
    sigs <- expand.df( sigs )
    
    ## summarize multiple probes within a study
    out <- summ.sigs.within(sigs)
  })
  
  tmp <- multimerge(study.sigs)
  sigs.up   <- tmp[ , grep("\\.up$", colnames(tmp), v=TRUE), drop=FALSE]
  sigs.down <- tmp[ , grep("\\.down$", colnames(tmp), v=TRUE), drop=FALSE]
  
  return(list(sigs.up=sigs.up, sigs.down=sigs.down))
}

sum.of.logs <- function(list.of.sigs){
  
  combsigs  <- combine.significances( list.of.sigs )
  
  sigs.up   <- combsigs$sigs.up
  valid.up  <- rowSums( !is.na(sigs.up) )
  F.stat.up <- -2*rowSums( log(sigs.up), na.rm=TRUE )
  F.pval.up <- pchisq( F.stat.up, 2*valid.up, lower.tail=FALSE )
  
  sigs.down   <- combsigs$sigs.down
  valid.down  <- rowSums( !is.na(sigs.down) )
  F.stat.down <- -2*rowSums( log(sigs.down), na.rm=TRUE )
  F.pval.down <- pchisq( F.stat.down, 2*valid.down, lower.tail=FALSE )
  
  out <- cbind(F.stat.up, F.pval.up, F.stat.down, F.pval.down)
  return(out)
}

summ.sigs.within <- function(sigs){
  pval.cols <- setdiff( 1:ncol(sigs), grep("keys", colnames(sigs)) )
  out <- sapply( sigs[ , pval.cols], function(x)
    tapply(x, sigs$keys, min) )
  return(out)
}


se <- function(x) sd(x)/length(x)

cleanNA <- function(x) return( x[!is.na(x) & is.finite(x) ] )

expand.df <- function (df, key.name = "keys", keys.sep = ",") {
  keys <- as.character(df[ ,key.name])
  skey <- strsplit( keys, split = keys.sep)
  df   <- df[rep(1:nrow(df), sapply(skey, length)), ]
  df[, key.name] <- unlist(skey)
  return(df)
}


multimerge <- function(mylist){
  unames <- unique( unlist( lapply( mylist, rownames ) ) )
  n      <- length(unames)
  
  out <- lapply( mylist, function(df){
    tmp <- matrix( nr=n, nc=ncol(df),
                   dimnames=list( unames, colnames(df) ) )
    tmp[ rownames(df), ] <- as.matrix(df)
    return(tmp)
  })
  
  bigout <- do.call( cbind, out )
  colnames(bigout) <- paste(rep( names(mylist), sapply(mylist, ncol) ),
                            sapply(mylist, colnames), sep="_")
  return(bigout)
}


pairwise.apply <- function(x, FUN, ...){
  n <- nrow(x)
  r <- rownames(x)
  output <- matrix(NA, nc=n, nr=n, dimnames=list(r, r))
  
  
  for(i in 1:n){
    for(j in 1:n){
      if(i >= j) next()
      output[i, j] <- FUN( x[i,], x[j,], ... )
    }
  }
  return(output)
}


check.metaGEM.study <- function(study){
  stopifnot( all( c("expr", "class", "keys") %in% names(study) ) )
  
  if( !all(levels(as.factor(study$class)) == c("0", "1")) )
    stop("study$class must be coded as 0 or 1")
  
  if( !is.character(study$keys) )
    stop("The keys must be stored as a character vector")
  
}

createAnnTable <- function(gem.data) {
  annTable = cbind(rownames(gem.data$expr), gem.data$keys)
  colnames(annTable) = c("probeid", "symbol")
  return(annTable)
}


filterCombinedES <- function(pooled.ES, summary=NULL, fdr=0.05, studies=0) {
  w <- which(pooled.ES$p.fdr < fdr & pooled.ES$n.studies >= studies)
  pooled.ES = pooled.ES[w,]
  if(!is.null(summary) && summary > 0) {
    pooled.ES = pooled.ES[which(pooled.ES$summary > summary),]
  } else if(!is.null(summary) && summary < 0) {
    pooled.ES = pooled.ES[which(pooled.ES$summary < summary),]
  }
  return(pooled.ES)
}


extractExprData <- function(genes, gems) {
  exprs = NULL
  for(i in 1:length(gems)) {
    cat("Processing data set ", i, "...", sep="")
    tempExprs = NULL
    junk = apply(as.matrix(genes), 1, function(x, keys) which(keys == x), keys=gems[[i]]$keys)
    for(j in 1:length(junk)) {
      if(length(junk[[j]]) == 0) {
        next
      }
      temp = gems[[i]]$expr[junk[[j]],]
      if(!is.vector(temp)) {
        temp = t(as.matrix(colMeans(temp) ))
      } else {
        temp = t(as.matrix(temp))
      }
      rownames(temp) = genes[j]
      tempExprs = rbind(tempExprs, temp)
    }
    tempExprs = data.frame(tempExprs)
    tempExprs$ID = rownames(tempExprs)
    if(i == 1) {
      exprs = tempExprs
    } else {
      exprs = merge(exprs, tempExprs, by="ID",all.x=T)
    }
    cat("Done.\n")
  }
  rownames(exprs) = exprs[,1]
  exprs = exprs[,2:dim(exprs)[2]]
  return(exprs)
}


extractDataFromGEM <- function(gem, genes) {
  tempExprs = NULL
  junk = lapply(as.matrix(genes), function(x, keys) which(keys == x), keys=gem$keys)
  for(j in 1:length(junk)) {
    if(length(junk[[j]]) == 0) {
      next
    }
    temp = gem$expr[junk[[j]],]
    if(!is.vector(temp)) {
      temp = t(as.matrix(colMeans(temp) ))
    } else {
      temp = t(as.matrix(temp))
    }
    rownames(temp) = genes[j]
    tempExprs = rbind(tempExprs, temp)
  }
  tempExprs = data.frame(tempExprs)
  
  return(tempExprs)
}

extractExprDataAsList <- function(genes, gems) {
  exprsList = NULL
  for(i in 1:length(gems)) {
    cat("Processing data set ", i, "...", sep="")
    tempExprs = NULL
    junk = apply(as.matrix(genes), 1, function(x, keys) which(keys == x), keys=gems[[i]]$keys)
    for(j in 1:length(junk)) {
      if(length(junk[[j]]) == 0) {
        next
      }
      temp = gems[[i]]$expr[junk[[j]],]
      if(!is.vector(temp)) {
        temp = t(as.matrix(colMeans(temp) ))
      } else {
        temp = t(as.matrix(temp))
      }
      rownames(temp) = genes[j]
      tempExprs = rbind(tempExprs, temp)
    }
    tempExprs = data.frame(tempExprs)
    exprsList[[i]] = tempExprs
    cat("Done.\n")
  }
  return(exprsList)
}

extractClassData <- function(gems) {
  classData = NULL
  for(i in 1:length(gems)) {
    temp = gems[[i]]$class
    classData = c(classData, temp)
  }
  return(classData)
}

createExprDataFor <- function(gem, genes) {
  tempExprs = NULL
  junk = apply(as.matrix(genes), 1, function(x, keys) which(keys == x), keys=gem$keys)
  for(j in 1:length(junk)) {
    if(length(junk[[j]]) == 0) {
      next
    }
    temp = gem$expr[junk[[j]],]
    if(!is.vector(temp)) {
      temp = t(as.matrix(colMeans(temp) ))
    } else {
      temp = t(as.matrix(temp))
    }
    rownames(temp) = genes[j]
    tempExprs = rbind(tempExprs, temp)
  }
  classData = gem$class
  
  exprData = tempExprs
  return(list(expr=exprData, y=classData))
}

leaveOneOutMetaAnalysisWrapper <- function(gems) {
  looResults = list()
  for(i in 1:length(gems)) {
    cat("Iteration: ", i, "...\n", sep="")
    looResults[[i]] = leaveOneOutMetaAnalysis(gems[-i])
    cat("Done.\n\n")
  }
  return(looResults)
}

runMetaAnalysis <- function(gems) {
  return(leaveOneOutMetaAnalysis(gems))
}

leaveOneOutMetaAnalysis <- function(gems) {
  annDB = createAnnTable(gems[[1]])
  
  if(length(gems) > 1) {
    for(i in 2:length(gems)) {
      tempAnnTable = createAnnTable(gems[[i]])
      
      commonProbes = match(annDB[,1], tempAnnTable[,1])
      commonProbes = commonProbes[!is.na(commonProbes)]
      if(length(commonProbes) > 0) {
        cat("Found common probes in", i, "\n", sep=" ")
        tempAnnTable = tempAnnTable[-commonProbes,]
      }
      annDB = rbind(annDB, tempAnnTable)
    }
  }
  rownames(annDB) = annDB[,1]
  annDB=as.matrix(annDB[,2])
  colnames(annDB) = c("symbol")
  
  cat("Computing effect sizes...")
  all.ES <- lapply( gems, effect.sizes ) 
  output.REM <- combine.effect.sizes( all.ES )
  
  cat("\nComputing pooled.ES...")
  pooled.ES <- output.REM$pooled.estimates
  pooled.ES$p.fdr <- p.adjust( pooled.ES$p.value, method="fdr" )
  pooled.ES <- pooled.ES[ order(pooled.ES$p.fdr), ]
  
  cat("\nComputing Fisher's output, then FDR-adjusting")
  all.Pvals <- lapply(gems, ttest.Pvalues)
  output.Fisher <- sum.of.logs(all.Pvals)
  output.Fisher <- adjust.fisher(output.Fisher=output.Fisher)
  cat("\nCounting votes...")
  cat("\n")
  return(list(gems = gems, annDB = annDB, all.ES = all.ES, output.REM = output.REM,
              pooled.ES = pooled.ES, all.Pvals = all.Pvals, output.Fisher = output.Fisher))
}

adjust.fisher <- function(output.Fisher, method="fdr"){
  F.Qval.up   <- p.adjust(output.Fisher[, "F.pval.up"], method=method)
  F.Qval.down <- p.adjust(output.Fisher[, "F.pval.down"], method=method)
  return(cbind(output.Fisher, F.Qval.up, F.Qval.down))
}


plotES <- function(allESs, col, names, minimum.x=NULL, maximum.x=NULL, legend=T, ...) {
  max.y = 0
  min.x = 0
  max.x = 0
  for(i in 1:length(allESs)) {
    if(max.y < max(density(allESs[[i]][,9], na.rm=TRUE)$y, na.rm=TRUE))
      max.y = max(density(allESs[[i]][,9], na.rm=TRUE)$y, na.rm=TRUE)
    if(min.x > min(density(allESs[[i]][,9], na.rm=TRUE)$x, na.rm=TRUE))
      min.x = min(density(allESs[[i]][,9], na.rm=TRUE)$x, na.rm=TRUE)
    if(max.x < max(density(allESs[[i]][,9], na.rm=TRUE)$x, na.rm=TRUE))
      max.x = max(density(allESs[[i]][,9], na.rm=TRUE)$x, na.rm=TRUE)
  }
  
  if(!is.null(minimum.x)) {
    min.x = minimum.x
  }
  if(!is.null(maximum.x)) {
    max.x = maximum.x
  }
  
  i = 1
  d = density(allESs[[i]][,9], na.rm=TRUE)
  plot(d, col=col[i], xlim=c(min.x, max.x), ylim=c(0, max.y), ...)
  for(i in 2:length(allESs)) {
    d = density(allESs[[i]][,9], na.rm=TRUE)
    lines(d, col=col[i], ...)
  }
  if(legend) {
    legend("topright",names,lty=1,col=col)
  }
  abline(v=0)
}


createMatrixFromGEM <- function(gem) {
  out = rbind(gem$class, gem$expr)
  rownames(out)[1] = "group"
  out = cbind(c("", gem$keys), out)
  return(out)
}

countSamples <- function(gems) {
  totalSamples = 0
  totalCases = 0
  for(i in 1:length(gems)) {
    totalSamples = totalSamples + length(gems[[i]]$class)
    totalCases = totalCases + sum(gems[[i]]$class)
  }
  totalControls = totalSamples - totalCases
  return(list(samples=totalSamples, cases=totalCases, controls=totalControls))
}


################  thresholding functions   ##########################################
examineLOOResults <- function(loo.results, threshold.ES, threshold.Studies = 0, threshold.Fisher) {
  
  pooled.ES = loo.results$pooled.ES
  output.Fisher = loo.results$output.Fisher
  
  commonRows = intersect(rownames(pooled.ES), rownames(output.Fisher))
  
  out <- cbind( pooled.ES[commonRows,], output.Fisher[commonRows, ] )
  out.pos <- out[ which(out$summary > 0), ]
  out.neg <- out[ which(out$summary < 0), ]
  
  junk.pos = filterCombinedES(out.pos, summary = 0.00001, fdr = threshold.ES, studies = threshold.Studies)
  junk2.pos = junk.pos[which(junk.pos$p.fdr <= threshold.ES & junk.pos$F.pval.up <= threshold.Fisher),]
  junk2.pos = junk2.pos[order(junk2.pos$summary, decreasing=T),]
  
  junk.neg = filterCombinedES(out.neg, summary = -0.00001, fdr = threshold.ES, studies = threshold.Studies)
  junk2.neg = junk.neg[which(junk.neg$p.fdr <= threshold.ES & junk.neg$F.pval.down <= threshold.Fisher),]
  junk2.neg = junk2.neg[order(junk2.neg$summary, decreasing=F),]
  
  return (list( loo.results=loo.results, siggenes.pos.ES = junk.pos, siggenes.pos.ESFisher = junk2.pos,
                siggenes.neg.ES = junk.neg, siggenes.neg.ESFisher = junk2.neg))
}


replaceValues <- function(x, thresholdValue = 0, replaceValue = 1) {
  for(i in 1:dim(x)[2]) {
    indices = which(x[,i] <= thresholdValue)
    if(length(indices) == 0)
      next
    x[indices,i] = replaceValue
  }
  return(x)
}

replaceNaNs <- function(x, replaceValue=1) {
  for(i in 1:dim(x)[2]) {
    indices = which(is.nan(x[,i]) == TRUE)
    if(length(indices) == 0)
      next
    x[indices,i] = replaceValue
  }
  return(x)  
}

replaceNAs <- function(x, replaceValue=1) {
  for(i in 1:dim(x)[2]) {
    indices = which(is.na(x[,i]) == TRUE)
    if(length(indices) == 0)
      next
    x[indices,i] = replaceValue
  }
  return(x)  
}


filterLOOResults <- function(discoveryResults, discoveryLOOResults, esThreshold, fisherThreshold, 
                             thresholdStudies = length(discoveryLOOResults)-1) {
  
  #Filter each LOO given significance levels
  discoveryLOOResultsFiltered = list()
  for(i in 1:length(discoveryLOOResults)) {
    cat("Filtering ", i, "...", sep="")
    discoveryLOOResultsFiltered[[i]] = examineLOOResults(discoveryLOOResults[[i]], 
                                                         threshold.ES = esThreshold, 
                                                         threshold.Studies = thresholdStudies, 
                                                         threshold.Fisher = fisherThreshold)
    cat("Done.\n")
  }
  return(discoveryLOOResultsFiltered)
}

getNonHetGenesLOO <- function(discoveryResults, discoveryLOOResultsFiltered, pValHet=0.05){
  # Find genes that are expressed across all leave-one-out meta analyses 
  pos.genes = rownames(discoveryLOOResultsFiltered[[1]]$siggenes.pos.ESFisher)
  for(i in 2:length(discoveryLOOResultsFiltered)) {
    pos.genes <- intersect(rownames(discoveryLOOResultsFiltered[[i]]$siggenes.pos.ESFisher), pos.genes)
  }
  pos.genes = sort(pos.genes)
  neg.genes = rownames(discoveryLOOResultsFiltered[[1]]$siggenes.neg.ESFisher)
  for(i in 2:length(discoveryLOOResultsFiltered)) {
    neg.genes <- intersect(rownames(discoveryLOOResultsFiltered[[i]]$siggenes.neg.ESFisher), neg.genes)
  }
  
  neg.genes = sort(neg.genes)
  
  ## Now remove genes with significant heterogeneity
  pos.genes = rownames(subset(discoveryResults$pooled.ES[pos.genes, ], pval.het > pValHet))
  neg.genes = rownames(subset(discoveryResults$pooled.ES[neg.genes, ], pval.het > pValHet))
  
  return(list(pos=pos.genes, neg=neg.genes))
}


forest.plot <- function(db, key, sort.names=T, print.labels=T, main=key, summlabel="Summary",
                        boxCol="blue", linesCol="lightblue", 
                        summaryCol="orange", textCol="red", boxsize=1, ...){
  stopifnot(identical(rownames(db$g), rownames(db$se.g)))
  stopifnot(identical(rownames(db$g), rownames(db$pooled.estimates)))
  
  g    <- cleanNA( db$g[ key, ] )
  if(sort.names) {
    g    <- g[ sort(names(g)) ]
  }
  names(g) = gsub("_g", "", names(g))
  
  se.g <- cleanNA( db$se.g[ key, ] )
  if(sort.names) {
    se.g <- se.g[ sort(names(se.g)) ]
  }
  names(se.g) = gsub("_se.g", "", names(se.g))
  
  ## stopifnot( identical( names(g), names(se.g) ) )
  
  pool    <- db$pooled.estimates[ key, "summary" ]
  se.pool <- db$pooled.estimates[ key, "se.summary" ]
  
  x.label <- "Standardized Mean Difference (log2 scale)"
  if(!print.labels) {
    names(g) <- NULL
    x.label <- ""
    summlabel <- ""
  } 
  
  metaplot( g, se.g, labels=names(g),
            summn=pool, sumse=se.pool, sumnn=1/se.pool^2, summlabel=summlabel,
            xlab=x.label, ylab="", main=main,
            colors=meta.colors(box=boxCol, lines=linesCol,
                               zero="black", summary=summaryCol, text=textCol), boxsize=boxsize, ... )
}

ESthresh <- function(discoveryResults, pos.genes, neg.genes, ESfold=1.5){
  ESthr <- log2(ESfold)
  pos <- discoveryResults$pooled.ES[pos.genes, ]
  neg <- discoveryResults$pooled.ES[neg.genes, ]
  pos <- pos[abs(pos$summary)>=ESthr, ]
  neg <- neg[abs(neg$summary)>=ESthr, ]
  return(list(pos=rownames(pos), neg=rownames(neg)))
}

geomMean <- function (x, na.rm = FALSE) 
{
  if (!is.numeric(x) && !is.complex(x) && !is.logical(x)) {
    warning("argument is not numeric or logical: returning NA")
    return(as.numeric(NA))
  }
  if (na.rm) 
    x <- x[!is.na(x)]
  if (any(x < 0)) 
    stop("'x' contains negative value(s)")
  return(prod(x)^(1/length(x)))
}

################   forward search functions   #############################################
getGenesData <- function(genes.mtx, genes){
  tmp <- genes.mtx[genes, ]
  rownames(tmp) <- genes
  tmp[is.na(tmp)] <- 1
  return(tmp)
}

## Here is the main function
forwardSearchWeighted <- function(discoveryResults, pos.genes, neg.genes, 
                                  yes.pos=NULL, yes.neg=NULL, forwardThresh=0.5){
  
  discovery.genes <- convertDiscoveryListToGenes(discoveryResults, pos.genes, neg.genes)  
  
  matchAll(discovery.genes, pos.genes, neg.genes, yes.pos, yes.neg)
  
  weights <- unlist(lapply(discovery.genes, function(GEM) length(GEM$class)))
  
  if(is.null(yes.pos) & is.null(yes.neg)){
    auc.weighted <- 0
  } else {
    pos.genes <- setdiff(pos.genes, yes.pos)
    neg.genes <- setdiff(neg.genes, yes.neg)
    auc.weighted <- getWeightedAUCsGenesList(discovery.genes, yes.pos, yes.neg, print=F)
  }
  
  while(TRUE){
    geneSearch <- forwardSearchInner(discovery.genes, pos.genes, neg.genes, 
                                     yes.pos, yes.neg, forwardThresh, 
                                     auc.weighted, weights)
    
    #Keep the gene that increases AUC the most, as long as it's > forwardThresh
    best <- max(geneSearch$diff, na.rm=T)
    cat("next best=", best, "\n")
    if(best >= forwardThresh){
      ##If multiple "best" genes, remove the first only
      bestGene <- as.character(geneSearch[geneSearch$diff==best, 1][1]) 
      bestAUC <- geneSearch[geneSearch$diff==best, 3][1]
      cat("Adding ", as.character(bestGene), bestAUC, "\n")
      if(bestGene %in% pos.genes){
        pos.genes <- setdiff(pos.genes, bestGene)
        yes.pos <- c(yes.pos, bestGene)
      } else {
        neg.genes <- setdiff(neg.genes, bestGene)
        yes.neg <- c(yes.neg, bestGene)
      }
    } else {break} ## stop when best<forwardThresh
    auc.weighted <- bestAUC
  }
  return(list(yes.pos, yes.neg))
}


## Function for inner loop in forwardSearchWeighted
## Calls the helper functions forwardSearchPos and forwardSearchNeg
forwardSearchInner <- function(discovery.genes, pos.genes, neg.genes, 
                               yes.pos, yes.neg, forwardThresh, 
                               auc.weighted, weights){
  if(length(pos.genes>0)){
    geneSearchPos <- data.frame(genes=pos.genes, orig=auc.weighted, search=0)
    auc <- forwardSearchPos(discovery.genes=discovery.genes, 
                            pos.test=pos.genes, yes.pos=yes.pos, yes.neg=yes.neg)
    geneSearchPos[ ,3] <- colSums(t(auc)*weights)
  } else {geneSearchPos <- data.frame(genes=NULL, orig=NULL, search=NULL)}   
  
  if(length(neg.genes>0)){
    geneSearchNeg <- data.frame(genes=neg.genes, orig=auc.weighted, search=0)
    auc <- forwardSearchNeg(discovery.genes=discovery.genes, 
                            neg.test=neg.genes, yes.pos=yes.pos, yes.neg=yes.neg) 
    geneSearchNeg[ ,3] <- colSums(t(auc)*weights)
  } else {geneSearchNeg <- data.frame(genes=NULL, orig=NULL, search=NULL)}  
  
  geneSearch <- rbind(geneSearchPos, geneSearchNeg)
  if(dim(geneSearch)[1]==0) break;
  geneSearch <- data.frame(geneSearch, diff=geneSearch$search - geneSearch$orig)
  
  return(geneSearch)
}

forwardSearchPos <- function(discovery.genes, pos.test, yes.pos, yes.neg){
  auc <- lapply(discovery.genes, function(GEM) {
    #posGenes <- GEM$genes[c(pos.test, yes.pos), ]
    posGenes <- getGenesData(GEM$genes, c(pos.test, yes.pos))
    
    negScore <- 0 #in case no neg genes
    if (sum(!is.na(match(yes.neg, rownames(GEM$genes)))) > 0){
      #negGenes <- GEM$genes[yes.neg, ]
      negGenes <- getGenesData(GEM$genes, yes.neg)
      negScore <- apply(negGenes, 2, geomMean)
    }
    
    #Since often fewer neg genes, can make them relatively less important
    ratio <- length(yes.neg)/(length(yes.pos)+1)
    
    auc <- pos.test #preallocate
    for(j in 1:length(pos.test)){
      #take score without each gene (LOO)
      tmp <- posGenes[c(pos.test[j], yes.pos), , drop=F]
      posScore <- apply(tmp, 2, geomMean)
      totalScore <- scale(posScore - ratio*negScore)          
      auc[j] <- efficientAUC(GEM$class, totalScore)
    }
    return(auc)
  }) 
  
  auc <- matrix(as.numeric(unlist(auc)), nrow=length(pos.test), 
                ncol=length(discovery.genes), byrow=F)
  auc[is.na(auc)] <- 0
  return(auc)
}


forwardSearchNeg <- function(discovery.genes, neg.test, yes.pos, yes.neg){
  auc <- lapply(discovery.genes, function(GEM) {
    posScore <- 0 #in case no pos genes
    if (sum(!is.na(match(yes.pos, rownames(GEM$genes)))) > 0){
      #posGenes <- GEM$genes[yes.pos, ]
      posGenes <- getGenesData(GEM$genes, yes.pos)
      posScore <- apply(posGenes, 2, geomMean)
    }
    
    #negGenes <- GEM$genes[c(neg.test, yes.neg), ]
    negGenes <- getGenesData(GEM$genes, c(neg.test, yes.neg))
    #Since often fewer neg genes, can make them relatively less important
    if(is.null(yes.pos)){ 
      ratio <- 1
    } else {
      ratio <- (length(yes.neg)+1)/length(yes.pos)
    }
    
    auc <- neg.test #preallocate
    for(j in 1:length(neg.test)){  
      #take score without each gene (LOO)
      tmp <- negGenes[c(neg.test[j], yes.neg), , drop=F]
      negScore <- apply(tmp, 2, geomMean)
      totalScore <- scale(posScore - ratio*negScore)
      auc[j] <- efficientAUC(GEM$class, totalScore)
    }    
    return(auc)
  }) 
  auc <- matrix(as.numeric(unlist(auc)), nrow=length(neg.test), 
                ncol=length(discovery.genes),  byrow=F)
  auc[is.na(auc)] <- 0
  return(auc)
}


getWeightedAUCs <- function(discovery.genes, pos.genes, neg.genes){
  counter <- 1
  weights <- integer()
  auc <- lapply(discovery.genes$gems, function(GEM) {
    GEMmin <- min(GEM$expr, na.rm=T)
    if (GEMmin <0) {GEM$expr <- GEM$expr + abs(GEMmin)}
    sink("trash")
    totalScore <- ScoreNegPos(GEM, pos.genes, neg.genes)
    sink()
    auc <- efficientAUC(GEM$class, totalScore)
    cat(discovery.genes$GSEnames[counter], auc, "\n")
    weights <<- c(weights, length(GEM$class))
    counter <<- counter + 1
    return(auc)
  }); 
  AUC <- as.numeric(unlist(auc))
  AUC <- AUC * weights
  cat("weighted AUC:", sum(AUC), "\n")
  return(sum(AUC))
}

efficientAUC <- function(labels, predictions) {
  #Taken from ROCR package
  levels <- sort(unique(labels))
  labels <- ordered(labels, levels = levels)
  
  n.pos <- sum(labels == levels[2])
  n.neg <- sum(labels == levels[1])
  
  pred.order <- order(predictions, decreasing = TRUE)
  predictions.sorted <- predictions[pred.order]
  
  tp <- cumsum(labels[pred.order] == levels[2])
  fp <- cumsum(labels[pred.order] == levels[1])
  dups <- rev(duplicated(rev(predictions.sorted)))
  tp <- c(0, tp[!dups])
  fp <- c(0, fp[!dups])
  fn <- n.pos - tp
  tn <- n.neg - fp
  
  tpr <- tp / (tp + fn)
  fpr <- fp / (fp + tn)
  roc = data.frame(x = fpr, y = tpr)
  roc <- roc[order(roc$x, roc$y),]
  
  i <- 2:nrow(roc)
  auc <- (roc$x[i] - roc$x[i - 1]) %*% (roc$y[i] + roc$y[i - 1])/2
  #auc <- signif(auc, 4)
  return(auc)
}

convertDiscoveryListToGenes <- function(discoveryResults, pos.genes, neg.genes){
  cat("Converting probes:genes for all GEMs")
  discovery.genes <- list()
  for (i in 1:length(discoveryResults$gems)){
    #Convert to genes and subset
    GEM.genes <- extractCleanFromGEM(discoveryResults$gems[[i]], c(pos.genes, neg.genes) )
    
    #scale positive
    GEMmin <- min(GEM.genes, na.rm=T)
    if (GEMmin <0) {GEM.genes <- GEM.genes + abs(GEMmin)}
    
    discovery.genes[[discoveryResults$GSEnames[i]]] <- list(genes = GEM.genes, class = discoveryResults$gems[[i]]$class)
    cat(".")
  }
  rm(GEM.genes, i, GEMmin)
  cat("   done. \n")
  return(discovery.genes)
}

ScoreGenesMtx <- function(GeneMtx, pos.genes, neg.genes){
  posScore <- 0 
  if (sum(pos.genes %in% rownames(GeneMtx)) > 0){
    posMatch <- getGenesData(GeneMtx, pos.genes)
    posScore <- apply(posMatch, 2, geomMean)
  }
  negScore <- 0
  if (sum(neg.genes %in% rownames(GeneMtx)) > 0){
    negMatch <- getGenesData(GeneMtx, neg.genes)
    negScore <- apply(negMatch, 2, geomMean)
  }
  
  ## Weight pos and neg by how many genes are being called
  if(length(posScore)>1){
    ratio <- length(neg.genes)/length(pos.genes); 
    negScore <- ratio*negScore
    totalScore <- scale(posScore - negScore)
  } else {
    totalScore <- scale(-negScore)
  }
  
  return(totalScore)
}

getAUCsGenesList <- function(genes.list, yes.pos, yes.neg, print=T){
  if(length(c(yes.pos,yes.neg))==0) return(0.5)
  counter <- 1
  auc <- lapply(genes.list, function(GEM) {
    score <- ScoreGenesMtx(GEM$genes, yes.pos, yes.neg)
    auc <- signif(efficientAUC(GEM$class, score), 3)
    if(print) cat(names(discovery.genes)[counter], auc, "\n")
    counter <<- counter + 1
    return(auc)
  }); 
  mean.AUC <- mean(as.numeric(unlist(auc)))
  if(print) cat("mean AUC:", mean.AUC, "\n")
  return(mean.AUC)
}

extractCleanFromGEM <- function(GEM, genes){
  GenesMtx <- extractDataFromGEM(GEM, genes)
  GenesMtx <- replaceValues(GenesMtx, 0, 1)
  GenesMtx <- replaceNaNs(GenesMtx, 1)
  GenesMtx <- replaceNAs(GenesMtx, 1)
  return(GenesMtx)
}

getWeightedAUCsGenesList <- function(genes.list, yes.pos, yes.neg, print=T){
  if(length(c(yes.pos,yes.neg))==0) return(0.5)
  weights <- unlist( lapply(genes.list, function(GEM) length(GEM$class)) )
  counter <- 1
  auc <- lapply(genes.list, function(GEM) {
    score <- ScoreGenesMtx(GEM$genes, yes.pos, yes.neg)
    auc <- efficientAUC(GEM$class, score)
    counter <<- counter + 1
    return(auc)
  }); 
  auc <- as.numeric(unlist(auc))
  weightedAUC <- auc * weights
  if(print) cat("sum weighted AUC:", sum(weightedAUC), "\n")
  return(sum(weightedAUC))
}

matchAll <- function(discovery.genes, pos.genes, neg.genes, yes.pos, yes.neg){
  discGenes <- unique(unlist(lapply(discovery.genes, function(x) rownames(x$genes))))
  if( !(all(yes.pos %in% pos.genes) )) stop("not all yes.pos %in% pos.genes")
  if( !(all(yes.neg %in% neg.genes) )) stop("not all yes.neg %in% neg.genes")
  if( !(all(neg.genes %in% discGenes) )) stop("not all neg.genes %in% discovery.genes")
  if( !(all(pos.genes %in% discGenes) )) stop("not all pos.genes %in% discovery.genes")
  if( !(all(discGenes %in% c(pos.genes, neg.genes)) )) stop("not all discovery.genes %in% c(pos.genes, neg.genes)")
}
