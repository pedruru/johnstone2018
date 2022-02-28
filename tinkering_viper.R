#######
#' Null model by sample permutation testing
#'
#' This function performs sample permutation and DESeq2 to generate a null model
#'
#' @param x Matrix containing the test dataset
#' @param y Matrix containing the reference dataset
#' @param colData (from DESeq2) for matrix input: a \code{DataFrame} or \code{data.frame} with at least a single column.
#' Rows of colData correspond to columns of countData
#' @param design (from DESeq2) a \code{formula} or \code{matrix}.
#' the \code{formula} expresses how the counts for each gene
#' depend on the variables in \code{colData}. Many R \code{formula} are valid,
#' including designs with multiple variables, e.g., \code{~ group + condition},
#' and designs with interactions, e.g., \code{~ genotype + treatment + genotype:treatment}.
#' See \code{\link{results}} for a variety of designs and how to extract results tables.
#' By default, the functions in this package will use 
#' the last variable in the formula for building results tables and plotting.
#' \code{~ 1} can be used for no design, although users need to remember
#' to switch to another design for differential testing
#' @param refcol (from DESeq2) the reference column in \code{colData}
#' @param reflvl (from DESeq2) the reference level¨
#' @param resname (from DESeq2) the name of the individual effect (coefficient) for
#' building a results table. Use this argument rather than \code{contrast}
#' for continuous variables, individual effects or for individual interaction terms.
#' The value provided to \code{name} must be an element of \code{resultsNames(object)}.
#' @param contrast (from DESeq2) this argument specifies what comparison to extract from
#' the \code{object} to build a results table. one of either:
#' \itemize{
#'  \item a character vector with exactly three elements:
#' the name of a factor in the design formula,
#' the name of the numerator level for the fold change,
#' and the name of the denominator level for the fold change
#' (simplest case)
#'  \item a list of 2 character vectors: the names of the fold changes
#' for the numerator, and the names of the fold changes
#' for the denominator.
#' these names should be elements of \code{resultsNames(object)}.
#' if the list is length 1, a second element is added which is the
#' empty character vector, \code{character()}.
#' (more general case, can be to combine interaction terms and main effects)
#'  \item a numeric contrast vector with one element
#' for each element in \code{resultsNames(object)} (most general case)
#' }
#' If specified, the \code{name} argument is ignored.
#' @param per Integer indicating the number of permutations
#' @param repos Logical, whether the permutations should be performed with reposition
#' @param seed Integer indicating the seed for the permutations, 0 for disable it
#' @param cores Integer indicating the number of cores to use (set to 1 in windows systems)
#' @param verbose Logical, whether progression messages should be printed in the terminal
#' @param ... Additional parameters added to keep compatibility
#' @return List of matrices of statistics and p-values with genes in rows and permutations in columns
#' @method DESeq2Null
#' @export
#' @seealso \code{\link{msviper}}, \code{\link{viper}}

DESeq2Null <- function(x, y, colData, design, refcol=NULL, reflvl=NULL, resname=NULL, contrast=NULL, per=1000, repos=TRUE, seed=1, cores=1, verbose=TRUE) {

  if (!requireNamespace("DESeq2", quietly = TRUE))
    BiocManager::install("DESeq2")
  
  if (seed>0) set.seed(round(seed))
  pb <- NULL
  if (verbose) {
    message(date(), "\nComputing the null model distribution by ", per, " permutations.")
  }
  if (cores>1) {
    res <- mclapply(1:per, function(i, x, y, repos) {
      expset <- cbind(x, y)
      repeat{
        sorder <- sample(ncol(expset), replace=repos)
        if (length(unique(sorder[1:ncol(x)]))>1 & length(unique(sorder[-(1:ncol(x))]))>1) break
        if (verbose) message("-", appendLF=FALSE)
      }
      x1 <- filterColMatrix(expset, sorder[1:ncol(x)])
      y1 <- filterColMatrix(expset, sorder[-(1:ncol(x))])
      ## largo <- rowSums(!is.na(x1))
      ## largoy <- rowSums(!is.na(y1))
      ## t <- ((rowMeans(x1, na.rm=TRUE) - rowMeans(y1, na.rm=TRUE))/sqrt(((largo - 1) *  rowVars(x1) + (largoy - 1) * rowVars(y1))/(largo + largoy - 2))/sqrt(1/largo + 1/largoy))[, 1]
      countData <- cbind(x1, y1)
      colData <- colData[colnames(countData), ]
      dds <- DESeq2::DESeqDataSetFromMatrix(countData, colData, design)
      if (!is.null(refcol) & !is.null(reflvl)) {
        colData(dds)[, refcol] <- relevel(colData(dds)[, refcol], reflvl)
      }
      dds <- DESeq2::DESeq(dds)
      if (!is.null(resname)) {
        dds <- data.frame(DESeq2::results(dds, name = resname)[, c("stat", "pvalue")])
      } else if (!is.null(contrast)) {
        dds <- data.frame(DESeq2::results(dds, contrast = contrast)[, c("stat", "pvalue")])
      } else {
        dds <- data.frame(DESeq2::results(dds)[, c("stat", "pvalue")])
      }
      return(dds)
      # z <- (qnorm(dds$pvalue/2, lower.tail = FALSE) * sign(dds$stat))
      # names(z) <- rownames(x)
      # return(z)
      ## t <- qnorm(pt(abs(t), largo + largoy - 2, lower.tail = FALSE), lower.tail=FALSE)*sign(t)
      ## names(t) <- rownames(x)
      ## return(t)
    }, x=x, y=y, repos=repos, mc.cores=cores)
    res <- lapply(res, function(x) x)
  }
  else {
    if (verbose) pb <- txtProgressBar(max=per, style=3)
    res <- lapply(1:per, function(i, x, y, repos, pb, verbose) {
      if (verbose) setTxtProgressBar(pb, i)
      expset <- cbind(x, y)
      repeat{
        sorder <- sample(ncol(expset), replace=repos)
        if (length(unique(sorder[1:ncol(x)]))>1 & length(unique(sorder[-(1:ncol(x))]))>1) break
        if (verbose) message("-", appendLF=FALSE)
      }
      x1 <- filterColMatrix(expset, sorder[1:ncol(x)])
      y1 <- filterColMatrix(expset, sorder[-(1:ncol(x))])
      ## largo <- rowSums(!is.na(x1))
      ## largoy <- rowSums(!is.na(y1))
      ## t <- ((rowMeans(x1, na.rm=TRUE) - rowMeans(y1, na.rm=TRUE))/sqrt(((largo - 1) *  rowVars(x1) + (largoy - 1) * rowVars(y1))/(largo + largoy - 2))/sqrt(1/largo + 1/largoy))[, 1]
      countData <- cbind(x1, y1)
      colData <- colData[colnames(countData), ]
      dds <- DESeq2::DESeqDataSetFromMatrix(countData, colData, design)
      if (!is.null(refcol) & !is.null(reflvl)) {
        colData(dds)[, refcol] <- relevel(colData(dds)[, refcol], reflvl)
      }
      dds <- DESeq2::DESeq(dds)
      if (!is.null(resname)) {
        dds <- data.frame(DESeq2::results(dds, name = resname)[, c("stat", "pvalue")])
      } else if (!is.null(contrast)) {
        dds <- data.frame(DESeq2::results(dds, contrast = contrast)[, c("stat", "pvalue")])
      } else {
        dds <- data.frame(DESeq2::results(dds)[, c("stat", "pvalue")])
      }
      return(dds)
      # z <- (qnorm(dds$pvalue/2, lower.tail = FALSE) * sign(dds$stat))
      # names(z) <- rownames(x)
      # return(z)
      ## t <- qnorm(pt(abs(t), largo + largoy - 2, lower.tail = FALSE), lower.tail=FALSE)*sign(t)
      ## names(t) <- rownames(x)
      ## return(t)
    }, x=x, y=y, repos=repos, pb=pb, verbose=verbose)
  }
  ## colnames(res) <- 1:per
  stat <- lapply(res, "[", "stat")
  stat <- do.call(cbind, stat)
  stat <- data.matrix(stat)
  colnames(stat) <- 1:per
  pvalue <- lapply(res, "[", "pvalue")
  pvalue <- do.call(cbind, pvalue)
  pvalue <- data.matrix(pvalue)
  colnames(pvalue) <- 1:per
  list <- list("stat" = stat, "pvalue" = pvalue)
  if (verbose) message("\n", date())
  return(list)
  # return(res)
}

#######
#' Plot msviper results
#'
#' This function generate a plot for msviper results showing the enrichment of the target genes for each significant master regulator on the gene expression signature
#'
#' @param x msviper object produced by \code{msviper} function
#' @param mrs Either an integer indicating the number of master regulators to include in the plot, or a character vector containing the names of the master regulators to include in the plot
#' @param color Vector of two components indicating the colors for the negative and positive parts of the regulon
#' @param FDR Optional matrix of FDR to include in the plot
#' @param bins Number of bins to split the vector of scores in order to compute the density color of the bars
#' @param cex Number indicating the text size scaling, 0 indicates automatic scaling
#' @param density Integrer indicating the number of steps for the kernel density. Zero for not ploting it
#' @param smooth Number indicating the proportion of point for smoothing the density distribution. Zero for not using the smoother
#' @param sep Number indicating the separation from figure and text
#' @param hybrid Logical, whether the 3-tail approach used for computingthe enrichment should be reflected in the plot
#' @param include Vector indicating the information to include as heatmap to the right of the msviper plot: expression and activity
#' @param gama Positive number indicating the exponential transformation for the activity and expression color scale
#' @param ... Given for compatibility to the plot generic function
#' @return Nothing, a plot is generated in the default output device
#' @method plot_msviper
#' @export
#' @seealso \code{\link{msviper}}

plot_msviper <- function(x, mrs=10, color=c("cornflowerblue","salmon"), FDR=NULL, bins=500, cex=0, density=0, smooth=0, sep=.2, hybrid=TRUE, include=c("expression", "activity"), gama=2, ...) {
  maobject <- x
  rm(x)
  marg <- par("mai")
  rlist <- maobject$signature
  if (ncol(rlist)>0) rlist <- rowMeans(rlist)
  if (length(mrs)==1 & is.numeric(mrs[1])) {
    mrs <- names(maobject$es$nes)[order(maobject$es$p.value)[1:round(mrs)]]
    mrs <- mrs[order(maobject$es$nes[match(mrs, names(maobject$es$nes))], decreasing=TRUE)]
  }
  groups <- maobject$regulon[match(mrs, names(maobject$regulon))]
  if (is.null(FDR)) FDR <- p.adjust(maobject$es$p.value, "fdr")[match(mrs, names(maobject$es$nes))]
  if (is.data.frame(FDR)) FDR <- as.matrix(FDR)
  if (is.matrix(rlist)) rlist <- rlist[,1]
  if (min(rlist)<0) rlist <- sort(rlist)
  else rlist <- sort(-rlist)
  groups <- groups[length(groups):1]
  color1 <- color
  layout(matrix(1:2, 1, 2), widths=c(10-length(include), length(include)))
  color <- rgb2hsv(col2rgb(color))
  satval <- color[3,]
  color <- color[1,]
  preg <- as.numeric(any(sapply(groups, function(x) any(x$tfmode<0))))+1
  textsize <- 1
  xlimit <- c(0, length(rlist)*(1+.2*max(nchar(names(groups)))/8))
  if (!is.null(FDR)) {
    if (is.matrix(FDR)) {
      FDR <- FDR[nrow(FDR):1, ]
      FDR <- FDR[, ncol(FDR):1]
      xlimit <- c(-length(rlist)*(sep*ncol(FDR)+.02), length(rlist)*(1+.2*max(nchar(names(groups)))/9))
    }
    else {
      FDR <- FDR[length(FDR):1]
      xlimit <- c(-length(rlist)*.12, length(rlist)*(1+.2*max(nchar(names(groups)))/9))
      textsize=.8
    }
  }
  if (length(include)>0 & include[1] != "") par(mai=c(marg[1:3], .05))
  plot(1, type="n", ylim=c(0, length(groups)), xlim=xlimit, axes=FALSE, ylab="", xlab="", yaxs="i")
  if (cex>0) textsize <- (length(groups)<=20)*cex+(length(groups)>20)*(20/length(groups))*cex
  switch(preg,
         {
           for (i in 1:length(groups)) {
             densi <- rep(0, length(rlist))
             x <- which(names(rlist) %in% names(groups[[i]]$tfmode))
             if (length(x)>0) {
               densi[x] <- 1
               denStep <- round(length(densi)/bins)
               x1 <- x[x<denStep]
               x2 <- x[x>=denStep & x <= (length(rlist)-denStep)]
               x3 <- x[x>(length(rlist)-denStep)]
               densiRes <- sapply(x2, function(i, densi, denStep)
                 sum(densi[(i-denStep):(i+denStep)]), densi=densi, denStep=denStep)
               densiRes <- densiRes/max(densiRes)
               temp <- rlist[x]
               if (satval[1]==0) temp <- hsv((temp<0)*color[1] + (temp>0)*color[2], satval, 1-densiRes)
               else temp <- hsv((sign(temp)<0)*color[1] + (sign(temp)>0)*color[2], densiRes, satval)
               for (ii in order(densiRes)) lines(c(x[ii], x[ii]),c(i-1, i), col=temp[ii])
               if (density>0) {
                 denStep <- round(length(densi)/density)
                 xpos <- seq(denStep, length(rlist)-denStep, length=density)
                 densiRes <- sapply(xpos, function(i, densi, denStep) {
                   sum(densi[(i-denStep):(i+denStep)])
                 }, densi=densi, denStep=denStep)
                 densiRes <- densiRes/max(densiRes)
                 if (smooth>0) densiRes <- smooth.spline(xpos, densiRes, spar=smooth)$y
                 lines(xpos, i+densiRes-1)
               }
             }
           }
           text(rep(length(rlist)*1.02, length(groups)), 1:length(groups)-.5, names(groups), adj=0, cex=textsize)
           if (!is.null(FDR)) {
             if (is.matrix(FDR)) {
               for (i in 1:ncol(FDR)) {
                 text(rep(-length(rlist)*(sep*(i-1)+.02), length(groups)), 1:length(groups)-.5, FDR[, i], adj=1, cex=.85*textsize)
                 text(-length(rlist)*(sep*(i-1)+.02), length(groups)+.5, colnames(FDR)[i], adj=1, cex=1)
               }
             }
             else {
               text(rep(-length(rlist)*.02, length(groups)), 1:length(groups)-.5, signif(FDR, 3), adj=1, cex=.85*textsize)
               text(0, length(groups)+.5, ifelse(max(FDR)>1, "oddsR", "FDR"), adj=1, cex=1.2)
             }
           }
           text(length(rlist)*1.02, length(groups)+.5, "MR", adj=0, cex=1.2)
         },
         {
           for (i in 1:length(groups)) {
             for (ii in 1:2) {
               tset <- groups[[i]]$tfmode
               tset <- tset[(tset<0 & ii==1)|(!(tset<0) & ii==2)]
               tset1 <- names(tset)[abs(tset)>.5]
               tset2 <- names(tset)[abs(tset)<.5]
               if (length(tset)>1) {
                 densi <- rep(0, length(rlist))
                 x <- match(names(tset), names(rlist))
                 tw1 <- rep(1, length(x))
                 if (hybrid) {
                   x <- match(tset1, names(rlist))
                   if (ii==1) {
                     x <- c(x, match(tset2, names(sort(-abs(rlist)*sign(maobject$es$nes[names(maobject$es$nes) == names(groups)[i]])))))
                   }
                   else {
                     x <- c(x, match(tset2, names(sort(abs(rlist)*sign(maobject$es$nes[names(maobject$es$nes) == names(groups)[i]])))))
                   }
                   tw1 <- groups[[i]]$likelihood[match(c(tset1, tset2), names(groups[[i]]$tfmode))]
                   tw1 <- tw1/max(tw1)
                 }
                 densi[x] <- 1
                 denStep <- round(length(densi)/bins)
                 x1 <- x[x<denStep]
                 x2 <- x[x>=denStep & x <= (length(rlist)-denStep)]
                 x3 <- x[x>(length(rlist)-denStep)]
                 densiRes <- sapply(x2, function(i, densi, denStep) {
                   sum(densi[(i-denStep):(i+denStep)])
                 }, densi=densi, denStep=denStep)
                 densiRes <- densiRes*(tw1[x>=denStep & x <= (length(rlist)-denStep)])
                 densiRes <- densiRes/max(densiRes)
                 temp <- rlist[x]
                 temp <- hsv(color[ii], densiRes, satval)
                 for (iii in order(densiRes)) {
                   lines(c(x[iii], x[iii]),c(i-1+(ii-1)/2, i-1+ii/2), col=temp[iii])
                 }
                 if (density>0) {
                   denStep <- round(length(densi)/density)
                   xpos <- seq(denStep, length(rlist)-denStep, length=density)
                   densiRes <- sapply(xpos, function(i, densi, denStep)
                     sum(densi[(i-denStep):(i+denStep)]), densi=densi, denStep=denStep)
                   densiRes <- densiRes/max(densiRes)
                   if (smooth>0) densiRes <- smooth.spline(xpos, densiRes, spar=smooth)$y
                   lines(xpos, i-1+densiRes/2+(ii-1)/2)
                 }
               }
             }
           }
           text(rep(length(rlist)*1.02, length(groups)), 1:length(groups)-.5, names(groups), adj=0, cex=textsize)
           if (!is.null(FDR)) {
             if (is.matrix(FDR)) {
               for (i in 1:ncol(FDR)) {
                 text(rep(-length(rlist)*(sep*(i-1)+.02), length(groups)), 1:length(groups)-.5, FDR[, i], adj=1, cex=.85*textsize)
                 text(-length(rlist)*(sep*(i-1)+.02), length(groups)+.5, colnames(FDR)[i], adj=1, cex=1)
               }
             }
             else {
               text(rep(-length(rlist)*.02, length(groups)), 1:length(groups)-.5, signif(FDR, 3), adj=1, cex=.85*textsize)
               axis(3, -.05*length(rlist), ifelse(max(FDR)>1, "oddsR", "FDR"), adj=1, cex=1.2, tick=FALSE, line=-.5)
             }
           }
           axis(3, length(rlist)*1.05, "MR", adj=0, cex=1.2, line=-.5, tick=FALSE)
         })
  abline(h=0:length(groups))
  lines(c(0, 0), c(0, length(groups)))
  lines(c(length(rlist), length(rlist)), c( 0, length(groups)))
  ss <- maobject$signature
  if (!is.null(dim(ss))) ss <- ss[, 1]
  x <- NULL
  xn <- NULL
  xpos <- NULL
  if("activity" %in% include) {
    x <- maobject$es$nes[match(mrs, names(maobject$es$nes))]/max(abs(maobject$es$nes))
    xn <- "NES"
  }
  if ("expression" %in% include) {
    x <- cbind(x, ss[match(mrs, names(ss))]/max(abs(ss)))
    xn <- c(xn, "EXP")
    xpos <- rank(-abs(ss))[match(mrs, names(ss))]
  }
  if (!is.null(x)) {
    if(is.null(dim(x))) dim(x) <- c(length(x), 1)
    rownames(x) <- colnames(x) <- NULL
    marg[2] <- .05
    # marg[4] <- .5
    par(mai=marg)
    scmax <- max(abs(x), na.rm=TRUE)
    x <- abs(x/scmax)^gama*sign(x)
    x <- filterRowMatrix(x, nrow(x):1)
    x1 <- x
    x1[is.na(x1)] <- 0
    coli <- hsv(ifelse(x1<0, color[1], color[2]), abs(x1), 1)
    coli[is.na(x)] <- hsv(0, 0, .5)
    image(1:ncol(x), 1:nrow(x), t(matrix(1:(ncol(x)*nrow(x)), nrow(x), ncol(x))), col=coli, ylab="", xlab="", axes=FALSE, yaxs="i")
    box()
    grid(ncol(x), nrow(x), col="black", lty=1)
    axis(3, 1:length(xn), xn, las=1, line=-.5, tick=FALSE)
    # axis(4, length(xpos):1, xpos, las=1, tick=FALSE, line=-.8, cex.axis=.85*textsize)
  }    
  par(mai=marg)
}