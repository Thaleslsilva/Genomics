#exemplo de uso
#pro.man<-data.frame(CHR=map2$Chromosome,BP=map2$Position,P=var.expl)
#manhattan(pro.man)     

manhattan<-function(dataframe, colors=c("darkblue","darkred","darkgreen","gray15"),
                    cex.x.axis=1, limitchromosomes=1:29,annotate=NULL,
                    ...)
{
  d=dataframe
  genomewideline<-5*sd(d$P)
  suggestiveline<-4*sd(d$P)
  ymax<-max(d$P)
  ymin<-min(abs(d$P))
  #throws error if you don't have columns named CHR, BP, and P in your data frame.
  if (!("CHR" %in% names(d) & "BP" %in% names(d) & "P" %in% names(d))) 
    stop("Make sure your data frame contains columns CHR, BP, and P")
  
  # limits chromosomes to plot
  if (any(limitchromosomes)) d=d[d$CHR %in% limitchromosomes, ]
  
  # remove na's, sort by CHR and BP
  d=subset(na.omit(d[order(d$CHR, d$BP), ]))
  
  # abs(sol)
  d$logp = abs(d$P) 
  
  # sets colors based on colors argument.
  colors <- rep(colors,max(d$CHR))[1:max(d$CHR)]
  
  # creates continuous position markers for x axis for entire chromosome. also creates tick points.
  d$pos=NA
  ticks=NULL
  lastbase=0
  numchroms=length(unique(d$CHR))
  if (numchroms==1) {
    d$pos=d$BP
    ticks=floor(length(d$pos))/2+1
  } else {
    for (i in unique(d$CHR)) {
      if (i==1) {
        d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP
      } else {
        lastbase=lastbase+tail(subset(d,CHR==i-1)$BP, 1)
        d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP+lastbase
      }
      ticks=c(ticks, d[d$CHR==i, ]$pos[floor(length(d[d$CHR==i, ]$pos)/2)+1])
    }
  }
  
  # create the plot
  if (numchroms==1) {
    # if you only have a single chromosome, the x axis is the chromosomal position
    with(d, plot(pos, logp, ylim=c(ymin,ymax), ylab="abs(SNP effect)", 
                 xlab=paste("Chromosome",unique(d$CHR),"position"), ...))
  } else {
    # if you have multiple chromosomes, first make the plot with no x-axis (xaxt="n")
    with(d, plot(pos, logp, ylim=c(ymin,ymax), ylab="abs(SNP effect)", 
                 xlab="Chromosome", xaxt="n", type="n", ...))
    # then make an axis that has chromosome number instead of position
    axis(1, at=ticks, lab=unique(d$CHR), cex.axis=cex.x.axis)
    icol=1
    for (i in unique(d$CHR)) {
      with(d[d$CHR==i, ],points(pos, logp, col=colors[icol], ...))
      icol=icol+1
    }
  }
  
  # create a new data frame with rows from the original data frame where SNP is in annotate character vector.
  # then plot those points over the original graph, but with a larger point size and a different color.
  if (!is.null(annotate)) {
    d.annotate=d[which(d$SNP %in% annotate), ]
    with(d.annotate, points(pos, logp, col="red", cex=2.5, ...)) 
  }
  
  # add threshold lines
  #if (suggestiveline) abline(h=suggestiveline, col="blue")
  #if (genomewideline) abline(h=genomewideline, col="red")
}
