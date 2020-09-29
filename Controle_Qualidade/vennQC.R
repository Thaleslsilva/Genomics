circle <- function(x, y, r, ...) {
 ang <- seq(0, 2*pi, length = 100)
 xx <- x + r * cos(ang)
 yy <- y + r * sin(ang)
 polygon(xx, yy, ...)
}

venndia <- function(nA, nB, nC,nAB,nAC,nBC,nABC,nD=NULL,getdata=FALSE, ...){
 
 par(mar=c(2, 2, 0, 0))
 plot(-10, -10, ylim=c(0, 12), xlim=c(0, 12), axes=FALSE, ...)
 circle(x=3, y=6, r=3, col=rgb(1,0,0,.5), border=NA)
 circle(x=6, y=6, r=3, col=rgb(0,.5,.1,.5), border=NA)
 circle(x=4.5, y=3, r=3, col=rgb(0,0,1,.5), border=NA)
 if (!is.null(nD)){circle(x=9,y=3,r=1,col=rgb(0.5,0.5,0.5,.5),border=NA)
    text( x=c(1.2, 7.7, 4.5,9), y=c(7.8, 7.8, 0.8,1.5), 
          c("MAF", "HWE", "Call rate","r2"), cex=3, col="gray30" )} else {
 text( x=c(1.2, 7.7, 4.5), y=c(7.8, 7.8, 0.8), c("MAF", "HWE", "Call rate"), cex=3, col="gray30" )}
 
 text(
  x=c(2, 7, 4.5, 4.5, 3, 6, 4.5,9), 
  y=c(7, 7, 2  , 7  , 4, 4, 5,2.5), 
  c(nA, nB, nC, nAB, nAC, nBC, nABC,nD), 
  cex=2
 )
 
 if(getdata){
  list(A=uniqueA, B=uniqueB, C=uniqueC, 
       AB=intersAB , AC=intersAC , BC=intersBC , 
       ABC=intersABC
  )
 }
}