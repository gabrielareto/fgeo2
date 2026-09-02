# <function>
# <name>
# skewness
# </name>
# <description>
# Sample skewness. The biased portion is the population skewness; correction is for finite sample.  
# D. N. Joanes and C. A. Gill. “Comparing Measures of Sample Skewness and Kurtosis”. The Statistician 47(1):183–189
# </description>
# <arguments>
#
# </arguments>
# <sample>
# 
# </sample>
# <source>
skewness=function(x)
{
 x=x[!is.na(x)]
 n=length(x)
 if(length(x)<=2) return(0)

 mn=mean(x)
 sumcube=sum((x-mn)^3)
 sumsq=sum((x-mn)^2)
 biased=sqrt(n)*sumcube/(sumsq^1.5)
 correction=sqrt(n)*sqrt((n-1)/(n-2))

 return(biased*correction)
}
# </source>
# </function>
# 
#
#
# <function>
# <name>
# skewness
# </name>
# <description>
# Standard error of skewness. Depends only on sample size. 
# </description>
# <arguments>
#
# </arguments>
# <sample>
# 
# </sample>
# <source>
se.skewness=function(x)
{
 x=x[!is.na(x)]
 n=length(x)
 if(length(x)<=2) return(0)

 part1=6/(n-2)
 part2=n/(n+1)
 part3=(n-1)/(n+3)
 return(sqrt(part1*part2*part3))
}
# </source>
# </function>
# 
# 
# 
# <function>
# <name>
# skewness
# </name>
# <description>
# Sample kurtosis. The biased portion is the population kurtosis; corrected is for finite sample.  
# </description>
# <arguments>
#
# </arguments>
# <sample>
# 
# </sample>
# <source>
kurtosis=function(x)
{
 x=x[!is.na(x)]
 n=length(x)
 if(length(x)<=3) return(0)

 mn=mean(x)
 moment4=sum((x-mn)^4)/n
 moment2=sum((x-mn)^2)/n
 biased=moment4/(moment2^2)
 
 part1=(biased*(n+1)+6)/(n-2)
 part2=(n-1)/(n-3)

 return(part1*part2)
}
# </source>
# </function>
# 
# 
# <function>
# <name>
# skewness
# </name>
# <description>
# Standard error of kurtosis. Depends only on sample size. 
# </description>
# <arguments>
#
# </arguments>
# <sample>
# 
# </sample>
# <source>
se.kurtosis=function(x)
{
 x=x[!is.na(x)]
 n=length(x)
 if(length(x)<=3) return(0)

 SES=se.skewness(x)
 part=n/(n-3)
 part1=n*part-1/(n-3)
 part2=n+5
 
 return(2*SES*sqrt(part1/part2))
}
# </source>
# </function>
# 
# 
# <function>
# <name>
# regslope
# </name>
# <description>
# Returns slope of regression as single scalar (for use with apply).
# </description>
# <arguments>
# 
# </arguments>
# <sample>
# 
# </sample>
# <source>
regslope=function(y,x)
{
 if(length(y[is.infinite(y)])>0) return(NA)
 if(length(y[is.infinite(x)])>0) return(NA)

 fit=lm(y~x,na.action="na.exclude")
 return(summary(fit)$coef[2,1])
}



# </source>
# </function>
# 
# 
# 
# <function>
# <name>
# regslope.noint
# </name>
# <description>
#  Returns slope of regression with no intercept as single scalar (for use with apply).

# </description>
# <arguments>
# 
# </arguments>
# <sample>
# 
# </sample>
# <source>
regslope.noint=function(y,x)
{
 fit=lm(y~x+0)
 return(summary(fit)$coef[2,1])
}



# </source>
# </function>
# 
# 
# 
# <function>
# <name>
# regress.plot
# </name>
# <description>
#  Performs regression in convenient way and returns coefficients and
# probabilities in a single vector, and plots a graph.

# </description>
# <arguments>
# 
# </arguments>
# <sample>
# 
# </sample>
# <source>
regress.plot=function(x,y,title="",graphit=T,add=F,pts=T,clr="blue",ptnames=NULL,xrange=NULL,yrange=NULL,xtitle=NULL,ytitle=NULL)
{
 fit=lm(y~x)
 regcoef=summary(fit)$coef[,1]
 prob=summary(fit)$coef[,4]
 rsq=cor(x,y)^2

 if(is.null(xrange)) xrange=range(x)
 if(is.null(yrange)) yrange=range(y)
 if(is.null(xtitle)) xtitle="x"
 if(is.null(ytitle)) ytitle="y"
 if(graphit)
  {
   if(add & pts) points(x,y,pch=16)
   if(!add & pts) plot(x,y,pch=16,main=title)
   abline(fit,col=clr)

   if(!is.null(ptnames)) identify(x,y,ptnames)
  }

 return(list(coef=c(regcoef,prob,rsq),full=summary(fit)))
}



# </source>
# </function>
# 
# 
# 
# <function>
# <name>
# regress.loglog
# </name>
# <description>
#  Performs regression and graphs in a convenient way: with or without log-transforming x and y variables (the option addone
# can be included to handle zeros for log-transformation), with or
# without manual point labelling, without or without the best-fit line added, and with many options for colors and points. 
# add can be a vector of length 2, a constant to be added to every value
# of x, y to remove zeroes.

# </description>
# <arguments>
# 
# </arguments>
# <sample>
# 
# </sample>
# <source>
regress.loglog=function(x,y,xlog=TRUE,ylog=TRUE,addone=NULL,graphit=TRUE,xrange=NULL,yrange=NULL,add=FALSE,pts=16,lwidth=1,
                        drawline="solid",title="",xtitle=NULL,ytitle=NULL,ptnames=NULL,ptsize=1,clr="blue",lineclr="red",includeaxis=TRUE)
{
 exist = !is.na(x) & !is.na(y)
 x=x[exist]
 y=y[exist]
 ptnames=ptnames[exist]

 if(!is.null(addone)) 
  {
   x[x==0]=x[x==0]+addone[1]
   y[y==0]=y[y==0]+addone[2]
  }

 pos = !is.na(x)
 if(xlog) pos = pos & x>0
 if(ylog) pos = pos & y>0
 x=x[pos]
 if(length(x)==0)  return(list(coef=NULL,prob=NULL,rsq=NULL,pred=NULL,full=NULL))  ## Added March 2010

 y=y[pos]
 ptnames=ptnames[pos]

 if(xlog) xreg=log(x)
 else xreg=x
 if(ylog) yreg=log(y)
 else yreg=y

 fit=lm(yreg~xreg)
 
 regcoef=summary(fit)$coef[,1]
 prob=summary(fit)$coef[,4]
 rsq=cor(xreg,yreg)^2

 if(is.null(xrange)) xrange=range(x)
 if(is.null(yrange)) yrange=range(y)
 if(is.null(xtitle)) xtitle="x"
 if(is.null(ytitle)) ytitle="y"

 if(ylog) predy=exp(fit$fitted)
 else predy=fit$fitted

 if(graphit)
  {
   logaxs=""
   if(xlog) logaxs=pst(logaxs,"x")
   if(ylog) logaxs=pst(logaxs,"y")

   if(add & !is.null(pts)) points(x,y,pch=pts,col=clr,cex=ptsize,cex.lab=ptsize,cex.axis=ptsize)

   if(!add & !is.null(pts))
      plot(x,y,pch=pts,main=title,log=logaxs,xlim=xrange,ylim=yrange,xlab=xtitle,ylab=ytitle,col=clr,
	       cex=ptsize,cex.lab=ptsize,cex.axis=ptsize,axes=includeaxis)
   if(!includeaxis) box()

   ord=order(x)
   if(!is.null(drawline)) lines(x[ord],predy[ord],col=lineclr,lty=drawline,lwd=lwidth)

   if(!is.null(ptnames)) identify(x,y,ptnames)
  }

 return(list(coef=regcoef,prob=prob,rsq=rsq,pred=predy,full=summary(fit)))
}



# </source>
# </function>
# 
# 
# 
# <function>
# <name>
# majoraxisreg
# </name>
# <description>
#  A major axis regression with parameters fitted by optim. The regression
# is the line which minimizes perpendicular distance summed over all points
# (and squared).

# </description>
# <arguments>
# 
# </arguments>
# <sample>
# 
# </sample>
# <source>
majoraxisreg=function(x,y,title="",graphit=F,add=F,pts=T,clr="blue",xtitle="x",ytitle="y",ptsize=1,labsize=1)
{
 start.param=c(1,1)

 fit=optim(start.param,minum.perpdist,x=x,y=y)
 m=fit$par[2]
 b=fit$par[1]

 if(graphit)
  {
   if(add & pts) points(x,y,pch=16)
   if(!add & pts) plot(x,y,pch=16,main=title,xlab=xtitle,ylab=ytitle,cex=ptsize,cex.lab=labsize,cex.axis=labsize)
   abline(b,m,col=clr)
  }

 return(fit$par)
}


# </source>
# </function>
# 
# 
# 
# <function>
# <name>
# minum.perpdist
# </name>
# <description>
#  The sum of squares used by majoraxisreg.

# </description>
# <arguments>
# 
# </arguments>
# <sample>
# 
# </sample>
# <source>
minum.perpdist=function(param,x,y)
    return(sumsq(perpendicular.distance(param[1],param[2],x,y)))



# </source>
# </function>
# 
# 
# 
# <function>
# <name>
# majoraxisreg.no.int
# </name>
# <description>
#  Major axis regression with no intercept. Only a slope
# is returned. Below is the same for standard regression.

# </description>
# <arguments>
# 
# </arguments>
# <sample>
# 
# </sample>
# <source>
majoraxisreg.no.int=function(x,y)
{
 inc=(!is.na(x)&!is.na(y))
 x=x[inc]
 y=y[inc]

 a=sum(x*y)
 b=sum(x^2)-sum(y^2)
 c=(-a)

 answer=numeric()
 answer[1]=(-b+sqrt(b^2-4*a*c))/(2*a)
 answer[2]=(-b-sqrt(b^2-4*a*c))/(2*a)
 return(max(answer,na.rm=T))
}
# </source>
# </function>
# 
# 
# <function>
# <name>
# standardreg.no.int
# </name>
# <description>
#  Standard regression with no intercept.

# </description>
# <arguments>
# 
# </arguments>
# <sample>
# 
# </sample>
# <source>
standardreg.no.int=function(x,y)
{
 inc=!is.na(x)&!is.na(y)
 x=x[inc]
 y=y[inc]

 return(sum(x*y)/sum(x^2))
}
# </source>
# </function>
# 
# 
# 
# <function>
# <name>
# autoregression
# </name>
# <description>
# Autocorrelation with a given lag of a vector y.
# </description>
# <arguments>
# 
# </arguments>
# <sample>
# 
# </sample>
# <source>
autoregression=function(y,lag,xlog=FALSE,ylog=FALSE,graphit=TRUE)
{
 N=length(y)

 lagindex=(1+lag):N
 nonlagindex=1:(N-lag)

 ylag=y[lagindex]
 ynon=y[nonlagindex]

 return(regress.loglog(ynon,ylag,xlog=xlog,ylog=ylog,graphit=graphit))
}
# </source>
# </function>

# <function>
# <name>
# bootstrap.corr
# </name>
# <description>
#  Running bootstrap on a correlation. Any columsn can be chosen from the submitted dataset, by number or name.

# </description>
# <arguments>
# 
# </arguments>
# <sample>
# 
# </sample>
# <source>
bootstrap.corr=function(dataset,xcol=(-1),ycol=(-1),xcolname="no",ycolname="no",boot=100)
{
 norecords=dim(dataset)[1]

 if(xcol<0) xcol=which(colnames(dataset)==xcolname)
 else xcolname=colnames(dataset)[xcol]
 if(ycol<0) ycol=which(colnames(dataset)==ycolname)
 else ycolname=colnames(dataset)[ycol]

 x=dataset[,xcol]
 y=dataset[,ycol]
 plot(x,y,pch=16,xlab=xcolname,ylab=ycolname)

 inter=corrcoef=rsq=prob=numeric()

 for(i in 1:boot)
  {
   s=sample(1:norecords,norecords,replace=T)

   fit=lm(y[s]~x[s])

   inter[i]=summary(fit)$coef[1,1]
   corrcoef[i]=summary(fit)$coef[2,1]
   prob[i]=summary(fit)$coef[2,4]
   rsq[i]=summary(fit)$r.squared

   if(i<=100) abline(fit)
  }

 if(boot<200)
  {
   interci=bootconf(inter)
   corrci=bootconf(corrcoef)
   probci=bootconf(prob)
   rsqci=bootconf(rsq)
  }
 else
  {
   interci=quantile(inter,prob=c(.025,.975))
   corrci=quantile(corrcoef,prob=c(.025,.975))
   probci=quantile(prob,prob=c(.025,.975))
   rsqci=quantile(rsq,prob=c(.025,.975))
  }

 return(list(inter=interci,corr=corrci,prob=probci,rsq=rsqci))
}



# </source>
# </function>
# 
# 
# 
# <function>
# <name>
# bootconf
# </name>
# <description>
#  A simple calculation of confidence limits based on the SD of a vector.

# </description>
# <arguments>
# 
# </arguments>
# <sample>
# 
# </sample>
# <source>
bootconf=function(x)
{
 lower=mean(x)-1.96*sd(x)
 upper=mean(x)+1.96*sd(x)

 return(c(lower,upper))
}


# </source>
# </function>

# <function>
# <name>
# CI
# </name>
# <description>
# Confidence limits (quantiles) from a vector at specified probabilities. Default is 95% confidence interval. 

# </description>
# <arguments>
# 
# </arguments>
# <sample>
# 
# </sample>
# <source>
CI=function(x,prob=c(.025,.975),na.rm=FALSE) 
     return(quantile(x,prob=prob,na.rm=na.rm))



# </source>
# </function>
# 
# 
# 
# <function>
# <name>
# hist.compare
# </name>
# <description>
# Compares two histograms with a Kolmogorov approach.
# </description>
# <arguments>
# 
# </arguments>
# <sample>
# 
# </sample>
# <source>
hist.compare=function(x,y,div,breaks=NULL)
{
 y=y[!is.na(y)]
 x=x[!is.na(x)]
 obsrange=range(x)
 
 if(is.null(breaks)) breaks=seq(obsrange[1],obsrange[2]+div,by=div)

 xcat=cut(x,breaks=breaks,right=FALSE)
 # xcount=table(xcat)
 ycat=cut(y,breaks=breaks,right=FALSE)
 ycount=table(ycat)
 # value=tapply(y,ycat,mean)

 pred=ycount/sum(ycount)
 pred[pred==0]=1/(2*sum(ycount))
 pred=pred/sum(pred)

 m=match(xcat,names(pred))
 llike=log(pred[m])
  browser()

 return(sum(llike))
}



# </source>
# </function>
# 
# 
# 
# <function>
# <name>
# harmonic.mean
# </name>
# <description>
#  Harmonic mean of a vector x. NAs and nonzero values can be ignored, and a constant can be added to every x.

# </description>
# <arguments>
# 
# </arguments>
# <sample>
# 
# </sample>
# <source>
harmonic.mean=function(x,add=0,na.rm=TRUE)
{
 x=x+add

 missing=which(x<=0 | is.na(x))
 if(length(missing)>0)
  {
   if(!na.rm) return(NA)
   x=x[-missing]
  }

 logx=log(x)
 return(exp(mean(logx)))
}



# </source>
# </function>
# 
# 
# 
# <function>
# <name>
# cumul.above
# </name>
# <description>
#  Given y values as a function of x, this seeks the x at which the curve passes through a given y. It sets
# a variable whichabove to 0 for all cases where y>cutoff, otherwise 0, then fits a logistic regression.
# The midpoint of the logistic regression is a good estimate. 

# </description>
# <arguments>
# 
# </arguments>
# <sample>
# 
# </sample>
# <source>
cumul.above=function(x,y,cutoff,logaxs="",graphit=TRUE,returnfull=FALSE)
{
 whichabove=y>cutoff
 if(logaxs=="x" | logaxs=="xy") runx=log(x)
 else runx=x

 fit=glm(whichabove~runx,family=binomial)

 if(logaxs=="x" | logaxs=="xy") mid=exp(-fit$coef[1]/fit$coef[2])
 else mid=(-fit$coef[1]/fit$coef[2])

 xord=x[order(x)]
 yord=y[order(x)]
 if(graphit)
  {
   plot(xord,yord,pch=16,log=logaxs)
   abline(v=mid)
   abline(h=cutoff)
  }

 return(mid)
}



# </source>
# </function>
# 
# 
# 
# <function>
# <name>
# sumsq
# </name>
# <description>
#  A trivial function used in minimizing sums of squares.

# </description>
# <arguments>
# 
# </arguments>
# <sample>
# 
# </sample>
# <source>
sumsq=function(x) 
  return(sum(x^2,na.rm=T))



# </source>
# </function>
# 
# 
# 
# <function>
# <name>
# is.odd
# </name>
# <description>
#  A trivial function to test whether numbers (scalar or vector) are odd. 

# </description>
# <arguments>
# 
# </arguments>
# <sample>
# 
# </sample>
# <source>
is.odd=function(x)
{
 answer=rep(FALSE,length(x))
 y=which(x%%2!=0)
 answer[y]=TRUE
 return(answer)
}


# </source>
# </function>
# 
# 
# 
# <function>
# <name>
# border.distance
# </name>
# <description>
#  Returns distance from a point to the nearest boundary of a rectangle (plot). Accepts either separate
# x-y coordinates, or an object where x is first column, y is second. The lower left corner of the plot is
# assumed to be 0,0. 

# </description>
# <arguments>
# 
# </arguments>
# <sample>
# 
# </sample>
# <source>
border.distance=function(x,y,plotdim=c(1000,500),pt=NULL)
{
 if(!is.null(pt)) 
  {
   x=pt[,1]
   y=pt[,2]
  }
  
 missingx = is.na(x) | x<0 | x>=plotdim[1]
 missingy = is.na(y) | y<0 | y>=plotdim[2]

 closestx=closesty=closest=numeric()

 lefthalf = x<plotdim[1]/2 & !is.na(x)
 bottomhalf = y<plotdim[2]/2  & !is.na(y)

 closestx[lefthalf]=x[lefthalf]
 closestx[!lefthalf]=plotdim[1]-x[!lefthalf]
 closesty[bottomhalf]=y[bottomhalf]
 closesty[!bottomhalf]=plotdim[2]-y[!bottomhalf]

 xcloser=closestx<closesty

 closest[xcloser]=closestx[xcloser]
 closest[!xcloser]=closesty[!xcloser]
 closest[missingx | missingy]=NA

 return(closest)
}



# </source>
# </function>
# 
# 
# 
# <function>
# <name>
# regsum
# </name>
# <description>
#  This carries out either first or second order polynomial regression,
# finds the x- and y-values at y's peak if its second order,
# otherwise the x-intercept.

# </description>
# <arguments>
# 
# </arguments>
# <sample>
# 
# </sample>
# <source>
regsum=function(x,y,poly=1,graphit=F,yrange=c(-1,-1),yaxslab="multiplier",
                newgraph=F,add=F,pts=T,ptype=16,ltype="solid")
{
 x2=x^2
 x3=x^3
 if(poly==1) fit=lm(y~x)
 else if(poly==2) fit=lm(y~x+x2)
 else if(poly==3) fit=lm(y~x+x2+x3)

 a=summary(fit)$coef[1,1]
 b=summary(fit)$coef[2,1]
 c=d=0
 if(poly>1) c=summary(fit)$coef[3,1]
 if(poly>2) d=summary(fit)$coef[4,1]

 pred=a+b*x+c*x2+d*x3

 if(graphit)
  {
   if(yrange[2]<0) yrange=c(1,max(y))
   if(newgraph) win.graph(height=4,width=6)
   if(add)
    {
     if(pts) points(x-31,y,pch=ptype)
     lines(x-31,pred,lty=ltype)
    }
   else
    {
     if(pts)
      {
       plot(x-31,y,pch=ptype,ylim=yrange,xlab="date",ylab=yaxslab)
       lines(x-31,pred,lty=ltype)
      }
     else plot(x-31,pred,pch=ptype,ylim=yrange,xlab="date",ylab=yaxslab,type="l",lty=ltype)
    }
  }

 if(poly==1)
  {
   peakday=(-a/b)
   peak=NA
  }
 else if(poly==2)
  {
   peakday=(-b/(2*c))
   peak=a-b^2/(4*c)
  }

 return(list(pred=pred,peak=peak,peakday=peakday))
}



# </source>
# </function>
# 
# 
# 
# <function>
# <name>
# colMedians
# </name>
# <description>
#  For convenient medians, like colMeans.

# </description>
# <arguments>
# 
# </arguments>
# <sample>
# 
# </sample>
# <source>
colMedians=function(mat,na.rm=TRUE)
  return(apply(mat,2,median,na.rm=na.rm))
# </source>
# </function>
# 
# 
# <function>
# <name>
# midPoint
# </name>
# <description>
#  Midpoint of any vector.

# </description>
# <arguments>
# 
# </arguments>
# <sample>
# 
# </sample>
# <source>
midPoint=function(x,na.rm=TRUE)
  return(min(x)+0.5*diff(range(x,na.rm=na.rm)))
# </source>
# </function>
# 
# 
