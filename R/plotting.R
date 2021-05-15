#***************************************************************************----
# Main exported functions ----

#' Plot summary of a trimmed distribution: pdf + cdf + quantile + random
#'
#' @param d distributions3 object
#' @param trim Numeric vector of size 2, trimming bounds.
#' @param doTrunc Logical, do truncation? default FALSE, i.e. rectification is used.
#' @param xgrid Numeric vector, x-values where pdf and cdf are evaluated
#' @param pgrid Numeric vector in (0,1), probabilities where quantile function is evaluated
#' @param nsim Integer, number of simulated values
#' @param line.size,line.colour,area.fill,area.alpha,point.colour,point.size,lab.size,txt.size,
#'     graphical parameters
#' @return a ggplot object.
#' @examples
#' # Define Normal distribution
#' norm <- Normal(mu=0.75,sigma=0.5)
#' plot_trim(norm)
#' plot_trim(norm,trim=c(0,1))
#' plot_trim(norm,trim=c(0,1),doTrunc=TRUE)
#' @import ggplot2
#' @export

plot_trim <- function(d,trim=c(-Inf,Inf),doTrunc=FALSE,
                      xgrid=quantile_trim(d,seq(0.001,0.999,length.out=1000)),
                      pgrid=seq(0.001,0.999,length.out=1000),
                      nsim=200,
                      line.size=1,line.colour='black',
                      area.fill='salmon',area.alpha=0.5,
                      point.colour='black',point.size=2,
                      lab.size=12,txt.size=4){
  g1=plot_trim_pdf(d,trim,doTrunc,xgrid,
                   line.size,line.colour,
                   area.fill,area.alpha,
                   point.colour,point.size,
                   lab.size,txt.size)
  g2=plot_trim_cdf(d,trim,doTrunc,xgrid,
                   line.size,line.colour,lab.size)
  g3=plot_trim_quantile(d,trim,doTrunc,pgrid,
                        line.size,line.colour,lab.size)
  g4=plot_trim_random(d,trim,doTrunc,nsim,
                      point.colour,point.size,lab.size)
  g=gridExtra::grid.arrange(g1,g3,g2,g4,ncol=2)
  return(g)
}

#' Plot pdf of a trimmed distribution
#'
#' @param d distributions3 object
#' @param trim Numeric vector of size 2, trimming bounds.
#' @param doTrunc Logical, do truncation? default FALSE, i.e. rectification is used.
#' @param xgrid Numeric vector, x-values where pdf is evaluated
#' @param line.size,line.colour,area.fill,area.alpha,point.colour,point.size,lab.size,txt.size,
#'     graphical parameters
#' @return a ggplot object.
#' @examples
#' # Define Normal distribution
#' norm <- Normal(mu=0.75,sigma=0.5)
#' plot_trim_pdf(norm)
#' plot_trim_pdf(norm,trim=c(0,1))
#' plot_trim_pdf(norm,trim=c(0,1),doTrunc=TRUE)
#' @import ggplot2
#' @export

plot_trim_pdf <- function(d,trim=c(-Inf,Inf),doTrunc=FALSE,
                          xgrid=quantile_trim(d,seq(0.001,0.999,length.out=1000)),
                          line.size=1,line.colour='black',
                          area.fill='salmon',area.alpha=0.5,
                          point.colour='black',point.size=3,
                          lab.size=14,txt.size=5){
  x=xgrid
  # Add trimming bounds if they are in view
  if(trim[1] >= xgrid[1]){
    x=c(x,trim[1])
  }
  if(trim[2] <= xgrid[length(xgrid)]){
    x=c(x,trim[2])
  }
  # Compute pdf
  y=pdf_trim(d,x,trim,doTrunc)
  # Put NA at trimming bounds to make a disconnected line
  mask= x %in% trim
  y[mask]=NA
  # Start ggplot
  DF=data.frame(x=x,y=y)
  g=ggplot()+geom_line(data=DF,aes(x=x,y=y),size=line.size,colour=line.colour)
  g=g+xlab('x')+ylab('pdf f(x)')
  g=g+theme_bw()+theme(axis.title=element_text(size=lab.size))
  # Handle area + stem + text for rectified distribution
  p=getPs(d,trim)
  if(!doTrunc){
    if(trim[1] >= xgrid[1]){
      # Area before bound
      mask= (x<=trim[1])
      y=pdf_trim(d,x[mask])
      DF=data.frame(x=x[mask],y=y)
      g=g+geom_area(data=DF,aes(x=x,y=y),fill=area.fill,alpha=area.alpha)
      # stem
      DF=DF[DF$x==trim[1],]
      g=g+geom_segment(data=DF,aes(x=x,y=0,xend=x,yend=y))
      g=g+geom_point(data=DF,aes(x=x,y=y),colour=point.colour,size=point.size)
      # text
      label=paste0('p=',round(p$left,2),'  ')
      DF=cbind(DF,label=label)
      g=g+geom_text(data=DF,aes(x=x,y=y,label=label),hjust='right',size=txt.size)
    }
    if(trim[2] <= xgrid[length(xgrid)]){
      # Area after bound
      mask= (x>=trim[2])
      y=pdf_trim(d,x[mask])
      DF=data.frame(x=x[mask],y=y)
      g=g+geom_area(data=DF,aes(x=x,y=y),fill=area.fill,alpha=area.alpha)
      # stem
      DF=DF[DF$x==trim[2],]
      g=g+geom_segment(data=DF,aes(x=x,y=0,xend=x,yend=y))
      g=g+geom_point(data=DF,aes(x=x,y=y),colour=point.colour,size=point.size)
      # text
      label=paste0('  p=',round(p$right,2))
      DF=cbind(DF,label=label)
      g=g+geom_text(data=DF,aes(x=x,y=y,label=label),hjust='left',size=txt.size)
    }
  }
  return(g)
}

#' Plot cdf of a trimmed distribution
#'
#' @param d distributions3 object
#' @param trim Numeric vector of size 2, trimming bounds.
#' @param doTrunc Logical, do truncation? default FALSE, i.e. rectification is used.
#' @param xgrid Numeric vector, x-values where cdf is evaluated
#' @param line.size,line.colour,lab.size, graphical parameters
#' @return a ggplot object.
#' @examples
#' # Define Normal distribution
#' norm <- Normal(mu=0.75,sigma=0.5)
#' plot_trim_cdf(norm)
#' plot_trim_cdf(norm,trim=c(0,1))
#' plot_trim_cdf(norm,trim=c(0,1),doTrunc=TRUE)
#' @import ggplot2
#' @export

plot_trim_cdf <- function(d,trim=c(-Inf,Inf),doTrunc=FALSE,
                          xgrid=quantile_trim(d,seq(0.001,0.999,length.out=1000)),
                          line.size=1,line.colour='black',
                          lab.size=14){
  x=xgrid
  # Add trimming bounds if they are in view
  if(trim[1] >= xgrid[1]){
    x=c(x,trim[1])
  }
  if(trim[2] <= xgrid[length(xgrid)]){
    x=c(x,trim[2])
  }
  # Compute cdf
  y=cdf_trim(d,x,trim,doTrunc)
  # Put NA at trimming bounds to make a disconnected line
  if(!doTrunc){
    mask= x %in% trim
    y[mask]=NA
  }
  # Start ggplot
  DF=data.frame(x=x,y=y)
  g=ggplot()+geom_line(data=DF,aes(x=x,y=y),size=line.size,colour=line.colour)
  g=g+xlab('x')+ylab('cdf F(x)')
  g=g+theme_bw()+theme(axis.title=element_text(size=lab.size))
  return(g)
}

#' Plot quantile function of a trimmed distribution
#'
#' @param d distributions3 object
#' @param trim Numeric vector of size 2, trimming bounds.
#' @param doTrunc Logical, do truncation? default FALSE, i.e. rectification is used.
#' @param pgrid Numeric vector in (0,1), probabilities where quantile function is evaluated
#' @param line.size,line.colour,lab.size, graphical parameters
#' @return a ggplot object.
#' @examples
#' # Define Normal distribution
#' norm <- Normal(mu=0.75,sigma=0.5)
#' plot_trim_quantile(norm)
#' plot_trim_quantile(norm,trim=c(0,1))
#' plot_trim_quantile(norm,trim=c(0,1),doTrunc=TRUE)
#' @import ggplot2
#' @export

plot_trim_quantile <- function(d,trim=c(-Inf,Inf),doTrunc=FALSE,
                          pgrid=seq(0.001,0.999,length.out=1000),
                          line.size=1,line.colour='black',
                          lab.size=14){
  # Compute quantile
  x=pgrid
  y=quantile_trim(d,pgrid,trim,doTrunc)
  # Start ggplot
  DF=data.frame(x=x,y=y)
  g=ggplot()+geom_line(data=DF,aes(x=x,y=y),size=line.size,colour=line.colour)
  g=g+xlab('p')+ylab('quantile Q(p)')
  g=g+theme_bw()+theme(axis.title=element_text(size=lab.size))
  return(g)
}

#' Illustrate random function of a trimmed distribution
#'
#' @param d distributions3 object
#' @param trim Numeric vector of size 2, trimming bounds.
#' @param doTrunc Logical, do truncation? default FALSE, i.e. rectification is used.
#' @param nsim Integer, number of simulated values
#' @param point.size,point.colour,lab.size, graphical parameters
#' @return a ggplot object.
#' @examples
#' # Define Normal distribution
#' norm <- Normal(mu=0.75,sigma=0.5)
#' plot_trim_random(norm)
#' plot_trim_random(norm,trim=c(0,1))
#' plot_trim_random(norm,trim=c(0,1),doTrunc=TRUE)
#' @import ggplot2
#' @export

plot_trim_random <- function(d,trim=c(-Inf,Inf),doTrunc=FALSE,nsim=200,
                             point.colour='black',point.size=3,lab.size=14){
  # generate values
  x=1:nsim
  y=random_trim(d,nsim,trim,doTrunc)
  # Start ggplot
  DF=data.frame(x=x,y=y)
  g=ggplot()+geom_point(data=DF,aes(x=x,y=y),size=point.size,colour=point.colour)
  g=g+xlab('Index')+ylab('Realizations')
  g=g+theme_bw()+theme(axis.title=element_text(size=lab.size))
  return(g)
}
