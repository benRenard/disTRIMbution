f=file.path('data-raw','CowhouseCreek.txt')
Q=read.table(f,header=TRUE)
years=min(Q$year):max(Q$year)
n=length(years)
CowhouseCreek=data.frame(year=integer(n),
                         min=numeric(n),
                         duration=numeric(n))
for(i in 1:n){
  mask=Q$year==years[i]
  CowhouseCreek$year[i]=years[i]
  CowhouseCreek$min[i]=min(Q$streamflow[mask])
  CowhouseCreek$duration[i]=sum(Q$streamflow[mask]<=0.01)/sum(Q$streamflow[mask]>=0)
}
save(CowhouseCreek,file=file.path('data','CowhouseCreek.RData'))
