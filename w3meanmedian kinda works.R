newt3=function(f,g,h,x){ #newt with BLS and eigenvalue transformation
  d=length(x)
  maxit=50
  b=matrix(rep(NA,d*maxit),nrow=d)
  b[,1]=x #x_0
  x1=x
  it=2
  while( log(norm(g(x1),"2"))>-5 && it<maxit ){
    if(mean(eigen(h(b[,it-1]))$values>0)==1){
      Bk=h(b[,it-1])
    }
    else{
      Ek=diag(eigen(h(b[,it-1]))$values)
      Ek[Ek <= 0] = (10^{-4}) #our machine precision, replacing all negative eigen values with it.
      Bk=Ek
    }
    b[,it] = b[,it-1]-bls2(f,g,Bk,b[,it-1])*solve(Bk)%*%g(b[,it-1])    #replace h in bls2 with Bk probably
    x1=b[,it]
    it=it+1
  }
  b[,(it-1):maxit]=b[,it-1]
  return(b)
}
plottern2=function(f,g,h,xs){  #for newtons
  zs=matrix(0,nrow=(dim(xs)[1]),ncol=50)
  for(j in 1:(dim(xs)[1])){
    x=unlist(xs[j,],use.names=F)
    ys=newt3(f,g,h,x)
    for(i in 1:50)
    {
      zs[j,i]=(norm(g(ys[,i]),"2"))
    }
  }
  zs=log(colMeans(zs))#[1:200]#/dim(xs)[1]
  xes=(1:length(zs)-1)#[1:200]
  zs1=data.frame(xes,zs)
  zs=melt(zs1, id.vars='xes')
  ggplot(zs,aes(xes,value))+geom_line()+labs(x="Iteration",y="log(norm(grad f))" , title=paste("Newton on",as.character(substitute(f))), subtitle="with x1,x2 in (-10,10)")
}
plottern2(f2,gf2,hf2,xs) #means
#change means, medians. perhabs maxit

plottern2(f2,gf2,hf2,xs) #medians
