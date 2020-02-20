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
  medianf=log(colMedians(zs)) #prevoius zs=log(colmed(zs))
  meanf=log(colMeans(zs)) #  in r*2
  xes=(1:length(meanf)-1)#[1:200] #in r^200
  zs1=data.frame(xes,medianf,meanf)
  zs=melt(zs1, id.vars='xes',variable.name='Function')
  ggplot(zs,aes(x=xes,y=value,colour=Function))+geom_line()+labs(x="Iteration",y="log(norm(grad f))" , title=paste("Newton on",as.character(substitute(f))), subtitle="with x1,x2 in (-10,10)")
}
plottern2(f1,gf1,hf1,xs) 
plottern2(f2,gf2,hf2,xs) 
plottern2(f3,gf3,hf3,xs) 
plottern2(f4,gf4,hf4,xs) 
plottern2(f5,gf5,hf5,xs) 
xs


#############################################################################
#####################################
#graddesc
graddesc2=function(f,g,x){
  d=length(x)
  maxit=20000
  b=matrix(rep(NA,d*maxit),nrow=d)
  b[,1]=x
  x1=x
  it=2
  while( log(norm(g(x1),"2"))>-5 && it<maxit ){
    #  for (i in 2:n){
    b[,it] = b[,it-1]+bls1(f,g,b[,it-1])*(-g(b[,it-1]))
    x1 = b[,it-1]+bls1(f,g,b[,it-1])*(-g(b[,it-1]))
    it=it+1
  }
  #print(it)
  b[,it:maxit]=b[,it-1]
  return(b)
  #return( list(x = log(norm(x1,"2"))),b[,1:it])  #change to b instead of b[,n] for all numbers
}

plotterg2=function(f,g,xs){  #for newtons
  zs=matrix(0,nrow=(dim(xs)[1]),ncol=20000)
  for(j in 1:(dim(xs)[1])){
    x=unlist(xs[j,],use.names=F)
    ys=graddesc2(f,g,x)
    for(i in 1:20000)
    {
      zs[j,i]=(norm(g(ys[,i]),"2"))
    }
  }
  medianf=log(colMedians(zs)) #prevoius zs=log(colmed(zs))
  meanf=log(colMeans(zs)) #  in r*2
  xes=(1:length(meanf)-1)#[1:200] #in r^200
  zs1=data.frame(xes,medianf,meanf)
  zs=melt(zs1, id.vars='xes',variable.name='Function')
  ggplot(zs,aes(x=xes,y=value,colour=Function))+geom_line()+labs(x="Iteration",y="log(norm(grad f))" , title=paste("Steepest descent on",as.character(substitute(f))), subtitle="with x1,x2 in (-10,10)")
}
plotterg2(f1,gf1,xs)
plotterg2(f2,gf2,xs)
plotterg2(f3,gf3,xs)
plotterg2(f4,gf4,xs)
plotterg2(f5,gf5,xs)
