#refer to w2.R for functions f,g,h 1-5
x0=c(-10,10)

bls2=function(f,g,h,x){ #for newt
  c=10e-4  #c in (0,1) #c_1=10^-4 page 33 - wolf cond. c less than 0.5 goldstein
  rho=0.5 #rho in (0,1). why 0.01?  #close to 1 many iterations, close to 0 few i suppose.
  a=1#starting value
  B=h
  p=-solve(B)%*%g(x) 
  #it=1                  #if we want to count iterations
  if(f(x+a*p) <= f(x)+c*a*t(g(x)%*%p)){ 
    #   print(it)
    return(a)
  }
  repeat{
    #    it=it+1
    a=rho*a
    if(f(x+a*p) <= f(x)+c*a*t(g(x)%*%p)){
      break
    }
  }
  #  print(it)
  return(a)
}

newt2=function(f,g,h,x){ #newt with BLS and eigenvalue transformation
  d=length(x)
  maxit=20000
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
  print(it)
  return(b[,1:it-1])
}
test=newt2(f3,gf3,hf3,x0)
test[,dim(test)[2]]

newt2(f1,gf1,hf1,x0)
newt2(f2,gf2,hf2,x0)
newt2(f3,gf3,hf3,x0)
newt2(f2,gf2,hf2,x0)

log(norm(gf1(x0),"2"))
###############################################################################################################################
###########for grad desc
bls1=function(f,g,x){ #for graddesc
  c=10e-4 #c in (0,1) #c_1=10^-4 page 33 - wolf cond. c less than 0.5 goldstein
  rho=0.25 #rho in (0,1). why 0.01?  #close to 1 many iterations, close to 0 few i suppose.
  #  ah=1         #set =1 for newton and quasi newton
  a=0.1 #for graddesc try small a 
  p=-g(x) 
  #  it=1
  if(f(x+a*p) <= f(x)+c*a*t(g(x)%*%p)){ 
    #    print(it)
    return(a)
  }
  repeat{
    #   it=it+1
    a=rho*a
    if(f(x+a*p) <= f(x)+c*a*t(g(x)%*%p)){
      break
    }
  }
  # print(it)
  return(a)
}
test1=graddesc(f2,gf2,x0)
test1[,dim(test1)[2]]
dim(test1)

graddesc=function(f,g,x){
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
  return(b[,1:it-1])
  #return( list(x = log(norm(x1,"2"))),b[,1:it])  #change to b instead of b[,n] for all numbers
}
##########################################
gf3(c( 3.326605e-27,1.321139e-17))
gf3(c(-8.354678e-20,-1.852140e-17))
#plots
plotterg=function(f,g,x){ #for gradient descent             #perhaps let c,rho,a be changable values here.
  ys1=graddesc(f,g,x)
  zs=rep(NA,dim(ys1)[2])
  for(i in 1:length(zs)){
    zs[i]=log(norm(g(ys1[,i]),"2")) #to avoid -infs
  }
  xs=1:length(zs)-1
  zs1=data.frame(zs)
  zs1=data.frame(cbind(xs,zs))
  zs=melt(zs1, id.vars='xs')

    print(ys1[,dim(ys1)[2]]) #gives final value
  ggplot(zs,aes(xs,value))+geom_line()+labs(x="Iteration",y="log(norm(grad f))" , title=paste("Steepest descent on",as.character(substitute(f))), subtitle="with x0=(-10,10)")
#  plot(1:length(zs),zs,xlab = "iteration",ylab="lognorm of gradient",type = "l")
}

plotterg(f1,gf1,x0) 
plotterg(f2,gf2,x0) #rho 0.2 #or 25
plotterg(f3,gf3,x0) #0.25
plotterg(f4,gf4,x0)
plotterg(f5,gf5,x0)


c(1,2,3,4)[-1]

log(norm(gf1(c(-9,7)),"2"))

plotterg(f2,gf2,x0)
plotterg(f3,gf3,x0)
plotterg(f4,gf4,x0)
plotterg(f5,gf5,x0)

c(1,2,3)[,1]


######################
plottern=function(f,g,h,x){ #for newtons
  ys1=newt2(f,g,h,x)
  zs=rep(NA,dim(ys1)[2])
  for(i in 1:length(zs)){
    zs[i]=log(norm(g(ys1[,i]),"2")) #to avoid -infs
  }
  xs=1:length(zs)-1
  zs1=data.frame(zs)
  zs1=data.frame(cbind(xs,zs))
  zs=melt(zs1, id.vars='xs')
  #  plot(1:length(zs),zs,xlab = "iteration",ylab="lognorm of gradient",type = "l")
  ggplot(zs,aes(xs,value))+geom_line()+labs(x="Iteration",y="log(norm(grad f))" , title=paste("Newton on",as.character(substitute(f))), subtitle="with x0=(-10,10)")
}
plottern(f1,gf1,hf1,x0)
plottern(f2,gf2,hf2,x0)
plottern(f3,gf3,hf3,x0*10e-10)

plottern(f4,gf4,hf4,x0*10e-12)
plottern(f5,gf5,hf5,x0*10e-8)

hf3(x0)

hf4(x0*10e-8)
hf4(x0*10e-10)
hf4(x0*10e-12)

hf5(x0*10e-6)
hf5(x0*10e-10)
hf4(x0*10e-10)

eigen(hf3(x0*10e-10))$values
x0*10e-10

toString(123)
String(123)

log(norm(gf3(c(-2510,512)),"2"))
