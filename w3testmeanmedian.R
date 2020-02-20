?title




x1=-5:5
y1=x1
xs=data.frame(expand.grid(x1,y1))
dim(xs)

plottern(f2,gf2,hf2,unlist(xs[,],use.names = F)) #change title every time..

unlist(xs[1,],use.names = F)
xs
#lol



xs[1,]
dim(xs)[1]
newt3=function(f,g,h,x){ #newt with BLS and eigenvalue transformation
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
  b[,(it-1):maxit]=b[,it-1]
  return(b)
}
test1=newt3(f1,gf1,hf1,x0)
test1[,2]
norm(test1[,2],"2")

if(norm(test1[,4],"2")==0) {
 print("a")
}

# log(norm(test1[,2],"2"))
dim(xs)
xs[1,]

test=matrix(c(1,2,3,4),nrow=2)
test
dim(xs)
colMeans(test)
colMedians(test)
dim(xs)[1]
test=colMeans(test)

dim(newt2(f1,gf1,hf1,x0))

plottern2=function(f,g,h,xs){  #for newtons
  zs=matrix(0,nrow=(dim(xs)[1]),ncol=20000)
  for(j in 1:(dim(xs)[1])){
      x=unlist(xs[j,],use.names=F)
      ys=newt3(f,g,h,x)
    for(i in 1:dim(ys)[1])
      {
      zs[j][i]=log(norm(g(ys[,i]),"2"))
      }
  }
  zs=rowMedians(zs)
  #  zs=log(zs)/(dim(xs)[1])
  xes=1:length((zs))-1
  zs1=data.frame(xes,zs)
#  zs1=data.frame(cbind(xes,zs))
  zs=melt(zs1, id.vars='xes')
  #  plot(1:length(zs),zs,xlab = "iteration",ylab="lognorm of gradient",type = "l")
  ggplot(zs,aes(xes,value))+geom_line()+labs(x="Iteration",y="log(norm(grad f))" , title="Newton on f2", subtitle="with xs")
}
plottern2(f2,gf2,hf2,xs)
plottern2(f1,gf1,hf1,xs)

install.packages(matrixStats)
library(matrixStats)
test=matrix(c(1,23,4,5,6,6),nrow=3)
colMedians(test)
length(rowMeans(test))
test


log(c(2,2))/2
plottern2(f1,gf1,hf1,xs)

plottern(f2,gf2,hf2,x0)
#















#########################################################################
plottern2=function(x){ #for newtons
  ys1=newt2(f1,gf1,hf1,x)
  ys2=newt2(f2,gf2,hf2,x)
  zs1=rep(NA,dim(ys1)[2])
  zs2=rep(NA,dim(ys2)[2])
  for(i in 1:length(zs)){
    zs1[i]=log(norm(gf1(ys1[,i]),"2")) #to avoid -infs
    zs2[i]=log(norm(gf2(ys1[,i]),"2")) #to avoid -infs
  }
  xs=1:max(dim(zs)[2],dim(zs)[2])-1
  zs=data.frame(xs,zs1,zs2)
  zs=melt(zs1, id.vars='xs',variable.name = "Function")
  #  plot(1:length(zs),zs,xlab = "iteration",ylab="lognorm of gradient",type = "l")
  ggplot(zs,aes(xs,value))+geom_line(aes(colour=Function))+labs(x="Iteration",y="log(norm(grad f))" , title="Newton on f2", subtitle="with x0=(-10,10)")
}
plottern2(x0)

(hf1(x0))


ys1=newt2(f1,gf1,hf1,x0)
ys2=newt2(f2,gf2,hf2,x0)
zs1=rep(NA,dim(ys1)[2])
zs2=rep(NA,dim(ys2)[2])
length(zs2)
ys2

zs2=log(norm(gf2(ys2),"2"))

for(i in 1:length(zs2)){
  zs1[i]=log(norm(gf1(ys1[,i]),"2")) #to avoid -infs
  zs2[i]=log(norm(gf2(ys1[,i]),"2")) #to avoid -infs
}
zs2
length(zs2)
plot(1:50,na.omit(zs2))
length(na.omit(zs2))

xs=1:max(dim(zs)[2],dim(zs)[2])-1
zs=data.frame(xs,zs1,zs2)
zs=melt(zs1, id.vars='xs',variable.name = "Function")
#  plot(1:length(zs),zs,xlab = "iteration",ylab="lognorm of gradient",type = "l")

ggplot(zs,aes(xs,value))+geom_line(aes(colour=Function))+labs(x="Iteration",y="log(norm(grad f))" , title="Newton on f2", subtitle="with x0=(-10,10)")





test=function(x){
  print(as.character(x))
}

test(f1)
as.character(substitute(f1))
f1$
test(123)


dat=melt(dat, id.vars='xs',variable.name="Function")
ggplot(dat,aes(xs,value))+geom_line(aes(colour=Function))+labs(x='Iteration',y='value',col='function',title='CG for x=(-10,10)')+
  scale_color_manual(labels = c("f1", "f2",'f3','f4','f5'), values = c(1,2,3,4,5))

xs=1:length(zs)
zs1=data.frame(zs)
zs1=data.frame(cbind(xs,zs))
zs=melt(zs1, id.vars='xs')
zs[1,]
library(reshape2)

str(f1)

character(substitute(toString(f1)

ggplot(zs,aes(xs,value))+geom_line()+labs(title=paste('NEWT',"on",2))

ggplot(zs,aes(xs,value))+geom_line()+labs(x="Iteration",y="log(norm(grad f))" , title=("Newton on",(f)), subtitle="with x0=(-10,10)")


plottern(f2,gf2,hf2,x0)



###############
ys1=graddesc(f1,gf1,x0)
dim(ys1)
log(norm(ys1[,1],"2"))
#ys=newt2(f2,gf2,hf2,c(-10,-20),100) #newt2 uses bls

ys1[,1]#x0
ys1[,100] #(1,1) is min for rosenbeck
ys1[,4422] #(1,1) is min for rosenbeck

zs=rep(NA,dim(ys1)[2])

ys1[,4423]

norm(gf1(ys1[,4422]),"2")

for(i in 1:length(zs)){
  zs[i]=log(norm(gf1(ys1[,i]),"2")) #to avoid -infs
}

plot(1:length(zs),zs,xlab = "iteration",ylab="lognorm of gradient",type = "l")


xs=1:length(zs)
zs1=data.frame(zs)
zs1=data.frame(cbind(xs,zs))
zs=melt(zs1, id.vars='xs')
zs[1,]
library(reshape2)
ggplot(zs,aes(xs,value))+geom_line()
#+geom_line


#+labs(x='Iteration',y='value',col='function',title='CG for x=(-10,10)')+
#  scale_color_manual(labels = c("f1", "f2",'f3','f4','f5'), values = c(1,2,3,4,5))



library(ggplot2)
opt3=function(x,n){ #CG
  print(x)
  xs=1:n
  ys=rep(NA,n)
  ys2=rep(NA,n)
  ys3=rep(NA,n)
  ys4=rep(NA,n)
  ys5=rep(NA,n)
  for (i in 1:n){
    rev=optim(x,f1,method="CG",control=list(maxit=i))
    ys[i]=log(norm(gf1(rev$par),"2"))
    rev2=optim(x,f2,method="CG",control=list(maxit=i))
    ys2[i]=log(norm(gf2(rev2$par),"2"))
    rev3=optim(x,f3,method="CG",control=list(maxit=i))
    ys3[i]=log(norm(gf3(rev3$par),"2"))
    rev4=optim(x,f4,method="CG",control=list(maxit=i))
    ys4[i]=log(norm(gf4(rev4$par),"2"))
    rev5=optim(x,f5,method="CG",control=list(maxit=i))
    ys5[i]=log(norm(gf5(rev5$par),"2"))
  }
  dat=data.frame(xs,ys,ys2,ys3,ys4,ys5)
  #  y=as.data.frame(cbind(ys,ys2,ys3,ys4,ys5))
  #  colnames(dat)=c("Iteration","Values")
  dat=melt(dat, id.vars='xs',variable.name="Function")
  ggplot(dat,aes(xs,value))+geom_line(aes(colour=Function))+labs(x='Iteration',y='value',col='function',title='CG for x=(-10,10)')+
    scale_color_manual(labels = c("f1", "f2",'f3','f4','f5'), values = c(1,2,3,4,5))
}
###########################################################

bls2=function(f,g,h,x){ #for newt
  c=10e-4  #c in (0,1) #c_1=10^-4 page 33 - wolf cond. c less than 0.5 goldstein
  rho=0.9 #rho in (0,1). why 0.01?  #close to 1 many iterations, close to 0 few i suppose.
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
#rho 0.8, 8.25e-18
bls2(f3,gf3,Bk1,x0)
#rho=0.1, 1e-18

hf3(x0)
x0
x0+bls2(f3,gf3,Bk1,x0)*(solve(Bk1))%*%gf3(x0)

x0+bls2(f3,gf3,Bk1,x0)*(solve(Bk1))%*%gf3(x0)
x0+bls2(f3,gf3,Bk1,x0)*t(solve(Bk1))*gf3(x0)

eigen(hf3(x0))$values


Ek1=diag(eigen(hf3(x0))$values)
Ek1[Ek1 <= 0] = (10^{-4}) #our machine precision, replacing all negative eigen values with it.
Bk1=Ek1
Bk1

newt3=function(f,g,h,x){ #newt with BLS and eigenvalue transformation
  d=length(x)
  maxit=20000
  b=matrix(rep(NA,d*maxit),nrow=d)
  b[,1]=x #x_0
  x1=x
  it=2
  while( log(norm(g(x1),"2"))>-10 && it<maxit ){
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
test=newt3(f3,gf3,hf3,x0)
test[,dim(test)[2]]
