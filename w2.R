#ellipsoid:
#f1=function(x,a){
#  sum(a^{(i-1)/(length(x)-1)}*x[i]^2)  
#}
#with forloops:

#ellipsoid:
f1=function(x){
  alpha=1000
  d=length(x)
  c=0
  for (i in 1:d){
    c=c+alpha^((i-1)/(d-1))*x[i]^2
    #    print(c)
  }
  c
}
gf1=function(x){
  alpha=1000
  d=length(x)
  g=rep(NA,d)
  for (i in 1:d){
    g[i]=alpha^((i-1)/(d-1))*2*x[i]
  }
  g
}
hf1=function(x){
  alpha=1000
  d=length(x)
  h=matrix(0,nrow=d,ncol=d)
  for (i in 1:d){           #will be non-zero only on the diagonal.
    h[i,i]=2*alpha^((i-1)/(d-1)) 
  }
  h
}
hf1(x1)
2*
#rosenbeck banana:
f2=function(x){
  if(length(x)==2){
    (1-x[1])^2+100*(x[2]-x[1]^2)^2}
  else NULL
}
gf2=function(x){ 
  g=rep(NA,2)
  if(length(x)==2){
    g[1]=2*x[1]-2+100*(4*x[1]^3-4*x[2]*x[1])
    g[2]=100*(2*x[2]-2*x[1]^2)
  }
  else NULL
  g
}

hf2=function(x){ 
  h=matrix(NA,nrow=2,ncol=2)
  if(length(x)==2){
    h[1,1]=2+100*(12*x[1]^2-4*x[2])
    h[2,2]=100*2
    h[1,2]=100*(-4*x[1])
    h[2,1]=h[1,2]
  }
  else NULL
  h
}

#############
#log-ellipsoid:
f3=function(x){
  epsilon=10^{-16}
  return(log(epsilon+f1(x)))
}
gf3=function(x){
  epsilon=10^{-16}
  return(gf1(x)/(epsilon+f1(x)))
}

#hf3=function(x){
#  epsilon=10^{-16}
#  alpha=1000
#  d=length(x)
#  H=matrix(0,nrow=d,ncol=d)
#  k2=0
#  for (k in 1:d){
#    k2= k2+alpha^((k-1)/(d-1))*x[k]^2
#  }
#  for(i in 1:d){
#      for(j in 1:d){
#        if(j!=i){
#        H[i,j]= (-4*x[i]*x[j]*alpha^((j+i-2)/(d-1)))/(epsilon+k2)^2
#          }}
#      H[i,i]=(2*alpha^((i-1)/(d-1))*(epsilon+k2)-(2*alpha^((i-1)/(d-1))*x[i])^2)/((epsilon+k2)^2)
#  }
#  H
#}
########this is nicer v
hf3=function(x){
  epsilon=10^{-16}
  alpha=1000
  return(hf1(x)/(epsilon+f1(x))-gf1(x)%*%t(gf1(x))/(epsilon+f1(x))^2)
}
hf4(x0)
hf4
hf3( x0)

eigen(hf3( x0))$values
te
te[1,2]
#much better
h=function(z){ #h function
  q=10^8
  return((log(1+exp(-abs(q*z)))+max(q*z,0))/q)
  }

hd=function(z){ #logistic function
  q=10^8
  return(plogis(q*z))
}
10^8


#attractive sector functions:
f4=function(x){
  q=10^8
 # h(x)=log(1+exp(q*x))/q
  return(
    sum(h(x)+100*h(-x))
  )
}
gf4=function(x){
  q=10^8
  plogis(q*x)-100*plogis(-q*x) #logistic function instead of  h, could also just do plogis(x)=hd(x)
}

diag(te)
hf4=function(x){
 q=10^8
 d=length(x)
 H=matrix(0,nrow=d,ncol=d)
 for (i in 1:d){
  if(x[i]<0){
       H[i,i]=101*q*exp(q*x[i])/(1+exp(q*x[i]))^2
    }
   else{
     H[i,i]=101*q*exp(-q*x[i])/(1+exp(-q*x[i]))^2
   }
  }
 H
}

f5=function(x){
  q=10^8
#  h(x)=log(1+exp(q*x))/q
  return(
    sum(h(x)^2+100*h(-x)^2)
  )
}


gf5=function(x){ #without forloops:
  q=10^8
  #  h(x)=log(1+exp(q*x))/q
  2*h(x)*hd(x)-200*h(-x)*hd(-x)
}

hf5=function(x){
  q=10^8
  d=length(x)
  H=matrix(0,nrow=d,ncol=d)
  for (i in 1:d){
   H[i,i]=(-200*exp(q*x[i])*log(1+exp(-q*x[i]))+2*exp(q*x[i])*log(1+exp(q*x[i]))+2*exp(2*q*x[i])-200)/(1+exp(q*x[i]))^2 
  }
  H
}

#gf5=function(x){
#  q=10^8
#  d=length(x)
#  g=rep(NA,d)
#  #  h(x)=log(1+exp(q*x))/q
#  for (i in 1:d){
#    g[i]=2*h(x[i])*hd(x[i])-200*h(-x[i])*hd(-x[i])
#    
#  }
#  g
#}
hf5=function()

library(ggplot2)
opt=function(x,f,g,n){
  xs=seq(1,n,1)
  ys=rep(NA,n)
  for (i in 1:n){
    rev=optim(x,f,method="BFGS",control=list(maxit=i))
    ys[i]=log(norm(g(rev$par),"2"))
    if (i==n){
      print(rev$par)
    }
  }
  print(x0)
  dat=cbind(xs,ys)
  dat=as.data.frame(dat)
  colnames(dat)=c("Iteration","Value")
  ggplot(dat,aes(Iteration,Value))+geom_line()
}
x0=rnorm(2,0,10)
x0=c(-10^{-17},10^{-17})
x0=c(-10^{-17})
norm(x0,"2")


y0=c(5.521419e-16,2.538079e-16)
gf3(y0)
y0
gf3(c(0,0))
sqrt(10^8)
optim(x0,f3,method="BFGS")
?optim
x0=c(-10,10)
x0=c(10,10)
opt(x0,f3,gf3,100)

gf3(c(-10,10))
gf3(y0)

plogis(-10) == 1-plogis(10)
plogis(-10)
1-plogis(10)

gf4(c(0,0))

f3
gf3(c(5.521419e-16, 2.538079e-16))
gf3(c(0,0))
opt(x0,f1,gf1,100)
opt(x0,f2,gf2,100)
opt(x0,f3,gf3,100)
opt(x0,f4,gf4,100)
opt(x0,f5,gf5,100)

length(seq(0,10,1))

library(reshape2)
opt2=function(x,n){ #BFGS
  xs=1:n
  ys=rep(NA,n)
  ys2=rep(NA,n)
  ys3=rep(NA,n)
  ys4=rep(NA,n)
  ys5=rep(NA,n)
  for (i in 1:n){
    rev=optim(x,f1,method="BFGS",control=list(maxit=i))
    ys[i]=log(norm(gf1(rev$par),"2"))
    rev2=optim(x,f2,method="BFGS",control=list(maxit=i))
    ys2[i]=log(norm(gf2(rev2$par),"2"))
    rev3=optim(x,f3,method="BFGS",control=list(maxit=i))
    ys3[i]=log(norm(gf3(rev3$par),"2"))
    rev4=optim(x,f4,method="BFGS",control=list(maxit=i))
    ys4[i]=log(norm(gf4(rev4$par),"2"))
    rev5=optim(x,f5,method="BFGS",control=list(maxit=i))
    ys5[i]=log(norm(gf5(rev5$par),"2"))
  }
  dat=data.frame(xs,ys,ys2,ys3,ys4,ys5)
#  y=as.data.frame(cbind(ys,ys2,ys3,ys4,ys5))
#  colnames(dat)=c("Iteration","Values")
  dat=melt(dat, id.vars='xs',variable.name="Function")
  ggplot(dat,aes(xs,value))+geom_line(aes(colour=Function))+labs(x='Iteration',y='value',col='function',title='BFGS for x=(-10,10)')+
  scale_color_manual(labels = c("f1", "f2",'f3','f4','f5'), values = c(1,2,3,4,5))
}

x0=c(-10,10)
opt2(x0,115)
x0

#####################################
t0=1:10
t1=rnorm(10,0,5)
t2=rnorm(10,0,10)
df=data.frame(t0,t1,t2)
df=melt(df,id.vars = 't0',variable.name = 'series')
ggplot(df,aes(t0,value))+geom_line(aes(colour=series))+labs(x='LOL',y='rofl',title='tittles')

data.frame(cbind(t0,t1,t2))
##########################################
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
?optim
x0=c(10,10)
optim(x0,f3,method="BFGS")$par

test=optim(x0,f3,method="BFGS")$par
test
gf3(test)
norm(gf3(test),"2")
optim(x0,f2,method="BFGS",control=list(maxit=80))

optim(x0,f5,method="CG",control=list(maxit=1))

opt3(x0,55)

