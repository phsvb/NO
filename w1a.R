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

f1(c(1,2,3))
gf1(c(1,2,3))

f3(c(1,2,3))
gf3(c(1,2,3))
#gradient
gf1=function(x){
  alpha=1000
  d=length(x)
  g=rep(NA,d)
  for (i in 1:d){
    g[i]=alpha^((i-1)/(d-1))*2*x[i]
  }
  g
}
#hessian:
hf1=function(x){
  alpha=1000
  d=length(x)
  h=matrix(0,nrow=d,ncol=d)
  for (i in 1:d){           #will be non-zero only on the diagonal.
    h[i,i]=alpha^((i-1)/(d-1))*2 
  }
  h
}

gf1(c(0,0))#hence 0,0 is an extremum
hf1(c(0,0))#p.d., so (0,0) is a strict local minimizer of f.
#####################
#rosenbeck banana:
f2=function(x){
  if(length(x)==2){
    (1-x[1])^2+100*(x[2]-x[1]^2)^2}
  else NULL
}


optim(c(-2,5),f2)

?optim
optim(c(-2,5),f2,method="BFGS") #default is 
optim(c(-2,5),f2,method="CG") #default is 
?optim




#gradient of f2
gf2=function(x){ 
  g=rep(NA,2)
  if(length(x)==2){
    g[1]=2*x[1]-2+100*(4*x[1]^3-4*x[2]*x[1])
    g[2]=100*(2*x[2]-2*x[1]^2)
  }
  else NULL
  g
}

#hessian
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

gf2(c(1,1))#gradient is zero, hence this is a extremum
hf2(c(1,1)) # is a symmetric matrix. trace is positive, so is determinant 
det(hf2(c(1,1))) #hence is p.d.. So by 2.4 (1,1) is a strict local minimizer. 




#############
#log-ellipsoid:
f3=function(x){
  epsilon=10^{-16}
  return(log(epsilon+f1(x)))
}
#gradient #maybe?
gf3=function(x){
  epsilon=10^{-16}
  return(gf1(x)/(epsilon+f1(x)))
}
gf3(c(1,2,3))

#hessian
hf3=function(x){

}


#attractive sector functions:
f4=function(x){
  q=10^8
  h(x)=log(1+exp(q*x))/q
  return(
    sum(h(x)+100*h(-x))
  )
}

f5=function(x){
  q=10^8
  h(x)=log(1+exp(q*x))/q
  return(
    sum(h(x)^2+100*h(-x)^2)
  )
}

#########################
#to plot functions we adjust them from the general case to the 2case
pf1=function(x,y){
  alpha=1
  d=2
  c=0
  c=alpha^((1-1)/(2-1))*x^2+alpha^((2-1)/(2-1))*y^2
  return(c)
}
x=seq(-1,1,0.1)
y=seq(-1,1,0.1) 
z=outer(x,y,pf1)

persp(x, y, z,
      main="Plot of funtion 1-Ellipsoid function, with alpha=1",
      xlab="x in [-1,1]",
      ylab="y in [-1,1]",
      zlab = "f(x,y)",
      theta = 30, phi = 15,
      col = "springgreen", shade = 0.5)


##
##plot of rosenbeck f2
pf2=function(x,y){
  (1-x)^2+100*(y-x^2)^2
}
x=seq(-2,2,0.1)
y=seq(-1,3,0.1)
z=outer(x,y,pf2)
persp(x, y, z,
      main="Plot of funtion 2-Rosenbeck banana",
      xlab="x in [-2,2]",
      ylab="y in [-1,3]",
      zlab = "f(x,y)",
      theta = 120, phi = 15,
      col = "springgreen", shade = 0.5)



########plot of f3

pf3=function(x,y){
  epsilon=10^{-16}
  return(log(epsilon+pf1(x,y)))
}
x=seq(-0.01,0.01,0.001)
y=seq(-0.01,0.01,0.001)
z=outer(x,y,pf3)
persp(x, y, z,
      main="Plot of funtion 3",
      xlab="x in [-1,1]",
      ylab="y in [-1,1]",
      zlab = "f(x,y)",
      theta = 120, phi = 30,
      col = "springgreen", shade = 0.5)

######################plot of f4
pf3=function(x,y){
  epsilon=10^{-16}
  return(log(epsilon+pf1(x,y)))
}
x=seq(-1,1,0.1)
y=seq(-1,1,0.1)
z=outer(x,y,pf3)
persp(x, y, z,
      main="Plot of funtion 3",
      xlab="x in [-1,1]",
      ylab="y in [-1,1]",
      zlab = "f(x,y)",
      theta = 120, phi = 30,
      col = "springgreen", shade = 0.5)

############################
#this is so desperate
h=function(z){
  q=10^8
  return(log(1+exp(q*z))/q)
}
#much better
h=function(z){
  q=10^8
  return((log(1+exp(-abs(q*z)))+max(q*z,0))/q)
}
h(-1)
h(10)

pf4=function(x,y){
#  q=10^8
#  h=function(z){log(1+exp(q*z))/q}
  return(
    (h(x)+100*h(-x))+(h(y)+100*h(-y))
  )
}
x=seq(-10,10,1)
y=x
z=outer(x,y,pf4)
persp(x, y, z,
      main="Plot of funtion 4",
      xlab="x in [-1,1]",
      ylab="y in [-1,1]",
      zlab = "f(x,y)",
      theta = 120, phi = 30,
      col = "springgreen", shade = 0.5)
h(-1)
########################
#5
pf5=function(x,y){
  #  q=10^8
  #  h=function(z){log(1+exp(q*z))/q}
  return(
    (h(x)^2+100*h(-x)^2)+(h(y)^2+100*h(-y)^2)
  )
}
x=seq(-1*10^{-6},10^{-6},10^{-7})
y=x
z=outer(x,y,pf5)
persp(x, y, z,
      main="Plot of funtion 5",
      xlab="x in [-1,1]",
      ylab="y in [-1,1]",
      zlab = "f(x,y)",
      theta = 120, phi = 30,
      col = "springgreen", shade = 0.5)


exp(0.1)
exp(-0.1)
