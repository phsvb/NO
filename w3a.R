bls1=function(f,g,x,ah){
  c=0.5 #c in (0,1) #c_1=10^-4 page 33 - wolf cond. c less than 0.5 goldstein
  rho=0.1 #rho in (0,1). why 0.01?
#  ah=1         #set =1 for newton and quasi newton
  a=ah
  p=-g(x)#/norm(g(x),"2") #pk with grad. desc. 
  if(f(x+a*p) <= f(x)+c*a*t(g(x)%*%p)){return(a)}
  repeat{
    a=rho*a
    if(f(x+a*p) <= f(x)+c*a*t(g(x)%*%p)){
      break
    }
  }
  return(a)
}

f1(x0-gf1(x0))
f1(x0)+


bls1(f1,gf1,x0)

floor(runif(1, 0, 10))

graddesc=function(f,g,x,n){
  d=length(x)
  b=matrix(rep(NA,2*n),nrow=d)
  b[,1]=x
  for (i in 2:n){
    b[,i] = b[,i-1]+bls1(f,g,b[,i-1],0.001)*(-gf1(b[,i-1]))
  }
  b[,n]  #change to b instead of b[,n] for all numbers
}


graddesc(f1,gf1,x0,2)
graddesc(f1,gf1,x0,10)
#x1 is very slowly...



x0
matrix(I,2)
diag(2)

eigen(hf1(x0))$values
diag(2)%*%diag(eigen(hf1(x0))$values)%*%diag(2)

eigen(hf1(x0))$values


solve(diag(2)%*%diag(2))


newt=function(f,g,h,x,n){
  d=length(x)
  b=matrix(rep(NA,d*n),nrow=d)
  b[,1]=x #x_0
  for (i in 2:n){
      Bk=h(b[,i-1])
#      Bk=diag(d)%*%diag(eigen(h(b[,i-1]))$values)%*%diag(d)
#      Pk=b[,i-1]-bls1(f,g,b[,i-1],1)*g(b[,i-1])%*%solve(Bk)
      b[,i] = b[,i-1]-g(b[,i-1])%*%solve(Bk)
        #Pk #b[,i-1]+bls1(f,g,b[,i-1],0.001)*(-gf1(b[,i-1]))
  }
  b[,n]  #change to b instead of b[,n] for all numbers
}
x0
#one part of it is converging... the other is not
newt(f1,gf1,hf1,x0,2)
newt(f1,gf1,hf1,x0,4)
newt(f1,gf1,hf1,x0,10)

f1
hf1
gf1


x0
gf1(x0)

x0-gf1(x0)%*%solve(hf1(x0))

gf1(x0)
hf1(x0)


-solve(hf1(x0))%*%gf1(x0)
#################################for f1: convergence after 1 iteratation
newt=function(f,g,h,x,n){ #newt with BLS. no transformation of hessian
  d=length(x)
  b=matrix(rep(NA,d*n),nrow=d)
  b[,1]=x #x_0
  for (i in 2:n){
    Bk=h(b[,i-1])
    #      Bk=diag(d)%*%diag(eigen(h(b[,i-1]))$values)%*%diag(d)
    #      Pk=b[,i-1]-bls1(f,g,b[,i-1],1)*g(b[,i-1])%*%solve(Bk)
    b[,i] = b[,i-1]-bls2(f,g,h,b[,i-1])*solve(Bk)%*%g(b[,i-1])
    #Pk #b[,i-1]+bls1(f,g,b[,i-1],0.001)*(-gf1(b[,i-1]))
  }
  b[,n]  #change to b instead of b[,n] for all numbers
}
#for f1 conv after 1 iteration
newt(f1,gf1,hf1,x0,5)
newt(f2,gf2,hf2,x0,1000)
newt(f3,gf3,hf3,x0,4)

eigen(hf2(x0))$values>0
eigen(hf2(x0))$values[1]

if(eigen(hf2(x0))$values>0){
  print("yes")
}


if(c(eigen(hf2(x0))$values)>c(0,0)){
  print("yes")
}

c(eigen(hf2(x0))$values)>c(0,0)

mean(eigen(hf2(x0))$values>0)

##
newt2=function(f,g,h,x,n){#newt without BLS
  d=length(x)
  b=matrix(rep(NA,d*n),nrow=d)
  b[,1]=x #x_0
  for (i in 2:n){
    Bk=h(b[,i-1])
    #      Bk=diag(d)%*%diag(eigen(h(b[,i-1]))$values)%*%diag(d)
    #      Pk=b[,i-1]-bls1(f,g,b[,i-1],1)*g(b[,i-1])%*%solve(Bk)
    b[,i] = b[,i-1]-solve(Bk)%*%g(b[,i-1])
    #Pk #b[,i-1]+bls1(f,g,b[,i-1],0.001)*(-gf1(b[,i-1]))
  }
  b  #change to b instead of b[,n] for all numbers
}

newt2(f2,gf2,hf2,x0,5)
newt(f2,gf2,hf2,x0,5)

bls2()
solve(hf1(x0))

solve(hf1(x0))%*%gf1(x0)
te[,1]-1*solve(hf1(x0))%*%gf1(x0)

x0

newt(f2,gf2,hf2,x0,5)

newt(f3,gf3,hf3,x0,10)

bls1(f1,gf1,x0,1)

x0-solve(hf1(x0))%*%gf1(x0)
x0-bls2(f1,gf1,hf1,x0)*solve(hf1(x0))%*%gf1(x0)

bls2=function(f,g,h,x){ #for newt
  c=0.5 #c in (0,1) #c_1=10^-4 page 33 - wolf cond. c less than 0.5 goldstein
  rho=0.1 #rho in (0,1). why 0.01?
  #  ah=1         #set =1 for newton and quasi newton
  a=1
  p=-solve(h(x))%*%g(x) #/norm(g(x),"2") #pk with grad. desc. 
  if(f(x+a*p) <= f(x)+c*a*t(g(x)%*%p)){return(a)}
  repeat{
    a=rho*a
    if(f(x+a*p) <= f(x)+c*a*t(g(x)%*%p)){
      break
    }
  }
  return(a)
}

bls2(f1,gf1,hf1,x0)

x0-bls2(f1,gf1,hf1,x0)*gf1(x0)%*%solve(hf1(x0))

x0-bls2(f1,gf1,hf1,x0)*solve(hf1(x0))%*%gf1(x0)
