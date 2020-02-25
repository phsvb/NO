x0


hf3(x0)


B=diag(eigen(hf3(x0))$values)         #spectral decomp w. machine presicion
for (i in 1:dim(B)[1]){
  if (B[i,i]<=0){
    B[i,i]=1e-4
  }
}
B

#p(lambda)=
-(solve(B)+1e-04*diag(dim(B)[1]))%*%

diag(3)

eigen(B+lambda*diag(dim(B)[1]))$values>=0

if( mean(eigen(B+lambda*diag(dim(B)[1]))$values>=0)==1) 
  
  
#check thm 4.1

p=function(L){
  return(-solve(B+L*diag(2))%*%gf3(x0))
}

l0=1e-4
l1=10
p(l0)
p(l1)

norm(p(l0),"2")
norm(p(l1),"2")
#pretend delta is 10
delta=10
(norm(p(l0),"2")>delta && delta>norm(p(l1),"2"))
norm(p(10/2),"2")
norm(norm(p(10/2),"2") - 10,"2")
#################################################
########working implementation of bijection algorithm
if(norm(p(l0),"2")>delta && delta>norm(p(l1),"2")){
  Lt=(l0+l1)/2
  while (norm(norm(p(Lt),"2")-delta,"2")>1e-2){
    if (norm(p(Lt),"2")>delta){
      l0=Lt
    }
    else{
      l1=Lt
    }
    Lt=(l0+l1)/2
  }
  Lt
}
Lt

Ps=-solve(B+Lt*diag(2))%*%gf3(x0)
norm(Ps,"2")
Lt*(delta-norm(Ps,"2"))

delta
Ps
#rho
B
mk=function(x,f,g,B,p){
  return(f(x)+t(g(x))%*%p+1/2*t(p)%*%B%*%p)
}
mk(x0,f3,gf3,B,Ps)
mk(x0,f3,gf3,B,0)
gf3(x0)

(f3(x0)-f3(x0+Ps))/(mk(x0,f3,gf3,B,c(0,0))-mk(x0,f3,gf3,B,Ps)) #3.47 big!
 

p(Lt)
norm(p(Lt),"2")

p*=-solve(B+lamba*diag(dim(B)[1]))%*%g

-solve(B+1*diag(dim(B)[1]))%*%gf3(x0)

norm(-solve(B+1*diag(dim(B)[1]))%*%gf3(x0),"2")
1*()

(B+diag(2))


p=function(la){
  return(-solve(B+l*diag(2))%*%g)
}
norm(p(l0),"2")

#change diag(2) to diag(dim(B)[1])

temp1=function(f,g,h,x,lamax=10){            #lamax is probably delta max, delta_0 chosen somehow
  mk=function(x,f,g,B,p){
    return(f(x)+t(g(x))%*%p+1/2*t(p)%*%B%*%p)
  }
  p=function(la){
    return(-solve(B+la*diag(2))%*%g(x)) 
  }
  d=length(x)
  n=100  
  xs=matrix(NA,nrow=d,ncol=d)
  xs[,1]=x
  delta=5       #max range of cicle. Right?
  for (i in 2:n){
    if(mean(eigen(h(xs[,i-1]))$values>0)==1){
      B=h(xs[,i-1])
    }
    else{
      Ek=diag(eigen(h(xs[,i-1]))$values)
      for (j in 1:dim(Ek)[1]){                #DO LIKE THIS FAM
        if(Ek[j,j] <= 0)
        { 
          Ek[j,j]= (10^{-4})
        } #our machine precision, replacing all negative eigen values with it.
        B=Ek
      } 
    #then now we got a B
    #to find a lambda we use the bijection algorithm
    }#found B now 
    #now to find lambda 
    L0=1e-4
    L1=lamax
    if(norm(p(L0),"2")>delta && delta>norm(p(L1),"2")){
      Lt=(L0+L1)/2
      while (norm(norm(p(Lt),"2")-delta,"2")>1e-6){ #threshold. make smaller
        if (norm(p(Lt),"2")>delta){
          l0=Lt
        }
        else{
          l1=Lt
        }
        Lt=(l0+l1)/2
      }
      Lt
    }
    #found lambda... what do now
    print(Lt) #by thm 4.1 we can now find p^*
    Pstar=-solve(B+Lt*diag(2))%*%g(xs[,i-1]) #found P_k now. can change to k-valued
    rho= f(xs[,i-1]-f(xs[,i-1]+Ps))/(mk(xs[,i-1],f,g,B,c(0,0))-mk(xs[,i-1],f,g,B,Pstar))
    if (rho <1/4){
      delta[k+1]=1/4*delta[k]
    }
    else{
      if (rho>3/4 && norm(norm(Ps,"2")-delta[k],"2")<1e-6){ #threshold can be changed
        delta[k+1]=min(2*delta[k],deltahat)
      }
      else{
        delta[k+1]=delta[k]
      }
    }
    if (rho>eta){
      xs[,i]=xs[,i-1]+Pstar
    }
    else{
      xs[,i]=xs[,i-1]
    }
  }
#><
}
#delta hat is overal bound on step lengths...

l0=1e-4
l1=lamax

eigen(hf3(x0))$values
################
temp1=function(f,g,h,x,deltamax,delta0){            #lamax is probably delta max, delta_0 chosen somehow
  eta=1/8
  mk=function(x,f,g,B,p){
    return(f(x)+t(g(x))%*%p+1/2*t(p)%*%B%*%p)
  }
  p=function(la){
    return(-solve(B+la*diag(dim(B)[1]))%*%g(x)) 
  }
  d=length(x)
  n=10 
  xs=matrix(NA,nrow=d,ncol=n)
  xs[,1]=x
  deltas=rep(NA,n)
  deltas[1]=delta0
  for (i in 2:n){
    if(mean(eigen(h(xs[,i-1]))$values>0)==1){
      B=h(xs[,i-1])
    }
    else{
      Ek=diag(eigen(h(xs[,i-1]))$values)
      for (j in 1:dim(Ek)[1]){                #DO LIKE THIS FAM
        if(Ek[j,j] <= 0)
        { 
          Ek[j,j]= 1e-4
        } #our machine precision, replacing all negative eigen values with it.
        B=Ek
      } 
      #then now we got a B
      #to find a lambda we use the bijection algorithm
    }#found B now 
    #now to find lambda 
    L0=1e-4       #lambda 0
    L1=deltamax   #deltamax
    if(norm(p(L0),"2") > deltas[i-1] && norm(p(L1),"2") < deltas[i-1]){         #using the bijective algorithm to find lambda
      Lt=(L0+L1)/2
      while (norm(norm(p(Lt),"2")-deltas[i-1],"2")>1e-6){ #threshold. make smaller
        if (norm(p(Lt),"2")>delta[i-1]){
          l0=Lt
        }
        else{
          l1=Lt
        }
        Lt=(l0+l1)/2
      }
      Lt
    }
    #found lambda... what do now
    #    print(Lt) #by thm 4.1 we can now find p^*
    Pstar=-solve(B+Lt*diag(dim(B)[1]))%*%g(xs[,i-1]) #found P_k now. can change to k-valued
    rho= (f(xs[,i-1])-f(xs[,i-1]+Pstar))/(mk(xs[,i-1],f,g,B,c(0,0))-mk(xs[,i-1],f,g,B,Pstar))
    if (rho <1/4){
      deltas[i]=1/4*deltas[i-1]
    }
    else{
      if (rho>3/4 && norm(norm(Ps,"2")-deltas[i-1],"2")<1e-6){ #threshold can be changed
        deltas[i]=min(2*deltas[i-1],deltamax) 
      }
      else{
        deltas[i]=deltas[i-1]
      }
    }
    if (rho>eta){
      xs[,i]=xs[,i-1]+Pstar
    }
    else{
      xs[,i]=xs[,i-1]
    }
  }
  print(xs)
  print(deltas)#><
#  xs
#  deltas
}
#temp1(f1,gf1,hf1,x0,0,0) #not too bad
temp1(f1,gf1,hf1,x0,1,0.1) #not too bad

temp1(f1,gf1,hf1,x0,100,10) #not too bad

temp1(f2,gf2,hf2,x0,10,11) #lol yikes

temp1(f4,gf4,hf4,x0,10,1)#lol

