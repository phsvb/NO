bls2=function(f,g,h,x){ #for newt
  c=0.5 #c in (0,1) #c_1=10^-4 page 33 - wolf cond. c less than 0.5 goldstein
  rho=0.1 #rho in (0,1). why 0.01?
  #  ah=1         #set =1 for newton and quasi newton
  a=1
  B=h
  p=-solve(B)%*%g(x) 
  if(f(x+a*p) <= f(x)+c*a*t(g(x)%*%p)){return(a)}
  repeat{
    a=rho*a
    if(f(x+a*p) <= f(x)+c*a*t(g(x)%*%p)){
      break
    }
  }
  return(a)
}

newt=function(f,g,h,x,n){ #newt with BLS and eigenvalue transformation
  d=length(x)
  b=matrix(rep(NA,d*n),nrow=d)
  b[,1]=x #x_0
  for (i in 2:n){
    if(mean(eigen(h(b[,i-1]))$values>0)==1){
      Bk=h(b[,i-1])
    }
    else{
      Ek=diag(eigen(h(b[,i-1]))$values)
      Ek[Ek < 0] = sqrt(10^{-16}) #our machine precision, replacing all negative eigen values with it.
      Bk=Ek
      }
    #pk=-B_k^-1 g(x)
    #      Bk=diag(d)%*%diag(eigen(h(b[,i-1]))$values)%*%diag(d)
    #      Pk=b[,i-1]-bls1(f,g,b[,i-1],1)*g(b[,i-1])%*%solve(Bk)
    #Pk #b[,i-1]+bls1(f,g,b[,i-1],0.001)*(-gf1(b[,i-1]))
    b[,i] = b[,i-1]-bls2(f,g,Bk,b[,i-1])*solve(Bk)%*%g(b[,i-1])    #replace h in bls2 with Bk probably
  }
  b  #change to b instead of b[,n] for all numbers
}


log(norm(gf1(ys[,1]),"2"))
log(norm(gf1(ys[,1]),"2"))
length(gf1(ys))
plot(1:20,log(norm(gf1(ys),"2")))
##########################################

ys=newt(f2,gf2,hf2,x0,1000)#good
ys[,1]#x0
ys[,400] #(1,1) is min for rosenbeck
zs=rep(NA,dim(ys)[2])

for(i in 1:length(zs)){
  zs[i]=log(norm(gf1(ys[,i]),"2")) #to avoid -infs
}
plot(1:length(zs),zs)
#shows convergence^
#######################################
gf1(c(0,0))
gf2(c(1,1))
newt(f2,gf2,hf2,x0,2)#so slow..

newt(f2,gf2,hf2,x0,200)#so slow..
hf2(x0)

newt(f3,gf3,hf3,c(1,10),100)

newt(f4,gf4,hf4,c(1,10),2)

gf3(c(0,0))
x0
newt(f3,gf3,hf3,x0,10000)
mean(eigen(hf2(x0))$values>0)




mean(eigen(hf1(x0))$values>0)==1


diag(eigen(hf2(x0))$values)


te2=diag(c(1:2,-3))
te2
#te2<0
te2[te2 < 0] <- 10^{-8}
te2
eigen(te2)$values
te3=diag(eigen(te2)$values)
te3[te3 < 0] <- 10^{-8}

te2
te[te2 < 0] <- 10^{-8}
Ek=diag(eigen(hf2(x0))$values)

diag(eigen(hf2(x0))$values)[diag(eigen(hf2(x0))$values)<0]=10^{-8}
Ek
Ek[2,2]=-1
Ek
Ek[Ek < 0] = 10^{-8}

Ek


library(tidyr)
x=-5:5
y=x
expand.grid(x,y)
