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
    #      Bk=diag(d)%*%diag(eigen(h(b[,i-1]))$values)%*%diag(d)
    #      Pk=b[,i-1]-bls1(f,g,b[,i-1],1)*g(b[,i-1])%*%solve(Bk)
    b[,i] = b[,i-1]-bls2(f,g,h,b[,i-1])*solve(Bk)%*%g(b[,i-1])    #replace h in bls2 with Bk probably
    #Pk #b[,i-1]+bls1(f,g,b[,i-1],0.001)*(-gf1(b[,i-1]))
  }
  b[,n]  #change to b instead of b[,n] for all numbers
}
newt(f1,gf1,hf1,x0,2)

x0
newt(f3,gf3,hf3,x0,10000)
mean(eigen(hf2(x0))$values>0)




mean(eigen(hf1(x0))$values>0)==1


diag(eigen(hf2(x0))$values)


te2=diag(c(1:2,-3))
te2
te2<0
te2[te2 < 0] <- 10^{-8}
te2

te2
te[te2 < 0] <- 10^{-8}
Ek=diag(eigen(hf2(x0))$values)

diag(eigen(hf2(x0))$values)[diag(eigen(hf2(x0))$values)<0]=10^{-8}
Ek
Ek[2,2]=-1
Ek
Ek[Ek < 0] = 10^{-8}

Ek
