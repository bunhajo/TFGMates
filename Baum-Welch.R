#ALGORITMO BAUM-WELCH
BaumWelch = function(returns, mu, sigma, n_states=3, Tolerance=7*10^{-2}, maxstep=1000){

  change_likelihood=c(rep(Inf, n_states))
  likelihood=data.frame()
  
  returns=as.data.frame(returns)
  mu = as.data.frame(mu)
  A=data.frame(rep(1/n_states,n_states))
  A[1:n_states]=rep(1/n_states,n_states)
  # Nuestra Matriz de transción incial
  # A=data.frame(c(0.9,0.1),c(0.1,0.9))
  p=rep(1/n_states,n_states)
  # Nuestras probabilidades iniciales
  # p = c(0.6, 0.4)
  k=ncol(returns)
  if(k==1){sigma=as.data.frame(sigma)}
  L=nrow(returns)
  
  B=data.frame(c(rep(0,L)))
  B[1:n_states]=c(rep(0,L))
  forward=B
  backward=B
  smoothed=B
  xi=vector("list", L-1)
  xi[1:L-1]=list(data.frame(c(rep(0,n_states)),c(rep(0,n_states))))
  iteration=1
  while(change_likelihood[1] > Tolerance & change_likelihood[2] > Tolerance & iteration<=maxstep){
  
  #SECCION A
  for(i in 1:n_states){
    if(k!=1){
    R=returns
    for(j in 1:nrow(returns)){
      R[j,]=returns[j,]-mu[,i]
    }
        B[,i] = exp(-.5*apply((as.matrix(R)%*%solve(as.matrix(sigma[[i]])))*as.matrix(R), 1, function(x)sum(x)))/((2*pi)^(k/2)*sqrt(abs(det(as.matrix(sigma[[i]])))))
      }else{
        B[,i] = exp(-.5*((returns-mu[i,])*sigma[i,]^(-1)*(returns-mu[i,])))/(sqrt(2*pi*sigma[i,]))
      }
    }
    
  
  # SECCION B
  forward[1,]=p*B[1,]
  forward[1,] = forward[1,]/sum(forward[1,])
  
  for (t in 2:L){
   aux=c(rep(0,n_states))
   for(i in 1:n_states){
    aux[i] =  aux[i] + sum(forward[t-1,]*A[,i])
   }
   forward[t,]=aux*B[t,]
   forward[t,] = forward[t,]/sum(forward[t,])
  }
  
  #SECCION C
  
  backward[L,]=B[L,]
  backward[L,]=backward[L,]/sum(backward[L,])
  
  t=L-1
  while (t>=1){
      aux=c(rep(0,n_states))
      for(i in 1:n_states){
        aux[i] =  aux[i] + sum(A[,i]*backward[t+1,])
      }
  backward[t,]=aux*B[t+1,]
  backward[t,]=backward[t,]/sum(backward[t,])
  t=t-1
  }               
                 
  #SECCION D
  for (t in 1:L){
    smoothed[t,]= forward[t,]*backward[t,]
    smoothed[t,]= smoothed[t,]/sum(smoothed[t,])
  }
  
  t=1
  while(t<=L-1){
  aux = diag(0,n_states,n_states)
  
  for(j in 1:n_states){
    aux[j,]=as.matrix(forward[t,j]*backward[t+1,]*B[t+1,])
  } 
  
  xi[[t]] = A*aux
  xi[[t]]=xi[[t]]/sum(xi[[t]])
  t=t+1
  }
  
  #SECCION E
  p=smoothed[1,]
  
  exp_num_transitions=data.frame(c(rep(0,n_states)),c(rep(0,n_states)))
  for(i in 1:n_states){
    for(j in 1:n_states){
     exp_num_transitions[i,j]=sum(sapply(lapply(xi, function(x)sum(x[i,j])), function(x)(sum(x))))  
    }
  }
  
  for(i in 1:n_states){
    if(sum(sapply(lapply(xi, function(x)sum(x[i,])), function(x)(sum(x))))!=0){
    A[i,] = exp_num_transitions[i,]/sum(sapply(lapply(xi, function(x)sum(x[i,])), function(x)(sum(x))))
    }else{A[i,]=c(rep(0,n_states))}
  if(k!=1){
    for(j in 1:k){
      if(sum(smoothed[,i])!=0){mu[j,i]=sum(smoothed[,i]*returns[,j])/sum(smoothed[,i])}else{mu[i,j]=0}
    }
  }else{
      if(sum(smoothed[,i])!=0){
        mu[i,]=sum(smoothed[,i]*returns)/sum(smoothed[,i])
        sigma[i,]=apply((smoothed[,i]*(returns-mu[i,])^2),2, function(x)sum(x)/sum(smoothed[,i]))
      }
  }
  
    #TOLERANCE
    if(k!=1){
      R=returns
      for(j in 1:nrow(returns)){
        R[j,]=returns[j,]-mu[,i]
      }
    likelihood[iteration, i] = log(sum(smoothed[,i]/sum(smoothed[,i])*exp(-.5*apply((as.matrix(R)%*%solve(as.matrix(sigma[[i]])))*as.matrix(R), 1, function(x)sum(x)))/((2*pi)^(k/2)*sqrt(abs(det(as.matrix(sigma[[i]])))))))
    }else{
    likelihood[iteration, i] = log(sum(smoothed[,i]/sum(smoothed[,i])*exp(-.5*((returns-mu[i,])*sigma^(-1)*(returns-mu[i,])))/((2*pi)^(k/2)*sqrt(cov(as.data.frame(returns))))))
    }
  }
  if(iteration>1){change_likelihood= abs(likelihood[iteration,]-likelihood[iteration-1,])}
  iteration=iteration+1
  }

return(list(mu=mu, sigma=sigma, p = p, A=A, forward=forward, backward=backward, smoothed=smoothed, exp_num_transitions))

}







