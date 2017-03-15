BaumWelch = function(returns, w, p, n_states=2, Tolerance=7*10^{-2}){
  
  change_likelihood=c(rep(Inf, n_states))
  likelihood=data.frame()
  
  returns=as.data.frame(returns)
  w = as.data.frame(w)
  A=data.frame(rep(1/n_states,n_states))
  A[1:n_states]=rep(1/n_states,n_states)
  #A=data.frame(c(0.9,0.1),c(0.1,0.9))
  p=rep(1/n_states,n_states)
  k=ncol(returns)
  L=nrow(returns)
  
  B=data.frame(c(rep(0,L)))
  B[1:n_states]=c(rep(0,L))
  forward=B
  backward=B
  smoothed=B
  xi=vector("list", L-1)
  xi[1:L-1]=list(data.frame(c(rep(0,n_states)),c(rep(0,n_states))))
  iteration=1
  while(change_likelihood[1] > Tolerance & change_likelihood[2] > Tolerance){
    
    #SECCION A
    
    
    for(i in 1:n_states){
        B[,i] = exp(-.5*apply(as.matrix(w)*as.matrix(retornos), 1, function(x)sum(x)))/((2*pi)^(k/2))
    }
    
    # % for t=1:T
    # %     B(t,1) =  exp(-.5*((returns(t)-w(1))/sigma(1)).^2)/(sqrt(2*pi)*sigma(1));
    # %     B(t,2) =  exp(-.5*((returns(t)-w(2))/sigma(2)).^2)/(sqrt(2*pi)*sigma(2));
    # % end
    # 
    
    
    # SSECCION B
    
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
          #Mirar
          if(sum(smoothed[,i])!=0){w[j,i]=sum(smoothed[,i]*returns[,j])/sum(smoothed[,i])}else{w[i,j]=0}
          
          #sigma[i,]=sqrt(sum(smoothed[,i]*(returns-as.data.frame(w)[i,])^2)/sum(smoothed[,i]))
          #sigma[i,]=sqrt(apply((smoothed[,i]%*%(returns-as.data.frame(w)[i,])^2),1, function(x) sum(x))/sum(smoothed[,i]))
        }
      }else{
        w[i,]=sum(smoothed[,i]*returns)/sum(smoothed[,i])
        sigma[i,]=apply((smoothed[,i]*(returns-w[i,])^2),2, function(x)sum(x)/sum(smoothed[,i]))
      }
      
      
      
      #TOLERANCE
      likelihood[iteration, i] = log(sum(smoothed[,i]/sum(smoothed[,i])*exp(-.5*apply(as.matrix(w)*as.matrix(retornos), 1, function(x)sum(x)))/((2*pi)^(k/2))))

    }
    if(iteration>1){change_likelihood= abs(likelihood[iteration,]-likelihood[iteration-1,])}
    iteration=iteration+1
  }
  
  return(list(w=w, p = p, A=A, forward=forward, backward=backward, smoothed=smoothed, esp_num_transitions=exp_num_transitions))
  
}