BaumWelch_ANN = function(returns, w, p, n_states=2, Tolerance=7*10^{-2}){
  
  change_likelihood=c(rep(Inf, n_states))
  likelihood=data.frame()
  
  returns=as.data.frame(returns)
  w = as.data.frame(w)
  A=data.frame(rep(1/n_states,n_states))
  A[1:n_states]=rep(1/n_states,n_states)
  #A=data.frame(c(0.9,0.1),c(0.1,0.9))
  p=rep(1/n_neurons, n_neurons)
  
  n_inputs=ncol(returns)
  L=nrow(returns)
  m = matrix(0, ncol = n_neurons, nrow = L)
  B=list(as.data.frame(m), as.data.frame(m))
  B[1:n_states]=list(data.frame(c(rep(0,L))))
  #B[1:n_states]=c(rep(0,L))
  forward=B
  backward=B
  smoothed=B
  xi=vector("list", L-1)
  xi[1:L-1]=list(data.frame(c(rep(0,n_states)),c(rep(0,n_states))))
  iteration=1
  
  #Generate Poisson Train
  n_trains=n_inputs*n_states
  high_freq=n_trains/2
  low_freq=n_trains/10
  df=1/n_trains
  m = matrix(ncol = L, nrow = n_inputs)
  for(i in 1:L){
    m[,i] = as.double(runif(n_trains)<ifelse(i %in% stdp[[1]], high_freq*df, low_freq*df))
  }
  
  
  while(change_likelihood[1] > Tolerance & change_likelihood[2] > Tolerance){
    
    #SECCION A
    for(i in 1:n_states){
      for(j in 1:L){
        B[[i]][j,]=exp(apply(w_new[[i]][seq((i-1)*n_inputs+1:i*n_inputs),]*as.double(returns[j,]),2,function(x)sum(x))+p[i,])
      }
      R[[i]][j,] = B[[i]][j,]/sum(B[[i]][j,])
    }
    
    rbind(apply(B[[1]], 1, function(x)max(x)), apply(B[[2]], 1, function(x)max(x)))
    

    #Repensarlo, hay que almacenar 5*20*730 datos, no sé bien con qué formato --- epsp depende también de la entrada y, no solo de t.
    w_change = function(w,L, n_neurons, n_states, stdp){
      aux=list(data.frame())
      for(i_1 in 1:L){
        aux[[i_1]]=w
        if(i_1 %in% stdp[[1]]){poisson_train=m[,i_1]}else{poisson_train=rev(m[,i_1])}
          for(i_2 in 1:n_trains){
            if(poisson_train[i_2]==1){aux[[i_1]][i_2,]=exp(-aux[[i_1]][i_2,]+1)-1}else{aux[[i_1]][i_2,]=rep(-1,n_neurons)}
          }
        }
      return(aux)
    }
    w_new=list(data.frame())
    for(i in 1:L){
      w_new[[i]]=w+0.1*w_change[[i]]
    }
    delta=w_change(w, L, n_neurons, n_states, stdp)
    
    
   
    
    
    
    
    
    
    
    
    #SECCION E
    
    p[1]=sum(R[,1])/length(R[,1])
    p[2]=sum(R[,2])/length(R[,2])
    # 
    # exp_num_transitions=data.frame(c(rep(0,n_states)),c(rep(0,n_states)))
    # for(i in 1:n_states){
    #   for(j in 1:n_states){
    #     exp_num_transitions[i,j]=sum(sapply(lapply(xi, function(x)sum(x[i,j])), function(x)(sum(x))))  
    #   }
    # }
    
    for(i in 1:n_states){
      if(k!=1){
        for(j in 1:k){
          #Mirar
          w[j,i]=sum(R[,i]*returns[,j])/sum(R[,i])
          
          #sigma[i,]=sqrt(sum(smoothed[,i]*(returns-as.data.frame(w)[i,])^2)/sum(smoothed[,i]))
          #sigma[i,]=sqrt(apply((smoothed[,i]%*%(returns-as.data.frame(w)[i,])^2),1, function(x) sum(x))/sum(smoothed[,i]))
        }
      }else{
        w[i,]=sum(B[,i]*returns)/sum(B[,i])
        # sigma[i,]=apply((smoothed[,i]*(returns-w[i,])^2),2, function(x)sum(x)/sum(smoothed[,i]))
      }
      #TOLERANCE
      # likelihood[iteration, i] = log(sum(smoothed[,i]/sum(smoothed[,i])*exp(t(w[,i]%*%t(returns)))))
      
    }
    if(iteration>1){change_likelihood= abs(likelihood[iteration,]-likelihood[iteration-1,])}
    iteration=iteration+1
  }
  
  return(list(w=w, p = p, R=R))
  
}
     
     # for t=1:T
     #     B(t,1) =  exp(-.5*((returns(t)-w(1))/sigma(1)).^2)/(sqrt(2*pi)*sigma(1));
     #     B(t,2) =  exp(-.5*((returns(t)-w(2))/sigma(2)).^2)/(sqrt(2*pi)*sigma(2));
     # end

    
    
    # SSECCION B
    
    # forward[1,]=p*B[1,]
    # forward[1,] = forward[1,]/sum(forward[1,])
    # 
    # for (t in 2:L){
    #   aux=c(rep(0,n_states))
    #   for(i in 1:n_states){
    #     aux[i] =  aux[i] + sum(forward[t-1,]*A[,i])
    #   }
    #   forward[t,]=aux*B[t,]
    #   forward[t,] = forward[t,]/sum(forward[t,])
    # }
    # 
    # #SECCION C
    # 
    # backward[L,]=B[L,]
    # backward[L,]=backward[L,]/sum(backward[L,])
    # 
    # t=L-1
    # while (t>=1){
    #   aux=c(rep(0,n_states))
    #   for(i in 1:n_states){
    #     aux[i] =  aux[i] + sum(A[,i]*backward[t+1,])
    #   }
    #   backward[t,]=aux*B[t+1,]
    #   backward[t,]=backward[t,]/sum(backward[t,])
    #   t=t-1
    # }
    # 
    # #SECCION D
    # 
    # 
    # for (t in 1:L){
    #   smoothed[t,]= forward[t,]*backward[t,]
    #   smoothed[t,]= smoothed[t,]/sum(smoothed[t,])
    # }
    # 
    # t=1
    # while(t<=L-1){
    #   aux = diag(0,n_states,n_states)
    # 
    #   for(j in 1:n_states){
    #     aux[j,]=as.matrix(forward[t,j]*backward[t+1,]*B[t+1,])
    #   }
    # 
    #   xi[[t]] = A*aux
    #   xi[[t]]=xi[[t]]/sum(xi[[t]])
    #   t=t+1
    # }

