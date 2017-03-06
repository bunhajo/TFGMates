#ALGORITMO BAUM-WELCH


change_likelihood=data.frame(c(Inf, Inf))



Baum-Welch = function(returns, mu, sigma, p, n_states=2)


returns=as.data.frame(returns)
sigma=c(cov(retornos),cov(retornos))
A=data.frame(c(0.5,0.5),c(0.5,0.5))
p=c(0.5,0.5)
k=ncol(returns)
L=nrow(returns)

B=data.frame(c(rep(0,L)))
B[1:n_states]=c(rep(0,L))
forward=B
backward=B
smoothed=B
xi=vector("list", L-1)
xi[1:L-1]=list(data.frame(c(rep(0,n_states)),c(rep(0,n_states))))
while(change_likelihood[iteration-1,1] > Tolerance && change_likelihood[iteration-1,2] > Tolerance){

#SECCION A

  #Arreglar
for(i in 1:n_states){
  if(k!=1){
    B[,i] = exp(-.5*apply(t((returns-mu[,i]))%*%sigma^(-1)*t((returns-mu[,i])), 1, function(x)sum(x)))/((2*pi)^(k/2)*sqrt(abs(det(sigma))))
  }else{
    B[,i] = exp(-.5*((returns-mu[i,])*sigma^(-1)*(returns-mu[i,])))/((2*pi)^(k/2)*sqrt(cov(as.data.frame(returns))))
  }
}
# % for t=1:T
# %     B(t,1) =  exp(-.5*((returns(t)-mu(1))/sigma(1)).^2)/(sqrt(2*pi)*sigma(1));
# %     B(t,2) =  exp(-.5*((returns(t)-mu(2))/sigma(2)).^2)/(sqrt(2*pi)*sigma(2));
# % end
# 
  
  
# SSECCION B

forward[1,]=p*B[1,]
forward[1,] = forward[1,]/sum(forward[1,])

for (t in 2:L){
 aux=c(rep(0,n_states))
 for(i in 1:length(A)){
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
    for(i in 1:length(A)){
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
aux = data.frame(c(rep(0,n_states)),c(rep(0,n_states)))

for(j in 1:n_states){
  aux[j,]=forward[t,j]*backward[t+1,]*B[t+1,]
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
A[i,] = exp_num_transitions[i,]/sum(sapply(lapply(xi, function(x)sum(x[i,])), function(x)(sum(x))))
if(k!=1){
  for(j in 1:k)
  mu[i,k]=sum(smoothed[,i]*returns[,j])/sum(smoothed[,i])
  #sigma[i,]=sqrt(sum(smoothed[,i]*(returns-as.data.frame(mu)[i,])^2)/sum(smoothed[,i]))
  #sigma[i,]=sqrt(apply((smoothed[,i]%*%(returns-as.data.frame(mu)[i,])^2),1, function(x) sum(x))/sum(smoothed[,i]))
  }else{
  mu[i]=sum(smoothed[,i]*returns)/sum(smoothed[,i])
  sigma[i]=sqrt(apply((smoothed[,i]*(returns-as.data.frame(mu)[i,])^2),1, function(x) sum(x))/sum(smoothed[,i]))
  }
}


#TOLERANCE

# for i=1:2
# A(i,:)=exp_num_transitions(i,:)/sum(sum(xi(i,:,:),2),3);%(2)
#        sigma(i)=sqrt(sum(smoothed(:,i).*(returns-mu(i)).^2)/sum(smoothed(:,i)));%(3)
#        end
#        
#        %TOLERANCE
#        %%
#          
#          likelihood(iteration,1)=log(...
#                                      (((exp(-.5*((returns(t)-mu(1)*ones(T,1))/sigma(1)).^2)/(sqrt(2*pi)*sigma(1)))'*smoothed(:,1))/sum(smoothed(:,1)))...
#                                      );
#                                        likelihood(iteration,2)=log(...
#                                        (((exp(-.5*((returns(t)-mu(2)*ones(T,1))/sigma(2)).^2)/(sqrt(2*pi)*sigma(2)))'*smoothed(:,2))/sum(smoothed(:,2)))...
#                                        );
#        
#        change_likelihood(iteration,:)= abs(likelihood(iteration,:)-likelihood(iteration-1,:));
#        iteration=iteration+1;
#        
# }               
                              
                              
                              