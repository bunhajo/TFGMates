####
source(paste(getwd(), "/","Baum-Welch.R", sep=""))

mu=c(mean(returns$returns[which(returns$returns>0)]), mean(returns$returns[which(returns$returns<0)]))
sigma=c(var(returns), var(returns))
BaumWelch(returns = data, mu=mu, sigma=sigma, n_states = 2)

#### 
aux=aux %>% group_by(table1_visitnumber, table1_hits_eventInfo_eventAction, table1_hits_product_v2ProductCategory) %>% mutate(obs= n())
aux=unique(aux)
mu=matrix(1,ncol =  30, nrow=15)
for(i in 1:15){mu[i,2*i-1]=10
mu[i,2*i]=5}

for(i in min(aux$table1_visitnumber):max(aux$table1_visitnumber)){
  aux2=aux[which(aux$table1_visitnumber==i),]  
  for(j in 1:length(states)){
    observ[i,2*j-1]=sum(aux2$obs[which((aux2$table1_hits_eventInfo_eventAction=="Product Detail" | aux2$table1_hits_eventInfo_eventAction=="Products Impressions On List") & aux2$table1_hits_product_v2ProductCategory==states[j])])
    observ[i,2*j]=sum(aux2$obs[which(aux2$table1_hits_eventInfo_eventAction=="Add to cart" & aux2$table1_hits_product_v2ProductCategory==states[j])])
  }
}
####

ggplot(aux, aes(x=FECHA, y=APERTURA, group=1, colour=smooth))+geom_line() + scale_x_discrete(limits=aux$FECHA)
