#### Ejemplo HMM
source(paste(getwd(), "/","Baum-Welch.R", sep=""))

df=read.csv(paste0(getwd(), "/Ibex35.csv"), header = T)[1:200,]
returns=double()
for (i in 1:199){returns[i]=(df$APERTURA[i+1]-df$APERTURA[i])/df$APERTURA[i]}
mu=c(mean(returns[which(returns>0)]), mean(returns[which(returns<0)]))
sigma=c(var(returns[which(returns>0)]), var(returns[which(returns<0)]))
var=BaumWelch(returns = returns, mu=mu, sigma=sigma, n_states = 2)
smoothed=var$smoothed[,1]
ggplot(df[-1,], aes(x=FECHA, y=APERTURA, group=1, colour=smoothed))+geom_line(size=0.2) + 
  scale_x_discrete(limits=df[-1,]$FECHA,breaks=df[-1,]$FECHA[seq(1,191,10)]) +scale_color_gradient(low="red",high="green") +
  metri_style



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

