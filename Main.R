#libraries
library(ggplot2)
library(dplyr)
library(neuralnet)


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
  scale_x_discrete(limits=df[-1,]$FECHA,breaks=df[-1,]$FECHA[seq(1,191,10)]) +scale_color_gradient(low="red",high="green")
  

#### Pruebas Desigual ####

#HMM

Prueb_HMM =  Prueb_HMM %>% group_by(table1_fullvisitorId) %>% 
             mutate(n_session=max(table1_visitnumber)-min(table1_visitnumber)+1)

max_sess=max(Prueb_HMM$n_session)

ID=unique(Prueb_HMM$table1_fullvisitorId[which(Prueb_HMM$n_session==max_sess)])
prova=Prueb_HMM[which(Prueb_HMM$table1_fullvisitorId==ID),]
prova=prova[-which(prova$table1_hits_product_v2ProductCategory %in% c("", "(not set)")),]
aux=prova %>% group_by(table1_visitnumber, table1_hits_eventInfo_eventAction, table1_hits_product_v2ProductCategory) %>% summarise(obs= n())
aux=unique(aux)
states=unique(aux$table1_hits_product_v2ProductCategory)
mu=matrix(1,ncol = 2*length(states), nrow=length(states))
for(i in 1:length(states)){mu[i,2*i-1]=5
               mu[i,2*i]=5}

observ=data.frame()
for(i in min(aux$table1_visitnumber):max(aux$table1_visitnumber)){
  aux2=aux[which(aux$table1_visitnumber==i),]  
  for(j in 1:length(states)){
    observ[i,2*j-1]=sum(aux2$obs[which((aux2$table1_hits_eventInfo_eventAction=="Product Detail" | aux2$table1_hits_eventInfo_eventAction=="Products Impressions On List") & aux2$table1_hits_product_v2ProductCategory==states[j])])
    observ[i,2*j]=sum(aux2$obs[which(aux2$table1_hits_eventInfo_eventAction=="Add to cart" & aux2$table1_hits_product_v2ProductCategory==states[j])])
  }
}
observ=observ[-which(apply(observ, 1, function(x)sum(x))==0),]
rem_col=as.double(which(apply(observ, 2, function(x)sum(x))==0))
mu=mu[,-rem_col]
observ=observ[,-rem_col]

var=BaumWelch(returns = observ, mu=a_mu, sigma=diag(1,ncol(observ)), n_states = length(states), Tolerance = 100)


#### Red Neuronal

Prueba_ANN  = read.delim(paste0(getwd(), "/Prueba_TFG.txt")) 
df=DesigualB2C_ANN
net=neuralnet(transactions ~ sessions + users + bounces + newusers, DesigualB2C_ANN[1:300,], 
              hidden = 300, threshold = 0.001, stepmax = 500000)
copy_net=net
test=data.frame(net$net.result, Prueba_ANN$transactions)


# Diario
df_retornos=df[-nrow(df),]
for(i in 1:(nrow(df_retornos))){
  for(j in 2:(ncol(df_retornos))){
    df_retornos[i, j]=(df[i+1,j]-df[i,j])/df[i,j]
  }
}


retornos=df_retornos[,c(2:5)]
# mu=data.frame(good=as.double(sapply(retornos[which(df_retornos$transactions>0),], function(x)median(x))),
#               meh= as.double(sapply(retornos[which(df_retornos$transactions>-0.05 & df_retornos$transactions<0.05),], function(x)median(x))),
#               bad= as.double(sapply(retornos[which(df_retornos$transactions<(-0)),], function(x)median(x))))

ggplot(data=retornos)+
  geom_density(aes(x=sessions), fill="blue", alpha=.2)+
  geom_density(aes(x=users), fill="green", alpha=0.2)+
  geom_density(aes(x=bounces), fill="red", alpha=.2)+
  geom_density(aes(x=newusers), fill="yellow",alpha=0.2)+
  stat_function(fun=dnorm, args = list(mean=0.07, sd=0.1))+
  scale_x_continuous(breaks = round(seq(min(retornos$sessions), max(retornos$sessions), by = 0.1),1)) 


mu=data.frame(good=rep(0.1,4),
              meh=rep(0,4),
              bad=rep(-.1,4))

sigma=list()
for(i in 1:4){sigma[[i]]=var(retornos)}

var=BaumWelch(returns = retornos, mu=mu, sigma= sigma, n_states = 2)
#smoothed=c(0,as.double(var$smoothed[,1]>var$smoothed[,2]))


smoothed = var$smoothed
names(smoothed)=c("1","2")
A=double()
for (i in 1:nrow(smoothed)) {
  A[i]=as.double(names(smoothed)[which(t(smoothed[i,])==max(t(smoothed[i,])))])
}


DesigualB2C_ANN$smoothed=as.double(c(2, A))

ggplot(DesigualB2C_ANN, aes(x=date, y=sessions, group=1, colour= factor(smoothed)))+geom_line(size=0.2) + 
  scale_x_discrete(limits=df$date,breaks=df$date[seq(1,366,10)])+scale_colour_manual(labels=c("1","2","3"), values=c("green", "yellow", "red"))



#Semanal
df_wretornos=df_w[-nrow(df_w),]
for(i in 1:(nrow(df_wretornos))){
  for(j in 2:(ncol(df_wretornos))){
    df_wretornos[i, j]=(df_w[i+1,j]-df_w[i,j])/df_w[i,j]
  }
}
wretornos=df_wretornos[,c(2:5)]

ggplot(data=wretornos)+
  geom_density(aes(x=w_sessions), fill="blue", alpha=.2)+
  geom_density(aes(x=w_users), fill="green", alpha=0.2)+
  geom_density(aes(x=w_bounces), fill="red", alpha=.2)+
  geom_density(aes(x=w_newusers), fill="yellow",alpha=0.2)+
  #stat_function(fun=dnorm, args = list(mean=0.07, sd=0.1))+
  scale_x_continuous(breaks = round(seq(min(retornos$sessions), max(retornos$sessions), by = 0.1),1)) 

mu=data.frame(good=c(mean(wretornos$w_sessions[which(wretornos$w_sessions>0)]),
                     mean(wretornos$w_users[which(wretornos$w_users>0)]),
                     mean(wretornos$w_bounces[which(wretornos$w_bounces>0)]),
                     mean(wretornos$w_newusers[which(wretornos$w_newusers>0)])),
              bad=c(mean(wretornos$w_sessions[which(wretornos$w_sessions<0)]),
                    mean(wretornos$w_users[which(wretornos$w_users<0)]),
                    mean(wretornos$w_bounces[which(wretornos$w_bounces<0)]),
                    mean(wretornos$w_newusers[which(wretornos$w_newusers<0)])))


# mu=data.frame(good=as.double(sapply(wretornos[which(df_wretornos$w_transactions>0),], function(x)median(x))),
#               meh= c(0,0,0,0),
#               bad= as.double(sapply(wretornos[which(df_wretornos$w_transactions<0),], function(x)median(x))))

sigma=list()
sigma[[1]]=var(wretornos[which(as.double(apply(wretornos, 1, function(x)mean(x)))>0),])
sigma[[2]]=var(wretornos[which(as.double(apply(wretornos, 1, function(x)mean(x)))<0),])

var=BaumWelch(returns = wretornos, mu=mu, sigma= sigma, n_states = 2, Tolerance = 0.001)

ggplot(prov1, aes(x=w_sessions))+
  geom_density(aes(x=w_sessions), fill="blue", alpha=.2)+
  geom_density(aes(x=w_users), fill="green", alpha=0.2)+
  geom_density(aes(x=w_bounces), fill="red", alpha=.2)+
  geom_density(aes(x=w_newusers), fill="yellow",alpha=0.2)+
  stat_function(fun=dnorm, args = list(mean=0.14, sd=0.1))+
  stat_function(fun=dnorm, args = list(mean=var$mu[2,2], sd=0.1))+
  stat_function(fun=dnorm, args = list(mean=var$mu[3,2], sd=0.1))+
  stat_function(fun=dnorm, args = list(mean=var$mu[4,2], sd=0.1))+
  stat_function(fun=dnorm, args = list(mean=mu[1,1], sd=1))+
  stat_function(fun=dnorm, args = list(mean=mu[2,1], sd=1))+
  stat_function(fun=dnorm, args = list(mean=mu[3,1], sd=1))+
  stat_function(fun=dnorm, args = list(mean=mu[4,1], sd=1))
  scale_x_continuous(breaks = round(seq(min(retornos$sessions), max(retornos$sessions), by = 0.1),1)) 
  
  smoothed = var$smoothed
  names(smoothed)=c("1","2","3")
  A=double()
  for (i in 1:nrow(smoothed)) {
    A[i]=as.double(names(smoothed)[which(t(smoothed[i,])==max(t(smoothed[i,])))])
  }
  
  
  df_w$smoothed=as.double(c(2, A))

  ggplot(df_w[,c(1:2, ncol(df_w))], aes(x=week, y=w_sessions,  group=1,colour= factor(smoothed)))+geom_line(size=0.2) + 
    scale_x_discrete(limits=df_w$week,breaks=df_w$week[seq(1,53,10)])+scale_colour_manual(labels=c("1","2","3"), values=c("green", "yellow", "red"))
  