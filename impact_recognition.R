#libraries
library(ggplot2)
library(dplyr)
library(neuralnet)


source(paste(getwd(), "/","Baum-Welch.R", sep=""))
source(paste(getwd(), "/","GHMM-ANN.R", sep=""))

df  = read.delim("DesigualB2C_ANN.txt")
n_inputs=1 # múltiplo de n_var
n_neurons=10
n_states=2

df_1=df[1:365,]
df_2=df[c(-1:-365),]

df_retornos=df[-nrow(df),]
for(i in 1:(nrow(df_retornos))){
  for(j in 2:(ncol(df_retornos))){
    df_retornos[i, j]=(df[i+1,j]-df[i,j])/df[i,j]
  }
}

net_1=neuralnet( ~ round(sessions), df, 
              hidden = 10, threshold = 10 , stepmax = 10000)


#### ANN
  
  df_retornos$smoothed=0
  up=which(df_retornos$transactions>0)
  St_4 = df_retornos[which(df_retornos$transactions>0),]
  df_retornos$smoothed=1
  down=which(df_retornos$transactions<0)
  St_3 = df_retornos[which(df_retornos$transactions<0),]
  stdp=list(up, down)
  
  # St_2 = df_retornos
  # df_retornos$smoothed=0
  # df_retornos$smoothed[which(df_retornos$sessions<(-0.1))]=1
  # St_1 = df_retornos
  
  
  # net_st_1=neuralnet(smoothed ~ sessions + users + bounceRate, St_1, 
  #                        hidden = 4, threshold = 0.001, stepmax = 100000)
  # net_st_2=neuralnet(smoothed ~ sessions + users + bounceRate + transactions, St_2, 
  #                    hidden = 4, threshold = 0.001, stepmax = 100000)
  net_st_3=neuralnet(smoothed ~ sessions + users + newusers + pageviews + avgTimeOnPage, St_3, 
                     hidden = 10, threshold = 0.01, stepmax = 1000000)
  net_st_4=neuralnet(smoothed ~ sessions + users + newusers + pageviews + avgTimeOnPage, St_4, 
                     hidden = 10, threshold = 0.01, stepmax = 100000)
  
  test=data.frame(St_3$smoothed, round(as.data.frame(net_st_3$net.result)[,1]))
  
  w=rbind(as.data.frame(net_st_3$weights[[1]][[1]])[-1,] , as.data.frame(net_st_4$weights[[1]][[1]])[-1,])
  p=rbind(as.data.frame(net_st_3$weights[[1]][[1]])[1,] , as.data.frame(net_st_4$weights[[1]][[1]])[1,])
  
  
  
  
  #################################################
  w=matrix(ncol = n_inputs, nrow = n_neurons*n_states)
  for(i in 1:nrow(w)){w[i,]=runif(ncol(w))}
  p=matrix(nrow = n_neurons*n_states, ncol = 1)
  for(i in 1:nrow(p)){p[i,]=runif(ncol(p))}
  source(paste(getwd(), "/","GHMM-ANN.R", sep=""))
  
  results=BaumWelch_ANN(returns = df_retornos[,9], w=w,p=p,n_states = 2, n_neurons=10, stdp=stdp)
  
  net_st_4$weights[[1]][[1]]=data.frame(c(p[2], w[,2] ))
  net_st_3$weights[[1]][[1]]=data.frame(c(p[1], w[,1] ))
  
  compute_w <-function (x, covariate, rep = 1) 
  {
    nn <- x
    linear.output <- nn$linear.output
    weights <- nn$weights[[rep]]
    length.weights <- length(weights)
    covariate <- as.matrix(cbind(1, covariate))
    act.fct <- function(x){x}
    neurons <- list(covariate)
    if (length.weights > 1) 
      for (i in 1:(length.weights - 1)) {
        temp <- neurons[[i]] %*% weights[[i]]
        act.temp <- act.fct(temp)
        neurons[[i + 1]] <- cbind(1, act.temp)
      }
    temp <- neurons[[length.weights]] %*% weights[[length.weights]]
    if (linear.output) 
      net.result <- temp
    else net.result <- act.fct(temp)
    list(neurons = neurons, net.result = net.result)
  }
  
  C=compute_w(net_st_3, df_retornos[,c(2,3,4,6)])$net.result
  D=compute_w(net_st_4, df_retornos[,c(2,3,4,6)])$net.result
  
  names(smoothed)=c("1","2")
  D=double()
  for (i in 1:nrow(smoothed)) {
    D[i]=as.double(names(smoothed)[which(t(smoothed[i,])==max(t(smoothed[i,])))])
  }
  
  df_retornos$smoothed=B
  df$smoothed=c(1,B)
  ggplot(df, aes(x=date, y=sessions, group=1, colour= factor(smoothed)))+geom_point(size=3)+geom_line(color="black") + geom_smooth()+theme(axis.text.x = element_text(angle=30))+
    scale_x_discrete(limits=df$date,breaks=df$date[seq(1,366,10)])+scale_colour_manual(labels=c("1","2","3"), values=c("red", "green", "yellow"))
  
  
  ggplot(df_retornos, aes(x=date, y=sessions, group=1, colour= factor(smoothed)))+geom_line(size=0.2) + geom_smooth()+theme(axis.text.x = element_text(angle=30))+
scale_x_discrete(limits=df$date,breaks=df$date[seq(1,366,10)])+scale_colour_manual(labels=c("1","2","3"), values=c("red", "green", "yellow"))
  
  
  
  
  # 
  # retornos_sess=df_retornos[,2]
  # retornos_users=df_retornos[,3]
  # retornos_transactions=df_retornos[,6]
  # retornos_newusers=df_retornos[,5]
  # retornos_bounces=df_retornos[,4]
  # 
  # var_sess=BaumWelch(returns = retornos_sess, #mu=c(median(retornos_sess[which(retornos_sess>0)]),mean(retornos_sess[which(retornos_sess>0)]), mean(retornos_sess[which(retornos_sess<0)]), median(retornos_sess[which(retornos_sess<0)])),
  #                    mu=c(-0.2,0.2),
  #                    sigma= rep(var(retornos_sess),2), n_states = 2)
  # var_users=BaumWelch(returns = retornos_users, #mu=c(median(retornos_users[which(retornos_users>0)]),0,median(retornos_users[which(retornos_users<0)])), 
  #                     mu=c(-0.4,-0.2,0.2,0.3),                   
  #                     sigma= rep(var(retornos_users),4), n_states = 4)
  # var_transactions=BaumWelch(returns = retornos_transactions, #mu=c(median(retornos_transactions[which(retornos_transactions>0)]),0,median(retornos_transactions[which(retornos_transactions<0)])),
  #                            mu=c(-0.4,-0.2,0.2,0.3),
  #                            sigma= rep(var(retornos_transactions),4), n_states = 4)
  # var_newusers=BaumWelch(returns = retornos_newusers, #mu=c(median(retornos_newusers[which(retornos_newusers>0)]),0,median(retornos_newusers[which(retornos_newusers<0)])),
  #                        mu=c(-0.4,-0.2,0.2,0.3),
  #                        sigma= rep(var(retornos_newusers),4), n_states = 4)
  # var_bounces=BaumWelch(returns = retornos_bounces, #mu=c(median(retornos_bounces[which(retornos_bounces>0)]),0,median(retornos_bounces[which(retornos_bounces<0)])),
  #                       mu=c(-0.4,-0.2,0.2,0.3),
  #                       sigma= rep(var(retornos_bounces),4), n_states = 4)
  # 
  # smoothed_users = var_users$smoothed
  # names(smoothed_users)=c("1","2","3","4")
  # A=double()
  # for (i in 1:nrow(smoothed_users)) {
  #   A[i]=as.double(names(smoothed_users)[which(t(smoothed_users[i,])==max(t(smoothed_users[i,])))])
  # }
  # 
  # smoothed_sessions= smoothed
  # names(smoothed_sessions)=c("1","2")
  # B=double()
  # for (i in 1:nrow(smoothed_sessions)) {
  #   B[i]=as.double(names(smoothed_sessions)[which(t(smoothed_sessions[i,])==max(t(smoothed_sessions[i,])))])
  # }
  # 
  # 
  # smoothed_transactions= var_transactions$smoothed
  # names(smoothed_transactions)=c("1","2","3","4")
  # C=double()
  # for (i in 1:nrow(smoothed_transactions)) {
  #   C[i]=as.double(names(smoothed_transactions)[which(t(smoothed_transactions[i,])==max(t(smoothed_transactions[i,])))])
  # }
  # 
  # 
  # smoothed_newusers= var_newusers$smoothed
  # names(smoothed_newusers)=c("1","2","3","4")
  # D=double()
  # for (i in 1:nrow(smoothed_newusers)) {
  #   D[i]=as.double(names(smoothed_newusers)[which(t(smoothed_newusers[i,])==max(t(smoothed_newusers[i,])))])
  # }
  # 
  # 
  # smoothed_bounces= var_bounces$smoothed
  # names(smoothed_bounces)=c("1","2","3", "4")
  # E=double()
  # for (i in 1:nrow(smoothed_bounces)) {
  #   E[i]=as.double(names(smoothed_bounces)[which(t(smoothed_bounces[i,])==max(t(smoothed_bounces[i,])))])
  # }
  # 
  # 
  # df_retornos$smoothed=as.double(E)
  # df$smoothed=as.double(c(2,E))
  # ggplot(df_retornos, aes(x=date, y=sessions, group=1, colour= factor(round(smoothed,0))))+geom_line(size=0.2) + geom_smooth()+theme(axis.text.x = element_text(angle=30))+
  #   scale_x_discrete(limits=df$date,breaks=df$date[seq(1,366,10)])+scale_colour_discrete(h=c(0,180), l=60, c=180)
  # 
  # 
  # df_retornos$smoothed=as.double(A)
  # df$smoothed=as.double(c(2,A))
  # p=ggplot(df, aes(x=date, y=sessions, group=1, colour= factor(smoothed)))+geom_line(size=0.2) + geom_smooth()+theme(axis.text.x = element_text(angle=30))+
  #   scale_x_discrete(limits=df$date,breaks=df$date[seq(1,400,10)])+scale_colour_manual(labels=c("1","2","3"), values=c("green", "yellow", "red"))
  # 
  # 
  # 
  # df_retornos$smoothed=as.double(B)
  # df$smoothed=as.double(c(2,B))
  # p=ggplot(df_retornos, aes(x=date, y=sessions, group=1, colour= factor(smoothed)))+geom_point(size=3)+geom_line(color="black") + geom_smooth()+theme(axis.text.x = element_text(angle=30))+
  #   scale_x_discrete(limits=df$date,breaks=df$date[seq(1,366,10)])+scale_colour_discrete(h=c(0,180), l=50)
  # q=ggplot(df[-1,], aes(x=date, y=sessions, group=1, colour= factor(smoothed)))+geom_point(size=3)+geom_line(color="black") + geom_smooth()+theme(axis.text.x = element_text(angle=30))+
  #   scale_x_discrete(limits=df$date,breaks=df$date[seq(1,366,10)])+scale_colour_discrete(h=c(0,180), l=50)
  # 
  # multiplot(p,q)
  # 
  # 
  # ggplot(data=df_retornos, aes(x=sessions))+
  #   geom_density(alpha=.2, fill="blue")+
  #   stat_function(fun=dnorm, args = list(mean=var_sess$mu[1,], sd=sqrt(var_sess$sigma[1,])), colour="green")+
  #   stat_function(fun=dnorm, args = list(mean=var_sess$mu[2,], sd=sqrt(var_sess$sigma[2,])), colour="red")+
  #   stat_function(fun=dnorm, args = list(mean=var_sess$mu[4,], sd=sqrt(var_sess$sigma[3,])), colour="orange")+
  #   stat_function(fun=dnorm, args = list(mean=var_sess$mu[3,], sd=sqrt(var_sess$sigma[4,])), colour="yellow")+
  #   scale_x_continuous(breaks = round(seq(min(df_retornos$sessions), max(df_retornos$sessions), by = 0.1),1)) 
  # geom_density(aes(x=users), fill="green", alpha=0.2)+
  #   geom_density(aes(x=bounces), fill="red", alpha=.2)+
  #   geom_density(aes(x=newusers), fill="yellow",alpha=0.2)+  
  
  