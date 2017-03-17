#libraries
library(ggplot2)
library(dplyr)
library(neuralnet)


#### Ejemplo HMM
source(paste(getwd(), "/","Baum-Welch.R", sep=""))


df  = read.delim("DesigualB2C_ANN.txt")

df_retornos=df[-nrow(df),]
for(i in 1:(nrow(df_retornos))){
  for(j in 2:(ncol(df_retornos))){
    df_retornos[i, j]=(df[i+1,j]-df[i,j])/df[i,j]
  }
}

retornos_sess=df_retornos[,2]
retornos_users=df_retornos[,3]
retornos_transactions=df_retornos[,6]
retornos_newusers=df_retornos[,5]
retornos_bounces=df_retornos[,4]

var_sess=BaumWelch(returns = retornos_sess, mu=c(median(retornos_sess[which(retornos_sess>0)]),0,median(retornos_sess[which(retornos_sess<0)])),
                   sigma= rep(var(retornos_sess),3), n_states = 3)
var_users=BaumWelch(returns = retornos_users, mu=c(median(retornos_users[which(retornos_users>0)]),0,median(retornos_users[which(retornos_users<0)])), 
                    sigma= rep(var(retornos_users),3), n_states = 3)
var_transactions=BaumWelch(returns = retornos_transactions, mu=c(median(retornos_transactions[which(retornos_transactions>0)]),0,median(retornos_transactions[which(retornos_transactions<0)])),
                           sigma= rep(var(retornos_transactions),3), n_states = 3)
var_newusers=BaumWelch(returns = retornos_newusers, mu=c(median(retornos_newusers[which(retornos_newusers>0)]),0,median(retornos_newusers[which(retornos_newusers<0)])),
                           sigma= rep(var(retornos_newusers),3), n_states = 3)
var_bounces=BaumWelch(returns = retornos_bounces, mu=c(median(retornos_bounces[which(retornos_bounces>0)]),0,median(retornos_bounces[which(retornos_bounces<0)])),
                       sigma= rep(var(retornos_bounces),3), n_states = 3)

smoothed_users = var_users$smoothed
names(smoothed_users)=c("1","2","3")
A=double()
for (i in 1:nrow(smoothed_users)) {
  A[i]=as.double(names(smoothed_users)[which(t(smoothed_users[i,])==max(t(smoothed_users[i,])))])
}

smoothed_sessions= var_sess$smoothed
names(smoothed_sessions)=c("1","2","3")
B=double()
for (i in 1:nrow(smoothed_sessions)) {
  B[i]=as.double(names(smoothed_sessions)[which(t(smoothed_sessions[i,])==max(t(smoothed_sessions[i,])))])
}


smoothed_transactions= var_transactions$smoothed
names(smoothed_transactions)=c("1","2","3")
C=double()
for (i in 1:nrow(smoothed_transactions)) {
  C[i]=as.double(names(smoothed_transactions)[which(t(smoothed_transactions[i,])==max(t(smoothed_transactions[i,])))])
}


smoothed_newusers= var_newusers$smoothed
names(smoothed_newusers)=c("1","2","3")
D=double()
for (i in 1:nrow(smoothed_newusers)) {
  D[i]=as.double(names(smoothed_newusers)[which(t(smoothed_newusers[i,])==max(t(smoothed_newusers[i,])))])
}


smoothed_bounces= var_bounces$smoothed
names(smoothed_bounces)=c("1","2","3")
E=double()
for (i in 1:nrow(smoothed_bounces)) {
  E[i]=as.double(names(smoothed_bounces)[which(t(smoothed_bounces[i,])==max(t(smoothed_bounces[i,])))])
}


df_retornos$smoothed=as.double(C)
df$smoothed=as.double(c(2,C))
r=ggplot(df_retornos, aes(x=date, y=transactions, group=1, colour= factor(smoothed)))+geom_line(size=0.2) + geom_smooth()+theme(axis.text.x = element_text(angle=30))+
  scale_x_discrete(limits=df$date,breaks=df$date[seq(1,366,10)])+scale_colour_manual(labels=c("1","2","3"), values=c("green", "yellow", "red"))


df_retornos$smoothed=as.double(A)
df$smoothed=as.double(c(2,A))
p=ggplot(df, aes(x=date, y=sessions, group=1, colour= factor(smoothed)))+geom_line(size=0.2) + geom_smooth()+theme(axis.text.x = element_text(angle=30))+
  scale_x_discrete(limits=df$date,breaks=df$date[seq(1,366,10)])+scale_colour_manual(labels=c("1","2","3"), values=c("green", "yellow", "red"))



df_retornos$smoothed=as.double(B)
df$smoothed=as.double(c(2,B))
q=ggplot(df, aes(x=date, y=users, group=1, colour= factor(smoothed)))+geom_line(size=0.2) + geom_smooth()+theme(axis.text.x = element_text(angle=30))+
  scale_x_discrete(limits=df$date,breaks=df$date[seq(1,366,10)])+scale_colour_manual(labels=c("1","2","3"), values=c("green", "yellow", "red"))

multiplot(p,q,r)
