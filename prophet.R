##Prophet
library(prophet)
library(dplyr)

df  = read.delim("DesigualB2C_ANN.txt")
df_1=df[1:365,]
df_2=df[c(-1:-365),]
names(df_1)[1:2]=c("ds", "y")

m = prophet(df_1)
future = make_future_dataframe(m, periods = 365)
forecast = predict(m, future)

forecast = full_join(forecast, df[,1:2])

ggplot(forecast, aes(x=ds,y=yhat , ymin = yhat_lower, ymax= yhat_upper))+
  geom_ribbon(alpha=0.2)+geom_point()+
  geom_point(aes(y=y),color="green")

