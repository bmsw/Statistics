with(mydata, lm(shannon))

with(mydata, plot(shannon, rich))

attach(mydata)
Counts<-rich
Time<-shannon
exponential.model <- lm(log(Counts)~ Time)
summary(exponential.model)
timevalues <- seq(6, 9, 0.01)
Counts.exponential2 <- exp(predict(exponential.model,list(Time=timevalues)))
plot(Time, Counts,pch=16)
lines(timevalues, Counts.exponential2,lwd=2, col = "red", xlab = "Time (s)", ylab = "Counts")
