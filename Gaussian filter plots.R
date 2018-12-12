





x<-seq(-4000,4000,length.out=10000)
y1<-dunif(x,-250,250)
# plot 1
plot(x, y1, type="n", lty=1, xlab="distance",ylab="Density", main="")
abline(v=-125,lty=2)
abline(v=125,lty=2)
lines(x,dunif(x,-250,250),lwd=2)
lines(x,dunif(x,-1000,1000),lwd=2, col=2)
legend("topright", inset=.05,  legend=c("250m cell","500m buffer", "1000m buffer"), lwd=c(1,2,2), lty=c(2, 1,1), col=c(1,1,2))

# plot 2
plot(x, y1, type="n", lty=1, xlab="distance",ylab="Density", main="")
abline(v=-125,lty=2)
abline(v=125,lty=2)
lines(x,dunif(x,-250,250),lwd=2)
mix<-0.5*dunif(x,250,1000)+0.5*dunif(x,-1000,-250)
lines(x,mix,col=4,lwd=2)
legend("topright", inset=.05,  legend=c("250m cell","500m buffer", "500-1000m 'doughnut'"), lwd=c(1,2,2), lty=c(2, 1,1), col=c(1,1,4))

# plot 3
plot(x, y1, type="n", lty=1, xlab="distance",ylab="Density", main="")
abline(v=-125,lty=2)
abline(v=125,lty=2)

lines(x, dnorm(x,0,sd=250),lwd=2, col=1)
lines(x, dnorm(x,0,500), lwd=2, col=4)
lines(x, dnorm(x,0,750), lwd=2, col="darkgreen")
legend("topright", inset=.05,  legend=c("250m cell","Gaussian filter (mean=0, sigma=250m)","GF (mean=0, sigma=500m)","GF (mean=0, sigma=750m)"), lwd=c(1,2,2,2), lty=c(2,1,1,1), col=c(1,1,4,"darkgreen"))

# Plot 4
plot(x, y1, type="n", lty=1, xlab="distance",ylab="Density", main="")
abline(v=-125,lty=2)
abline(v=125,lty=2)
mix1<- 0.33*dnorm(x,-1000,500) +0.33*dnorm(x,0,300)+ 0.33*dnorm(x,1000,500)
lines(x, mix1, col=2,lwd=2)
mix2<- 0.5*dnorm(x,-1000,400) +0.5*dnorm(x,1000,400) 
lines(x,mix2, col=4, lwd=2)

legend("topright", inset=.05,  legend=c("250m cell","Gaussian mixture filter (3 modes at -1000, 0, 1000)","Gaussian mixture filter (2 modes at -1000, 1000)"), lwd=c(1,2,2), lty=c(2,1,1), col=c(1,2,4))




# plot 1
plot(x, y1, type="n", lty=1, xlab="distance",ylab="Density", main="")
abline(v=-125,lty=2)
abline(v=125,lty=2)
lines(x,dunif(x,-250,250),lwd=2)
legend("topright", inset=.05,  legend=c("250m cell","500m buffer"), lwd=c(1,2), lty=c(2, 1), col=c(1,1))

# plot 2
plot(x, y1, type="n", lty=1, xlab="distance",ylab="Density", main="")
abline(v=-125,lty=2)
abline(v=125,lty=2)
lines(x,dunif(x,-250,250),lwd=2, lty=2)
lines(x,dunif(x,-1000,1000),lwd=2, col=2)
legend("topright", inset=.05,  legend=c("250m cell","500m buffer", "1000m buffer"), lwd=c(1,2,2), lty=c(2, 2,1), col=c(1,1,2))

# plot 3
plot(x, y1, type="n", lty=1, xlab="distance",ylab="Density", main="")
abline(v=-125,lty=2)
abline(v=125,lty=2)
lines(x,dunif(x,-250,250),lwd=2, lty=2)
lines(x,dunif(x,-1000,1000),lwd=2, col=2, lty=2)
mix<-0.5*dunif(x,300,1500)+0.5*dunif(x,-1500,-300)
lines(x,mix,col=4,lwd=2)
legend("topright", inset=.05,  legend=c("250m cell","500m buffer","1000m buffer", "300-1500m 'doughnut'"), lwd=c(1,2,2,2), lty=c(2,2,2,1), col=c(1,1,2,4))

# plot 4
plot(x, y1, type="n", lty=1, xlab="distance",ylab="Density", main="")
abline(v=-125,lty=2)
abline(v=125,lty=2)
lines(x, dnorm(x,0,sd=250),lwd=2, col=1)
legend("topright", inset=.05,  legend=c("250m cell","Gaussian filter (mean=0, sigma=250m)"), lwd=c(1,2), lty=c(2,1), col=c(1,1))

# plot 5
plot(x, y1, type="n", lty=1, xlab="distance",ylab="Density", main="")
abline(v=-125,lty=2)
abline(v=125,lty=2)

lines(x, dnorm(x,0,sd=250),lwd=2, col=1)
lines(x, dnorm(x,0,500), lwd=2, col=4)
lines(x, dnorm(x,0,750), lwd=2, col="darkgreen")
legend("topright", inset=.05,  legend=c("250m cell","GF (mean=0, sigma=250m)","GF (mean=0, sigma=500m)","GF (mean=0, sigma=750m)"), lwd=c(1,2,2,2), lty=c(2,1,1,1), col=c(1,1,4,"darkgreen"))

# Plot 6
plot(x, y1, type="n", lty=1, xlab="distance",ylab="Density", main="")
abline(v=-125,lty=2)
abline(v=125,lty=2)
mix2<- 0.5*dnorm(x,-1000,400) +0.5*dnorm(x,1000,400) 
lines(x,mix2, col=4, lwd=2)
legend("topright", inset=.05,  legend=c("250m cell","Gaussian mixture 'doughtnut' (2 modes at -1000, 1000)"), lwd=c(1,2), lty=c(2,1), col=c(1,4))

# Plot 7
plot(x, y1, type="n", lty=1, xlab="distance",ylab="Density", main="")
abline(v=-125,lty=2)
abline(v=125,lty=2)
mix2<- 0.5*dnorm(x,-1000,400) +0.5*dnorm(x,1000,400) 
lines(x,mix2, col=4, lwd=2)
mix1<- 0.33*dnorm(x,-1000,400) +0.33*dnorm(x,0,350)+ 0.33*dnorm(x,1000,400)
lines(x, mix1, col=2,lwd=2)
legend("topright", inset=.05,  legend=c("250m cell","Gaussian mixture 'doughnut' (2 modes at -1000, 1000)","Gaussian mixture 'cinnamon roll' (3 modes at -1000, 0, 1000)"), lwd=c(1,2,2), lty=c(2,1,1), col=c(1,4,2))



