### This script plots data output from the parcel model

graphics.off()

dat <- read.table('parcelout.txt')

names(dat) = c('Pressure','TempK','Theta','eTheta','Wl','Wv','Wtot','RH')

plot(dat$TempK,dat$Pressure/100,xlim=c(220,400),ylim=c(1000,200),type='l',
		xlab='Temperature [K]',ylab='Pressure [mb]',lwd=2)

lines(dat$Theta,dat$Pressure/100,col='red',lwd=2)
lines(dat$eTheta,dat$Pressure/100,col='blue',lwd=2)

legend('bottomleft',c('T',expression(theta),
						expression(theta[e])),col=c('black','red','blue'),lty=1,lwd=2)
title('Parcel Model Output')



dev.new()

plot(log10(dat$Wv),dat$Pressure/100,type='l',ylim=c(1000,250),xlim=c(-6,-1),
					xlab=bquote(paste(log[10],' Mixing Ratio [kg/kg]')),ylab='Pressure [mb]',lwd=2)

lines(log10(dat$Wl),dat$Pressure/100,col='red',lwd=2)
lines(log10(dat$Wtot),dat$Pressure/100,col='blue',lwd=2)

legend('topright',c('Water Vapor','Liquid Water','Total Water'),
							col=c('black','red','blue'),lty=1,lwd=2)

title('Parcel Model Output')

#Output from part (c 1): used as initial conditions for part (c 2)
#        Pressure   TempK    Theta     eTheta    Wl            Wv            Wtot
#5000    24985 		227.6361 338.3360  339.3153  0.0005434013  0.0002642887  0.0008076900
