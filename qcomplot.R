## This script plots the output from the qcom model

qout <- read.table('qout.dat')
plot(seq(286.5,288,length.out=20),seq(1,22,length.out=20)*25,type='n',xlab='Temperature (K)',ylab='Height')

for(i in 1:20){
  if(i==1){
    lcol <- 'blue'
    lwid <- 2
  }else if(i==20){
    lcol <- 'green'
    lwid <- 2
  }else{
    lcol <- 'grey35'
    lwid <- 1
  }
  
lines(qout[i,],(seq(1,22,length.out=22)-1)*25,col=lcol,lwd=lwid)
}

x <- ((qout[1,22]-qout[1,1])/500*seq(0,20,length.out=20)*25) + qout[1,1]
lines(x,seq(1,20,length.out=20)*25,col='red',lwd=2)

legend('topright',c('Initial','Final','Analytical'),col=c('blue','green','red'),lty=1,lwd=2)
title('1-Dimensional QCOM output')