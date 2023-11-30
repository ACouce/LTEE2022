# .. This script resamples datasets to plot cumulative DFEs with 90% confidence intervals

# .. make condensed version of data sets with only fitness effects and standard errors
Rsimple<-Rtable[Rtable$fitted1>-0.25 & Rtable$fitted1<0.15, c('fitted1', 'sterr1')]
Tsimple<-Ttable[Ttable$fitted1>-0.25 & Ttable$fitted1<0.15,c('fitted1', 'sterr1')]
Fsimple<-Ftable[Ftable$fitted1>-0.25 & Ftable$fitted1<0.15,c('fitted1', 'sterr1')]

# .. retrieve required dimensions and build objects to store the confidence intervals around the cumulative DFEs
qq<-qqplot(Rsimple$fitted1, Tsimple$fitted1, plot=FALSE)
CIS_RT<-matrix(0,200,length(qq$x))
CIS_RF<-matrix(0,200,length(qq$x))

# .. conduct resampling with replacement using the logarithm of the regression standard errors as weights
for (b in 1:200) { 

	Rr<-sample(Rsimple$fitted1, nrow(Rsimple), replace = TRUE, prob= log(1/Rsimple$sterr1))
	Tr<-sample(Tsimple$fitted1, nrow(Tsimple), replace = TRUE, prob= log(1/Tsimple$sterr1))
	Fr<-sample(Fsimple$fitted1, nrow(Fsimple), replace = TRUE, prob= log(1/Fsimple$sterr1))

	RT<-qqplot(Rr, Tr, cex=1.5, lwd=0.5, col=rgb(0.1,0.56,1,0.8), xlim=c(-0.25,0.1), ylim=c(-0.25,0.1), xaxt = 'n', yaxt = 'n', type='l', plot=FALSE)
	RF<-qqplot(Rr, Fr, cex=1.5, lwd=0.5, col=rgb(0.1,0.56,1,0.8), xlim=c(-0.25,0.1), ylim=c(-0.25,0.1), xaxt = 'n', yaxt = 'n', type='l', plot=FALSE)

	CIS_RT[b,1:length(RT$y)]<-RT$y-RT$x
	CIS_RF[b,1:length(RF$y)]<-RF$y-RF$x	
}

# .. initialize graphical output
png("cDFEs+CI.png", width = (15.63*0.9)*(2/2.22), height = (8.54)/2.22, units = 'in', res = 300)
par(mfrow=c(1,4), family="sans", cex=1)

# .. plot cumulative DFE for 2K
RT<-qqplot(Rsimple$fitted1, Tsimple$fitted1, xlim=c(-0.25,0.1), ylim=c(-0.25,0.1), cex=1, lwd=0.5, col='white', xlab='selection coeff. (Anc)', ylab='selection coeff. (2K)', xaxt = 'n', yaxt = 'n', cex.lab=1.5)

# .. compute upper and lower limits for confidence interval
var<-apply(CIS_RT, 2, function(x) quantile(x,0.975)-quantile(x,0.025))
UL<-RT$y+var[1:length(RT$x)]
LL<-RT$y-var[1:length(RT$x)]

# .. smooth these upper and lower limits
x<-RT$x
y<-UL
valUL <- loess(y ~ x, span=0.05)
y<-LL
valLL <- loess(y ~ x, span=0.05)

# .. plot confidence interval region
polygon(x = c(RT$x,rev(RT$x)), y = c(predict(valUL), rev(predict(valLL))), col= rgb(0.8,0.43,0.42,0.25), border=NA)

# .. plot actual points
points(RT$x, RT$y, col=rgb(0.8,0.43,0.42,0.9), lwd=0.75)

# .. custom axes
axis(side = 1, at = seq(-0.25, 0.1, by=0.05), labels=NA,  cex.axis=1.33)
axis(side = 1, at = seq(-0.2, 0.1, by=0.1), labels = seq(-0.2, 0.1, by=0.1),  cex.axis=1.33)
axis(side = 2, at = seq(-0.25, 0.1, by=0.05), labels=NA,  cex.axis=1.33)
axis(side = 2, at = seq(-0.2, 0.1, by=0.1), labels = seq(-0.2, 0.1, by=0.1),  cex.axis=1.33)

# .. add reference lines
abline(a = 0, b = 1, lty = 2, lwd=1)
abline(v=0, lty = 2, lwd=1)
abline(h=0, lty = 2, lwd=1)

# .. homogenize border lines
box(lwd=1.5)

# .. plot cumulative DFE for 15K
RF<-qqplot(Rsimple$fitted1, Fsimple$fitted1, xlim=c(-0.25,0.1), ylim=c(-0.25,0.1), cex=1, lwd=0.5, col='white', xlab='selection coeff. (Anc)', ylab='selection coeff. (2K)', xaxt = 'n', yaxt = 'n', cex.lab=1.5) 

# .. compute upper and lower limits for confidence interval
var<-apply(CIS_RF, 2, function(x) quantile(x,0.975)-quantile(x,0.025))
UL<-RF$y+var[1:length(RF$x)]
LL<-RF$y-var[1:length(RF$x)]

# .. smooth these upper and lower limits
x<-RF$x
y<-UL
valUL <- loess(y ~ x, span=0.05)
y<-LL
valLL <- loess(y ~ x, span=0.05)

# .. plot confidence interval region
polygon(x = c(RF$x,rev(RF$x)), y = c(predict(valUL), rev(predict(valLL))), col= rgb(0.42,0.57,0.8,0.2), border=NA)

# .. plot actual points
points(RF$x, RF$y, col=rgb(0.42,0.57,0.8,0.6), lwd=0.75)

# .. custom axes
axis(side = 1, at = seq(-0.25, 0.1, by=0.05), labels=NA,  cex.axis=1.33)
axis(side = 1, at = seq(-0.2, 0.1, by=0.1), labels = seq(-0.2, 0.1, by=0.1),  cex.axis=1.33)
axis(side = 2, at = seq(-0.25, 0.1, by=0.05), labels=NA,  cex.axis=1.33)
axis(side = 2, at = seq(-0.2, 0.1, by=0.1), labels = seq(-0.2, 0.1, by=0.1),  cex.axis=1.33)

# .. add reference lines
abline(a = 0, b = 1, lty = 2, lwd=1)
abline(v=0, lty = 2, lwd=1)
abline(h=0, lty = 2, lwd=1)

# .. homogenize border lines
box(lwd=1.5)

# .. close graphical output
dev.off()
