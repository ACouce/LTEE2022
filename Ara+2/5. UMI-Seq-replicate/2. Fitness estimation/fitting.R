#This script estimates fitness by fitting a linear regression model to the frequency trajectories of pooled insertions

# .. initialize graphical output (Linux)
x11(width = 10.4, height = 4)
layout(matrix(c(1,2,3,4),1,4))

# .. load library for nice colors
library("RColorBrewer")
palette(brewer.pal(n = 12, name = "Paired"))
colpal<-palette(brewer.pal(n = 12, name = "Paired"))

# .. load list of neutral alleles and prepare associated objects
neutrals<-read.table("neutral_compilation.dat", header=FALSE, comment.char='') # genes annotated as cryptic, ara operon
neutral_alleles<-ttable$alle[ttable$alle %in% neutrals[,1]]
list_neutral<-c(NA,NA,NA,NA,NA,NA)

# .. set number of generations along the time series
generations<-c(seq(log(100,2), log(100,2)*5,log(100,2))) 

# .. create an empty plot to fill with raw count trajectories
plot(generations, generations, ylim=c(1,1e4), log='y', col='white')

# .. plotting trajectories for neutrals and building a filtered list (t1 raw count > 10)
for (l in 1:length(neutral_alleles)) {
    neutral<-ttable[ttable$alle==as.character(neutral_alleles[l]), 5:9]
    if (neutral$t1<=30) {next}
    lines(generations, neutral[1,], col=colpal[floor(runif(10,1,11+1))], cex=0.5, lwd=0.5)
    list_neutral<-rbind(list_neutral, neutral)
}
list_neutral<-list_neutral[-1,]

# .. estimating and plotting consensus "wild-type"
cons_wt<-matrix(NA, 5,1)
mincov<-quantile(list_neutral$t1, 0.25)
for (i in 1:5) {
    cons_wt[i,1]<-quantile(list_neutral[list_neutral[,1]>mincov,i], 0.75, na.rm='TRUE')
}
lines(generations, cons_wt[,1], col='black', lwd=3, type='b')

# .. creating a version of table with frequencies instead of raw counts (normalized by the consensus "wild-type")
allelname<-unique(ttable$alle)
all<-ttable[ttable$alle%in%as.character(allelname[]), 5:9]
all<-all/cons_wt[col(all)]
ttable2<-ttable
ttable2[ttable2$alle%in%as.character(allelname[]), 5:9]<-ttable2[ttable2$alle%in%as.character(allelname[]), 5:9]/cons_wt[col(all)]

# .. create an empty plot to fill with the frequency trajectories of neutral alleles (should produce mostly straight horizontal lines)
plot(generations, generations, ylim=c(0.003,15), log='y', col='white')

# .. plotting frequency trajectories for neutrals
for (l in 1:length(neutral_alleles)) {
    neutral<-ttable2[ttable2$alle==as.character(neutral_alleles[l]), 5:9] 
    lines(generations, neutral[1,], col=colpal[floor(runif(10,1,11+1))], cex=0.5, lwd=0.5)
}

# .. prepare objects to conduct the fitting
fitted <- matrix(NA,nrow=alleles, ncol=12)
ratio <- as.matrix(all[1,1:5])

# .. create an empty plot to fill with the frequency trajectories of all alleles
plot(generations, log(ratio), ylim=c(-6,4), col='white')

# .. loop to estimate fitness as slopes of frequency trajectories relative to the consensus "wild-type", and record details about the fit to a linear regression model
for (i in 1:alleles) {
	ratio <- as.matrix(all[i,1:5])
	ratio[is.na(ratio)=='TRUE']<-0
	if (sum(ratio[1:5]==0)<3) {
		fit<-lm(log(ratio[ratio>0])~generations[ratio>0])
        	fit1<-lm(log(ratio[ratio>0])~generations[ratio>0], weights=ratio[ratio>0]^1)
        	fit2<-lm(log(ratio[ratio>0])~generations[ratio>0], weights=ratio[ratio>0]^2)

		fitted[i,1]<- coef(fit)["generations[ratio > 0]"]
        	fitted[i,2]<- coef(fit1)["generations[ratio > 0]"]
        	fitted[i,3]<- coef(fit2)["generations[ratio > 0]"]

		fitted[i,4]<- summary(fit)$coef[2,2]
		fitted[i,5]<- summary(fit1)$coef[2,2]
		fitted[i,6]<- summary(fit2)$coef[2,2]

		fitted[i,7]<- summary(fit)$coef[2,4]
		fitted[i,8]<- summary(fit1)$coef[2,4]
		fitted[i,9]<- summary(fit2)$coef[2,4]

		fitted[i,10]<- max(cooks.distance(fit))/min(cooks.distance(fit))	# root mean squared error 
        	fitted[i,11]<- max(cooks.distance(fit1))/min(cooks.distance(fit1))
        	fitted[i,12]<- max(cooks.distance(fit2))/min(cooks.distance(fit2))

       		if ((i %% 100)==0) {						# plot 1% of all data
       			points(generations, log(ratio), type='o', col=colpal[floor(runif(10,1,11+1))], cex=0.3, lwd=0.3)
       		}
	} else {
         	fitted[i,1:4]<-NA       
        }
}

# .. prepare objects to plot the distribution of fitness effects of the neutral set
wts<-ttable[ttable$alle%in%neutrals[,1], ]
wts<-fitted[rownames(ttable)%in%rownames(wts),2]
hist(wts, col='grey', breaks=20, xlim=c(-0.2,0.1), main=filename2)

# .. update table by adding the fit and its details
ttable$fitted<-fitted[,1]
ttable$fitted1<-fitted[,2]
ttable$fitted2<-fitted[,3]
ttable$sterr<-fitted[,4]
ttable$sterr1<-fitted[,5]
ttable$sterr2<-fitted[,6]
ttable$pval<-fitted[,7]
ttable$pval1<-fitted[,8]
ttable$pval2<-fitted[,9]
ttable$rmse<-fitted[,10]
ttable$rmse1<-fitted[,11]
ttable$rmse2<-fitted[,12]

# .. write the table for further reference
write.table(ttable, file = filename3, sep = "\t")

