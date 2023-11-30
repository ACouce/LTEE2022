# .. This script estimates and plots pairwise comparisons among segments from the same locus 

# .. load filtered dataset
ttable<-read.table(file = "Rfitted_fil.txt", sep="\t", header=TRUE, as.is=TRUE, comment.char = "")
x<-ttable[(ttable$fitted1>0 & ttable$t1>=10 & ttable$abn>1 & ttable$sterr1<0.01) | (ttable$fitted1<=0 & ttable$t1>=10 & ttable$abn>1 & ttable$sterr1<1) ,c(1,2,4,6,10,12,14,17)]

# .. create empty list to store pairwise comparisons
list<-c(0,0,0,0,0,0,0)

# .. retrieve loci for which more than one segment has observations (otherwise, pairwise comparisons ara not possible)
tab<-table(x$ORF)
loci<-names(tab[tab>1])

# .. loop for the pairwise comparisons among internal segments (N- or C-terminal segments are considered only if show no clear deviation)
n<-0
for (r in loci) {
	alls<-x[x$ORF%in%r,]											# load alleles
	if (length(grep(r,alls$alle)) != length(alls$alle)) {next}						# to deal with overlapping loci
	pwise<-t(combn(alls$fitted1,2))										# create matrix with all possible combinations
	var<-abs(pwise[,1]-pwise[,2])										# estimate absolute difference across all possible combinations 	
	dev<-as.numeric(which.max(colSums(as.matrix(dist(alls$fitted1, upper=TRUE)))))				# identify combination with highest difference (index)
	max<-alls[dev,]$alle											# identify combination with highest difference (name)

	if (abs(alls[dev,]$fitted1-mean(alls$fitted1))/abs(mean(alls$fitted1))>=0.25) {				# if deviation represents > 25% of locus average (as in filter_and_count.R)	
		if (length(grep('IR',r))==1) {next}								# skip loci for which functional similarity accross internal fragments is not guaranteed (i.e., intergenic regions, pseudogenes, etc.)
		if (length(grep('ESCRE',r))==1) {next}
		if (length(grep('ECB',r))==1) {next}
		if (length(grep('#',r))==1) {next}

		if (nrow(alls)==2) {										# if only two segments have observations
			alls<-alls[grep('-2|-3|-4',alls$alle),]							# retrieve internal segments
			if (nrow(alls)<2) {next}									# if one of the two is a terminal segment, skip
			pwise<-t(combn(alls$fitted1,2))								# otherwise, record pairwise difference
			var<-abs(pwise[,1]-pwise[,2])
			
		} else if (nrow(alls)>2) {									# if more than two segments have observations	
			if (length(grep('-1|-5', max)) == 1) {							# is maximum deviation due to a N- or C-terminal segment?
				alls<-alls[alls$alle!=max,]							# if true, remove it and recalculate pairwise differences
				pwise<-t(combn(alls$fitted1,2))
				var<-abs(pwise[,1]-pwise[,2])							# provisionally keeping the other N- or C-terminal segment
				dev<-as.numeric(which.max(colSums(as.matrix(dist(alls$fitted1, upper=TRUE))))) 
				max<-alls[dev,]$alle
				
				if (abs(alls[dev,]$fitted1-mean(alls$fitted1))/abs(mean(alls$fitted1))>=0.25) {	# is maximum deviation still > 25% of locus average?
					dev<-abs(alls$fitted1-mean(alls$fitted1))
					max<-alls[which(dev==max(dev)),]$alle
					if (sum(grepl('-1|-5',max))>0) {						# is maximum deviation due to the other N- or C-terminal segment?
						alls<-alls[alls$alle!=max,]					# if true, remove it and recalculate pairwise differences
						if (nrow(alls)<2) {next}						# if no other segments are left, skip
						pwise<-t(combn(alls$fitted1,2))					# otherwise, record pairwise difference
						var<-abs(pwise[,1]-pwise[,2])
					}
				}
			} else {			
				n<-n+1
				alls<-alls[grep('-2|-3|-4',alls$alle),]						# retrieve internal segments
				if (nrow(alls)<2) {next}								# if no other segments are left, skip
				pwise<-t(combn(alls$fitted1,2))							# record pairwise difference
				var<-abs(pwise[,1]-pwise[,2])
			}
		}
	}

	pwise<-cbind(r, t(combn(alls$fitted1,2)),var,mean(var),t(combn(alls$alle,2)))				# create element of general list
	list<-rbind(list, pwise)											# update list
}

# .. remove first element (empty)
list<-list[-1,]

# .. create a data frame for greater ease of manipulation
df<-data.frame(ORF=list[,1], V2=as.numeric(list[,2]), V3=as.numeric(list[,3]), var=as.numeric(list[,4]))

# .. initialize graphical output
png(file="bulk.png", width=1500, height=1800, res=135)
par(mfrow=c(2,2), family="sans", cex=2)
par(lwd=1.5)

# .. calculate global correlation
cr<-round(cor(df[,]$V2, df[,]$V3, use='complete.obs'),digits=3)

# .. add a column with average fitness for each pairwise comparison
df$ave<-(df$V2+df$V3)/2

# .. create an empty plot with desired details and parameters
plot(df[df$V3<0.015 | df$V2<0.015,]$V2, df[df$V3<0.015 | df$V2<0.015,]$V3, col='white', ylim=c(-0.24,0.1), xlim=c(-0.24,0.1), main=cr, xlab='Error (s)', ylab='Fitness (s)')

# .. add useful reference lines
abline(v=0,lty=2, col='black')
abline(h=0,lty=2, col='black')
abline(0,1, lty=2)

# .. add the actual points (deleterious)
points(df[df$V3<0.015 | df$V2<0.015,]$V2, df[df$V3<0.015 | df$V2<0.015,]$V3, lwd=abs(df[df$V3<0.015 | df$V2<0.015,]$ave*10)+0.1, pch=1, cex=0.4, col=grey(0,0.25), ylim=c(-0.24,0.1), xlim=c(-0.24,0.1))

# .. add the actual points (beneficial)
points(df[df$V3>0.015 | df$V2>0.015,]$V2, df[df$V3>0.015 | df$V2>0.015,]$V3, lwd=abs(df[df$V3>0.015 | df$V2>0.015,]$ave*10)+0.1, pch=1, cex=0.4, col=rgb(1,0.4,0.3,0.75))

# .. set a lower limit for visualization purposes
df[df$var<0.0001,]$var<-0.0001

# .. create custom palette
dr<-as.vector(c(col2rgb('#ca0020'),180))/255
lr<-as.vector(c(col2rgb('#f4a582'),180))/255
n<-as.vector(c(col2rgb('#f7f7f7'),180))/255
lb<-as.vector(c(col2rgb('#92c5de'),180))/255
db<-as.vector(c(col2rgb('#0571b0'),180))/255

# .. prepare objects for boxplot
xHL<-rep(1,nrow(df[df$ave<=-0.05,]))+rnorm(nrow(df[df$ave<=-0.05,]), 0,0.1)
xML<-rep(2,nrow(df[df$ave>-0.05 & df$ave<=-0.015,]))+rnorm(nrow(df[df$ave>-0.05 & df$ave<=-0.015,]), 0,0.1)
xN<-rep(3,nrow(df[df$ave>-0.015 & df$ave<=0.015,]))+rnorm(nrow(df[df$ave>-0.015 & df$ave<=0.015,]), 0,0.1)
xMB<-rep(4,  nrow(df[df$ave>0.015 & df$ave<=0.05,]))+rnorm(nrow(df[df$ave>0.015 & df$ave<=0.05,]), 0,0.1)
xHB<-rep(5, nrow(df[df$ave>0.05,]))+rnorm(nrow(df[df$ave>0.05,]), 0,0.1)

# .. plot individual points before boxplot
plot(c(xHL,xML,xMB,xHB), c(df[df$ave<=-0.05,]$var, df[df$ave>-0.05 & df$ave<=-0.015,]$var, df[df$ave>0.015 & df$ave<=0.05,]$var, df[df$ave>0.05,]$var), xlim=c(0.25,5.75), ylim=c(0.0001,0.1), log='y', xaxt='n', yaxt='n', xlab=NA, ylab=NA, lwd=0.75, cex=0.4, col=grey(0,0.3))
points(c(xN), c(df[df$ave>-0.015 & df$ave<=0.015,]$var), xlim=c(0,6), ylim=c(0.0005,0.1), lwd=0.3, cex=0.45, col=grey(0,0.25))

# .. add the boxplot
boxplot(df[df$ave<=-0.05,]$var, df[df$ave>-0.05 & df$ave<=-0.015,]$var, df[df$ave>-0.015 & df$ave<=0.015,]$var, df[df$ave>0.015 & df$ave<=0.05,]$var, df[df$ave>0.05,]$var, ylim=c(0.0001,0.1), col=c(rgb(db[1], db[2], db[3], db[4]), rgb(lb[1], lb[2], lb[3], lb[4]), rgb(n[1], n[2], n[3], n[4]), rgb(lr[1], lr[2], lr[3], lr[4]), rgb(dr[1], dr[2], dr[3], dr[4])), log='y', add=TRUE, outline=FALSE,  xaxt='n', yaxt='n', xlab=NA, ylab=NA)

# .. custom y-axis
ylab=seq(-4,-1,by=1)
yat=10^ylab
axis(2, at=yat, labels=c(0.0001, 0.001, 0.01, 0.1), las=3, lwd=1, tcl=-0.75, cex.axis=1, )

ylab=seq(2,9, by=1)
for (l in seq(-4,-1,by=1)) {
	yat=ylab*10^l
	axis(2,at=yat, labels=NA, las=1, lwd=0.5, tcl=-0.5)
}

# .. custom x-axis
axis(1, at=1:5, labels=c('delet.', 'mild del.', 'neutral', 'mild ben.', 'benef.'), las=2)

# .. graphical parameters for last plots (histograms)
line<-par(lwd=1)
breaks<-seq(-0.075,0.1, length.out=50)

# .. plot first histogram
hist(c(as.numeric(list[as.numeric(list[,3])>0.015,3]),(as.numeric(list[as.numeric(list[,2])>0.015,2]))), xlim=c(-0.1,0.1), ylim=c(0,200), breaks=breaks, col=rgb(1,0.4,0.3,0.75), xlab='Fitness (s)', ylab='Number of loci', xaxt='n', yaxt='n', main=NA)

# .. custom axis line width
axis(side = 1, lwd = 1.5)
axis(side = 2, lwd = 1.5)

# .. useful reference line
abline(v=0,lty=2, col='black')

# .. plot second histogram
hist(c(as.numeric(list[as.numeric(list[,3])>0.015,2]),(as.numeric(list[as.numeric(list[,2])>0.015,3]))), xlim=c(-0.1,0.1), ylim=c(0,200), breaks=breaks, col=rgb(1,0.4,0.3,0.75), xlab='Fitness (s)', ylab='Number of loci', xaxt='n', yaxt='n', main=NA)

# .. custom axis line width
axis(side = 1, lwd = 1.5)
axis(side = 2, lwd = 1.5)

# .. useful reference line
abline(v=0,lty=2, col='black')

# .. close graphical output
dev.off()
