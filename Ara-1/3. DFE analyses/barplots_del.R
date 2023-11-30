# ..  This script plots the most deleterious alleles in the ancestral background and their fitness effects in the 2K and 15K backgrounds

# .. load deleterious tail using thdel as threshold for approximate neutrality
thdel<--0.015
R_del<-Rtable[!is.na(Rtable$fitted1) & Rtable$fitted1>=-1 & Rtable$fitted1<=thdel, ]
R_del<-R_del[!duplicated(R_del$site),] 			# to impede overlapping regions being counted twice!

# .. preparing objects for a loop to create a table with effects of same allele in Anc, 2K and 15K
RORFs<-as.character(R_del$ORF)
Rlist<-as.character(R_del$alle)
Rtabs<-as.data.frame(matrix(0, length(Rlist), 3))

for (i in 1:length(Rlist)) {
	orf<-R_del[which(R_del$alle==Rlist[i]),]$ORF
	
	Rtabs[i,1]<-R_del[which(R_del$alle==Rlist[i]),]$fitted1
	
	if (length(Ttable[which(Ttable$alle==Rlist[i]),]$fitted1)>0) {		# if exact same allele is available, use it..
		Rtabs[i,2]<-Ttable[which(Ttable$alle==Rlist[i]),]$fitted1	
	} else {Rtabs[i,2]<-mean(Ttable[which(Ttable$ORF==orf),]$fitted1)}	# ..otherwise, take locus average
	
	if (length(Ftable[which(Ftable$alle==Rlist[i]),]$fitted1)>0) {		# if exact same allele is available, use it..
		Rtabs[i,3]<-Ftable[which(Ftable$alle==Rlist[i]),]$fitted1
	} else {Rtabs[i,3]<-mean(Ftable[which(Ftable$ORF==orf),]$fitted1)}	# ..otherwise, take locus average
}

# .. add useful columns to the table
Rtabs$name<-RORFs
Rtabs$site<-R_del$site
Rtabs$abn<-R_del$abn
Rtabs$block<-NA

# .. add column names
colnames(Rtabs) <- c('R', 'M', 'K', 'name', 'site','abn','block')

# .. sort by site (useful for clustering)
Rtabs<-Rtabs[order(Rtabs$site),]

# .. load local context information (uncomment in case 'barplots_ben.R' has not been run yet)
#locctx<-read.table(file = 'R606_localctx.txt', sep = "\t", header = FALSE)
#colnames(locctx) <- c('name', 'str', 'end', 'std', 'ctx')

# .. create a matrix to store block labels (clustered by transcriptional units, as per locctx)
g<-matrix(NA,dim(Rtabs)[1],1)
k<-0

# .. a loop to assign alleles to same block if they belong to the same TU
for (i in 1:dim(Rtabs)[1]) {

	exctname<-paste("\\b",Rtabs[i,]$name,"\\b", sep="")	# trick to ensure perfect match (grep is promiscuous)

	if (!is.na(g[i])) {					# if already assigned a block, keep it and assign it to its neighbours (de facto merging contiguous TUs)

		partners<-which(Rtabs$name%in%locctx[grep(exctname,locctx$ctx),]$name)
		g[partners]<-g[i]
		
	} else {
		if (length(Rtabs[Rtabs$name%in%locctx[grep(exctname,locctx$ctx),]$name,]$name)<=1) { 	# if is a singleton, assign it to a new block
			k<-k+1
			g[i]<-k
		} else {										# otherwise, add their neighbours too to a new block
			partners<-which(Rtabs$name%in%locctx[grep(exctname,locctx$ctx),]$name)
			if (any(!is.na(g[partners]))) {							# to deal with cases in which an allele has been already assigned to a contiguous TU
				g[partners]<-min(g[partners], na.rm = TRUE)
			} else {
				k<-k+1
				g[partners]<-k			
			}
		}		
	}
}

# .. add blocks to table
Rtabs$block<-g

#let's ensure we can compare datasets
Rtabs<-na.omit(Rtabs)

# .. lets sort the table by the fitness of each block 
# .. first, a loop to compute fitness of each block (upper quartile)
blck_fit<-matrix(0, nrow(Rtabs), 1)
for (i in 1:nrow(Rtabs)) {
	blck_fit[i]<-quantile(Rtabs[Rtabs$block%in%(Rtabs[i,]$block),]$R, 0.75)
}


# .. second, blck_fit to do the sorting 
Rtabs$val<-blck_fit
Rtabs<-Rtabs[order(Rtabs$val, decreasing=TRUE),]

# .. rename blocks in order to help with visualization
rename<-matrix(NA, nrow(Rtabs), 1)
r<-0
for (i in 1:nrow(Rtabs)) {
	if (!is.na(rename[i])) {		# if already assigned a block, keep it and assign it to its neighbours
		next }	
	r<-r+1
	rename[which(Rtabs$block%in%(Rtabs[i,]$block))]<-r
}
Rtabs$block<-rename

# .. prepare final table to plot by keeping the non singletons (i.e., "multiple")
multiple<-which(as.vector(table(Rtabs$block))>1)
toplot<-Rtabs[Rtabs$block%in%multiple,]

# .. let's only show the average fitness per locus (instead of individual values for each segment)
delorf<-unique(toplot$name)
tsingle<-as.data.frame(matrix(0, length(delorf), 8))
n<-0


for (i in delorf) {
		n<-n+1
		tsingle[n,1:3]<-apply(toplot[toplot$name%in%i,1:3], 2, mean,na.rm=TRUE)	# average fitness per locus 
		tsingle[n,4:7]<-toplot[toplot$name%in%i,4:7][1,]				# just copy details
		tsingle[n,8]<-mean(as.numeric(tsingle[n,1:3]), na.rm=TRUE)		# estimate average accross backgrounds
	}
toplot<-tsingle[,1:8]

# .. add column names
colnames(toplot) <- c('R', 'M', 'K', 'name', 'site','abn','block','bvalue')

# .. lets sort the again, now by the average fitness accross backgrounds
toplot<-toplot[order(toplot$bvalue),]

# .. to ensure we visualize the most deleterious accross backgrounds
toplot<-toplot[toplot$R<=-0.033 & toplot$M<=-0.033 & toplot$K<=-0.033,]

# .. to update block labels (we missed some along the way)
oldblocks<-toplot$block
newblocks<-1:nrow(toplot)
n=0
for (b in unique(toplot$block)) {
	n=n+1
	print(c(b,n))
	newblocks[which(oldblocks==b)]<-n
}
toplot$block<-newblocks
toplot<-toplot[order(toplot$block),]

# .. initialize graphical output
png(file="dels.png", width=1100, height=500*0.75, res=48)
par(family="sans", cex=2)

# .. show the most beneficial 55 loci in the ancestor (38 blocks, cut just for visualization purposes)
toplot<-toplot[toplot$block<toplot[55,]$block,]

# .. let's make a barplot and store coordinates to add white and gray shaded areas grouping loci from same block
coord<-barplot(t(toplot[,1:3]), beside = TRUE, space=c(0,1), axes = FALSE, axisnames=FALSE, ylim=c(-0.25, 0), col='white', las=2, cex.names=1.5, cex.axis=1.5, horiz = FALSE, add = FALSE, border = 'NA')
coord<-rbind(coord,t(toplot$block))

# .. adding white and gray shaded areas
for (r in seq(1, max(toplot$block), by=2)) {
	rect(min(coord[1:3,coord[4,]==r+1])-1,-3,max(coord[1:3,coord[4,]==r+1])+1,3,col="white",lty=1, lwd=0.25)
	rect(min(coord[1:3,coord[4,]==r])-1,-3,max(coord[1:3,coord[4,]==r])+1,3,col="grey85",lty=1, lwd=0.25)
	}

# .. the actual barplot
barplot(t(toplot[,1:3]), beside = TRUE, space=c(0,1), yaxt='n', names.arg=toplot$name, font=3, ylim=c(-0.25, 0), col=c('black', col2K, col15K), las=2, cex.names=1.1, cex.axis=1.2, horiz = FALSE, add = TRUE, lwd=1.5)

# ..  create matrix to store average standard errors
ster<-matrix(NA,nrow(toplot),3)
for (s in 1:nrow(toplot)) {
	name<-toplot[s,]$name
	ster[s,1]<-mean(Rtable[Rtable$ORF==name,]$sterr1, na.rm=TRUE)
	ster[s,2]<-mean(Ttable[Ttable$ORF==name,]$sterr1, na.rm=TRUE)
	ster[s,3]<-mean(Ftable[Ftable$ORF==name,]$sterr1, na.rm=TRUE)
}

# .. use standard errors to add errors bars indicate ~95% confidence intervals
segments(coord[1,], toplot$R-2*ster[,1], coord[1,], toplot$R+2*ster[,1])
segments(coord[2,], toplot$M-2*ster[,2], coord[2,], toplot$M+2*ster[,2], col=col2K)
segments(coord[3,], toplot$K-2*ster[,3], coord[3,], toplot$K+2*ster[,3], col=col15K)

# .. custom y-axis
axis(2, at = seq(-0.25, 0.1, 0.025), labels=NA, lwd=2, cex.axis = 1, tck=-0.05)
axis(2, at = seq(-0.25, 0.1, 0.05),  labels=seq(-0.25, 0.1, 0.05), lwd=2, cex.axis = 1.2, tck=-0.05)
abline(h=0, lty=2)

# .. close graphical output
dev.off()

# .. initialize graphical output
png(file="dels2.png", width=300*0.8, height=460*0.85, res=48)
par(family="sans", cex=2)

# .. plot summary barplots
boxplot(toplot[toplot$R<=-0.033 & toplot$M<=-0.033 & toplot$K<=-0.033,1:3], ylim=c(-0.25, 0), col=c('grey24', col2K, col15K), cex.names=1, cex=0.8, names=c('Anc', '2K', '15K'), ylab='selection coeff. (s)', cex.lab=1.5,  cex.axis=1.2, lwd=1.2,  yaxt='n', las=2)

# .. custom y-axis
axis(2, at = seq(-0.25, 0.1, 0.025), labels=NA, lwd=2, cex.axis = 1, tck=-0.05)
axis(2, at = seq(-0.25, 0.1, 0.05),  labels=seq(-0.25, 0.1, 0.05), lwd=2, cex.axis = 1.2, tck=-0.05)
abline(h=0, lty=2)

# .. homogenize border lines
box(lwd=1.5)

# .. close graphical output
dev.off()
