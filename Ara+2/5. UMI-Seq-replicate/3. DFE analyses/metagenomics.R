# ..  This script analyzes the determinants of the observed frequency of beneficial mutations in the LTEE. It plots the relationship bewtween these frequencies with both our fitness estimates and the mutational target size.

# .. load list with start and end positions, and sense of all loci
ORFs<-read.table('parsed_R606genoscope_IU.txt', sep="\t", header=FALSE, as.is=TRUE, comment.char = "")
colnames(ORFs)<-c('ORF','str','end','S')

# .. load metagenomic time-series data for the LTEE (only alleles appearing in at least two independent populations) 
Mikemap<-read.table(file = 'Mikemap_par.dat', sep = "\t", header = TRUE)
Mikemap[Mikemap$reanno%in%c('rbsR','hsrA') & Mikemap$ann=='sv',]$reanno<-'rbsB' #... manually correct table to standardize annotation of mutations in the rbs operon and facilitate identification of parallelism (see Fig. 3 from Cooper (2001) J Bacteriol.)

# .. load table with most beneficial alleles in the ancestral background (produced by 'barplots_ben.R')
Rtabs<-read.table(file = 'Rtoanal.dat', sep = "\t", header = TRUE)

# .. create table to record statistics to understand the prevalence of alleles in the metagenomic data
Rblocks<-as.data.frame(sort(unique(Rtabs$block)))
colnames(Rblocks) <- c('block')
Rblocks$name<-NA
Rblocks$R<-NA
Rblocks$M<-NA
Rblocks$K<-NA
Rblocks$med<-NA
Rblocks$max<-NA
Rblocks$per<-NA
Rblocks$nobs1<-NA
Rblocks$nobs5<-NA
Rblocks$nobs50<-NA
Rblocks$nobs90<-NA
Rblocks$nobs100<-NA
Rblocks$nlin<-NA
Rblocks$tmin<-NA
Rblocks$size<-NA
Rblocks$lof<-NA

for (i in 1:dim(Rblocks)[1]) {
	b<-Rblocks[i,]$block														# retrieve block
	names<-as.character(unique(Rtabs[Rtabs$block==b,]$name))										# retrieve loci names in block
	loc<-unique(unlist(strsplit(as.character(locctx[locctx$name%in%names,]$ctx), split=",")))						# retrieve loci names of block's local neighborhood (to detect polar effects)
	Rblocks[i,]$R<-max(Rtabs[Rtabs$block==b,]$R, na.rm = TRUE)									# largest fitness in block (Anc)
	Rblocks[i,]$M<-max(Rtabs[Rtabs$block==b,]$M, na.rm = TRUE)									# largest fitness in block (2K)
	Rblocks[i,]$K<-max(Rtabs[Rtabs$block==b,]$K, na.rm = TRUE)									# largest fitness in block (15K)	
	Rblocks[i,]$med<-median(c(Rtabs[Rtabs$block==b,]$R,Rtabs[Rtabs$block==b,]$M,Rtabs[Rtabs$block==b,]$K), na.rm = TRUE)			# median fitness accross backgrounds	
	Rblocks[i,]$max<-max(Mikemap[Mikemap$reanno%in%names & Mikemap$par>2,]$max)							# maximum frequency in the LTEE by any allele in block 
	Rblocks[i,]$per<-sum(Mikemap[Mikemap$reanno%in%names,]$per)									# persistence alleles in block at detectable frequency in the LTEE
	Rblocks[i,]$name<-as.character(Rtabs[which(Rtabs$R==max(Rtabs[Rtabs$block==b,]$R)),]$name)						# name of locus with largest fitness in block
	Rblocks[i,]$nobs1<-sum(Mikemap[Mikemap$reanno%in%names,]$max>=0)									# total number of alleles in block ever reaching a frequency >0% in the LTEE
	Rblocks[i,]$nobs5<-sum(Mikemap[Mikemap$reanno%in%names,]$max>=0.05)								# total number of alleles in block ever reaching a frequency >5% in the LTEE					
	Rblocks[i,]$nobs50<-sum(Mikemap[Mikemap$reanno%in%names,]$max>=0.5)								# total number of alleles in block ever reaching a frequency >50% in the LTEE
	Rblocks[i,]$nobs90<-sum(Mikemap[Mikemap$reanno%in%names,]$max>=0.9)								# total number of alleles in block ever reaching a frequency >90% in the LTEE
	Rblocks[i,]$nobs100<-sum(Mikemap[Mikemap$reanno%in%names,]$max>=0.99)								# total number of alleles in block ever reaching a frequency >99% in the LTEE
	Rblocks[i,]$nlin<-length(unique(Mikemap[Mikemap$reanno%in%names,]$ara))								# total number of parallel populations in which they were detected
	Rblocks[i,]$tmin<-min(Mikemap[Mikemap$reanno%in%names,]$tmax)									# minumum time to reach maximum frequency 
	Rblocks[i,]$size<-sum(abs(ORFs[ORFs$ORF%in%names,]$str-ORFs[ORFs$ORF%in%names,]$end))						# target size of block as a whole
	Rblocks[i,]$lof<-sum(Mikemap[Mikemap$reanno%in%names,]$ann%in%c('frameshift','sv','nonsense'))/dim(Mikemap[Mikemap$reanno%in%names,])[1]	# frequency of LOF mutations among mutations in block
}

# .. sort table by fitness in ancestor
Rblocks<-Rblocks[order(Rblocks$R, decreasing = FALSE),]

# .. initialize graphical output
png(file="corr.png", width = (15.63*0.9)*(2/2.22), height = (8.54)/2.22, units = 'in', res = 300)
par(mfrow=c(1,4), family="sans", cex=1)

# .. create custom palette
values<-Rblocks$R/diff(range(Rblocks$R))
values<-round(values/min(values))^2
cols = colorRampPalette(c("lightblue", "#ff3f0f"))(max(values))

# .. plot correlation between observed beneficial mutations in the LTEE and mutational target size 
plot(Rblocks$size,Rblocks$nobs50, main=round(cor(Rblocks$size,Rblocks$nobs50),2), col=cols[values], pch = 19, lwd=2, cex=Rblocks$R/0.04, xlab='target size (bp)', ylab='observed alleles', cex.lab=1.5, cex.axis=1.33)

# .. homogenize border lines
box(lwd=1.5)

# .. create custom palette
values<-Rblocks$size/diff(range(Rblocks$size))
values<-round((values/min(values))^1)
cols = colorRampPalette(c("#7b3294", "#b8e186"))(max(values))

# .. plot correlation between observed beneficial mutations in the LTEE and fitness estimates
plot(Rblocks$R,Rblocks$nobs50, main=round(cor(Rblocks$R,Rblocks$nobs50),2), col=cols[values], pch = 19, lwd=2, cex=sqrt(Rblocks$size/1500), xlab='selection coeff. (s)', ylab='observed alleles', cex.lab=1.5, cex.axis=1.33)

# .. homogenize border lines
box(lwd=1.5)

# .. close graphical output
dev.off()
