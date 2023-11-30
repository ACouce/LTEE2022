# This script compares fitness estimates obtained for the same libraries with the HiSeq versus the UMI-Seq approaches

# .. load data sets
UMI<-read.table(file = "Rfitted_logHS.txt", sep="\t", header=TRUE, as.is=TRUE, comment.char = "")	# ancestor for Ara+2 from UMI-Seq
HS<-read.table(file = "RHSfitted_log_n3.txt", sep="\t", header=TRUE, as.is=TRUE, comment.char = "")	# ancestor for Ara+2 from HiSeq, but with loci divided in just 3 segments (for comparative purposes)

# .. filter data sets to ensure comparison among reliable loci
UMI<-UMI[which(!is.na(UMI$fitted1)),]
HS<-HS[which(!is.na(HS$fitted1)),]
UMI<-UMI[(UMI$fitted1<=-0 & UMI$sterr1<1 & UMI$t1>=5) | (UMI$fitted1>0 & UMI$sterr1<0.012 & UMI$t1>=5),]
HS<-HS[(HS$fitted1<=-0 & HS$sterr1<1 & HS$t1>=20) | (HS$fitted1>0 & HS$sterr1<0.008 & HS$t1>=20) & HS$abn>5,]

# .. load genome coordinates with loci divided in just 3 segments
alls<-read.table('parsed_R606genoscope_IUD_n3_ORF.txt')


# .. retrieve common loci 
comm<-intersect(unique(UMI$ORF), unique(HS$ORF))
comm<-comm[which(!is.na(comm))]

# .. loop to make the comparisons
comp<-matrix(NA,1,2)
for (i in 1:length(comm)) {

	alleles<-as.character(alls[grep(paste0("^",comm[i],'-'), alls$V1),1])		# retrieve alleles

	for (a in alleles) {
	
		H<-UMI[UMI$alle%in%a,]
		if (nrow(H)==0) {next}							# if exact allele is not present, skip
		H<-H[,]$fitted1
		
		M<-HS[HS$alle%in%a,]
		if (nrow(HS[HS$alle%in%a,])==0) {next}					# if exact allele is not present, skip
		M<-M[,]$fitted1
				
		if ((abs(H-M)/abs(median(c(H,M)))>0.75) & (max(abs(H-M))>0.01)) {		# if deviation is greater than 1%, and represetns more than 75% of median fitness values for locus across datasets

			if (((UMI[UMI$alle%in%a,]$t1/HS[HS$alle%in%a,]$t1)<0.4 | (UMI[UMI$alle%in%a,]$t1/HS[HS$alle%in%a,]$t1)>2.5)) {	 # if deviation due to large discrepancies in coverage, skip		
				next
			}				
			if (sum(grepl('-1|-3',a))>0) {					# if deviation involves C- or N-terminal fragments, skip
				next	
			}
		}		
		test<-cbind(H,M)
		comp<-rbind(comp,test)
	}
}

# .. initialize graphical output
png(file="corr_UMISeq.png", width=614*1.04, height=378*1, res=48)
par(mfrow=c(1,2), family="sans", cex=2)

# .. plot version 1
plot(comp[,1], comp[,2], xlim=c( -0.25,0.11), ylim=c( -0.25,0.11), cex=0.5, lwd=1.5, xlab='Fitness (s)', ylab='Fitness (s)', cex.lab=1.3, main=round(cor(comp[,1], comp[,2], use="complete.obs"),3), col='dodgerblue')
abline(0,1, lty=2, lwd=0.5)
abline(v=0, lty=1)
abline(h=0, lty=1)
box(lwd=2)

# .. plot version 2
plot(comp[,1], comp[,2], xlim=c(-0.25,0.11), ylim=c(-0.25,0.11),cex=0.8, lwd=0.5, xlab='selection coeff. (s)', ylab='selection coeff. (s)', cex.lab=1.3, main=cor(comp[,1], comp[,2], use="complete.obs"))
abline(0,1, lty=2)
abline(v=0, lty=1)
abline(h=0, lty=1)
box(lwd=2)

# .. close graphical output
dev.off()
