# ..  This script analyzes how well the estimated DFEs capture the beneficial mutations becoming numerically dominant in the LTEE.

# .. a useful function
'%!in%' <- function(x,y)!('%in%'(x,y))

# .. load local context information
locctx<-read.table(file = 'R606_localctx.txt', sep = "\t", header = FALSE)
colnames(locctx) <- c('name', 'str', 'end', 'std', 'ctx')

# .. load metagenomic time-series data for the LTEE (only alleles appearing in at least two independent populations) 
Mikemap<-read.table(file = 'Mikemap_par.dat', sep = "\t", header = TRUE)
Mikemap[Mikemap$reanno%in%c('rbsR','hsrA') & Mikemap$ann=='sv',]$reanno<-'rbsB' #... manually correct table to standardize annotation of mutations in the rbs operon and facilitate identification of parallelism (see Fig. 3 from Cooper (2001) J Bacteriol.)

# .. define global parameters
replicates<-1									# by definition, just one for the real data (in contrast to simulated data )
gens<-c(seq(1000,50000,500))							# generations at which the meta-genomic data was sampled
label<-'real'									# flag to enable using 'repre.R' both with the real and the simulated data

# .. load the table with the beneficial mutations from the DFE of the ancestor
Ntabs<-read.table(file = 'Rtoanal.dat', sep = "\t", header = TRUE)

# .. define parameters for ancestor
parthres<-2									# paralellism threshold (number of independent populations in which a mutation is observed)
genth<-1000									# first time point to consider
lins<-c('m1', 'm2', 'm3', 'm4', 'm5', 'm6', 'p1', 'p2', 'p3', 'p4', 'p5', 'p6')	# population names (all for the ancestor)
fixed<-c('araA','recD')								# already fixed in the ancestor

# .. run the script that estimates representativity with the above details
source('repre.R')
write.table(rep, 'Rrep.dat', row.names = FALSE, sep = "\t",  col.names = TRUE)	# write table with representativity results (fraction of drivers captured by the DFE)
write.table(LOF, 'LOF.dat', row.names = FALSE, sep = "\t",  col.names = TRUE)	# uncomment here and within 'repre.R' for analyses fo prevalence of loss-of-function mutations

# .. load the table with beneficial mutations from the DFE at 2K
Ntabs<-read.table(file = 'Ttoanal.dat', sep = "\t", header = TRUE)

# .. define parameters for 2K
parthres<-1
genth<-2000
lins<-c('p2')
fixed<-as.character(Mikemap[!is.na(Mikemap$max) & Mikemap$max>=0.5 & Mikemap$ann%!in%c('synonymous') & Mikemap$tmax<=2000 & Mikemap$ara%in%lins,]$reanno)

# .. run the script that estimates representativity with the above details
source('repre.R')
write.table(rep, 'Trep.dat', row.names = FALSE, sep = "\t",  col.names = TRUE)	# write table with representativity results (fraction of drivers captured by the DFE)

# .. load the table with beneficial mutations from the DFE at 15K
Ntabs<-read.table(file = 'Ftoanal.dat', sep = "\t", header = TRUE)

# .. define parameters for 15K
parthres<-1
genth<-15000
lins<-c('p2')
fixed<-as.character(Mikemap[!is.na(Mikemap$max) & Mikemap$max>=0.5 & Mikemap$ann%!in%c('synonymous') & Mikemap$tmax<=15000 & Mikemap$ara%in%lins,]$reanno)

# .. run the script that estimates representativity with the above details
source('repre.R')
write.table(rep, 'Frep.dat', row.names = FALSE, sep = "\t",  col.names = TRUE)	# write table with representativity results (fraction of drivers captured by the DFE)
