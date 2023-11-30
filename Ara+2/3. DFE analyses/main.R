#This master script handles the pipeline for analyzing different aspects of the distribution of fitness effects estimated in folder "2. Fitness estimation"

# .. remove any previous content in the session
rm(list=ls())

# .. start the stopwatch
clock <- Sys.time()

# .. purge out poor-quality fits and outliers presumed to be hitchhikers
source('cleansing.R')

# .. define color scheme for the plots
col2K<-'#cd6e6c'
col15K<-'#6c91cd'

# .. build and plot the DFE
source('DFE.R')

# .. load filtered datasets
Rtable<-read.table(file = "Rfitted_fil.txt", sep="\t", header=TRUE, as.is=TRUE, comment.char = "")
Ttable<-read.table(file = "2Kfitted_fil.txt", sep="\t", header=TRUE, as.is=TRUE, comment.char = "")
Ftable<-read.table(file = "15Kfitted_fil.txt", sep="\t", header=TRUE, as.is=TRUE, comment.char = "")

# .. remove NAs
Rtable<-Rtable[!is.na(Rtable$fitted1),]
Ttable<-Ttable[!is.na(Ttable$fitted1),]
Ftable<-Ftable[!is.na(Ftable$fitted1),]

# .. resample datasets to plot cumulative DFEs with 90% confidence intervals
source('resampling_DFE.R')

# .. zoom-in plots and exponential fit of beneficial tails
source('ben_tail.R')

# .. fits to other several common distributions
source('fittingdis.R')

# .. resample datasets to assess the role of noise in analyses of the shape of DFE
source('resampling.R')

# .. visualize the top beneficial mutations in the ancestor and their effects in the 2K and 15K backgrounds
source('barplots_ben.R')

# .. visualize the change of effects of beneficial mutations accross all backgrounds
source('segments_ben.R')

# .. visualize DFEs for beneficial mutations, both in ancestor and later generations, and compare with the resulting DFE produced by these mutations in the other backgrounds
source('overlaps.R')

# visualize deleterious mutations accross backgrounds
source('barplots_del.R')

# are beneficials found in the actual experiment?
source('metagenomics.R')

# best model to capture realized substitutions
source('lm_3d.R')

# reproducibility among segments from the same locus 
source('comparison_pwise.R')

# .. stop the stopwatch
print(Sys.time()-clock)
