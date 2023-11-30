# .. this script estimates how well either real or simulated DFEs capture the drivers of adaptation in the LTEE, as well as the fraction of drivers that can be understood as loss-of-function mutations

rep<-matrix(0, length(gens), replicates)
LOF<-matrix(0, length(gens), replicates)	# uncomment here and within 'representativeness_real.R' for analyses fo prevalence of loss-of-function mutations 


for (r in 1:replicates) {
	print(r)
	t<-0
	for (time in gens) {
		t<-t+1
		if (time<genth) {next}
		# identify the main drivers of adaptation up to this generation (non-synonymous mutations that reached >50% in frequency)
		drivfull<-Mikemap[!is.na(Mikemap$max) & Mikemap$max>=0.5 & Mikemap$ann%!in%c('synonymous') & Mikemap$reanno%!in%fixed & Mikemap$tmax<=time & Mikemap$ara%in%lins,]
		drivers<-as.character(drivfull$reanno)		
		tab<-table(drivers)				# count in how many independent populations these drivers arose
		fix<-names(tab[tab>=parthres])			# keep those appearing in at least 2 populations (only for the ancestor)
		fix<-paste(as.character(fix), collapse=",")	# manipulations to facilitate matching
		fix<-unique(unlist(strsplit(as.character(fix), split=",")))
		nfix<-length(fix)				# total number of drivers

		if (label=='sims') {					# if processing the simulated data
			bens<-unique(Ntabs[Ntabs$rep==r,]$name)		# load beneficial mutations from the DFE of replicate r
			bens<-paste0('\\b', bens, '\\b', collapse = "|")	# manipulations to facilitate matching via grepl			
		} else {						# if processing the real data
			bens<-unique(Ntabs$name)				# load beneficial mutations from the real DFEs
			bens<-paste0('\\b', bens, '\\b', collapse = "|")	# manipulations to facilitate matching via grepl
		}
		
		count<-sum(grepl(bens,fix))			# counting matches
		rep[t,r]<-count/nfix				# store results
		
		if (time %in% c(1500, 2000,3000,4000,6000,7000,seq(5000,50000,5000))) {	# to monitor progress
			print(time)		
		}
				
		# uncomment all below and within 'representativeness_real.R' for analyses fo prevalence of loss-of-function mutations 
		countLOF<-0	
		drivfull<-Mikemap[!is.na(Mikemap$max) & Mikemap$max>=0.5 & Mikemap$ann%!in%c('synonymous') & Mikemap$reanno%!in%fixed & Mikemap$tmax<=time+500 & Mikemap$ara%in%lins,]
		for (d in fix) {
			if (sum(drivfull[drivfull$reanno%in%d,]$ann%in%c('frameshift','nonsense','sv'))) {	#presumed loss-of-function mutations	
				countLOF<-countLOF+1
			}
		}
		LOF[t,r]<-countLOF/nfix	

		}
}
