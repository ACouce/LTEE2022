#this script handles issues arising from n-plicated ORFs

# .. load genome maps with and withough n-plicates
index_old<-read.table('parsed_R606genoscope_I.txt', sep="\t", header=FALSE, as.is=TRUE)
index_new<-read.table('parsed_R606genoscope_IU.txt', sep="\t", header=FALSE, as.is=TRUE, comment.char = "") #The default is comment.char = "#"!!

# .. load table of alleles with their abundance accross time points
file<-read.table('Rtable(fuzz).txt', sep="\t", header=TRUE, as.is=TRUE)

newtable<-file
for (i in 1:dim(file)[1]) {			# algorithm to do the cross-correction
    if (!(newtable[i,1]%in%index_new$V1)) {
        for (n in 1:dim(index_new)[1]) {
            if (index_new[n,2]<=newtable[i,2] && index_new[n,3]>=newtable[i,2]) {
                newtable[i,1]<-index_new[n,1]
                break
            }
        }
    }
}

write.table(newtable, file = "Rtable_updated(fuzz).txt", sep = "\t")
