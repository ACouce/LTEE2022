#!/bin/bash
start=`date +%s`			#to mesure execution time by subtracting the start from the end time

#export PATH=$PATH:/path/to/folder/bwa*	#change accordingly, we used bwa-0.7.17
export PATH=$PATH:/home/alejandro/Software/bwa-0.7.17

prefix="ara2-2K_" 		#change accordingly

#DATA EXTRACTION AND PROCESSING PIPELINE	#works as a module, so it can be plugged into a for loop or just copied and pasted as needed

#READ FILE NAMES
nameF=$sample1F				
nameR=$sample1R

gunzip $nameF'.gz'
gunzip $nameR'.gz'

#EXTRACT INSERTION REGION
python3 extract14N_F_fuzzy.py $nameF		#script to extract neighbouring regions from forward reads (beware, Python 3)
python3 extract14N_R_fuzzy.py $nameR		#script to extract neighbouring regions from reverse reads

#PREPARE FILE NAMES
IFS='_' read -a array <<< "$nameF"		#trick to pass to bash the name of the output files from Python
samplename=$prefix${array[1]}
extension="_14Ntoblast_F.fasta"
filenameF=$samplename$extension

IFS='_' read -a array <<< "$nameR"		#trick to pass to bash the name of the output files from Python
samplename=$prefix${array[1]}
extension="_14Ntoblast_R.fasta"
filenameR=$samplename$extension

#LOCAL ALIGNMENT USING THE BWA ALGORITHM
bwa index R606_genoscope.fasta			#reference genome from https://mage.genoscope.cns.fr/microscope/
bwa aln R606_genoscope.fasta $filenameF > bashaln.sai
bwa samse R606_genoscope.fasta bashaln.sai $filenameF > F_bash.sam
bwa aln R606_genoscope.fasta $filenameR > bashaln.sai
bwa samse R606_genoscope.fasta bashaln.sai $filenameR > R_bash.sam

#CLEANING UP THE .SAM FILE (see http://bio-bwa.sourceforge.net/bwa.shtml)
grep "XT:A:U" F_bash.sam > $samplename'_cleanF.sam' 
grep "XT:A:U" R_bash.sam > $samplename'_cleanR.sam'
grep "NM:i:0\|NM:i:1" $samplename'_cleanF.sam' > $samplename'_clean_uniqueF.sam'
grep "NM:i:0\|NM:i:1" $samplename'_cleanR.sam' > $samplename'_clean_uniqueR.sam'

#MERGE INFORMATION FROM BOTH READS
###first, just labels
grep -Po '^@A\s*\S+' $samplename'_clean_uniqueF.sam' > labelsF.sam
grep -Po '^@A\s*\S+' $samplename'_clean_uniqueR.sam' > labelsR.sam
##second, compile index of cases with two reads
grep -Fxf labelsF.sam labelsR.sam > both.sam
##third, compile index of cases with just one read
grep -Fxv -f both.sam labelsF.sam > UlabelsF.sam
grep -Fxv -f both.sam labelsR.sam > UlabelsR.sam
#https://stackoverflow.com/questions/1016244/unix-command-to-find-string-set-intersections-or-outliers

python3 new_merging.py merged.sam $samplename'_clean_uniqueF.sam' $samplename'_clean_uniqueR.sam' #this increases quality in the overlapped dataset, and reduces issues with counting twice the same

##MAP TO ORF (#countperallele best handled now via R)
finalname=$prefix${array[2]}
python3 new_map.py merged.sam $finalname'_merged.txt'

#REMOVE INTERMEDIATE FILES
rm *.{sam,pac,sa,bwt,ann,amb,sai}
rm *14Ntoblast_F.fasta
rm *14Ntoblast_R.fasta

gzip *.fastq

end=`date +%s`
runtime=$((end-start))			#to mesure execution time
echo $start, $end, $runtime
