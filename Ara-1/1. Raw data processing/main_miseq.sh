#!/bin/bash
start=`date +%s`			#to mesure execution time by subtracting the start from the end time
#export PATH=$PATH:/path/to/folder/bwa*	#change accordingly, we used bwa-0.7.17
export PATH=$PATH:/home/alejandro/Software/bwa-0.7.17

sample1F=R0_R1.fastq	#file with forward reads, change accordingly
sample1R=R0_R2.fastq	#file with reverse reads, change accordingly

#DATA EXTRACTION AND PROCESSING PIPELINE	#works as a module, so it can be plugged into a for loop or just copied and pasted as needed

#READ FILE NAMES
nameF=$sample1F				
nameR=$sample1R

#EXTRACT INSERTION REGION
python extract14N_Fold.py $nameF		#script to extract neighbouring regions from forward reads (beware, Python 2.7)
python extract14N_Rold.py $nameR		#script to extract neighbouring regions from reverse reads

#PREPARE FILE NAMES
IFS='_' read -a array <<< "$nameF"		#trick to pass to bash the name of the output files from Python
prefix=${array[0]}
extension="_14Ntoblast_F.fasta"
filenameF=$prefix$extension

IFS='_' read -a array <<< "$nameR"		#trick to pass to bash the name of the output files from Python
extension="_14Ntoblast_R.fasta"
filenameR=$prefix$extension

#LOCAL ALIGNMENT USING THE BWA ALGORITHM
bwa index R606_genoscope.fasta			#reference genome from https://mage.genoscope.cns.fr/microscope/
bwa aln R606_genoscope.fasta $filenameF > bashaln.sai
bwa samse R606_genoscope.fasta bashaln.sai $filenameF > F_bash.sam
bwa aln R606_genoscope.fasta $filenameR > bashaln.sai
bwa samse R606_genoscope.fasta bashaln.sai $filenameR > R_bash.sam

#CLEANING UP THE .SAM FILE (see http://bio-bwa.sourceforge.net/bwa.shtml)
grep "XT:A:U" F_bash.sam > $prefix'_cleanF.sam'
grep "XT:A:U" R_bash.sam > $prefix'_cleanR.sam'
grep "NM:i:0\|NM:i:1" $prefix'_cleanF.sam' > $prefix'_clean_uniqueF.sam'
grep "NM:i:0\|NM:i:1" $prefix'_cleanR.sam' > $prefix'_clean_uniqueR.sam'

#MERGE INFORMATION FROM BOTH READS
python merging.py $prefix $prefix'_clean_uniqueF.sam' $prefix'_clean_uniqueR.sam'
 
#MAP TO ORF
python map.py $prefix $prefix'_merged_results.sam'

#SUMMARY PER ORF
python countperallele.py $prefix $prefix'_results.txt'

#REMOVE INTERMEDIATE FILES
rm *.{sam,pac,sa,bwt,ann,amb,sai}
rm *14Ntoblast_F.fasta
rm *14Ntoblast_R.fasta

end=`date +%s`
runtime=$((end-start))			#to mesure execution time
echo $start, $end, $runtime
