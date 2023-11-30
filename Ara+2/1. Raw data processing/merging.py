#!/usr/bin/python
#this script merges .sam files from forward and reverse reads. It needs to be called by main.sh.
import re
import sys
print 'executing', sys.argv[0], 'with input files:', sys.argv[2], sys.argv[3]

outputname='_'.join([str(sys.argv[1]),'merged_results.sam'])
inputnameF=str(sys.argv[2])
inputnameR=str(sys.argv[3])

def main():

	 new=open(outputname,"w")
	 file1=open(inputnameF)
	 file2=open(inputnameR)
	 tot=0
	 exception=0
	 nlines1=0
	 nlines2=0
	 while (1):	#retrieve the highest cardinal number (F reads are numbered in sequential order s2, s6, s10, s14, etc...)
		line1 = file1.readline()
		if not line1:
			file1.close()
			break
		spl1=line1.split()
		max1=spl1[0]
	 while (1):	#retrieve the highest cardinal number (R reads are too numbered in sequential order s2, s6, s10, s14, etc...)
		line2 = file2.readline()
		if not line2:
			file2.close()
			break
		spl2=line2.split()
		max2=spl2[0]

	 print max1, max2	

	 flag1=0
	 flag2=0
	 flag_end1=0
	 flag_end2=0
	 s_ini=-2
	
	 #ALGORITHM TO COMPARE BOTH FILES SEARCHING IN BOTH THE SERIES s2, s6, s10, s14, ...
	 file1=open(inputnameF)
	 file2=open(inputnameR)
   	 while (1):
		s_ini=s_ini+4
		query = 's'+str(s_ini)
	
		while (1): #find the line with the appropriate query label (s2, s6, etc) in file F

			if flag_end1==1:
				spl1[3]=0
				break
			if flag1==1:
				spl1=line1.split()
				if (spl1[0] == query):
					flag1=0
					break
				else:
					spl1[3]=0		#if query is absent flag as 0
					break

			line1 = file1.readline() #read line
			if not line1:
				file1.close()
				file1=open(inputnameF)	
				spl1[3]=0		#if query is absent flag as 0
				flag_end=1
				break
				
			spl1=line1.split() #read query label
			if (int(re.search(r'\d+', spl1[0]).group()) > s_ini): # check if we are ahead
				flag1=1
				spl1[3]=0		#if query is absent flag as 0
				break
			nlines1+=1
			if (spl1[0] == query):
				break

		while (1): #find the line with the appropriate query label (s2, s6, etc) in file R

			if flag_end2==1:
				spl2[3]=0
				break
			if flag2==1:
				spl2=line2.split()
				if (spl2[0] == query):
					flag2=0
					break
				else:
					spl2[3]=0		#if query is absent flag as 0
					break
			line2 = file2.readline()
			if not line2:
				file2.close()
				file2=open(inputnameR)
				spl2[3]=0		#if query is absent flag as 0
				flag_end2=1
				break
			spl2=line2.split()
			if (int(re.search(r'\d+', spl2[0]).group()) > s_ini):
				flag2=1
				spl2[3]=0		#if query is absent flag as 0
				break
			nlines2+=1
			if (spl2[0] == query):
				break

		#COMPARE BOTH READS and KEEP ONE
		if (int(re.search(r'\d+', max1).group()) < s_ini) and (int(re.search(r'\d+', max2).group()) < s_ini):	
			break
		if (spl1[3] == spl2[3] and spl1[3] == 0): #clausule 1: both reads are null/absent
			continue
		elif  (spl1[3] == spl2[3] and spl1[0] == query): #clausule 2: both reads are present and coincide
			new.write(str(line1))
			tot+=1
			if tot%5000==0:
				print tot, s_ini
			continue
		elif (spl1[3] == 0 and spl2[0] == query): #clausule 3: F read is null/absent and read R is present
			new.write(str(line2))
			tot+=1
			if tot%5000==0:
				print tot, s_ini
			continue
		elif (spl2[3] == 0 and spl1[0] == query): #clausule 4: R read is null/absent and read F is present
			new.write(str(line1))
			tot+=1
			if tot%5000==0:
				print tot, s_ini
			continue
		elif (abs(int(spl1[3])-int(spl2[3]))<=1): #clausule 5: F and R reads differ in a +1 or -1 frameshift, ans both give unique mapping.
			new.write(str(line1))
			tot+=1
			if tot%5000==0:
				print tot, s_ini
			continue
		else:
			print 'exception', exception, tot, s_ini
			print line1
			print line2
			exception+=1
			continue
	 file1.close()
	 file2.close()
# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':
    main()
