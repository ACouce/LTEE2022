#!/usr/bin/python
#this script maps insertion sites to the list of loci of the reference genome. It needs to be called by main.sh.
import re
import sys
print 'executing', sys.argv[0], 'with input file:', sys.argv[2]

outputname='_'.join([str(sys.argv[1]),'results.txt']) 
inputname=str(sys.argv[2])

cutoff=0
def main():
	 file=open(inputname)
	 ref=open("parsed_R606genoscope_IU.txt")	#the list of loci of the reference genome
	 new=open(outputname,"w")
	 cee=0
	 tot=0
	 itg=0
   	 while (1):				#algorithm to read the whole file
		line = file.readline()
		tot+=1
		flag=0
		overl_flag=0
		counter=0
		if not line:
			print 'endoffile!'
			break
		spl=line.split()

		if str(spl[5])!='14M':		#not reliable
			print '14m'
			continue
		if int(spl[1])==4:		#not aligned
			print '4'
			continue
		if int(spl[1])==0:		#forward stranded
			realpos=int(spl[3])+14
		if int(spl[1])==16:		#reverse stranded
			realpos=int(spl[3])
		while (1):
			refline = ref.readline()
			cee+=1
			if not refline:
				ref.close()
				ref=open("parsed_R606genoscope_IU.txt")
				if flag==0:
					itg+=1
				break
			refspl=refline.split()
			if refspl[3]=='F':
				if 'IR' in refspl[0]:  #in the parsing of the .gbk file from GenBank, IRs were arbitrarily assigned a F orientation
					cutoff=0
				genlen=int(refspl[2])-int(refspl[1])
				if realpos>=int(refspl[1]) and realpos<=(int(refspl[2])-int(cutoff*genlen)):
					new.write(str(refspl[1])+"\t"+ str(refspl[0])+"\t"+ str(realpos))
					new.write("\n")
					flag=1
					overl_flag=1
					if tot%5000==0:
						print tot, itg
				cutoff=0
			elif refspl[3]=='C':
				genlen=int(refspl[2])-int(refspl[1]) #the parsed genome file displays the positions in 5->3!
				if realpos>=(int(refspl[1])+int(cutoff*genlen)) and realpos<=int(refspl[2]):
					new.write(str(refspl[1])+"\t"+ str(refspl[0])+"\t"+ str(realpos))
					new.write("\n")
					flag=1
					overl_flag=1
					if tot%5000==0:
						print tot, itg
			if overl_flag==1: #this trick would allow recording overlapping genes by ensuring that the next two ORFs are always tested (because there are even triple overlaps)
				counter+=1
				if counter>2:
					ref.close()
					ref=open("parsed_R606genoscope_IU.txt")
					break 
# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':
    main()
