#!/usr/bin/python
#this script extracts the region adjacent to insertion sites. It needs to be called by main_umiseq.sh.
import re
import sys
import mmap	# important library to enable memory mapping
import pandas as pd
import bisect

print('executing', sys.argv[0], 'with input file:', sys.argv[2])

outputname=str(sys.argv[2])

cutoff=0
def main():
	file1=open('merged.sam')
	ref = pd.read_csv("parsed_R606genoscope_I.txt", sep='\t', header=None)
	new=open(outputname,"w")

	tot=0

	mm_query = mmap.mmap(file1.fileno(), length=0,access=mmap.ACCESS_READ,offset=0)	# trick to speed up file I/O via memory mapping
	
	while(1):
	
		line1 = mm_query.readline()    
		if not line1:
			file1.close()
			break
		dec=line1.decode()
		spl=dec.split()
		tot+=1
		if tot%100000==0:
			print(tot)
		if int(spl[1])==0: 							#forward stranded
			realpos=int(spl[3])
		if int(spl[1])==16: 							#reverse stranded
			realpos=int(spl[3])+15		
		
		umi = spl[0].split("-UMI:")[1]

		orf=ref.iloc[bisect.bisect_left(ref[2], realpos)][0]			 #bisect_left returns the left-most index to insert the given element, while keeping order
		new.write(str(orf)+"\t"+ str(realpos)+"\t"+ str(umi))
		new.write("\n")				

	file1.close()

# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':
    main()
