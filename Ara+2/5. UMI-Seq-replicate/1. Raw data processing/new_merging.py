#!/usr/bin/python
#this script extracts the region adjacent to insertion sites. It needs to be called by main_umiseq.sh.
import re
import sys
import mmap	# important library to enable memory mapping
print ('executing', sys.argv[0], 'with input files:', sys.argv[1], sys.argv[2], sys.argv[3])

outputname=str(sys.argv[1]) # this is the 'merged.fasta'
inputnameF=str(sys.argv[2]) # this is one .sam
inputnameR=str(sys.argv[3]) # this is the other .sam

def main():

	new=open(outputname,"w")
	file1= open('both.sam', "r")
	fileF=open(inputnameF, "r")
	fileR=open(inputnameR, "r") 
	tot=0
	exception=0
	 
	mm_query = mmap.mmap(file1.fileno(), length=0,access=mmap.ACCESS_READ,offset=0)	# trick to speed up file I/O via memory mapping
	mm_objectF = mmap.mmap(fileF.fileno(), length=0,access=mmap.ACCESS_READ,offset=0)
	mm_objectR = mmap.mmap(fileR.fileno(), length=0,access=mmap.ACCESS_READ,offset=0)
	 
	while (1):		
		#CHECK INDEX WITH ONLY READS PRESENT IN BOTH FILES
		line1 = mm_query.readline()    
		if not line1:
			file1.close()
			break
		dec=line1.decode()
		spl=dec.split()
		label=spl[0]
	
		# LOOK UP FIRST FILE
		while (1):
			
			line=mm_objectF.readline()
			if not line:
				mm_objectF = mmap.mmap(fileF.fileno(), length=0,access=mmap.ACCESS_READ,offset=0)
				line=mm_objectF.readline()
			if label in line.decode():
				lineF=line.decode()
				spl1=lineF.split()	             
				break

		# LOOK UP SECOND FILE			
		while (1):
			line=mm_objectR.readline()
			if not line:
				mm_objectR = mmap.mmap(fileR.fileno(), length=0,access=mmap.ACCESS_READ,offset=0)
				line=mm_objectR.readline()		
			if label in line.decode():
				lineR=line.decode()
				spl2=lineR.split()	             
				break
		
		if spl1[0]!=spl2[0]:
			print(tot,label)
			print(spl1)
			print(spl2)			

		#COMPARE BOTH READS and KEEP ONE
		if (spl1[3] == spl2[3] and spl1[5] != '15M'): 	#clausule 1: both reads are null
			continue
		elif  (spl1[3] == spl2[3]): 			#clausule 2: both sites coincide
			new.write(str(lineF))
			tot+=1
			if tot%100000==0:
				print ('both',tot, exception)
			continue
		elif (spl1[5] != '15M' and spl2[5] == '15M'):	 #clausule 3: F read is bad
			spl2[1]=str(16-int(spl2[1]))		#"convert" it to a F strand
			lineR='\t'.join(spl2) 			#regenerate edited line
			new.write(str(lineR))
			new.write("\n")			
			tot+=1
			continue
		elif (spl2[5] != '15M' and spl1[5] == '15M'): 	#clausule 4: R read is bad
			new.write(str(lineF))
			tot+=1
			continue
		elif (spl2[12] != 'NM:i:0' and spl1[12] == 'NM:i:0'): #a mutation in one read happens to also find a unique 15M elsewhere! (but decidable on the basis of NM:i:0)
			new.write(str(lineF))
			tot+=1
			continue
		elif (spl2[12] == 'NM:i:0' and spl1[12] != 'NM:i:0'): #a mutation in one read happens to also find a unique 15M elsewhere! (but decidable on the basis of NM:i:0)
			spl2[1]=str(16-int(spl2[1]))		#"convert" it to a F strand
			lineR='\t'.join(spl2) 			#regenerate edited line
			new.write(str(lineR))
			new.write("\n")	
			tot+=1
			continue

		else:
			exception+=1
	file1.close()
	fileF.close()
	fileR.close()
	
#LET's ADD THE SINGLETONS
	file1=open('UlabelsF.sam', "r")
	fileF=open(inputnameF, "r")

	exception=0
	 
	mm_query = mmap.mmap(file1.fileno(), length=0,access=mmap.ACCESS_READ,offset=0)	# trick to speed up file I/O via memory mapping
	mm_objectF = mmap.mmap(fileF.fileno(), length=0,access=mmap.ACCESS_READ,offset=0)
	
	while(1):
	
		line1 = mm_query.readline()    
		if not line1:
			print ('endofsearch!')
			file1.close()
			break
		dec=line1.decode()
		spl=dec.split()
		label=spl[0]
	
		while (1):
			line=mm_objectF.readline()
			if label in line.decode():
				lineF=line.decode()
				spl1=lineF.split()	             
				break

		if (spl1[5] == '15M'):
			new.write(str(lineF))
			tot+=1
			if tot%100000==0:
				print ('F', tot, exception)		

	file1.close()
	fileF.close()
	
	#REVERSE ONES
	file1=open('UlabelsR.sam', "r")
	fileR=open(inputnameR, "r")

	exception=0
	 
	mm_query = mmap.mmap(file1.fileno(), length=0,access=mmap.ACCESS_READ,offset=0)
	mm_objectR = mmap.mmap(fileR.fileno(), length=0,access=mmap.ACCESS_READ,offset=0)
	
	while(1):
	
		line1 = mm_query.readline()    
		if not line1:
			print ('endofsearch!')
			file1.close()
			break
		dec=line1.decode()
		spl=dec.split()
		label=spl[0]
	
		while (1):
			line=mm_objectR.readline()
			if label in line.decode():
				lineR=line.decode()
				spl1=lineR.split()	             
				break

		if (spl1[5] == '15M'):			
			spl1[1]=str(16-int(spl1[1]))	#"convert" it to a F strand
			lineR='\t'.join(spl1)		#regenerate edited line
			new.write(str(lineR))
			new.write("\n")				
			tot+=1
			if tot%100000==0:
				print ('R', tot, exception)	
	file1.close()
	fileR.close()

# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':
    main()
