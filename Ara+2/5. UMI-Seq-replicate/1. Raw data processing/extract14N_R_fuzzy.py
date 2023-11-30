#!/usr/bin/python
#this script extracts the region adjacent to insertion sites. It needs to be called by main_umiseq.sh.
import re
import regex
import sys
print ('executing', sys.argv[0])
print ('input file:', sys.argv[1])

def rev_comp(st):					#useful function to compute the reverse complement
    nn = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(nn[n] for n in reversed(st))

inputname=str(sys.argv[1])
spl=inputname.split('_')
outputname='_'.join([spl[0],spl[1],'14Ntoblast_R.fasta'])

def main():
	file=open(inputname)
	new=open(outputname,"w")
	cee=0
	tot=0

	while (1):			#algorithm to read the whole file
		line = file.readline()
		if not line:
			print ('endoffile!')
			break
		if re.match(r'@A0', line):
			label = line.split()[0] #if no argument splits on whitespace
			line = file.readline()          
			find1=regex.findall("(TAACAGGTTGGCTGATAAGTCCCCCAA){s<=1}", str.strip(line))	#str.strip to remove \n
			
			if len(find1)>0:

				if len(line[line.index(find1[0])-15:line.index(find1[0])])==0: 		#we dont have to worry here whether read is incomplete
					continue
				
				line=line.split('\n')[0] 						#to prevent jumping line in cases of short umis
				
				umi=line[line.index(find1[0])+27:line.index(find1[0])+27+5]		#retrieve UMI
				
				if re.search("N", umi): 							#re.match looks only at the begining of the string
					continue
				if(len(umi)<5):
					continue
				
				umi=rev_comp(umi)							#get the reverse complement
				
				tot+=1
				if tot%10050==0:
					print (tot, line[line.index(find1[0])-15:line.index(find1[0])])
				
				new.write('>'+str(label)+'-UMI:'+umi)
				new.write("\n")
				new.write(line[line.index(find1[0])-15:line.index(find1[0])])
				new.write("\n")
			
# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':
    main()
