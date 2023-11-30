#!/usr/bin/python
#this script extracts the region adjacent to insertion sites. It needs to be called by main.sh.
import re
import regex
import sys
print 'executing', sys.argv[0]
print 'input file:', sys.argv[1]

inputname=str(sys.argv[1])
spl=inputname.split('_')
outputname='_'.join([spl[0],spl[1],'14Ntoblast_R.fasta'])

endF='CTGTTA' #if this one AND secondF are seen, get the 14bp just afterwards

length1=15
length2=len(endF)

def main():
	 file=open(inputname)
	 new=open(outputname,"w")
	 cee=0
	 tot=0
	 #algorithm to read the whole file
   	 while (1):
		line = file.readline()
		tot+=1
		if not line:
			print 'endoffile!'
			break
		if re.match(r'@K0', line):
			line = file.readline()
			tot+=1
			find1=regex.findall("(TCATCCAACCTGTTA){s<=1}", line)	#if match to this one, get the 14bp just afterwards
			find2=regex.findall("(CTGTCTCTTATACAC){s<=1}", line)	#if match to this one, get the 25bp just before
			
			if len(find1)>0:
				#if (line.index(firstF)+length1+14)<line.index(firstF): #we dont have to worry here whether read is incomplete
				#	continue
				if tot%10050==0:
					print tot, line[line.index(find1[0])+length1:line.index(find1[0])+length1+14] 
				new.write('>s'+str(tot))
				new.write("\n")
				new.write(line[line.index(find1[0])+length1:line.index(find1[0])+length1+14])
				new.write("\n")
			elif len(find2)>0:
				if (line.index(find2[0])-25)<0: 
					continue
				window = line[line.index(find2[0])-25:line.index(find2[0])]
				if endF in window:
					if tot%10050==0:
						print tot, window[window.index(endF)+length2:window.index(endF)+length2+14]
					new.write('>s'+str(tot))
					new.write("\n")
					new.write(window[window.index(endF)+length2:window.index(endF)+length2+14])
					new.write("\n")
			
# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':
    main()
