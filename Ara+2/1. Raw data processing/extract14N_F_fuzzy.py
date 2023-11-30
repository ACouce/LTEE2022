#!/usr/bin/python
#this script extracts the region adjacent to insertion sites. It needs to be called by main.sh.
import re
import regex
import sys
print 'executing', sys.argv[0]
print 'input file:', sys.argv[1]

inputname=str(sys.argv[1])
spl=inputname.split('_')
outputname='_'.join([spl[0],spl[1],'14Ntoblast_F.fasta'])

def main():
	 file=open(inputname)
	 new=open(outputname,"w")
	 cee=0
	 tot=0
	
   	 while (1):	 #algorithm to read the whole file
		line = file.readline()
		tot+=1
		if not line:
			print 'endoffile!'
			break
		if re.match(r'@K0', line):
			line = file.readline()
			tot+=1
			find1=regex.findall("(TAACAGGTTGGATGA){s<=1}", line)	#if match with this one, get the 14bp just before
			find2=regex.findall("(CTGTCTCTTATACAC){s<=1}", line)	#if match with this one, get the 49bp-35bp just before
			
			if len(find1)>0:
				if (line.index(find1[0])-14)<0: #if read is incomplete
					continue
				if tot%10050==0:
					print tot, line[line.index(find1[0])-14:line.index(find1[0])] 
				new.write('>s'+str(tot))
				new.write("\n")
				new.write(line[line.index(find1[0])-14:line.index(find1[0])])
				new.write("\n")
			elif len(find2)>0:
				if (line.index(find2[0])-49)<0: #if read is incomplete
					continue
				if tot%10050==0:
					print tot, line[line.index(find2[0])-49:line.index(find2[0])-35] 
				new.write('>s'+str(tot))
				new.write("\n")
				new.write(line[line.index(find2[0])-49:line.index(find2[0])-35])
				new.write("\n")
			
# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':
    main()
