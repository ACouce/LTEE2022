#!/usr/bin/python
#this script counts the number of insertions observed in the same site. It needs to be called by main.sh.
import sys
import numpy
import scipy.stats as ss
print 'executing', sys.argv[0], 'with input files:', sys.argv[1], sys.argv[2]

inputname='_'.join([str(sys.argv[1]),'results.txt'])
outputname=str(sys.argv[2])

def main():

    collection1=[]
    collection2=[]
    n=0
    count=0
    file=open(inputname)
    while (1):		#algorithm to create iterable collections
	    line = file.readline()
	    if not line:
                print 'endoffile'
                break
	    spl=line.split()
	    gene=spl[1]
	    pos=spl[2]
	    n+=1
	    collection1.append(gene)
	    collection2.append(pos)	    
    
    from collections import Counter
    freqs = Counter(collection2)
    print len(collection2)
    print freqs
    new=open(outputname,"w")

    for p,m in freqs.iteritems():
	indexes=[]
	for i,x in enumerate(collection2):
		if x == p:
			indexes.append(i)

	T=[collection1[i] for i in indexes] 
	count+=1
	if count%10050==0:
	    print count
	new.write(str(collection1[i])+"\t"+str(p)+"\t"+str(m)+"\t")
	new.write("\n")
    new.close()

# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':
    main()
