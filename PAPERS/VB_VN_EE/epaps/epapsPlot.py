from numpy import *
import pylab

# ----------------------------------------------------------------------
def readLinesFromFile(file_name):
    ''' Try to open the file, and read in all the lines, if it doesn't 
    exist, exits'''
    try:
        in_file = open(file_name)
        lines = in_file.readlines()
        in_file.close()
    except:
        print 'File ' + str(file_name) + ' not found'
        sys.exit()

    return lines
# ----------------------------------------------------------------------
numberOfStates=zeros(0,float)
partialWeight=zeros(0,float)
partialEntropy=zeros(0,float)

lines=readLinesFromFile('epapsN7.dat') 
for line in lines:
	splittedLine=line.split()
	numberOfStates=append(numberOfStates, float(splittedLine[0]))
	partialWeight=append(partialWeight, float(splittedLine[1]))
	partialEntropy=append(partialEntropy, float(splittedLine[2]))

errorInTheEntropy=partialEntropy[len(partialEntropy)-1]-partialEntropy

truncation=int(180)
truncNumberOfStates=numberOfStates[:truncation]
truncPartialWeight=partialWeight[:truncation]
truncPartialEntropy=partialEntropy[:truncation]
truncErrorInTheEntropy=errorInTheEntropy[:truncation]

pylab.figure()
pylab.semilogy(truncNumberOfStates, truncErrorInTheEntropy, 'bo-')
pylab.semilogy(truncNumberOfStates, truncPartialWeight, 'ro-')
pylab.legend(('$S^{vN}_{m}$','truncation error'), 'upper right')
#pylab.plot(numberOfStates, partialEntropy, 'bo-')
#pylab.plot(numberOfStates, partialWeight, 'bo-')
pylab.xlabel('$m$')
#pylab.ylabel('$S_{m}$')
pylab.axis([10, max(truncNumberOfStates), min(truncPartialWeight), 
	max(truncErrorInTheEntropy)])
#pylab.semilogx(numberOfStates, partialEntropy, 'bo-')
#pylab.loglog(numberOfStates, partialWeight, 'bo-')
pylab.savefig('temp.png')

