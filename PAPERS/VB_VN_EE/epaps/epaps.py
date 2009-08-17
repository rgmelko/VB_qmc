from numpy import *
from math import log
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
class Evals:
	'''
	A class to hold eigenvalues of the DM
	'''
	def __init__(self):
		'''
		Constructor
		'''
		self.values=zeros(0, float)
	
	def append(self, value):
		'''
		Appends a value
		'''
		self.values=append(self.values, float(value))

	def orderFromHigherToLower(self):
		'''
		Orders the evals in descending order
		'''
		self.values.sort()
		self.values=self.values[::-1]

	def weight(self):
		'''
		Sum all the evals
		'''
		return 1.0-self.values.sum()

	def normalize(self):
		'''
		Normalizes the evals so they become probs
		'''
		norm=sum(self.values)
		self.values=self.values/norm
	
	def entropy(self):
		'''
		Calculates the entropy
		'''
		result=0.0
		for i in range(0,len(self.values)):
			if self.values[i]>1.0e-20:
				result-=self.values[i]*math.log(self.values[i])
		return result

	def __str__(self):
		'''
		Prints out the thing
		'''
		return str(self.values)
	
	def __len__(self):
		'''
		Length of the array, i.e. number of states kept
		'''
		return len(self.values)
# ----------------------------------------------------------------------
# Reads the file into a instance of the Evals class
#line=readLinesFromFile('test.txt') # This line is just for debug
line=readLinesFromFile('tmp.txt') # This line is for the real thing
data=line[0].split() 
#print data # This line is just for debug

allTheEvals=Evals()
for i in range(0,len(data)):
	allTheEvals.append(data[i])

#print allTheEvals # This line is just for debug
allTheEvals.orderFromHigherToLower()
#print allTheEvals # This line is just for debug

# Creates a new instance of Evals with just the largest eigenvalues and
# then prints out in a single line the number of states, the weight, and the
# entropy corresponding to these evals.
highestEvals=Evals()
# You calculate a value ofr entropy, etc... each step state
#step= 1 # This line is just for debug
step=10 # This line is for the real thing
for i in range(0, len(allTheEvals), step):
	highestEvals.values=allTheEvals.values[0:i]
	partialWeight=highestEvals.weight()
	highestEvals.normalize()
	partialEntropy=highestEvals.entropy()
	print str(len(highestEvals))+' '+str(partialWeight)\
			+' '+str(partialEntropy)
allTheWeight=allTheEvals.weight()
allTheEvals.normalize()
print str(len(allTheEvals))+' '+str(allTheWeight)+' '\
		+str(allTheEvals.entropy())
# ----------------------------------------------------------------------
