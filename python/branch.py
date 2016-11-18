from array import array

#Branch class
class Branch(object) :

	def __init__(self,ttree_name,ptree_name,typestring,initial_value) :
		#A branch has a name each for its use in ttrees and process trees
		self.__ttree_name = ttree_name
		self.__ptree_name = ptree_name
		#Based on the initial value and variable type/array length, create arrays for ttrees and ptree use
		self.__arraylength = 1
		self.__pythonvtype = typestring
		tss = typestring.split('/')
		if len(tss)==2 :
			self.__arraylength = int(tss[0])
			self.__pythonvtype = tss[1]
		self.__rootvtype = convertToRootType(self.__pythonvtype)
		self.__initial_value = initial_value
		self.__ttree_array = None
		self.__ptree_array = None
		if ttree_name!=None :
			self.__ttree_array = array(self.__pythonvtype,self.__arraylength*[initial_value])
		if ptree_name!=None :
			self.__ptree_array = array(self.__pythonvtype,self.__arraylength*[initial_value])

	#Getters/Setters/Adders
	def getTTreeName(self) :
		return self.__ttree_name
	def getTTreeArray(self) :
		return self.__ttree_array
	def getPTreeName(self) :
		return self.__ptree_name
	def getPTreeArray(self) :
		return self.__ptree_array
	def getTTreeValue(self) :
		if self.__arraylength!=1 :
			print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
			print '!!!!!!   WARNING, METHOD SHOULD NOT BE USED TO GET VALUE FOR BRANCH '+self.__ptree_name+'   !!!!!!'
			print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
			return self.__initial_value
		else :
			return self.__ttree_array[0]
	def getPTreeValue(self) :
		if self.__arraylength!=1 :
			print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
			print '!!!!!!   WARNING, METHOD SHOULD NOT BE USED TO GET VALUE FOR BRANCH '+self.__ptree_name+'   !!!!!!'
			print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
			return self.__initial_value
		else :
			return self.__ptree_array[0]
	def getPTreeValueAtIndex(self,index) :
		return self.__ptree_array[index]
	#return the "leaf list" string of the branch for process trees
	def getPTreeLeafList(self) :
		s = self.__ptree_name
		if self.__arraylength>1 :
			s+='['+str(self.__arraylength)+']'
		s+='/'+self.__rootvtype
		return s

	def __del__(self) :
		pass

	def __str__(self) :
		pass

#return the ROOT version of the variable type given the python name
def convertToRootType(ptype) :
	if ptype==None :
		return ptype
	elif ptype=='d' :
		return 'D'
	elif ptype=='i' :
		return 'I'
	elif ptype=='I' :
		return 'i'
	elif ptype=='f' :
		return 'F'
	else :
		print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
		print '!!!!!!   WARNING, VARIABLE TYPE NOT RECOGNIZED FROM PYTHON TYPE '+ptype+'   !!!!!!'
		print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'