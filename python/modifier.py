#Modifier class
class Modifier(object) :

	def __init__(self,name) :
		self.__name = name

	def getName(self) :
		return self.__name
	def isJECModifier(self) :
		return isinstance(self,JEC_Modifier)
	def isSSModifier(self) :
		return isinstance(self,Simple_Systematic)
	def isPDFModifier(self) :
		return isinstance(self,PDF_Modifier)

	def __del__(self) :
		pass

	def __str__(self) :
		s = self.__name
		return s

#Fit_Parameter subclass
class Fit_Parameter(Modifier) :

	def __init__(self,name,value,sigma,fix) :
		#Fit parameters have a name, a value, a symmetric uncertainty, and a boolean for whether or not it's fixed
		Modifier.__init__(self,name)
		self.__value = value
		self.__sigma = sigma
		self.__fix = fix

	#Getters/Setters/Adders
	def isFixed(self) :
		return self.__fix
	def getUpValue(self) :
		return self.__value+self.__sigma
	def getDownValue(self) :
		return self.__value-self.__sigma
	def getNomValue(self) :
		return self.__value

	def __del__(self) :
		pass

	def __str__(self) :
		s = 'Fit_Parameter: (Modifier:%s, value:%.2f, sigma:%.2f, fix:%s)'%(Modifier.__str__(self),self.__value,self.__sigma,str(self.__fix))
		return s

#PDF_Modifier subclass
class PDF_Modifier(Modifier) :

	def __init__(self,PDF_branch_name,ipdflambda) :
		#PDF Modifiers have a name (including the PDF reweight branch name and lambda number)
		name = PDF_branch_name+'_lambda_'+str(ipdflambda)
		Modifier.__init__(self,name)
		self.__branch_name = PDF_branch_name
		#PDF Modifiers also have an index in the given branch to find their nominal, up, and down indicies
		self.__nom_index  = 0
		self.__up_index   = 2*ipdflambda-1
		self.__down_index = 2*ipdflambda 

	def getUpIndex(self) :
		return self.__up_index
	def getDownIndex(self) :
		return self.__down_index

	def __del__(self) :
		pass

	def __str__(self) :
		pass

#JEC_Modifier subclass
class JEC_Modifier(Modifier) :

	def __init__(self,name) :
		Modifier.__init__(self,name)

	def __del__(self) :
		pass

	def __str__(self) :
		pass

#Simple_Systematic subclass
class Simple_Systematic(Modifier) :

	def __init__(self,ttree_name,ptree_name) :
		#A Simple systematic has a name for its branch in both the ttrees and the process trees, as well as their +/-1 sigma names
		Modifier.__init__(self,ptree_name)
		self.__ttree_name = ttree_name
		self.__tree_name_up = None; self.__ttree_name_down = None
		if ttree_name!= None :
			self.__tree_name_up = ttree_name+'_hi'
			self.__ttree_name_down = ttree_name+'_low'
		self.__ptree_name = ptree_name
		self.__ptree_name_up = ptree_name+'_up'
		self.__ptree_name_down = ptree_name+'_down'

	def getPTreeNameUp(self) :
		return self.__ptree_name_up
	def getPTreeNameDown(self) :
		return self.__ptree_name_down

	def __del__(self) :
		pass

	def __str__(self) :
		pass