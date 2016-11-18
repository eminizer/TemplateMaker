class Parameter(object) :

	def __init__(self,name) :
		self.__name = name

	def __del__(self) :
		pass

	def __str__(self) :
		s = self.__name
		return s

class Fit_Parameter(Parameter) :

	def __init__(self,name,value,sigma,fix) :
		Parameter.__init__(self,name)
		self.__value = value
		self.__sigma = sigma
		self.__fix = fix

	def __del__(self) :
		pass

	def __str__(self) :
		s = 'Fit_Parameter: (Parameter:%s, value:%.2f, sigma:%.2f, fix:%s)'%(Parameter.__str__(self),self.__value,self.__sigma,str(self.__fix))
		return s