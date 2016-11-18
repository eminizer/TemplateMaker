import copy
from process import MC_Process, NTMJ_Process, Data_Process

#Lepton_Type class
class Lepton_Type(object) :
	
	def __init__(self,name) :
		self.__name = name #lepton flavor name

	#Getters/Setters/Adders
	def getLepType(self) :
		return self.__name

	def __del__(self) :
		pass

	def __str__(self) :
		s = self.__name
		return s

#Channel subclass
class Channel(Lepton_Type) :

	#process names
	__mc_process_names = ['fgg','fqq','fbck']
	__ntmj_process_names = ['fntmj']
	__data_process_names = ['DATA']

	def __init__(self,name,fit_parameter_tuple,include_PDF,include_JEC,include_sss) :
		#find the lepton flavor in this channel
		leptontypename = autoset_lepton_type_name(name)
		#initialize the Lepton_Type
		Lepton_Type.__init__(self,leptontypename)
		#Set the channel name and charge from the name
		self.__name = name
		self.__charge = autoset_charge(name)
		#Make the list of processes that will be in this channel, with all of the wiggles necessary
		self.__process_list = self.__make_process_list__(fit_parameter_tuple,include_PDF,include_JEC,include_sss)

	#Getters/Setters/Adders
	def getName(self) :
		return self.__name
	def getCharge(self) :
		return self.__charge
	def getProcessList(self) :
		return self.__process_list

	#Internal class functions
	#Make the list of processes in this channel
	def __make_process_list__(self,fit_parameter_tuple,include_PDF,include_JEC,include_sss) :
		channel_plist = []
		#first add the MC processes
		for mpn in self.__mc_process_names :
			channel_plist.append(MC_Process(self.__name+'__'+mpn,fit_parameter_tuple,include_PDF,include_JEC,include_sss))
		#copy the list so we can associate it with the NTMJ process
		mc_plist = copy.copy(channel_plist)
		#add the NTMJ processes
		for npn in self.__ntmj_process_names :
			channel_plist.append(NTMJ_Process(self.__name+'__'+npn,fit_parameter_tuple,mc_plist,include_PDF,include_JEC,include_sss))
		#add the DATA processes
		for dpn in self.__data_process_names :
			channel_plist.append(Data_Process(self.__name+'__'+dpn))
		return channel_plist

	def __del__(self) :
		pass

	def __str__(self) :
		s = 'Channel: (name:%s, leptype:%s, charge:%i, \n'%(self.__name,Lepton_Type.__str__(self),self.__charge)
		s+= '			process_list = ['
		for i in range(len(self.__process_list)) :
			proc = self.__process_list[i]
			if i==0 :
				s+=proc.__str__()
			else :
				s+='				'+proc.__str__()
			if i==len(self.__process_list)-1 :
				s+='])\n'
			else :
				s+=',\n'
		return s

#Automatically returns the lepton type name based on the channel name
def autoset_lepton_type_name(name) :
	leptontypename = ''
	if name.startswith('mu') :
		leptontypename = 'mu'
	elif name.startswith('el') :
		leptontypename = 'el'
	else :
		print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
		print '!!!!!!   WARNING, LEPTON TYPE NOT RECOGNIZED FROM CHANNEL '+name+'   !!!!!!'
		print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	return leptontypename

#returns the integer lepton charge based on the channel name (0 if charge summed, +/-1 otherwise)
def autoset_charge(name) :
	charge = 0
	if name.endswith('plus') :
		charge = 1
	elif name.endswith('minus') :
		charge = -1
	return charge