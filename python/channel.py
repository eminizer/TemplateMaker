import copy
from process import MC_Process, Data_Process

#Channel subclass
class Channel(object) :

	#process names
	__mc_process_names = ['fgg','fqq','fbck','fwjets']
	__QCD_process_names = ['fqcd']
	__data_process_names = ['DATA']

	def __init__(self,name,fit_parameter_tuple,include_JEC,include_sss) :
		#find the event topology in this channel
		self.__topology = autoset_topology(name)
		#find the lepton flavor in this channel
		self.__lep_type = autoset_lepton_type_name(name)
		#Set the channel name and charge from the name
		self.__name = name
		self.__charge = autoset_charge(name)
		#Set the channel's region based on the name
		self.__region = autoset_region(name)
		#Make the list of processes that will be in this channel, with all of the wiggles necessary
		self.__process_list = self.__make_process_list__(fit_parameter_tuple,include_JEC,include_sss)

	#Getters/Setters/Adders
	def getName(self) :
		return self.__name
	def getTopology(self) :
		return self.__topology
	def getLepType(self) :
		return self.__lep_type
	def getCharge(self) :
		return self.__charge
	def getRegion(self) :
		return self.__region
	def getProcessList(self) :
		return self.__process_list

	#Internal class functions
	#Make the list of processes in this channel
	def __make_process_list__(self,fit_parameter_tuple,include_JEC,include_sss) :
		channel_plist = []
		#first add the MC processes
		for mpn in self.__mc_process_names :
			channel_plist.append(MC_Process(self.__name+'__'+mpn,fit_parameter_tuple,include_JEC,include_sss))
		#next add the QCD processes
		for qcdpn in self.__QCD_process_names :
			channel_plist.append(QCD_Process(self.__name+'__'+qcdpn,channel_plist,include_JEC,include_sss))
		#add the DATA processes
		for dpn in self.__data_process_names :
			channel_plist.append(Data_Process(self.__name+'__'+dpn))
		return channel_plist

	def __del__(self) :
		pass

	def __str__(self) :
		s = 'Channel: (name:%s, leptype:%s, charge:%i, \n'%(self.__name,self.__lep_type,self.__charge)
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

#returns the topology name based on the channel name
def autoset_topology(name) :
	return name[:2]

#Automatically returns the lepton type name based on the channel name
def autoset_lepton_type_name(name) :
	leptontypename = ''
	if name.split('_')[1].startswith('mu') :
		leptontypename = 'mu'
	elif name.split('_')[1].startswith('el') :
		leptontypename = 'el'
	else :
		print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
		print '!!!!!!   WARNING, LEPTON TYPE NOT RECOGNIZED FROM CHANNEL '+name+'   !!!!!!'
		print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	return leptontypename

#returns the integer lepton charge based on the channel name (0 if charge summed, +/-1 otherwise)
def autoset_charge(name) :
	charge = 0
	if name.split('_')[1].endswith('plus') :
		charge = 1
	elif name.split('_')[1].endswith('minus') :
		charge = -1
	return charge

#returns the region name based on the channel name
def autoset_region(name) :
	ns = name.split('_')
	if len(ns)==3 :
		return ns[2]
	elif len(ns)==4 :
		return ns[2]+'_'+ns[3]
