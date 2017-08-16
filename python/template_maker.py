import copy
from math import *
from array import array
from ROOT import TFile, TTree, TF1, TCanvas, TLegend, TGraphErrors, TVirtualFitter, kBlue
from channel import Channel
from branch import Branch

#Template_Group class
class Template_Group(object) :

	#Luminosities
	__LUMINOSITY_BTOF_MU = 19690.184
	__LUMINOSITY_BTOF_EL = 19171.010
	__LUMINOSITY_GH_MU   = 16226.452
	__LUMINOSITY_GH_EL   = 16214.862
	__LUMI_MIN_BIAS = 0.026

	def __init__(self,fit_parameter_tuple,sum_charges,include_mu,include_el,topology_list,include_JEC,include_sss) :
		#switches for which templates to auto-create
		self.__include_JEC = include_JEC
		self.__include_mu = include_mu
		self.__include_el = include_el
		self.__include_sss = include_sss
		#make the list of channels, which will have their processes and templates automatically added
		self.__channel_list = make_channel_list(fit_parameter_tuple,include_mu,include_el,sum_charges,topology_list,include_JEC,include_sss)
		self.__aux_obj_list = []
		#print 'Created new template group:\n'+self.__str__()

	#given a file path, add the events to the file to the processes
	def add_file_to_processes(self,ttree_file_path) :
		#skip the JEC files if we don't need them
		if (not self.__include_JEC) and (ttree_file_path.find('JES')!=-1 or ttree_file_path.find('JER')!=-1) :
			return
		#skip the lepton-separated data files if we don't need them
		if (not self.__include_mu and ttree_file_path.find('SingleMu')!=-1) or (not self.__include_el and ttree_file_path.find('SingleEl')!=-1) :
			return
		#get the list of contributing processes
		ps_for_file = []
		for channel in self.__channel_list :
			for process in channel.getProcessList() :
				if get_contrib_weight(ttree_file_path,process)!=0. :
					ps_for_file.append(process)
		if len(ps_for_file)==0 :
			return
#		s = 'ps_for_file = [' #DEBUG
#		for p in ps_for_file : #DEBUG
#			s+=p.getName()+',' #DEBUG
#		s+=']' #DEBUG
#		print s #DEBUG
		#open the file
		filep = TFile(ttree_file_path)
		#get the tree from the file
		tree = filep.Get('tree')
		nEntries = tree.GetEntries()
		#make lists of observables that we'll need from the original ttree
		self.__ttree_cut_branches = make_dict_of_ttree_cut_branches()
		#for every event in the file
		for entry in range(nEntries) :
			percentdone = 100.*entry/nEntries
			if percentdone%1.<100./nEntries and percentdone%1.<((100.*(entry-1)/nEntries)%1.) :
				print '		%d%% done (%d out of %d)'%(percentdone,entry,nEntries)
			#set branch addresses for cuts
			for branch in self.__ttree_cut_branches.values() :
				tree.SetBranchAddress(branch.getTTreeName(),branch.getTTreeArray())
			tree.GetEntry(entry)
			#Check basic selection cuts
			if self.__ttree_cut_branches['fullselection'].getTTreeValue()!=1 :
				continue
			#get the event type and lepton flavor once per event
			etype = self.__ttree_cut_branches['eventType'].getTTreeValue()
			lepflav = self.__ttree_cut_branches['lepflavor'].getTTreeValue()
			#for each channel
			for channel in self.__channel_list :
				cleptype = channel.getLepType()
				#check that the event topology and lepton flavors/charges agree
				if int(channel.getTopology().split('t')[1])!=self.__ttree_cut_branches['eventTopology'].getTTreeValue() :
					continue
				if cleptype=='mu' and lepflav!=1 or cleptype=='el' and lepflav!=2 :
					continue
				if channel.getCharge()!=0 and self.__ttree_cut_branches['lep_Q'].getTTreeValue()!=channel.getCharge() and self.__ttree_cut_branches['addTwice'].getTTreeValue()!=1 :
					continue
				#for each process in the channel
				for process in channel.getProcessList() :
					if not process in ps_for_file :
						continue
					#for ttbar files, cut on the event type
					if ttree_file_path.find('_TT')!=-1 :
						pname = process.getName()
						if (pname.find('fqq')!=-1 and etype!=0) or (pname.find('fgg')!=-1 and etype!=1) or (pname.find('fbck')!=-1 and etype!=2 and etype!=3) :
							continue
					#get the contribution weight
					contrib_weight = get_contrib_weight(ttree_file_path,process)
					#get the branch dictionary
					bdict = process.getBranchDict()
					#explicitly fill the contribution weight and luminosity branches for this process
					bdict['contrib_weight'].getPTreeArray()[0] = contrib_weight
					if 'lumi_BtoF' in bdict.keys() :
						if process.isMCProcess() :
							if lepflav==1 :
								bdict['lumi_BtoF'].getPTreeArray()[0] = self.__LUMINOSITY_BTOF_MU
								bdict['lumi_BtoF_up'].getPTreeArray()[0] = (1.0+self.__LUMI_MIN_BIAS)*self.__LUMINOSITY_BTOF_MU
								bdict['lumi_BtoF_down'].getPTreeArray()[0] = (1.0-self.__LUMI_MIN_BIAS)*self.__LUMINOSITY_BTOF_MU
							elif lepflav==2 :
								bdict['lumi_BtoF'].getPTreeArray()[0] = self.__LUMINOSITY_BTOF_EL
								bdict['lumi_BtoF_up'].getPTreeArray()[0] = (1.0+self.__LUMI_MIN_BIAS)*self.__LUMINOSITY_BTOF_EL
								bdict['lumi_BtoF_down'].getPTreeArray()[0] = (1.0-self.__LUMI_MIN_BIAS)*self.__LUMINOSITY_BTOF_EL
						else :
							bdict['lumi_BtoF'].getPTreeArray()[0] = 1.0
							bdict['lumi_BtoF_up'].getPTreeArray()[0] = 1.0
							bdict['lumi_BtoF_down'].getPTreeArray()[0] = 1.0
					if 'lumi_GH' in bdict.keys() :
						if process.isMCProcess() :
							if lepflav==1 :
								bdict['lumi_GH'].getPTreeArray()[0] = self.__LUMINOSITY_GH_MU
								bdict['lumi_GH_up'].getPTreeArray()[0] = (1.0+self.__LUMI_MIN_BIAS)*self.__LUMINOSITY_GH_MU
								bdict['lumi_GH_down'].getPTreeArray()[0] = (1.0-self.__LUMI_MIN_BIAS)*self.__LUMINOSITY_GH_MU
							elif lepflav==2 :
								bdict['lumi_GH'].getPTreeArray()[0] = self.__LUMINOSITY_GH_EL
								bdict['lumi_GH_up'].getPTreeArray()[0] = (1.0+self.__LUMI_MIN_BIAS)*self.__LUMINOSITY_GH_EL
								bdict['lumi_GH_down'].getPTreeArray()[0] = (1.0-self.__LUMI_MIN_BIAS)*self.__LUMINOSITY_GH_EL
						else :
							bdict['lumi_GH'].getPTreeArray()[0] = 1.0
							bdict['lumi_GH_up'].getPTreeArray()[0] = 1.0
							bdict['lumi_GH_down'].getPTreeArray()[0] = 1.0
					#set the branch addresses and re-get the event
					for branch in process.getBranchDict().values() :
						bttreename = branch.getTTreeName()
						if bttreename!=None :
							tree.SetBranchAddress(bttreename,branch.getPTreeArray())
					tree.GetEntry(entry)
					#fill the tree for this process
					process.fillTree(ttree_file_path)

	def build_templates(self) :
		#Start of with the weightsums
		self.__build_weightsums__()
		#then build all of the templates
		for c in self.__channel_list :
			for p in c.getProcessList() :
				if p.isMCProcess() or p.isDataProcess() :
					p.buildTemplates(c.getCharge())

	def get_list_of_process_trees(self) :
		treelist = []
		for chan in self.__channel_list :
			for p in chan.getProcessList() :
				for t in p.getTreeDict().values() :
					treelist.append(t)
		return treelist

	def get_list_of_auxiliary_objects(self) :
		for c in self.__channel_list :
			for p in c.getProcessList() :
				for t in p.getTemplateList() :
					self.__aux_obj_list.append(t.getHisto3D())
					self.__aux_obj_list.append(t.getHistoX())
					self.__aux_obj_list.append(t.getHistoY())
					self.__aux_obj_list.append(t.getHistoZ())
		return self.__aux_obj_list

	def get_list_of_1D_histos(self) :
		hist_list = []
		for c in self.__channel_list :
			for p in c.getProcessList() :
				for t in p.getTemplateList() :
					hist_list.append(t.convertTo1D())
		return hist_list

#	#return true if the event passes the muon hadronic pretagging criteria
#	def __event_passes_mu_pretag__(self) :
#		branches = self.__ttree_cut_branches
#		flav_ind_preselection = branches['scaled_hadt_pt'].getTTreeValue()>300. and branches['scaled_hadt_M'].getTTreeValue()>50.
#		muon_preselection = flav_ind_preselection and branches['muTrig'].getTTreeValue()==1 and branches['muon1_pt'].getTTreeValue()>branches['ele1_pt'].getTTreeValue() 
#		muon_kinematics   = branches['muon1_pt'].getTTreeValue()>40. and abs(branches['muon1_eta'].getTTreeValue())<2.4
#		muon_ID = branches['muon1_ID'].getTTreeValue()==1
#		muon_2D = branches['muon1_relPt'].getTTreeValue()>25. or branches['muon1_dR'].getTTreeValue()>0.5
#		ltm = branches['scaled_lept_M'].getTTreeValue()
#		lep_top_mass = ltm>140. and ltm<900. #CHANGED
#		muon_full_leptonic = muon_preselection and muon_kinematics and muon_ID and muon_2D and lep_top_mass 
#		muon_hadronic_pretag = muon_full_leptonic and branches['hadt_tau21'].getTTreeValue()>0.1
#		return muon_hadronic_pretag
#
#	#return true if the event passes the electron hadronic pretagging criteria
#	def __event_passes_el_pretag__(self) :
#		branches = self.__ttree_cut_branches
#		flav_ind_preselection = branches['scaled_hadt_pt'].getTTreeValue()>300. and branches['scaled_hadt_M'].getTTreeValue()>50.
#		ele_preselection = flav_ind_preselection and branches['elTrig'].getTTreeValue()==1 and branches['ele1_pt'].getTTreeValue()>branches['muon1_pt'].getTTreeValue()
#		ele_kinematics = branches['ele1_pt'].getTTreeValue()>40. and abs(branches['ele1_eta'].getTTreeValue())<2.4
#		ele_ID = branches['ele1_ID'].getTTreeValue()==1
#		ele_2D = branches['ele1_relPt'].getTTreeValue()>25. or branches['ele1_dR'].getTTreeValue()>0.5
#		ltm = branches['scaled_lept_M'].getTTreeValue()
#		lep_top_mass = ltm>140. and ltm<900. #CHANGED
#		ele_full_leptonic = ele_preselection and ele_kinematics and ele_ID and ele_2D and lep_top_mass
#		ele_hadronic_pretag = ele_full_leptonic and branches['hadt_tau21'].getTTreeValue()>0.1
#		return ele_hadronic_pretag

	#Build the weightsums for function replacement in the MC processes
	def __build_weightsums__(self) :
		#Make a dictionary of all the weightsums for each possible template type
		all_weightsums = {}
		for p in self.__channel_list[0].getProcessList() :
			if not p.isMCProcess() :
				continue
			for t in p.getTemplateList() :
				ttype = t.getType()
				if not ttype in all_weightsums.keys() :
					all_weightsums[ttype]=copy.deepcopy(t.getWeightsumDict())
		#Get weightsums for individual channels and processes
		for channel in self.__channel_list :
			for p in channel.getProcessList() :
				if p.isMCProcess() :
					print '	Getting weightsums for MC process '+p.getName()
					p.buildWeightsums(channel.getCharge())
		#sum all of the channels/processes together by template type
		for ttype in all_weightsums :
			for c in self.__channel_list :
				for p in c.getProcessList() :
					if not p.isMCProcess() :
						continue
					for t in p.getTemplateList() :
						if t.getType()==ttype or (t.getType()=='nominal' and ttype not in p.getListOfTTypes()):
							for ws in all_weightsums[ttype] :
								for ti in all_weightsums[ttype][ws] :
									if ti=='wname' :
										continue
									all_weightsums[ttype][ws][ti]+=t.getWeightsumDict()[ws][ti]
		#Copy over the weightsum dictionaries based on template types
		for c in self.__channel_list :
			for p in c.getProcessList() :
				for t in p.getTemplateList() :
					ttype = t.getType()
					if ttype=='fit__up' or ttype=='fit__down' :
						ttype='nominal'
					t.setWeightsumDict(copy.deepcopy(all_weightsums[ttype]))
	#	print '	FINAL WEIGHTSUMS: '
	#	for c in self.__channel_list :
	#		for p in c.getProcessList() :
	#			if not p.isMCProcess() :
	#				continue
	#			for t in p.getTemplateList() :
	#				print '		Template '+t.getName()+': '
	#				wsd = t.getWeightsumDict()
	#				for weightsum_name in wsd :
	#					s = '			'+weightsum_name+': {'
	#					for ttree_identifier in wsd[weightsum_name] :
	#						if ttree_identifier=='wname' :
	#							continue
	#						s+='%s:%.2f,'%(ttree_identifier,wsd[weightsum_name][ttree_identifier])
	#					s+='}'
	#					print s

	def __del__(self) :
		pass

	def __str__(self) :
		s = 'Template_Group object:\n'
		s+= '	channel_list = ['
		for i in range(len(self.__channel_list)) :
			chan = self.__channel_list[i]
			if i==0 :
				s+=chan.__str__()
			else :
				s+='		'+chan.__str__()
			if i==len(self.__channel_list)-1 :
				s+=']\n'
			else :
				s+=',\n'
		return s

#make the list of channels depending on the types of leptons and whether or not to sum their charges
def make_channel_list(fit_parameter_tuple,include_mu,include_el,sum_charges,topology_list,include_JEC,include_sss) :
	channel_list = []
	for topology in topology_list :
		if include_mu :
			if sum_charges :
				channel_list.append(Channel(topology+'_mu',fit_parameter_tuple,include_JEC,include_sss))
			else :
				channel_list.append(Channel(topology+'_muplus',fit_parameter_tuple,include_JEC,include_sss))
				channel_list.append(Channel(topology+'_muminus',fit_parameter_tuple,include_JEC,include_sss))
		if include_el :
			if sum_charges :
				channel_list.append(Channel(topology+'_el',fit_parameter_tuple,include_JEC,include_sss))
			else :
				channel_list.append(Channel(topology+'_elplus',fit_parameter_tuple,include_JEC,include_sss))
				channel_list.append(Channel(topology+'_elminus',fit_parameter_tuple,include_JEC,include_sss))
	return channel_list

#make the list of branches that are only needed for checking selection cuts and nothing else
def make_dict_of_ttree_cut_branches() :
	branches = {}
	#branches needed only for cuts (not anything else)
	branches['fullselection']=Branch('fullselection',None,'I',2)
	branches['eventTopology']=Branch('eventTopology',None,'i',0)
	branches['eventType']=Branch('eventType',None,'i',5)
	branches['addTwice']=Branch('addTwice',None,'I',0)
	branches['lep_Q']=Branch('lep_Q',None,'i',0)
	branches['lepflavor']=Branch('lepflavor',None,'I',0)
	return branches

#return the weight that should be applied to this event when adding to the given process
def get_contrib_weight(filepath,process) :
	pname = process.getName()
	#if it's a data process return 1 only if it's a data file
	if process.isDataProcess() :
		if filepath.find('SingleMu')!=-1 and pname[3:].startswith('mu')!=-1 :
			return 1.0
		elif filepath.find('SingleEl')!=-1 and pname[3:].startswith('el')!=-1 :
			return 1.0
		else :
			return 0.
	#otherwise if it's a MC process
	elif process.isMCProcess() :
		#check the ttbar processes
		if pname.find('fqq')!=-1 or pname.find('fgg')!=-1 or pname.find('fbck')!=-1 :
			if filepath.find('_TT')!=-1 :
				return 1.0
			return 0.
		#check the background processes
		elif pname.find('fbck')!=-1 :
			bck_stems = ['ST_s-c','ST_t-c_top','ST_tW-c_top','ST_t-c_antitop','ST_tW-c_antitop']
			for bckstem in bck_stems :
				if filepath.find(bckstem)!=-1 :
					return 1.0
			return 0.
		else :
			print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
			print '!!!!!!   WARNING, MC PROCESS TYPE NOT RECOGNIZED FROM NAME '+pname+'   !!!!!!'
			print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
			return 0.
	#otherwise I'm super confused about what type of process it is in the first place haha
	else :
		print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
		print '!!!!!!   WARNING, PROCESS TYPE NOT RECOGNIZED FROM NAME '+pname+'   !!!!!!'
		print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
		return 0.
