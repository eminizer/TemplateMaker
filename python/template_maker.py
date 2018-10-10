import copy
from math import *
from array import array
from ROOT import TFile, TTree, Double, TH1D
from channel import Channel
from branch import Branch
import multiprocessing
import os

#Luminosities
LUMINOSITY_BTOF_MU = 19690.184
LUMINOSITY_BTOF_EL = 19171.010
LUMINOSITY_GH_MU   = 16226.452
LUMINOSITY_GH_EL   = 16214.862
LUMI_MIN_BIAS = 0.025

#Template_Group class
class Template_Group(object) :

	def __init__(self,fit_parameter_tuple,outname,sum_charges,include_mu,include_el,topology_list,include_JEC,include_sss,n_procs) :
		#switches for which templates to auto-create
		self.__include_JEC = include_JEC
		self.__include_mu = include_mu
		self.__include_el = include_el
		self.__include_sss = include_sss
		self.__n_procs = n_procs
		#make the list of channels, which will have their processes and templates automatically added
		self.__channel_list = make_channel_list(fit_parameter_tuple,include_mu,include_el,sum_charges,topology_list,include_JEC,include_sss)
		self.__outname = outname 
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
		#make lists of observables that we'll need from the original ttree
		self.__ttree_cut_branches = make_dict_of_ttree_cut_branches()
		#copy every branch we need from the tree (for cuts and for adding to processes) to a dictionary
		branches_to_copy = {}
		for tcb in self.__ttree_cut_branches.values() :
			branches_to_copy[tcb.getTTreeName()] = copy.deepcopy(tcb)
		for chan in self.__channel_list :
			for proc in chan.getProcessList() :
				bdict = proc.getBranchDict()
				for b in bdict.values() :
					if b.getTTreeName()!=None and b.getTTreeName() not in branches_to_copy.keys() :
						branches_to_copy[b.getTTreeName()] = copy.deepcopy(b)
		#spawn processes for adding events to appropriate processes
		procs = []
		#first segment the tree
		print '		Segmenting trees for %d parallel processes...'%(self.__n_procs)
		for i in range(self.__n_procs) :
			p = multiprocessing.Process(target=self.copyTreeSegment, args=(ttree_file_path,i,copy.deepcopy(branches_to_copy)))
			procs.append(p)
			p.start()
		for p in procs :
			p.join()
		print '		Done.'
		#now actually add the events from the segmented trees to the total process
		print '		Adding events from trees to process trees...'
		procs = []
		for i in range(self.__n_procs) :
			p = multiprocessing.Process(target=self.filterCopyTree, args=(i==self.__n_procs-1,i,ttree_file_path,ps_for_file))
			#p = multiprocessing.Process(target=self.filterCopyTree, args=(True,i,ttree_file_path,ps_for_file))
			procs.append(p)
			p.start()
		#make sure all processes completed 
		for p in procs :
			p.join()
		#remove the garbage files that held the segmented trees
		print '		Removing segmented tree files'
		os.system('rm -rf segment_*.root')
		print '		Done.'
		#hadd together the separate process tree files
		print '		Aggregating process tree files'
		filetype = ttree_file_path.split('/')[-1].rstrip('skim_all.root')
		cmd = 'hadd -f '+filetype+'_'+self.__outname+'_process_trees_all.root '
		for i in range(self.__n_procs) :
			cmd+=filetype+'_'+self.__outname+'_process_trees_'+str(i+1)+'.root '
		os.system(cmd)
		os.system('rm -rf '+filetype+'_'+self.__outname+'_process_trees_?.root '+filetype+'_'+self.__outname+'_process_trees_??.root')
		print '		Done'

	def copyTreeSegment(self,ttree_file_path,i,bdict) :
		#open the file
		filep = TFile(ttree_file_path)
		#get the tree from the file
		tree = filep.Get('tree')
		nEntries = tree.GetEntries()
		#set branches to copy
		tree.SetBranchStatus('*',0)
		for b in bdict.values() :
			bname = b.getTTreeName()
			tree.SetBranchStatus(bname,1)
			tree.SetBranchAddress(bname,b.getTTreeArray())
		#figure out how many events to copy starting from where
		thisstart = i*(nEntries/self.__n_procs)
		thisend = thisstart+nEntries/self.__n_procs
		if i==self.__n_procs-1 :
			thisend = nEntries
		#open a garbage file to hold the new tree
		garbageFile = TFile('segment_'+str(i+1)+'.root','recreate')
		#print '		Begin copying tree segment %d of %d (entries %d to %d, %d of %d total entries)'%(i+1,self.__n_procs,thisstart,thisend,thisend-thisstart,nEntries) #DEBUG
		#thistree = tree.CopyTree('','',thisnentries,thisstart)
		thistree = tree.CloneTree(0)
		for ientry in range(thisstart,thisend) :
			tree.GetEntry(ientry)
			thistree.Fill()
		thistree.Write()
		garbageFile.Close()
		filep.Close()
		#print '		Done copying tree segment %d of %d'%(i+1,self.__n_procs) #DEBUG

	def filterCopyTree(self,printbool,iproc,ttree_file_path,ps_for_file) :
		#get the segmented tree
		filep = TFile('segment_'+str(iproc+1)+'.root')
		tree = filep.Get('tree')
		nentries = tree.GetEntries()
		#print 'number of entries for tree in %s = %d'%(multiprocessing.current_process().name,nentries) #DEBUG
		#print 'process: %s, tree = %s, tree_cut_branches = %s, chanlist = %s'%(multiprocessing.current_process().name,tree,tree_cut_branches,chanlist) #DEBUG
		tree.SetMakeClass(1)
		#open up the file that the process trees will go in
		ptreesfilep = TFile(ttree_file_path.split('/')[-1].rstrip('skim_all.root')+'_'+self.__outname+'_process_trees_'+str(iproc+1)+'.root','recreate')
		#add trees and keep track of them in a dictionary of dictionaries
		treedicts = {}
		for chan in self.__channel_list :
			for proc in chan.getProcessList() :
				procname = proc.getName()
				procbranches = copy.deepcopy(proc.getBranchDict())
				treedicts[procname] = {'branches':procbranches}
				for key, value in proc.getTreeDict().iteritems() :
					newtree = TTree(value.GetName(),'recreate')
					for branch in procbranches.values() :
						newtree.Branch(branch.getPTreeName(),branch.getPTreeArray(),branch.getPTreeLeafList())
					treedicts[procname][key] = newtree
		#get a copy of the branches for ttree cuts
		tree_cut_branches = copy.deepcopy(self.__ttree_cut_branches)
		#set branch addresses for cuts
		for branch in tree_cut_branches.values() :
			tree.SetBranchAddress(branch.getTTreeName(),branch.getTTreeArray())
		#just real quick see if we need to check the event type for every one of the events
		checkeventtype = ttree_file_path.find('powheg_TT')!=-1 or ttree_file_path.find('mcatnlo_TT')!=-1 
		#do selections and copy relevant events from reconstructor ttrees to process ttrees
		for entry in range(nentries) :
			if printbool :
				percentdone = 100.*entry/nentries
				if percentdone%1.<100./nentries and percentdone%10.<((100.*(entry-1)/nentries)%10.) :
					print '		%d%% done'%(percentdone)
			check = tree.GetEntry(entry)
			#Check basic selection cuts
			basiccuts = {}
			basiccuts['SR'] = getBasicCutDict(tree_cut_branches,'SR')
			basiccuts['WJets_CR'] = getBasicCutDict(tree_cut_branches,'WJets_CR')
			if basiccuts['SR'].values().count(True)==0 and basiccuts['WJets_CR'].values().count(True)==0 :
				continue
			#get the event type and lepton flavor once per event
			etype = tree_cut_branches['eventType'].getTTreeValue()
			lepflav = tree_cut_branches['lepflavor'].getTTreeValue()
			#look through each channel
			for channel in self.__channel_list :
				#first check that this event belongs in this channel's region of phase space
				cregion = channel.getRegion()
				if basiccuts[cregion].values().count(True)==0 :
					continue
				#check that the event topology and lepton flavors/charges agree (these are the same for every process in this channel)
				cleptype = channel.getLepType()
				if int(channel.getTopology().split('t')[1])!=tree_cut_branches['eventTopology'].getTTreeValue() :
					continue
				if (cleptype=='mu' and lepflav!=1) or (cleptype=='el' and lepflav!=2) :
					continue
				if channel.getCharge()!=0 and tree_cut_branches['lep_Q'].getTTreeValue()!=channel.getCharge() and tree_cut_branches['addTwice'].getTTreeValue()!=1 :
					continue
				#now look through each process in the channel
				for process in channel.getProcessList() :
					#first off make sure this file contributes to processes of this type
					if not process in ps_for_file :
						continue
					#get the name of this process
					pname = process.getName()
					#for ttbar files, cut on the event type
					if checkeventtype :
						if ((pname.find('fqq')!=-1 or pname.find('fqp')!=-1 or pname.find('fqm')!=-1) and etype!=0) or (pname.find('fg')!=-1 and etype!=1) or (pname.find('fbck')!=-1 and etype!=2 and etype!=3) :
							continue
					#if we made it to this point, this event should be copied into the process's tree(s)!!
					#get the contribution weight
					contrib_weight = get_contrib_weight(ttree_file_path,process)
					#get the branch dictionary for this process
					bdict = treedicts[pname]['branches']
					#set the branch addresses and re-get the event
					for branch in bdict.values() :
						bttreename = branch.getTTreeName()
						if bttreename!=None :
							tree.SetBranchAddress(bttreename,branch.getPTreeArray())
					tree.GetEntry(entry)
					#explicitly fill the contribution weight and luminosity branches for this process
					bdict['contrib_weight'].getPTreeArray()[0] = contrib_weight
					if 'lumi_BtoF' in bdict.keys() :
						if process.isMCProcess() or process.isFitProcess():
							if lepflav==1 :
								bdict['lumi_BtoF'].getPTreeArray()[0] = LUMINOSITY_BTOF_MU
								bdict['lumi_BtoF_up'].getPTreeArray()[0] = (1.0+LUMI_MIN_BIAS)*LUMINOSITY_BTOF_MU
								bdict['lumi_BtoF_down'].getPTreeArray()[0] = (1.0-LUMI_MIN_BIAS)*LUMINOSITY_BTOF_MU
							elif lepflav==2 :
								bdict['lumi_BtoF'].getPTreeArray()[0] = LUMINOSITY_BTOF_EL
								bdict['lumi_BtoF_up'].getPTreeArray()[0] = (1.0+LUMI_MIN_BIAS)*LUMINOSITY_BTOF_EL
								bdict['lumi_BtoF_down'].getPTreeArray()[0] = (1.0-LUMI_MIN_BIAS)*LUMINOSITY_BTOF_EL
						else :
							bdict['lumi_BtoF'].getPTreeArray()[0] = 1.0
							bdict['lumi_BtoF_up'].getPTreeArray()[0] = 1.0
							bdict['lumi_BtoF_down'].getPTreeArray()[0] = 1.0
					if 'lumi_GH' in bdict.keys() :
						if process.isMCProcess() or process.isFitProcess() :
							if lepflav==1 :
								bdict['lumi_GH'].getPTreeArray()[0] = LUMINOSITY_GH_MU
								bdict['lumi_GH_up'].getPTreeArray()[0] = (1.0+LUMI_MIN_BIAS)*LUMINOSITY_GH_MU
								bdict['lumi_GH_down'].getPTreeArray()[0] = (1.0-LUMI_MIN_BIAS)*LUMINOSITY_GH_MU
							elif lepflav==2 :
								bdict['lumi_GH'].getPTreeArray()[0] = LUMINOSITY_GH_EL
								bdict['lumi_GH_up'].getPTreeArray()[0] = (1.0+LUMI_MIN_BIAS)*LUMINOSITY_GH_EL
								bdict['lumi_GH_down'].getPTreeArray()[0] = (1.0-LUMI_MIN_BIAS)*LUMINOSITY_GH_EL
						else :
							bdict['lumi_GH'].getPTreeArray()[0] = 1.0
							bdict['lumi_GH_up'].getPTreeArray()[0] = 1.0
							bdict['lumi_GH_down'].getPTreeArray()[0] = 1.0
					#fill the tree for this process
					fillTree(ttree_file_path,process.getJECModList(),treedicts[pname],process.isMCProcess(),process.isFitProcess(),basiccuts[cregion])
					#reset the branch addresses
					for branch in tree_cut_branches.values() :
						tree.SetBranchAddress(branch.getTTreeName(),branch.getTTreeArray())
		#write all the trees
		for p in treedicts.keys() :
			for key, value in treedicts[p].iteritems() :
				if key!='branches' :
					treedicts[p][key].Write()
		#close the files
		ptreesfilep.Close()
		filep.Close()

	def build_QCD_templates(self,ptree_filename) :
		self.__ptree_filename = ptree_filename
		#Start of with the weightsums to put in the fit functions everywhere
		self.__build_weightsums__()
		#Next build each template's conversion factor
		self.__calculate_conversion_factors__()
		#then build all of the templates
		for c in self.__channel_list :
			for p in c.getProcessList() :
				if p.isQCDProcess() :
					p.buildTemplates(c.getCharge(),self.__ptree_filename,self.__n_procs)
		#finally set the NQCD numbers now that the templates are actually built
		self.__set_NQCD_values__()

	def build_templates(self) :
		#build all of the templates
		for c in self.__channel_list :
			for p in c.getProcessList() :
				if p.isMCProcess() or p.isDataProcess() or p.isFitProcess() :
					p.buildTemplates(c.getCharge(),self.__ptree_filename,self.__n_procs)

	def get_list_of_process_trees(self) :
		treelist = []
		for chan in self.__channel_list :
			for p in chan.getProcessList() :
				for t in p.getTreeDict().values() :
					treelist.append(t)
		return treelist

	def get_list_of_auxiliary_objects(self) :
		#add histogram projections
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
		#convert templates to 1D
		for c in self.__channel_list :
			#do observed data templates first
			for p in c.getProcessList() :
				if p.isDataProcess() :
					for t in p.getTemplateList() :
						temp1D,binlist=t.convertTo1D()
						hist_list.append(temp1D)
			#next do all the nominal templates, making a dictionary of the bins to zero for each process along the way
			bins_to_zero = {}
			nom_1D_histos = {}
			for p in c.getProcessList() :
				if p.isDataProcess() :
					continue
				for t in p.getTemplateList() :
					if t.getModifier()==None :
						temp1D,binlist=t.convertTo1D()
						bins_to_zero[p.getName()]=binlist
						nom_1D_histos[p.getName()]=temp1D
						hist_list.append(temp1D)
			#finally do all the systematics wiggle templates
			for p in c.getProcessList() :
				if p.isDataProcess() :
					continue
				for t in p.getTemplateList() :
					if t.getModifier()==None :
						continue
					temp1D,binlist=t.convertTo1D(bins_to_zero[p.getName()],nom_1D_histos[p.getName()])
					hist_list.append(temp1D)
		return hist_list

	#Build the weightsums for function replacement in the MC processes
	def __build_weightsums__(self) :
		#Make a dictionary of all the weightsums for each possible template type
		all_weightsums = {}
		#first key is topology_region
		for c in self.__channel_list :
			cregion = c.getRegion(); ctopology = c.getTopology()
			cidentifier = ctopology+'_'+cregion
			if not cidentifier in all_weightsums.keys() :
				all_weightsums[cidentifier] = {}
			#second key is the type of template (parameter wiggles)
			for p in c.getProcessList() :
				if not p.isMCProcess() :
					continue
				for t in p.getTemplateList() :
					ttype = t.getType()
					if not ttype in all_weightsums[cidentifier].keys() :
						all_weightsums[cidentifier][ttype]=copy.deepcopy(t.getWeightsumDict())
		#Get weightsums for individual channels and processes
		for channel in self.__channel_list :
			for p in channel.getProcessList() :
				if p.isMCProcess() :
					print '	Getting weightsums for MC process '+p.getName()
					p.buildWeightsums(channel.getCharge(),self.__ptree_filename,self.__n_procs)
		#sum all of the channels/processes together by template type
		for cidentifier in all_weightsums :
			for ttype in all_weightsums[cidentifier] :
				for c in self.__channel_list :
					if not c.getTopology()+'_'+c.getRegion()==cidentifier :
						continue
					for p in c.getProcessList() :
						if not p.isMCProcess() :
							continue
						for t in p.getTemplateList() :
							if t.getType()==ttype or (t.getType()=='nominal' and ttype not in p.getListOfTTypes()):
								for ws in all_weightsums[cidentifier][ttype] :
									for ti in all_weightsums[cidentifier][ttype][ws] :
										if ti=='wname' :
											continue
										all_weightsums[cidentifier][ttype][ws][ti]+=t.getWeightsumDict()[ws][ti]
		#Copy over the weightsum dictionaries based on template types
		for c in self.__channel_list :
			cregion = c.getRegion(); ctopology = c.getTopology()
			cidentifier = ctopology+'_'+cregion
			for p in c.getProcessList() :
				for t in p.getTemplateList() :
					t.setWeightsumDict(copy.deepcopy(all_weightsums[cidentifier][t.getType()]))
		#save nominal values for everything except QCD in the 'sig' tree to the auxiliary file
		for cidentifier in all_weightsums.keys() :
			for ws in all_weightsums[cidentifier]['nominal'].keys() :
				if ws!='NQCD' :
					newweightsumhist = TH1D(cidentifier+'_'+ws,cidentifier+'_'+ws,1,0.,1.)
					newweightsumhist.Fill(0.5,all_weightsums[cidentifier]['nominal'][ws]['sig'])
					self.__aux_obj_list.append(newweightsumhist)
		print '	FINAL WEIGHTSUMS: '
		for c in self.__channel_list :
			cregion = c.getRegion(); ctopology = c.getTopology()
			cidentifier = ctopology+'_'+cregion
			print '	Topology_Region '+cidentifier+': '
			for p in c.getProcessList() :
				if not p.isMCProcess() :
					continue
				for t in p.getTemplateList() :
					print '		Template '+t.getName()+': '
					wsd = t.getWeightsumDict()
					for weightsum_name in wsd :
						s = '			'+weightsum_name+': {'
						for ttree_identifier in wsd[weightsum_name] :
							if ttree_identifier=='wname' :
								continue
							s+='%s:%.2f,'%(ttree_identifier,wsd[weightsum_name][ttree_identifier])
						s+='}'
						print s

	#calculate the conversion factors to apply to each QCD template's shape in the sideband
	def __calculate_conversion_factors__(self) :
		print '	Calculating conversion factors'
		#make a dictionary of event numbers
		event_numbers = {}
		for c in self.__channel_list :
			event_numbers[c.getName()] = {}
			charge = c.getCharge()
			for p in c.getProcessList() :
				if p.isQCDProcess() :
					event_numbers[c.getName()] = p.getEventNumbers(charge,self.__ptree_filename,self.__n_procs)
		#Make the conversion factors dictionary per channel/template type
		conv_fac_dict = {}
		for c in self.__channel_list :
			cname = c.getName()
			if cname not in conv_fac_dict.keys() :
				conv_fac_dict[cname] = {}
			for p in c.getProcessList() :
				if not p.isQCDProcess() :
					continue
				for t in p.getTemplateList() :
					ttype = t.getType()
					if not ttype in conv_fac_dict[cname].keys() :
						print '		Finding conversion factor for template '+t.getName()
						#make the conversion function
						n_qcd_a = event_numbers[cname][ttype]['qcd_a']
						n_qcd_b = event_numbers[cname][ttype]['qcd_b']
						conv_fac_dict[cname][ttype] = n_qcd_a/n_qcd_b if n_qcd_b!=0. else 0.
						print '			n_qcd_a=%.2f, n_qcd_b=%.2f, factor=%.5f'%(n_qcd_a,n_qcd_b,conv_fac_dict[cname][ttype])
				for t in p.getTemplateList() :
					ttype = t.getType()
					t.setConversionFactor(conv_fac_dict[cname][ttype])

	def __set_NQCD_values__(self) :
		print '	Getting and resetting total NQCD values'
		#get the individual values
		nqcd_values = {}
		for c in self.__channel_list :
			nqcd_values[c.getName()] = {}
			for p in c.getProcessList() :
				if not p.isQCDProcess() :
					continue
				for t in p.getTemplateList() :
					ttype = t.getType()
					nqcd_values[c.getName()][ttype] = t.fixNQCDValues()
		#sum them over channels per region and NTMJ process type
		nqcd_values2 = {}
		for c in self.__channel_list :
			cregion = c.getRegion(); ctopology = c.getTopology()
			cidentifier = ctopology+'_'+cregion
			if not cidentifier in nqcd_values2.keys() :
				nqcd_values2[cidentifier] = {}
			for p in c.getProcessList() :
				if not p.isQCDProcess() :
					continue
				for t in p.getTemplateList() :
					ttype = t.getType()
					if not ttype in nqcd_values2[cidentifier].keys() :
						nqcd_values2[cidentifier][ttype]=0.
					nqcd_values2[cidentifier][ttype]+=nqcd_values[c.getName()][ttype]
		#print 'nqcd_values2=%s'%(nqcd_values2) #DEBUG
		#set the values in the templates
		for c in self.__channel_list :
			if c.getTopology()=='t1' or (c.getTopology()=='t2' and c.getRegion()!='WJets_CR') :
				continue
			for p in c.getProcessList() :
				for t in p.getTemplateList() :
					print '		Resetting NQCD value for template %s (channel %s, process %s, template type %s'%(t.getName(),c.getName(),p.getName(),t.getType())
					t.setNQCDValue(nqcd_values2[c.getTopology()+'_'+c.getRegion()][t.getType()])
		#add the nominal NQCD values in each region to the auxiliary object list
		for cidentifier in nqcd_values2.keys() :
			if len(nqcd_values2[cidentifier].keys())<1 :
				continue
			#print 'Getting new histogram for cidentifier %s with keys:'%(cidentifier) #DEBUG
			#for k in nqcd_values2[cidentifier].keys() : #DEBUG
			#	print '	'+str(k) #DEBUG
			newweightsumhist = TH1D(cidentifier+'_NQCD',cidentifier+'_NQCD',1,0.,1.)
			newweightsumhist.Fill(0.5,nqcd_values2[cidentifier]['nominal'])
			self.__aux_obj_list.append(newweightsumhist)

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
		region_names = ['SR','WJets_CR'] if topology!='t3' else ['SR']
		for region in region_names :
			if include_mu :
				if sum_charges :
					channel_list.append(Channel(topology+'_mu_'+region,fit_parameter_tuple,include_JEC,include_sss))
				else :
					channel_list.append(Channel(topology+'_muplus_'+region,fit_parameter_tuple,include_JEC,include_sss))
					channel_list.append(Channel(topology+'_muminus_'+region,fit_parameter_tuple,include_JEC,include_sss))
			if include_el :
				if sum_charges :
					channel_list.append(Channel(topology+'_el_'+region,fit_parameter_tuple,include_JEC,include_sss))
				else :
					channel_list.append(Channel(topology+'_elplus_'+region,fit_parameter_tuple,include_JEC,include_sss))
					channel_list.append(Channel(topology+'_elminus_'+region,fit_parameter_tuple,include_JEC,include_sss))
	return channel_list

#make the list of branches that are only needed for checking selection cuts and nothing else
def make_dict_of_ttree_cut_branches() :
	branches = {}
	#branches needed only for cuts (not anything else)
	branches['fullselection']=Branch('fullselection',None,'I',2)
	branches['wjets_cr_selection']=Branch('wjets_cr_selection',None,'I',2)
	branches['qcd_A_SR_selection']=Branch('qcd_A_SR_selection',None,'I',2)
	branches['qcd_B_SR_selection']=Branch('qcd_B_SR_selection',None,'I',2)
	branches['qcd_C_SR_selection']=Branch('qcd_C_SR_selection',None,'I',2)
	branches['qcd_A_CR_selection']=Branch('qcd_A_CR_selection',None,'I',2)
	branches['qcd_B_CR_selection']=Branch('qcd_B_CR_selection',None,'I',2)
	branches['qcd_C_CR_selection']=Branch('qcd_C_CR_selection',None,'I',2)
	branches['eventTopology']=Branch('eventTopology',None,'i',0)
	branches['eventType']=Branch('eventType',None,'i',5)
	branches['addTwice']=Branch('addTwice',None,'I',0)
	branches['lep_Q']=Branch('lep_Q',None,'i',0)
	branches['lepflavor']=Branch('lepflavor',None,'I',0)
	return branches

#return the weight that should be applied to this event when adding to the given process
def get_contrib_weight(filepath,process) :
	pname = process.getName()
	#if it's a data or QCD process return 1 only if it's a data file
	if process.isDataProcess() or process.isQCDProcess() :
		if filepath.find('SingleMu')!=-1 and pname[3:].startswith('mu') :
			return 1.0
		elif filepath.find('SingleEl')!=-1 and pname[3:].startswith('el') :
			return 1.0
		else :
			return 0.
	#otherwise if it's a MC or Fit process
	elif process.isMCProcess() or process.isFitProcess() :
		#make sure it's not a data file
		if filepath.find('SingleMu')!=-1 or filepath.find('SingleEl')!=-1 :
			return 0.
		#check the ttbar processes
		if pname.split('__')[1].startswith('fq') or pname.split('__')[1].startswith('fg') :
			if filepath.find('powheg_TT')!=-1 or filepath.find('mcatnlo_TT')!=-1 :
				return 1.0
			return 0.
		#check the WJets processes
		if pname.find('fwjets')!=-1 :
			if filepath.find('WJets_')!=-1 :
				return 1.0
			return 0.
		#check the background processes (non-semilep ttbar, single top, and DYJets)
		elif pname.find('fbck')!=-1 :
			if filepath.find('powheg_TT')!=-1 or filepath.find('mcatnlo_TT')!=-1 :
				return 1.0
			else :
				bck_stems = ['ST_s-c','ST_t-c_top','ST_tW-c_top','ST_t-c_antitop','ST_tW-c_antitop','DYJets_M-50_HT-']
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

def getBasicCutDict(bs,region) :
	returndict = {}
	if region=='SR' :
		returndict['fullselection'] = bs['fullselection'].getTTreeValue()==1
		returndict['qcd_A'] = bs['qcd_A_SR_selection'].getTTreeValue()==1
		returndict['qcd_B'] = bs['qcd_B_SR_selection'].getTTreeValue()==1
		returndict['qcd_C'] = bs['qcd_C_SR_selection'].getTTreeValue()==1
	elif region=='WJets_CR' :
		returndict['fullselection'] = bs['wjets_cr_selection'].getTTreeValue()==1
		returndict['qcd_A'] = bs['qcd_A_CR_selection'].getTTreeValue()==1
		returndict['qcd_B'] = bs['qcd_B_CR_selection'].getTTreeValue()==1
		returndict['qcd_C'] = bs['qcd_C_CR_selection'].getTTreeValue()==1
	return returndict

def fillTree(ttree_file_path,jecmodlist,ptreedict,isMCProcess,isFitProcess,regioncutdict) :
	#print 'regioncutdict = %s'%(regioncutdict) #DEBUG
	treestring = 'sig'
	if regioncutdict['qcd_A'] :
		treestring='qcd_a'
	elif regioncutdict['qcd_B'] :
		treestring='qcd_b'
	elif regioncutdict['qcd_C'] :
		treestring='qcd_c'
	#If it's a MC or Fit distribution, add the JEC part of the name to the tree identifier if necessary
	if (isMCProcess or isFitProcess) and (ttree_file_path.find('_up')!=-1 or ttree_file_path.find('_down')!=-1) :
		JEC_ID = ''
		for jecmod in jecmodlist :
			if ttree_file_path.find(jecmod.getName()+'_up') != -1 :
				JEC_ID = jecmod.getName()+'Up'
				break
			elif ttree_file_path.find(jecmod.getName()+'_down') != -1 :
				JEC_ID = jecmod.getName()+'Down'
				break
		treestring+='_'+JEC_ID
	#if treestring!='sig' : #DEBUG
	#	print 'treestring = %s'%(treestring) #DEBUG
	ptreedict[treestring].Fill()
