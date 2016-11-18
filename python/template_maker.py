import copy
from math import *
from array import array
from ROOT import TFile, TTree, TF1, TCanvas, TLegend, TGraphErrors, TVirtualFitter, kBlue
from lepton_type import Channel
from branch import Branch

#Template_Group class
class Template_Group(object) :

	#Luminosity
	__LUMINOSITY = 12917.
	#Luminosity min bias
	__LUMI_MIN_BIAS = 0.062

	def __init__(self,fit_parameter_tuple,sum_charges,include_mu,include_el,include_PDF,include_JEC,include_sss) :
		self.__include_JEC = include_JEC
		self.__include_mu = include_mu
		self.__include_el = include_el
		self.__include_sss = include_sss
		#make the list of channels, which will have their processes and templates automatically added
		self.__channel_list = make_channel_list(fit_parameter_tuple,include_mu,include_el,sum_charges,include_PDF,include_JEC,include_sss)
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
			if percentdone%10.<100./nEntries and percentdone%10.<((100.*(entry-1)/nEntries)%10.) :
				print '		%d%% done (%d out of %d)'%(percentdone,entry,nEntries)
			#set branch addresses for cuts
			for branch in self.__ttree_cut_branches.values() :
				tree.SetBranchAddress(branch.getTTreeName(),branch.getTTreeArray())
			tree.GetEntry(entry)
			#Check basic selection cuts
			if self.__ttree_cut_branches['fullselection'].getTTreeValue()!=1 or self.__ttree_cut_branches['hadt_SDM'].getTTreeValue()<25. or self.__ttree_cut_branches['hadt_SDM'].getTTreeValue()>750. or self.__ttree_cut_branches['M'].getTTreeValue()<750. :
				continue
			#for each channel
			for channel in self.__channel_list :
				#for each process in the channel
				for process in channel.getProcessList() :
					if not process in ps_for_file :
						continue
					#get the contribution weight
					contrib_weight = get_contrib_weight(ttree_file_path,process)
					#set the branch addresses and re-get the event
					for branch in process.getBranchDict().values() :
						if branch.getTTreeName()!=None :
							tree.SetBranchAddress(branch.getTTreeName(),branch.getPTreeArray())
					tree.GetEntry(entry)
					bdict = process.getBranchDict()
					#check the charge separation
					if channel.getCharge()!=0 and bdict['Q_l'].getPTreeValue()!=channel.getCharge() and bdict['addTwice'].getPTreeValue()!=1 :
						continue
					#explicitly fill the contribution weight and luminosity branches for this process
					bdict['contrib_weight'].getPTreeArray()[0] = contrib_weight
					if 'luminosity' in bdict.keys() :
						if process.isMCProcess() :
							bdict['luminosity'].getPTreeArray()[0] = self.__LUMINOSITY
							bdict['luminosity_up'].getPTreeArray()[0] = (1.0+self.__LUMI_MIN_BIAS)*self.__LUMINOSITY
							bdict['luminosity_down'].getPTreeArray()[0] = (1.0-self.__LUMI_MIN_BIAS)*self.__LUMINOSITY
						else :
							bdict['luminosity'].getPTreeArray()[0] = 1.0
							bdict['luminosity_up'].getPTreeArray()[0] = 1.0
							bdict['luminosity_down'].getPTreeArray()[0] = 1.0
					#fill the tree for this process
					process.fillTree(ttree_file_path)

	def build_NTMJ_templates(self) :
		#build the weightsums for using in applying the function reweights
		self.__build_weightsums__()
		#build the conversion functions to apply 
		self.__build_conversion_functions__()
		#build the templates by applying the conversion functions and so on
		self.__build_NTMJ_templates__()
		#set the NNTMJ numbers
		self.__set_NNTMJ_values__()

	def build_templates(self) :
		for c in self.__channel_list :
			for p in c.getProcessList() :
				if not (p.isMCProcess() or p.isDataProcess()) :
					continue
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

	#build the conversion functions that the NTMJ templates have
	def __build_conversion_functions__(self) :
		print '	Building conversion functions'
		#make a dictionary of event numbers
		event_numbers = {}
		for c in self.__channel_list :
			event_numbers[c.getName()] = {}
			charge = c.getCharge()
			for p in c.getProcessList() :
				if p.isNTMJProcess() :
					pname = p.getName().split('__')[1]
					event_numbers[c.getName()][pname] = {}
					p.getEventNumbers(charge,event_numbers[c.getName()][pname])
		#sum the event numbers over lepton type per template type
		event_numbers2 = copy.deepcopy(event_numbers)
		for c in self.__channel_list :
			leptype = c.getLepType()
			for p in c.getProcessList() :
				if not p.isNTMJProcess() :
					continue
				process_name = p.getName().split('__')[1]
				for t in p.getTemplateList() :
					if t.getType()=='fit__up' or t.getType()=='fit__down' :
						continue
					for channame in event_numbers2 :
						if channame==c.getName() or not channame.startswith(leptype) :
							continue
						for pname in event_numbers2[channame] :
							if not pname == process_name :
								continue
							for ttype in event_numbers2[channame][pname] :
								if ttype==t.getType() :
									for ttid in event_numbers2[channame][pname][ttype] :
										event_numbers[c.getName()][process_name][ttype][ttid]+=event_numbers2[channame][pname][ttype][ttid]
									break
		for c in self.__channel_list :
			cname = c.getName()
			for p in c.getProcessList() :
				if not p.isNTMJProcess() :
					continue
				pname = p.getName().split('__')[1]
				for t in p.getTemplateList() :
					ttype = t.getType()
					if ttype.find('fit__')!=-1 :
						continue
					s = '		Event numbers for template '+t.getName()+' = {'
					for ttid in event_numbers[cname][pname][ttype] :
						s+="'%s':%.2f,"%(ttid,event_numbers[cname][pname][ttype][ttid])
					s+='}'
					print s
		#Make the conversion functions dictionary per lepton/template type
		conv_func_dict = {}
		for c in self.__channel_list :
			channame = c.getName()
			ltype = c.getLepType()
			if ltype not in conv_func_dict.keys() :
				conv_func_dict[ltype] = {}
			for p in c.getProcessList() :
				if not p.isNTMJProcess() :
					continue
				pname = p.getName()
				ptype = p.getName().split('__')[1]
				if not ptype in conv_func_dict[ltype].keys() :
					conv_func_dict[ltype][ptype] = {}
				for t in p.getTemplateList() :
					ttype = t.getType()
					if ttype=='fit__up' or ttype=='fit__down' :
						continue
					if not ttype in conv_func_dict[ltype][ptype].keys() :
						print '		Adding conversion function for template '+t.getName()
						#make the conversion function
						n_lm1_ps = event_numbers[channame][ptype][ttype]['lm1_ps']
						n_lm1_fs = event_numbers[channame][ptype][ttype]['lm1_fs']
						n_lm2_ps = event_numbers[channame][ptype][ttype]['lm2_ps']
						n_lm2_fs = event_numbers[channame][ptype][ttype]['lm2_fs']
						n_hm_ps  = event_numbers[channame][ptype][ttype]['hm_ps']
						n_hm_fs  = event_numbers[channame][ptype][ttype]['hm_fs']
#						n_lm1_ps = 154.38 
#						n_lm1_fs = 10679.57 
#						n_lm2_ps = 298.36 
#						n_lm2_fs = 7083.99 
#						n_hm_ps  = 4.56 
#						n_hm_fs  = 42.18 
#						if ltype == 'el' :
#							n_lm1_ps = 38.66
#							n_lm1_fs = 2755.10
#							n_lm2_ps = 68.83
#							n_lm2_fs = 1909.33
#							n_hm_ps  = 0.50
#							n_hm_fs  = 24.39
						print '			lm1_ps=%.2f, lm1_fs=%.2f, lm2_ps=%.2f, lm2_fs=%.2f, hm_ps=%.2f, hm_fs=%.2f'%(n_lm1_ps,n_lm1_fs,n_lm2_ps,n_lm2_fs,n_hm_ps,n_hm_fs)
						minn = min(n_lm1_ps,min(n_lm1_fs,min(n_lm2_ps,min(n_lm2_fs,min(n_hm_ps,n_hm_fs)))))
						if minn<0. :
							n_lm1_ps+=-1.*minn+0.5
							n_lm1_fs+=-1.*minn+0.5
							n_lm2_ps+=-1.*minn+0.5
							n_lm2_fs+=-1.*minn+0.5
							n_hm_ps+=-1.*minn+0.5
							n_hm_fs+=-1.*minn+0.5
							print '			modified lm1_ps=%.2f, lm1_fs=%.2f, lm2_ps=%.2f, lm2_fs=%.2f, hm_ps=%.2f, hm_fs=%.2f'%(n_lm1_ps,n_lm1_fs,n_lm2_ps,n_lm2_fs,n_hm_ps,n_hm_fs)	
						x1 = p.getX1(); x2 = p.getX2(); x3 = p.getX3()
						x1_err = p.getX1_err(); x2_err = p.getX2_err(); x3_err = p.getX3_err()
						y1 = n_lm1_ps/n_lm1_fs; y2 = n_lm2_ps/n_lm2_fs; y3 = n_hm_ps/n_hm_fs
						y1_err = y1*sqrt(1./n_lm1_ps+1./n_lm1_fs)
						y2_err = y2*sqrt(1./n_lm2_ps+1./n_lm2_fs)
						y3_err = y3*sqrt(1./n_hm_ps+1./n_hm_fs)
						#Build the TGraph to fit with the conversion function
						n=3
						xs  = array('d',[x1,x2,x3])
						xes = array('d',[x1_err,x2_err,x3_err])
						ys  = array('d',[y1,y2,y3])
						yes = array('d',[y1_err,y2_err,y3_err])
						gr = TGraphErrors(n,xs,ys,xes,yes)
						#Define the conversion function
						conv_func = TF1('conv_func','[0]*x+[1]',x1-x1_err,x3+x3_err)
						#Fit the graph
						gr.Fit('conv_func')
						#get the fitter
						fitter = TVirtualFitter.GetFitter()
						#Build the functions
						nom_func = TF1(ltype+'_'+ttype+'_nom_func','[0]*x+[1]',x1-x1_err,x3+x3_err) 
						nom_func.SetParameter(0,conv_func.GetParameter(0)) 
						nom_func.SetParameter(1,conv_func.GetParameter(1))
						conv_func_dict[ltype][ptype][ttype] = nom_func
						#Make the plot of the fit
						canv = TCanvas(ltype+'_'+ttype+'_conv_func_canv',ltype+' '+ttype+' conversion function canvas',900,900)
						gr.SetTitle(ltype+' '+ttype+' conversion function fit')
						gr.GetXaxis().SetTitle('hadronic top candidate softdrop mass (GeV)')
						gr.GetYaxis().SetTitle('conversion rate (N_{passed}/N_{failed})')
						gr.GetXaxis().SetRangeUser(x1-x1_err-20.,x3+x3_err+20.)
						gr.GetYaxis().SetRangeUser(0.,1.1*nom_func.Eval(x3+x3_err+20.))
						gr.SetMarkerStyle(21)
						gr.Draw('AP')
						nom_func.SetLineWidth(3)
						nom_func.SetLineColor(kBlue)
						nom_func.Draw('L SAME')
						leg = TLegend(0.62,0.67,0.9,0.9)
						leg.AddEntry(gr,'measured rates','PE')
						leg.AddEntry(nom_func,'linear fit','L')
						if ttype == 'nominal' and self.__include_sss :
							fit_up_func = TF1(ltype+'_'+ttype+'_fit_up_func','[0]*x+[1]+sqrt([2]*[2]+2*x*[3]*[3]+x*x*[4]*[4])',x1-x1_err,x3+x3_err) 
							fit_up_func.SetParameter(0,conv_func.GetParameter(0)) 
							fit_up_func.SetParameter(1,conv_func.GetParameter(1))
							fit_up_func.SetParameter(2,conv_func.GetParError(1))
							fit_up_func.SetParameter(3,fitter.GetCovarianceMatrixElement(0,1))
							fit_up_func.SetParameter(4,conv_func.GetParError(0))
							fit_down_func = TF1(ltype+'_'+ttype+'_fit_down_func','[0]*x+[1]-sqrt([2]*[2]+2*x*[3]*[3]+x*x*[4]*[4])',x1-x1_err,x3+x3_err) 
							fit_down_func.SetParameter(0,conv_func.GetParameter(0)) 
							fit_down_func.SetParameter(1,conv_func.GetParameter(1))
							fit_down_func.SetParameter(2,conv_func.GetParError(1))
							fit_down_func.SetParameter(3,fitter.GetCovarianceMatrixElement(0,1))
							fit_down_func.SetParameter(4,conv_func.GetParError(0))
							conv_func_dict[ltype][ptype]['fit__up'] = fit_up_func
							conv_func_dict[ltype][ptype]['fit__down'] = fit_down_func
							fit_up_func.SetLineWidth(3)
							fit_up_func.SetLineStyle(2)
							fit_up_func.SetLineColor(kBlue)
							fit_up_func.Draw('L SAME')
							fit_down_func.SetLineWidth(3)
							fit_down_func.SetLineStyle(2)
							fit_down_func.SetLineColor(kBlue)
							fit_down_func.Draw('L SAME')
							leg.AddEntry(fit_up_func,'fit error','L')
						leg.Draw()
						self.__aux_obj_list.append(copy.deepcopy(canv))
				for t in p.getTemplateList() :
					ttype = t.getType()
					t.setConversionFunction(conv_func_dict[ltype][ptype][ttype])
	
	#build the templates from the NTMJ processes using the conversion functions
	def __build_NTMJ_templates__(self) :
		'	Building NTMJ templates (actually)'
		for c in self.__channel_list :
			for p in c.getProcessList() :
				if p.isNTMJProcess() :
					p.buildTemplates(c.getCharge())

	#set the NNTMJ values after creating the NTMJ templates
	def __set_NNTMJ_values__(self) :
		print '	Getting and resetting NNTMJ values'
		#get the individual values
		nntmj_values = {}
		for c in self.__channel_list :
			nntmj_values[c.getName()] = {}
			for p in c.getProcessList() :
				if not p.isNTMJProcess() :
					continue
				ptype = p.getName().split('__')[1]
				nntmj_values[c.getName()][ptype] = {}
				for t in p.getTemplateList() :
					ttype = t.getType()
					nntmj_values[c.getName()][ptype][ttype] = t.fixNNTMJValues()
		#sum them over channels and NTMJ process types
		nntmj_values2 = copy.deepcopy(nntmj_values)
		for c in self.__channel_list :
			for p in c.getProcessList() :
				if not p.isNTMJProcess() :
					continue
				ptype = p.getName().split('__')[1]
				for t in p.getTemplateList() :
					ttype = t.getType()
					for oc in nntmj_values2 :
						for op in nntmj_values2[oc] :
							if op!=ptype :
								continue
							for ot in nntmj_values2[oc][op] :
								if ot == ttype :
									nntmj_values[c.getName()][ptype][ttype]+=nntmj_values2[oc][op][ot]
									break
		#set the values in the MC templates
		for c in self.__channel_list :
			for p in c.getProcessList() :
				if p.isNTMJProcess() :
					continue		
				for t in p.getTemplateList() :
					ttype = t.getType()
					ptypes = []
					for p in c.getProcessList() :
						if not p.isNTMJProcess() :
							continue
						ptypes.append(p.getName().split('__')[1])
					value = 0.
					for ptype in ptypes :
						value += nntmj_values[c.getName()][ptype][ttype]
					t.setNNTMJValue(value)

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
def make_channel_list(fit_parameter_tuple,include_mu,include_el,sum_charges,include_PDF,include_JEC,include_sss) :
	channel_list = []
	if include_mu :
		if sum_charges :
			channel_list.append(Channel('mu',fit_parameter_tuple,include_PDF,include_JEC,include_sss))
		else :
			channel_list.append(Channel('muplus',fit_parameter_tuple,include_PDF,include_JEC,include_sss))
			channel_list.append(Channel('muminus',fit_parameter_tuple,include_PDF,include_JEC,include_sss))
	if include_el :
		if sum_charges :
			channel_list.append(Channel('el',fit_parameter_tuple,include_PDF,include_JEC,include_sss))
		else :
			channel_list.append(Channel('elplus',fit_parameter_tuple,include_PDF,include_JEC,include_sss))
			channel_list.append(Channel('elminus',fit_parameter_tuple,include_PDF,include_JEC,include_sss))
	return channel_list

#make the list of branches that are only needed for checking selection cuts and nothing else
def make_dict_of_ttree_cut_branches() :
	branches = {}
	#branches needed only for cuts (not anything else)
	branches['muTrig']=Branch('muTrig',None,'i',2)
	branches['muon1_pt']=Branch('muon1_pt',None,'f',-1.0)
	branches['ele1_pt']=Branch('ele1_pt',None,'f',-1.0)
	branches['scaled_lepW_pt']=Branch('scaled_lepW_pt',None,'f',-1.0)
	branches['scaled_hadt_pt']=Branch('scaled_hadt_pt',None,'f',-1.0)
	branches['scaled_hadt_M']=Branch('scaled_hadt_M',None,'f',-1.0)
	branches['muon1_eta']=Branch('muon1_eta',None,'f',100.)
	branches['muon1_ID']=Branch('muon1_ID',None,'I',2)
	branches['muon1_relPt']=Branch('muon1_relPt',None,'f',-1.0)
	branches['muon1_dR']=Branch('muon1_dR',None,'f',-1.0)
	branches['ele1_eta']=Branch('ele1_eta',None,'f',100.)
	branches['ele1_ID']=Branch('ele1_ID',None,'I',2)
	branches['ele1_relPt']=Branch('ele1_relPt',None,'f',-1.0)
	branches['ele1_dR']=Branch('ele1_dR',None,'f',-1.0)
	branches['elTrig']=Branch('elTrig',None,'i',2)
	branches['scaled_lept_M']=Branch('scaled_lept_M',None,'f',-1.0)
	branches['hadt_tau21']=Branch('hadt_tau21',None,'f',-1.0)
	branches['hadt_SDM']=Branch('hadt_SDM',None,'f',-1.0)
	branches['M']=Branch('M',None,'f',-1.0)
	branches['fullselection']=Branch('fullselection',None,'I',2)
	return branches

#return the weight that should be applied to this event when adding to the given process
def get_contrib_weight(filepath,process) :
	#if it's a data or NTMJ process return 1 only if it's a data file
	if process.isDataProcess() or process.isNTMJProcess() :
		if filepath.find('SingleMu')!=-1 and process.getName().startswith('mu')!=-1 :
			return 1.0
		elif filepath.find('SingleEl')!=-1 and process.getName().startswith('el')!=-1 :
			return 1.0
		else :
			return 0.
	#otherwise if it's a MC process
	elif process.isMCProcess() :
		#if it's a qqbar/gg/bck file/process return 1
		if process.getName().find('fqq')!=-1 :
			if filepath.find('qq_semilep_TT')!=-1 :
				return 1.0
			else :
				return 0.
		elif process.getName().find('fgg')!=-1 :
			if filepath.find('gg_semilep_TT')!=-1 :
				return 1.0
			else :
				return 0.
		elif process.getName().find('fbck')!=-1 :
			bck_stems = ['dilep_TT','had_TT','ST_s-c','ST_t-c_top','ST_tW-c_top','ST_t-c_antitop','ST_tW-c_antitop']
			for bckstem in bck_stems :
				if filepath.find(bckstem)!=-1 :
					return 1.0
			return 0.
		else :
			print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
			print '!!!!!!   WARNING, MC PROCESS TYPE NOT RECOGNIZED FROM NAME '+process.getName()+'   !!!!!!'
			print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
			return 0.
	#otherwise I'm super confused about what type of process it is haha
	else :
		print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
		print '!!!!!!   WARNING, PROCESS TYPE NOT RECOGNIZED FROM NAME '+process.getName()+'   !!!!!!'
		print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
		return 0.
