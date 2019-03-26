import copy
from ROOT import TTree, TFile
import multiprocessing
from modifier import Fit_Parameter, JEC_Modifier, Simple_Systematic
from branch import Branch
from template import Template

#Jet Energy Correction Systematic names
JEC_NAMES = ['JES','JER',
			 #'AK4JESPU',
			 #'AK4JESEta',
			 #'AK4JESPt',
			 #'AK4JESScale',
			 #'AK4JESTime',
			 #'AK4JESFlav',
			 #'AK4JERStat',
			 #'AK4JERSys',
			 #'AK8JESPU',
			 #'AK8JESEta',
			 #'AK8JESPt',
			 #'AK8JESScale',
			 #'AK8JESTime',
			 #'AK8JESFlav',
			 #'AK8JERStat',
			 #'AK8JERSys',
			 ]

#Functions
PREFAC_1 = '( @NTOT@ - @NBCK@ * #Rbck# - @NWJETS@ * #Rwjets# - @NQCD@ * #Rqcd# ) * ( 1. / @NTTBAR@ )'
PREFAC_2 = '( @NTTBAR@ - @NQQBAR@ * #Rqqbar# ) * ( 1. / ( @NTTBAR@ - @NQQBAR@ ) )'
FGG  = '( 1. + #mu# * ( 1. + #mu# ) * ( @NG1@ / ( @NTTBAR@ - @NQQBAR@ ) )'
FGG += ' + ( #mu# * #mu# + #d# * #d# ) * ( ( 1. + #mu# ) * ( @NG2@ / ( @NTTBAR@ - @NQQBAR@ ) )'
FGG += ' + (1. - 5. * #mu# ) * ( @NG3@ / ( @NTTBAR@ - @NQQBAR@ ) )'
FGG += ' + ( #mu# * #mu# + #d# * #d# ) * ( @NG4@ / ( @NTTBAR@ - @NQQBAR@ ) ) ) )'
FQQ  = '( 1. + ( 2. * #mu# + #mu# * #mu# - #d# * #d# ) * ( @NQ1@ / @NQQBAR@ ) + ( #mu# * #mu# + #d# * #d# ) * ( @NQ2@ / @NQQBAR@ ) )'
fgg_func  = PREFAC_1+' * '+PREFAC_2+' * (1. / '+FGG+' ) * ( 1. + #mu# * (1. + #mu# ) * $wg1$'
fgg_func += ' + ( #mu# * #mu# + #d# * #d# ) * ( ( 1. + #mu# ) * $wg2$ + ( 1. - 5. * #mu# ) * $wg3$'
fgg_func += ' + ( #mu# * #mu# + #d# * #d# ) * $wg4$ ) )'
fqq_func  = PREFAC_1+' * ( 1. / '+FQQ+' ) * #Rqqbar# * ( 1. + #Afb# * $wqa0$ + ( 2. * #mu# + #mu# * #mu# - #d# * #d# ) * ( $wqs1$ + #Afb# * $wqa1$ )'
fqq_func += ' + ( #mu# * #mu# + #d# * #d# ) * ( $wqs2$ + #Afb# * $wqa2$ ) )'
fbck_func = '#Rbck#'
fwjets_func = '#Rwjets#'
fqcd_func = '#Rqcd#'
#-------------------------------
fqp0_func = '( 0.5 * (1 + $wqa0$) )'
fqm0_func = '( 0.5 * (1 - $wqa0$) )'
fq0_func  = '( 1. + ( 0.1 * $wqa0$ ) )'
fq1_func  = '( $wqs1$ + ( 0.1 * $wqa1$ ) )'
fq2_func  = '( $wqs2$ + ( 0.1 * $wqa2$ ) )'
fg0_func  = '( 1. )'
fg1_func  = '( $wg1$ )'
fg2_func  = '( $wg2$ )'
fg3_func  = '( $wg3$ )'
fg4_func  = '( $wg4$ )'

#Process class
class Process(object) :

	def __init__(self,name,fit_parameter_tuple) :
		#Set the name and base function based on the process name
		self.__name = name
		self.__leptype = autoset_lepton_type_name(name)
		self.__trig_reg = autoset_trigger_region(name)
		self.__topology = autoset_topology(name)
		self.__base_function = autoset_base_function(name.split('__')[1])
		#Automatically add fit parameters from the base function
		self.__fit_parameter_list = make_fit_parameter_list(self.__base_function,fit_parameter_tuple)
		#modifier lists
		self.__jec_modifier_list = []
		self.__ss_modifier_list  = []
		#template list
		self.__template_list = [Template(name,name+' process nominal template',None)]
		#branch dict
		self.__branch_dict = self.__initialize_branch_dict__()
		#tree dict
		self.__tree_dict = {'sig':TTree(name+'__sig_ptree',name+'__sig_ptree'),
							'qcd_a':TTree(name+'__qcda_ptree',name+'__qcda_ptree'),
							'qcd_b':TTree(name+'__qcdb_ptree',name+'__qcdb_ptree'),
							'qcd_c':TTree(name+'__qcdc_ptree',name+'__qcdc_ptree')}

	#Getters/Setters/Adders
	def getName(self) :
		return self.__name
	def getLepType(self) :
		return self.__leptype
	def getTrigReg(self) :
		return self.__trig_reg
	def getTopology(self) :
		return self.__topology
	def getFitParameterList(self) :
		return self.__fit_parameter_list
	def getFitParameterNameList(self) :
		namelist = []
		for par in self.__fit_parameter_list :
			namelist.append(par.getName())
		return namelist
	def getJECModifierList(self) :
		return self.__jec_modifier_list
	def getJECModifierNameList(self) :
		namelist = []
		for par in self.__jec_modifier_list :
			namelist.append(par.getName())
		return namelist
	def getSSModifierList(self) :
		return self.__ss_modifier_list
	def getSSModifierNameList(self) :
		namelist = []
		for par in self.__ss_modifier_list :
			namelist.append(par.getName())
		return namelist
	def getBranchDict(self) :
		return self.__branch_dict
	def getTreeDict(self) :
		return self.__tree_dict
	def setTree(self,ttid,thistree) :
		self.__tree_dict[ttid] = thistree
	def getJECModList(self) :
		return self.__jec_modifier_list
	def getModifierList(self) :
		return self.__fit_parameter_list+self.__jec_modifier_list+self.__ss_modifier_list
	def getTemplateList(self) :
		return self.__template_list
	def getListOfTTypes(self) :
		ttype_list = []
		for t in self.__template_list :
			ttype_list.append(t.getType())
		return ttype_list
	def getBaseFunction(self) :
		return self.__base_function
	def isMCProcess(self) :
		return isinstance(self,MC_Process)	
	def isDataProcess(self) :
		return isinstance(self,Data_Process)
	def isQCDProcess(self) :
		return isinstance(self,QCD_Process)
	def isFitProcess(self) :
		return isinstance(self,Fit_Process)
	def addFitParameter(self,par) :
		self.__fit_parameter_list.append(par)
	def addJECModifier(self,jecmod) :
		self.__jec_modifier_list.append(jecmod)
	def addSSModifier(self,ssmod) :
		self.__ss_modifier_list.append(ssmod)
	def addBranch(self,branch) :
		self.__branch_dict[branch.getPTreeName()]=branch

	#Private class functions
	#Add templates required by the current list of fit parameters
	def __add_fit_parameter_templates__(self) :
		for fitpar in self.__fit_parameter_list :
			if fitpar.isFixed() :
				continue
			parname = fitpar.getName()
			tempnames = ((self.__name+'__par_'+parname+'Up',self.__name+' process par '+parname+' up template'),
						 (self.__name+'__par_'+parname+'Down',self.__name+' process par '+parname+' down template'))
			for tempname in tempnames :
				self.__template_list.append(Template(tempname[0],tempname[1],fitpar))
	#Make the list of simple systematic modifiers based on the modifier names in the ttrees and process trees
	def __make_ss_modifier_list__(self,ss_list) :
		for i in range(len(ss_list)) :
			ttree_name = ss_list[i][0]
			ptree_name = ss_list[i][1]
			issplit = ss_list[i][2]
			self.__ss_modifier_list.append(Simple_Systematic(ttree_name,ptree_name,issplit))
	#make the list of JEC modifiers based on the names of the wiggled JEC factors
	def __make_JEC_modifier_list__(self,JEC_names) :
		for JEC_name in JEC_names :
			self.__jec_modifier_list.append(JEC_Modifier(JEC_name))
	#add the templates required by the current list of JEC modifiers
	def __add_JEC_templates__(self) :
		for jecmod in self.__jec_modifier_list :
			parname = jecmod.getName()
			tempnames = ((self.__name+'__'+parname+'Up',self.__name+' process '+parname+' up template'),
						 (self.__name+'__'+parname+'Down',self.__name+' process '+parname+' down template'))
			for tempname in tempnames :
				self.__template_list.append(Template(tempname[0],tempname[1],jecmod))
	#add the templates required by the current list of simple systematic modifiers
	def __add_ss_templates__(self) :
		for ssmod in self.__ss_modifier_list :
			parname = ssmod.getName()
			tempnames = ((self.__name+'__'+parname+'Up',self.__name+' process '+parname+' up template'),
						 (self.__name+'__'+parname+'Down',self.__name+' process '+parname+' down template'))
			for tempname in tempnames :
				self.__template_list.append(Template(tempname[0],tempname[1],ssmod))
	#adds branches for observables 
	def __initialize_branch_dict__(self) :
		branch_dict = {}
		branch_dict['cstar'] 		  = Branch('cstar','cstar','f',100.)
		branch_dict['x_F'] 			  = Branch('x_F','x_F','f',100.)
		branch_dict['M'] 			  = Branch('M','M','f',-1.0)
		branch_dict['lep_Q'] 		  = Branch('lep_Q','lep_Q','i',0)
		branch_dict['addTwice'] 	  = Branch('addTwice','addTwice','I',2)
		branch_dict['contrib_weight'] = Branch(None,'contrib_weight','f',1.0)
		return branch_dict
	#add the trees that will be needed by the JEC wiggles to the dictionary
	def __add_JEC_trees__(self) :
		nominal_t_names = self.__tree_dict.keys()
		for jecmod in self.__jec_modifier_list :
			for tname in nominal_t_names :
				newtname = tname+'_'+jecmod.getName()+'Up'
				self.__tree_dict[newtname] = TTree(self.__name+'_'+newtname+'_tree',self.__name+'_'+newtname+'_tree')
				newtname = tname+'_'+jecmod.getName()+'Down'
				self.__tree_dict[newtname] = TTree(self.__name+'_'+newtname+'_tree',self.__name+'_'+newtname+'_tree')
	#Set up the tree by adding the branches
	def __initialize_trees__(self) :
		print '	Initializing trees for process '+self.__name
		for branch in self.__branch_dict.values() :
			print '		Adding branch '+branch.getPTreeName()
			for t in self.__tree_dict.values() :
				t.Branch(branch.getPTreeName(),branch.getPTreeArray(),branch.getPTreeLeafList())

	def __del__(self) :
		pass

	def __str__(self) :
		s = '(name: %s, base_function: %s, \n'%(self.__name,self.__base_function.replace(' ',''))
		s+='					template_list: ['
		for i in range(len(self.__template_list)) :
			temp = self.__template_list[i]
			if i==0 :
				s+=temp.__str__()
			else :
				s+='						'+temp.__str__()
			if i==len(self.__template_list)-1 :
				s+=']'
			else :
				s+=',\n'
		s+=')'
		return s

#MC_Process subclass
class MC_Process(Process) :

	#JEC names
	__JEC_names = JEC_NAMES
	#list of constant reweights
	__const_reweights_ttrees = ['weight']
	__const_reweights_ptrees = ['cs_weight']
	#__const_reweights_ttrees = ['weight','sf_top_pt_rw_hi']
	#__const_reweights_ptrees = ['cs_weight','top_pt_re_weight']
	#list of event reweighting factors for qqbar and gg distributions
	__qqbar_rws = ['wqs1','wqs2','wqa0','wqa1','wqa2','wqs1_opp','wqs2_opp','wqa0_opp','wqa1_opp','wqa2_opp']
	__gg_rws = ['wg1','wg2','wg3','wg4','wg1_opp','wg2_opp','wg3_opp','wg4_opp']
	__other_rws	= ['wega','wegc']

	def __init__(self,name,fit_parameter_tuple,include_JEC,include_sss) :
		#Set up the Process for this MC_Process
		Process.__init__(self,name,fit_parameter_tuple)
		#define the simple systematics
		self.__simple_systematics = [('sf_pileup','pileup_weight',False),
							('sf_trig_eff',self.getLepType()+'_trig_eff_weight_'+self.getTrigReg(),True),
							('sf_lep_ID',self.getLepType()+'_ID_weight',True),
							('sf_lep_iso',self.getLepType()+'_iso_weight',True),
							('sf_btag_eff_flavb','btag_eff_weight_flavb_'+self.getTrigReg(),False),
							('sf_btag_eff_flavc','btag_eff_weight_flavc_'+self.getTrigReg(),False),
							#('sf_btag_eff_heavy','btag_eff_weight_heavy'+self.getTrigReg(),False),
							('sf_btag_eff_light','btag_eff_weight_light_'+self.getTrigReg(),False),
							('sf_ttag_eff_merged','ttag_eff_weight_merged',False),
							('sf_ttag_eff_semimerged','ttag_eff_weight_semimerged',False),
							('sf_ttag_eff_notmerged','ttag_eff_weight_notmerged',False),
							('sf_mu_R','ren_scale_weight',False),
							('sf_mu_F','fact_scale_weight',False),
							('sf_scale_comb','comb_scale_weight',False),
							('sf_pdf_alphas','pdfas_weight',False),
							('sf_B_frag','B_frag_weight',False),
							('sf_B_br','B_br_weight',False),
							('sf_top_pt_rw_v2','top_pt_re_weight',False),
							(None,'lumi',True),
							]
		self.__add_fit_parameter_templates__()
		if include_JEC :
			self.__make_JEC_modifier_list__(self.__JEC_names)
			self.__add_JEC_templates__()
		self.__make_ss_modifier_list__(self.__simple_systematics)
		if include_sss :
			self.__add_ss_templates__()
		#Add to the list of branches
		self.__add_to_dict_of_branches__()
		#Add the JEC trees to the dictionary
		self.__add_JEC_trees__()

	#Public functions
	def buildWeightsums(self,channelcharge,ptfn,n_procs) :
		procs = []
		pipe_ends = []
		for t in self.getTemplateList() :
			if len(procs)>=n_procs :
				for item in pipe_ends :
					item[1].setWeightsumDict(item[0].recv())
				for p in procs :
					p.join()
				procs = []; pipe_ends = []
			p_wsd, c_wsd = multiprocessing.Pipe()
			p = multiprocessing.Process(target=buildWeightsumsParallel, args=(c_wsd,t,copy.deepcopy(self.getTreeDict()),copy.deepcopy(self.getBranchDict()),
																			  self.__const_reweights_ptrees,self.getSSModifierList(),channelcharge,ptfn))
			p.start()
			procs.append(p)
			pipe_ends.append((p_wsd,t))
		for item in pipe_ends :
			item[1].setWeightsumDict(item[0].recv())
		for p in procs :
			p.join()

	def buildTemplates(self,channelcharge,ptfn,n_procs) :
		print '	Building templates for process %s...'%(self.getName())
		procs = []
		pipe_ends = []
		#for each template
		tlist = self.getTemplateList()
		for i in range(len(tlist)) :
			if len(procs)>=n_procs :
				for item in pipe_ends :
					item[1].setHistos(item[0].recv())
				for p in procs :
					p.join()
				procs = []; pipe_ends = []
			t = tlist[i]
			p_hl, c_hl = multiprocessing.Pipe()
			#make the JEC append
			JEC_append = ''
			mod = t.getModifier()
			if mod!=None and mod.isJECModifier() :
				JEC_append = '_'+t.getType()
			p = multiprocessing.Process(target=buildTemplatesParallel, args=(c_hl,t,channelcharge,copy.deepcopy(self.getTreeDict()['sig'+JEC_append]),
																			 copy.deepcopy(self.getBranchDict()),self.getConstantReweightsPTreesList(),
								 											 self.getSSModifierList(),self.getBaseFunction(),self.getFitParameterList(),ptfn,
								 											 i==len(tlist)-1))
			p.start()
			procs.append(p)
			pipe_ends.append((p_hl,t))
		for item in pipe_ends :
			item[1].setHistos(item[0].recv())
		for p in procs :
			p.join()

	#Getters/Setters/Adders
	def getConstantReweightsTTreesList(self) :
		return self.__const_reweights_ttrees
	def getConstantReweightsPTreesList(self) :
		return self.__const_reweights_ptrees
	def getPDFTTreeName(self) :
		return self.__PDF_branch_name
	def getPDFBranchLength(self) :
		return self.__PDF_set_length

	#private class methods
	#add to the branch dict based on a bunch of reweights and observables
	def __add_to_dict_of_branches__(self) :
		#constant reweights
		for i in range(len(self.__const_reweights_ttrees)) :
			self.addBranch(Branch(self.__const_reweights_ttrees[i],self.__const_reweights_ptrees[i],'f',1.0))
		#systematic reweights
		for i in range(len(self.__simple_systematics)) :
			ttree_name = self.__simple_systematics[i][0]
			ptree_name = self.__simple_systematics[i][1]
			if self.__simple_systematics[i][2] :
				if ttree_name!=None :
					self.addBranch(Branch(ttree_name+'_BtoF',ptree_name+'_BtoF','f',1.0))
					self.addBranch(Branch(ttree_name+'_BtoF_hi',ptree_name+'_BtoF_up','f',1.0))
					self.addBranch(Branch(ttree_name+'_BtoF_low',ptree_name+'_BtoF_down','f',1.0))
					self.addBranch(Branch(ttree_name+'_GH',ptree_name+'_GH','f',1.0))
					self.addBranch(Branch(ttree_name+'_GH_hi',ptree_name+'_GH_up','f',1.0))
					self.addBranch(Branch(ttree_name+'_GH_low',ptree_name+'_GH_down','f',1.0))
				else :
					self.addBranch(Branch(None,ptree_name+'_BtoF','f',1.0))
					self.addBranch(Branch(None,ptree_name+'_BtoF_up','f',1.0))
					self.addBranch(Branch(None,ptree_name+'_BtoF_down','f',1.0))
					self.addBranch(Branch(None,ptree_name+'_GH','f',1.0))
					self.addBranch(Branch(None,ptree_name+'_GH_up','f',1.0))
					self.addBranch(Branch(None,ptree_name+'_GH_down','f',1.0))
			else :
				if ttree_name!=None :
					self.addBranch(Branch(ttree_name,ptree_name,'f',1.0))
					self.addBranch(Branch(ttree_name+'_hi',ptree_name+'_up','f',1.0))
					self.addBranch(Branch(ttree_name+'_low',ptree_name+'_down','f',1.0))
				else :
					self.addBranch(Branch(None,ptree_name,'f',1.0))
					self.addBranch(Branch(None,ptree_name+'_up','f',1.0))
					self.addBranch(Branch(None,ptree_name+'_down','f',1.0))
		#event reweights
		if self.getName().find('fqq')!=-1 :
			for i in range(len(self.__qqbar_rws)) :
				rw = self.__qqbar_rws[i]
				self.addBranch(Branch(rw,rw,'f',1.0))	
		elif self.getName().find('fgg')!=-1 :
			for i in range(len(self.__gg_rws)) :
				rw = self.__gg_rws[i]
				self.addBranch(Branch(rw,rw,'f',1.0))

	def __del__(self) :
		pass

	def __str__(self) :
		s = 'MC_Process: (Process: %s)'%(self.__str__())
		return s

#QCD_Process subclass
class QCD_Process(Process) :

	def __init__(self,name,fit_parameter_tuple,mc_plist,include_JEC,include_sss) :
		#Start up the Process
		Process.__init__(self,name,fit_parameter_tuple)
		#Set the list of MC_Processes to subtract from the data events in making this template eventually
		self.__mc_process_list = copy.deepcopy(mc_plist)
		#Add modifiers and templates required by the MC_Processes
		self.__add_mc_fit_parameters__()
		#self.__add_fit_parameter_templates__()
		self.__add_jec_modifiers__()
		#if include_JEC :
		#	self.__add_JEC_templates__()
		self.__add_ss_modifiers__()
		#if include_sss :
		#	self.__add_ss_templates__()
		#Initialize the tree
		self.__initialize_trees__()

	#get the MC-subtracted numbers of data events in the sideband regions for each necessary template type
	def getEventNumbers(self,charge,ptfn,n_procs) :
		print '		Getting event numbers for process '+self.getName()
		manager = multiprocessing.Manager()
		eventnumberssubdict = {}; eventnumberssubdict_p = manager.dict()
		#for each template holding only data events
		procs = []
		for t in self.getTemplateList() :
			if len(procs)>=n_procs :
				for p in procs :
					p.join()
				procs = []
			p = multiprocessing.Process(target=self.__getEventNumbersParallel__, args=(t,charge,ptfn,eventnumberssubdict_p))
			p.start()
			procs.append(p)
		for p in procs :
			p.join()
		for k in eventnumberssubdict_p.keys() :
			eventnumberssubdict[k] = eventnumberssubdict_p[k]
		for t in self.getTemplateList() :
			s= '			Event numbers for template '+t.getName()+' = {'
			for ttid in eventnumberssubdict[t.getType()] :
				s+="'%s':%.2f,"%(ttid,eventnumberssubdict[t.getType()][ttid])
			s+='}'
			print s
		return eventnumberssubdict
	def __getEventNumbersParallel__(self,t,charge,ptfn,eventnumberssubdict) :
		print '			Adding numbers from template '+t.getName()
		#get the event numbers
		eventnumberssubdict[t.getType()] = {'qcd_a':0.,'qcd_b':0.}
		#from the data events in the distribution
		eventnumberssubdict[t.getType()] = t.getEventNumbers(charge,None,self.getTreeDict(),self.getBranchDict(),None,None,eventnumberssubdict[t.getType()],
															 self.getBaseFunction(),self.getFitParameterList(),1.0,ptfn)
		print '				After adding data: %s'%(eventnumberssubdict[t.getType()]) #DEBUG
		#subtract off the MC events in each MC process
		for mcp in self.__mc_process_list :
			print '			Subtracting numbers from MC Process '+mcp.getName()+' for template '+t.getName()
			#make the JEC append
			JEC_append = ''
			mod = t.getModifier()
			if mod!=None and mod.isJECModifier() :
				JEC_append = '_'+t.getType()
			eventnumberssubdict[t.getType()] = t.getEventNumbers(charge,JEC_append,mcp.getTreeDict(),mcp.getBranchDict(),mcp.getConstantReweightsPTreesList(),
																  mcp.getSSModifierList(),eventnumberssubdict[t.getType()],
																  '('+self.getBaseFunction()+')*('+mcp.getBaseFunction()+')',
							  									  self.getFitParameterList()+mcp.getFitParameterList(),-1.0,ptfn)
			print '				After subtracting %s MC: %s'%(mcp.getName(),eventnumberssubdict[t.getType()]) #DEBUG


	#build the templates from the data and MC events
	def buildTemplates(self,channelcharge,ptfn,n_procs) :
		print '	Building templates for QCD process %s...'%(self.getName())
		procs = []
		pipe_ends = []
		#for each template
		tlist = self.getTemplateList()
		for i in range(len(tlist)) :
			if len(procs)>=n_procs :
				for item in pipe_ends :
					item[1].setHistos(item[0].recv())
				for p in procs :
					p.join()
				procs = []; pipe_ends = []
			t = tlist[i]
			p_hl, c_hl = multiprocessing.Pipe()
			p = multiprocessing.Process(target=self.__buildQCDTemplatesParallel__, args=(c_hl,t,channelcharge,ptfn,i==len(tlist)-1))
			p.start()
			procs.append(p)
			pipe_ends.append((p_hl,t))
		for item in pipe_ends :
			item[1].setHistos(item[0].recv())
		for p in procs :
			p.join()
	def __buildQCDTemplatesParallel__(self,child_histolist,t,channelcharge,ptfn,pb) :
		print '			Building template '+t.getName()
		#add the data events
		t.addTreeToTemplates(channelcharge,'qcd_c',self.getTreeDict()['qcd_c'],self.getBranchDict(),None,
							 None,self.getBaseFunction(),self.getFitParameterList(),1.0,ptfn,pb)
		#make the JEC append
		JEC_append = ''
		mod = t.getModifier()
		if mod!=None and mod.isJECModifier() :
			JEC_append = '_'+t.getType()
		#subtract all of the MC events
		for mcp in self.__mc_process_list :
			t.addTreeToTemplates(channelcharge,'qcd_c',mcp.getTreeDict()['qcd_c'+JEC_append],mcp.getBranchDict(),mcp.getConstantReweightsPTreesList(),
								 mcp.getSSModifierList(),'('+self.getBaseFunction()+')*('+mcp.getBaseFunction()+')',
								 self.getFitParameterList()+mcp.getFitParameterList(),-1.0,ptfn,pb)
		child_histolist.send(t.getHistos())
		child_histolist.close()

	#Getters/Setters/Adders
	def setMCProcessList(self,mcpl) :
		self.__mc_process_list = mcpl
	def getMCProcessList(self) :
		return self.__mc_process_list

	#Functions to add modifiers that apply to the MC Processes
	def __add_mc_fit_parameters__(self) :
		for mc_p in self.__mc_process_list :
			for par in mc_p.getFitParameterList() :
				parname = par.getName()
				if not parname in self.getFitParameterNameList() :
					self.addFitParameter(copy.deepcopy(par))
	def __add_jec_modifiers__(self) :
		for mc_p in self.__mc_process_list :
			for jecmod in mc_p.getJECModifierList() :
				jecmodname = jecmod.getName()
				if not jecmodname in self.getJECModifierNameList() :
					self.addJECModifier(copy.deepcopy(jecmod))
	def __add_ss_modifiers__(self) :
		for mc_p in self.__mc_process_list :
			for ssmod in mc_p.getSSModifierList() :
				ssmodname = ssmod.getName()
				if not ssmodname in self.getSSModifierNameList() :
					self.addSSModifier(copy.deepcopy(ssmod))

#Data_Process subclass
class Data_Process(Process) :

	def __init__(self,name) :
		Process.__init__(self,name,None)
		#Initialize the tree
		Process.__initialize_trees__(self)

	def buildTemplates(self,channelcharge,ptfn,n_procs) :
		print '	Building templates for process %s...'%(self.getName())
		procs = []
		pipe_ends = []
		#for each template
		tlist = self.getTemplateList()
		for i in range(len(tlist)) :
			if len(procs)>=n_procs :
				for item in pipe_ends :
					item[1].setHistos(item[0].recv())
				for p in procs :
					p.join()
				procs = []; pipe_ends = []
			t = tlist[i]
			p_hl, c_hl = multiprocessing.Pipe()
			p = multiprocessing.Process(target=buildTemplatesParallel,args=(c_hl,t,channelcharge,copy.deepcopy(self.getTreeDict()['sig']),
																			copy.deepcopy(self.getBranchDict()),None,None,self.getBaseFunction(),
																			self.getFitParameterList(),ptfn,i==len(tlist)-1))
			p.start()
			procs.append(p)
			pipe_ends.append((p_hl,t))
		for item in pipe_ends :
			item[1].setHistos(item[0].recv())
		for p in procs :
			p.join()

	def __del__(self) :
		pass

	def __str__(self) :
		s = 'Data_Process: (Process: %s)'%(self.__str__())
		return s

#Fit_Process subclass
class Fit_Process(Process) :

	#JEC names
	__JEC_names = JEC_NAMES
	#list of constant reweights
	__const_reweights_ttrees = ['weight']
	__const_reweights_ptrees = ['cs_weight']
	#__const_reweights_ttrees = ['weight','sf_top_pt_rw_hi']
	#__const_reweights_ptrees = ['cs_weight','top_pt_re_weight']
	#list of event reweighting factors for qqbar and gg distributions
	__qqbar_rws = ['wqs1','wqs2','wqa0','wqa1','wqa2','wqs1_opp','wqs2_opp','wqa0_opp','wqa1_opp','wqa2_opp']
	__gg_rws = ['wg1','wg2','wg3','wg4','wg1_opp','wg2_opp','wg3_opp','wg4_opp']
	__other_rws	= ['wega','wegc']

	def __init__(self,name,include_JEC,include_sss) :
		#Set up the Process for this MC_Process
		Process.__init__(self,name,None)
		#define the simple systematics
		self.__simple_systematics = [('sf_pileup','pileup_weight',False),
							('sf_trig_eff',self.getLepType()+'_trig_eff_weight_'+self.getTrigReg(),True),
							('sf_lep_ID',self.getLepType()+'_ID_weight',True),
							('sf_lep_iso',self.getLepType()+'_iso_weight',True),
							('sf_btag_eff_flavb','btag_eff_weight_flavb_'+self.getTrigReg(),False),
							('sf_btag_eff_flavc','btag_eff_weight_flavc_'+self.getTrigReg(),False),
							#('sf_btag_eff_heavy','btag_eff_weight_heavy'+self.getTrigReg(),False),
							('sf_btag_eff_light','btag_eff_weight_light_'+self.getTrigReg(),False),
							('sf_ttag_eff_merged','ttag_eff_weight_merged',False),
							('sf_ttag_eff_semimerged','ttag_eff_weight_semimerged',False),
							('sf_ttag_eff_notmerged','ttag_eff_weight_notmerged',False),
							('sf_mu_R','ren_scale_weight',False),
							('sf_mu_F','fact_scale_weight',False),
							('sf_scale_comb','comb_scale_weight',False),
							('sf_pdf_alphas','pdfas_weight',False),
							('sf_B_frag','B_frag_weight',False),
							('sf_B_br','B_br_weight',False),
							('sf_top_pt_rw_v2','top_pt_re_weight',False),
							(None,'lumi',True),
							]
		if include_JEC :
			self.__make_JEC_modifier_list__(self.__JEC_names)
			self.__add_JEC_templates__()
		self.__make_ss_modifier_list__(self.__simple_systematics)
		if include_sss :
			self.__add_ss_templates__()
		#Add to the list of branches
		self.__add_to_dict_of_branches__()
		#Add the JEC trees to the dictionary
		self.__add_JEC_trees__()

	def buildTemplates(self,channelcharge,ptfn,n_procs) :
		print '	Building templates for fit process %s...'%(self.getName())
		procs = []
		pipe_ends = []
		#for each template
		tlist = self.getTemplateList()
		for i in range(len(tlist)) :
			if len(procs)>=n_procs :
				for item in pipe_ends :
					item[1].setHistos(item[0].recv())
				for p in procs :
					p.join()
				procs = []; pipe_ends = []
			t = tlist[i]
			p_hl, c_hl = multiprocessing.Pipe()
			#make the JEC append
			JEC_append = ''
			mod = t.getModifier()
			if mod!=None and mod.isJECModifier() :
				JEC_append = '_'+t.getType()
			p = multiprocessing.Process(target=buildTemplatesParallel, args=(c_hl,t,channelcharge,copy.deepcopy(self.getTreeDict()['sig'+JEC_append]),
																			 copy.deepcopy(self.getBranchDict()),self.getConstantReweightsPTreesList(),
								 											 self.getSSModifierList(),self.getBaseFunction(),self.getFitParameterList(),ptfn,
								 											 i==len(tlist)-1))
			p.start()
			procs.append(p)
			pipe_ends.append((p_hl,t))
		for item in pipe_ends :
			item[1].setHistos(item[0].recv())
		for p in procs :
			p.join()

	#Getters/Setters/Adders
	def getConstantReweightsTTreesList(self) :
		return self.__const_reweights_ttrees
	def getConstantReweightsPTreesList(self) :
		return self.__const_reweights_ptrees
	def getPDFTTreeName(self) :
		return self.__PDF_branch_name
	def getPDFBranchLength(self) :
		return self.__PDF_set_length

	#private class methods
	#add to the branch dict based on a bunch of reweights and observables
	def __add_to_dict_of_branches__(self) :
		#constant reweights
		for i in range(len(self.__const_reweights_ttrees)) :
			self.addBranch(Branch(self.__const_reweights_ttrees[i],self.__const_reweights_ptrees[i],'f',1.0))
		#systematic reweights
		for i in range(len(self.__simple_systematics)) :
			ttree_name = self.__simple_systematics[i][0]
			ptree_name = self.__simple_systematics[i][1]
			if self.__simple_systematics[i][2] :
				if ttree_name!=None :
					self.addBranch(Branch(ttree_name+'_BtoF',ptree_name+'_BtoF','f',1.0))
					self.addBranch(Branch(ttree_name+'_BtoF_hi',ptree_name+'_BtoF_up','f',1.0))
					self.addBranch(Branch(ttree_name+'_BtoF_low',ptree_name+'_BtoF_down','f',1.0))
					self.addBranch(Branch(ttree_name+'_GH',ptree_name+'_GH','f',1.0))
					self.addBranch(Branch(ttree_name+'_GH_hi',ptree_name+'_GH_up','f',1.0))
					self.addBranch(Branch(ttree_name+'_GH_low',ptree_name+'_GH_down','f',1.0))
				else :
					self.addBranch(Branch(None,ptree_name+'_BtoF','f',1.0))
					self.addBranch(Branch(None,ptree_name+'_BtoF_up','f',1.0))
					self.addBranch(Branch(None,ptree_name+'_BtoF_down','f',1.0))
					self.addBranch(Branch(None,ptree_name+'_GH','f',1.0))
					self.addBranch(Branch(None,ptree_name+'_GH_up','f',1.0))
					self.addBranch(Branch(None,ptree_name+'_GH_down','f',1.0))
			else :
				if ttree_name!=None :
					self.addBranch(Branch(ttree_name,ptree_name,'f',1.0))
					self.addBranch(Branch(ttree_name+'_hi',ptree_name+'_up','f',1.0))
					self.addBranch(Branch(ttree_name+'_low',ptree_name+'_down','f',1.0))
				else :
					self.addBranch(Branch(None,ptree_name,'f',1.0))
					self.addBranch(Branch(None,ptree_name+'_up','f',1.0))
					self.addBranch(Branch(None,ptree_name+'_down','f',1.0))
		#event reweights
		if self.getName().find('fq')!=-1 :
			for i in range(len(self.__qqbar_rws)) :
				rw = self.__qqbar_rws[i]
				self.addBranch(Branch(rw,rw,'f',1.0))	
		elif self.getName().find('fg')!=-1 :
			for i in range(len(self.__gg_rws)) :
				rw = self.__gg_rws[i]
				self.addBranch(Branch(rw,rw,'f',1.0))

	def __del__(self) :
		pass

	def __str__(self) :
		s = 'Fit_Process: (Process: %s)'%(self.__str__())
		return s

#Return the choice of base function according to process name
def autoset_base_function(name) :
	base_function = None
	if name == 'fgg' :
		base_function = fgg_func
	elif name == 'fqq' :
		base_function = fqq_func
	elif name == 'fbck' :
		base_function = fbck_func
	elif name == 'fwjets' :
		base_function = fwjets_func
	elif name == 'fqcd' :
		base_function = fqcd_func
	elif name == 'data_obs' :
		base_function = ''
	elif name == 'fqp0' :
		base_function = fqp0_func
	elif name == 'fqm0' :
		base_function = fqm0_func
	elif name == 'fq0' :
		base_function = fq0_func
	elif name == 'fq1' :
		base_function = fq1_func
	elif name == 'fq2' :
		base_function = fq2_func
	elif name == 'fg0' :
		base_function = fg0_func
	elif name == 'fg1' :
		base_function = fg1_func
	elif name == 'fg2' :
		base_function = fg2_func
	elif name == 'fg3' :
		base_function = fg3_func
	elif name == 'fg4' :
		base_function = fg4_func
	else :
		print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
		print '!!!!!!   WARNING, PROCESS TYPE NOT RECOGNIZED FROM NAME: '+name+'   !!!!!!'
		print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	return base_function

#make the list of Fit Parameter objects based on the tuple created by input options
def make_fit_parameter_list(func,fit_parameter_tuple) :
	fit_parameter_list = []
	if fit_parameter_tuple!=None :
		fs = func.split('#')
		for fit_par in fit_parameter_tuple :
			if fit_par[0] in fs :
				fit_parameter_list.append(Fit_Parameter(fit_par[0],fit_par[1],fit_par[2],fit_par[3]))
	return fit_parameter_list

def buildWeightsumsParallel(child_weightsumdict,t,ttd,bd,crws,ssrws,channelcharge,ptfn) :
	JEC_append = ''
	mod = t.getModifier()
	if mod!=None and mod.isJECModifier() :
		JEC_append = '_'+t.getType()
	print '		Doing template '+t.getName()
	t.buildWeightsums(ttd,bd,crws,ssrws,JEC_append,channelcharge,ptfn)
	child_weightsumdict.send(t.getWeightsumDict())
	child_weightsumdict.close()

def buildTemplatesParallel(child_histolist,t,channelcharge,thistree,bd,crws,ssrws,bf,fpl,ptfn,pb) :
	print '		Building template '+t.getName()
	#add all of the MC events
	t.addTreeToTemplates(channelcharge,'sig',thistree,bd,crws,ssrws,bf,fpl,1.0,ptfn,pb)
	child_histolist.send(t.getHistos())
	child_histolist.close()

#Automatically returns the lepton type name based on the process name
def autoset_lepton_type_name(name) :
	leptontypename = ''
	if name.split('__')[0].split('_')[1].startswith('mu') :
		leptontypename = 'mu'
	elif name.split('__')[0].split('_')[1].startswith('el') :
		leptontypename = 'el'
	else :
		print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
		print '!!!!!!   WARNING, LEPTON TYPE NOT RECOGNIZED FROM PROCESS '+name+'   !!!!!!'
		print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	return leptontypename

#Automatically returns the trigger region based on the process name
def autoset_trigger_region(name) :
	trigreg = ''
	if name.split('__')[0].split('_')[0] in ['t1','t2'] :
		trigreg = 'b'
	elif name.split('__')[0].split('_')[0] in ['t3'] :
		trigreg = 'r'
	else :
		print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
		print '!!!!!!   WARNING, TRIGGER REGION NOT RECOGNIZED FROM PROCESS '+name+'   !!!!!!'
		print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	return trigreg

#Returns the decay topology based on the process name
def autoset_topology(name) :
	top = ''
	if name.split('__')[0].split('_')[0] == 't1' :
		top = 't1'
	elif name.split('__')[0].split('_')[0] == 't2' :
		top = 't2'
	elif name.split('__')[0].split('_')[0] == 't3' :
		top = 't3'
	else :
		print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
		print '!!!!!!     WARNING, topology NOT RECOGNIZED FROM PROCESS '+name+'    !!!!!!'
		print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	return top


