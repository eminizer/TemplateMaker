import copy
from ROOT import TTree
import multiprocessing
from modifier import Fit_Parameter, JEC_Modifier, Simple_Systematic
from branch import Branch
from template import Template

#Functions
PREFAC_1 = '( @NTOT@ - @NBCK@ * #Rbck# - @NWJETS@ * #Rwjets# ) * ( 1. / @NTTBAR@ )'
PREFAC_2 = '( @NTTBAR@ - @NQQBAR@ * #Rqqbar# ) * ( 1. / ( @NTTBAR@ - @NQQBAR@ ) )'
FGG  = '( 1. + #mu# * ( 1. - #mu# ) * ( @NG1@ / ( @NTTBAR@ - @NQQBAR@ ) )'
FGG += ' + ( #mu# * #mu# + #d# * #d# ) * ( ( 1. + #mu# ) * ( @NG2@ / ( @NTTBAR@ - @NQQBAR@ ) )'
FGG += ' + (1. - 5. * #mu# ) * ( @NG3@ / ( @NTTBAR@ - @NQQBAR@ ) )'
FGG += ' + ( #mu# * #mu# + #d# * #d# ) * ( @NG4@ / ( @NTTBAR@ - @NQQBAR@ ) ) ) )'
FQQ  = '( 1. + ( 2. * #mu# + #mu# * #mu# - #d# * #d# ) * ( @NQ1@ / @NQQBAR@ ) + ( #mu# * #mu# + #d# * #d# ) * ( @NQ2@ / @NQQBAR@ ) )'
fgg_func  = PREFAC_1+' * '+PREFAC_2+' * (1. / '+FGG+' ) * ( 1. + #mu# * (1. - #mu# ) * $wg1$'
fgg_func += ' + ( #mu# * #mu# + #d# * #d# ) * ( ( 1. + #mu# ) * $wg2$ + ( 1. - 5. * #mu# ) * $wg3$'
fgg_func += ' + ( #mu# * #mu# + #d# * #d# ) * $wg4$ ) )'
fqq_func  = PREFAC_1+' * ( 1. / '+FQQ+' ) * #Rqqbar# * ( 1. + #Afb# * $wqa0$ + ( 2. * #mu# + #mu# * #mu# - #d# * #d# ) * ( $wqs1$ + #Afb# * $wqa1$ )'
fqq_func += ' + ( #mu# * #mu# + #d# * #d# ) * ( $wqs2$ + #Afb# * $wqa2$ ) )'
fbck_func = '#Rbck#'
fwjets_func = '#Rwjets#'

#Process class
class Process(object) :

	def __init__(self,name,fit_parameter_tuple) :
		#Set the name and base function based on the process name
		self.__name = name
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
		self.__tree_dict = {'sig':TTree(name+'_sig_ptree',name+'_sig_ptree')}

	#Getters/Setters/Adders
	def getName(self) :
		return self.__name
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
			tempnames = ((self.__name+'__par_'+parname+'__up',self.__name+' process par '+parname+' up template'),
						 (self.__name+'__par_'+parname+'__down',self.__name+' process par '+parname+' down template'))
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
			tempnames = ((self.__name+'__'+parname+'__up',self.__name+' process '+parname+' up template'),
						 (self.__name+'__'+parname+'__down',self.__name+' process '+parname+' down template'))
			for tempname in tempnames :
				self.__template_list.append(Template(tempname[0],tempname[1],jecmod))
	#add the templates required by the current list of simple systematic modifiers
	def __add_ss_templates__(self) :
		for ssmod in self.__ss_modifier_list :
			parname = ssmod.getName()
			tempnames = ((self.__name+'__'+parname+'__up',self.__name+' process '+parname+' up template'),
						 (self.__name+'__'+parname+'__down',self.__name+' process '+parname+' down template'))
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
				newtname = tname+'_'+jecmod.getName()+'__up'
				self.__tree_dict[newtname] = TTree(self.__name+'_'+newtname+'_tree',self.__name+'_'+newtname+'_tree')
				newtname = tname+'_'+jecmod.getName()+'__down'
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
	__JEC_names = ['JES','JER']
	#list of constant reweights
	__const_reweights_ttrees = ['weight']
	__const_reweights_ptrees = ['cs_weight']
	#list of systematic reweights
	__simple_systematics = [('sf_pileup','pileup_weight',False),
							('sf_trig_eff','trig_eff_weight',True),
							('sf_lep_ID','lep_ID_weight',True),
							('sf_lep_iso','lep_iso_weight',True),
							('sf_lep_trk','lep_trk_weight',False),
							('sf_btag_eff','btag_eff_weight',False),
							('sf_mu_R','ren_scale_weight',False),
							('sf_mu_F','fact_scale_weight',False),
							('sf_scale_comb','comb_scale_weight',False),
							('sf_pdf_alphas','pdfas_weight',False),
							(None,'lumi',True)]
	#list of event reweighting factors for qqbar and gg distributions
	__qqbar_rws = ['wqs1','wqs2','wqa0','wqa1','wqa2','wqs1_opp','wqs2_opp','wqa0_opp','wqa1_opp','wqa2_opp']
	__gg_rws = ['wg1','wg2','wg3','wg4','wg1_opp','wg2_opp','wg3_opp','wg4_opp']
	__other_rws	= ['wega','wegc']

	def __init__(self,name,fit_parameter_tuple,include_JEC,include_sss) :
		#Set up the Process for this MC_Process
		Process.__init__(self,name,fit_parameter_tuple)
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
	def buildWeightsums(self,channelcharge,ptfn) :
		procs = []
		pipe_ends = []
		for t in self.getTemplateList() :
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

	def buildTemplates(self,channelcharge,ptfn) :
		print '	Building templates for process %s...'%(self.getName())
		procs = []
		pipe_ends = []
		#for each template
		tlist = self.getTemplateList()
		for i in range(len(tlist)) :
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

	#Getters/Setters/Adders
	def getMCProcessList(self) :
		return self.__mc_process_list

	#Functions to add modifiers that apply to the MC Processes
	def __add_mc_fit_parameters__(self) :
		for mc_p in self.__mc_process_list :
			for par in mc_p.getFitParameterList() :
				parname = par.getName()
				if not parname in self.getFitParameterNameList() :
					self.addFitParameter(self,copy.deepcopy(par))
	def __add_pdf_modifiers__(self) :
		for mc_p in self.__mc_process_list :
			for pdfmod in mc_p.getPDFModifierList() :
				pdfmodname = pdfmod.getName()
				if not pdfmodname in self.getPDFModifierNameList() :
					self.addPDFModifier(self,copy.deepcopy(pdfmod))
	def __add_jec_modifiers__(self) :
		for mc_p in self.__mc_process_list :
			for jecmod in mc_p.getJECModifierList() :
				jecmodname = jecmod.getName()
				if not jecmodname in self.getJECModifierNameList() :
					self.addJECModifier(self,copy.deepcopy(jecmod))
	def __add_ss_modifiers__(self) :
		for mc_p in self.__mc_process_list :
			for ssmod in mc_p.getSSModifierList() :
				ssmodname = ssmod.getName()
				if not ssmodname in self.getSSModifierNameList() :
					self.addSSModifier(self,copy.deepcopy(ssmod))

#Data_Process subclass
class Data_Process(Process) :

	def __init__(self,name) :
		Process.__init__(self,name,None)
		#Initialize the tree
		Process.__initialize_trees__(self)

	def buildTemplates(self,channelcharge,ptfn) :
		print '	Building templates for process %s...'%(self.getName())
		procs = []
		pipe_ends = []
		#for each template
		tlist = self.getTemplateList()
		for i in range(len(tlist)) :
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
	elif name == 'DATA' :
		base_function = ''
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


