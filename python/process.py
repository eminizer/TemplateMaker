import copy
from ROOT import TTree
from modifier import Fit_Parameter, PDF_Modifier, JEC_Modifier, Simple_Systematic
from branch import Branch
from template import Template

#Functions
PREFAC_1 = '( @NTOT@ - @NBCK@ * #Rbck# - @NNTMJ@ * #Rntmj# ) * ( 1. / @NTTBAR@ )'
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
fntmj_func = '#Rntmj#'

#Cuts for conversion function calculation
LM1_LOW = 25.
LM1_HI  = 65.
LM2_LOW = 65.
LM2_HI  = 105.
SM_LOW  = 105.
SM_HI   = 220.
HM_LOW  = 220.
HM_HI   = 750.
SUBS_CUT = 0.69

#Process class
class Process(object) :

	def __init__(self,name,fit_parameter_tuple) :
		#Set the name and base function based on the process name
		self.__name = name
		self.__base_function = autoset_base_function(name.split('__')[1])
		#Automatically add fit parameters from the base function
		self.__fit_parameter_list = make_fit_parameter_list(self.__base_function,fit_parameter_tuple)
		#modifier lists
		self.__pdf_modifier_list = []
		self.__jec_modifier_list = []
		self.__ss_modifier_list  = []
		#template list
		self.__template_list = [Template(name,name+' process nominal template',None)]
		#branch dict
		self.__branch_dict = self.__initialize_branch_dict__()
		#tree dict
		self.__tree_dict = {'lm1_ps':TTree(name+'_lm1_ps_tree',name+'_lm1_ps_tree'),
							'lm1_fs':TTree(name+'_lm1_fs_tree',name+'_lm1_fs_tree'),
							'lm2_ps':TTree(name+'_lm2_ps_tree',name+'_lm2_ps_tree'),
							'lm2_fs':TTree(name+'_lm2_fs_tree',name+'_lm2_fs_tree'),
							'sm_ps':TTree(name+'_sm_ps_tree',name+'_sm_ps_tree'),
							'sm_fs':TTree(name+'_sm_fs_tree',name+'_sm_fs_tree'),
							'hm_ps':TTree(name+'_hm_ps_tree',name+'_hm_ps_tree'),
							'hm_fs':TTree(name+'_hm_fs_tree',name+'_hm_fs_tree')}

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
	def getPDFModifierList(self) :
		return self.__pdf_modifier_list
	def getPDFModifierNameList(self) :
		namelist = []
		for par in self.__pdf_modifier_list :
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
		return self.__fit_parameter_list+self.__pdf_modifier_list+self.__jec_modifier_list+self.__ss_modifier_list
	def getTemplateList(self) :
		return self.__template_list
	def getListOfTTypes(self) :
		ttype_list = []
		for t in self.__template_list :
			ttype_list.append(t.getType())
		return ttype_list
	def getBaseFunction(self) :
		return self.__base_function
	def getX1(self) :
		return (LM1_HI+LM1_LOW)/2.
	def getX2(self) :
		return (LM2_HI+LM2_LOW)/2.
	def getX3(self) :
		return (HM_HI+HM_LOW)/2.
	def getX1_err(self) :
		return (LM1_HI-LM1_LOW)/2.
	def getX2_err(self) :
		return (LM2_HI-LM2_LOW)/2.
	def getX3_err(self) :
		return (HM_HI-HM_LOW)/2.
	def isMCProcess(self) :
		return isinstance(self,MC_Process)
	def isNTMJProcess(self) :
		return isinstance(self,NTMJ_Process)	
	def isDataProcess(self) :
		return isinstance(self,Data_Process)
	def addFitParameter(self,par) :
		self.__fit_parameter_list.append(par)
	def addPDFModifier(self,pdfmod) :
		self.__pdf_modifier_list.append(pdfmod)
	def addJECModifier(self,jecmod) :
		self.__jec_modifier_list.append(jecmod)
	def addSSModifier(self,ssmod) :
		self.__ss_modifier_list.append(ssmod)
	def addBranch(self,branch) :
		self.__branch_dict[branch.getPTreeName()]=branch
	def fillTree(self,ttree_file_path) :
		treestring = ''
		hadtSDM = self.__branch_dict['hadt_SDM'].getPTreeValue()
		if (hadtSDM>50. and hadtSDM<LM1_LOW) or hadtSDM>HM_HI :
			return
		elif hadtSDM>LM1_LOW and hadtSDM<=LM1_HI :
			treestring+='lm1'
		elif hadtSDM>LM2_LOW and hadtSDM<=LM2_HI :
			treestring+='lm2'
		elif hadtSDM>SM_LOW and hadtSDM<=SM_HI :
			treestring+='sm'
		elif hadtSDM>HM_LOW and hadtSDM<=HM_HI :
			treestring+='hm'
		else :
			print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
			print '!!!!!!   MASS REGION NOT IDENTIFIABLE, HADT_M = '+str(hadtSDM)+'   !!!!!!'
			print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
		treestring+='_'
		hadttau32 = self.__branch_dict['hadt_tau32'].getPTreeValue()
		if hadttau32>SUBS_CUT :
			treestring+='fs'
		else :
			treestring+='ps'
		#If it's a MC distribution, add the JEC part of the name to the tree identifier if necessary
		if self.isMCProcess() and (ttree_file_path.find('_up')!=-1 or ttree_file_path.find('_down')!=-1) :
			JEC_ID = ''
			for jecmod in self.__jec_modifier_list :
				if ttree_file_path.find(jecmod.getName()+'_up') != -1 :
					JEC_ID = jecmod.getName()+'__up'
					break
				elif ttree_file_path.find(jecmod.getName()+'_down') != -1 :
					JEC_ID = jecmod.getName()+'__down'
					break
			treestring+='_'+JEC_ID
		self.__tree_dict[treestring].Fill()


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
	#Make the list of pdf modifiers based on the name and size of the pdf reweight branch in the ttrees
	def __make_PDF_modifier_list__(self,PDF_branch_name,PDF_set_length) :
		for i in range(PDF_set_length) :
			ipdflambda = i+1
			self.__pdf_modifier_list.append(PDF_Modifier(PDF_branch_name,ipdflambda))
	#Make the list of simple systematic modifiers based on the modifier names in the ttrees and process trees
	def __make_ss_modifier_list__(self,ss_name_list_ttrees,ss_name_list_ptrees) :
		for i in range(len(ss_name_list_ttrees)) :
			ttree_name = ss_name_list_ttrees[i]
			ptree_name = ss_name_list_ptrees[i]
			self.__ss_modifier_list.append(Simple_Systematic(ttree_name,ptree_name))
	#add the templates required by the current list of pdf modifiers
	def __add_PDF_templates__(self) :
		for pdfmod in self.__pdf_modifier_list :
			parname = pdfmod.getName()
			tempnames = ((self.__name+'__'+parname+'__up',self.__name+' process '+parname+' up template'),
						 (self.__name+'__'+parname+'__down',self.__name+' process '+parname+' down template'))
			for tempname in tempnames :
				self.__template_list.append(Template(tempname[0],tempname[1],pdfmod))
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
		branch_dict['hadt_M'] 		  = Branch('scaled_hadt_M','hadt_M','f',-1.0)
		branch_dict['hadt_SDM'] 	  = Branch('hadt_SDM','hadt_SDM','f',-1.0)
		branch_dict['hadt_tau32'] 	  = Branch('hadt_tau32','hadt_tau32','f',-1.0)
		branch_dict['Q_l'] 			  = Branch('Q_l','Q_l','i',0)
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

	#PDF name
	__PDF_branch_name = 'sf_pdf_alphas'
	__PDF_set_length  = 1
	#JEC names
	__JEC_names = ['JES','JER']
	#list of constant reweights
	__const_reweights_ttrees = ['weight']
	__const_reweights_ptrees = ['cs_weight']
	#list of systematic reweights
	__simple_systematics_ttrees = ['sf_pileup',     'sf_lep_ID',     'sf_mu_R', 		 'sf_mu_F',			  'sf_scale_comb', 	   'sf_pdf_alphas', None]
	__simple_systematics_ptrees = ['pileup_weight', 'lep_ID_weight', 'ren_scale_weight', 'fact_scale_weight', 'comb_scale_weight', 'pdfas_weight',  'luminosity']
	#list of event reweighting factors
	__event_reweights = ['wg1','wg2','wg3','wg4','wqs1','wqs2','wqa0','wqa1','wqa2',
						 'wg1_opp','wg2_opp','wg3_opp','wg4_opp','wqs1_opp','wqs2_opp','wqa0_opp','wqa1_opp','wqa2_opp',
						 'wega','wegc']

	def __init__(self,name,fit_parameter_tuple,include_PDF,include_JEC,include_sss) :
		#Set up the Process for this MC_Process
		Process.__init__(self,name,fit_parameter_tuple)
		Process.__add_fit_parameter_templates__(self)
		Process.__make_PDF_modifier_list__(self,self.__PDF_branch_name,self.__PDF_set_length)
		if include_PDF :
			Process.__add_PDF_templates__(self)
		Process.__make_JEC_modifier_list__(self,self.__JEC_names)
		if include_JEC :
			Process.__add_JEC_templates__(self)
		Process.__make_ss_modifier_list__(self,self.__simple_systematics_ttrees,self.__simple_systematics_ptrees)
		if include_sss :
			Process.__add_ss_templates__(self)
		#Add to the list of branches
		self.__add_to_dict_of_branches__()
		#Add the JEC trees to the dictionary
		Process.__add_JEC_trees__(self)
		#Initialize the tree
		Process.__initialize_trees__(self)

	#Public functions
	def buildWeightsums(self,channelcharge) :
		for t in Process.getTemplateList(self) :
			JEC_append = ''
			mod = t.getModifier()
			if mod!=None and mod.isJECModifier() :
				JEC_append = '_'+t.getType()
			print '		Doing template '+t.getName()
			t.buildWeightsums(Process.getTreeDict(self),Process.getBranchDict(self),self.__const_reweights_ptrees,Process.getSSModifierNameList(self),JEC_append,channelcharge)

	def buildTemplates(self,channelcharge) :
		#for each template
		for t in Process.getTemplateList(self) :
			print '	Building template '+t.getName()
			#make the JEC append
			JEC_append = ''
			mod = t.getModifier()
			if mod!=None and mod.isJECModifier() :
				JEC_append = '_'+t.getType()
			#add all of the MC events
			t.addTreeToTemplates(channelcharge,'sm_ps',Process.getTreeDict(self)['sm_ps'+JEC_append],Process.getBranchDict(self),self.getConstantReweightsPTreesList(),
								 Process.getSSModifierNameList(self),Process.getBaseFunction(self),Process.getFitParameterList(self),1.0)


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
		#PDF reweights
		Process.addBranch(self,Branch(self.__PDF_branch_name,'pdf_weights',str(self.__PDF_set_length*2+1)+'/f',1.0))
		#constant reweights
		for i in range(len(self.__const_reweights_ttrees)) :
			Process.addBranch(self,Branch(self.__const_reweights_ttrees[i],self.__const_reweights_ptrees[i],'f',1.0))
		#systematic reweights
		for i in range(len(self.__simple_systematics_ttrees)) :
			ttree_name_nominal = self.__simple_systematics_ttrees[i]
			ttree_name_up = None; ttree_name_down = None;
			if ttree_name_nominal!=None :
				ttree_name_up=ttree_name_nominal+'_hi'; ttree_name_down=ttree_name_nominal+'_low'
			ptree_name_nominal = self.__simple_systematics_ptrees[i]
			ptree_name_up = None; ptree_name_down = None
			if ptree_name_nominal!=None :
				ptree_name_up = ptree_name_nominal+'_up'; ptree_name_down = ptree_name_nominal+'_down'
			Process.addBranch(self,Branch(ttree_name_nominal,ptree_name_nominal,'f',1.0))
			Process.addBranch(self,Branch(ttree_name_up,ptree_name_up,'f',1.0))
			Process.addBranch(self,Branch(ttree_name_down,ptree_name_down,'f',1.0))
		#event reweights
		for i in range(len(self.__event_reweights)) :
			rw = self.__event_reweights[i]
			Process.addBranch(self,Branch(rw,rw,'f',1.0))

	def __del__(self) :
		pass

	def __str__(self) :
		s = 'MC_Process: (Process: %s)'%(Process.__str__(self))
		return s

#NTMJ_Process subclass
class NTMJ_Process(Process) :

	def __init__(self,name,fit_parameter_tuple,mc_plist,include_PDF,include_JEC,include_sss) :
		#Start up the Process
		Process.__init__(self,name,fit_parameter_tuple)
		#Set the list of MC_Processes to subtract from the data events in making this template eventually
		self.__mc_process_list = mc_plist
		#Add modifiers and templates required by the MC_Processes
		self.__add_mc_fit_parameters__()
		Process.__add_fit_parameter_templates__(self)
		self.__add_pdf_modifiers__()
		if include_PDF :
			Process.__add_PDF_templates__(self)
		self.__add_jec_modifiers__()
		if include_JEC :
			Process.__add_JEC_templates__(self)
		self.__add_ss_modifiers__()
		Process.addSSModifier(self,Simple_Systematic(None,'fit')) #adding for NTMJ fit function systematics
		if include_sss :
			Process.__add_ss_templates__(self)
		#Initialize the tree
		Process.__initialize_trees__(self)

	#get the MC-subtracted numbers of data events in the sideband regions for each necessary template type
	def getEventNumbers(self,charge,eventnumberssubdict) :
		print '		Getting event numbers for process '+Process.getName(self)
		#for each template holding only data events
		for t in Process.getTemplateList(self) :
			if t.getType()=='fit__up' or t.getType()=='fit__down' :
				continue
			print '			Adding numbers from template '+t.getName()
			#get the event numbers
			eventnumberssubdict[t.getType()] = {'lm1_ps':0.,'lm1_fs':0.,'lm2_ps':0.,'lm2_fs':0.,'hm_ps':0.,'hm_fs':0.,}
			#from the data events in the distribution
			t.getEventNumbers(charge,None,Process.getTreeDict(self),Process.getBranchDict(self),None,None,
							  eventnumberssubdict[t.getType()],Process.getBaseFunction(self),Process.getFitParameterList(self),1.0)
		#subtract off the MC events in each MC process
		for mcp in self.__mc_process_list :
			for t in mcp.getTemplateList() :
				print '			Adding numbers from template '+t.getName()
				#make the JEC append
				JEC_append = ''
				mod = t.getModifier()
				if mod!=None and mod.isJECModifier() :
					JEC_append = '_'+t.getType()
				t.getEventNumbers(charge,JEC_append,mcp.getTreeDict(),mcp.getBranchDict(),mcp.getConstantReweightsPTreesList(),
								  mcp.getSSModifierNameList(),eventnumberssubdict[t.getType()],'('+Process.getBaseFunction(self)+')*('+mcp.getBaseFunction()+')',
								  Process.getFitParameterList(self)+mcp.getFitParameterList(),-1.0)
		for t in Process.getTemplateList(self) :
			if t.getType().find('fit__')!=-1 :
				continue
			s= '			Event numbers for template '+t.getName()+' = {'
			for ttid in eventnumberssubdict[t.getType()] :
				s+="'%s':%.2f,"%(ttid,eventnumberssubdict[t.getType()][ttid])
			s+='}'
			print s

	#build the NTMJ templates from the data and MC events
	def buildTemplates(self,channelcharge) :
		#for each template
		for t in Process.getTemplateList(self) :
			print '			Building template '+t.getName()
			#add the data events
			t.addTreeToTemplates(channelcharge,'sm_fs',Process.getTreeDict(self)['sm_fs'],Process.getBranchDict(self),None,
								 None,Process.getBaseFunction(self),Process.getFitParameterList(self),1.0)
			#make the JEC append
			JEC_append = ''
			mod = t.getModifier()
			if mod!=None and mod.isJECModifier() :
				JEC_append = '_'+t.getType()
			#subtract all of the MC events
			for mcp in self.__mc_process_list :
				t.addTreeToTemplates(channelcharge,'sm_fs',mcp.getTreeDict()['sm_fs'+JEC_append],mcp.getBranchDict(),mcp.getConstantReweightsPTreesList(),
									 mcp.getSSModifierNameList(),'('+Process.getBaseFunction(self)+')*('+mcp.getBaseFunction()+')',
									 Process.getFitParameterList(self)+mcp.getFitParameterList(),-1.0)

	#Getters/Setters/Adders
	def getMCProcessList(self) :
		return self.__mc_process_list

	#Functions to add modifiers that apply to the MC Processes
	def __add_mc_fit_parameters__(self) :
		for mc_p in self.__mc_process_list :
			for par in mc_p.getFitParameterList() :
				parname = par.getName()
				if not parname in Process.getFitParameterNameList(self) :
					Process.addFitParameter(self,copy.deepcopy(par))
	def __add_pdf_modifiers__(self) :
		for mc_p in self.__mc_process_list :
			for pdfmod in mc_p.getPDFModifierList() :
				pdfmodname = pdfmod.getName()
				if not pdfmodname in Process.getPDFModifierNameList(self) :
					Process.addPDFModifier(self,copy.deepcopy(pdfmod))
	def __add_jec_modifiers__(self) :
		for mc_p in self.__mc_process_list :
			for jecmod in mc_p.getJECModifierList() :
				jecmodname = jecmod.getName()
				if not jecmodname in Process.getJECModifierNameList(self) :
					Process.addJECModifier(self,copy.deepcopy(jecmod))
	def __add_ss_modifiers__(self) :
		for mc_p in self.__mc_process_list :
			for ssmod in mc_p.getSSModifierList() :
				ssmodname = ssmod.getName()
				if not ssmodname in Process.getSSModifierNameList(self) :
					Process.addSSModifier(self,copy.deepcopy(ssmod))

	def __del__(self) :
		pass

	def __str__(self) :
		s = 'NTMJ_Process: (Process: %s, MC_process_list (names): ['%(Process.__str__(self))
		for i in range(len(self.__mc_process_list)) :
			mcp = self.__mc_process_list[i]
			s+=mcp.getName()
			if i!=len(self.__mc_process_list)-1 :
				s+=', '
			else :
				s+=']'
		s+=')'
		return s

#Data_Process subclass
class Data_Process(Process) :

	def __init__(self,name) :
		Process.__init__(self,name,None)
		#Initialize the tree
		Process.__initialize_trees__(self)

	def buildTemplates(self,channelcharge) :
		#for each template
		for t in Process.getTemplateList(self) :
			print '	Building template '+t.getName()
			#add all of the MC events
			t.addTreeToTemplates(channelcharge,'sm_ps',Process.getTreeDict(self)['sm_ps'],Process.getBranchDict(self),None,
								 None,Process.getBaseFunction(self),Process.getFitParameterList(self),1.0)

	def __del__(self) :
		pass

	def __str__(self) :
		s = 'Data_Process: (Process: %s)'%(Process.__str__(self))
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
	elif name == 'fntmj' :
		base_function = fntmj_func
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


