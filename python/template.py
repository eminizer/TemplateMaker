from ROOT import TH3D, TH1D
from array import array
from math import *

#Template class
class Template(object) :

	#histogram limits
	__XBINS = array('d',[-1.0,-0.8,-0.6,-0.4,-0.2,0.,0.2,0.4,0.6,0.8,1.0])
	__YBINS = array('d',[0.,0.15,0.45,1.0])
	__ZBINS = array('d',[750.,1000.,1250.,1500.,1750.,2000.,2500.,3000.,4000.])

	def __init__(self,name,formatted_name,modifier) :
		#A Template has a name and a formatted name
		print '	Adding template with name '+name
		self.__name = name
		self.__type = 'nominal'
		if len(name.split('__'))==4 :
			self.__type = name.split('__')[2]+'__'+name.split('__')[3]
		self.__formatted_name = formatted_name
		#Templates have modifiers associated with them
		self.__modifier = modifier
		#Dictionary of weightsums for function replacements
		self.__weightsum_dict = {'NQQBAR':{'wname':None,'sig':0.},
								 'NQ1':{'wname':'wqs1','sig':0.},
								 'NQ2':{'wname':'wqs2','sig':0.},
								 'NGG':{'wname':None,'sig':0.},
								 'NG1':{'wname':'wg1','sig':0.},
								 'NG2':{'wname':'wg2','sig':0.},
								 'NG3':{'wname':'wg3','sig':0.},
								 'NG4':{'wname':'wg4','sig':0.},
								 'NBCK':{'wname':None,'sig':0.}}
		#Templates will eventually need a conversion function
		self.__conversion_function = None
		#Templates have 3D and 1D projection histograms
		self.__histo_3D = TH3D(name,	   formatted_name+'; c*; |x_{F}|; M (GeV)',len(self.__XBINS)-1,self.__XBINS,len(self.__YBINS)-1,self.__YBINS,len(self.__ZBINS)-1,self.__ZBINS)
		self.__histo_x  = TH1D(name+'_x',formatted_name+' X Projection; c*',len(self.__XBINS)-1,self.__XBINS)
		self.__histo_y  = TH1D(name+'_y',formatted_name+' Y Projection; |x_{F}|',len(self.__YBINS)-1,self.__YBINS)
		self.__histo_z  = TH1D(name+'_z',formatted_name+' Z Projection; M (GeV)',len(self.__ZBINS)-1,self.__ZBINS)
		#Set the directories of the newly created histograms
		self.__histo_3D.SetDirectory(0); self.__histo_x.SetDirectory(0); self.__histo_y.SetDirectory(0); self.__histo_z.SetDirectory(0)

	#get the event numbers for this process
	def buildWeightsums(self,ttree_dict,branch_dict,const_rw_name_list,ss_rw_name_list,JEC_append,channelcharge) :
		#Make a list of weightsum names that we'll be building
		weightsum_names = []
		for name in self.__weightsum_dict :
			if self.__name.find('fqq')!=-1 :
				if name.find('NQ')!=-1 :
					weightsum_names.append(name)
			elif self.__name.find('fgg')!=-1 :
				if name.find('NG')!=-1 :
					weightsum_names.append(name)
			elif self.__name.find('fbck')!=-1 :
				if name.find('NBCK')!=-1 :
					weightsum_names.append(name)
			elif self.__name.find('fntmj')!=-1 :
				if name.find('NNTMJ')!=-1 :
					weightsum_names.append(name)
		#Make a list of ttree identifiers
		ttree_identifiers = []
		for name in self.__weightsum_dict[weightsum_names[0]] :
			if name!='wname' :
				ttree_identifiers.append(name)
		#open each tree
		for identifier in ttree_identifiers :
			realidentifier = identifier+JEC_append
			tree = ttree_dict[realidentifier]
			#set branch addresses
			for branch in branch_dict.values() :
				tree.SetBranchAddress(branch.getPTreeName(),branch.getPTreeArray())
			#loop over entries
			nEntries = tree.GetEntries()
			for entry in range(nEntries) :
				tree.GetEntry(entry)
				#get event weight
				eweight = self.__get_event_weight__(branch_dict,const_rw_name_list,ss_rw_name_list)
				#for each weightsum we're making, apply the last remaining weight and then increment the value
				for weightsumname in weightsum_names :
					wname = self.__weightsum_dict[weightsumname]['wname']
					if channelcharge==0 :
						if wname==None :
							self.__weightsum_dict[weightsumname][identifier]+=2*eweight
						else :
							self.__weightsum_dict[weightsumname][identifier]+=eweight*branch_dict[wname].getPTreeValue()
							self.__weightsum_dict[weightsumname][identifier]+=eweight*branch_dict[wname+'_opp'].getPTreeValue()
					else :
						if wname==None :
							self.__weightsum_dict[weightsumname][identifier]+=eweight
						elif branch_dict['Q_l'].getPTreeValue()==channelcharge :
							self.__weightsum_dict[weightsumname][identifier]+=eweight*branch_dict[wname].getPTreeValue()
						elif branch_dict['Q_l'].getPTreeValue()!=channelcharge :
							self.__weightsum_dict[weightsumname][identifier]+=eweight*branch_dict[wname+'_opp'].getPTreeValue()

	#get the numbers of events in the regions specified by the keys of the event_numbers dictionary supplied
	def getEventNumbers(self,channelcharge,jecappend,ttree_dict,branch_dict,const_rw_name_list,ss_rw_list,event_numbers,function,fit_par_list,extra_weight) :
		#for each region (and tree)
		for ttree_identifier in event_numbers :
			#First replace the total event numbers in the function
			functionstring = self.__replace_function_string__(ttree_identifier,function,fit_par_list)
			#loop over the tree
			realidentifier = ttree_identifier
			if jecappend!=None :
				realidentifier+=jecappend
			tree = ttree_dict[realidentifier]
			for branch in branch_dict.values() :
				tree.SetBranchAddress(branch.getPTreeName(),branch.getPTreeArray())
			nEntries = tree.GetEntries()
			for entry in range(nEntries) :
				tree.GetEntry(entry)
				#get the event weight
				eweight = self.__get_event_weight__(branch_dict,const_rw_name_list,ss_rw_list)
				eweight_opp = eweight
				#multiply by the function weights if necessary
				if functionstring!=None and functionstring!='' :
					fweight, fweight_opp = self.__get_function_weight__(branch_dict,functionstring)
					eweight*=fweight; eweight_opp*=fweight_opp
				#multiply by the last +/-1 factor or whatever
				eweight*=extra_weight; eweight_opp*=extra_weight
				#add to the event numbers
				if channelcharge==0 or branch_dict['Q_l'].getPTreeValue()==channelcharge:
					event_numbers[ttree_identifier]+=eweight
				if branch_dict['addTwice'].getPTreeValue()==1 and (channelcharge==0 or branch_dict['Q_l'].getPTreeValue()!=channelcharge) :
					event_numbers[ttree_identifier]+=eweight_opp

	#add the given tree to the template
	def addTreeToTemplates(self,channelcharge,ttree_identifier,tree,branch_dict,const_rw_name_list,ss_rw_list,function,fit_par_list,extra_weight) :
		#replace the total event numbers in the function
		functionstring = self.__replace_function_string__(ttree_identifier,function,fit_par_list)
		#loop over the tree
		for branch in branch_dict.values() :
			tree.SetBranchAddress(branch.getPTreeName(),branch.getPTreeArray())
		nEntries = tree.GetEntries()
#		added = 0 #DEBUG
#		n=0 #DEBUG
		for entry in range(nEntries) :
#			n+=1 #DEBUG
			tree.GetEntry(entry)
			#get the event weight
			eweight = self.__get_event_weight__(branch_dict,const_rw_name_list,ss_rw_list)
#			s = 'event weight = '+str(eweight) #DEBUG
			#multiply by the function weights
			eweight_opp = eweight
			if functionstring!=None and functionstring!='' :
				fweight, fweight_opp = self.__get_function_weight__(branch_dict,functionstring)
				eweight*=fweight; eweight_opp*=fweight_opp
#				s+=', after func. weight = '+str(eweight) #DEBUG
			#multiply by the last +/-1 factor from the function call
			eweight*=extra_weight; eweight_opp*=extra_weight
#			s+=', after extra weight = '+str(eweight) #DEBUG
			#add to the templates
			x = branch_dict['cstar'].getPTreeValue()
			y = abs(branch_dict['x_F'].getPTreeValue())
			z = branch_dict['M'].getPTreeValue()
			if channelcharge==0 or branch_dict['Q_l'].getPTreeValue()==channelcharge:
				self.Fill(x,y,z,eweight)
#				added+=eweight #DEBUG
			if branch_dict['addTwice'].getPTreeValue()==1 and (channelcharge==0 or branch_dict['Q_l'].getPTreeValue()==(-1*channelcharge)) :
				self.Fill(-1.0*x,y,z,eweight_opp)
#				added+=eweight_opp #DEBUG
#			s+=' (weightsum = '+str(added)+', n = '+str(n)+')' #DEBUG
#			print s #DEBUG
#		print 'TEMPLATE '+self.__name+' INTEGRAL = '+str(self.__histo_3D.Integral())+' (total weight added = '+str(added)+')' #DEBUG

	#Fill the histograms given x, y, and z values
	def Fill(self,x,y,z,w) :
		inxbounds = x>=self.__XBINS[0] and x<self.__XBINS[len(self.__XBINS)-1]
		inybounds = y>=self.__YBINS[0] and y<self.__YBINS[len(self.__YBINS)-1]
		inzbounds = z>=self.__ZBINS[0] and z<self.__ZBINS[len(self.__ZBINS)-1]
		if inxbounds and inybounds and inzbounds :
			self.__histo_3D.Fill(x,y,z,w)
			self.__histo_x.Fill(x,w)
			self.__histo_y.Fill(y,w)
			self.__histo_z.Fill(z,w)

	#convertTo1D takes a 3D distribution and makes it 1D for use with theta
	def convertTo1D(self) :
		nBins = self.__histo_3D.GetNbinsX()*self.__histo_3D.GetNbinsY()*self.__histo_3D.GetNbinsZ()
		newHisto = TH1D(self.__histo_3D.GetName(),self.__histo_3D.GetTitle(),nBins,0.,nBins-1.)
		newHisto.SetDirectory(0)
		realbincounter = 1
		nglobalbins = self.__histo_3D.GetSize()
		for k in range(nglobalbins) :
			if not self.__histo_3D.IsBinOverflow(k) and not self.__histo_3D.IsBinUnderflow(k) :
				if not self.__histo_3D.GetBinContent(k) < 0. :
					newHisto.SetBinContent(realbincounter,self.__histo_3D.GetBinContent(k))
					if self.__histo_3D.GetBinContent(k)-self.__histo_3D.GetBinError(k)<0. :
						newHisto.SetBinError(realbincounter,self.__histo_3D.GetBinContent(k))
					else :
						newHisto.SetBinError(realbincounter,self.__histo_3D.GetBinError(k))
				realbincounter+=1
		return newHisto

	#make_from_1D_histo takes a 1D distribution and makes a template out of it!
	def make_from_1D_histo(self,histo_1D) :
		nglobalbins = self.__histo_3D.GetSize()
		global1Dbincounter = 1
		for k in range(nglobalbins) :
			if not self.__histo_3D.IsBinOverflow(k) and not self.__histo_3D.IsBinUnderflow(k) :
				content = histo_1D.GetBinContent(global1Dbincounter)
				error   = histo_1D.GetBinError(global1Dbincounter)
				self.__histo_3D.SetBinContent(k,content)
				self.__histo_3D.SetBinError(k,error)
				binx = array('i',[0]); biny = array('i',[0]); binz = array('i',[0])
				self.__histo_3D.GetBinXYZ(k,binx,biny,binz)
				self.__histo_x.SetBinContent(binx[0],self.__histo_x.GetBinContent(binx[0])+content)
				self.__histo_y.SetBinContent(biny[0],self.__histo_y.GetBinContent(biny[0])+content)
				self.__histo_z.SetBinContent(binz[0],self.__histo_z.GetBinContent(binz[0])+content)
				self.__histo_x.SetBinError(binx[0],self.__histo_x.GetBinError(binx[0])+error*error)
				self.__histo_y.SetBinError(biny[0],self.__histo_y.GetBinError(biny[0])+error*error)
				self.__histo_z.SetBinError(binz[0],self.__histo_z.GetBinError(binz[0])+error*error)
				global1Dbincounter+=1
		hs = [self.__histo_x,self.__histo_y,self.__histo_z]
		for h in hs :
			for k in range(h.GetSize()) :
				if not h.IsBinUnderflow(k) and not h.IsBinOverflow(k) :
					h.SetBinError(k,sqrt(h.GetBinError(k)))

	#Getters/Setters/Adders
	def getName(self) :
		return self.__name
	def getModifier(self) :
		return self.__modifier
	def getWeightsumDict(self) :
		return self.__weightsum_dict
	def setWeightsumDict(self,newdict) :
		self.__weightsum_dict = newdict
	def getType(self) :
		return self.__type
	def setConversionFunction(self,func) :
		self.__conversion_function = func
	def getHisto3D(self) :
		return self.__histo_3D
	def getHistoX(self) :
		return self.__histo_x
	def getHistoY(self) :
		return self.__histo_y
	def getHistoZ(self) :
		return self.__histo_z
	#Private methods
	#get the weight for this event based on reweights
	def __get_event_weight__(self,branch_dict,const_rw_name_list,ss_rw_list) :
		eweight = 1.0
		mod = self.__modifier
		#constant reweights
		#s = '' #DEBUG
		if const_rw_name_list!=None :
			for constrw in const_rw_name_list :
				eweight*=branch_dict[constrw].getPTreeValue()
		#		s+= 'after const rw'+constrw+' = '+str(eweight)+', ' #DEBUG
		#simple systematics
		if ss_rw_list!=None :
			ss_weights = [1.0,1.0] #BtoF first, then GH
			for ss in ss_rw_list :
				ssname = ss.getName()
				thisValueBtoF = 1.0; thisValueGH = 1.0
				if mod!=None and mod.isSSModifier() and mod.getName()==ssname :
					if self.__name.find('__up')!=-1 :
						thisValueBtoF=branch_dict[ssname+'_BtoF_up'].getPTreeValue()
						thisValueGH=branch_dict[ssname+'_GH_up'].getPTreeValue()
					elif self.__name.find('__down')!=-1 :
						thisValueBtoF=branch_dict[ssname+'_BtoF_down'].getPTreeValue()
						thisValueGH=branch_dict[ssname+'_GH_down'].getPTreeValue()
				else :
					if ss.isSplit() :
						thisValueBtoF = branch_dict[ssname+'_BtoF'].getPTreeValue()
						thisValueGH   = branch_dict[ssname+'_GH'].getPTreeValue()
					else :
						thisValueBtoF = thisValueGH = branch_dict[ssname].getPTreeValue()
				ss_weights[0]*=thisValueBtoF; ss_weights[1]*=thisValueGH
			eweight*=ss_weights[0]+ss_weights[1]
		#		s+='after systematics rw = '+str(eweight)+', ' #DEBUG
		#half the weight if we're adding it twice
		if branch_dict['addTwice'].getPTreeValue() == 1 :
			eweight*=0.5
		#s+='final = '+str(eweight) #DEBUG
		#print s #DEBUG
		return eweight
	#return a new function string with the total event numbers and fit parameters replaced as necessary
	def __replace_function_string__(self,ttree_identifier,function,fit_par_list) :
		#print 'OLD FUNCTION FOR TTREE IDENTIFIER '+ttree_identifier+' AND TEMPLATE '+self.__name+' = '+function.replace(' ','')
		newfunction1 = ''
		fssen = function.split('@')
		ngg = self.__weightsum_dict['NGG'][ttree_identifier]
		nqq = self.__weightsum_dict['NQQBAR'][ttree_identifier]
		nbck = self.__weightsum_dict['NBCK'][ttree_identifier]
		others_list = ['NG1','NG2','NG3','NG4','NQ1','NQ2']
		for s in fssen :
			if s == 'NTOT' :
				newfunction1+='('+str(ngg+nqq+nbck+nntmj)+')'
			elif s == 'NBCK' :
				newfunction1+='('+str(nbck)+')'
			elif s == 'NTTBAR' :
				newfunction1+='('+str(ngg+nqq)+')'
			elif s == 'NQQBAR' :
				newfunction1+='('+str(nqq)+')'
			elif s == 'NGG' :
				newfunction1+='('+str(ngg)+')'
			elif s in others_list :
				newfunction1+='('+str(self.__weightsum_dict[s][ttree_identifier])+')'
			else :
				newfunction1+=s
		newfunction2 = ''
		fssfp = newfunction1.split('#')
		for s in fssfp :
			replaced = False
			for fitpar in fit_par_list :
				if fitpar.getName() == s :
					if self.__type == 'par_'+fitpar.getName()+'__up' :
						newfunction2+='('+str(fitpar.getUpValue())+')'
					elif self.__type == 'par_'+fitpar.getName()+'__down' :
						newfunction2+='('+str(fitpar.getDownValue())+')'
					else :
						newfunction2+='('+str(fitpar.getNomValue())+')'
					replaced = True
					break
			if not replaced :
				newfunction2+=s
		#print 'NEW FUNCTION FOR TTREE IDENTIFIER '+ttree_identifier+' AND TEMPLATE '+self.__name+' = '+newfunction2.replace(' ','')
		return newfunction2
	#return the weight and opposite weight after applying the function reweights
	def __get_function_weight__(self,branch_dict,functionstring) :
		#replace the individual weights in the function string
		newfunction = ''
		newfunction_opp = ''
		fsser = functionstring.split('$')
		for s in fsser :
			if s in branch_dict.keys() :
				newfunction+='('+str(branch_dict[s].getPTreeValue())+')'
				newfunction_opp+='('+str(branch_dict[s+'_opp'].getPTreeValue())+')'
			else :
				newfunction+=s
				newfunction_opp+=s
		#print '	NEW FUNCTION WITH WEIGHTS = '+str(eval(newfunction))+' = '+newfunction.replace(' ','')
		return eval(newfunction), eval(newfunction_opp)

	def __del__(self) :
		pass

	def __str__(self) :
		s = self.__name
		return s