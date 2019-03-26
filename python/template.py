from ROOT import TH2D, TH1D, TFile
from array import array
import copy
import random
from math import *

#value with which to fill zeroed bins
ZEROED_BIN_CONTENT=0.000001
#bin arrays
#cstar
cstar_in_one = array('d',[-1.00,1.00])
cstar_in_two = array('d',[-1.00,0.00,1.00])
cstar_in_four = array('d',[-1.00,-0.50,0.00,0.50,1.00])
cstar_in_six = array('d',[-1.00,-0.50,-0.25,0.00,0.25,0.50,1.00])
cstar_in_eight = array('d',[-1.00,-0.60,-0.40,-0.20,0.00,0.20,0.40,0.60,1.00])
cstar_in_ten = array('d',[-1.00,-0.80,-0.60,-0.40,-0.20,0.00,0.20,0.40,0.60,0.80,1.00])
cstar_in_sixteen = array('d',[-1.00,-0.875,-0.75,-0.625,-0.50,-0.375,-0.25,-0.125,0.00,0.125,0.25,0.375,0.50,0.625,0.75,0.875,1.00])
cstar_in_twenty = array('d',[-1.00,-0.90,-0.80,-0.70,-0.60,-0.50,-0.40,-0.30,-0.20,-0.10,0.00,0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.00])
#M
t1_M_in_one = array('d',[350.,5550.])
t1_M_in_two = array('d',[350.,1100.,5550.])
t1_M_in_three = array('d',[350.,900.,1250.,5550.])
t1_M_in_six = array('d',[350.,700.,900.,1075.,1250.,1600.,5550.])
t2_M_in_one = array('d',[350.,4500.])
t2_M_in_two = array('d',[350.,730.,5650.])
t2_M_in_three = array('d',[350.,730.,860.,5650.])
t2_M_in_six = array('d',[350.,550.,650.,750.,875.,1050.,5650.])
t3_M_in_one = array('d',[350.,2650.])
t3_M_in_two = array('d',[350.,475.,2650.])
t3_M_in_three = array('d',[350.,425.,500.,2650.])
t3_M_in_four = array('d',[350.,420.,475.,550.,2650.])
t3_M_in_five = array('d',[350.,400.,450.,500.,600.,2650.])
#dictionary of decided binnings
BINS = {
'SR':{
	't1':{
		'mu':[{'ylo':0.00,'yhi':0.10,'xbins':cstar_in_six,'zbins':t1_M_in_one},
			  {'ylo':0.10,'yhi':0.24,'xbins':cstar_in_two,'zbins':t1_M_in_two},
			  {'ylo':0.24,'yhi':1.00,'xbins':cstar_in_two,'zbins':t1_M_in_one},
			 ],
		'el':[{'ylo':0.00,'yhi':0.10,'xbins':cstar_in_four,'zbins':t1_M_in_one},
			  {'ylo':0.10,'yhi':0.24,'xbins':cstar_in_two,'zbins':t1_M_in_one},
			  {'ylo':0.24,'yhi':1.00,'xbins':cstar_in_two,'zbins':t1_M_in_one},
			 ],
		},
	't2':{
		'mu':[{'ylo':0.00,'yhi':0.04,'xbins':cstar_in_six,'zbins':t2_M_in_three},
			  {'ylo':0.04,'yhi':0.09,'xbins':cstar_in_six,'zbins':t2_M_in_two},
			  {'ylo':0.09,'yhi':0.15,'xbins':cstar_in_six,'zbins':t2_M_in_one},
			  {'ylo':0.15,'yhi':0.24,'xbins':cstar_in_two,'zbins':t2_M_in_three},
			  {'ylo':0.24,'yhi':1.00,'xbins':cstar_in_two,'zbins':t2_M_in_two},
			],
		'el':[{'ylo':0.00,'yhi':0.15,'xbins':cstar_in_four,'zbins':t2_M_in_one},
			  {'ylo':0.15,'yhi':1.00,'xbins':cstar_in_two,'zbins':t2_M_in_one},
			],
		},
	't3':{
		'mu':[{'ylo':0.00,'yhi':0.02,'xbins':cstar_in_twenty,'zbins':t3_M_in_three},
			  {'ylo':0.02,'yhi':0.04,'xbins':cstar_in_ten,'zbins':t3_M_in_five},
			  {'ylo':0.04,'yhi':0.06,'xbins':cstar_in_sixteen,'zbins':t3_M_in_two},
			  {'ylo':0.06,'yhi':0.08,'xbins':cstar_in_sixteen,'zbins':t3_M_in_two},
			  {'ylo':0.08,'yhi':0.10,'xbins':cstar_in_eight,'zbins':t3_M_in_two},
			  {'ylo':0.10,'yhi':0.115,'xbins':cstar_in_four,'zbins':t3_M_in_two},
			  {'ylo':0.115,'yhi':0.13,'xbins':cstar_in_four,'zbins':t3_M_in_two},
			  {'ylo':0.13,'yhi':0.15,'xbins':cstar_in_six,'zbins':t3_M_in_one},
			  {'ylo':0.15,'yhi':0.18,'xbins':cstar_in_six,'zbins':t3_M_in_one},
			  {'ylo':0.18,'yhi':1.00,'xbins':cstar_in_four,'zbins':t3_M_in_one},
			],
		'el':[{'ylo':0.00,'yhi':0.025,'xbins':cstar_in_ten,'zbins':t3_M_in_five},
			  {'ylo':0.025,'yhi':0.05,'xbins':cstar_in_twenty,'zbins':t3_M_in_two},
			  {'ylo':0.05,'yhi':0.075,'xbins':cstar_in_sixteen,'zbins':t3_M_in_two},
			  {'ylo':0.075,'yhi':0.10,'xbins':cstar_in_six,'zbins':t3_M_in_four},
			  {'ylo':0.10,'yhi':0.125,'xbins':cstar_in_four,'zbins':t3_M_in_two},
			  {'ylo':0.125,'yhi':0.155,'xbins':cstar_in_six,'zbins':t3_M_in_one},
			  {'ylo':0.155,'yhi':1.00,'xbins':cstar_in_four,'zbins':t3_M_in_one},
			],
		},
	},
'WJets_CR':{
	't1':{
		'mu':[{'ylo':0.00,'yhi':1.00,'xbins':cstar_in_one,'zbins':t1_M_in_six}],
		'el':[{'ylo':0.00,'yhi':1.00,'xbins':cstar_in_one,'zbins':t1_M_in_three}],
		},
	't2':{
		'mu':[{'ylo':0.00,'yhi':1.00,'xbins':cstar_in_one,'zbins':t2_M_in_six}],
		'el':[{'ylo':0.00,'yhi':1.00,'xbins':cstar_in_one,'zbins':t2_M_in_one}],
		},	
	},
}

#Template class
class Template(object) :

	def __init__(self,name,formatted_name,modifier) :
		#A Template has a name and a formatted name
		#print '	Adding template with name '+name
		self._name = name
		namesplit = name.split('__')
		self._topology = namesplit[0].split('_')[0]
		self._leptype = namesplit[0].split('_')[1]
		self._region = namesplit[0][len(self._topology+'_'+self._leptype+'_'):]
		self._process = namesplit[1]
		self._type = 'nominal'
		if len(namesplit)==3 :
			self._type = namesplit[2]
		self._formatted_name = formatted_name
		#Templates have modifiers associated with them
		self._modifier = modifier
		#Dictionary of weightsums for function replacements
		self._weightsum_dict = {'NQQBAR':{'wname':None,'sig':0.,'qcd_a':0.,'qcd_b':0.,'qcd_c':0.},
								 'NQ1':{'wname':'wqs1','sig':0.,'qcd_a':0.,'qcd_b':0.,'qcd_c':0.},
								 'NQ2':{'wname':'wqs2','sig':0.,'qcd_a':0.,'qcd_b':0.,'qcd_c':0.},
								 'NGG':{'wname':None,'sig':0.,'qcd_a':0.,'qcd_b':0.,'qcd_c':0.},
								 'NG1':{'wname':'wg1','sig':0.,'qcd_a':0.,'qcd_b':0.,'qcd_c':0.},
								 'NG2':{'wname':'wg2','sig':0.,'qcd_a':0.,'qcd_b':0.,'qcd_c':0.},
								 'NG3':{'wname':'wg3','sig':0.,'qcd_a':0.,'qcd_b':0.,'qcd_c':0.},
								 'NG4':{'wname':'wg4','sig':0.,'qcd_a':0.,'qcd_b':0.,'qcd_c':0.},
								 'NBCK':{'wname':None,'sig':0.,'qcd_a':0.,'qcd_b':0.,'qcd_c':0.},
								 'NWJETS':{'wname':None,'sig':0.,'qcd_a':0.,'qcd_b':0.,'qcd_c':0.},
								 'NQCD':{'wname':None,'sig':0.,'qcd_a':0.,'qcd_b':0.,'qcd_c':0.}}
		#set the binning based on the channel name
		self._bins = BINS[self._region][self._topology][self._leptype[:2]]
		#QCD templates will eventually need a conversion factor
		self._conversion_factor = None
		#Templates are dictionaries of 2D (cstar,M) histograms keyed by tuples of |x_F| ranges
		self._histos = {}
		self._ntotalbins=0 #number of bins for the unrolled 1D histogram
		for d in self._bins :
			histoname = name+'_y='+str(d['ylo'])+'to'+str(d['yhi'])
			histoformattedname = formatted_name+' y from '+str(d['ylo'])+' to '+str(d['yhi'])
			newHisto=TH2D(histoname,histoformattedname,len(d['xbins'])-1,d['xbins'],len(d['zbins'])-1,d['zbins'])
			newHisto.SetDirectory(0)
			self._histos[(d['ylo'],d['yhi'])]=newHisto
			self._ntotalbins+=(len(d['xbins'])-1)*(len(d['zbins'])-1)

	#get the event numbers for this process
	def buildWeightsums(self,ttree_dict,branch_dict,const_rw_name_list,ss_rw_name_list,JEC_append,channelcharge,ptfn) :
		#Make a list of weightsum names that we'll be building
		weightsum_names = []
		for name in self._weightsum_dict :
			if self._process=='fqq' :
				if name.find('NQ')!=-1 and name!='NQCD' :
					weightsum_names.append(name)
			elif self._process=='fgg' :
				if name.find('NG')!=-1 :
					weightsum_names.append(name)
			elif self._process=='fbck' :
				if name.find('NBCK')!=-1 :
					weightsum_names.append(name)
			elif self._process=='fwjets' :
				if name.find('NWJETS')!=-1 :
					weightsum_names.append(name)
			elif self._process=='fqcd' :
				if name.find('NQCD')!=-1 :
					weightsum_names.append(name)
		#Make a list of ttree identifiers
		ttree_identifiers = []
		for name in self._weightsum_dict[weightsum_names[0]] :
			if name!='wname' :
				ttree_identifiers.append(name)
		#open each tree
		for identifier in ttree_identifiers :
			realidentifier = identifier+JEC_append
			filep = TFile(ptfn)
			tree = filep.Get(ttree_dict[realidentifier].GetName())
			#set branch addresses
			for branch in branch_dict.values() :
				tree.SetBranchAddress(branch.getPTreeName(),branch.getPTreeArray())
			#loop over entries
			nEntries = tree.GetEntries()
			for entry in range(nEntries) :
				tree.GetEntry(entry)
				#get event weight
				eweight = self.__get_event_weight__(branch_dict,const_rw_name_list,ss_rw_name_list,(identifier.find('qcd_b')!=-1 or identifier.find('qcd_c')!=-1))
				#for each weightsum we're making, apply the last remaining weight and then increment the value
				for weightsumname in weightsum_names :
					wname = self._weightsum_dict[weightsumname]['wname']
					if channelcharge==0 :
						if wname==None :
							self._weightsum_dict[weightsumname][identifier]+=2*eweight
						else :
							self._weightsum_dict[weightsumname][identifier]+=eweight*branch_dict[wname].getPTreeValue()
							self._weightsum_dict[weightsumname][identifier]+=eweight*branch_dict[wname+'_opp'].getPTreeValue()
					else :
						if wname==None :
							self._weightsum_dict[weightsumname][identifier]+=eweight
						elif branch_dict['lep_Q'].getPTreeValue()==channelcharge :
							self._weightsum_dict[weightsumname][identifier]+=eweight*branch_dict[wname].getPTreeValue()
						elif branch_dict['lep_Q'].getPTreeValue()!=channelcharge :
							self._weightsum_dict[weightsumname][identifier]+=eweight*branch_dict[wname+'_opp'].getPTreeValue()
			filep.Close()

	#get the numbers of events in the regions specified by the keys of the event_numbers dictionary supplied
	def getEventNumbers(self,channelcharge,jecappend,ttree_dict,branch_dict,const_rw_name_list,ss_rw_list,event_numbers,function,fit_par_list,extra_weight,ptfn) :
		#for each region (and tree)
		for ttree_identifier in event_numbers :
			#First replace the total event numbers in the function
			functionstring = self.__replace_function_string__(ttree_identifier,function,fit_par_list)
			#loop over the tree
			realidentifier = ttree_identifier
			if jecappend!=None :
				realidentifier+=jecappend
			filep = TFile(ptfn)
			tree = filep.Get(ttree_dict[realidentifier].GetName())
			for branch in branch_dict.values() :
				tree.SetBranchAddress(branch.getPTreeName(),branch.getPTreeArray())
			nEntries = tree.GetEntries()
			for entry in range(nEntries) :
				#if entry==1 : print 'doing entry %d of %d for template with name %s'%(entry,nEntries,self.getName()) #DEBUG
				tree.GetEntry(entry)
				#get the event weight
				eweight = self.__get_event_weight__(branch_dict,const_rw_name_list,ss_rw_list,(ttree_identifier.find('qcd_b')!=-1 or ttree_identifier.find('qcd_c')!=-1),entry==1)
				eweight_opp = eweight
				#if entry==1 : print '	initial = %.4f'%(eweight) #DEBUG
				#multiply by the function weights if necessary
				if functionstring!=None and functionstring!='' :
					fweight, fweight_opp = self.__get_function_weight__(branch_dict,functionstring)
					eweight*=fweight; eweight_opp*=fweight_opp
				#if entry==1 : print '	after function weight = %.4f'%(eweight) #DEBUG
				#multiply by the last +/-1 factor or whatever
				eweight*=extra_weight; eweight_opp*=extra_weight
				#if entry==1 : print '	after extra weight = %.4f'%(eweight) #DEBUG
				#add to the event numbers
				if channelcharge==0 or branch_dict['lep_Q'].getPTreeValue()==channelcharge:
					event_numbers[ttree_identifier]+=eweight
				if branch_dict['addTwice'].getPTreeValue()==1 and (channelcharge==0 or branch_dict['lep_Q'].getPTreeValue()!=channelcharge) :
					event_numbers[ttree_identifier]+=eweight_opp
				#if entry==1 : print 'for template %s with realidentifier %s eweight for entry %d = %.6f'%(self.getName(),realidentifier,entry,eweight) #DEBUG
			filep.Close()
		return event_numbers

	#add the given tree to the template
	def addTreeToTemplates(self,channelcharge,ttree_identifier,tree,branch_dict,const_rw_name_list,ss_rw_list,function,fit_par_list,extra_weight,ptfn,pb) :
		#replace the total event numbers in the function
		functionstring = self.__replace_function_string__(ttree_identifier,function,fit_par_list)
		filep = TFile(ptfn)
		tree = filep.Get(tree.GetName())
		#loop over the tree
		for branch in branch_dict.values() :
			tree.SetBranchAddress(branch.getPTreeName(),branch.getPTreeArray())
		nEntries = tree.GetEntries()
#		added = 0 #DEBUG
#		n=0 #DEBUG
		for entry in range(nEntries) :
			if pb :
				percentdone = 100.*entry/nEntries
				if percentdone%1.<100./nEntries and percentdone%10.<((100.*(entry-1)/nEntries)%10.) :
					print '		%d%% done'%(percentdone)
#			n+=1 #DEBUG
			tree.GetEntry(entry)
			#get the event weight
			eweight = self.__get_event_weight__(branch_dict,const_rw_name_list,ss_rw_list,(ttree_identifier.find('qcd_b')!=-1 or ttree_identifier.find('qcd_c')!=-1))
#			s = 'event weight = '+str(eweight) #DEBUG
			#multiply by the function weights
			eweight_opp = eweight
			if functionstring!=None and functionstring!='' :
				fweight, fweight_opp = self.__get_function_weight__(branch_dict,functionstring)
				eweight*=fweight; eweight_opp*=fweight_opp
#				s+=', after func. weight = '+str(eweight) #DEBUG
			#if it has a conversion factor, apply it.
			if self._conversion_factor!=None :
				eweight*=self._conversion_factor; eweight_opp*=self._conversion_factor
			#multiply by the last +/-1 factor from the function call
			eweight*=extra_weight; eweight_opp*=extra_weight
#			s+=', after extra weight = '+str(eweight) #DEBUG
			#add to the templates
			x = branch_dict['cstar'].getPTreeValue()
			y = abs(branch_dict['x_F'].getPTreeValue())
			z = branch_dict['M'].getPTreeValue()
			if channelcharge==0 or branch_dict['lep_Q'].getPTreeValue()==channelcharge:
				self.Fill(x,y,z,eweight)
#				added+=eweight #DEBUG
			if branch_dict['addTwice'].getPTreeValue()==1 and (channelcharge==0 or branch_dict['lep_Q'].getPTreeValue()==(-1*channelcharge)) :
				self.Fill(-1.0*x,y,z,eweight_opp)
#				added+=eweight_opp #DEBUG
#			s+=' (weightsum = '+str(added)+', n = '+str(n)+')' #DEBUG
#			print s #DEBUG
#		print 'TEMPLATE '+self._name+' INTEGRAL = '+str(self._histo_3D.Integral())+' (total weight added = '+str(added)+')' #DEBUG

	def fixNQCDValues(self) :
		value = 0
		for h in self._histos.values() :
			value+=h.Integral()
		print '		Signal region NQCD for template '+self._name+' = '+str(value)
		return value
	def setNQCDValue(self,value) :
		self._weightsum_dict['NQCD']['sig']=value

	#Fill the histograms given x, y, and z values
	def Fill(self,x,y,z,w) :
		for yt,h in self._histos.items() :
			if y>=yt[0] and y<yt[1] :
				h.Fill(x,z,w)
				break

	#convertTo1D takes the list of 2D histograms and makes it into one unrolled 1D for use with combine
	def convertTo1D(self) :
		newHisto = TH1D(self._name,self._formatted_name,self._ntotalbins,0.,self._ntotalbins-1.)
		newHisto.SetDirectory(0)
		realbincounter = 1
		sortedytlist = sorted(self._histos)
		for yt in sortedytlist :
			h = self._histos[yt]
			nglobalbins=h.GetSize()
			for k in range(nglobalbins) :
				if not h.IsBinOverflow(k) and not h.IsBinUnderflow(k) :
					cont=max(h.GetBinContent(k),ZEROED_BIN_CONTENT)
					err = cont if cont==ZEROED_BIN_CONTENT else h.GetBinError(k)
					newHisto.SetBinContent(realbincounter,cont)
					newHisto.SetBinError(realbincounter,err)
					#print('name = {}, realbincounter = {}, content = {}, error = {}'.format(self._name,realbincounter,cont,err))
					realbincounter+=1
		return newHisto

	#make_from_1D_histo takes an unrolled 1D distribution and makes a template out of it!
	def make_from_1D_histo(self,histo1D) :
		global1Dbincounter = 1
		sortedytlist = sorted(self._histos)
		for yt in sortedytlist :
			h=self._histos[yt]
			nglobalbins=h.GetSize()
			for k in range(nglobalbins) :
				if not h.IsBinOverflow(k) and not h.IsBinUnderflow(k) :
					cont=histo1D.GetBinContent(global1Dbincounter)
					err=histo1D.GetBinError(global1Dbincounter)
					h.SetBinContent(k,cont)
					h.SetBinError(k,err)
					global1Dbincounter+=1

	#Getters/Setters/Adders
	def getName(self) :
		return self._name
	def getModifier(self) :
		return self._modifier
	def getWeightsumDict(self) :
		return self._weightsum_dict
	def setWeightsumDict(self,newdict) :
		self._weightsum_dict = newdict
	def getType(self) :
		return self._type
	def setConversionFactor(self,fac) :
		self._conversion_factor = fac
	def getHistos(self) :
		return self._histos
	def setHistos(self,histos) :
		self._histos = histos
	def getYbins(self) :
		ybinslist = []
		sortedytlist = sorted(self._histos)
		for i in range(len(sortedytlist)) :
			if i==0 :
				ybinslist.append(sortedytlist[i][0])
			ybinslist.append(sortedytlist[i][1])
		return array('d',ybinslist)
	def getBins(self) :
		return self._bins

	#Private methods
	#get the weight for this event based on reweights
	def __get_event_weight__(self,branch_dict,const_rw_name_list,ss_rw_list,skipiso=False,print_ind_weights=False) :
		eweight = 1.0
		mod = self._modifier
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
		#		if isinstance(ss,str) : #DEBUG
		#			print 'ss = %s'%(ss) #DEBUG
				ssname = ss.getName()
				if ssname.find('_iso_weight')!=-1 and skipiso :
					continue
				thisValueBtoF = 1.0; thisValueGH = 1.0
				if mod!=None and mod.isSSModifier() and mod.getName()==ssname :
					if self._name.endswith('Up') :
						if ss.isSplit() :
							thisValueBtoF=branch_dict[ssname+'_BtoF_up'].getPTreeValue()
							thisValueGH=branch_dict[ssname+'_GH_up'].getPTreeValue()
						else :
							thisValueBtoF = thisValueGH = branch_dict[ssname+'_up'].getPTreeValue()
					elif self._name.endswith('Down') :
						if ss.isSplit() :
							thisValueBtoF=branch_dict[ssname+'_BtoF_down'].getPTreeValue()
							thisValueGH=branch_dict[ssname+'_GH_down'].getPTreeValue()
						else :
							thisValueBtoF = thisValueGH = branch_dict[ssname+'_down'].getPTreeValue()
				else :
					if ss.isSplit() :
						thisValueBtoF = branch_dict[ssname+'_BtoF'].getPTreeValue()
						thisValueGH   = branch_dict[ssname+'_GH'].getPTreeValue()
					else :
						thisValueBtoF = thisValueGH = branch_dict[ssname].getPTreeValue()
				ss_weights[0]*=thisValueBtoF; ss_weights[1]*=thisValueGH
		#		s+='after sys. %s rw: %s '%(ssname,ss_weights) #DEBUG
			eweight*=ss_weights[0]+ss_weights[1]
		#	s+='after systematics rw = '+str(eweight)+', ' #DEBUG
		#half the weight if we're adding it twice
		if branch_dict['addTwice'].getPTreeValue() == 1 :
			eweight*=0.5
		#s+='final = '+str(eweight) #DEBUG
		#if print_ind_weights : #DEBUG
		#	print s #DEBUG
		return eweight
	#return a new function string with the total event numbers and fit parameters replaced as necessary
	def __replace_function_string__(self,ttree_identifier,function,fit_par_list) :
		#print 'OLD FUNCTION FOR TTREE IDENTIFIER '+ttree_identifier+' AND TEMPLATE '+self._name+' = '+function.replace(' ','')
		newfunction1 = ''
		fssen = function.split('@')
		ngg = self._weightsum_dict['NGG'][ttree_identifier]
		nqq = self._weightsum_dict['NQQBAR'][ttree_identifier]
		nbck = self._weightsum_dict['NBCK'][ttree_identifier]
		nwjets = self._weightsum_dict['NWJETS'][ttree_identifier]
		nqcd = self._weightsum_dict['NQCD'][ttree_identifier]
		others_list = ['NG1','NG2','NG3','NG4','NQ1','NQ2']
		for s in fssen :
			if s == 'NTOT' :
				newfunction1+='('+str(ngg+nqq+nbck+nwjets+nqcd)+')'
			elif s == 'NBCK' :
				newfunction1+='('+str(nbck)+')'
			elif s == 'NWJETS' :
				newfunction1+='('+str(nwjets)+')'
			elif s == 'NTTBAR' :
				newfunction1+='('+str(ngg+nqq)+')'
			elif s == 'NQQBAR' :
				newfunction1+='('+str(nqq)+')'
			elif s == 'NGG' :
				newfunction1+='('+str(ngg)+')'
			elif s == 'NQCD' :
				newfunction1+='('+str(nqcd)+')'
			elif s in others_list :
				newfunction1+='('+str(self._weightsum_dict[s][ttree_identifier])+')'
			else :
				newfunction1+=s
		newfunction2 = ''
		fssfp = newfunction1.split('#')
		for s in fssfp :
			replaced = False
			for fitpar in fit_par_list :
				if fitpar.getName() == s :
					if self._type == 'par_'+fitpar.getName()+'Up' :
						newfunction2+='('+str(fitpar.getUpValue())+')'
					elif self._type == 'par_'+fitpar.getName()+'Down' :
						newfunction2+='('+str(fitpar.getDownValue())+')'
					else :
						newfunction2+='('+str(fitpar.getNomValue())+')'
					replaced = True
					break
			if not replaced :
				newfunction2+=s
		#print 'NEW FUNCTION FOR TTREE IDENTIFIER '+ttree_identifier+' AND TEMPLATE '+self._name+' = '+newfunction2.replace(' ','') #DEBUG
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
		#print '	NEW FUNCTION TO EVALUATE = %s'%(newfunction.replace(' ','')) #DEBUG
		#print '	NEW FUNCTION WITH WEIGHTS = '+str(eval(newfunction))+' = '+newfunction.replace(' ','') #DEBUG
		return eval(newfunction), eval(newfunction_opp)

	def __del__(self) :
		pass

	def __str__(self) :
		s = self._name
		return s