from ROOT import TH3D, TH1D, TFile
from array import array
import copy
from math import *

#binning dictionary: first key = channel name, second key = "x","y", or "z"
BINS_FINE = {
#type-3
't3_muplus_SR':{'x':array('d',[-1.000000,-0.900000,-0.800000,-0.700000,-0.600000,-0.500000,-0.400000,-0.300000,-0.200000,-0.100000,0.000000,0.100000,0.200000,0.300000,0.400000,0.500000,0.600000,0.700000,0.800000,0.900000,1.000000]),
'y':array('d',[0.000000,0.026134,0.046960,0.062639,0.078182,0.096256,0.114598,0.134829,0.154853,0.180000,0.600000]),
'z':array('d',[300.000000,378.889343,408.826233,434.843140,459.773987,485.331299,513.404663,547.011841,592.698303,671.100586,2000.000000])},
't3_muminus_SR':{'x':array('d',[-1.000000,-0.900000,-0.800000,-0.700000,-0.600000,-0.500000,-0.400000,-0.300000,-0.200000,-0.100000,0.000000,0.100000,0.200000,0.300000,0.400000,0.500000,0.600000,0.700000,0.800000,0.900000,1.000000]),
'y':array('d',[0.000000,0.022668,0.039310,0.054528,0.069643,0.083554,0.103297,0.126237,0.149162,0.180000,0.600000]),
'z':array('d',[300.000000,379.389374,408.797333,434.966095,459.401093,484.967529,513.116272,546.592285,591.992432,670.087830,2000.000000])},
't3_elplus_SR':{'x':array('d',[-1.000000,-0.900000,-0.800000,-0.700000,-0.600000,-0.500000,-0.400000,-0.300000,-0.200000,-0.100000,0.000000,0.100000,0.200000,0.300000,0.400000,0.500000,0.600000,0.700000,0.800000,0.900000,1.000000]),
'y':array('d',[0.000000,0.029700,0.048474,0.065726,0.080974,0.095914,0.117348,0.138890,0.162571,0.180000,0.600000]),
'z':array('d',[300.000000,381.137756,411.353455,437.247009,462.342651,487.694885,516.518982,550.481812,596.201660,675.163330,2000.000000])},
't3_elminus_SR':{'x':array('d',[-1.000000,-0.900000,-0.800000,-0.700000,-0.600000,-0.500000,-0.400000,-0.300000,-0.200000,-0.100000,0.000000,0.100000,0.200000,0.300000,0.400000,0.500000,0.600000,0.700000,0.800000,0.900000,1.000000]),
'y':array('d',[0.000000,0.017768,0.036985,0.053670,0.070198,0.086368,0.104584,0.128377,0.151648,0.180000,0.600000]),
'z':array('d',[300.000000,380.972839,411.131287,437.039062,462.022888,487.937286,516.535889,550.934631,596.562866,674.151123,2000.000000])},
#type-2 signal region
't2_muplus_SR':{'x':array('d',[-1.000000,-0.600000,-0.400000,-0.200000,0.000000,0.200000,0.400000,0.600000,1.000000]),
'y':array('d',[0.000000,0.069440,0.120388,0.154646,0.200000,0.600000]),
'z':array('d',[300.000000,563.186035,658.825012,753.855286,887.899414,3000.000000])},
't2_muminus_SR':{'x':array('d',[-1.000000,-0.600000,-0.400000,-0.200000,0.000000,0.200000,0.400000,0.600000,1.000000]),
'y':array('d',[0.000000,0.068279,0.121150,0.166912,0.200000,0.600000]),
'z':array('d',[300.000000,566.299805,659.969666,753.882568,882.469238,3000.000000])},
't2_elplus_SR':{'x':array('d',[-1.000000,-0.600000,-0.400000,-0.200000,0.000000,0.200000,0.400000,0.600000,1.000000]),
'y':array('d',[0.000000,0.047222,0.114373,0.160741,0.200000,0.600000]),
'z':array('d',[300.000000,592.787537,754.956665,885.332397,1057.502441,3000.000000])},
't2_elminus_SR':{'x':array('d',[-1.000000,-0.600000,-0.400000,-0.200000,0.000000,0.200000,0.400000,0.600000,1.000000]),
'y':array('d',[0.000000,0.030225,0.084181,0.153692,0.200000,0.600000]),
'z':array('d',[300.000000,591.313721,767.121765,903.888550,1068.874512,3000.000000])},
#type-2 control region
't2_muplus_WJets_CR':{'x':array('d',[-1.000000,-0.600000,-0.400000,-0.200000,0.000000,0.200000,0.400000,0.600000,1.000000]),
'y':array('d',[0.000000,0.067848,0.107569,0.150669,0.200000,0.600000]),
'z':array('d',[300.000000,607.664185,690.382812,788.903931,947.918518,3000.000000])},
't2_muminus_WJets_CR':{'x':array('d',[-1.000000,-0.600000,-0.400000,-0.200000,0.000000,0.200000,0.400000,0.600000,1.000000]),
'y':array('d',[0.000000,0.075112,0.121483,0.160047,0.200000,0.600000]),
'z':array('d',[300.000000,604.437317,689.684387,787.984741,945.892212,3000.000000])},
't2_elplus_WJets_CR':{'x':array('d',[-1.000000,-0.600000,-0.400000,-0.200000,0.000000,0.200000,0.400000,0.600000,1.000000]),
'y':array('d',[0.000000,0.096725,0.129873,0.161719,0.200000,0.600000]),
'z':array('d',[300.000000,692.641968,818.512756,957.834045,1210.355957,3000.000000])},
't2_elminus_WJets_CR':{'x':array('d',[-1.000000,-0.600000,-0.400000,-0.200000,0.000000,0.200000,0.400000,0.600000,1.000000]),
'y':array('d',[0.000000,0.085480,0.135928,0.166615,0.200000,0.600000]),
'z':array('d',[300.000000,685.536072,807.605530,956.955139,1212.451660,3000.000000])},
#type-1 signal region
't1_muplus_SR':{'x':array('d',[-1.000000,-0.500000,-0.250000,0.000000,0.250000,0.500000,1.000000]),
'y':array('d',[0.000000,0.132467,0.300000,0.600000]),
'z':array('d',[500.000000,1108.365967,3000.000000])},
't1_muminus_SR':{'x':array('d',[-1.000000,-0.500000,-0.250000,0.000000,0.250000,0.500000,1.000000]),
'y':array('d',[0.000000,0.150009,0.300000,0.600000]),
'z':array('d',[500.000000,1106.829102,3000.000000])},
't1_elplus_SR':{'x':array('d',[-1.000000,-0.500000,-0.250000,0.000000,0.250000,0.500000,1.000000]),
'y':array('d',[0.000000,0.161069,0.300000,0.600000]),
'z':array('d',[500.000000,1176.718018,3000.000000])},
't1_elminus_SR':{'x':array('d',[-1.000000,-0.500000,-0.250000,0.000000,0.250000,0.500000,1.000000]),
'y':array('d',[0.000000,0.162364,0.300000,0.600000]),
'z':array('d',[500.000000,1180.208252,3000.000000])},
#type-1 control region
't1_muplus_WJets_CR':{'x':array('d',[-1.000000,-0.500000,-0.250000,0.000000,0.250000,0.500000,1.000000]),
'y':array('d',[0.000000,0.139960,0.300000,0.600000]),
'z':array('d',[500.000000,941.432434,3000.000000])},
't1_muminus_WJets_CR':{'x':array('d',[-1.000000,-0.500000,-0.250000,0.000000,0.250000,0.500000,1.000000]),
'y':array('d',[0.000000,0.121112,0.300000,0.600000]),
'z':array('d',[500.000000,940.678223,3000.000000])},
't1_elplus_WJets_CR':{'x':array('d',[-1.000000,-0.500000,-0.250000,0.000000,0.250000,0.500000,1.000000]),
'y':array('d',[0.000000,0.141358,0.300000,0.600000]),
'z':array('d',[500.000000,1074.641968,3000.000000])},
't1_elminus_WJets_CR':{'x':array('d',[-1.000000,-0.500000,-0.250000,0.000000,0.250000,0.500000,1.000000]),
'y':array('d',[0.000000,0.145268,0.300000,0.600000]),
'z':array('d',[500.000000,1064.129150,3000.000000])},
}

BINS_COARSE = {}

#Template class
class Template(object) :

	def __init__(self,name,formatted_name,modifier,binning='fine') :
		#A Template has a name and a formatted name
		#print '	Adding template with name '+name
		self.__name = name
		self.__type = 'nominal'
		if len(name.split('__'))==3 :
			self.__type = name.split('__')[2]
		self.__formatted_name = formatted_name
		#Templates have modifiers associated with them
		self.__modifier = modifier
		#Dictionary of weightsums for function replacements
		self.__weightsum_dict = {'NQQBAR':{'wname':None,'sig':0.,'qcd_a':0.,'qcd_b':0.,'qcd_c':0.},
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
		#set the binning arrays based on the channel name
		if binning=='fine' :
			self._XBINS = copy.copy(BINS_FINE[self.__name.split('__')[0]]['x'])
			self._YBINS = copy.copy(BINS_FINE[self.__name.split('__')[0]]['y'])
			self._ZBINS = copy.copy(BINS_FINE[self.__name.split('__')[0]]['z'])
		elif binning=='coarse' :
			self._XBINS = copy.copy(BINS_COARSE[self.__name.split('__')[0]]['x'])
			self._YBINS = copy.copy(BINS_COARSE[self.__name.split('__')[0]]['y'])
			self._ZBINS = copy.copy(BINS_COARSE[self.__name.split('__')[0]]['z'])
		else :
			print 'ERROR: cannot determine which binning to use based on binning option %s'%(binning)
			return
		#QCD templates will eventually need a conversion factor
		self.__conversion_factor = None
		#Templates have 3D and 1D projection histograms
		self.__histo_3D = TH3D(name,formatted_name+'; c*; |x_{F}|; M (GeV)',len(self._XBINS)-1,self._XBINS,len(self._YBINS)-1,self._YBINS,len(self._ZBINS)-1,self._ZBINS)
		self.__histo_x  = TH1D(name+'_x',formatted_name+' X Projection; c*',len(self._XBINS)-1,self._XBINS)
		self.__histo_y  = TH1D(name+'_y',formatted_name+' Y Projection; |x_{F}|',len(self._YBINS)-1,self._YBINS)
		self.__histo_z  = TH1D(name+'_z',formatted_name+' Z Projection; M (GeV)',len(self._ZBINS)-1,self._ZBINS)
		#Set the directories of the newly created histograms
		self.__histo_3D.SetDirectory(0); self.__histo_x.SetDirectory(0); self.__histo_y.SetDirectory(0); self.__histo_z.SetDirectory(0)

	#get the event numbers for this process
	def buildWeightsums(self,ttree_dict,branch_dict,const_rw_name_list,ss_rw_name_list,JEC_append,channelcharge,ptfn) :
		#Make a list of weightsum names that we'll be building
		weightsum_names = []
		for name in self.__weightsum_dict :
			if self.__name.find('fqq')!=-1 :
				if name.find('NQ')!=-1 and name!='NQCD' :
					weightsum_names.append(name)
			elif self.__name.find('fgg')!=-1 :
				if name.find('NG')!=-1 :
					weightsum_names.append(name)
			elif self.__name.find('fbck')!=-1 :
				if name.find('NBCK')!=-1 :
					weightsum_names.append(name)
			elif self.__name.find('fwjets')!=-1 :
				if name.find('NWJETS')!=-1 :
					weightsum_names.append(name)
			elif self.__name.find('fqcd')!=-1 :
				if name.find('NQCD')!=-1 :
					weightsum_names.append(name)
		#Make a list of ttree identifiers
		ttree_identifiers = []
		for name in self.__weightsum_dict[weightsum_names[0]] :
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
						elif branch_dict['lep_Q'].getPTreeValue()==channelcharge :
							self.__weightsum_dict[weightsumname][identifier]+=eweight*branch_dict[wname].getPTreeValue()
						elif branch_dict['lep_Q'].getPTreeValue()!=channelcharge :
							self.__weightsum_dict[weightsumname][identifier]+=eweight*branch_dict[wname+'_opp'].getPTreeValue()
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
			if self.__conversion_factor!=None :
				eweight*=self.__conversion_factor; eweight_opp*=self.__conversion_factor
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
#		print 'TEMPLATE '+self.__name+' INTEGRAL = '+str(self.__histo_3D.Integral())+' (total weight added = '+str(added)+')' #DEBUG

	def fixNQCDValues(self) :
		value = self.__histo_3D.Integral()
		print '		Signal region NQCD for template '+self.__name+' = '+str(value)
		return value
	def setNQCDValue(self,value) :
		self.__weightsum_dict['NQCD']['sig']=value

	#Fill the histograms given x, y, and z values
	def Fill(self,x,y,z,w) :
		inxbounds = x>=self._XBINS[0] and x<self._XBINS[len(self._XBINS)-1]
		inybounds = y>=self._YBINS[0] and y<self._YBINS[len(self._YBINS)-1]
		inzbounds = z>=self._ZBINS[0] and z<self._ZBINS[len(self._ZBINS)-1]
		if inxbounds and inybounds and inzbounds :
			self.__histo_3D.Fill(x,y,z,w)
			self.__histo_x.Fill(x,w)
			self.__histo_y.Fill(y,w)
			self.__histo_z.Fill(z,w)

	#convertTo1D takes a 3D distribution and makes it 1D for use with combine
	def convertTo1D(self,bins_to_zero=None) :
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
		#correct any zero bins not in data
		zeroed_bins = []
		if self.__name.find('data_obs')==-1 :
			#fill zero or less bins with a very small value
			fillervalue = 0.00010
			if self.__type.endswith('Up') : 
				fillervalue = 0.0001005
			elif self.__type.endswith('Down') :
				fillervalue=0.0000995
			realbincounter = 1
			for k in range(nglobalbins) :
				if not self.__histo_3D.IsBinOverflow(k) and not self.__histo_3D.IsBinUnderflow(k) :
					if (bins_to_zero!=None and realbincounter in bins_to_zero) and not self.__modifier.getName()=='top_pt_re_weight' :
						newHisto.SetBinContent(realbincounter,fillervalue)
						zeroed_bins.append(realbincounter)
					elif not self.__histo_3D.GetBinContent(k) > 0. :
						newHisto.SetBinContent(realbincounter,0.00010)
						zeroed_bins.append(realbincounter)
					realbincounter+=1
		return newHisto,zeroed_bins

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
	def setConversionFactor(self,fac) :
		self.__conversion_factor = fac
	def getHisto3D(self) :
		return self.__histo_3D
	def getHistoX(self) :
		return self.__histo_x
	def getHistoY(self) :
		return self.__histo_y
	def getHistoZ(self) :
		return self.__histo_z
	def getHistos(self) :
		return [self.__histo_3D,self.__histo_x,self.__histo_y,self.__histo_z]
	def setHistos(self,hlist) :
		self.__histo_3D = hlist[0]
		self.__histo_x  = hlist[1]
		self.__histo_y  = hlist[2]
		self.__histo_z  = hlist[3]
	#Private methods
	#get the weight for this event based on reweights
	def __get_event_weight__(self,branch_dict,const_rw_name_list,ss_rw_list,skipiso=False,print_ind_weights=False) :
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
			#if isinstance(ss,str) : #DEBUG
			#	print 'ss = %s'%(ss) #DEBUG
				ssname = ss.getName()
				if ssname.find('_iso_weight')!=-1 and skipiso :
					continue
				thisValueBtoF = 1.0; thisValueGH = 1.0
				if mod!=None and mod.isSSModifier() and mod.getName()==ssname :
					if self.__name.endswith('Up') :
						if ss.isSplit() :
							thisValueBtoF=branch_dict[ssname+'_BtoF_up'].getPTreeValue()
							thisValueGH=branch_dict[ssname+'_GH_up'].getPTreeValue()
						else :
							thisValueBtoF = thisValueGH = branch_dict[ssname+'_up'].getPTreeValue()
					elif self.__name.endswith('Down') :
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
		#print 'OLD FUNCTION FOR TTREE IDENTIFIER '+ttree_identifier+' AND TEMPLATE '+self.__name+' = '+function.replace(' ','')
		newfunction1 = ''
		fssen = function.split('@')
		ngg = self.__weightsum_dict['NGG'][ttree_identifier]
		nqq = self.__weightsum_dict['NQQBAR'][ttree_identifier]
		nbck = self.__weightsum_dict['NBCK'][ttree_identifier]
		nwjets = self.__weightsum_dict['NWJETS'][ttree_identifier]
		nqcd = self.__weightsum_dict['NQCD'][ttree_identifier]
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
				newfunction1+='('+str(self.__weightsum_dict[s][ttree_identifier])+')'
			else :
				newfunction1+=s
		newfunction2 = ''
		fssfp = newfunction1.split('#')
		for s in fssfp :
			replaced = False
			for fitpar in fit_par_list :
				if fitpar.getName() == s :
					if self.__type == 'par_'+fitpar.getName()+'Up' :
						newfunction2+='('+str(fitpar.getUpValue())+')'
					elif self.__type == 'par_'+fitpar.getName()+'Down' :
						newfunction2+='('+str(fitpar.getDownValue())+')'
					else :
						newfunction2+='('+str(fitpar.getNomValue())+')'
					replaced = True
					break
			if not replaced :
				newfunction2+=s
		#print 'NEW FUNCTION FOR TTREE IDENTIFIER '+ttree_identifier+' AND TEMPLATE '+self.__name+' = '+newfunction2.replace(' ','') #DEBUG
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
		s = self.__name
		return s