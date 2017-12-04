from ROOT import TH3D, TH1D, TFile
from array import array
import copy
from math import *

#binning dictionary: first key = channel name, second key = "x","y", or "z"
BINS = {'t3_muplus_SR':{'x':array('d',[-1.00000,-0.90000,-0.80000,-0.70000,-0.60000,-0.50000,-0.40000,-0.30000,-0.20000,-0.10000,0.00000,0.10000,0.20000,0.30000,0.40000,0.50000,0.60000,0.70000,0.80000,0.90000,1.00000]),
'y':array('d',[0.00000,0.01691,0.02982,0.03981,0.04977,0.06245,0.07368,0.08181,0.09279,0.10543,0.11968,0.13544,0.14612,0.16075,0.17759,0.19505,0.22127,0.25083,0.30000]),
'z':array('d',[300.00000,305.74890,314.63672,320.27063,327.75574,329.87878,332.00891,335.20145,338.67850,341.83798,346.91391,350.31372,352.16711,353.78729,356.11069,359.73636,361.44217,367.68130,371.54974,372.93878,374.93176,377.00763,379.59964,382.62524,1500.00000])},
't2_muplus_WJets_CR':{'x':array('d',[-1.00000,-0.85714,-0.71429,-0.57143,-0.42857,-0.28571,-0.14286,0.00000,0.14286,0.28571,0.42857,0.57143,0.71429,0.85714,1.00000]),
'y':array('d',[0.00000,0.07091,0.10823,0.14652,0.18219,0.21799,0.24799,0.31768,0.40000]),
'z':array('d',[300.00000,426.94275,475.51422,543.91333,741.89246,856.23193,987.26947,1144.23633,1413.81812,1784.56128,2300.00000])},
't3_elplus_SR':{'x':array('d',[-1.00000,-0.90000,-0.80000,-0.70000,-0.60000,-0.50000,-0.40000,-0.30000,-0.20000,-0.10000,0.00000,0.10000,0.20000,0.30000,0.40000,0.50000,0.60000,0.70000,0.80000,0.90000,1.00000]),
'y':array('d',[0.00000,0.01646,0.03289,0.04244,0.05527,0.06861,0.08054,0.09000,0.10235,0.11674,0.12633,0.14073,0.16257,0.18123,0.19298,0.21306,0.23045,0.26233,0.30000]),
'z':array('d',[300.00000,313.54385,325.86469,334.70886,340.63721,343.48889,349.80060,356.94043,362.16412,369.42419,378.51459,382.12711,388.96729,393.23828,397.69586,401.48309,405.23770,408.16177,416.81284,424.52151,429.00797,432.15646,440.15027,443.03619,1500.00000])},
't1_elplus_SR':{'x':array('d',[-1.00000,-0.75000,-0.50000,-0.25000,0.00000,0.25000,0.50000,0.75000,1.00000]),
'y':array('d',[0.00000,0.16104,0.35349,0.60000]),
'z':array('d',[500.00000,859.18774,1377.63208,1832.34106,3000.00000])},
't1_muplus_WJets_CR':{'x':array('d',[-1.00000,-0.75000,-0.50000,-0.25000,0.00000,0.25000,0.50000,0.75000,1.00000]),
'y':array('d',[0.00000,0.13706,0.30875,0.60000]),
'z':array('d',[500.00000,1167.25098,1489.60571,1958.93555,3000.00000])},
't2_elminus_WJets_CR':{'x':array('d',[-1.00000,-0.85714,-0.71429,-0.57143,-0.42857,-0.28571,-0.14286,0.00000,0.14286,0.28571,0.42857,0.57143,0.71429,0.85714,1.00000]),
'y':array('d',[0.00000,0.09084,0.13573,0.16373,0.19464,0.23349,0.30292,0.37648,0.40000]),
'z':array('d',[300.00000,418.84286,485.12701,547.82996,843.03516,996.79987,1172.23120,1416.18774,1642.80640,1965.79321,2300.00000])},
't3_elminus_SR':{'x':array('d',[-1.00000,-0.90000,-0.80000,-0.70000,-0.60000,-0.50000,-0.40000,-0.30000,-0.20000,-0.10000,0.00000,0.10000,0.20000,0.30000,0.40000,0.50000,0.60000,0.70000,0.80000,0.90000,1.00000]),
'y':array('d',[0.00000,0.00964,0.02309,0.03907,0.04869,0.05982,0.06931,0.08539,0.09711,0.11145,0.11921,0.13547,0.14895,0.16756,0.19033,0.21060,0.22591,0.25472,0.30000]),
'z':array('d',[300.00000,315.32001,331.91760,345.42426,363.11716,381.12296,400.23334,436.25262,507.77448,536.13550,565.18896,598.25305,623.06470,655.78296,677.03448,715.84436,764.75879,832.10455,893.42664,962.03351,1039.26746,1113.81006,1223.83203,1378.43091,1500.00000])},
't3_muminus_SR':{'x':array('d',[-1.00000,-0.90000,-0.80000,-0.70000,-0.60000,-0.50000,-0.40000,-0.30000,-0.20000,-0.10000,0.00000,0.10000,0.20000,0.30000,0.40000,0.50000,0.60000,0.70000,0.80000,0.90000,1.00000]),
'y':array('d',[0.00000,0.01502,0.02833,0.04705,0.05490,0.06722,0.08072,0.09239,0.10211,0.11414,0.12633,0.13711,0.15585,0.16870,0.18954,0.20861,0.23880,0.27297,0.30000]),
'z':array('d',[300.00000,310.29221,321.04169,329.12308,334.27814,343.80292,351.38385,362.77777,370.65735,379.72992,389.29608,401.16022,417.14508,427.45251,439.25708,475.72272,506.96964,525.69232,547.81824,563.84766,577.55090,590.02356,610.75031,626.92969,1500.00000])},
't1_elminus_SR':{'x':array('d',[-1.00000,-0.75000,-0.50000,-0.25000,0.00000,0.25000,0.50000,0.75000,1.00000]),
'y':array('d',[0.00000,0.17717,0.33967,0.60000]),
'z':array('d',[500.00000,865.94940,1378.33411,1685.97620,3000.00000])},
't2_muminus_WJets_CR':{'x':array('d',[-1.00000,-0.85714,-0.71429,-0.57143,-0.42857,-0.28571,-0.14286,0.00000,0.14286,0.28571,0.42857,0.57143,0.71429,0.85714,1.00000]),
'y':array('d',[0.00000,0.07373,0.10418,0.13765,0.18569,0.23076,0.28350,0.33609,0.40000]),
'z':array('d',[300.00000,404.05518,468.19135,531.49347,784.38855,940.58496,1061.29749,1233.77393,1456.16846,1783.42810,2300.00000])},
't2_elplus_SR':{'x':array('d',[-1.00000,-0.85714,-0.71429,-0.57143,-0.42857,-0.28571,-0.14286,0.00000,0.14286,0.28571,0.42857,0.57143,0.71429,0.85714,1.00000]),
'y':array('d',[0.00000,0.04453,0.07876,0.12355,0.16248,0.21025,0.25656,0.32892,0.40000]),
'z':array('d',[300.00000,357.47830,405.58008,653.82422,905.52307,1022.43335,1099.89819,1244.01453,1362.01538,1670.84204,2300.00000])},
't1_muplus_SR':{'x':array('d',[-1.00000,-0.75000,-0.50000,-0.25000,0.00000,0.25000,0.50000,0.75000,1.00000]),
'y':array('d',[0.00000,0.15641,0.30757,0.60000]),
'z':array('d',[500.00000,810.53833,1335.70203,1871.76111,3000.00000])},
't2_elminus_SR':{'x':array('d',[-1.00000,-0.85714,-0.71429,-0.57143,-0.42857,-0.28571,-0.14286,0.00000,0.14286,0.28571,0.42857,0.57143,0.71429,0.85714,1.00000]),
'y':array('d',[0.00000,0.04917,0.09219,0.13562,0.16830,0.22114,0.26963,0.31276,0.40000]),
'z':array('d',[300.00000,354.16663,380.25574,428.74789,639.57947,887.92340,1102.52686,1297.17114,1530.81116,1755.49316,2300.00000])},
't1_elminus_WJets_CR':{'x':array('d',[-1.00000,-0.75000,-0.50000,-0.25000,0.00000,0.25000,0.50000,0.75000,1.00000]),
'y':array('d',[0.00000,0.17992,0.37258,0.60000]),
'z':array('d',[500.00000,603.58478,1367.24341,2016.21497,3000.00000])},
't1_elplus_WJets_CR':{'x':array('d',[-1.00000,-0.75000,-0.50000,-0.25000,0.00000,0.25000,0.50000,0.75000,1.00000]),
'y':array('d',[0.00000,0.16925,0.35644,0.60000]),
'z':array('d',[500.00000,659.11121,1310.55078,1736.47241,3000.00000])},
't2_muminus_SR':{'x':array('d',[-1.00000,-0.85714,-0.71429,-0.57143,-0.42857,-0.28571,-0.14286,0.00000,0.14286,0.28571,0.42857,0.57143,0.71429,0.85714,1.00000]),
'y':array('d',[0.00000,0.03182,0.06724,0.10608,0.14920,0.19500,0.24296,0.29955,0.40000]),
'z':array('d',[300.00000,346.63495,392.24628,485.06415,760.74268,876.91016,957.26440,1090.60034,1315.81152,1687.08618,2300.00000])},
't1_muminus_SR':{'x':array('d',[-1.00000,-0.75000,-0.50000,-0.25000,0.00000,0.25000,0.50000,0.75000,1.00000]),
'y':array('d',[0.00000,0.15334,0.30814,0.60000]),
'z':array('d',[500.00000,755.27478,1371.70532,1821.45142,3000.00000])},
't2_muplus_SR':{'x':array('d',[-1.00000,-0.85714,-0.71429,-0.57143,-0.42857,-0.28571,-0.14286,0.00000,0.14286,0.28571,0.42857,0.57143,0.71429,0.85714,1.00000]),
'y':array('d',[0.00000,0.03598,0.06943,0.10161,0.14247,0.18224,0.21544,0.28274,0.40000]),
'z':array('d',[300.00000,335.67422,375.10474,450.65985,773.23352,859.74432,974.70831,1079.76880,1335.55078,1717.48975,2300.00000])},
't2_elplus_WJets_CR':{'x':array('d',[-1.00000,-0.85714,-0.71429,-0.57143,-0.42857,-0.28571,-0.14286,0.00000,0.14286,0.28571,0.42857,0.57143,0.71429,0.85714,1.00000]),
'y':array('d',[0.00000,0.08746,0.12448,0.16991,0.22235,0.26212,0.29875,0.33969,0.40000]),
'z':array('d',[300.00000,469.08356,548.33405,622.89789,891.49719,999.66022,1089.36560,1285.86804,1509.19348,1777.44641,2300.00000])},
't1_muminus_WJets_CR':{'x':array('d',[-1.00000,-0.75000,-0.50000,-0.25000,0.00000,0.25000,0.50000,0.75000,1.00000]),
'y':array('d',[0.00000,0.13933,0.28781,0.60000]),
'z':array('d',[500.00000,1226.66968,1505.61023,1796.29919,3000.00000])},
}

#Template class
class Template(object) :

	def __init__(self,name,formatted_name,modifier) :
		#A Template has a name and a formatted name
		print '	Adding template with name '+name
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
		self._XBINS = copy.copy(BINS[self.__name.split('__')[0]]['x'])
		self._YBINS = copy.copy(BINS[self.__name.split('__')[0]]['y'])
		self._ZBINS = copy.copy(BINS[self.__name.split('__')[0]]['z'])
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
				eweight = self.__get_event_weight__(branch_dict,const_rw_name_list,ss_rw_name_list,identifier.find('qcd')!=-1)
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
				#print 'doing entry %d of %d for template with name %s'%(entry,nEntries,self.getName()) #DEBUG
				tree.GetEntry(entry)
				#get the event weight
				eweight = self.__get_event_weight__(branch_dict,const_rw_name_list,ss_rw_list,ttree_identifier.find('qcd')!=-1)
				eweight_opp = eweight
				#print '	initial = %.4f'%(eweight) #DEBUG
				#multiply by the function weights if necessary
				if functionstring!=None and functionstring!='' :
					fweight, fweight_opp = self.__get_function_weight__(branch_dict,functionstring)
					eweight*=fweight; eweight_opp*=fweight_opp
				#print '	after function weight = %.4f'%(eweight) #DEBUG
				#multiply by the last +/-1 factor or whatever
				eweight*=extra_weight; eweight_opp*=extra_weight
				#print '	after extra weight = %.4f'%(eweight) #DEBUG
				#add to the event numbers
				if channelcharge==0 or branch_dict['lep_Q'].getPTreeValue()==channelcharge:
					event_numbers[realidentifier]+=eweight
				if branch_dict['addTwice'].getPTreeValue()==1 and (channelcharge==0 or branch_dict['lep_Q'].getPTreeValue()!=channelcharge) :
					event_numbers[realidentifier]+=eweight_opp
				#print 'for template %s with realidentifier %s eweight for entry %d = %.6f'%(self.getName(),realidentifier,entry,eweight) #DEBUG
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
			eweight = self.__get_event_weight__(branch_dict,const_rw_name_list,ss_rw_list,ttree_identifier.find('qcd')!=-1)
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
	def __get_event_weight__(self,branch_dict,const_rw_name_list,ss_rw_list,skipiso=False) :
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
				if ssname=='lep_iso_weight' and skipiso :
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