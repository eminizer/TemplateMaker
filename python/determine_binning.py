#imports
from ROOT import TChain, TMinuit, Long, Double, TFile, TH3D
from array import array
from datetime import date
from math import *
import multiprocessing
import numpy as np

#global variables
total_ttree_dir_path = '/uscms_data/d3/eminizer/ttbar_13TeV/CMSSW_8_1_0/src/Analysis/Reconstructor/test/total_ttree_files/'
BINS_IN_USE = {
#type-3
't3_muplus_SR':{'x':array('d',[-1.00000,-0.93333,-0.86667,-0.80000,-0.73333,-0.66667,-0.60000,-0.53333,-0.46667,-0.40000,-0.33333,-0.26667,-0.20000,-0.13333,-0.06667,0.00000,0.06667,0.13333,0.20000,0.26667,0.33333,0.40000,0.46667,0.53333,0.60000,0.66667,0.73333,0.80000,0.86667,0.93333,1.00000]),
'y':array('d',[0.00000,0.02165,0.04189,0.06237,0.07682,0.09849,0.11826,0.13516,0.15447,0.17486,0.20695,0.23776,0.30000]),
'z':array('d',[300.00000,372.87604,399.71985,423.32367,445.06909,466.55261,488.74268,513.13080,540.92096,575.17975,623.84082,707.41034,1500.00000])},
't3_muminus_SR':{'x':array('d',[-1.00000,-0.93333,-0.86667,-0.80000,-0.73333,-0.66667,-0.60000,-0.53333,-0.46667,-0.40000,-0.33333,-0.26667,-0.20000,-0.13333,-0.06667,0.00000,0.06667,0.13333,0.20000,0.26667,0.33333,0.40000,0.46667,0.53333,0.60000,0.66667,0.73333,0.80000,0.86667,0.93333,1.00000]),
'y':array('d',[0.00000,0.01695,0.03934,0.06121,0.07990,0.09236,0.11050,0.13230,0.15060,0.18037,0.21003,0.24133,0.30000]),
'z':array('d',[300.00000,372.99780,400.05176,423.24091,444.72675,466.15851,488.15814,512.28015,540.11108,574.47864,622.88477,707.37244,1500.00000])},
't3_elplus_SR':{'x':array('d',[-1.00000,-0.93333,-0.86667,-0.80000,-0.73333,-0.66667,-0.60000,-0.53333,-0.46667,-0.40000,-0.33333,-0.26667,-0.20000,-0.13333,-0.06667,0.00000,0.06667,0.13333,0.20000,0.26667,0.33333,0.40000,0.46667,0.53333,0.60000,0.66667,0.73333,0.80000,0.86667,0.93333,1.00000]),
'y':array('d',[0.00000,0.02983,0.04811,0.06513,0.08087,0.09450,0.10802,0.12414,0.14259,0.16794,0.19613,0.23064,0.30000]),
'z':array('d',[300.00000,374.85919,402.40796,425.89838,447.90344,469.81821,492.34473,517.20758,545.59766,580.57910,630.14075,716.79761,1500.00000])},
't3_elminus_SR':{'x':array('d',[-1.00000,-0.93333,-0.86667,-0.80000,-0.73333,-0.66667,-0.60000,-0.53333,-0.46667,-0.40000,-0.33333,-0.26667,-0.20000,-0.13333,-0.06667,0.00000,0.06667,0.13333,0.20000,0.26667,0.33333,0.40000,0.46667,0.53333,0.60000,0.66667,0.73333,0.80000,0.86667,0.93333,1.00000]),
'y':array('d',[0.00000,0.01904,0.03919,0.05701,0.07045,0.08302,0.09895,0.11783,0.14302,0.16337,0.18438,0.21541,0.30000]),
'z':array('d',[300.00000,374.47806,401.98401,425.61566,447.90192,469.48907,492.22876,517.14563,545.34753,581.00000,630.30835,716.62079,1500.00000])},
#type-2 signal region
't2_muplus_SR':{'x':array('d',[-1.00000,-0.85714,-0.71429,-0.57143,-0.42857,-0.28571,-0.14286,0.00000,0.14286,0.28571,0.42857,0.57143,0.71429,0.85714,1.00000]),
'y':array('d',[0.00000,0.03742,0.07834,0.11990,0.14789,0.18838,0.24034,0.30220,0.40000]),
'z':array('d',[300.00000,603.25977,765.57068,2300.00000])},
't2_muminus_SR':{'x':array('d',[-1.00000,-0.85714,-0.71429,-0.57143,-0.42857,-0.28571,-0.14286,0.00000,0.14286,0.28571,0.42857,0.57143,0.71429,0.85714,1.00000]),
'y':array('d',[0.00000,0.06020,0.09284,0.12212,0.16685,0.20663,0.25314,0.30168,0.40000]),
'z':array('d',[300.00000,601.98999,764.78754,2300.00000])},
't2_elplus_SR':{'x':array('d',[-1.00000,-0.85714,-0.71429,-0.57143,-0.42857,-0.28571,-0.14286,0.00000,0.14286,0.28571,0.42857,0.57143,0.71429,0.85714,1.00000]),
'y':array('d',[0.00000,0.04497,0.09096,0.12185,0.16669,0.22265,0.26073,0.34057,0.40000]),
'z':array('d',[300.00000,612.18353,862.03961,2300.00000])},
't2_elminus_SR':{'x':array('d',[-1.00000,-0.85714,-0.71429,-0.57143,-0.42857,-0.28571,-0.14286,0.00000,0.14286,0.28571,0.42857,0.57143,0.71429,0.85714,1.00000]),
'y':array('d',[0.00000,0.04086,0.08652,0.15359,0.20935,0.25600,0.29741,0.35061,0.40000]),
'z':array('d',[300.00000,607.87329,860.54993,2300.00000])},
#type-2 control region
't2_muplus_WJets_CR':{'x':array('d',[-1.00000,-0.85714,-0.71429,-0.57143,-0.42857,-0.28571,-0.14286,0.00000,0.14286,0.28571,0.42857,0.57143,0.71429,0.85714,1.00000]),
'y':array('d',[0.00000,0.07027,0.10914,0.14975,0.18860,0.23231,0.27742,0.32712,0.40000]),
'z':array('d',[300.00000,653.97552,820.31079,2300.00000])},
't2_muminus_WJets_CR':{'x':array('d',[-1.00000,-0.85714,-0.71429,-0.57143,-0.42857,-0.28571,-0.14286,0.00000,0.14286,0.28571,0.42857,0.57143,0.71429,0.85714,1.00000]),
'y':array('d',[0.00000,0.07377,0.11958,0.15177,0.18145,0.22464,0.27217,0.32715,0.40000]),
'z':array('d',[300.00000,654.21356,818.80200,2300.00000])},
't2_elplus_WJets_CR':{'x':array('d',[-1.00000,-0.85714,-0.71429,-0.57143,-0.42857,-0.28571,-0.14286,0.00000,0.14286,0.28571,0.42857,0.57143,0.71429,0.85714,1.00000]),
'y':array('d',[0.00000,0.09437,0.13660,0.18305,0.22429,0.26389,0.30472,0.35505,0.40000]),
'z':array('d',[300.00000,743.87195,982.92511,2300.00000])},
't2_elminus_WJets_CR':{'x':array('d',[-1.00000,-0.85714,-0.71429,-0.57143,-0.42857,-0.28571,-0.14286,0.00000,0.14286,0.28571,0.42857,0.57143,0.71429,0.85714,1.00000]),
'y':array('d',[0.00000,0.08645,0.12456,0.15790,0.18272,0.22561,0.28721,0.34180,0.40000]),
'z':array('d',[300.00000,742.55719,985.60760,2300.00000])},
#type-1 signal region
't1_muplus_SR':{'x':array('d',[-1.00000,-0.75000,-0.50000,-0.25000,0.00000,0.25000,0.50000,0.75000,1.00000]),
'y':array('d',[0.00000,0.12528,0.28831,0.60000]),
'z':array('d',[500.00000,3000.00000])},
't1_muminus_SR':{'x':array('d',[-1.00000,-0.75000,-0.50000,-0.25000,0.00000,0.25000,0.50000,0.75000,1.00000]),
'y':array('d',[0.00000,0.15001,0.30706,0.60000]),
'z':array('d',[500.00000,3000.00000])},
't1_elplus_SR':{'x':array('d',[-1.00000,-0.75000,-0.50000,-0.25000,0.00000,0.25000,0.50000,0.75000,1.00000]),
'y':array('d',[0.00000,0.19086,0.35487,0.60000]),
'z':array('d',[500.00000,3000.00000])},
't1_elminus_SR':{'x':array('d',[-1.00000,-0.75000,-0.50000,-0.25000,0.00000,0.25000,0.50000,0.75000,1.00000]),
'y':array('d',[0.00000,0.16227,0.31716,0.60000]),
'z':array('d',[500.00000,3000.00000])},
#type-1 control region
't1_muplus_WJets_CR':{'x':array('d',[-1.00000,-0.75000,-0.50000,-0.25000,0.00000,0.25000,0.50000,0.75000,1.00000]),
'y':array('d',[0.00000,0.11865,0.27983,0.60000]),
'z':array('d',[500.00000,3000.00000])},
't1_muminus_WJets_CR':{'x':array('d',[-1.00000,-0.75000,-0.50000,-0.25000,0.00000,0.25000,0.50000,0.75000,1.00000]),
'y':array('d',[0.00000,0.12102,0.28791,0.60000]),
'z':array('d',[500.00000,3000.00000])},
't1_elplus_WJets_CR':{'x':array('d',[-1.00000,-0.75000,-0.50000,-0.25000,0.00000,0.25000,0.50000,0.75000,1.00000]),
'y':array('d',[0.00000,0.16833,0.34323,0.60000]),
'z':array('d',[500.00000,3000.00000])},
't1_elminus_WJets_CR':{'x':array('d',[-1.00000,-0.75000,-0.50000,-0.25000,0.00000,0.25000,0.50000,0.75000,1.00000]),
'y':array('d',[0.00000,0.14534,0.31970,0.60000]),
'z':array('d',[500.00000,3000.00000])},
}


######			FILTER BINS FUNCTION			#####
def filter_bins(a) :
	#list of bin edge tuples
	bin_edge_tuples = []
	#first find the average and std. deviation of the bin widths
	bin_widths = [(a[i+1]-a[i]) for i in range(len(a)-1)]
	bin_widths_sq = [w**2 for w in bin_widths]
	avg_width = sum(bin_widths)/len(bin_widths)
	#avg_width = (a[-1]-a[0])/(len(a)-1)
	std_dev_bin_width = sqrt(abs(sum(bin_widths_sq)/len(bin_widths)-avg_width**2))
	#loop over bins
	for i in range(len(a)-1) :
		#get this bin's width
		width = a[i+1]-a[i]
		#add a tuple of this bin's edges if it's wide enough to keep
		#if width>(avg_width-2*std_dev_bin_width) :
		if width>(0.1*avg_width) :
			bin_edge_tuples.append((a[i],a[i+1]))
	#array of filtered bins to return
	filtered_bins = array('f',(len(bin_edge_tuples)+1)*[0.])
	filtered_bins[0]=a[0] #first LHS should be the minimum even if we're losing the first bin
	#now loop over the bin edge tuples and adjust them to be continuous
	for i in range(len(bin_edge_tuples)-1) :
		#push the right hand side of this bin, splitting the difference between its given RHS and the LHD of the next bin
		thisright = bin_edge_tuples[i][1]
		nextleft  = bin_edge_tuples[i+1][0]
		filtered_bins[i+1]=thisright+0.5*(nextleft-thisright)
	filtered_bins[len(bin_edge_tuples)]=a[-1] #last value again should be the end of the dataset even if we lose the last bin
	return filtered_bins

######			BAYESIAN BLOCKS FUNCTION			#####
def bayesian_blocks(t,const=4.) :

	# copy and sort the array
	t = np.sort(t)
	N = t.size

	# create length-(N + 1) array of cell edges
	edges = np.concatenate([t[:1],
							0.5 * (t[1:] + t[:-1]),
							t[-1:]])
	block_length = t[-1] - edges

	# arrays needed for the iteration
	nn_vec = np.ones(N)
	best = np.zeros(N, dtype=float)
	last = np.zeros(N, dtype=int)

	#-----------------------------------------------------------------
	# Start with first data cell; add one cell at each iteration
	#-----------------------------------------------------------------
	for K in range(N):
		# Compute the width and count of the final bin for all possible
		# locations of the K^th changepoint
		width = block_length[:K + 1] - block_length[K + 1]
		count_vec = np.cumsum(nn_vec[:K + 1][::-1])[::-1]

		# evaluate fitness function for these possibilities
		fit_vec = count_vec * (np.log(count_vec) - np.log(width))
		fit_vec -= const  # 4 comes from the prior on the number of changepoints
		fit_vec[1:] += best[:K]

		# find the max of the fitness: this is the K^th changepoint
		i_max = np.argmax(fit_vec)
		last[K] = i_max
		best[K] = fit_vec[i_max]

	#-----------------------------------------------------------------
	# Recover changepoints by iteratively peeling off the last block
	#-----------------------------------------------------------------
	change_points =  np.zeros(N, dtype=int)
	i_cp = N
	ind = N
	while True:
		i_cp -= 1
		change_points[i_cp] = ind
		if ind == 0:
			break
		ind = last[ind - 1]
	change_points = change_points[i_cp:]

	return filter_bins(edges[change_points])

#channel class
class Channel(object) :

	#init
	def __init__(self,name,cutstring,nxbins,xmin,xmax,nybins,ymin,ymax,nzbins,zmin,zmax,inityconst=1.,initzconst=1.) :
		self._name = name
		self._cutstring = cutstring
		self._nxbins = nxbins
		self._xmin = xmin
		self._xmax = xmax
		self._nybins = nybins
		self._ymin = ymin
		self._ymax = ymax
		self._nzbins = nzbins
		self._zmin = zmin
		self._zmax = zmax
		self._inityconst = inityconst
		self._initzconst = initzconst
		self._data_chain = None
		self._skimmed_tree = None
		self._x_data = None
		self._y_data = None
		self._z_data = None
		self._best_fit_x_bins = None
		self._best_fit_y_bins = None
		self._best_fit_z_bins = None

	######			PUBLIC FUNCTIONS			#####
	def skimDataChain(self) :
		self._skimmed_tree = self._data_chain.CopyTree(self._cutstring)

	def readSkimmedTreeIntoArrays(self) :
		print 'reading skimmed tree into data arrays in channel %s'%(self._name)
		#number of data points is number of entries in the tree
		nEntries = self._skimmed_tree.GetEntries()
		x_data = np.zeros(nEntries,dtype=np.dtype(float))
		y_data = np.zeros(nEntries,dtype=np.dtype(float))
		z_data = np.zeros(nEntries,dtype=np.dtype(float))
		#declare and set branch addresses
		cstar = array('f',[0.])
		x_F = array('f',[0.])
		M = array('f',[0.])
		self._skimmed_tree.SetBranchAddress('cstar',cstar)
		self._skimmed_tree.SetBranchAddress('x_F',x_F)
		self._skimmed_tree.SetBranchAddress('M',M) 
		#fill data arrays
		tossed = 0
		realcount = 0
		for i in range(nEntries) :
			self._skimmed_tree.GetEntry(i)
			if cstar[0]>self._xmin and cstar[0]<self._xmax and abs(x_F[0])>self._ymin and abs(x_F[0])<self._ymax and M[0]>self._zmin and M[0]<self._zmax :
				x_data[realcount]=cstar[0]
				y_data[realcount]=abs(x_F[0])
				z_data[realcount]=abs(M[0])
				realcount+=1
			else :
				#print '	tossing data point with (cstar, |x_F|, M) = (%.4f, %.4f, %.4f)'%(cstar[0],abs(x_F[0]),M[0]) #DEBUG
				tossed+=1
		print 'removed %d out of %d events (%.3f%%) in channel %s'%(tossed,nEntries,100.*(tossed/nEntries),self._name)
		#copy into properly-sized arrays 
		self._x_data = x_data[:realcount]
		self._y_data = y_data[:realcount]
		self._z_data = z_data[:realcount]

	def doBayesianBlocks(self) :
		print 'doing bayesian_blocks in channel %s'%(self._name)
		#use even-width bins in the cstar dimension
		xbinwidth = (self._xmax-self._xmin)/(self._nxbins)
		bb_x_bins = array('f',(self._nxbins+1)*[0.])
		for i in range(len(bb_x_bins)) :
			bb_x_bins[i] = self._xmin+i*xbinwidth
		print 'bb_x_bins = %s'%(bb_x_bins)
	#	#x_F dimension uses the Bayesian Blocks algorithm
	#	const=self._inityconst; lastlen=500; switches=0
	#	while True :
	#		bb_y_bins = bayesian_blocks(self._y_data,const)
	#		#bb_y_bins = bayesian_blocks([self._y_data[i] for i in range(0,len(self._y_data),10)],const) #shortened bayesian blocks array to reduce runtime
	#		print '	y const = %.6f, nybins = %d'%(const,len(bb_y_bins)-1)
	#		if len(bb_y_bins)==self._nybins+1 or switches>4 :
	#			break
	#		elif len(bb_y_bins)>self._nybins+1 :
	#			if lastlen<self._nybins+1 :
	#				switches+=1
	#			const+=(10./(10**switches))
	#		elif len(bb_y_bins)<self._nybins+1 :
	#			if lastlen>self._nybins+1 :
	#				switches+=1
	#			const-=(10./(10**switches))
	#		lastlen=len(bb_y_bins)
	#	bb_y_bins[0] = self._ymin; bb_y_bins[-1] = self._ymax
	#	print 'bb_y_bins = %s'%(bb_y_bins)
		#dummy y-dimension returning bins in use
		bb_y_bins = BINS_IN_USE[self._name]['y']
		print 'bb_y_bins = %s'%(bb_y_bins)
		#M dimension returns equal-content bins
		bb_z_bins = array('f',(self._nzbins+1)*[0.])
		self._z_data = np.sort(self._z_data)
		nEventsPerZBin = int(round(len(self._z_data)/self._nzbins))
		bb_z_bins[0] = self._zmin
		for i in range(1,self._nzbins) :
			bb_z_bins[i]=0.5*(self._z_data[i*nEventsPerZBin-1]+self._z_data[i*nEventsPerZBin])
		bb_z_bins[-1] = self._zmax
		print 'bb_z_bins = %s'%(bb_z_bins)
		#copy results into the binning arrays
		self._best_fit_x_bins = array('d',len(bb_x_bins)*[0.])
		for i in range(len(bb_x_bins)) :
			self._best_fit_x_bins[i]=bb_x_bins[i]
		self._best_fit_y_bins = array('d',len(bb_y_bins)*[0.])
		for i in range(len(bb_y_bins)) :
			self._best_fit_y_bins[i]=bb_y_bins[i]
		self._best_fit_z_bins = array('d',len(bb_z_bins)*[0.])
		for i in range(len(bb_z_bins)) :
			self._best_fit_z_bins[i]=bb_z_bins[i]

	def getBestFitHistoList(self) :
		#declare 3D histogram
		histoname = '%s_best_fit_3D_histo'%(self._name)
		self._best_fit_histo = TH3D(histoname,'%s channel best fit 3D histogram; c*; |x_{F}|; M [GeV]'%(self._name),
									len(self._best_fit_x_bins)-1,self._best_fit_x_bins,len(self._best_fit_y_bins)-1,self._best_fit_y_bins,
									len(self._best_fit_z_bins)-1,self._best_fit_z_bins)
		#draw skimmed tree into this histogram
		cstar = array('f',[-1.0]);  self._skimmed_tree.SetBranchAddress('cstar',cstar)
		x_F   = array('f',[-1.0]);  self._skimmed_tree.SetBranchAddress('x_F',x_F)
		M 	  = array('f',[-1.0]);  self._skimmed_tree.SetBranchAddress('M',M) 
		for i in range(self._skimmed_tree.GetEntries()) :
			self._skimmed_tree.GetEntry(i)
			self._best_fit_histo.Fill(cstar[0],abs(x_F[0]),M[0])
		#return the histogram and its projections
		return self._best_fit_histo, self._best_fit_histo.ProjectionX(), self._best_fit_histo.ProjectionY(), self._best_fit_histo.ProjectionZ()

	def getCodeLineList(self) :
		#x bins
		x_bins_string = "'x':array('d',["
		for i in range(len(self._best_fit_x_bins)) :
			x_bins_string+='%.5f'%(self._best_fit_x_bins[i])
			if i<self._nxbins :
				x_bins_string+=','
			else :
				x_bins_string+='])'
		#y bins
		y_bins_string = "'y':array('d',["
		for i in range(len(self._best_fit_y_bins)) :
			y_bins_string+='%.5f'%(self._best_fit_y_bins[i])
			if i<self._nybins :
				y_bins_string+=','
			else :
				y_bins_string+='])'
		#z bins
		z_bins_string = "'z':array('d',["
		for i in range(len(self._best_fit_z_bins)) :
			z_bins_string+='%.5f'%(self._best_fit_z_bins[i])
			if i<self._nzbins :
				z_bins_string+=','
			else :
				z_bins_string+='])'
		return x_bins_string, y_bins_string, z_bins_string


	######			GETTERS/SETTERS			#####
	def setSkimmedTree(self,tree) :
		self._skimmed_tree = tree
	def getSkimmedTree(self) :
		return self._skimmed_tree
	def setDataChain(self,chain) :
		self._data_chain = chain
	def getName(self) :
		return self._name

#add channels
all_channels = {}
#signal region
#charge separated
all_channels['t1_muplus_SR']  = Channel('t1_muplus_SR','eventTopology==1 && lepflavor==1 && lep_Q>0 && fullselection==1',8,-1.,1.,3,0.,0.6,1,500.,3000.,100.,45.)
all_channels['t1_muminus_SR'] = Channel('t1_muminus_SR','eventTopology==1 && lepflavor==1 && lep_Q<0 && fullselection==1',8,-1.,1.,3,0.,0.6,1,500.,3000.,55.,45.)
all_channels['t1_elplus_SR']  = Channel('t1_elplus_SR','eventTopology==1 && lepflavor==2 && lep_Q>0 && fullselection==1',8,-1.,1.,3,0.,0.6,1,500.,3000.,45.,35.)
all_channels['t1_elminus_SR'] = Channel('t1_elminus_SR','eventTopology==1 && lepflavor==2 && lep_Q<0 && fullselection==1',8,-1.,1.,3,0.,0.6,1,500.,3000.,25.,25.)
all_channels['t2_muplus_SR']  = Channel('t2_muplus_SR','eventTopology==2 && lepflavor==1 && lep_Q>0 && fullselection==1',14,-1.,1.,8,0.,0.4,3,300.,2300.,50.)
all_channels['t2_muminus_SR'] = Channel('t2_muminus_SR','eventTopology==2 && lepflavor==1 && lep_Q<0 && fullselection==1',14,-1.,1.,8,0.,0.4,3,300.,2300.,50.)
all_channels['t2_elplus_SR']  = Channel('t2_elplus_SR','eventTopology==2 && lepflavor==2 && lep_Q>0 && fullselection==1',14,-1.,1.,8,0.,0.4,3,300.,2300.,5.)
all_channels['t2_elminus_SR'] = Channel('t2_elminus_SR','eventTopology==2 && lepflavor==2 && lep_Q<0 && fullselection==1',14,-1.,1.,8,0.,0.4,3,300.,2300.,4.)
all_channels['t3_muplus_SR']  = Channel('t3_muplus_SR','eventTopology==3 && lepflavor==1 && lep_Q>0 && fullselection==1',30,-1.,1.,12,0.,0.3,12,300.,1500.,6.1)
all_channels['t3_muminus_SR'] = Channel('t3_muminus_SR','eventTopology==3 && lepflavor==1 && lep_Q<0 && fullselection==1',30,-1.,1.,12,0.,0.3,12,300.,1500.,6.1)
all_channels['t3_elplus_SR']  = Channel('t3_elplus_SR','eventTopology==3 && lepflavor==2 && lep_Q>0 && fullselection==1',30,-1.,1.,12,0.,0.3,12,300.,1500.,110.)
all_channels['t3_elminus_SR'] = Channel('t3_elminus_SR','eventTopology==3 && lepflavor==2 && lep_Q<0 && fullselection==1',30,-1.,1.,12,0.,0.3,12,300.,1500.,6.7)#118.)
#charge summed
#all_channels['t1_mu_SR']  = Channel('t1_mu_SR','eventTopology==1 && lepflavor==1 && fullselection==1',8,-1.,1.,3,0.,0.6,1,500.,3000.,100.,45.)
#all_channels['t1_el_SR']  = Channel('t1_el_SR','eventTopology==1 && lepflavor==2 && fullselection==1',8,-1.,1.,3,0.,0.6,1,500.,3000.,45.,35.)
#all_channels['t2_mu_SR']  = Channel('t2_mu_SR','eventTopology==2 && lepflavor==1 && fullselection==1',14,-1.,1.,8,0.,0.4,3,300.,2300.,140.,5.)
#all_channels['t2_el_SR']  = Channel('t2_el_SR','eventTopology==2 && lepflavor==2 && fullselection==1',14,-1.,1.,8,0.,0.4,3,300.,2300.,10.,10.)
#all_channels['t3_mu_SR']  = Channel('t3_mu_SR','eventTopology==3 && lepflavor==1 && fullselection==1',20,-1.,1.,9,0.,0.3,6,300.,1500.,5.,2.)
#all_channels['t3_el_SR']  = Channel('t3_el_SR','eventTopology==3 && lepflavor==2 && fullselection==1',20,-1.,1.,9,0.,0.3,6,300.,1500.,400.,50.)
##boosted W+Jets control regions
#charge separated
all_channels['t1_muplus_WJets_CR']  = Channel('t1_muplus_WJets_CR','eventTopology==1 && lepflavor==1 && lep_Q>0 && wjets_cr_selection==1',8,-1.,1.,3,0.,0.6,1,500.,3000.,35.,25.)
all_channels['t1_muminus_WJets_CR'] = Channel('t1_muminus_WJets_CR','eventTopology==1 && lepflavor==1 && lep_Q<0 && wjets_cr_selection==1',8,-1.,1.,3,0.,0.6,1,500.,3000.,25.,15.)
all_channels['t1_elplus_WJets_CR']  = Channel('t1_elplus_WJets_CR','eventTopology==1 && lepflavor==2 && lep_Q>0 && wjets_cr_selection==1',8,-1.,1.,3,0.,0.6,1,500.,3000.,15.,25.)
all_channels['t1_elminus_WJets_CR'] = Channel('t1_elminus_WJets_CR','eventTopology==1 && lepflavor==2 && lep_Q<0 && wjets_cr_selection==1',8,-1.,1.,3,0.,0.6,1,500.,3000.,25.,8.)
all_channels['t2_muplus_WJets_CR']  = Channel('t2_muplus_WJets_CR','eventTopology==2 && lepflavor==1 && lep_Q>0 && wjets_cr_selection==1',14,-1.,1.,8,0.,0.4,3,300.,2300.,50.)
all_channels['t2_muminus_WJets_CR'] = Channel('t2_muminus_WJets_CR','eventTopology==2 && lepflavor==1 && lep_Q<0 && wjets_cr_selection==1',14,-1.,1.,8,0.,0.4,3,300.,2300.,25.)
all_channels['t2_elplus_WJets_CR']  = Channel('t2_elplus_WJets_CR','eventTopology==2 && lepflavor==2 && lep_Q>0 && wjets_cr_selection==1',14,-1.,1.,8,0.,0.4,3,300.,2300.,5.)
all_channels['t2_elminus_WJets_CR'] = Channel('t2_elminus_WJets_CR','eventTopology==2 && lepflavor==2 && lep_Q<0 && wjets_cr_selection==1',14,-1.,1.,8,0.,0.4,3,300.,2300.,4.)
#charge summed
#all_channels['t1_mu_WJets_CR']  = Channel('t1_mu_WJets_CR','eventTopology==1 && lepflavor==1 && wjets_cr_selection==1',8,-1.,1.,3,0.,0.6,4,500.,3000.,35.,25.)
#all_channels['t1_el_WJets_CR']  = Channel('t1_el_WJets_CR','eventTopology==1 && lepflavor==2 && wjets_cr_selection==1',8,-1.,1.,3,0.,0.6,4,500.,3000.,15.,25.)
#all_channels['t2_mu_WJets_CR']  = Channel('t2_mu_WJets_CR','eventTopology==2 && lepflavor==1 && wjets_cr_selection==1',14,-1.,1.,8,0.,0.4,10,300.,2300.,25.,65.)
#all_channels['t2_el_WJets_CR']  = Channel('t2_el_WJets_CR','eventTopology==2 && lepflavor==2 && wjets_cr_selection==1',14,-1.,1.,8,0.,0.4,10,300.,2300.,6.,18.)

#data sample names
mu_data_samplenames = []
mu_data_samplenames.append('SingleMu_Run2016Bv2')
mu_data_samplenames.append('SingleMu_Run2016C')
mu_data_samplenames.append('SingleMu_Run2016D')
mu_data_samplenames.append('SingleMu_Run2016E')
mu_data_samplenames.append('SingleMu_Run2016F')
mu_data_samplenames.append('SingleMu_Run2016G')
mu_data_samplenames.append('SingleMu_Run2016Hv2')
mu_data_samplenames.append('SingleMu_Run2016Hv3')
el_data_samplenames = []
el_data_samplenames.append('SingleEl_Run2016Bv2')
el_data_samplenames.append('SingleEl_Run2016C')
el_data_samplenames.append('SingleEl_Run2016D')
el_data_samplenames.append('SingleEl_Run2016E')
el_data_samplenames.append('SingleEl_Run2016F')
el_data_samplenames.append('SingleEl_Run2016G')
el_data_samplenames.append('SingleEl_Run2016Hv2')
el_data_samplenames.append('SingleEl_Run2016Hv3')

#data chains
mu_data_chain = TChain('tree')
for n in mu_data_samplenames :
	mu_data_chain.Add(total_ttree_dir_path+n+'_skim_all.root')
el_data_chain = TChain('tree')
for n in el_data_samplenames :
	el_data_chain.Add(total_ttree_dir_path+n+'_skim_all.root')

#set data chains
for channame in all_channels.keys() :
	if channame.find('_mu')!=-1 :
		all_channels[channame].setDataChain(mu_data_chain)
	elif channame.find('_el')!=-1 :
		all_channels[channame].setDataChain(el_data_chain)
	else :
		print 'ERROR: failed to determine leptype for chain with name '+channame

#skim data chains in parallel
def skimChainsParallel(channel,child_tree) :
	print 'skimming chain for channel '+channel.getName()+'...'
	channel.skimDataChain()
	print 'done skimming chain for channel '+channel.getName()+', sending tree'
	child_tree.send(channel.getSkimmedTree())
	print 'done sending skimmed tree for channel '+channel.getName()
	child_tree.close()
procs = []; pipe_ends = []
for channel in all_channels.values() :
	p_t, c_t = multiprocessing.Pipe()
	pipe_ends.append((p_t,channel))
	p = multiprocessing.Process(target=skimChainsParallel,args=(channel,c_t))
	p.start()
	procs.append(p)
for pipe_end in pipe_ends :
	pipe_end[1].setSkimmedTree(pipe_end[0].recv())
for proc in procs :
	proc.join()

#use chains to build data arrays
for channel in all_channels.values() :
	channel.readSkimmedTreeIntoArrays()

#use the data arrays to determine the optimal binning
for channel in all_channels.values() :
	channel.doBayesianBlocks()

#write results to files
outrootfile = TFile('binning_test_histos_'+str(date.today())+'.root','recreate')
outtextfile = open('binning_test_codelines_'+str(date.today())+'.txt','w')
outtextfile.write('BINS = {')
for channel in all_channels.values() :
	print 'saving results for channel %s'%(channel.getName())
	#save data histograms with the given binning
	histolist = channel.getBestFitHistoList()
	for histo in histolist :
		histo.Write()
	#save lines of code to copy-paste into template.py
	outtextfile.write("'%s':{"%(channel.getName()))
	codelinelist = channel.getCodeLineList()
	for i in range(len(codelinelist)) :
		outtextfile.write("%s"%(codelinelist[i]))
		if i==len(codelinelist)-1 :
			outtextfile.write('}')
		outtextfile.write(',\n')
outtextfile.write('\n}\n')
outrootfile.Close()
outtextfile.close()
print 'All done!!'
