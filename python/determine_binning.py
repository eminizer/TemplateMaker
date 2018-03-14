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
'y':array('d',[0.00000,0.02758,0.04819,0.06722,0.08191,0.10083,0.11415,0.13403,0.15670,0.17709,0.20166,0.24187,0.30000]),
'z':array('d',[300.00000,403.26733,442.16650,475.80878,509.84674,548.18164,598.65088,683.80573,1500.00000])},
't3_muminus_SR':{'x':array('d',[-1.00000,-0.93333,-0.86667,-0.80000,-0.73333,-0.66667,-0.60000,-0.53333,-0.46667,-0.40000,-0.33333,-0.26667,-0.20000,-0.13333,-0.06667,0.00000,0.06667,0.13333,0.20000,0.26667,0.33333,0.40000,0.46667,0.53333,0.60000,0.66667,0.73333,0.80000,0.86667,0.93333,1.00000]),
'y':array('d',[0.00000,0.02758,0.04819,0.06722,0.08191,0.10083,0.11415,0.13403,0.15670,0.17709,0.20166,0.24187,0.30000]),
'z':array('d',[300.00000,402.91043,440.96478,474.33536,508.21732,546.35938,596.99371,682.79004,1500.00000])},
't3_elplus_SR':{'x':array('d',[-1.00000,-0.93333,-0.86667,-0.80000,-0.73333,-0.66667,-0.60000,-0.53333,-0.46667,-0.40000,-0.33333,-0.26667,-0.20000,-0.13333,-0.06667,0.00000,0.06667,0.13333,0.20000,0.26667,0.33333,0.40000,0.46667,0.53333,0.60000,0.66667,0.73333,0.80000,0.86667,0.93333,1.00000]),
'y':array('d',[0.00000,0.02601,0.04264,0.06316,0.08051,0.10235,0.11958,0.14073,0.16257,0.18803,0.21321,0.24119,0.30000]),
'z':array('d',[300.00000,404.68713,443.48633,477.28915,512.00159,551.66168,603.56793,693.08807,1500.00000])},
't3_elminus_SR':{'x':array('d',[-1.00000,-0.93333,-0.86667,-0.80000,-0.73333,-0.66667,-0.60000,-0.53333,-0.46667,-0.40000,-0.33333,-0.26667,-0.20000,-0.13333,-0.06667,0.00000,0.06667,0.13333,0.20000,0.26667,0.33333,0.40000,0.46667,0.53333,0.60000,0.66667,0.73333,0.80000,0.86667,0.93333,1.00000]),
'y':array('d',[0.00000,0.01780,0.03907,0.05983,0.07829,0.09679,0.11921,0.14017,0.16716,0.19435,0.21266,0.25430,0.30000]),
'z':array('d',[300.00000,404.60132,443.99030,478.13617,512.50482,551.60370,604.63599,694.63232,1500.00000])},
#type-2 signal region
't2_muplus_SR':{'x':array('d',[-1.00000,-0.85714,-0.71429,-0.57143,-0.42857,-0.28571,-0.14286,0.00000,0.14286,0.28571,0.42857,0.57143,0.71429,0.85714,1.00000]),
'y':array('d',[0.00000,0.03626,0.07458,0.11004,0.14555,0.18790,0.24134,0.30371,0.40000]),
'z':array('d',[300.00000,580.57172,751.94049,2300.00000])},
't2_muminus_SR':{'x':array('d',[-1.00000,-0.85714,-0.71429,-0.57143,-0.42857,-0.28571,-0.14286,0.00000,0.14286,0.28571,0.42857,0.57143,0.71429,0.85714,1.00000]),
'y':array('d',[0.00000,0.02928,0.06229,0.09681,0.12930,0.16976,0.21690,0.29217,0.40000]),
'z':array('d',[300.00000,579.56177,750.71912,2300.00000])},
't2_elplus_SR':{'x':array('d',[-1.00000,-0.85714,-0.71429,-0.57143,-0.42857,-0.28571,-0.14286,0.00000,0.14286,0.28571,0.42857,0.57143,0.71429,0.85714,1.00000]),
'y':array('d',[0.00000,0.04146,0.08171,0.11894,0.16248,0.21010,0.26052,0.32987,0.40000]),
'z':array('d',[300.00000,582.40131,824.45508,2300.00000])},
't2_elminus_SR':{'x':array('d',[-1.00000,-0.85714,-0.71429,-0.57143,-0.42857,-0.28571,-0.14286,0.00000,0.14286,0.28571,0.42857,0.57143,0.71429,0.85714,1.00000]),
'y':array('d',[0.00000,0.02095,0.06528,0.10061,0.15173,0.21084,0.26978,0.31276,0.40000]),
'z':array('d',[300.00000,571.34265,820.16565,2300.00000])},
#type-2 control region
't2_muplus_WJets_CR':{'x':array('d',[-1.00000,-0.85714,-0.71429,-0.57143,-0.42857,-0.28571,-0.14286,0.00000,0.14286,0.28571,0.42857,0.57143,0.71429,0.85714,1.00000]),
'y':array('d',[0.00000,0.07091,0.11474,0.14652,0.18781,0.21799,0.24784,0.31771,0.40000]),
'z':array('d',[300.00000,647.65326,814.10364,2300.00000])},
't2_muminus_WJets_CR':{'x':array('d',[-1.00000,-0.85714,-0.71429,-0.57143,-0.42857,-0.28571,-0.14286,0.00000,0.14286,0.28571,0.42857,0.57143,0.71429,0.85714,1.00000]),
'y':array('d',[0.00000,0.07390,0.10424,0.13915,0.18314,0.23072,0.28350,0.33034,0.40000]),
'z':array('d',[300.00000,647.70184,813.19189,2300.00000])},
't2_elplus_WJets_CR':{'x':array('d',[-1.00000,-0.85714,-0.71429,-0.57143,-0.42857,-0.28571,-0.14286,0.00000,0.14286,0.28571,0.42857,0.57143,0.71429,0.85714,1.00000]),
'y':array('d',[0.00000,0.08117,0.12181,0.16974,0.22103,0.25909,0.29891,0.35423,0.40000]),
'z':array('d',[300.00000,726.63721,949.88147,2300.00000])},
't2_elminus_WJets_CR':{'x':array('d',[-1.00000,-0.85714,-0.71429,-0.57143,-0.42857,-0.28571,-0.14286,0.00000,0.14286,0.28571,0.42857,0.57143,0.71429,0.85714,1.00000]),
'y':array('d',[0.00000,0.09406,0.13573,0.17199,0.20903,0.24973,0.30305,0.37694,0.40000]),
'z':array('d',[300.00000,726.50244,964.18347,2300.00000])},
#type-1 signal region
't1_muplus_SR':{'x':array('d',[-1.00000,-0.75000,-0.50000,-0.25000,0.00000,0.25000,0.50000,0.75000,1.00000]),
'y':array('d',[0.00000,0.13030,0.30785,0.60000]),
'z':array('d',[500.00000,3000.00000])},
't1_muminus_SR':{'x':array('d',[-1.00000,-0.75000,-0.50000,-0.25000,0.00000,0.25000,0.50000,0.75000,1.00000]),
'y':array('d',[0.00000,0.15324,0.30814,0.60000]),
'z':array('d',[500.00000,3000.00000])},
't1_elplus_SR':{'x':array('d',[-1.00000,-0.75000,-0.50000,-0.25000,0.00000,0.25000,0.50000,0.75000,1.00000]),
'y':array('d',[0.00000,0.16864,0.35377,0.60000]),
'z':array('d',[500.00000,3000.00000])},
't1_elminus_SR':{'x':array('d',[-1.00000,-0.75000,-0.50000,-0.25000,0.00000,0.25000,0.50000,0.75000,1.00000]),
'y':array('d',[0.00000,0.17003,0.33829,0.60000]),
'z':array('d',[500.00000,3000.00000])},
#type-1 control region
't1_muplus_WJets_CR':{'x':array('d',[-1.00000,-0.75000,-0.50000,-0.25000,0.00000,0.25000,0.50000,0.75000,1.00000]),
'y':array('d',[0.00000,0.13706,0.30875,0.60000]),
'z':array('d',[500.00000,3000.00000])},
't1_muminus_WJets_CR':{'x':array('d',[-1.00000,-0.75000,-0.50000,-0.25000,0.00000,0.25000,0.50000,0.75000,1.00000]),
'y':array('d',[0.00000,0.13939,0.30110,0.60000]),
'z':array('d',[500.00000,3000.00000])},
't1_elplus_WJets_CR':{'x':array('d',[-1.00000,-0.75000,-0.50000,-0.25000,0.00000,0.25000,0.50000,0.75000,1.00000]),
'y':array('d',[0.00000,0.14230,0.35644,0.60000]),
'z':array('d',[500.00000,3000.00000])},
't1_elminus_WJets_CR':{'x':array('d',[-1.00000,-0.75000,-0.50000,-0.25000,0.00000,0.25000,0.50000,0.75000,1.00000]),
'y':array('d',[0.00000,0.14625,0.31357,0.60000]),
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
		##x_F dimension uses the Bayesian Blocks algorithm
		#const=self._inityconst; lastlen=500; switches=0
		#while True :
		#	bb_y_bins = bayesian_blocks(self._y_data,const)
		#	print '	y const = %.6f, nybins = %d'%(const,len(bb_y_bins)-1)
		#	if len(bb_y_bins)==self._nybins+1 or switches>4 :
		#		break
		#	elif len(bb_y_bins)>self._nybins+1 :
		#		if lastlen<self._nybins+1 :
		#			switches+=1
		#		const+=(10./(10**switches))
		#	elif len(bb_y_bins)<self._nybins+1 :
		#		if lastlen>self._nybins+1 :
		#			switches+=1
		#		const-=(10./(10**switches))
		#	lastlen=len(bb_y_bins)
		#bb_y_bins[0] = self._ymin; bb_y_bins[-1] = self._ymax
		#print 'bb_y_bins = %s'%(bb_y_bins)
		#dummy y-dimension returning bins in use
		bb_y_bins = BINS_IN_USE[self._name]['y']
		print 'bb_y_bins = %s'%(bb_y_bins)
		#M dimension returns equal-content bins
		bb_z_bins = array('f',(self._nzbins+1)*[0.])
		self._z_data = np.sort(self._z_data)
		nEventsPerZBin = int(round(len(self._z_data)/self._nzbins))
		bb_z_bins[0] = self._zmin
		for i in range(1,self._nzbins) :
			bb_z_bins[i]=0.5*(self._z_data[i*nEventsPerZBin]+self._z_data[i*nEventsPerZBin+1])
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
#all_channels['t1_muplus_SR']  = Channel('t1_muplus_SR','eventTopology==1 && lepflavor==1 && lep_Q>0 && fullselection==1',8,-1.,1.,3,0.,0.6,1,500.,3000.,100.,45.)
#all_channels['t1_muminus_SR'] = Channel('t1_muminus_SR','eventTopology==1 && lepflavor==1 && lep_Q<0 && fullselection==1',8,-1.,1.,3,0.,0.6,1,500.,3000.,55.,45.)
#all_channels['t1_elplus_SR']  = Channel('t1_elplus_SR','eventTopology==1 && lepflavor==2 && lep_Q>0 && fullselection==1',8,-1.,1.,3,0.,0.6,1,500.,3000.,45.,35.)
#all_channels['t1_elminus_SR'] = Channel('t1_elminus_SR','eventTopology==1 && lepflavor==2 && lep_Q<0 && fullselection==1',8,-1.,1.,3,0.,0.6,1,500.,3000.,25.,25.)
#all_channels['t2_muplus_SR']  = Channel('t2_muplus_SR','eventTopology==2 && lepflavor==1 && lep_Q>0 && fullselection==1',14,-1.,1.,8,0.,0.4,3,300.,2300.,50.)
#all_channels['t2_muminus_SR'] = Channel('t2_muminus_SR','eventTopology==2 && lepflavor==1 && lep_Q<0 && fullselection==1',14,-1.,1.,8,0.,0.4,3,300.,2300.,50.)
#all_channels['t2_elplus_SR']  = Channel('t2_elplus_SR','eventTopology==2 && lepflavor==2 && lep_Q>0 && fullselection==1',14,-1.,1.,8,0.,0.4,3,300.,2300.,5.)
#all_channels['t2_elminus_SR'] = Channel('t2_elminus_SR','eventTopology==2 && lepflavor==2 && lep_Q<0 && fullselection==1',14,-1.,1.,8,0.,0.4,3,300.,2300.,4.)
all_channels['t3_muplus_SR']  = Channel('t3_muplus_SR','eventTopology==3 && lepflavor==1 && lep_Q>0 && fullselection==1',30,-1.,1.,12,0.,0.3,8,300.,1500.,2.)
all_channels['t3_muminus_SR'] = Channel('t3_muminus_SR','eventTopology==3 && lepflavor==1 && lep_Q<0 && fullselection==1',30,-1.,1.,12,0.,0.3,8,300.,1500.,71.)
all_channels['t3_elplus_SR']  = Channel('t3_elplus_SR','eventTopology==3 && lepflavor==2 && lep_Q>0 && fullselection==1',30,-1.,1.,12,0.,0.3,8,300.,1500.,58.)
all_channels['t3_elminus_SR'] = Channel('t3_elminus_SR','eventTopology==3 && lepflavor==2 && lep_Q<0 && fullselection==1',30,-1.,1.,12,0.,0.3,8,300.,1500.,61.)
#charge summed
#all_channels['t1_mu_SR']  = Channel('t1_mu_SR','eventTopology==1 && lepflavor==1 && fullselection==1',8,-1.,1.,3,0.,0.6,1,500.,3000.,100.,45.)
#all_channels['t1_el_SR']  = Channel('t1_el_SR','eventTopology==1 && lepflavor==2 && fullselection==1',8,-1.,1.,3,0.,0.6,1,500.,3000.,45.,35.)
#all_channels['t2_mu_SR']  = Channel('t2_mu_SR','eventTopology==2 && lepflavor==1 && fullselection==1',14,-1.,1.,8,0.,0.4,3,300.,2300.,140.,5.)
#all_channels['t2_el_SR']  = Channel('t2_el_SR','eventTopology==2 && lepflavor==2 && fullselection==1',14,-1.,1.,8,0.,0.4,3,300.,2300.,10.,10.)
#all_channels['t3_mu_SR']  = Channel('t3_mu_SR','eventTopology==3 && lepflavor==1 && fullselection==1',20,-1.,1.,9,0.,0.3,6,300.,1500.,5.,2.)
#all_channels['t3_el_SR']  = Channel('t3_el_SR','eventTopology==3 && lepflavor==2 && fullselection==1',20,-1.,1.,9,0.,0.3,6,300.,1500.,400.,50.)
##boosted W+Jets control regions
#charge separated
#all_channels['t1_muplus_WJets_CR']  = Channel('t1_muplus_WJets_CR','eventTopology==1 && lepflavor==1 && lep_Q>0 && wjets_cr_selection==1',8,-1.,1.,3,0.,0.6,1,500.,3000.,35.,25.)
#all_channels['t1_muminus_WJets_CR'] = Channel('t1_muminus_WJets_CR','eventTopology==1 && lepflavor==1 && lep_Q<0 && wjets_cr_selection==1',8,-1.,1.,3,0.,0.6,1,500.,3000.,25.,15.)
#all_channels['t1_elplus_WJets_CR']  = Channel('t1_elplus_WJets_CR','eventTopology==1 && lepflavor==2 && lep_Q>0 && wjets_cr_selection==1',8,-1.,1.,3,0.,0.6,1,500.,3000.,15.,25.)
#all_channels['t1_elminus_WJets_CR'] = Channel('t1_elminus_WJets_CR','eventTopology==1 && lepflavor==2 && lep_Q<0 && wjets_cr_selection==1',8,-1.,1.,3,0.,0.6,1,500.,3000.,25.,8.)
#all_channels['t2_muplus_WJets_CR']  = Channel('t2_muplus_WJets_CR','eventTopology==2 && lepflavor==1 && lep_Q>0 && wjets_cr_selection==1',14,-1.,1.,8,0.,0.4,3,300.,2300.,50.)
#all_channels['t2_muminus_WJets_CR'] = Channel('t2_muminus_WJets_CR','eventTopology==2 && lepflavor==1 && lep_Q<0 && wjets_cr_selection==1',14,-1.,1.,8,0.,0.4,3,300.,2300.,25.)
#all_channels['t2_elplus_WJets_CR']  = Channel('t2_elplus_WJets_CR','eventTopology==2 && lepflavor==2 && lep_Q>0 && wjets_cr_selection==1',14,-1.,1.,8,0.,0.4,3,300.,2300.,5.)
#all_channels['t2_elminus_WJets_CR'] = Channel('t2_elminus_WJets_CR','eventTopology==2 && lepflavor==2 && lep_Q<0 && wjets_cr_selection==1',14,-1.,1.,8,0.,0.4,3,300.,2300.,4.)
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
