#imports
from ROOT import TChain, TMinuit, Long, Double, TFile, TH3D
from array import array
from datetime import date
from math import *
import multiprocessing

#global variables
KEEP_X_AXIS_FIXED = False
total_ttree_dir_path = '/uscms_data/d3/eminizer/ttbar_13TeV/CMSSW_8_0_26_patch1/src/Analysis/Reconstructor/test/total_ttree_files/'
global_tree = None
global_cstar = array('f',[-1.0])
global_x_F   = array('f',[-1.0])
global_M 	 = array('f',[-1.0])
global_nXbins = 0; global_xMin = 0.; global_xMax = 0.
global_nYbins = 0; global_yMin = 0.; global_yMax = 0.
global_nZbins = 0; global_zMin = 0.; global_zMax = 0.

######			MINUIT MINIMIZATION FUNCTION			#####
def fcnt(npar, deriv, f, par, flag) :
	#take the absolute value of all the bin fractions
	for i in range(global_nXbins+global_nYbins+global_nZbins) :
		par[i] = abs(par[i])
	#first build the bin arrays
	#x bins
	xbins_i = array('d',(global_nXbins+1)*[0.])
	xbins_i[0] = global_xMin
	for i in range(global_nXbins) :
		xbins_i[i+1]=xbins_i[i]+(global_xMax-global_xMin)*par[i]
	#y bins
	ybins_i = array('d',(global_nYbins+1)*[0.])
	ybins_i[0] = global_yMin
	for i in range(global_nYbins) :
		ybins_i[i+1]=ybins_i[i]+(global_yMax-global_yMin)*par[(global_nXbins-1)+i]
	#z bins
	zbins_i = array('d',(global_nZbins+1)*[0.])
	zbins_i[0] = global_zMin
	for i in range(global_nZbins) :
		zbins_i[i+1]=zbins_i[i]+(global_zMax-global_zMin)*par[(global_nXbins+global_nYbins-2)+i]
	#then make a histogram to fill
	histo_i = TH3D('histo_for_fitting','histo_for_fitting',global_nXbins,xbins_i,global_nYbins,ybins_i,global_nZbins,zbins_i)
	#fill the histogram
	global_tree.SetBranchAddress('cstar',global_cstar)
	global_tree.SetBranchAddress('x_F',global_x_F)
	global_tree.SetBranchAddress('M',global_M) 
	for i in range(global_tree.GetEntries()) :
		global_tree.GetEntry(i)
		#print 'filling with cstar=%.2f, |x_F|=%.3f, M=%.2f'%(global_cstar[0],abs(global_x_F[0]),global_M[0])
		histo_i.Fill(global_cstar[0],abs(global_x_F[0]),global_M[0])
	#finally calculate the likelihood in pieces
	#norm is sum over all bins of number of events in the bin + s
	s = par[global_nXbins+global_nYbins+global_nZbins-3]
	norm = 0.
	nglobalbins = histo_i.GetSize()
	for k in range(nglobalbins) :
		if histo_i.IsBinOverflow(k) or histo_i.IsBinUnderflow(k) :
			continue
		norm+=histo_i.GetBinContent(k)+s
	#assemble the full likelihood
	L=0.
	for k in range(nglobalbins) :
		if histo_i.IsBinOverflow(k) or histo_i.IsBinUnderflow(k) :
			continue
		binx = array('i',[0]); biny = array('i',[0]); binz = array('i',[0])
		histo_i.GetBinXYZ(k,binx,biny,binz)
		xwidth = histo_i.GetXaxis().GetBinWidth(binx[0])/(global_xMax-global_xMin)
		ywidth = histo_i.GetYaxis().GetBinWidth(biny[0])/(global_yMax-global_yMin)
		zwidth = histo_i.GetZaxis().GetBinWidth(binz[0])/(global_zMax-global_zMin)
		cont = histo_i.GetBinContent(k)
		#print 'content in bin %d is %.2f'%(k,cont)
		if cont!=0 :
			num = cont+s-1.
			denom = (xwidth*ywidth*zwidth)*(norm-1.)
			if denom<=0. :
				print 'WARNING; DENOMINATOR IS %.4f (xw=%.2f, yw=%.2f, zw=%.2f, norm=%.2f)'%(denom,xwidth,ywidth,zwidth,norm)
				f[0]=100000000000000.
				return
			if (num/denom)<=0. :
				print 'WARNING; cont*num/denom<=0.! (cont=%d, num=%.2f, denom=%.2f)'%(cont,num,denom)
				f[0]=100000000000000.
				return
			L+=cont*log(num/denom)
	sumxweights = sum([par[i] for i in range(global_nXbins)])
	sumyweights = sum([par[global_nXbins+i] for i in range(global_nYbins)])
	sumzweights = sum([par[global_nXbins+global_nYbins+i] for i in range(global_nZbins)])
	f[0] = Double(0.)
	if L<=0. :
		print 'WARNING: L<=0. L=%.4f'%(L)
		f[0]+=10000.
	else :
		f[0] += -2.*log(L)
	f[0]+=100000.*(sumxweights-1.)**2+100000.*(sumyweights-1.)**2+100000.*(sumzweights-1.)**2+100.*(s-1.)**2
	print 'current f[0]=%.6f'%(f[0]) #DEBUG
	#print '	xbins=%s'%(xbins_i) #DEBUG
	#print '	ybins=%s'%(ybins_i) #DEBUG
	#print '	zbins=%s'%(zbins_i) #DEBUG

#channel class
class Channel(object) :

	#init
	def __init__(self,name,cutstring,nxbins,xmin,xmax,nybins,ymin,ymax,nzbins,zmin,zmax) :
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
		self._initial_x = array('d',(self._nxbins)*[(self._xmax-self._xmin)/self._nxbins])
		self._initial_y = array('d',(self._nybins)*[(self._ymax-self._ymin)/self._nybins])
		self._initial_z = array('d',(self._nzbins)*[(self._zmax-self._zmin)/self._nzbins])
		if self._name.find('t1_')!=-1 :
			self._initial_y = array('d',[0.25,0.25,0.5])
			self._initial_z = array('d',[0.25,0.1,0.15,0.5])
		self._data_chain = None
		self._skimmed_tree = None
		self._best_fit_x_bins = None
		self._best_fit_y_bins = None
		self._best_fit_z_bins = None

	######			PUBLIC FUNCTIONS			#####
	def skimDataChain(self) :
		self._skimmed_tree = self._data_chain.CopyTree(self._cutstring)

	def prepareFit(self) :
		#number of parameters in the fit will be equal to number of undetermined bin edges nxbins+nybins+nzbins plus one for the s smoothing parameter
		self._minimizer = TMinuit(self._nxbins+self._nybins+self._nzbins+1)
		ierflag = Long(1)
		#add parameters to the fitter
		#fractional bin widths
		for i in range(self._nxbins) : #x bin widths
			self._minimizer.mnparm(i,'x_'+str(i),self._initial_x[i],0.35*self._initial_x[i],0.,1.,ierflag)
		for i in range(self._nybins) : #y bin widths
			self._minimizer.mnparm((self._nxbins)+i,'y_'+str(i),self._initial_y[i],0.35*self._initial_y[i],0.,1.,ierflag)
		for i in range(self._nzbins) : #z bin widths
			self._minimizer.mnparm((self._nxbins+self._nybins)+i,'z_'+str(i),self._initial_z[i],0.35*self._initial_z[i],0.,1.,ierflag)
		#'s' smoothing parameter
		self._minimizer.mnparm((self._nxbins+self._nybins+self._nzbins),'s',1.,0.2,0.,0.,ierflag)
		#set the function to use
		self._minimizer.SetFCN(fcnt)

	def doFit(self) :
		print '--------------------------------------------------------------------------------'
		print 'starting fit for channel '+self._name+'...'
		#start by setting the global variables
		global global_tree
		global global_nXbins 
		global global_xMin
		global global_xMax
		global global_nYbins 
		global global_yMin
		global global_yMax
		global global_nZbins 
		global global_zMin
		global global_zMax
		global_tree = self._skimmed_tree
		global_nXbins = self._nxbins; global_xMin = self._xmin; global_xMax = self._xmax
		global_nYbins = self._nybins; global_yMin = self._ymin; global_yMax = self._ymax
		global_nZbins = self._nzbins; global_zMin = self._zmin; global_zMax = self._zmax
		#and some Minuit-specific stuff
		ierflag = Long(1)
		arglist = array( 'd', [100000.] )
		if KEEP_X_AXIS_FIXED : #fix the binning on the x-axis
			for i in range(self._nxbins) :
				print 'fixing parameter %d'%(i)
				self._minimizer.mnfixp(0,ierflag)
		#minimize and get the error flag
		self._minimizer.mnexcm('MIGRAD', arglist, 1, ierflag)
		errflag = int(ierflag)
		#Check fit Chi2 
		tmp1 = Double(1.0); tmp2 = Double(1.0); tmp3 = Double(1.0)
		self._minimizer.mnstat(tmp1,tmp2,tmp3,Long(1),Long(1),Long(1))
		chi2value = float(tmp1)
		print 'COMPLETED fit for channel %s with error flag %s and final chi^2=%.5f'%(self._name,errflag,chi2value)
		if KEEP_X_AXIS_FIXED : #unfix the x-axis binning
			self._minimizer.mnfree(0)
		#Get the bestfit parameters back from minuit
		#x bins
		self._best_fit_x_bins = array('d',(self._nxbins+1)*[0.])
		self._best_fit_x_bins[0]=self._xmin
		for i in range(self._nxbins) :
			tmp = Double(1.0); tmp2 = Double(1.0)
			self._minimizer.GetParameter(i,tmp,tmp2)
			self._best_fit_x_bins[i+1] = self._best_fit_x_bins[i]+(self._xmax-self._xmin)*abs(float(tmp))
		#y bins
		self._best_fit_y_bins = array('d',(self._nybins+1)*[0.])
		self._best_fit_y_bins[0]=self._ymin
		for i in range(self._nybins) :
			tmp = Double(1.0); tmp2 = Double(1.0)
			self._minimizer.GetParameter((self._nxbins)+i,tmp,tmp2)
			self._best_fit_y_bins[i+1] = self._best_fit_y_bins[i]+(self._ymax-self._ymin)*abs(float(tmp))
		#z bins
		self._best_fit_z_bins = array('d',(self._nzbins+1)*[0.])
		self._best_fit_z_bins[0]=self._zmin
		for i in range(self._nzbins) :
			tmp = Double(1.0); tmp2 = Double(1.0)
			self._minimizer.GetParameter((self._nxbins+self._nybins)+i,tmp,tmp2)
			self._best_fit_z_bins[i+1] = self._best_fit_z_bins[i]+(self._zmax-self._zmin)*abs(float(tmp))
		print 'FINAL BIN ARRAYS:'
		print 'xbins = %s'%(self._best_fit_x_bins)
		print 'ybins = %s'%(self._best_fit_y_bins)
		print 'zbins = %s'%(self._best_fit_z_bins)
		print '--------------------------------------------------------------------------------'
		#reset the global variables
		global_tree = None
		global_cstar[0]=-1; global_x_F[0]=-1; global_M[0]=-1
		global_nXbins = 0; global_xMin = 0.; global_xMax = 0.
		global_nYbins = 0; global_yMin = 0.; global_yMax = 0.
		global_nZbins = 0; global_zMin = 0.; global_zMax = 0.

	def getBestFitHistoList(self) :
		#declare 3D histogram
		histoname = '%s_best_fit_3D_histo'%(self._name)
		self._best_fit_histo = TH3D(histoname,'%s channel best fit 3D histogram; c*; |x_{F}|; M [GeV]'%(self._name),
									self._nxbins,self._best_fit_x_bins,self._nybins,self._best_fit_y_bins,self._nzbins,self._best_fit_z_bins)
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
		x_bins_string = "_XBINS = array('d',["
		for i in range(self._nxbins+1) :
			x_bins_string+=str(self._best_fit_x_bins[i])
			if i<self._nxbins-2 :
				x_bins_string+=','
			else :
				x_bins_string+='])'
		#y bins
		y_bins_string = "_YBINS = array('d',["
		for i in range(self._nybins+1) :
			y_bins_string+=str(self._best_fit_y_bins[i])
			if i<self._nybins-2 :
				y_bins_string+=','
			else :
				y_bins_string+='])'
		#z bins
		z_bins_string = "_ZBINS = array('d',["
		for i in range(self._nzbins+1) :
			z_bins_string+=str(self._best_fit_z_bins[i])
			if i<self._nzbins-2 :
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
#all_channels['t1_muplus_SR']  = Channel('t1_muplus_SR','eventTopology==1 && lepflavor==1 && lep_Q>0 && fullselection==1',6,-1.,1.,3,0.,0.6,4,500.,3000.)
#all_channels['t1_muminus_SR'] = Channel('t1_muminus_SR','eventTopology==1 && lepflavor==1 && lep_Q<0 && fullselection==1',20,-1.,1.,4,0.,0.6,5,500.,3000.)
#all_channels['t1_elplus_SR']  = Channel('t1_elplus_SR','eventTopology==1 && lepflavor==2 && lep_Q>0 && fullselection==1',20,-1.,1.,4,0.,0.6,5,500.,3000.)
#all_channels['t1_elminus_SR'] = Channel('t1_elminus_SR','eventTopology==1 && lepflavor==2 && lep_Q<0 && fullselection==1',20,-1.,1.,4,0.,0.6,5,500.,3000.)
#all_channels['t2_muplus_SR']  = Channel('t2_muplus_SR','eventTopology==2 && lepflavor==1 && lep_Q>0 && fullselection==1',20,-1.,1.,3,0.,0.4,6,300.,2300.)
#all_channels['t2_muminus_SR'] = Channel('t2_muminus_SR','eventTopology==2 && lepflavor==1 && lep_Q<0 && fullselection==1',20,-1.,1.,3,0.,0.4,6,300.,2300.)
#all_channels['t2_elplus_SR']  = Channel('t2_elplus_SR','eventTopology==2 && lepflavor==2 && lep_Q>0 && fullselection==1',20,-1.,1.,3,0.,0.4,6,300.,2300.)
#all_channels['t2_elminus_SR'] = Channel('t2_elminus_SR','eventTopology==2 && lepflavor==2 && lep_Q<0 && fullselection==1',20,-1.,1.,3,0.,0.4,6,300.,2300.)
all_channels['t3_muplus_SR']  = Channel('t3_muplus_SR','eventTopology==3 && lepflavor==1 && lep_Q>0 && fullselection==1',20,-1.,1.,3,0.,0.3,6,300.,1500.)
#all_channels['t3_muminus_SR'] = Channel('t3_muminus_SR','eventTopology==3 && lepflavor==1 && lep_Q<0 && fullselection==1',20,-1.,1.,3,0.,0.3,6,300.,1500.)
#all_channels['t3_elplus_SR']  = Channel('t3_elplus_SR','eventTopology==3 && lepflavor==2 && lep_Q>0 && fullselection==1',20,-1.,1.,3,0.,0.3,6,300.,1500.)
#all_channels['t3_elminus_SR'] = Channel('t3_elminus_SR','eventTopology==3 && lepflavor==2 && lep_Q<0 && fullselection==1',20,-1.,1.,3,0.,0.3,6,300.,1500.)
##boosted W+Jets control regions
#all_channels['t1_muplus_WJets_CR']  = Channel('t1_muplus_WJets_CR','eventTopology==1 && lepflavor==1 && lep_Q>0 && wjets_cr_selection==1',20,-1.,1.,4,0.,0.6,5,500.,3000.)
#all_channels['t1_muminus_WJets_CR'] = Channel('t1_muminus_WJets_CR','eventTopology==1 && lepflavor==1 && lep_Q<0 && wjets_cr_selection==1',20,-1.,1.,4,0.,0.6,5,500.,3000.)
#all_channels['t1_elplus_WJets_CR']  = Channel('t1_elplus_WJets_CR','eventTopology==1 && lepflavor==2 && lep_Q>0 && wjets_cr_selection==1',20,-1.,1.,4,0.,0.6,5,500.,3000.)
#all_channels['t1_elminus_WJets_CR'] = Channel('t1_elminus_WJets_CR','eventTopology==1 && lepflavor==2 && lep_Q<0 && wjets_cr_selection==1',20,-1.,1.,4,0.,0.6,5,500.,3000.)
#all_channels['t2_muplus_WJets_CR']  = Channel('t2_muplus_WJets_CR','eventTopology==2 && lepflavor==1 && lep_Q>0 && wjets_cr_selection==1',20,-1.,1.,3,0.,0.4,6,300.,2300.)
#all_channels['t2_muminus_WJets_CR'] = Channel('t2_muminus_WJets_CR','eventTopology==2 && lepflavor==1 && lep_Q<0 && wjets_cr_selection==1',20,-1.,1.,3,0.,0.4,6,300.,2300.)
#all_channels['t2_elplus_WJets_CR']  = Channel('t2_elplus_WJets_CR','eventTopology==2 && lepflavor==2 && lep_Q>0 && wjets_cr_selection==1',20,-1.,1.,3,0.,0.4,6,300.,2300.)
#all_channels['t2_elminus_WJets_CR'] = Channel('t2_elminus_WJets_CR','eventTopology==2 && lepflavor==2 && lep_Q<0 && wjets_cr_selection==1',20,-1.,1.,3,0.,0.4,6,300.,2300.)

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

#set up fits for each channel
print 'setting up fits...'
for channel in all_channels.values() :
	channel.prepareFit()
print 'done'

#do the fits for each channel
for channel in all_channels.values() :
	channel.doFit()

#write results to files
outrootfile = TFile('binning_test_histos_'+str(date.today())+'.root','recreate')
outtextfile = open('binning_test_codelines_'+str(date.today())+'.txt','w')
for channel in all_channels.values() :
	print 'saving results for channel %s'%(channel.getName())
	#save data histograms with the given binning
	histolist = channel.getBestFitHistoList()
	for histo in histolist :
		histo.Write()
	#save line of code to copy-past into template.py
	outtextfile.write('--------- binning for channel %s ---------'%(channel.getName()))
	codelinelist = channel.getCodeLineList()
	for codeline in codelinelist :
		outtextfile.write(codeline)
	outtextfile.write('------- end binning for channel %s -------\n'%(channel.getName()))
outrootfile.Close()
outtextfile.close()
print 'All done!!'
