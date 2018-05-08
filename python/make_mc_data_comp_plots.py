from ROOT import *
import CMS_lumi, tdrstyle
from math import *
from array import array
from optparse import OptionParser
from template import Template

gROOT.SetBatch()

#dictionary of the process names and their drawing colors
procs = ['fqcd','fwjets','fbck','fg0','fg1','fg2','fg3','fg4','fqm0','fqm1','fqm2','fqp0','fqp1','fqp2'] #this is ordered to always stack the MC histograms in the same way
prefit_procs = ['fqcd','fwjets','fbck','fgg','fqq']
proc_colors = {'fqcd':kYellow,
				'fwjets':kGreen,
				'fbck':kMagenta,
		 		'fgg':kBlue,
		 		'fg0':kBlue,
		 		'fg1':kBlue,
		 		'fg2':kBlue,
		 		'fg3':kBlue,
		 		'fg4':kBlue,
		 		'fqq':kRed+2,
		 		'fqm0':kRed+2,
		 		'fqm1':kRed+2,
		 		'fqm2':kRed+2,
		 		'fqp0':kRed+2,
		 		'fqp1':kRed+2,
		 		'fqp2':kRed+2}
proc_leg_names = {'fqcd':'QCD',
				  'fwjets':'W+Jets',
				  'fbck':'Single top/DY Jets',
				  'fg0':'gg #rightarrow t#bar{t}',
				  'fgg':'gg #rightarrow t#bar{t}',
				  'fqp0':'q#bar{q} #rightarrow t#bar{t}',
				  'fqq':'q#bar{q} #rightarrow t#bar{t}'}

parser = OptionParser()
parser.add_option('-M','--mode',  type='choice', action='store', dest='mode', choices=['prefit','postfit'])
parser.add_option('--tfilename',   type='string', action='store', default=None, dest='tfilename')
parser.add_option('--cfilename',   type='string', action='store', default=None, dest='cfilename')
parser.add_option('--outfilename', type='string', action='store', dest='outfilename')
(options, args) = parser.parse_args()

#Set some TDR options
tdrstyle.setTDRStyle()
iPeriod = 4 #13TeV iPeriod = 1*(0/1 7 TeV) + 2*(0/1 8 TeV)  + 4*(0/1 13 TeV)
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Preliminary"

#Plot class :
class Plot(object) :

	def __init__(self,channame,dim,lPos=3,iPos=11) :
		self._dim = dim #dimension (x,y,z)
		plotvars = {'x':'c*','y':'|x_{F}|','z':'M'}
		self._plotvar = plotvars[self._dim] #string of plot variable
		self._lPos = lPos
		self._iPos = iPos
		self._channame = channame #name of the channel this plot belongs to
		#pieces of the plot
		#MC stack
		self._MC_stack = THStack(channame+'_'+self._dim+'_stack',';;Events/bin')
		#make a new template to clone histogram shapes from
		dummy_template_name = channame+'__'+self._dim+'_dummy_template'
		dummy_template = Template(dummy_template_name,dummy_template_name,None)
		#clone histogram for MC error histograms, residual plots, and MC err residuals
		self._MC_err_histo = None; self._resid = None; self._MC_err_resid = None
		if self._dim=='x' :
			self._MC_err_histo = dummy_template.getHistoX().Clone()
			self._resid = dummy_template.getHistoX().Clone()
			self._MC_err_resid = dummy_template.getHistoX().Clone()
		elif self._dim=='y' :
			self._MC_err_histo = dummy_template.getHistoY().Clone()
			self._resid = dummy_template.getHistoY().Clone()
			self._MC_err_resid = dummy_template.getHistoY().Clone()
		elif self._dim=='z' :
			self._MC_err_histo = dummy_template.getHistoZ().Clone()
			self._resid = dummy_template.getHistoZ().Clone()
			self._MC_err_resid = dummy_template.getHistoZ().Clone()
		if self._MC_err_histo==None :
			print 'ERROR: could not get histogram cloned for plot dimension %s in channel %s'%(self._dim,self._channame)
		#Set attributes and directories
		#MC err histograms
		self._MC_err_histo.SetFillColor(kBlack); self._MC_err_histo.SetMarkerStyle(1); 
		self._MC_err_histo.SetFillStyle(3013); self._MC_err_histo.SetStats(0); self._MC_err_histo.SetDirectory(0)
		#residual plots
		self._resid.SetTitle(';'+self._plotvar+';Data/MC'); self._resid.SetStats(0); self._resid.SetMarkerStyle(20); self._resid.SetDirectory(0)
		#MC err residuals
		self._MC_err_resid.SetTitle(';'+self._plotvar+';Data/MC'); self._MC_err_resid.SetFillColor(kBlack)
		self._MC_err_resid.SetStats(0); self._MC_err_resid.SetMarkerStyle(1); self._MC_err_resid.SetFillStyle(3013)
		self._MC_err_resid.SetDirectory(0)
		#"oneline" to go at y=1. on the residuals plot
		self._oneline = TLine(self._resid.GetXaxis().GetBinLowEdge(1),1.,self._resid.GetXaxis().GetBinUpEdge(self._resid.GetNbinsX()),1.)
		self._oneline.SetLineWidth(3); self._oneline.SetLineStyle(2)
		#initialize the legend
		legwidth = 0.25
		legheight = 0.1+len(procs)*0.05
		x2 = 0.9; y2 = 0.83
		if self._lPos==1 :
			x2 = 0.44
		elif self._lPos==2 :
			x2 = 0.7
		self._leg = TLegend(x2-legwidth,y2-legheight,x2,y2)
		self._lumi_obj = None

	def addMChisto(self,h,pname) :
		#add the histogram to the stack
		print '		Adding histogram %s to stack for channel %s and process %s (size=%d)'%(h,self._channame,pname,h.GetSize())
		self._MC_stack.Add(h,'hist')
		#increment MC error values
		for j in range(h.GetSize()) :
			if not h.IsBinOverflow(j) and not h.IsBinUnderflow(j) :
				oldcont = self._MC_err_histo.GetBinContent(j)
				newcont = oldcont+h.GetBinContent(j)
				self._MC_err_histo.SetBinContent(j,newcont)
				olderr = self._MC_err_histo.GetBinError(j)
				newerr = olderr+h.GetBinError(j)**2+h.GetBinContent(j)
				self._MC_err_histo.SetBinError(j,newerr)
#				if j==5 : #DEBUG
#					print '			bin %d content (%.2f -> %.2f) and error (%.2f -> %.2f)'%(j,oldcont,newcont,olderr,newerr) #DEBUG
		#add a Legend entry
		if pname in proc_leg_names.keys() :
			self._leg.AddEntry(h,proc_leg_names[pname],'F')

	def setMCerrors(self) :
		#in every bin take the sqrt of the MC error since the initial values were sums of err**2
		print '		Resetting MC stack errors in channel %s (size=%d)'%(self._channame,self._MC_err_histo.GetSize())
		for j in range(self._MC_err_histo.GetSize()) :
			if not self._MC_err_histo.IsBinOverflow(j) and not self._MC_err_histo.IsBinUnderflow(j) :
				self._MC_err_histo.SetBinError(j,sqrt(self._MC_err_histo.GetBinError(j)))
				self._MC_err_resid.SetBinContent(j,1.)
				#also set the fractional error for the MC error residual histogram
				mccont = self._MC_stack.GetStack().Last().GetBinContent(j)
				if mccont!=0. :
					self._MC_err_resid.SetBinError(j,self._MC_err_histo.GetBinError(j)/mccont)

	def setDataHisto(self,h) :
		self._data_histo = h
		self._data_histo.SetMarkerStyle(20)
		self._leg.AddEntry(self._data_histo,'data','PE')
#		print '		setting data histo in channel %s to %s'%(self._channame,self._data_histo) #DEBUG

	def calculateResiduals(self) :
		#loop over bins
#		print '		calculating residuals in channel %s using data histo %s with size=%d'%(self._channame,self._data_histo,self._data_histo.GetSize()) #DEBUG
#		print '		data histo integral = %.2f, MC stack integral = %.2f'%(self._data_histo.Integral(),sum(h.Integral() for h in self._MC_stack.GetHists())) #DEBUG
		for j in range(self._data_histo.GetSize()) :
			if self._data_histo.IsBinUnderflow(j) or self._data_histo.IsBinOverflow(j) :
				continue
			if self._MC_err_histo.GetBinContent(j)!=0. and self._data_histo.GetBinContent(j)!=0. :
				resid_cont = 0.; resid_err = 0;
				datacont = self._data_histo.GetBinContent(j)
				MC_cont  = self._MC_stack.GetStack().Last().GetBinContent(j)
				if MC_cont!=0. :
					resid_cont=datacont/abs(MC_cont)
					resid_err = sqrt(datacont)/abs(MC_cont)
				self._resid.SetBinContent(j,resid_cont)
				self._resid.SetBinError(j,resid_err)
#				if j==5 : #DEBUG
#					print '			bin %d residual = %.2f +/- %.2f'%(j,resid_cont,resid_err) #DEBUG

	def plotOnCanvas(self) :
		outfile.cd()
		#declare the canvas
		cname = self._channame+'_'+self._dim+'_canv'
		self._canv = TCanvas(cname,cname,1100,900)
		self._canv.cd()
		#reset the maximum of the stack to show the whole thing
		self._MC_stack.SetMaximum(1.10*max(self._MC_stack.GetMaximum()+sqrt(self._MC_stack.GetMaximum()),self._data_histo.GetMaximum()+sqrt(self._data_histo.GetMaximum())))
		#add the MC uncertainty to the legend
		self._leg.AddEntry(self._MC_err_histo,'MC uncertainty','F')
		#build the final plot
		self._canv.cd()
		gStyle.SetTitleFontSize(0.00001)
		padname = self._channame+'_'+self._dim
		self._histo_pad = TPad(padname+'_histo_pad',padname+'_histo_pad',0,0.25,1,1)
		self._resid_pad = TPad(padname+'_resid_pad',padname+'_resid_pad',0,0,1.,0.25)
		self._histo_pad.SetLeftMargin(0.16); self._histo_pad.SetRightMargin(0.05) 
		self._histo_pad.SetTopMargin(0.11);	 self._histo_pad.SetBottomMargin(0.02)
		self._histo_pad.SetBorderMode(0)
		self._resid_pad.SetLeftMargin(0.16); self._resid_pad.SetRightMargin(0.05)
		self._resid_pad.SetTopMargin(0.0);   self._resid_pad.SetBottomMargin(0.42)
		self._resid_pad.SetBorderMode(0)
		self._histo_pad.Draw(); self._resid_pad.Draw()
		#plot on the histogram pad
		self._histo_pad.cd()
		self._MC_stack.Draw(); self._MC_err_histo.Draw("SAME E2"); self._data_histo.Draw("SAME PE")
		if self._lPos!=0 :
			self._leg.Draw("SAME")
		#make the channel identifier text
		chanTxt = 'type '+self._channame.split('_')[0].split('t')[1]
		if self._channame.split('_')[1].startswith('mu') :
			chanTxt+=' #mu + jets'
		elif self._channame.split('_')[1].startswith('el') :
			chanTxt+=' e + jets'
		if self._channame.split('_')[1].endswith('plus') :
			chanTxt+=' (Q>0'
		elif self._channame.split('_')[1].endswith('minus') :
			chanTxt+=' (Q<0'
		if self._channame.find('SR')!=-1 :
			chanTxt+=' SR)'
		elif self._channame.find('WJets_CR')!=-1 :
			chanTxt+=' WJets CR)'
		chantxt = TLatex()
		chantxt.SetNDC()
		chantxt.SetTextAngle(0)
		chantxt.SetTextColor(kBlack)
		chantxt.SetTextFont(42)
		chantxt.SetTextAlign(11) 
		chantxt.SetTextSize(0.6*0.11)
		chantxt.DrawLatex(0.16,1-0.11+0.2*0.11,chanTxt)
		self._MC_stack.GetHistogram().GetXaxis().SetLabelSize(0.)
		#plot on the residuals pad
		self._resid_pad.cd()
		self._MC_err_resid.Draw("E2"); self._oneline.Draw("SAME"); self._resid.Draw("SAME PE")
		self._MC_err_resid.GetXaxis().SetLabelSize(0.15)
		self._MC_err_resid.GetYaxis().SetLabelSize(0.15)
		self._MC_err_resid.GetYaxis().SetTitleOffset(0.25)
		self._MC_err_resid.GetXaxis().SetTitleSize((0.75/0.25)*self._MC_err_resid.GetXaxis().GetTitleSize())
		self._MC_err_resid.GetYaxis().SetTitleSize((0.75/0.25)*self._MC_err_resid.GetYaxis().GetTitleSize())
		self._MC_err_resid.GetYaxis().SetRangeUser(0.7,1.3)
		self._MC_err_resid.GetYaxis().SetNdivisions(504)
		self._canv.Update()
		#plot the CMS_Lumi lines on the canvases
		self._lumi_obj=CMS_lumi.CMS_lumi(self._histo_pad, iPeriod, self._iPos)
		return self._canv


#Process class
class Process(object) :

	def __init__(self,pname,histo1D) :
		self._name = pname
		self._color = proc_colors[self._name]
		self._histo1D = histo1D
		self._histo1D.SetDirectory(0)
		print '			new process with name %s and histo1D %s'%(self._name, self._histo1D)

	def getHistoProjections(self,channame) :
		#make a new template
		newname = channame+'__'+self._name
		self._newtemp = Template(newname+'__POSTFIT',newname+'__POSTFIT',None)
		#make its histograms from the 1D histo from the file
		self._newtemp.make_from_1D_histo(self._histo1D)
		#get the three projection histograms
		self._x_histo = self._newtemp.getHistoX()
		self._y_histo = self._newtemp.getHistoY()
		self._z_histo = self._newtemp.getHistoZ()
		threehistos = [self._x_histo,self._y_histo,self._z_histo]
		#set attributes
		for h in threehistos :
			h.SetFillColor(self._color); h.SetLineColor(self._color); h.SetMarkerStyle(21); h.SetMarkerColor(self._color); 
			h.GetXaxis().SetLabelSize(0.); h.SetDirectory(0)
		#return histogram projections
		return threehistos

	############### Getters/Setters ###############
	def getName(self) :
		return self._name
	def getHisto1D(self) :
		return self._histo1D

#plotGroup class
class PlotGroup(object) :

	def __init__(self,channame) :
		self._channame = channame #channel name (one group per channel)
		self._processes = [] #list of contributing processes

	def addProcess(self,pname,histo1D) :
		self._processes.append(Process(pname,histo1D))

	def initializePlots(self) :
		#make new plots for each dimension
		self._x_plot = Plot(self._channame,'x',2)
		self._y_plot = Plot(self._channame,'y')
		self._z_plot = Plot(self._channame,'z')

	def addMCHistograms(self) :
		#for every process
		for proc in self._processes :
			print '		Getting projection histograms for process %s'%(proc.getName())
			#get the three 1-D projection histograms
			x_histo, y_histo, z_histo = proc.getHistoProjections(self._channame)
#			print '		Returned histograms:' #DEBUG
#			print '			x: %s'%(x_histo) #DEBUG
#			print '			y: %s'%(y_histo) #DEBUG
#			print '			z: %s'%(z_histo) #DEBUG
			#add them to the plots
			self._x_plot.addMChisto(x_histo,proc.getName())
			self._y_plot.addMChisto(y_histo,proc.getName())
			self._z_plot.addMChisto(z_histo,proc.getName())
		#after all the processes have been added, take the sqrt of the bin errors because I was incrementing by the err**2
		self._x_plot.setMCerrors()
		self._y_plot.setMCerrors()
		self._z_plot.setMCerrors()

	def buildResiduals(self) :
		#first get the 1D input data graph from the initial template file and copy it into a 1D histogram like in the MC processes
		newname = self._channame+'__data_obs'
		newdatatemp = Template(newname+'__POSTFIT',newname+'__POSTFIT',None)
		data1Dhisto=None
		if options.mode=='postfit' :
			combine_filep = TFile.Open(options.cfilename)
			data_graph = combine_filep.Get('shapes_fit_s/%s/data'%(self._channame))
			data1Dhisto = self._processes[0].getHisto1D().Clone()
			data1Dhisto.SetDirectory(0)
			data1Dhisto.Reset()
#			print '		ngraphpoints=%d data1Dhisto: nbins=%d, integral before filling=%.2f'%(data_graph.GetN(),data1Dhisto.GetSize()-2,data1Dhisto.Integral()) #DEBUG
			for i in range(data_graph.GetN()) :
				px = array('d',[0.]); py = array('d',[0.])
				data_graph.GetPoint(i,px,py)
				err = data_graph.GetErrorY(i)
#				if i==350 or i==475 : #DEBUG
#					print '		bin %d in data has px=%d, py=%.4f, err=%.4f'%(i,px[0],py[0],err) #DEBUG
				data1Dhisto.Fill(px[0],py[0])
				data1Dhisto.SetBinError(data1Dhisto.FindFixBin(px[0]),err)
			combine_filep.Close()
		elif options.mode=='prefit' :
			initial_templates_file = TFile(options.tfilename)
			data1Dhisto = initial_templates_file.Get(newname).Clone()
			data1Dhisto.SetDirectory(0)
			initial_templates_file.Close()
		newdatatemp.make_from_1D_histo(data1Dhisto)
		self._data_histo_x=newdatatemp.getHistoX()
		self._data_histo_y=newdatatemp.getHistoY()
		self._data_histo_z=newdatatemp.getHistoZ()
#		print '		plot group data histos: 1D=%s, x=%s, y=%s, z=%s (integrals: 1D=%.2f, x=%.2f, y=%.2f, z=%.2f)'%(data1Dhisto,self._data_histo_x,self._data_histo_y,self._data_histo_z,data1Dhisto.Integral(),self._data_histo_x.Integral(),self._data_histo_y.Integral(),self._data_histo_z.Integral()) #DEBUG
		#add them to the plot
		self._x_plot.setDataHisto(self._data_histo_x)
		self._y_plot.setDataHisto(self._data_histo_y)
		self._z_plot.setDataHisto(self._data_histo_z)
		#calculate residuals
		self._x_plot.calculateResiduals()
		self._y_plot.calculateResiduals()
		self._z_plot.calculateResiduals()

	def makePlots(self) :
		self._x_canv = self._x_plot.plotOnCanvas()
		self._y_canv = self._y_plot.plotOnCanvas()
		self._z_canv = self._z_plot.plotOnCanvas()
		outfile.cd()
		self._x_canv.Write()
		self._y_canv.Write()
		self._z_canv.Write()

	############### Getters/Setters ###############
	def getName(self) :
		return self._channame

############################################################################################
#################### 					Main Script 					####################
############################################################################################

#dict of all PlotGroups (keys = channel names, values = PlotGroup object)
all_plot_groups = {}

if options.mode=='postfit' :
	#initialize channels based on reading the input file from Combine
	#Open the input file with the postfit shapes from Combine
	combine_filep = TFile.Open(options.cfilename)
	#cd to the directory with fit results
	gDirectory.cd('shapes_fit_s')
	#each key in this directory's name is a channel name
	print 'Setting up plot groups...'
	for k1 in gDirectory.GetListOfKeys() :
		channel_name = k1.GetName()
		print '	Adding channel with name %s'%(channel_name)
		##skip the control region
		#if channel_name.find('WJets_CR')!=-1 :
		#	continue
		#initialize a plot group for this channel
		all_plot_groups[channel_name] = PlotGroup(channel_name)
		#cd to this channel's subdirectory
		combine_filep.cd('shapes_fit_s/%s'%(channel_name))
		#find the keys in this directory that are process names
		for procname in procs :
			for k2 in gDirectory.GetListOfKeys() :
				if not k2.GetName()==procname :
					continue
				#add processes to the PlotGroup
				print '		Adding process with name %s'%(k2.GetName())
				all_plot_groups[channel_name].addProcess(k2.GetName(),gDirectory.Get(k2.GetName()))
	combine_filep.Close()
elif options.mode=='prefit' :
	#initialize channels hardcoded
	all_channel_names = ['t1_muplus_SR','t1_muminus_SR','t1_elplus_SR','t1_elminus_SR',
						 't1_muplus_WJets_CR','t1_muminus_WJets_CR','t1_elplus_WJets_CR','t1_elminus_WJets_CR',
						 't2_muplus_SR','t2_muminus_SR','t2_elplus_SR','t2_elminus_SR',
						 't2_muplus_WJets_CR','t2_muminus_WJets_CR','t2_elplus_WJets_CR','t2_elminus_WJets_CR',
						 't3_muplus_SR','t3_muminus_SR','t3_elplus_SR','t3_elminus_SR',]
	#open the template file to get the 1D histograms
	template_filep = TFile.Open(options.tfilename)
	print 'Setting up plot groups...'
	#loop over the channel names
	for channel_name in all_channel_names :
		print '	Adding channel with name %s'%(channel_name)
		#initialize a plot group for this channel
		all_plot_groups[channel_name] = PlotGroup(channel_name)
		for procname in prefit_procs :
			if (channel_name.startswith('t1') or (channel_name.startswith('t2') and channel_name.find('SR')!=-1)) and procname=='fqcd' :
				continue
			#add processes to the PlotGroup
			print '		Adding process with name %s'%(procname)
			all_plot_groups[channel_name].addProcess(procname,template_filep.Get(channel_name+'__'+procname))
	template_filep.Close()

#start the output file
outfile = TFile(options.outfilename,'recreate')

#Initialize the plot objects in each channel
print 'Initializing plot objects...'
for pg in all_plot_groups.values() :
	print '	in channel %s'%(pg.getName())
	pg.initializePlots()
print 'Done'

#build 3D histograms out of 1D post-fit histograms from combine; add to histogram stacks and MC uncertainty graphs
print 'Building 3D histograms from 1D postfit histograms...'
for pg in all_plot_groups.values() :
	print '	in channel %s'%(pg.getName())
	pg.addMCHistograms()
print 'Done'

#build residuals plots
print 'Building residual plots...'
for pg in all_plot_groups.values() :
	print '	in channel %s'%(pg.getName())
	pg.buildResiduals()
print 'Done'

#plot and save all plots
print 'Plotting on canvases...'
for pg in all_plot_groups.values() :
	print '	in channel %s'%(pg.getName())
	pg.makePlots()
print 'Done'

outfile.Close()