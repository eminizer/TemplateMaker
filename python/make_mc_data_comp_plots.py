from ROOT import *
import os
import CMS_lumi, tdrstyle
from math import *
from array import array
from optparse import OptionParser
from template import Template

gROOT.SetBatch()

#dictionary of the process names and their drawing colors
procs = ['fqcd','fwjets','fbck','fg0','fg1','fg2','fg3','fg4','fq0','fqp0','fqm0','fq1','fq2'] #this is ordered to always stack the MC histograms in the same way
prefit_procs = ['fqcd','fwjets','fbck','fgg','fqq']
#prefit_procs = ['fqcd','fwjets','fbck','fg0','fq0']
proc_colors = {'fqcd':kYellow,
				'fwjets':kGreen,
				'fbck':kViolet-9,
		 		'fgg':kBlue-6,
		 		'fg0':kBlue-6,
		 		'fg1':kBlue-6,
		 		'fg2':kBlue-6,
		 		'fg3':kBlue-6,
		 		'fg4':kBlue-6,
		 		'fqq':kRed-9,
		 		'fq0':kRed-9,
		 		'fqp0':kRed-9,
		 		'fqm0':kRed-9,
		 		'fq1':kRed-9,
		 		'fq2':kRed-9,
				}
proc_leg_names = {'fqcd':'Multijet',
				  'fwjets':'W+jets',
				  'fbck':'Other t quark & Z/#gamma+jets',
				  'fg0':'gg #rightarrow t#bar{t}',
				  'fgg':'gg #rightarrow t#bar{t}',
				  'fq0':'q#bar{q} #rightarrow t#bar{t}',
				  'fqp0':'q#bar{q} #rightarrow t#bar{t}',
				  'fqq':'q#bar{q} #rightarrow t#bar{t}',
				}
x_axis_divs = {'t1_mu_SR':(12,True),
			   't1_el_SR':(10,True),
			   't1_mu_CR':(10,True),
			   't1_el_CR':(5,True),
			   't2_mu_SR':(510,False),
			   't2_el_SR':(10,True),
			   't2_mu_CR':(10,True),
			   't2_el_CR':(1,True),
			   't3_mu_SR':(1023,False),
			   't3_el_SR':(1017,False),
			}
ymax_vals = {'t1_mu_SR':1300.,
			 't1_el_SR':925.,
			 't1_mu_CR':2000.,
			 't1_el_CR':2000.,
			 't2_mu_SR':4100.,
			 't2_el_SR':1300.,
			 't2_mu_CR':15000.,
			 't2_el_CR':14000.,
			 't3_mu_SR':5600.,
			 't3_el_SR':4500.,
			 't1_muplus_SR':700.,
			 't1_elplus_SR':500.,
			 't1_muplus_CR':1300.,
			 't1_elplus_CR':925.,
			 't2_muplus_SR':2300.,
			 't2_elplus_SR':700.,
			 't2_muplus_CR':8000.,
			 't2_elplus_CR':7000.,
			 't3_muplus_SR':2750.,
			 't3_elplus_SR':2350.,
			 't1_muminus_SR':700.,
			 't1_elminus_SR':500.,
			 't1_muminus_CR':1000.,
			 't1_elminus_CR':925.,
			 't2_muminus_SR':2300.,
			 't2_elminus_SR':650.,
			 't2_muminus_CR':8000.,
			 't2_elminus_CR':6250.,
			 't3_muminus_SR':2750.,
			 't3_elminus_SR':2300.,
			}
x_F_line_bins = {'t1_mu_SR':[6.,10.],
				 't1_el_SR':[4.,6.],
				 't1_mu_CR':[],
				 't1_el_CR':[],
				 't2_mu_SR':[18.,30.,36.,42.],
				 't2_el_SR':[4.],
				 't2_mu_CR':[],
				 't2_el_CR':[],
				 't3_mu_SR':[60.,110.,142.,174.,190.,198.,206.,212.,218.],
				 't3_el_SR':[50.,90.,122.,146.,154.,160.],
			}
M_line_bins = {'t1_mu_SR':[8.],
			   't1_el_SR':[],
			   't1_mu_CR':[],
			   't1_el_CR':[],
			   't2_mu_SR':[6.,12.,24.,38.,40.,44.],
			   't2_el_SR':[],
			   't2_mu_CR':[],
			   't2_el_CR':[],
			   't3_mu_SR':[20.,40.,70.,80.,90.,100.,126.,158.,182.,194.,202.],
			   't3_el_SR':[10.,20.,30.,40.,70.,106.,128.,134.,140.,150.],
			}
cstar_mid_bins = {'t1_mu_SR':[3.,7.,9.,11.],
				  't1_el_SR':[2.,5.,7.],
				  't1_mu_CR':[],
				  't1_el_CR':[],
				  't2_mu_SR':[3.,9.,15.,21.,27.,33.,37.,39.,41.,43.,45.],
				  't2_el_SR':[2.,5.],
				  't2_mu_CR':[],
				  't2_el_CR':[],
				  't3_mu_SR':[10.,30.,50.,65.,75.,85.,95.,105.,118.,134.,150.,166.,178.,186.,192.,196.,200.,204.,209.,215.,221.],
				  't3_el_SR':[5.,15.,25.,35.,45.,60.,80.,98.,114.,125.,131.,137.,143.,148.,152.,157.,162.],
			}

parser = OptionParser()
parser.add_option('-M','--mode',  type='choice', action='store', dest='mode', choices=['prefit','postfit'])
parser.add_option('--tfilename',   type='string', action='store', default=None, dest='tfilename')
parser.add_option('--cfilename',   type='string', action='store', default=None, dest='cfilename')
parser.add_option('--outfilename', type='string', action='store', dest='outfilename')
parser.add_option('--use_top_pt_reweighting', type='string', action='store', default='no', dest='usetopptrw')
parser.add_option('--sumCharges',  action='store_true', dest='sumCharges')
(options, args) = parser.parse_args()

#Set some TDR options
tdrstyle.setTDRStyle()
iPeriod = 4 #13TeV iPeriod = 1*(0/1 7 TeV) + 2*(0/1 8 TeV)  + 4*(0/1 13 TeV)
#CMS_lumi.writeExtraText = 1
CMS_lumi.writeExtraText = False
CMS_lumi.extraText = "Preliminary"

#Plot class :
class Plot(object) :

	def __init__(self,channame,dim,lPos=3,iPos=11) :
		self._dim = dim #dimension (x,y,z)
		plotvars = {'x':'c*','y':'|x_{F}|','z':'M','1D':'Template bin'}
		self._plotvar = plotvars[self._dim] #string of plot variable
		self._lPos = lPos
		self._iPos = iPos
		self._channame = channame #name of the channel this plot belongs to
		self._canvsize = (1100,900)
		if self._channame.startswith('t3') :
			self._canvsize = (5500,900)
		#pieces of the plot
		#MC stack
		self._MC_stack = THStack(channame+'_'+self._dim+'_stack',';;Events/bin')
		#list of MC histograms for legend
		self._MC_histos = []
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
		elif self._dim=='1D' :
			self._MC_err_histo = dummy_template.convertTo1D().Clone()
			self._resid = dummy_template.convertTo1D().Clone()
			self._MC_err_resid = dummy_template.convertTo1D().Clone()
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
		legwidth = 0.622
		legheight = 0.272
		if self._channame.startswith('t3') :
			legwidth=0.561
			legheight=0.183
		x2 = 0.922; y2 = 0.857
		#if self._lPos==1 :
		#	x2 = 0.44
		#elif self._lPos==2 :
		#	x2 = 0.7
		self._leg = TLegend(x2-legwidth,y2-legheight,x2,y2)
		if self._channame.startswith('t3') :
			self._leg.SetNColumns(3)
		else :
			self._leg.SetNColumns(2)
		self._lumi_obj = None
		self._pnames_added_to_leg = []
		#make the channel identifier text (and the key for the plot info)
		self._chanTxt = 'Type-'+self._channame.split('_')[0].split('t')[1]
		self._ckey='t'+self._channame.split('_')[0].split('t')[1]+'_'
		if self._channame.split('_')[1].startswith('mu') :
			self._chanTxt+=' #mu+jets'
			self._ckey+='mu'
		elif self._channame.split('_')[1].startswith('el') :
			self._chanTxt+=' e+jets'
			self._ckey+='el'
		self._chanTxt+=' ('
		if self._channame.split('_')[1].endswith('plus') :
			self._chanTxt+='Q>0 '
			self._ckey+='plus'
		elif self._channame.split('_')[1].endswith('minus') :
			self._chanTxt+='Q<0 '
			self._ckey+='minus'
		if self._channame.find('SR')!=-1 :
			self._chanTxt+='SR)'
			self._ckey+='_SR'
		elif self._channame.find('WJets_CR')!=-1 :
			self._chanTxt+='W+jets CR)'
			self._ckey+='_CR'

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
		#add to the list of MC histos for the legend
		if pname in proc_leg_names.keys() and pname not in self._pnames_added_to_leg :
			self._MC_histos.append((h,pname))
			self._pnames_added_to_leg.append(pname)

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
		#add data legend entry
		self._leg.AddEntry(self._data_histo,'Data','PE')
		#add MC legend entries
		for i in reversed(range(len(self._MC_histos))) :
			self._leg.AddEntry(self._MC_histos[i][0],proc_leg_names[self._MC_histos[i][1]],'F')
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
		self._canv = TCanvas(cname,cname,self._canvsize[0],self._canvsize[1])
		self._canv.cd()
		#reset the maximum of the stack to show the whole thing
		maxfac = 1.50 if self._channame.startswith('t3') else 1.30
		self._MC_stack.SetMaximum(1.20*max(self._MC_stack.GetMaximum()+sqrt(self._MC_stack.GetMaximum()),self._data_histo.GetMaximum()+sqrt(self._data_histo.GetMaximum())))
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
		#make the lists of bin division lines
		self._x_F_histo_lines = []; self._M_histo_lines = []; self._cstar_histo_lines = []
		self._x_F_resid_lines = []; self._M_resid_lines = []; self._cstar_resid_lines = []
		for binnum in x_F_line_bins[self._ckey.replace('plus','').replace('minus','')] :
			newhistoline = TLine(binnum,0.,binnum,(3./20.)*ymax_vals[self._ckey])
			newresidline = TLine(binnum,0.7,binnum,1.3)
			for l in [newhistoline,newresidline] :
				l.SetLineWidth(2) if self._ckey.startswith('t3') else l.SetLineWidth(5)
				#l.SetLineStyle(2)
				l.SetLineColor(kRed+1)
			self._x_F_histo_lines.append(newhistoline)
			self._x_F_resid_lines.append(newresidline)
		if len(self._x_F_histo_lines)>0 :
			self._leg.AddEntry(self._x_F_histo_lines[0],'|x_{r}| bin edges','L')
		for binnum in M_line_bins[self._ckey.replace('plus','').replace('minus','')] :
			newhistoline = TLine(binnum,0.,binnum,(2./20.)*ymax_vals[self._ckey])
			newresidline = TLine(binnum,0.75,binnum,1.25)
			for l in [newhistoline,newresidline] :
				l.SetLineWidth(1) if self._ckey.startswith('t3') else l.SetLineWidth(4)
				l.SetLineStyle(2)
				l.SetLineColor(kBlue+1)
			self._M_histo_lines.append(newhistoline)
			self._M_resid_lines.append(newresidline)
		if len(self._M_histo_lines)>0 :
			self._leg.AddEntry(self._M_histo_lines[0],'m_{r} bin edges','L')
		for binnum in cstar_mid_bins[self._ckey.replace('plus','').replace('minus','')] :
			newhistoline = TLine(binnum,0.,binnum,(1./20.)*ymax_vals[self._ckey])
			newresidline = TLine(binnum,0.8,binnum,1.2)
			for l in [newhistoline,newresidline] :
				l.SetLineWidth(1) if self._ckey.startswith('t3') else l.SetLineWidth(4)
				l.SetLineStyle(5)
				l.SetLineColor(kMagenta+1)
			self._cstar_histo_lines.append(newhistoline)
			self._cstar_resid_lines.append(newresidline)
		if len(self._cstar_histo_lines)>0 :
			self._leg.AddEntry(self._cstar_histo_lines[0],'c_{r}*=0 points','L')
		#plot on the histogram pad
		self._histo_pad.cd()
		self._MC_stack.Draw()
		self._MC_stack.GetHistogram().GetXaxis().SetLabelSize(0.)
		self._MC_stack.GetHistogram().GetXaxis().SetNdivisions(x_axis_divs[self._ckey.replace('plus','').replace('minus','')][0])
		if x_axis_divs[self._ckey.replace('plus','').replace('minus','')][1] :
			self._MC_stack.GetHistogram().GetXaxis().CenterLabels()
		self._MC_stack.SetMaximum(ymax_vals[self._ckey])
		self._canv.Update()
		self._MC_err_histo.Draw("SAME E2"); self._data_histo.Draw("SAME PE")
		if self._lPos!=0 :
			self._leg.Draw("SAME")
		self._MC_stack.GetYaxis().SetTitleOffset(1.10)
		if self._channame.startswith('t3') :
			self._MC_stack.GetYaxis().SetTitleOffset(0.6)
			self._MC_stack.GetYaxis().SetTickLength(0.02)
		#plot the bin division lines
		for l in (self._x_F_histo_lines+self._M_histo_lines+self._cstar_histo_lines) :
			l.Draw('SAME')
		#plot the channel identifier text
		self._chantxt = TLatex()
		self._chantxt.SetNDC()
		self._chantxt.SetTextAngle(0)
		self._chantxt.SetTextColor(kBlack)
		self._chantxt.SetTextFont(42)
		self._chantxt.SetTextAlign(11) 
		self._chantxt.SetTextSize(0.6*0.11)
		self._chantxt.DrawLatex(0.16,1-0.11+0.2*0.11,self._chanTxt)
		#plot on the residuals pad
		self._resid_pad.cd()
		self._MC_err_resid.Draw("E2"); self._oneline.Draw("SAME"); self._resid.Draw("SAME PEX0")
		self._MC_err_resid.GetXaxis().SetLabelSize(0.15)
		self._MC_err_resid.GetYaxis().SetLabelSize(0.15)
		self._MC_err_resid.GetYaxis().SetTitleOffset(0.25)
		self._MC_err_resid.GetXaxis().SetTitleSize((0.75/0.25)*self._MC_err_resid.GetXaxis().GetTitleSize())
		self._MC_err_resid.GetYaxis().SetTitleSize((0.75/0.25)*self._MC_err_resid.GetYaxis().GetTitleSize())
		self._MC_err_resid.GetYaxis().SetRangeUser(0.7,1.3)
		self._MC_err_resid.GetYaxis().SetNdivisions(504)
		self._MC_err_resid.GetXaxis().SetTickLength(0.09)
		self._MC_err_resid.GetXaxis().SetNdivisions(x_axis_divs[self._ckey.replace('plus','').replace('minus','')][0])
		if x_axis_divs[self._ckey.replace('plus','').replace('minus','')][1] :
			self._MC_err_resid.GetXaxis().CenterLabels()
		if self._channame.startswith('t3') :
			self._MC_err_resid.GetYaxis().SetTitleOffset(0.20)
			self._MC_err_resid.GetYaxis().SetTickLength(0.02)
		#plot the bin division lines
		for l in (self._x_F_resid_lines+self._M_resid_lines+self._cstar_resid_lines) :
			l.Draw('SAME')
		self._canv.Update()
		#plot the CMS_Lumi lines on the canvases
		self._lumi_obj=CMS_lumi.CMS_lumi(self._histo_pad, iPeriod, self._iPos)
		##write the canvas as a .pdf
		self._canv.SaveAs('.pdf')
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
		self._histo1D.SetFillColor(self._color); self._histo1D.SetLineColor(self._color); self._histo1D.SetMarkerStyle(21); self._histo1D.SetMarkerColor(self._color); 
		self._histo1D.GetXaxis().SetLabelSize(0.); self._histo1D.SetDirectory(0)
		return self._histo1D

#plotGroup class
class PlotGroup(object) :

	def __init__(self,channame) :
		self._channame = channame #channel name (one group per channel)
		self._processes = [] #list of contributing processes

	def addProcess(self,pname,histo1D) :
		if options.mode=='prefit' :
			another_dummy_template_name = self._channame+'__fdummy'+'_dummy_template'
			another_dummy_template = Template(another_dummy_template_name,another_dummy_template_name,None)
			another_dummy_template.make_from_1D_histo(histo1D)
			histo1D = another_dummy_template.convertTo1D()
		self._processes.append(Process(pname,histo1D))

	def initializePlots(self) :
		#make new plots for each dimension
		#self._x_plot = Plot(self._channame,'x',2)
		#self._y_plot = Plot(self._channame,'y')
		#self._z_plot = Plot(self._channame,'z',0)
		self._1D_plot = Plot(self._channame,'1D')

	def addMCHistograms(self) :
		#for every process
		for proc in self._processes :
			print '		Getting projection histograms for process %s'%(proc.getName())
			##get the three 1-D projection histograms
			#x_histo, y_histo, z_histo = proc.getHistoProjections(self._channame)
			histo_1D = proc.getHisto1D()
#			print '		Returned histograms:' #DEBUG
#			print '			x: %s'%(x_histo) #DEBUG
#			print '			y: %s'%(y_histo) #DEBUG
#			print '			z: %s'%(z_histo) #DEBUG
			#add them to the plots
			#self._x_plot.addMChisto(x_histo,proc.getName())
			#self._y_plot.addMChisto(y_histo,proc.getName())
			#self._z_plot.addMChisto(z_histo,proc.getName())
			self._1D_plot.addMChisto(histo_1D,proc.getName())
		#after all the processes have been added, take the sqrt of the bin errors because I was incrementing by the err**2
		#self._x_plot.setMCerrors()
		#self._y_plot.setMCerrors()
		#self._z_plot.setMCerrors()
		self._1D_plot.setMCerrors()

	def buildResiduals(self) :
		#first get the 1D input data graph from the initial template file and copy it into a 1D histogram like in the MC processes
		cnameapps = ['plus','minus'] if options.sumCharges else ['']
		newname = self._channame+'__data_obs'
		newdatatemp = Template(newname+'__POSTFIT',newname+'__POSTFIT',None)
		data1Dhisto=None
		if options.mode=='postfit' :
			combine_filep = TFile.Open(options.cfilename)
			data1Dhisto = self._processes[0].getHisto1D().Clone()
			data1Dhisto.SetDirectory(0)
			data1Dhisto.Reset()
			cnamesplit = self._channame.split('_')
			for cnameapp in cnameapps :
				realchanname = cnamesplit[0]+'_'+cnamesplit[1]+cnameapp+'_'+cnamesplit[2]
				if len(cnamesplit)==4 :
					realchanname+='_'+cnamesplit[3]
				data_graph = combine_filep.Get('shapes_fit_s/%s/data'%(realchanname))
	#			print '		ngraphpoints=%d data1Dhisto: nbins=%d, integral before filling=%.2f'%(data_graph.GetN(),data1Dhisto.GetSize()-2,data1Dhisto.Integral()) #DEBUG
				for i in range(data_graph.GetN()) :
					px = array('d',[0.]); py = array('d',[0.])
					data_graph.GetPoint(i,px,py)
					err = data_graph.GetErrorY(i)
	#				if i==350 or i==475 : #DEBUG
	#					print '		bin %d in data has px=%d, py=%.4f, err=%.4f'%(i,px[0],py[0],err) #DEBUG
					data1Dhisto.Fill(px[0],py[0])
					#olderr = data1Dhisto.GetBinError(data1Dhisto.FindFixBin(px[0]))
					#data1Dhisto.SetBinError(data1Dhisto.FindFixBin(px[0]),sqrt(olderr**2+err**2))
					data1Dhisto.SetBinError(data1Dhisto.FindFixBin(px[0]),err)
			combine_filep.Close()
		elif options.mode=='prefit' :
			initial_templates_file = TFile(options.tfilename)
			data1Dhisto = initial_templates_file.Get(newname).Clone()
			data1Dhisto.SetDirectory(0)
			initial_templates_file.Close()
		newdatatemp.make_from_1D_histo(data1Dhisto)
		#self._data_histo_x=newdatatemp.getHistoX()
		#self._data_histo_y=newdatatemp.getHistoY()
		#self._data_histo_z=newdatatemp.getHistoZ()
		self._data_histo_1D=newdatatemp.convertTo1D()
#		print '		plot group data histos: 1D=%s, x=%s, y=%s, z=%s (integrals: 1D=%.2f, x=%.2f, y=%.2f, z=%.2f)'%(data1Dhisto,self._data_histo_x,self._data_histo_y,self._data_histo_z,data1Dhisto.Integral(),self._data_histo_x.Integral(),self._data_histo_y.Integral(),self._data_histo_z.Integral()) #DEBUG
		#add them to the plot
		#self._x_plot.setDataHisto(self._data_histo_x)
		#self._y_plot.setDataHisto(self._data_histo_y)
		#self._z_plot.setDataHisto(self._data_histo_z)
		self._1D_plot.setDataHisto(self._data_histo_1D)
		#calculate residuals
		#self._x_plot.calculateResiduals()
		#self._y_plot.calculateResiduals()
		#self._z_plot.calculateResiduals()
		self._1D_plot.calculateResiduals()

	def makePlots(self) :
		#self._x_canv = self._x_plot.plotOnCanvas()
		#self._y_canv = self._y_plot.plotOnCanvas()
		#self._z_canv = self._z_plot.plotOnCanvas()
		self._1D_canv = self._1D_plot.plotOnCanvas()
		outfile.cd()
		#self._x_canv.Write()
		#self._y_canv.Write()
		#self._z_canv.Write()
		self._1D_canv.Write()

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
	cnames_in_file = gDirectory.GetListOfKeys()
	#each key in this directory's name is a channel name
	print 'Setting up plot groups...'
	for procname in procs :
		for k1 in cnames_in_file :
			channel_name = k1.GetName()
			channel_id = channel_name
			print '	Adding channel with name %s'%(channel_name)
			##skip the control region
			#if channel_name.find('WJets_CR')!=-1 :
			#	continue
			#strip the "plus" or "minus" from the channel name if we want to sum over the two lepton charges
			if options.sumCharges :
				if channel_name.find('plus')!=-1 :
					cnamesplit = channel_name.split('plus')
					channel_id=cnamesplit[0]+cnamesplit[1]
				elif channel_name.find('minus')!=-1 :
					cnamesplit = channel_name.split('minus')
					channel_id=cnamesplit[0]+cnamesplit[1]
			#initialize a plot group for this channel
			if channel_id not in all_plot_groups.keys() :
				all_plot_groups[channel_id] = PlotGroup(channel_id)
			#cd to this channel's subdirectory
			combine_filep.cd('shapes_fit_s/%s'%(channel_name))
			#find the keys in this directory that are process names
			for k2 in gDirectory.GetListOfKeys() :
				if not k2.GetName()==procname :
					continue
				#add processes to the PlotGroup
				print '		Adding process with name %s'%(k2.GetName())
				all_plot_groups[channel_id].addProcess(k2.GetName(),gDirectory.Get(k2.GetName()))
			combine_filep.cd('..')
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
			#if channel_name.startswith('t1') and channel_name.split('_')[1].startswith('mu') and procname=='fqcd' :
			#	continue
			#add processes to the PlotGroup
			print '		Adding process with name %s'%(procname)
			if options.usetopptrw=='yes' :
				all_plot_groups[channel_name].addProcess(procname,template_filep.Get(channel_name+'__'+procname+'__top_pt_re_weightUp'))
			else :
				all_plot_groups[channel_name].addProcess(procname,template_filep.Get(channel_name+'__'+procname))
	template_filep.Close()

#start the output file
outname = options.outfilename if options.outfilename.endswith('.root') else options.outfilename+'.root'
outfile = TFile(outname,'recreate')

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

#put .pdfs in a folder
if not os.path.isdir(outname.replace('.root','')+'_pdfs') :
	os.system('mkdir '+outname.replace('.root','')+'_pdfs')
os.system('mv *_canv.pdf '+outname.replace('.root','')+'_pdfs')

outfile.Close()