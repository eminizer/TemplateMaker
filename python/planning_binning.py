from ROOT import *
import tdrstyle
from math import sqrt
from template import Template

#TDR plot style stuff
gROOT.SetBatch()
tdrstyle.setTDRStyle()

#output file name
outfilen = 'binning_investigation.root'

#open input file
infilep = TFile('../total_template_files/templates_powheg_all.root')

#lepton types
leptypes = ['elplus','elminus','muplus','muminus']
#processes
procs = ['fqp0','fqm0','fg0','fbck','fwjets','fqcd']
#analysis channels
templates = {'t1':{'SR':{},'WJets_CR':{}},'t2':{'SR':{},'WJets_CR':{}},'t3':{'SR':{}}}
for t in templates :
	for r in templates[t] :
		templates[t][r] = {}
		for lt in leptypes :
			templates[t][r][lt] = {}
			for p in procs :
				templates[t][r][lt][p] = {}
				n = t+'_'+lt+'_'+r+'__'+p
				templates[t][r][lt][p]['name']=n
				thishist = infilep.Get(n)
				templates[t][r][lt][p]['hist']=thishist
				newtemp = Template(n,n,None)
				newtemp.make_from_1D_histo(thishist)
				templates[t][r][lt][p]['temp']=newtemp

#open the output file
outfilep = TFile(outfilen,'recreate')

#hold on to all the plots
allplots = []
#declare the individual channel plots
ind_plots = {}
for t in templates :
	ind_plots[t]={}
	for r in templates[t] :
		ind_plots[t][r] = {}
		for lt in templates[t][r] :
			ind_plots[t][r][lt] = {}
			for p in templates[t][r][lt] :
				ind_plots[t][r][lt][p] = {}
				n = templates[t][r][lt][p]['name']+'_bad_bins_'
				this1Dhist = templates[t][r][lt][p]['temp'].convertTo1D()
				this3Dhist = templates[t][r][lt][p]['temp'].getHisto3D()
				ind_plots[t][r][lt][p]['bad_bins_x'] = TH1(n+'x',n+'x',this3Dhist.GetNbinsX(),0.5,this3Dhist.GetNbinsX()+0.5)
				ind_plots[t][r][lt][p]['bad_bins_y'] = TH1(n+'y',n+'y',this3Dhist.GetNbinsY(),0.5,this3Dhist.GetNbinsY()+0.5)
				ind_plots[t][r][lt][p]['bad_bins_z'] = TH1(n+'z',n+'z',this3Dhist.GetNbinsZ(),0.5,this3Dhist.GetNbinsZ()+0.5)
				allplots.append(ind_plots[t][r][lt][p]['bad_bins_x'])
				allplots.append(ind_plots[t][r][lt][p]['bad_bins_y'])
				allplots.append(ind_plots[t][r][lt][p]['bad_bins_z'])
				ind_plots[t][r][lt][p]['bad_bins_N'] = 0.
				#search template for zeroed bins




