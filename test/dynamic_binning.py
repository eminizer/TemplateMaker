#a script to help on the way to determining dynamic binning for templates

#imports
from ROOT import *
from array import array

#name of input file
ifn = '../total_ttree_files/powheg_TT_skim_all.root'

#name of output file
ofn = 'binning_investigations.root'

#weightstrings and other definitions
fqpws = '(0.5*(1.+wqa0))*((((19690.184*(lepflavor==1)+19171.010*(lepflavor==2))*sf_trig_eff_BtoF*sf_lep_ID_BtoF*sf_lep_iso_BtoF)+((16226.452*(lepflavor==1)+16214.862*(lepflavor==2))*sf_trig_eff_GH*sf_lep_ID_GH*sf_lep_iso_GH))*weight*sf_pileup*sf_ttag_eff_merged*sf_ttag_eff_semimerged*sf_ttag_eff_notmerged*sf_btag_eff_heavy*sf_btag_eff_light*sf_mu_R*sf_mu_F*sf_scale_comb*sf_pdf_alphas)'
fqmws = '(0.5*(1.-wqa0))*((((19690.184*(lepflavor==1)+19171.010*(lepflavor==2))*sf_trig_eff_BtoF*sf_lep_ID_BtoF*sf_lep_iso_BtoF)+((16226.452*(lepflavor==1)+16214.862*(lepflavor==2))*sf_trig_eff_GH*sf_lep_ID_GH*sf_lep_iso_GH))*weight*sf_pileup*sf_ttag_eff_merged*sf_ttag_eff_semimerged*sf_ttag_eff_notmerged*sf_btag_eff_heavy*sf_btag_eff_light*sf_mu_R*sf_mu_F*sf_scale_comb*sf_pdf_alphas)'
fqqws = '((((19690.184*(lepflavor==1)+19171.010*(lepflavor==2))*sf_trig_eff_BtoF*sf_lep_ID_BtoF*sf_lep_iso_BtoF)+((16226.452*(lepflavor==1)+16214.862*(lepflavor==2))*sf_trig_eff_GH*sf_lep_ID_GH*sf_lep_iso_GH))*weight*sf_pileup*sf_ttag_eff_merged*sf_ttag_eff_semimerged*sf_ttag_eff_notmerged*sf_btag_eff_heavy*sf_btag_eff_light*sf_mu_R*sf_mu_F*sf_scale_comb*sf_pdf_alphas)'
t1mup = '(fullselection==1 && M>350. && eventTopology==1 && eventType==0 && lepflavor==1 && lep_Q>0 {})*({})'
t1mum = '(fullselection==1 && M>350. && eventTopology==1 && eventType==0 && lepflavor==1 && lep_Q<0 {})*({})'
t1elp = '(fullselection==1 && M>350. && eventTopology==1 && eventType==0 && lepflavor==2 && lep_Q>0 {})*({})'
t1elm = '(fullselection==1 && M>350. && eventTopology==1 && eventType==0 && lepflavor==2 && lep_Q<0 {})*({})'
t2mup = '(fullselection==1 && M>350. && eventTopology==2 && eventType==0 && lepflavor==1 && lep_Q>0 {})*({})'
t2mum = '(fullselection==1 && M>350. && eventTopology==2 && eventType==0 && lepflavor==1 && lep_Q<0 {})*({})'
t2elp = '(fullselection==1 && M>350. && eventTopology==2 && eventType==0 && lepflavor==2 && lep_Q>0 {})*({})'
t2elm = '(fullselection==1 && M>350. && eventTopology==2 && eventType==0 && lepflavor==2 && lep_Q<0 {})*({})'
t3mup = '(fullselection==1 && M>350. && eventTopology==3 && eventType==0 && lepflavor==1 && lep_Q>0 {})*({})'
t3mum = '(fullselection==1 && M>350. && eventTopology==3 && eventType==0 && lepflavor==1 && lep_Q<0 {})*({})'
t3elp = '(fullselection==1 && M>350. && eventTopology==3 && eventType==0 && lepflavor==2 && lep_Q>0 {})*({})'
t3elm = '(fullselection==1 && M>350. && eventTopology==3 && eventType==0 && lepflavor==2 && lep_Q<0 {})*({})'
t1mupcr = '(wjets_cr_selection==1 && M>350. && eventTopology==1 && eventType==0 && lepflavor==1 && lep_Q>0 {})*({})'
t1mumcr = '(wjets_cr_selection==1 && M>350. && eventTopology==1 && eventType==0 && lepflavor==1 && lep_Q<0 {})*({})'
t1elpcr = '(wjets_cr_selection==1 && M>350. && eventTopology==1 && eventType==0 && lepflavor==2 && lep_Q>0 {})*({})'
t1elmcr = '(wjets_cr_selection==1 && M>350. && eventTopology==1 && eventType==0 && lepflavor==2 && lep_Q<0 {})*({})'
t2mupcr = '(wjets_cr_selection==1 && M>350. && eventTopology==2 && eventType==0 && lepflavor==1 && lep_Q>0 {})*({})'
t2mumcr = '(wjets_cr_selection==1 && M>350. && eventTopology==2 && eventType==0 && lepflavor==1 && lep_Q<0 {})*({})'
t2elpcr = '(wjets_cr_selection==1 && M>350. && eventTopology==2 && eventType==0 && lepflavor==2 && lep_Q>0 {})*({})'
t2elmcr = '(wjets_cr_selection==1 && M>350. && eventTopology==2 && eventType==0 && lepflavor==2 && lep_Q<0 {})*({})'
t3mupcr = '(wjets_cr_selection==1 && M>350. && eventTopology==3 && eventType==0 && lepflavor==1 && lep_Q>0 {})*({})'
t3mumcr = '(wjets_cr_selection==1 && M>350. && eventTopology==3 && eventType==0 && lepflavor==1 && lep_Q<0 {})*({})'
t3elpcr = '(wjets_cr_selection==1 && M>350. && eventTopology==3 && eventType==0 && lepflavor==2 && lep_Q>0 {})*({})'
t3elmcr = '(wjets_cr_selection==1 && M>350. && eventTopology==3 && eventType==0 && lepflavor==2 && lep_Q<0 {})*({})'

#open the input file and get the tree from it
ifp = TFile.Open(ifn,'r')
totaltree = ifp.Get('tree')

#open the output file and copy the tree
ofp = TFile.Open(ofn,'recreate')
t1fqtree = totaltree.CopyTree('fullselection==1 && eventTopology==1 && eventType==0')
t2fqtree = totaltree.CopyTree('fullselection==1 && eventTopology==2 && eventType==0')
t3fqtree = totaltree.CopyTree('fullselection==1 && eventTopology==3 && eventType==0')
t1fqcrtree = totaltree.CopyTree('wjets_cr_selection==1 && eventTopology==1 && eventType==0 && M>=350.')
t2fqcrtree = totaltree.CopyTree('wjets_cr_selection==1 && eventTopology==2 && eventType==0 && M>=350.')

#bin arrays
#cstar
cstar_in_two = array('d',[-1.00,0.00,1.00])
cstar_in_four = array('d',[-1.00,-0.50,0.00,0.50,1.00])
cstar_in_six = array('d',[-1.00,-0.50,-0.25,0.00,0.25,0.50,1.00])
cstar_in_eight = array('d',[-1.00,-0.60,-0.40,-0.20,0.00,0.20,0.40,0.60,1.00])
cstar_in_ten = array('d',[-1.00,-0.80,-0.60,-0.40,-0.20,0.00,0.20,0.40,0.60,0.80,1.00])
cstar_in_sixteen = array('d',[-1.00,-0.875,-0.75,-0.625,-0.50,-0.375,-0.25,-0.125,0.00,0.125,0.25,0.375,0.50,0.625,0.75,0.875,1.00])
cstar_in_eighteen = array('d',[-1.00,-0.80,-0.70,-0.60,-0.50,-0.40,-0.30,-0.20,-0.10,0.00,0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.80,1.00])
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

#test histograms
#xF
twobinxF = TH1D('twobinxF','twobinxF',2,array('d',[0.00,0.15,1.00]))
threebinxF = TH1D('threebinxF','threebinxF',3,array('d',[0.00,0.10,0.24,1.00]))
fourbinxF = TH1D('fourbinxF','fourbinxF',4,array('d',[0.00,0.06,0.14,0.24,1.00]))
fivebinxF = TH1D('fivebinxF','fivebinxF',5,array('d',[0.00,0.04,0.09,0.15,0.24,1.00]))
sevenbinxF = TH1D('sevenbinxF','sevenbinxF',7,array('d',[0.00,0.025,0.05,0.075,0.10,0.125,0.155,1.00]))
tenbinxF = TH1D('tenbinxF','tenbinxF',10,array('d',[0.00,0.02,0.04,0.06,0.08,0.10,0.115,0.13,0.15,0.18,1.00]))
#cstar
twobincstar = TH1D('twobincstar','twobincstar',2,-1.,1.)
fourbincstar = TH1D('fourbincstar','fourbincstar',4,-1.,1.)
sixbincstar = TH1D('sixbincstar','sixbincstar',6,cstar_in_six)
eightbincstar = TH1D('eightbincstar','eightbincstar',8,cstar_in_eight)
tenbincstar = TH1D('tenbincstar','tenbincstar',10,cstar_in_ten)
sixteenbincstar = TH1D('sixteenbincstar','sixteenbincstar',16,cstar_in_sixteen)
eighteenbincstar = TH1D('eighteenbincstar','eighteenbincstar',18,cstar_in_eighteen)
twentybincstar = TH1D('twentybincstar','twentybincstar',20,cstar_in_twenty)
#M
onebint1M = TH1D('onebint1M','onebint1M',1,t1_M_in_one)
twobint1M = TH1D('twobint1M','twobint1M',2,t1_M_in_two)
threebint1M = TH1D('threebint1M','threebint1M',3,t1_M_in_three)
sixbint1M = TH1D('sixbint1M','sixbint1M',6,t1_M_in_six)
onebint2M = TH1D('onebint2M','onebint2M',1,t2_M_in_one)
twobint2M = TH1D('twobint2M','twobint2M',2,t2_M_in_two)
threebint2M = TH1D('threebint2M','threebint2M',3,t2_M_in_three)
sixbint2M = TH1D('sixbint2M','sixbint2M',6,t2_M_in_six)

#some copied/pasted drawing commands
t2fqtree.Draw('abs(x_F)>>fivebinxF',t2mup.format('',fqpws))
t1fqcrtree.Draw('M>>twobint1M',t1mupcr.format('',fqqws))

#dictionary of decided binnings
SR_bins = {
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
}
CR_bins = {
't1':{
	'mu':t1_M_in_six,
	'el':t1_M_in_three,
	},
't2':{
	'mu':t2_M_in_six,
	'el':t2_M_in_one,
	},
}


##lists of files, cuts, and weights for each process we want to template-a-tize
#pre = '../total_ttree_files/'; post = '_skim_all.root'
#procdict = {'fqp0':{'files':[pre+'powheg_TT'+post],
#					'cuts':'eventType==0',
#					'weights':'(0.5*(1.+$wqa0$))*((((19690.184*($lepflavor$==1)+19171.010*($lepflavor$==2))*$sf_trig_eff_BtoF$*$sf_lep_ID_BtoF$*$sf_lep_iso_BtoF$)+((16226.452*($lepflavor$==1)+16214.862*($lepflavor$==2))*$sf_trig_eff_GH$*$sf_lep_ID_GH$*$sf_lep_iso_GH$))*$weight$*$sf_pileup$*$sf_ttag_eff_merged$*$sf_ttag_eff_semimerged$*$sf_ttag_eff_notmerged$*$sf_btag_eff_heavy$*$sf_btag_eff_light$*$sf_mu_R$*$sf_mu_F$*$sf_scale_comb$*$sf_pdf_alphas$)'
#					},
#		}
#
##dict of channel names with cuts and binning
#chandict = {'t1_muplus_SR':{'cuts':'eventTopology==1 && lepflavor==1 && lep_Q>0',
#							'bins':[{'ylo':0.,'yhi':0.8,
#									 'xbins':array('d',[-1.0,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,1.0]),
#									 'zbins':array('d',[500.,1350.,4000.])
#									 },
#									],
#		}
#
##list of all the objects we don't want to lose bc pyROOT : /
#objs = []
#
##quick template class
#class QuickTemplate(object) :
#	#initializer
#	def __init__(self,name,bindictlist) :
#		self._name = name
#		self._histos = {}
#		n1Dbins=0
#		for ib,bindict in self._bindictlist.enumerate() :
#			xbins=bindict['xbins']; zbins=bindict['zbins']
#			tempname=name+'_xz'+str(ib)
#			self._histos[ib]={'ylo':bindict['ylo'],'yhi':bindict['yhi'],'xzhist':TH2D(tempname,tempname,len(xbins)-1,xbins,len(zbins)-1,zbins)}
#			n1Dbins+=(len(xbins)-1)*(len(zbins)-1)
#		#empty 1-dimensional histogram
#		self._histo1D=TH1D(name,name,n1Dbins,0,n1Dbins-1)
#	#build the template from the file chain, the cutstring, and the weightstring
#	def build(self,chain,cutstring,weightstring) :
#		#first make the cuts
#		tree = chain.CopyTree(cutstring)
#		#then loop over events in the tree and add them to the template
#
##open the output file
#ofp = TFile.Open(outfilename,'recreate')
#
##loop over channels
#for cn in chandict.keys() :
#	#and over processes
#	for pn in procdict.keys() :
#		#initialize the template
#		template = QuickTemplate(cn+'__'+pn,chandict[cn]['xbins'],chandict[cn]['ybins'],chandict[cn]['zbins'])
#		#open the input files and chain them up
#		chain = TChain('tree')
#		for fn in procdict[pn]['files'] :
#			chain.Add(fn)
#		#build the template from the chain, the cuts, and the weights
#		allcuts = '(('+chandict[cn]['cuts']+') && ('+procdict[pn]['cuts']+'))'
#		weights = procdict[pn]['weights']
#		template.build(chain,allcuts,weights)

