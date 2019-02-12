from ROOT import *
#import tdrstyle
from math import sqrt

#TDR plot style stuff
gROOT.SetBatch()
#tdrstyle.setTDRStyle()

#open input file
infilep = TFile('../total_template_files/templates_powheg_aggregated_all.root')

#lepton types to sum over
#leptypes = ['elplus','elminus','muplus','muminus']
leptypes = ['muplus']

#lists of systematics names 
systematics = [('pileup_weight','pileup'),
			   ('top_pt_re_weight','top pT rw'),
			   ('B_frag_weight','B fragmentation'),
			   ('B_br_weight','B branching ratio'),
			   ('isr','ISR'),
			   ('fsr','FSR'),
			   ('hdamp','h damp'),
			   ('tune','underlying event'),
			   ('cr','color reconnection'),
			   ]
#dictionary of histogram lists
#first key: topology
#second key: region
#third key: type (qq,gg,bck,WJets,QCD)
#hists = {'t1':{'SR':{},'WJets_CR':{}},'t2':{'SR':{},'WJets_CR':{}},'t3':{'SR':{}}}
hists = {'t3':{'SR':{}}}
types = ['fqp0','fqm0','fg0','fbck']
for top in hists :
	for reg in hists[top] :
		for t in types :
			hists[top][reg][t]={}
			#get the nominal templates, starting with the first lepton type
			print 'looking for histogram with name %s'%(top+'_'+leptypes[0]+'_'+reg+'__'+t) #DEBUG
			hists[top][reg][t]['nominal']=infilep.Get(top+'_'+leptypes[0]+'_'+reg+'__'+t)
			#add from all the lepton types
			for i in range(1,len(leptypes)) :
				hists[top][reg][t]['nominal'].Add(infilep.Get(top+'_'+leptypes[i]+'_'+reg+'__'+t))
			#Get the systematics up/down templates
			for sys in systematics :
				basesysname = sys[0]
				#make the list of different names that this systematic has in each lepton flavor category
				sysnames = []
				#some systematics change names based on channel/topology
				if basesysname=='trig_eff_weight' :
					for j in range(len(leptypes)) :
						newsysname = 'mu_trig_eff_weight' if leptypes[j].startswith('mu') else 'el_trig_eff_weight'
						newsysname+='_r' if top=='t3' else '_b'
						sysnames.append(newsysname)
				elif basesysname=='lep_ID_weight' :
					for j in range(len(leptypes)) :
						newsysname = 'mu_ID_weight' if leptypes[j].startswith('mu') else 'el_ID_weight'
						sysnames.append(newsysname)
				elif basesysname=='lep_iso_weight' :
					for j in range(len(leptypes)) :
						newsysname = 'mu_iso_weight' if leptypes[j].startswith('mu') else 'el_iso_weight'
						sysnames.append(newsysname)
				elif basesysname.startswith('btag_eff_weight') :
					newsysname = basesysname + ('_r' if top=='t3' else '_b')
					for j in range(len(leptypes)) :
						sysnames.append(newsysname)
				else :
					for j in range(len(leptypes)) :
						sysnames.append(basesysname)
				#get the systematic up/down template for the first lepton type
				print 'looking for histograms with name %s'%(top+'_'+leptypes[0]+'_'+reg+'__'+t+'__'+sysnames[0]+'(Up/Down)') #DEBUG
				hists[top][reg][t][basesysname+'_up']=infilep.Get(top+'_'+leptypes[0]+'_'+reg+'__'+t+'__'+sysnames[0]+'Up')
				hists[top][reg][t][basesysname+'_down']=infilep.Get(top+'_'+leptypes[0]+'_'+reg+'__'+t+'__'+sysnames[0]+'Down')
				#Add to it from the rest of the lepton types
				for j in range(1,len(leptypes)) :
					hists[top][reg][t][basesysname+'_up'].Add(infilep.Get(top+'_'+leptypes[j]+'_'+reg+'__'+t+'__'+sysnames[j]+'Up'))
					hists[top][reg][t][basesysname+'_down'].Add(infilep.Get(top+'_'+leptypes[j]+'_'+reg+'__'+t+'__'+sysnames[j]+'Down'))

#open the output file
outfilep = TFile('nuisance_investigation_plots.root','recreate')

#Set histogram attributes for 1D comparison plots
for top in hists :
	for reg in hists[top] :
		for t in hists[top][reg] :
			thisnomhist = hists[top][reg][t]['nominal']
			thisnomhist.SetLineWidth(2); thisnomhist.SetLineColor(kBlack)
			thisnomhist.SetTitle(''); thisnomhist.SetStats(0)
			thisnomhist.GetYaxis().SetTitle('Events')
			thisnomhist.GetXaxis().SetTitle('Bin')
			for sys in systematics :
				hists[top][reg][t][sys[0]+'_up'].SetLineWidth(2); hists[top][reg][t][sys[0]+'_up'].SetLineColor(kRed+2); hists[top][reg][t][sys[0]+'_up'].SetStats(0)
				hists[top][reg][t][sys[0]+'_down'].SetLineWidth(2); hists[top][reg][t][sys[0]+'_down'].SetLineColor(kBlue+2)

#canvases
canvs = {}
for top in hists :
	canvs[top] = {}
	for reg in hists[top] :
		canvs[top][reg] = {}
		for sys in systematics :
			canvs[top][reg][sys[0]] = TCanvas(top+'_'+reg+'_'+sys[0]+'_canv',top+'_'+reg+'_'+sys[0]+'_canv',1500,900)

#plot plots
for top in canvs :
	for reg in canvs[top] :
		for s in canvs[top][reg] :
			#divide the canvas to plot every process
			thiscanv = canvs[top][reg][s]
			thiscanv.Divide(2,2)
			#for each process...
			for i in range(len(types)) :
				#get the nominal and up/down histograms, and the ratio histogram
				thisnomhist = hists[top][reg][types[i]]['nominal']
				thisuphist = hists[top][reg][types[i]][s+'_up']
				thisdnhist = hists[top][reg][types[i]][s+'_down']
				#cd to this portion of the canvas
				thiscanv.cd(i+1)
				#draw the template histograms
				thisnomhist.SetTitle(top+' '+reg+' '+s+' '+types[i])
				thisnomhist.Draw("HIST"); 
				thisuphist.Draw("HIST SAME"); thisdnhist.Draw("HIST SAME")
				thiscanv.Update()
			thiscanv.Write()

#close output file
outfilep.Close()

