from ROOT import *
import CMS_lumi, tdrstyle

#TDR plot style stuff
gROOT.SetBatch()
tdrstyle.setTDRStyle()
iPeriod = 4 #13TeV iPeriod = 1*(0/1 7 TeV) + 2*(0/1 8 TeV)  + 4*(0/1 13 TeV)
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Preliminary"

#dimension
dim = 'x'

#open input file
infilep = TFile('../total_template_files/templates_powheg_aggregated_all_aux.root')

#lepton types names to sum over
leptypes = ['elplus','elminus','muplus','muminus']
#leptypes = ['elplus']
#leptypes = ['elminus']
#leptypes = ['muplus']
#leptypes = ['muminus']

#lists of systematics names and colors
systematics = [('pileup_weight',			 kRed,'Pileup Weight'),
			   ('JES',						 kRed+3,'Jet Energy Scale'),
			   ('JER',						 kRed-9,'Jet Energy Resolution'),
			   ('trig_eff_weight',			 kBlue,'Trigger eff.'),
			   ('lep_ID_weight',			 kBlue+3,'Lepton ID eff.'),
			   ('lep_iso_weight',			 kBlue-9,'Lepton isolation eff.'),
			   ('btag_eff_weight_flavb',	 kGreen-5,'b-tagging eff. (b-flav.)'),
			   ('btag_eff_weight_flavc',	 kGreen+3,'b-tagging eff. (c-flav.)'),
			   ('btag_eff_weight_light',	 kGreen-8,'b-tagging eff. (light flavs.)'),
			   ('ttag_eff_weight_merged',	 kOrange,'merged top-tagging eff.'),
			   ('ttag_eff_weight_semimerged',kOrange+2,'semimerged top-tagging eff.'),
			   ('ttag_eff_weight_notmerged', kOrange-6,'notmerged top-tagging eff.'),
			   ('top_pt_re_weight',			 kCyan,'top p_{T} reweight'),
			   ('lumi',						 kViolet-3,'Luminosity'),
			   ('ren_scale_weight',			 kMagenta,'Renormalization scale'),
			   ('fact_scale_weight',		 kMagenta+3,'Factorization scale'),
			   ('comb_scale_weight',		 kMagenta-10,'Combined #mu_{R}/#mu_{F} scale'),
			   ('pdfas_weight',				 kPink+6,'PDF/#alpha_{s} weight'),
			   ('B_frag_weight',			 kYellow+3,'B fragmentation weight'),
			   ('B_br_weight',				 kYellow-8,'B BR weight'),
			   ('isr',				 		 kPink-8,'ISR'),
			   ('fsr',				 		 kPink-6,'FSR'),
			   ('hdamp',				 	 kOrange,'hDamp'),
			   ('tune',				 		 kOrange+7,'underlying event'),
			   ('cr', 						 kViolet+2,'color reconnection'),
			   ]
ttbar_only = ['top_pt_re_weight','ren_scale_weight','fact_scale_weight','comb_scale_weight','pdfas_weight','B_frag_weight','B_br_weight','isr','fsr','hdamp','tune','cr']
#dictionary of histogram lists
#first key: topology
#second key: region
#third key: type (qq,gg,bck,WJets,QCD)
hists = {'t1':{'SR':{},'WJets_CR':{}},'t2':{'SR':{},'WJets_CR':{}},'t3':{'SR':{}}}
types = ['fqq','fgg','fbck','fwjets','fqcd']
for top in hists :
	for reg in hists[top] :
		for t in types :
			hists[top][reg][t]=[]
			#get the nominal templates, starting with the first lepton type
			print 'looking for histogram with name %s'%(top+'_'+leptypes[0]+'_'+reg+'__'+t+'_'+dim) #DEBUG
			hists[top][reg][t].append(infilep.Get(top+'_'+leptypes[0]+'_'+reg+'__'+t+'_'+dim).Clone())
			#add from all the lepton types
			for i in range(1,len(leptypes)) :
				if top=='t1' and leptypes[i].startswith('mu') and t=='fqcd' :
						continue
				hists[top][reg][t][-1].Add(infilep.Get(top+'_'+leptypes[i]+'_'+reg+'__'+t+'_'+dim).Clone())
			#Get the systematics up/down templates
			for i in range(len(systematics)) :
				basesysname = systematics[i][0]
				if t in ['fwjets','fqcd'] and basesysname in ttbar_only :
					continue
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
				elif basesysname.startswith('btag_eff_weight') :#or basesysname.startswith('ttag_eff_weight') :
					newsysname = basesysname + ('_r' if top=='t3' else '_b')
					for j in range(len(leptypes)) :
						sysnames.append(newsysname)
				else :
					for j in range(len(leptypes)) :
						sysnames.append(basesysname)
				#get the first lepton type
				print 'looking for histogram with name %s'%(top+'_'+leptypes[0]+'_'+reg+'__'+t+'__'+sysnames[0]+'Up_'+dim) #DEBUG
				hists[top][reg][t].append(infilep.Get(top+'_'+leptypes[0]+'_'+reg+'__'+t+'__'+sysnames[0]+'Up_'+dim).Clone())
				hists[top][reg][t].append(infilep.Get(top+'_'+leptypes[0]+'_'+reg+'__'+t+'__'+sysnames[0]+'Down_'+dim).Clone())
				#Add from the rest of the lepton types
				for j in range(1,len(leptypes)) :
					if top=='t1' and leptypes[j].startswith('mu') and t=='fqcd' :
						continue
					hists[top][reg][t][-2].Add(infilep.Get(top+'_'+leptypes[j]+'_'+reg+'__'+t+'__'+sysnames[j]+'Up_'+dim).Clone())
					hists[top][reg][t][-1].Add(infilep.Get(top+'_'+leptypes[j]+'_'+reg+'__'+t+'__'+sysnames[j]+'Down_'+dim).Clone())

#open the output file
outfilep = TFile('template_comparison_plots_systematics_'+dim+'.root','recreate')

#Set histogram attributes
projtype = 'c*'
if dim == 'y' :
	projtype = '|x_{F}|'
if dim == 'z' :
	projtype = 'M'
for top in hists :
	for reg in hists[top] :
		for t in hists[top][reg] :
			hists[top][reg][t][0].SetLineWidth(4); hists[top][reg][t][0].SetLineColor(kBlack)
			hists[top][reg][t][0].SetTitle(''); hists[top][reg][t][0].SetStats(0)
			hists[top][reg][t][0].GetYaxis().SetTitle('Events')
			hists[top][reg][t][0].GetYaxis().SetTitleOffset(1.5)
			hists[top][reg][t][0].GetXaxis().SetTitle(projtype)
			thismax = hists[top][reg][t][0].GetMaximum()
			for j in range(1,len(hists[top][reg][t])) :
				if hists[top][reg][t][j].GetMaximum()>thismax :
					thismax=hists[top][reg][t][j].GetMaximum()
			hists[top][reg][t][0].GetYaxis().SetRangeUser(0.,1.1*thismax)
			for j in range(1,len(hists[top][reg][t])) :
				hists[top][reg][t][j].SetLineWidth(4); hists[top][reg][t][j].SetLineColor(systematics[(j-1)/2][1])
				if (j-1)%2 == 0 :
					hists[top][reg][t][j].SetLineStyle(7)
				else :
					hists[top][reg][t][j].SetLineStyle(3)

#make a legend
leg = TLegend(0.12,0.17,0.9,0.2)
leg.SetNColumns(4)
leg.AddEntry(hists['t1']['SR']['fqq'][0],'Nominal Template','L')
leg.AddEntry(None,''); leg.AddEntry(None,''); leg.AddEntry(None,'')
for i in range(len(systematics)) :
	leg.AddEntry(hists['t1']['SR']['fqq'][2*i+1],systematics[i][2]+' Up','L') 
	leg.AddEntry(hists['t1']['SR']['fqq'][2*i+2],systematics[i][2]+' Down','L') 

#canvases
canvs = {}
for top in hists :
	canvs[top] = {}
	for reg in hists[top] :
		canvs[top][reg] = {}
		for t in hists[top][reg] :
			canvs[top][reg][t]=TCanvas(top+'_'+reg+'_'+t+'_canv',top+'_'+reg+'_'+t+'_canv',1100,900)

all_cms_lumi_objs = []
#plot plots
for top in hists :
	for reg in hists[top] :
		for t in hists[top][reg] :
			canvs[top][reg][t].cd()
			hists[top][reg][t][0].Draw('HIST')
			for i in range(1,len(hists[top][reg][t])) :
				hists[top][reg][t][i].Draw('SAME HIST')
				hists[top][reg][t][i].Draw('SAME HIST')
			leg.Draw()
			all_cms_lumi_objs.append(CMS_lumi.CMS_lumi(canvs[top][reg][t], iPeriod, 0))

#save canvases
outfilep.cd()
for top in hists :
	for reg in hists[top] :
		for t in hists[top][reg] :
			canvs[top][reg][t].Write()

