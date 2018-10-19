from ROOT import *
#import tdrstyle
from math import sqrt

#TDR plot style stuff
gROOT.SetBatch()
#tdrstyle.setTDRStyle()

#open input file
infilep = TFile('../total_template_files/templates_powheg_all.root')

#lepton types to sum over
leptypes = ['elplus','elminus','muplus','muminus']

#lists of systematics names 
systematics = [('pileup_weight','Pileup Weight'),
			   ('JES','Jet Energy Scale'),
			   ('JER','Jet Energy Resolution'),
			   ('trig_eff_weight','Trigger eff.'),
			   ('lep_ID_weight','Lepton ID eff.'),
			   ('lep_iso_weight','Lepton isolation eff.'),
			   ('btag_eff_weight_flavb','b-tagging eff. (b flavor)'),
			   ('btag_eff_weight_flavc','b-tagging eff. (c flavor)'),
			   ('btag_eff_weight_light','b-tagging eff. (light flavors)'),
			   ('ttag_eff_weight','top-tagging eff.'),
			   ('top_pt_re_weight','top p_{T} reweight'),
			   ('lumi','Luminosity'),
			   ('ren_scale_weight','Renormalization scale'),
			   ('fact_scale_weight','Factorization scale'),
			   ('comb_scale_weight','Combined #mu_{R}/#mu_{F} scale'),
			   ('pdfas_weight','PDF/#alpha_{s} weight'),
			   ]
#dictionary of histogram lists
#first key: topology
#second key: region
#third key: type (qq,gg,bck,WJets,QCD)
hists = {'t1':{'SR':{},'WJets_CR':{}},'t2':{'SR':{},'WJets_CR':{}},'t3':{'SR':{}}}
types = ['fqp0','fqm0','fg0','fbck','fwjets','fqcd']
for top in hists :
	for reg in hists[top] :
		for t in types :
			hists[top][reg][t]={}
			#get the nominal templates, starting with the first lepton type
			print 'looking for histogram with name %s'%(top+'_'+leptypes[0]+'_'+reg+'__'+t) #DEBUG
			hists[top][reg][t]['nominal']=infilep.Get(top+'_'+leptypes[0]+'_'+reg+'__'+t)#.Clone()
			#add from all the lepton types
			for i in range(1,len(leptypes)) :
				if top=='t1' and leptypes[i].startswith('mu') and t=='fqcd' :
						continue
				hists[top][reg][t]['nominal'].Add(infilep.Get(top+'_'+leptypes[i]+'_'+reg+'__'+t))#.Clone())
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
				elif basesysname.startswith('btag_eff_weight') or basesysname.startswith('ttag_eff_weight') :
					newsysname = basesysname + ('_r' if top=='t3' else '_b')
					for j in range(len(leptypes)) :
						sysnames.append(newsysname)
				else :
					for j in range(len(leptypes)) :
						sysnames.append(basesysname)
				#get the systematic up/down template for the first lepton type
				print 'looking for histogram with name %s'%(top+'_'+leptypes[0]+'_'+reg+'__'+t+'__'+sysnames[0]+'Up') #DEBUG
				hists[top][reg][t][basesysname+'_up']=infilep.Get(top+'_'+leptypes[0]+'_'+reg+'__'+t+'__'+sysnames[0]+'Up')#.Clone()
				hists[top][reg][t][basesysname+'_down']=infilep.Get(top+'_'+leptypes[0]+'_'+reg+'__'+t+'__'+sysnames[0]+'Down')#.Clone()
				#Add to it from the rest of the lepton types
				for j in range(1,len(leptypes)) :
					hists[top][reg][t][basesysname+'_up'].Add(infilep.Get(top+'_'+leptypes[j]+'_'+reg+'__'+t+'__'+sysnames[j]+'Up'))#.Clone())
					hists[top][reg][t][basesysname+'_down'].Add(infilep.Get(top+'_'+leptypes[j]+'_'+reg+'__'+t+'__'+sysnames[j]+'Down'))#.Clone())

#Build the ratio plot dictionary
ratio_plots = {}
for top in hists :
	ratio_plots[top] = {}
	for reg in hists[top] :
		ratio_plots[top][reg] = {}
		for sys in systematics :
			ratio_plots[top][reg][sys[0]] = {}
			for t in types :
				#make the ratio histogram
				thisrhist = hists[top][reg][t]['nominal'].Clone()
				thisrhist.Reset()
				for i in range(1,thisrhist.GetNbinsX()+1) :
					nomcont = hists[top][reg][t]['nominal'].GetBinContent(i)
					upcont = hists[top][reg][t][sys[0]+'_up'].GetBinContent(i)
					downcont = hists[top][reg][t][sys[0]+'_down'].GetBinContent(i)
					total_shift = upcont-downcont
					total_shift_err = sqrt(upcont+downcont)
					cont = 0.; err = 0.
					if nomcont!=0. :
						cont=total_shift/nomcont
						if total_shift!=0. :
							err = cont*sqrt((total_shift_err/total_shift)**2+(1./nomcont))
					thisrhist.SetBinContent(i,cont); thisrhist.SetBinError(i,err)
				#set some attributes
				thisrhist.SetLineWidth(2); thisrhist.SetLineColor(kGreen+2)
				thisrhist.SetMarkerColor(kGreen+2)
				thisrhist.SetMarkerStyle(20)
				#thisrhist.SetMarkerSize(2)
				thisrhist.SetTitle(sys[0]+' '+t+' shifts in '+top+'_'+reg); thisrhist.SetStats(0)
				ratio_plots[top][reg][sys[0]][t] = thisrhist

#open the output file
outfilep = TFile('nuisance_investigation_plots.root','recreate')

#Set histogram attributes for 1D comparison plots
for top in hists :
	for reg in hists[top] :
		for t in hists[top][reg] :
			thisnomhist = hists[top][reg][t]['nominal']
			thisnomhist.SetLineWidth(4); thisnomhist.SetLineColor(kBlack)
			thisnomhist.SetTitle(''); thisnomhist.SetStats(0)
			thisnomhist.GetYaxis().SetTitle('Events')
			thisnomhist.GetXaxis().SetTitle('Bin')
			for sys in systematics :
				hists[top][reg][t][sys[0]+'_up'].SetLineWidth(4); hists[top][reg][t][sys[0]+'_up'].SetLineColor(kRed+2); hists[top][reg][t][sys[0]+'_up'].SetStats(0)
				hists[top][reg][t][sys[0]+'_down'].SetLineWidth(4); hists[top][reg][t][sys[0]+'_down'].SetLineColor(kBlue+2)

##make a legend
#leg = TLegend(0.62,0.67,0.9,0.9)
#leg.AddEntry(hists['t1']['SR']['fqq'][0],'Nominal Template','L')
#for i in range(len(systematics)) :
#	leg.AddEntry(hists['t1']['SR']['fqq'][2*i+1],systematics[i][2]+' Up','L') 
#	leg.AddEntry(hists['t1']['SR']['fqq'][2*i+2],systematics[i][2]+' Down','L') 

#canvases
canvs = {}
for top in ratio_plots :
	canvs[top] = {}
	for reg in ratio_plots[top] :
		canvs[top][reg] = {}
		for sys in systematics :
			canvs[top][reg][sys[0]] = TCanvas(top+'_'+reg+'_'+sys[0]+'_canv',top+'_'+reg+'_'+sys[0]+'_canv',1500,900)

#plot plots
avg_shifts = {}
for top in canvs :
	avg_shifts[top] = {}
	for reg in canvs[top] :
		avg_shifts[top][reg] = {}
		for s in canvs[top][reg] :
			avg_shifts[top][reg][s] = 0.
			#divide the canvas to plot every process
			thiscanv = canvs[top][reg][s]
			thiscanv.Divide(2,3)
			#calculate the average shift 
			#for each process...
			for i in range(len(types)) :
				#get the nominal and up/down histograms, and the ratio histogram
				thisnomhist = hists[top][reg][types[i]]['nominal']
				thisuphist = hists[top][reg][types[i]][s+'_up']
				thisdnhist = hists[top][reg][types[i]][s+'_down']
				thisratiohist = ratio_plots[top][reg][s][types[i]]
				avg_shifts[top][reg][s]+=thisratiohist.GetMean()
				#cd to this portion of the canvas
				thiscanv.cd(i+1)
				#draw the ratio histogram
				thisratiohist.Draw("HIST P"); thiscanv.Update()
				#get the scale for the histograms
				scale = (gPad.GetUymax())/(thisuphist.GetMaximum())
				#print 'scale for %s %s shifts in %s = %.4f/%.4f=%.4f'%(s,types[i],top+'_'+reg,gPad.GetUymax(),thisuphist.GetMaximum(),scale) #DEBUG
				#scale and draw the other histograms
				#thisnomhist.Draw("HIST SAME"); 
				#thisnomhist.Scale(scale); 
				thisuphist.Scale(scale); thisdnhist.Scale(scale)
				#do it all over again now that the objects are scaled
				#thisratiohist.Draw("HIST P")
				#thisnomhist.Draw("HIST SAME"); 
				#thisuphist.Draw("HIST SAME"); thisdnhist.Draw("HIST SAME")
				#set the max of this up histogram to show the shifts
				thisuphist.SetMinimum(min(0.,1.1*min(scale*thisnomhist.GetMinimum(),thisuphist.GetMinimum(),thisdnhist.GetMinimum(),thisratiohist.GetMinimum())))
				thisuphist.SetMaximum(1.1*max(scale*thisnomhist.GetMaximum(),thisuphist.GetMaximum(),thisdnhist.GetMaximum(),thisratiohist.GetMaximum()))
				thisuphist.SetTitle(s+' '+types[i]+' shifts in '+top+'_'+reg)
				thisuphist.Draw("HIST"); thisnomhist.DrawNormalized("HIST SAME",scale*thisnomhist.Integral()); thisdnhist.Draw("HIST SAME")
				thisratiohist.Draw("HIST P SAME")
				thiscanv.Update()
			thiscanv.Write()
			avg_shifts[top][reg][s]/=len(types)
			print top+'_'+reg+' '+s+' average shift = '+str(avg_shifts[top][reg][s])

#save canvases
#outfilep.cd()
#for top in canvs :
#	for reg in canvs[top] :
#		for s in canvs[top][reg] :
#			canvs[top][reg][s].Write()
outfilep.Close()

