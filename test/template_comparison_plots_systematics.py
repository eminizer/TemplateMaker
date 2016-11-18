from ROOT import *

#dimension
dim = 'x'

#open input file
infilep = TFile('../templates_aux.root')

#channel names to sum over
channels = ['muplus','muminus','elplus','elminus']

#lists of systematics names and colors
systematics = [('pileup_weight',kMagenta),('JES',kRed+2),('JER',kRed-7),('lep_ID_weight',kAzure),('luminosity',kOrange-3),('ren_scale_weight',kGreen+2),('fact_scale_weight',kGreen-7),('comb_scale_weight',kSpring-7),('pdfas_weight',kSpring+3)]
systematics_names = ['Pileup Weight','Jet Energy Scale','Jet Energy Resolution','Lepton ID weight','Luminosity','Renormalization scale','Factorization scale','Combined #mu_{R}/#mu_{F} scale','PDF/#alpha_{s} weight']
#systematics = [('fit',kRed)]
#systematics_names = ['Conv. func. fit']
#list of histograms
hists = []
#one list for each of the qq, gg, bck, and NTMJ
for i in range(4) :
	hists.append([])
#get the nominal templates
hists[0].append(infilep.Get(channels[0]+'__fqq_'+dim))
hists[1].append(infilep.Get(channels[0]+'__fgg_'+dim))
hists[2].append(infilep.Get(channels[0]+'__fbck_'+dim))
hists[3].append(infilep.Get(channels[0]+'__fntmj_'+dim))
#add from all the channels
for i in range(1,len(channels)) :
	hists[0][0].Add(infilep.Get(channels[i]+'__fqq_'+dim).Clone())
	hists[1][0].Add(infilep.Get(channels[i]+'__fgg_'+dim).Clone())
	hists[2][0].Add(infilep.Get(channels[i]+'__fbck_'+dim).Clone())
	hists[3][0].Add(infilep.Get(channels[i]+'__fntmj_'+dim).Clone())
#Get the systematics up/down templates
for i in range(len(systematics)) :
	hists[0].append(infilep.Get(channels[0]+'__fqq__'+systematics[i][0]+'__up_'+dim))
	hists[0].append(infilep.Get(channels[0]+'__fqq__'+systematics[i][0]+'__down_'+dim))
	hists[1].append(infilep.Get(channels[0]+'__fgg__'+systematics[i][0]+'__up_'+dim))
	hists[1].append(infilep.Get(channels[0]+'__fgg__'+systematics[i][0]+'__down_'+dim))
	hists[2].append(infilep.Get(channels[0]+'__fbck__'+systematics[i][0]+'__up_'+dim))
	hists[2].append(infilep.Get(channels[0]+'__fbck__'+systematics[i][0]+'__down_'+dim))
	hists[3].append(infilep.Get(channels[0]+'__fntmj__'+systematics[i][0]+'__up_'+dim))
	hists[3].append(infilep.Get(channels[0]+'__fntmj__'+systematics[i][0]+'__down_'+dim))
	#Add from the rest of the channels
	for j in range(1,len(channels)) :
		hists[0][2*i+1].Add(infilep.Get(channels[j]+'__fqq__'+systematics[i][0]+'__up_'+dim).Clone())
		hists[0][2*i+2].Add(infilep.Get(channels[j]+'__fqq__'+systematics[i][0]+'__down_'+dim).Clone())
		hists[1][2*i+1].Add(infilep.Get(channels[j]+'__fgg__'+systematics[i][0]+'__up_'+dim).Clone())
		hists[1][2*i+2].Add(infilep.Get(channels[j]+'__fgg__'+systematics[i][0]+'__down_'+dim).Clone())
		hists[2][2*i+1].Add(infilep.Get(channels[j]+'__fbck__'+systematics[i][0]+'__up_'+dim).Clone())
		hists[2][2*i+2].Add(infilep.Get(channels[j]+'__fbck__'+systematics[i][0]+'__down_'+dim).Clone())
		hists[3][2*i+1].Add(infilep.Get(channels[j]+'__fntmj__'+systematics[i][0]+'__up_'+dim).Clone())
		hists[3][2*i+2].Add(infilep.Get(channels[j]+'__fntmj__'+systematics[i][0]+'__down_'+dim).Clone())

#open the output file
outfilep = TFile('template_comparison_plots_systematics_'+dim+'.root','recreate')

#Set histogram attributes
for i in range(len(hists)) :
	hists[i][0].SetLineWidth(4); hists[i][0].SetLineColor(kBlack)
	for j in range(1,len(hists[i])) :
		hists[i][j].SetLineWidth(4); hists[i][j].SetLineColor(systematics[(j-1)/2][1])
		if (j-1)%2 == 0 :
			hists[i][j].SetLineStyle(7)
		else :
			hists[i][j].SetLineStyle(3)
projtype = 'c*'
if dim == 'y' :
	projtype = '|x_{F}|'
if dim == 'z' :
	projtype = 'M'
hists[0][0].SetTitle('Systematic variations in q#bar{q} template, '+projtype+' projection; '+projtype+'; Events')
hists[1][0].SetTitle('Systematic variations in gg template, '+projtype+' projection; '+projtype+'; Events')
hists[2][0].SetTitle('Systematic variations in simulated background template, '+projtype+' projection; '+projtype+'; Events')
hists[3][0].SetTitle('Systematic variations in data-driven NTMJ background template, '+projtype+' projection; '+projtype+'; Events')

#make a legend
leg = TLegend(0.62,0.67,0.9,0.9)
leg.AddEntry(hists[0][0],'Nominal Template','L')
for i in range(len(systematics)) :
	leg.AddEntry(hists[0][2*i+1],systematics_names[i]+' Up','L') 
	leg.AddEntry(hists[0][2*i+2],systematics_names[i]+' Down','L') 

#canvases
canvs = [TCanvas('qq_canv','qq_canv',1100,900),TCanvas('gg_canv','gg_canv',1100,900),TCanvas('bck_canv','bck_canv',1100,900),TCanvas('ntmj_canv','ntmj_canv',1100,900)]

#plot plots
for i in range(len(hists)) :
	canvs[i].cd()
	hists[i][0].Draw('HIST')
	for j in range(len(systematics)) :
		hists[i][2*j+1].Draw('SAME HIST')
		hists[i][2*j+2].Draw('SAME HIST')
	leg.Draw()

#save canvases
outfilep.cd()
for canv in canvs :
	canv.Write()

