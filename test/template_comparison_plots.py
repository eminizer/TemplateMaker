from ROOT import *

#open input file
infilep = TFile('../templates_aux.root')

#channel names to sum over
channels = ['muplus','muminus','elplus','elminus']

#histogram name stems
noms = 		['__fqq', 				 '__fgg', 				'__fbck', '__fntmj']
afb_ups = 	['__fqq__par_Afb__up', 	 None, 					None, 	  '__fntmj__par_Afb__up']
afb_downs = ['__fqq__par_Afb__down', None, 					None, 	  '__fntmj__par_Afb__down']
d_ups = 	['__fqq__par_d__up', 	 '__fgg__par_d__up', 	None, 	  '__fntmj__par_d__up']
d_downs = 	['__fqq__par_d__down', 	 '__fgg__par_d__down', 	None, 	  '__fntmj__par_d__down']
mu_ups = 	['__fqq__par_mu__up', 	 '__fgg__par_mu__up', 	None, 	  '__fntmj__par_mu__up']
mu_downs = 	['__fqq__par_mu__down',  '__fgg__par_mu__down', None, 	  '__fntmj__par_mu__down']
allnames = [noms,afb_ups,afb_downs,d_ups,d_downs,mu_ups,mu_downs]
colors_dists = [kRed,kBlue,kYellow,kGreen]
colors_types = [0,3,3,-3,-3,-9,-9]
leg_names_dists = ['q#bar{q} #rightarrow t#bar{t}','gg/qg etc. #rightarrow t#bar{t}','simulated background','NTMJ background']
leg_names_types = ['nominal','Afb up','Afb down','d up','d down','#mu up','#mu down']

#array that will eventually hold all of the things
hists = []

#open the output file
outfilep = TFile('template_comparison_plots.root','recreate')
#outfilep = TFile('template_comparison_plots_error_bars.root','recreate')

#make canvases and legends for each plot
canvs = []
legs = []
for nom in noms :
	name = nom.lstrip('__')+'_canv'
	canvs.append(TCanvas(name,name,1100,900))
	legs.append(TLegend(0.62,0.67,0.9,0.9))

#for each of the plot types
for i in range(len(noms)) :
	hists.append([])
	#get the nominal histogram in the first channel
	hists[i].append(infilep.Get(channels[0]+noms[i]+'_x'))
	for k in range(1,len(channels)) :
		newname = channels[k]+noms[i]+'_x'
		hists[i][0].Add(infilep.Get(newname).Clone())
	#also get the other histograms in the first channel
	for j in range(1,len(allnames)) :
		if allnames[j][i]!=None :
			hists[i].append(infilep.Get(channels[0]+allnames[j][i]+'_x').Clone())
			#and for each of those sum up from all of the channels
			for k in range(1,len(channels)) :
				newname = channels[k]+allnames[j][i]+'_x'
				hists[i][len(hists[i])-1].Add(infilep.Get(newname).Clone())

#Set plot attributes and add to legend
for i in range(len(leg_names_dists)) :
	histj = 0
	for j in range(len(leg_names_types)) :
		if allnames[j][i]==None :
			continue
		hists[i][histj].SetLineWidth(4)
		hists[i][histj].SetMarkerStyle(21)
		hists[i][histj].SetFillStyle(0)
		if (j-1)%2==0 :
			hists[i][histj].SetLineStyle(7)
		elif j>0 :
			hists[i][histj].SetLineStyle(3)
		hists[i][histj].SetLineColor(colors_dists[i]+colors_types[j])
		legs[i].AddEntry(hists[i][histj],leg_names_dists[i]+', '+leg_names_types[j],'L')
		histj+=1

#Plot plots
for i in range(len(hists)) :
	canvs[i].cd()
	hists[i][0].Draw('HIST')
#	hists[i][0].Draw('E')
	for j in range(1,len(hists[i])) :
		hists[i][j].Draw('SAME HIST')
#		hists[i][j].Draw('SAME E')
	legs[i].Draw()

#Save plots
outfilep.cd()
for canv in canvs :
	canv.Write()