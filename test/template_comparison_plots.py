from ROOT import *
import CMS_lumi, tdrstyle

#TDR plot style stuff
gROOT.SetBatch()
tdrstyle.setTDRStyle()
iPeriod = 4 #13TeV iPeriod = 1*(0/1 7 TeV) + 2*(0/1 8 TeV)  + 4*(0/1 13 TeV)
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Preliminary"

#function to get a histogram from any subdirectory of the file
def getHistFromFile(filep,histname) :
	#print 'looking for histogram %s...'%(histname)
	for k in filep.GetListOfKeys() :
		if k.GetName()==histname :
			#print 'returning based on key %s'%(k)
			newhist = filep.Get(histname).Clone()
			newhist.SetTitle('')
			newhist.GetYaxis().SetTitle('Events')
			newhist.SetStats(0)
			return newhist
	#print 'failed to find histogram %s'%(histname)
	return None


#open input file
infilep = TFile('templates_aux.root')

#lepton types to sum over
leptypes = ['muplus','muminus','elplus','elminus']

#histogram name stems
noms = 		['__fqq', 				 '__fgg', 				'__fbck', '__fwjets', '__fqcd']
afb_ups = 	['__fqq__par_Afb__up', 	 None, 					None, 	  None, 	  '__fqcd__par_Afb__up']
afb_downs = ['__fqq__par_Afb__down', None, 					None, 	  None, 	  '__fqcd__par_Afb__down']
d_ups = 	['__fqq__par_d__up', 	 '__fgg__par_d__up', 	None, 	  None, 	  '__fqcd__par_d__up']
d_downs = 	['__fqq__par_d__down', 	 '__fgg__par_d__down', 	None, 	  None, 	  '__fqcd__par_d__down']
mu_ups = 	['__fqq__par_mu__up', 	 '__fgg__par_mu__up', 	None, 	  None, 	  '__fqcd__par_mu__up']
mu_downs = 	['__fqq__par_mu__down',  '__fgg__par_mu__down', None, 	  None, 	  '__fqcd__par_mu__down']
allnames = [noms,afb_ups,afb_downs,d_ups,d_downs,mu_ups,mu_downs]
colors_dists = [kRed,kBlue,kMagenta,kGreen,kYellow]
colors_types = [0,3,3,-3,-3,-9,-9]
leg_names_dists = ['q#bar{q} #rightarrow t#bar{t}','gg/qg etc. #rightarrow t#bar{t}','other top background','W+Jets background','QCD background']
leg_names_types = ['nominal','Afb up','Afb down','d up','d down','#mu up','#mu down']

#dictionary that will eventually hold all of the histos organized by topology and region
hists = {'t1':{'SR':[],'WJets_CR':[]},'t2':{'SR':[],'WJets_CR':[]},'t3':{'SR':[]}}

#open the output file
#outfilep = TFile('template_comparison_plots.root','recreate')
outfilep = TFile('template_comparison_plots_error_bars.root','recreate')

#make canvases, legends, and CMS lumi objects for each plot
canvs 	  = {'t1':{'SR':[],'WJets_CR':[]},'t2':{'SR':[],'WJets_CR':[]},'t3':{'SR':[]}}
legs 	  = {'t1':{'SR':[],'WJets_CR':[]},'t2':{'SR':[],'WJets_CR':[]},'t3':{'SR':[]}}
lumi_objs = {'t1':{'SR':[],'WJets_CR':[]},'t2':{'SR':[],'WJets_CR':[]},'t3':{'SR':[]}}
for top in hists :
	for reg in hists[top] :
		for nom in noms :
			if top=='t1' and nom=='__fqcd' :
				continue
			name = top+'_'+reg+'_'+nom.lstrip('__')+'_canv'
			canvs[top][reg].append(TCanvas(name,name,1100,900))
			legs[top][reg].append(TLegend(0.62,0.67,0.9,0.9))

#for each of the topologies
for top in hists :
	#and each of the regions
	for reg in hists[top] :
		#for each of the plot types
		for i in range(len(noms)) :
			if top=='t1' and noms[i]=='__fqcd' :
				continue
			thishistlist = hists[top][reg]
			thishistlist.append([])
			#get the nominal histogram in the first channel
			thishistlist[i].append(getHistFromFile(infilep,top+'_'+leptypes[0]+'_'+reg+noms[i]+'_x'))
			#sum over the other lepton types
			for k in range(1,len(leptypes)) :
				newname = top+'_'+leptypes[k]+'_'+reg+noms[i]+'_x'
				thishistlist[i][0].Add(getHistFromFile(infilep,newname).Clone())
			#also get the other histograms in the first channel
			for j in range(1,len(allnames)) :
				if allnames[j][i]!=None :
					thishistlist[i].append(getHistFromFile(infilep,top+'_'+leptypes[0]+'_'+reg+allnames[j][i]+'_x').Clone())
					#and for each of those sum up from all of the channels
					for k in range(1,len(leptypes)) :
						newname = top+'_'+leptypes[k]+'_'+reg+allnames[j][i]+'_x'
						thishistlist[i][-1].Add(getHistFromFile(infilep,newname).Clone())

#Set plot attributes and add to legend
for top in hists :
	#and each of the regions
	for reg in hists[top] :
		#for each type of histogram
		for i in range(len(leg_names_dists)) :
			if top=='t1' and noms[i]=='__fqcd' :
				continue
			thishistlist = thishistlist = hists[top][reg]
			histj = 0
			for j in range(len(leg_names_types)) :
				if allnames[j][i]==None :
					continue
				thishistlist[i][histj].SetLineWidth(4)
				thishistlist[i][histj].SetMarkerStyle(21)
				thishistlist[i][histj].SetFillStyle(0)
				if (j-1)%2==0 :
					thishistlist[i][histj].SetLineStyle(7)
				elif j>0 :
					thishistlist[i][histj].SetLineStyle(3)
				thishistlist[i][histj].SetLineColor(colors_dists[i]+colors_types[j])
				legs[top][reg][i].AddEntry(thishistlist[i][histj],leg_names_dists[i]+', '+leg_names_types[j],'L')
				histj+=1

#Plot plots
for top in hists :
	for reg in hists[top] :
		for i in range(len(hists[top][reg])) :
			thismax = hists[top][reg][i][0].GetMaximum()
			for j in range(1,len(hists[top][reg][i])) :
				newmax = hists[top][reg][i][j].GetMaximum()
				if newmax>thismax :
					thismax = newmax
			hists[top][reg][i][0].GetYaxis().SetRangeUser(0.,1.1*thismax)
			canvs[top][reg][i].cd()
		#	hists[top][reg][i][0].Draw('HIST')
			hists[top][reg][i][0].Draw('E')
			for j in range(1,len(hists[top][reg][i])) :
		#		hists[top][reg][i][j].Draw('SAME HIST')
				hists[top][reg][i][j].Draw('SAME E')
			legs[top][reg][i].Draw()
			lumi_objs[top][reg].append(CMS_lumi.CMS_lumi(canvs[top][reg][i], iPeriod, 0))

#Save plots
outfilep.cd()
for top in hists :
	for reg in hists[top] :
		for canv in canvs[top][reg] :
			canv.Write()