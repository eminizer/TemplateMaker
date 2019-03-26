# makes plots to help understand what the different systematic shifts look like

#imports
from ROOT import *
from template import Template
from array import array
from math import sqrt

gROOT.SetBatch()

#name of file to pull 1D templates from
infilename = '../total_template_files/templates_powheg_dynamic_binning_aggregated_v5_all.root'
#name of output file
outfilename = 'systematics_visualizations_dynamic_binning_agg_v5.root'

#dict of systematics to make plots for with the channels/processes to which they apply
t1_mu_cnames = ['t1_muplus_SR','t1_muminus_SR']#,'t1_muplus_WJets_CR','t1_muminus_WJets_CR']
t1_el_cnames = ['t1_elplus_SR','t1_elminus_SR']#,'t1_elplus_WJets_CR','t1_elminus_WJets_CR']
t1_cnames = t1_mu_cnames+t1_el_cnames
t2_mu_cnames = ['t2_muplus_SR','t2_muminus_SR']#,'t2_muplus_WJets_CR','t2_muminus_WJets_CR']
t2_el_cnames = ['t2_elplus_SR','t2_elminus_SR']#,'t2_elplus_WJets_CR','t2_elminus_WJets_CR']
t2_cnames = t2_mu_cnames+t2_el_cnames
t3_mu_cnames = ['t3_muplus_SR','t3_muminus_SR']
t3_el_cnames = ['t3_elplus_SR','t3_elminus_SR']
t3_cnames = t3_mu_cnames+t3_el_cnames
boosted_mu_cnames = t1_mu_cnames+t2_mu_cnames
boosted_el_cnames = t1_el_cnames+t2_el_cnames
mu_cnames = boosted_mu_cnames+t3_mu_cnames
el_cnames = boosted_el_cnames+t3_el_cnames
boosted_cnames = t1_cnames+t2_cnames
all_cnames = boosted_cnames+t3_cnames
ttbar_pnames = ['fqp0','fqm0','fg0','fbck']
wjets_pname = ['fwjets']
all_pnames = ttbar_pnames+wjets_pname
sysnames = {'JES':{'cns':all_cnames,'pns':all_pnames},
			'JER':{'cns':all_cnames,'pns':all_pnames},
			'pileup_weight':{'cns':all_cnames,'pns':all_pnames},
			'mu_trig_eff_weight_b':{'cns':boosted_mu_cnames,'pns':all_pnames},
			'el_trig_eff_weight_b':{'cns':boosted_el_cnames,'pns':all_pnames},
			'mu_trig_eff_weight_r':{'cns':t3_mu_cnames,'pns':all_pnames},
			'el_trig_eff_weight_r':{'cns':t3_el_cnames,'pns':all_pnames},
			'mu_ID_weight':{'cns':mu_cnames,'pns':all_pnames},
			'el_ID_weight':{'cns':el_cnames,'pns':all_pnames},
			'mu_iso_weight':{'cns':mu_cnames,'pns':all_pnames},
			'el_iso_weight':{'cns':el_cnames,'pns':all_pnames},
			'btag_eff_weight_flavb_b':{'cns':boosted_cnames,'pns':all_pnames},
			'btag_eff_weight_flavb_r':{'cns':t3_cnames,'pns':all_pnames},
			'btag_eff_weight_flavc_b':{'cns':boosted_cnames,'pns':all_pnames},
			'btag_eff_weight_flavc_r':{'cns':t3_cnames,'pns':all_pnames},
			'btag_eff_weight_light_b':{'cns':boosted_cnames,'pns':all_pnames},
			'btag_eff_weight_light_r':{'cns':t3_cnames,'pns':all_pnames},
			'ttag_eff_weight_merged':{'cns':t1_cnames,'pns':all_pnames},
			'ttag_eff_weight_semimerged':{'cns':t1_cnames,'pns':all_pnames},
			'ttag_eff_weight_notmerged':{'cns':t1_cnames,'pns':all_pnames},
			'ren_scale_weight':{'cns':all_cnames,'pns':ttbar_pnames},
			'fact_scale_weight':{'cns':all_cnames,'pns':ttbar_pnames},
			'comb_scale_weight':{'cns':all_cnames,'pns':ttbar_pnames},
			'pdfas_weight':{'cns':all_cnames,'pns':ttbar_pnames},
			'B_frag_weight':{'cns':all_cnames,'pns':ttbar_pnames},
			'B_br_weight':{'cns':all_cnames,'pns':ttbar_pnames},
			'top_pt_re_weight':{'cns':all_cnames,'pns':ttbar_pnames},
			'isr':{'cns':all_cnames,'pns':ttbar_pnames},
			'fsr':{'cns':all_cnames,'pns':ttbar_pnames},
			'hdamp':{'cns':all_cnames,'pns':ttbar_pnames},
			'tune':{'cns':all_cnames,'pns':ttbar_pnames},
			'cr':{'cns':all_cnames,'pns':ttbar_pnames},
			}
#dictionary of the colors per channel
ccolors = {'t1_muplus_SR':kOrange+9,
		   't1_muminus_SR':kOrange+9,
		   't1_elplus_SR':kViolet-4,
		   't1_elminus_SR':kViolet-4,
		   't1_muplus_WJets_CR':kOrange-3,
		   't1_muminus_WJets_CR':kOrange-3,
		   't1_elplus_WJets_CR':kViolet-6,
		   't1_elminus_WJets_CR':kViolet-6,
		   't2_muplus_SR':kGreen+3,
		   't2_muminus_SR':kGreen+3,
		   't2_elplus_SR':kMagenta+3,
		   't2_elminus_SR':kMagenta+3,
		   't2_muplus_WJets_CR':kGreen-7,
		   't2_muminus_WJets_CR':kGreen-7,
		   't2_elplus_WJets_CR':kMagenta-7,
		   't2_elminus_WJets_CR':kMagenta-7,
		   't3_muplus_SR':kRed+2,
		   't3_muminus_SR':kRed+2,
		   't3_elplus_SR':kBlue+2,
		   't3_elminus_SR':kBlue+2,
}
#some plot parameters to change them easily
ndevplotbins = 600
devplotlowedge = -1.
devplothighedge = 2.
ncontentbins=166
maxcont=830

#list of all the canvases that are going somewhere so that they don't get deleted because ROOT
all_canvs = []
other_objs = []

#open the input file
ifp = TFile(infilename,'r')

#for each of the systematics
for sn in sysnames.keys() :
	print('Making plots for {} ({} of {})...'.format(sn,sysnames.keys().index(sn)+1,len(sysnames.keys())))
	#declare the histogram stacks for the 1D deviation plots
	upstack = THStack(sn+'_upstack',sn+' up shifts; observed frac. shift up in bin content; # of bins'); other_objs.append(upstack)
	dnstack = THStack(sn+'_dnstack',sn+' down shifts; observed frac. shift down in bin content; # of bins'); other_objs.append(dnstack)
	#loop over the channels
	for cn in sysnames[sn]['cns'] :
		print('	in channel {}'.format(cn))
		#declare the histograms to fill with deviations from this channel
		updevhist=TH1F(sn+'_'+cn+'_updevs',sn+'_'+cn+'_updevs; observed frac. shift up in bin content; # of bins',ndevplotbins,devplotlowedge,devplothighedge)
		updevhist.SetMarkerStyle(21)
		updevhist.SetMarkerColor(ccolors[cn])
		updevhist.SetLineColor(ccolors[cn])
		updevhist.SetFillColor(ccolors[cn])
		other_objs.append(updevhist)
		dndevhist=TH1F(sn+'_'+cn+'_dndevs',sn+'_'+cn+'_dndevs; observed frac. shift down in bin content; # of bins',ndevplotbins,devplotlowedge,devplothighedge)
		dndevhist.SetMarkerStyle(21)
		dndevhist.SetMarkerColor(ccolors[cn])
		dndevhist.SetLineColor(ccolors[cn])
		dndevhist.SetFillColor(ccolors[cn])
		other_objs.append(dndevhist)
		#declare the histograms for the 2D plots 
		upshifts_as_cont  = TH2F(sn+'_'+cn+'_upshifts_as_cont','bin content vs. '+sn+' up shifts in channel '+cn,ndevplotbins/10,devplotlowedge,devplothighedge,ncontentbins,0.,maxcont); other_objs.append(upshifts_as_cont)
		dnshifts_as_cont  = TH2F(sn+'_'+cn+'_dnshifts_as_cont','bin content vs. '+sn+' down shifts in channel '+cn,ndevplotbins/10,devplotlowedge,devplothighedge,ncontentbins,0.,maxcont); other_objs.append(dnshifts_as_cont)
		#histograms for the data comparison plots
		thisdata = ifp.Get(cn+'__data_obs')
		thisdata.SetMarkerStyle(20)
		thismcnom = TH1D(sn+'_'+cn+'_mc_nominal',sn+'_'+cn+'_mc_nominal',thisdata.GetNbinsX(),0,thisdata.GetNbinsX()-1.)
		thismcup  = TH1D(sn+'_'+cn+'_mc_up',sn+'_'+cn+'_mc_up',thisdata.GetNbinsX(),0,thisdata.GetNbinsX()-1.)
		thismcdn  = TH1D(sn+'_'+cn+'_mc_dn',sn+'_'+cn+'_mc_dn',thisdata.GetNbinsX(),0,thisdata.GetNbinsX()-1.)
		thismcnom.SetLineWidth(3); thismcup.SetLineWidth(3); thismcdn.SetLineWidth(3)
		thismcup.SetLineColor(kRed+2); thismcdn.SetLineColor(kBlue+2)
		other_objs.append(thisdata); other_objs.append(thismcnom); other_objs.append(thismcup); other_objs.append(thismcdn)
		##histogram dicts for the projected 1D shift plots (make a dummy template to get the dimensions and to build the plots to project)
		#dummytemplate = Template(cn+'__not_real',cn+'__not_real',None)
		#dummyhistos = dummytemplate.getHistos()
		#nomhists = {}; uphists = {}; dnhists = {}
		#for yt,h in dummyhistos.items() :
		#	newnamestem = sn+'_'+cn
		#	newnameend = '_y='+str(yt[0])+'to'+str(yt[1])
		#	newnomhist = h.Clone(newnamestem+'_nom'+newnameend); newnomhist.SetDirectory(0); other_objs.append(newnomhist); nomhists[yt]=newnomhist
		#	newuphist  = h.Clone(newnamestem+'_up'+newnameend);  newuphist.SetDirectory(0);  other_objs.append(newuphist);  uphists[yt]=newuphist
		#	newdnhist  = h.Clone(newnamestem+'_dn'+newnameend);  newdnhist.SetDirectory(0);  other_objs.append(newdnhist);  dnhists[yt]=newdnhist
		##and build the dummy template from the observed data to loop over its histograms later
		#dummytemplate.make_from_1D_histo(thisdata)
		#and over the processes
		for pn in sysnames[sn]['pns'] :
			#get the nominal and up/down templates
			nom = ifp.Get(cn+'__'+pn)
			up  = ifp.Get(cn+'__'+pn+'__'+sn+'Up')
			dn  = ifp.Get(cn+'__'+pn+'__'+sn+'Down')
			#add to the histograms for the data comparison plots
			thismcnom.Add(nom)
			thismcup.Add(up)
			thismcdn.Add(dn)
			#loop over the bins
			for i in range(1,nom.GetNbinsX()+1) :
				#find the fractional shifts
				nomcont=nom.GetBinContent(i)
				upcont=up.GetBinContent(i)
				dncont=dn.GetBinContent(i)
				upshift = (upcont-nomcont)/nomcont
				dnshift = (dncont-nomcont)/nomcont
				#correct shifts to put over/underflow in the last/first bins
				upshift = min(max(devplotlowedge,upshift),devplothighedge-0.0000001)
				dnshift = min(max(devplotlowedge,dnshift),devplothighedge-0.0000001)
				#fill the 1D deviation histograms
				updevhist.Fill(upshift)
				dndevhist.Fill(dnshift)
				#if abs((upcont-nomcont)/nomcont+1)<0.00005 :
				#	print('-1 up deviation with nom={}, up={}, sys={}, {}__{} bin {}'.format(nomcont,upcont,sn,cn,pn,i))
				#if (dncont-nomcont)/nomcont>3.5 :
				#	print('large dn deviation with nom={}, dn={}, sys={}, {}__{} bin {}'.format(nomcont,dncont,sn,cn,pn,i))
				#fill the 2D plots
				upshifts_as_cont.Fill(upshift,nomcont)
				dnshifts_as_cont.Fill(dnshift,nomcont)
		#	#build full template objects from the nominal and up/down template and add to the histogram dictionaries
		#	sis = str(sysnames.keys().index(sn))
		#	nomtemplate = Template(cn+'__'+pn+'_nom_'+sis,cn+'__'+pn+'_nom_'+sis,None); nomtemplate.make_from_1D_histo(nom)
		#	uptemplate  = Template(cn+'__'+pn+'_up_'+sis,cn+'__'+pn+'_up_'+sis,None);   uptemplate.make_from_1D_histo(up)
		#	dntemplate  = Template(cn+'__'+pn+'_dn_'+sis,cn+'__'+pn+'_dn_'+sis,None);   dntemplate.make_from_1D_histo(dn)
		#	thispnomhistos = nomtemplate.getHistos(); thispuphistos=uptemplate.getHistos(); thispdnhistos=dntemplate.getHistos()
		#	for yt in dummyhistos :
		#		nomhists[yt].Add(thispnomhistos[yt])
		#		uphists[yt].Add(thispuphistos[yt])
		#		dnhists[yt].Add(thispdnhistos[yt])
		#for the data comparison plots add from the processes not affected by the systematics
		for pn in  all_pnames :
			if not pn in sysnames[sn]['pns'] :
				newhist = ifp.Get(cn+'__'+pn)
				thismcnom.Add(newhist)
				thismcup.Add(newhist)
				thismcdn.Add(newhist)
		#		sis = str(sysnames.keys().index(sn))
		#		newtemplate = Template(cn+'__'+pn+'_new_'+sis,cn+'__'+pn+'_new_'+sis,None); newtemplate.make_from_1D_histo(thismcnom)
		#		thispnewhistos = newtemplate.getHistos()
		#		for yt in dummyhistos :
		#			nomhists[yt].Add(thispnewhistos[yt])
		#			uphists[yt].Add(thispnewhistos[yt])
		#			dnhists[yt].Add(thispnewhistos[yt])
		#add the 1D histograms to the stacks
		upstack.Add(updevhist); dnstack.Add(dndevhist)
		#make the canvas for the 2D plots in this channel
		canv2D = TCanvas(sn+'_'+cn+'_2Dplots',sn+'_'+cn+'_2Dplots',1100,900)
		canv2D.Divide(1,2)
		canv2D.cd(1)
		upshifts_as_cont.Draw("COLZ")
		canv2D.cd(2)
		dnshifts_as_cont.Draw("COLZ")
		all_canvs.append(canv2D)
		#make the canvas for the data comparison plots in this channel
		canvdatacomp = TCanvas(sn+'_'+cn+'_data_comp',sn+'_'+cn+'_data_comp',1100,900)
		thismcup.SetMinimum(0.)
		thismcup.Draw('HIST')
		thismcdn.Draw('HISTSAME')
		thismcnom.Draw('HISTSAME')
		thisdata.Draw('PE1 SAME')
		all_canvs.append(canvdatacomp)
		##make and plot the projected 1D shift plots
		#ybins = dummytemplate.getYbins()
		#data_as_y = TH1D(sn+'_'+cn+'_data_as_y','',len(ybins)-1,ybins); other_objs.append(data_as_y)
		#upshifts_as_y = TH1D(sn+'_'+cn+'_upshifts_as_y','',len(ybins)-1,ybins); other_objs.append(upshifts_as_y)
		#dnshifts_as_y = TH1D(sn+'_'+cn+'_dnshifts_as_y','',len(ybins)-1,ybins); other_objs.append(dnshifts_as_y)
		#for i,yt in enumerate(sorted(dummyhistos)) :
		#	erra = array('d',[0.])
		#	cont = dummyhistos[yt].IntegralAndError(0,-1,0,-1,erra)
		#	data_as_y.SetBinContent(i+1,cont); data_as_y.SetBinError(i+1,erra[0])
		#	nomcont = nomhists[yt].IntegralAndError(0,-1,0,-1,erra)
		#	nomerr = erra[0]
		#	upcont = uphists[yt].IntegralAndError(0,-1,0,-1,erra)
		#	uperr = erra[0]
		#	dncont = dnhists[yt].IntegralAndError(0,-1,0,-1,erra)
		#	dnerr = erra[0]
		#	upshifts_as_y.SetBinContent(i+1,(upcont-nomcont)/nomcont)
		#	upshifts_as_y.SetBinError(i+1,sqrt((uperr/nomcont)**2+(upcont*nomerr/(nomcont**2))**2))
		#	dnshifts_as_y.SetBinContent(i+1,(dncont-nomcont)/nomcont)
		#	dnshifts_as_y.SetBinError(i+1,sqrt((dnerr/nomcont)**2+(dncont*nomerr/(nomcont**2))**2))
		#data_as_y.Scale(1./data_as_y.Integral())
		#data_as_y.SetMarkerStyle(20)
		#upshifts_as_y.SetMarkerStyle(22); dnshifts_as_y.SetMarkerStyle(23)
		#upshifts_as_y.SetMarkerColor(kRed+2); dnshifts_as_y.SetMarkerColor(kBlue+2)
		#upshifts_as_y.SetLineColor(kRed+2); dnshifts_as_y.SetLineColor(kBlue+2)
		#upshifts_as_y.SetStats(0)
		#upshifts_as_y.SetMinimum(min(min(upshifts_as_y.GetMinimum(),dnshifts_as_y.GetMinimum()),data_as_y.GetMinimum())-0.02)
		#upshifts_as_y.SetMaximum(max(max(upshifts_as_y.GetMaximum(),dnshifts_as_y.GetMaximum()),data_as_y.GetMaximum())+0.02)
		#shifts_as_y_canv = TCanvas(sn+'_'+cn+'_shifts_as_y',sn+'_'+cn+'_shifts_as_y',1100,900)
		#upshifts_as_y.Draw('PE1X0')
		#dnshifts_as_y.Draw('PE1X0 SAME')
		#data_as_y.Draw('P SAME')
		#all_canvs.append(shifts_as_y_canv)
		##and repeat the above but for the other two dimensions
		#shifts_as_x_canv = TCanvas(sn+'_'+cn+'_shifts_as_x',sn+'_'+cn+'_shifts_as_x',1800,900)
		#shifts_as_z_canv = TCanvas(sn+'_'+cn+'_shifts_as_z',sn+'_'+cn+'_shifts_as_z',1800,900)
		#if len(dummyhistos) == 2 :
		#	shifts_as_x_canv.Divide(2,1); shifts_as_z_canv.Divide(2,1)
		#elif len(dummyhistos) == 3 :
		#	shifts_as_x_canv.Divide(3,1); shifts_as_z_canv.Divide(3,1)
		#elif len(dummyhistos) == 5 :
		#	shifts_as_x_canv.Divide(3,2); shifts_as_z_canv.Divide(3,2)
		#elif len(dummyhistos) == 7 :
		#	shifts_as_x_canv.Divide(4,2); shifts_as_z_canv.Divide(4,2)
		#elif len(dummyhistos) == 10 :
		#	shifts_as_x_canv.Divide(4,3); shifts_as_z_canv.Divide(4,3)
		#all_canvs.append(shifts_as_x_canv); all_canvs.append(shifts_as_z_canv);
		#allbins = dummytemplate.getBins()
		#data_as_xs = []; upshifts_as_xs = []; dnshifts_as_xs = []
		#data_as_zs = []; upshifts_as_zs = []; dnshifts_as_zs = []
		#for i,bins in enumerate(allbins) :
		#	yt = (bins['ylo'],bins['yhi'])
		#	data_as_xs.append(dummyhistos[yt].ProjectionX().Clone(sn+'_'+cn+'_data_as_x_'+str(i))); other_objs.append(data_as_xs[-1])
		#	nom_as_x  = nomhists[yt].ProjectionX().Clone(sn+'_'+cn+'_nom_as_x_'+str(i)); other_objs.append(nom_as_x)
		#	up_as_x  = uphists[yt].ProjectionX().Clone(sn+'_'+cn+'_up_as_x_'+str(i)); other_objs.append(up_as_x)
		#	dn_as_x  = dnhists[yt].ProjectionX().Clone(sn+'_'+cn+'_dn_as_x_'+str(i)); other_objs.append(dn_as_x)
		#	upshifts_as_xs.append(TH1D(sn+'_'+cn+'_upshifts_as_x_'+str(i),'',len(bins['xbins'])-1,bins['xbins'])); other_objs.append(upshifts_as_xs[-1])
		#	dnshifts_as_xs.append(TH1D(sn+'_'+cn+'_dnshifts_as_x_'+str(i),'',len(bins['xbins'])-1,bins['xbins'])); other_objs.append(dnshifts_as_xs[-1])
		#	for j in range(1,data_as_xs[-1].GetNbinsX()+1) :
		#		nomcont = nom_as_x.GetBinContent(j); nomerr=nom_as_x.GetBinError(j)
		#		upcont = up_as_x.GetBinContent(j);   uperr=up_as_x.GetBinError(j)
		#		dncont = dn_as_x.GetBinContent(j);   dnerr=dn_as_x.GetBinError(j)
		#		upshifts_as_xs[-1].SetBinContent(j,(upcont-nomcont)/nomcont)
		#		upshifts_as_xs[-1].SetBinError(j,sqrt((uperr/nomcont)**2+(upcont*nomerr/(nomcont**2))**2))
		#		dnshifts_as_xs[-1].SetBinContent(j,(dncont-nomcont)/nomcont)
		#		dnshifts_as_xs[-1].SetBinError(j,sqrt((dnerr/nomcont)**2+(dncont*nomerr/(nomcont**2))**2))
		#	data_as_xs[-1].Scale(1./data_as_xs[-1].Integral())
		#	data_as_xs[-1].SetMarkerStyle(20)
		#	upshifts_as_xs[-1].SetMarkerStyle(22); dnshifts_as_xs[-1].SetMarkerStyle(23)
		#	upshifts_as_xs[-1].SetMarkerColor(kRed+2); dnshifts_as_xs[-1].SetMarkerColor(kBlue+2)
		#	upshifts_as_xs[-1].SetLineColor(kRed+2); dnshifts_as_xs[-1].SetLineColor(kBlue+2)
		#	upshifts_as_xs[-1].SetStats(0)
		#	upshifts_as_xs[-1].SetMinimum(min(min(upshifts_as_xs[-1].GetMinimum(),dnshifts_as_xs[-1].GetMinimum()),data_as_xs[-1].GetMinimum())-0.02)
		#	upshifts_as_xs[-1].SetMaximum(max(max(upshifts_as_xs[-1].GetMaximum(),dnshifts_as_xs[-1].GetMaximum()),data_as_xs[-1].GetMaximum())+0.02)
		#	shifts_as_x_canv.cd(i+1)
		#	upshifts_as_xs[-1].Draw('PE1X0')
		#	dnshifts_as_xs[-1].Draw('PE1X0 SAME')
		#	data_as_xs[-1].Draw('P SAME')
		#	data_as_zs.append(dummyhistos[yt].ProjectionY().Clone(sn+'_'+cn+'_data_as_z_'+str(i))); other_objs.append(data_as_zs[-1])
		#	nom_as_z  = nomhists[yt].ProjectionY().Clone(sn+'_'+cn+'_nom_as_z_'+str(i)); other_objs.append(nom_as_z)
		#	up_as_z  = uphists[yt].ProjectionY().Clone(sn+'_'+cn+'_up_as_z_'+str(i)); other_objs.append(up_as_z)
		#	dn_as_z  = dnhists[yt].ProjectionY().Clone(sn+'_'+cn+'_dn_as_z_'+str(i)); other_objs.append(dn_as_z)
		#	upshifts_as_zs.append(TH1D(sn+'_'+cn+'_upshifts_as_z_'+str(i),'',len(bins['zbins'])-1,bins['zbins'])); other_objs.append(upshifts_as_zs[-1])
		#	dnshifts_as_zs.append(TH1D(sn+'_'+cn+'_dnshifts_as_z_'+str(i),'',len(bins['zbins'])-1,bins['zbins'])); other_objs.append(dnshifts_as_zs[-1])
		#	for j in range(1,data_as_zs[-1].GetNbinsX()+1) :
		#		nomcont = nom_as_z.GetBinContent(j); nomerr=nom_as_z.GetBinError(j)
		#		upcont = up_as_z.GetBinContent(j);   uperr=up_as_z.GetBinError(j)
		#		dncont = dn_as_z.GetBinContent(j);   dnerr=dn_as_z.GetBinError(j)
		#		upshifts_as_zs[-1].SetBinContent(j,(upcont-nomcont)/nomcont)
		#		upshifts_as_zs[-1].SetBinError(j,sqrt((uperr/nomcont)**2+(upcont*nomerr/(nomcont**2))**2))
		#		dnshifts_as_zs[-1].SetBinContent(j,(dncont-nomcont)/nomcont)
		#		dnshifts_as_zs[-1].SetBinError(j,sqrt((dnerr/nomcont)**2+(dncont*nomerr/(nomcont**2))**2))
		#	data_as_zs[-1].Scale(1./data_as_zs[-1].Integral())
		#	data_as_zs[-1].SetMarkerStyle(20)
		#	upshifts_as_zs[-1].SetMarkerStyle(22); dnshifts_as_zs[-1].SetMarkerStyle(23)
		#	upshifts_as_zs[-1].SetMarkerColor(kRed+2); dnshifts_as_zs[-1].SetMarkerColor(kBlue+2)
		#	upshifts_as_zs[-1].SetLineColor(kRed+2); dnshifts_as_zs[-1].SetLineColor(kBlue+2)
		#	upshifts_as_zs[-1].SetStats(0)
		#	upshifts_as_zs[-1].SetMinimum(min(min(upshifts_as_zs[-1].GetMinimum(),dnshifts_as_zs[-1].GetMinimum()),data_as_zs[-1].GetMinimum())-0.02)
		#	upshifts_as_zs[-1].SetMaximum(max(max(upshifts_as_zs[-1].GetMaximum(),dnshifts_as_zs[-1].GetMaximum()),data_as_zs[-1].GetMaximum())+0.02)
		#	shifts_as_z_canv.cd(i+1)
		#	upshifts_as_zs[-1].Draw('PE1X0')
		#	dnshifts_as_zs[-1].Draw('PE1X0 SAME')
		#	data_as_zs[-1].Draw('P SAME')
	#make the 1D deviation plots
	devcanv = TCanvas(sn+'_devplots',sn+'_devplots',1100,900)
	devcanv.Divide(1,2)
	devcanv.cd(1)
	upstack.Draw()
	devcanv.cd(2)
	dnstack.Draw()
	all_canvs.append(devcanv)

#write canvases to output file
ofp = TFile(outfilename,'recreate')
for c in all_canvs :
	c.Write()
ofp.Close()

#close the input file
ifp.Close()




