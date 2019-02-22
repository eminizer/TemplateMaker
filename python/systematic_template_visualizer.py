# makes plots to help understand what the different systematic shifts look like

#imports
from ROOT import *
from template import Template

gROOT.SetBatch()

#name of file to pull 1D templates from
infilename = '../total_template_files/templates_powheg_corrected_v3_all.root'
#name of output file
outfilename = 'systematics_visualizations.root'

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
		#declare the histograms for the 2D plots (make a dummy template to get the binning)
		dummytemplate = Template(cn+'__not_real',cn+'__not_real',None)
		cstarbinarray = dummytemplate.getXbinArray()
		x_Fbinarray = dummytemplate.getYbinArray()
		Mbinarray = dummytemplate.getZbinArray()
		upshifts_as_cstar = TH2F(sn+'_'+cn+'_upshifts_as_cstar','c* vs. '+sn+' up shifts in channel '+cn,ndevplotbins/10,devplotlowedge,devplothighedge,len(cstarbinarray)-1,cstarbinarray); other_objs.append(upshifts_as_cstar)
		upshifts_as_x_F   = TH2F(sn+'_'+cn+'_upshifts_as_x_F','|x_{F}| vs. '+sn+' up shifts in channel '+cn,ndevplotbins/10,devplotlowedge,devplothighedge,len(x_Fbinarray)-1,x_Fbinarray); other_objs.append(upshifts_as_x_F)
		upshifts_as_M     = TH2F(sn+'_'+cn+'_upshifts_as_M','M vs. '+sn+' up shifts in channel '+cn,ndevplotbins/10,devplotlowedge,devplothighedge,len(Mbinarray)-1,Mbinarray); other_objs.append(upshifts_as_M)
		upshifts_as_cont  = TH2F(sn+'_'+cn+'_upshifts_as_cont','bin content vs. '+sn+' up shifts in channel '+cn,ndevplotbins/10,devplotlowedge,devplothighedge,ncontentbins,0.,maxcont); other_objs.append(upshifts_as_cont)
		dnshifts_as_cstar = TH2F(sn+'_'+cn+'_dnshifts_as_cstar','c* vs. '+sn+' down shifts in channel '+cn,ndevplotbins/10,devplotlowedge,devplothighedge,len(cstarbinarray)-1,cstarbinarray); other_objs.append(dnshifts_as_cstar)
		dnshifts_as_x_F   = TH2F(sn+'_'+cn+'_dnshifts_as_x_F','|x_{F}| vs. '+sn+' down shifts in channel '+cn,ndevplotbins/10,devplotlowedge,devplothighedge,len(x_Fbinarray)-1,x_Fbinarray); other_objs.append(dnshifts_as_x_F)
		dnshifts_as_M     = TH2F(sn+'_'+cn+'_dnshifts_as_M','M vs. '+sn+' down shifts in channel '+cn,ndevplotbins/10,devplotlowedge,devplothighedge,len(Mbinarray)-1,Mbinarray); other_objs.append(dnshifts_as_M)
		dnshifts_as_cont  = TH2F(sn+'_'+cn+'_dnshifts_as_cont','bin content vs. '+sn+' down shifts in channel '+cn,ndevplotbins/10,devplotlowedge,devplothighedge,ncontentbins,0.,maxcont); other_objs.append(dnshifts_as_cont)
		#and over the processes
		for pn in sysnames[sn]['pns'] :
			#get the nominal and up/down templates
			nom = ifp.Get(cn+'__'+pn)
			up  = ifp.Get(cn+'__'+pn+'__'+sn+'Up')
			dn  = ifp.Get(cn+'__'+pn+'__'+sn+'Down')
			#build a full template object from the nominal template
			nomtemplate = Template(cn+'__'+pn+'_'+str(sysnames.keys().index(sn)),cn+'__'+pn+'_'+str(sysnames.keys().index(sn)),None)
			nomtemplate.make_from_1D_histo(nom)
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
				thiscstar, thisx_F, thisM = nomtemplate.getCoordsFrom1DBin(i)
				upshifts_as_cstar.Fill(upshift,thiscstar)
				upshifts_as_x_F.Fill(upshift,thisx_F)
				upshifts_as_M.Fill(upshift,thisM)
				upshifts_as_cont.Fill(upshift,nomcont)
				dnshifts_as_cstar.Fill(dnshift,thiscstar)
				dnshifts_as_x_F.Fill(dnshift,thisx_F)
				dnshifts_as_M.Fill(dnshift,thisM)
				dnshifts_as_cont.Fill(dnshift,nomcont)
		#add the 1D histograms to the stacks
		upstack.Add(updevhist); dnstack.Add(dndevhist)
		#make the canvas for the 2D plots in this channel
		canv2D = TCanvas(sn+'_'+cn+'_2Dplots',sn+'_'+cn+'_2Dplots',1800,900)
		canv2D.Divide(4,2)
		canv2D.cd(1)
		upshifts_as_cstar.Draw("COLZ")
		canv2D.cd(2)
		upshifts_as_x_F.Draw("COLZ")
		canv2D.cd(3)
		upshifts_as_M.Draw("COLZ")
		canv2D.cd(4)
		upshifts_as_cont.Draw("COLZ")
		canv2D.cd(5)
		dnshifts_as_cstar.Draw("COLZ")
		canv2D.cd(6)
		dnshifts_as_x_F.Draw("COLZ")
		canv2D.cd(7)
		dnshifts_as_M.Draw("COLZ")
		canv2D.cd(8)
		dnshifts_as_cont.Draw("COLZ")
		all_canvs.append(canv2D)
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




