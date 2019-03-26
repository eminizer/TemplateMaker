#Imports
from ROOT import *
from template import Template
from correct_templates import correct_template

gROOT.SetBatch()

#name of file with templates to smooth
ifn = 'templates_powheg_dynamic_binning_corrected_all.root'
#name of file to put new templates in
ofn = 'templates_powheg_dynamic_binning_smoothed_v12_all.root'
#name of file to put fit canvases in
fcfn = 'smoothing_plots_v12.root'
#name of the algorithm to  use for smoothing
smoothing_algo = 'k3a'

#list of channels to smooth in
cnames = ['t1_muplus_SR',
		  't1_muminus_SR',
		  't1_elplus_SR',
		  't1_elminus_SR',
		  't2_muplus_SR',
		  't2_muminus_SR',
		  't2_elplus_SR',
		  't2_elminus_SR',
		  't3_muplus_SR',
		  't3_muminus_SR',
		  't3_elplus_SR',
		  't3_elminus_SR',
]
#list of processes to smooth for
pnames = ['fqp0',
		  'fqm0',
		  'fq0',
		  'fq1',
		  'fq2',
		  'fg0',
		  'fg1',
		  'fg2',
		  'fg3',
		  'fg4',
		  'fbck',
]
#list of systematics whose templates we want to smooth
snames = ['isr',
		  'fsr',
		  'hdamp',
		  'tune',
		  'cr',
]

#smoothing functions
def smooth_templates_1D_slices(nom,up,dn) :
	#make template objects from the unrolled distributions
	uptemplate = Template(up.GetName(),up.GetTitle(),None)
	dntemplate = Template(dn.GetName(),dn.GetTitle(),None)
	uptemplate.make_from_1D_histo(up)
	dntemplate.make_from_1D_histo(dn)
	#loop over the 2D histos and smooth them
	olduphistos = uptemplate.getHistos(); olddnhistos = dntemplate.getHistos()
	newuphistos = {}; newdnhistos = {}
	for yt,h in olduphistos.items() :
		new_h = h.Clone(); new_h.Reset()
		for i in range(1,h.GetNbinsY()+1) :
			h.GetYaxis().SetRange(i,i)
			projhisto = h.ProjectionX()
			projhisto.Smooth(1,smoothing_algo)
			for j in range(1,projhisto.GetNbinsX()) :
				new_h.SetBinContent(j,i,projhisto.GetBinContent(j))
				new_h.SetBinError(j,i,projhisto.GetBinError(j))
		h.GetYaxis().SetRange()
		newuphistos[yt]=new_h
	for yt,h in olddnhistos.items() :
		new_h = h.Clone(); new_h.Reset()
		for i in range(1,h.GetNbinsY()+1) :
			h.GetYaxis().SetRange(i,i)
			projhisto = h.ProjectionX()
			projhisto.Smooth(1,smoothing_algo)
			for j in range(1,projhisto.GetNbinsX()) :
				new_h.SetBinContent(j,i,projhisto.GetBinContent(j))
				new_h.SetBinError(j,i,projhisto.GetBinError(j))
		h.GetYaxis().SetRange()
		newdnhistos[yt]=new_h
	#reset the template object histos, convert them 1D, correct, and return them
	uptemplate.setHistos(newuphistos)
	dntemplate.setHistos(newdnhistos)
	newup = uptemplate.convertTo1D()
	newdn = dntemplate.convertTo1D()
	return correct_template(nom,newup,newdn)
def smooth_templates_2D_slices(nom,up,dn) :
	#make template objects from the unrolled distributions
	uptemplate = Template(up.GetName(),up.GetTitle(),None)
	dntemplate = Template(dn.GetName(),dn.GetTitle(),None)
	uptemplate.make_from_1D_histo(up)
	dntemplate.make_from_1D_histo(dn)
	#loop over the 2D histos and smooth them
	olduphistos = uptemplate.getHistos(); olddnhistos = dntemplate.getHistos()
	newuphistos = {}; newdnhistos = {}
	for yt,h in olduphistos.items() :
		h.Smooth(1,smoothing_algo)
		newuphistos[yt]=h
	for yt,h in olddnhistos.items() :
		h.Smooth(1,smoothing_algo)
		newdnhistos[yt]=h
	#reset the template object histos, convert them 1D, correct, and return them
	uptemplate.setHistos(newuphistos)
	dntemplate.setHistos(newdnhistos)
	newup = uptemplate.convertTo1D()
	newdn = dntemplate.convertTo1D()
	return correct_template(nom,newup,newdn)
def smooth_templates_polynomial_fit(nom,up,dn) :
	#make template objects from the unrolled distributions
	nomtemplate = Template(nom.GetName(),nom.GetTitle(),None)
	uptemplate = Template(up.GetName(),up.GetTitle(),None)
	dntemplate = Template(dn.GetName(),dn.GetTitle(),None)
	nomtemplate.make_from_1D_histo(nom)
	uptemplate.make_from_1D_histo(up)
	dntemplate.make_from_1D_histo(dn)
	#loop over the 2D histos and smooth them
	oldnomhistos = nomtemplate.getHistos()
	olduphistos = uptemplate.getHistos(); olddnhistos = dntemplate.getHistos()
	newuphistos = {}; newdnhistos = {}
	func = TF1('f','pol3',-1.,1.,)
	for yt,h in olduphistos.items() :
		if h.GetNbinsX()>2 :
			new_h = h.Clone(); new_h.Reset()
			for i in range(1,h.GetNbinsY()+1) :
				fcfp.cd()
				projhisto = h.ProjectionX('_px',i,i,'e')
				nomprojhisto = oldnomhistos[yt].ProjectionX('_nompx',i,i,'e')
				for j in range(1,projhisto.GetNbinsX()+1) :
					projhisto.SetBinContent(j,projhisto.GetBinContent(j)/nomprojhisto.GetBinContent(j))
				projhisto.Fit('f')
				print('{} + {} * x + {} * x**2 + {} * x**3'.format(func.GetParameter(0),func.GetParameter(1),func.GetParameter(2),func.GetParameter(3)))
				#print('i={}, nbinsX={}'.format(i,projhisto.GetNbinsX()))
				for j in range(1,projhisto.GetNbinsX()+1) :
					x1 = projhisto.GetBinLowEdge(j)
					x2 = x1+projhisto.GetBinWidth(j)
					#new_h.SetBinContent(j,i,func.Integral(x1,x2))
					#new_h.SetBinContent(j,i,0.5*(func.Eval(projhisto.GetBinCenter(j))+projhisto.GetBinContent(j)))
					new_h.SetBinContent(j,i,nomprojhisto.GetBinContent(j)*func.Eval(projhisto.GetBinCenter(j)))
					new_h.SetBinError(j,i,projhisto.GetBinError(j))
					#print('i={}, j={}, x1 = {}, x2 = {}, oldcont = {}, cont = {}, err = {}'.format(i,j,x1,x2,projhisto.GetBinContent(j),func.Eval(projhisto.GetBinCenter(j)),projhisto.GetBinError(j)))
				fitcanv = TCanvas(up.GetName()+'_'+str(yt)+'_'+str(i)+'_fit',up.GetName()+'_'+str(yt)+'_'+str(i)+'_fit',1100,900)
				projhisto.Draw("HIST")
				func.Draw("SAME")
				fitcanv.Write()
			#new_h.Scale(h.Integral()/new_h.Integral())
			#new_h.Scale(oldnomhistos[yt].Integral()/new_h.Integral())
			newuphistos[yt]=new_h
		else :
			newuphistos[yt]=h
	for yt,h in olddnhistos.items() :
		if h.GetNbinsX()>2 :
			new_h = h.Clone(); new_h.Reset()
			for i in range(1,h.GetNbinsY()+1) :
				fcfp.cd()
				projhisto = h.ProjectionX('_px',i,i,'e')
				nomprojhisto = oldnomhistos[yt].ProjectionX('_nompx',i,i,'e')
				for j in range(1,projhisto.GetNbinsX()+1) :
					projhisto.SetBinContent(j,projhisto.GetBinContent(j)/nomprojhisto.GetBinContent(j))
				projhisto.Fit('f')
				print('{} + {} * x + {} * x**2 + {} * x**3'.format(func.GetParameter(0),func.GetParameter(1),func.GetParameter(2),func.GetParameter(3)))
				#print('i={}, nbinsX={}'.format(i,projhisto.GetNbinsX()))
				for j in range(1,projhisto.GetNbinsX()+1) :
					x1 = projhisto.GetBinLowEdge(j)
					x2 = x1+projhisto.GetBinWidth(j)
					#new_h.SetBinContent(j,i,func.Integral(x1,x2))
					#new_h.SetBinContent(j,i,0.5*(func.Eval(projhisto.GetBinCenter(j))+projhisto.GetBinContent(j)))
					new_h.SetBinContent(j,i,nomprojhisto.GetBinContent(j)*func.Eval(projhisto.GetBinCenter(j)))
					new_h.SetBinError(j,i,projhisto.GetBinError(j))
					#print('i={}, j={}, x1 = {}, x2 = {}, oldcont = {}, cont = {}, err = {}'.format(i,j,x1,x2,projhisto.GetBinContent(j),func.Eval(projhisto.GetBinCenter(j)),projhisto.GetBinError(j)))
				fitcanv = TCanvas(dn.GetName()+'_'+str(yt)+'_'+str(i)+'_fit',dn.GetName()+'_'+str(yt)+'_'+str(i)+'_fit',1100,900)
				projhisto.Draw("HIST")
				func.Draw("SAME")
				fitcanv.Write()
			#new_h.Scale(h.Integral()/new_h.Integral())
			#new_h.Scale(oldnomhistos[yt].Integral()/new_h.Integral())
			newdnhistos[yt]=new_h
		else :
			newdnhistos[yt]=h
	#reset the template object histos, convert them 1D, correct, and return them
	uptemplate.setHistos(newuphistos)
	dntemplate.setHistos(newdnhistos)
	newup = uptemplate.convertTo1D()
	newdn = dntemplate.convertTo1D()
	return correct_template(nom,newup,newdn)
def smooth_templates(nom,up,dn) :
	return smooth_templates_polynomial_fit(nom,up,dn)
	#return smooth_templates_1D_slices(nom,up,dn)
	#return smooth_templates_2D_slices(nom,up,dn)

#open the old file
ifp = TFile.Open(ifn,'r')

#open the new file
ofp = TFile.Open(ofn,'recreate')

#open the plot file
fcfp = TFile.Open(fcfn,'recreate')

#for each of the systematics
altered_keys = []
all_new_histos = []
for s in snames :
	#pnames = pnames+['fwjets'] if s in ['JES','JER'] else pnames
	#for each channel and process
	for c in cnames :
		for p in  pnames :
			#get the old templates
			nom=ifp.Get(c+'__'+p)
			up=ifp.Get(c+'__'+p+'__'+s+'Up'); altered_keys.append(c+'__'+p+'__'+s+'Up')
			dn=ifp.Get(c+'__'+p+'__'+s+'Down'); altered_keys.append(c+'__'+p+'__'+s+'Down')
			#get the smoothed templates
			print('doing {}Up/Down...'.format(c+'__'+p+'__'+s))
			smoothedup, smootheddown = smooth_templates(nom,up,dn)
			all_new_histos.append(smoothedup); all_new_histos.append(smootheddown)
			#print('smoothedup = {}, smootheddn = {}'.format(smoothedup,smootheddown))
			#set them to write to the new file
			smoothedup.SetDirectory(ofp)
			smootheddown.SetDirectory(ofp)

#copy over all the other stuff
print('copying over unsmoothed templates')
for k in ifp.GetListOfKeys() :
	n=k.GetName()
	if n not in altered_keys :
		obj=k.ReadObj()
		obj.SetDirectory(ofp)

#write the new file
print('writing new file')
ofp.Write()
ofp.Close()

ifp.Close()