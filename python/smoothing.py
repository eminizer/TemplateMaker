#Imports
from math import sqrt
from optparse import OptionParser
import os
import random
from ROOT import *
from template import Template
from lowess import lowess

gROOT.SetBatch()
gErrorIgnoreLevel=kError

##########								Parser Options								##########

parser = OptionParser()
#Input template filepath?
parser.add_option('--infile', type='string', action='store', dest='infilepath',
	help='Which template file is holding the fit histograms? Required.')
#run name? (used in filenames)
parser.add_option('--name', type='string', action='store', default='smoothed', dest='name',
	help='Name of this smoothing procedure.')
#channel/process/ masks
parser.add_option('--channels', type='string', action='store', default='*', dest='channels',
	help='Channels to smooth in separated by double underscores')
parser.add_option('--processes', type='string', action='store', default='*', dest='processes',
	help='Processes to smooth separated by double underscores')
parser.add_option('--nuisances', type='string', action='store', default='*', dest='nuisances',
	help='Nuisances whose corresponding templates should be smoothed separated by double underscores')
#LOWESS options
parser.add_option('--deg', type='int', action='store', default=1, dest='deg',	   	  
	help='Degree of polynomial for LOWESS smoothing function (0=constant, 1=linear, etc.)')
parser.add_option('-K','--kernel',  type='choice', action='store', default='epanechnikov', dest='kernel', choices=['epanechnikov', 'tri_cube', 'bi_square','k5a','k5b','k3a'],
	help='Kernel for LOWESS function ("epanechnikov", "tri_cube", "bi_square", "k5a", "k5b" or "k3a")')
parser.add_option('--robust', action='store_true', dest='robust',
	help='Apply robustification procedure?') 

(options, args) = parser.parse_args()

#process some options into lists
channels  = [x.lower() for x in options.channels.split('__')] if options.channels!='*' else '*'
processes = [x.lower() for x in options.processes.split('__')] if options.processes!='*' else '*'
nuisances = [x.lower() for x in options.nuisances.split('__')] if options.nuisances!='*' else '*'

##########						Function to Remove Zeroed Bins						##########

def removeZeroBins(histo,nomhisto=None) :
	nbins = histo.GetNbinsX()
	#loop over all the bins
	for i in range(1,nbins+1) :
		#if the content of the histogram is <=zero or the nominal content is <=zero/0.00001
		if histo.GetBinContent(i)<=0. or (nomhisto!=None and (nomhisto.GetBinContent(i)<=0. or nomhisto.GetBinContent(i)==0.00001)) :
			#if it's the nominal histogram, just set the content to 0.00001
			if nomhisto==None : 
				histo.SetBinContent(i,0.00001)
			else :
				#otherwise get the content and error shift from the nominal histogram's content
				cont = 0.00001 if nomhisto.GetBinContent(i)<=0. else nomhisto.GetBinContent(i)
				shift = 0.1*cont*(random.random()-0.5)
				if histo.GetName().endswith('Up') :
					histo.SetBinContent(i,cont+shift)
				elif histo.GetName().endswith('Down') :
					histo.SetBinContent(i,cont-shift)
	return histo

##########						Template Smoothing Function						##########

def smoothTemplate(old_histo1D,name,deg,kernel,robust,xrbin=1,yrbin=1) :
	#first make a new Template object built from the old 1D histogram
	temp = Template(name,name,None)
	if not (xrbin==1 and yrbin==1) :
		temp = rebin2DSlicesFromTemplate(temp,xrbin,yrbin)
	temp.make_from_1D_histo(old_histo1D)
	#get this template in three dimensions
	histo3D = temp.getHisto3D()
	#make a list of x-y projected 2D histograms, one for each Z bin
	proj_histos = {}
	nxbins = histo3D.GetXaxis().GetNbins()
	nybins = histo3D.GetYaxis().GetNbins()
	nzbins = histo3D.GetZaxis().GetNbins()
	for i in range(1,nzbins+1) :
		histo3D.GetZaxis().SetRange(i,i)
		#projstring='yxe' if i==1 else 'yx'
		projstring='yx'
		newprojhisto = histo3D.Project3D(projstring)
		proj_histos[i] = newprojhisto.Clone(newprojhisto.GetName()+'_'+str(i))
		#print 'projected nxbins=%d=?%d, nybins=%d=?%d'%(proj_histos[-1].GetXaxis().GetNbins(),nxbins,proj_histos[-1].GetYaxis().GetNbins(),nybins) #DEBUG
	histo3D.GetZaxis().SetRange(-1,nzbins+5)
	#smooth each projection
	smoothed_proj_histos = {}	
	for zbin,proj_histo in proj_histos.iteritems() :
		#make the numpy objects needed for the LOWESS algorithm
		#x=np.array([nxbins*[0.],nybins*[0.]]) #range of the problem: cstar/|x_F| coordinates at bin centers
		#y=np.zeros()
		#x0=np.zeros()
		##fill the objects from the projected histograms
		#for i in range(nglobalbins) :
		#	if not proj_histo.IsBinOverflow(i) and not proj_histo.IsBinUnderflow(i) :
		#		#put stuff here
		proj_histo.Smooth(1,kernel)
		smoothed_proj_histos[zbin]=proj_histo
	#build the total template back up from its projections
	histo3D.Reset()
	for k in range(1,nzbins+1) :
		thisprojhisto=smoothed_proj_histos[k]
		for i in range(1,nxbins+1) :
			for j in range(1,nybins+1) :
				histo3D.SetBinContent(histo3D.GetBin(i,j,k),thisprojhisto.GetBinContent(thisprojhisto.GetBin(i,j,k)))
	temp.setHistos([histo3D,histo3D.ProjectionX(),histo3D.ProjectionY(),histo3D.ProjectionZ()])
	return temp

##########						Comparison Plot Function						##########

def makeComparisonPlots(name,old_histo1D,new_template,nom_histo1D,afp) :
	#first we need to make template objects from the old (unsmoothed) and nominal 1D histograms
	old_template = Template(old_histo1D.GetName(),old_histo1D.GetName(),None,name.split('_')[-1])
	old_template.make_from_1D_histo(old_histo1D)
	nom_template = Template(nom_histo1D.GetName()+'_t',nom_histo1D.GetName()+'_t',None,name.split('_')[-1])
	nom_template.make_from_1D_histo(nom_histo1D)
	#cd to the auxiliary output file
	afp.cd()
	#get the histograms and set their attributes
	old_x_proj = old_template.getHistoX()
	old_y_proj = old_template.getHistoY()
	new_x_proj = new_template.getHistoX()
	new_y_proj = new_template.getHistoY()
	nom_x_proj = nom_template.getHistoX()
	nom_y_proj = nom_template.getHistoY()
	new_histo1D,garbage = new_template.convertTo1D()
	old_x_proj.SetLineWidth(4); old_x_proj.SetLineColor(kRed)
	old_y_proj.SetLineWidth(4); old_y_proj.SetLineColor(kRed)
	new_x_proj.SetLineWidth(4); new_x_proj.SetLineColor(kBlue)
	new_y_proj.SetLineWidth(4); new_y_proj.SetLineColor(kBlue)
	nom_x_proj.SetLineWidth(4); nom_x_proj.SetLineColor(kBlack)
	nom_y_proj.SetLineWidth(4); nom_y_proj.SetLineColor(kBlack)
	old_histo1D.SetLineWidth(4); old_histo1D.SetLineColor(kRed)
	new_histo1D.SetLineWidth(4); new_histo1D.SetLineColor(kBlue)
	nom_histo1D.SetLineWidth(4); nom_histo1D.SetLineColor(kBlack)
	nom_x_proj.GetYaxis().SetRangeUser(0.,1.02*max(nom_x_proj.GetMaximum(),old_x_proj.GetMaximum(),new_x_proj.GetMaximum()))
	nom_y_proj.GetYaxis().SetRangeUser(0.,1.02*max(nom_y_proj.GetMaximum(),old_y_proj.GetMaximum(),new_y_proj.GetMaximum()))
	nom_histo1D.GetYaxis().SetRangeUser(0.,1.02*max(nom_histo1D.GetMaximum(),old_histo1D.GetMaximum(),new_histo1D.GetMaximum()))
	#make the 1D ratio histogram and line at one
	ratio_histo = old_histo1D.Clone(old_histo1D.GetName()+'__ratios')
	ratio_histo.Reset()
	ratio_histo.SetLineWidth(1); ratio_histo.SetLineColor(kBlack)
	for i in range(1,ratio_histo.GetXaxis().GetNbins()+1) :
		oc = old_histo1D.GetBinContent(i)
		nc = new_histo1D.GetBinContent(i)
		oe = old_histo1D.GetBinError(i)
		ne = new_histo1D.GetBinError(i)
		if oc>0. :
			ratio_histo.SetBinContent(i,nc/oc)
			if nc>0. :
				ratio_histo.SetBinError(i,(nc/oc)*sqrt((oe/oc)**2+(ne/nc)**2))
	oneline = TLine(0.,1.,ratio_histo.GetXaxis().GetNbins(),1.)
	oneline.SetLineWidth(3); oneline.SetLineColor(kBlack)
	#create the canvas to hold the plots
	canv = TCanvas(name,name+' smoothing',1200,1200)
	canv.Divide(2,2)
	#make the x projection comparison plot
	canv.cd(1)
	nom_x_proj.Draw("HIST E")
	old_x_proj.Draw("HIST E SAME")
	new_x_proj.Draw("HIST E SAME")
	#make the y projection comparison plot
	canv.cd(2)
	nom_y_proj.Draw("HIST E")
	old_y_proj.Draw("HIST E SAME")
	new_y_proj.Draw("HIST E SAME")
	#make the 1D unrolled comparison plot
	canv.cd(3)
	nom_histo1D.Draw("HIST E")
	old_histo1D.Draw("HIST E SAME")
	new_histo1D.Draw("HIST E SAME")
	#make the ratio comparison plot
	canv.cd(4)
	ratio_histo.Draw("PE1")
	oneline.Draw("SAME")
	#write the canvas to the file
	canv.Write()
	return canv

##########						Rebinning Functions						##########

#given a 1D histogram, rebins and returns the resulting 1D histogram
def rebin2DSlicesFrom1D(hist,xrbin=5,yrbin=5) :
	#make a template object for the old 1D histogram
	name_fine = hist.GetName()+'_fine'; name_coarse = hist.GetName()+'_coarse'
	temp_fine = Template(name_fine,name_fine,None,binning='fine')
	temp_fine.make_from_1D_histo(hist)
	#get this template in three dimensions
	histo3D_fine = temp_fine.getHisto3D()
	#make a list of x-y projected 2D histograms, one for each Z bin
	proj_histos = {}
	nzbins = histo3D_fine.GetZaxis().GetNbins()
	for i in range(1,nzbins+1) :
		histo3D_fine.GetZaxis().SetRange(i,i)
		#projstring='yxe' if i==1 else 'yx'
		projstring='yx'
		newprojhisto = histo3D_fine.Project3D(projstring)
		proj_histos[i] = newprojhisto.Clone(newprojhisto.GetName()+'_'+str(i))
	histo3D_fine.GetZaxis().SetRange(-1,nzbins+5)
	#rebin each projection
	rebinned_proj_histos = {}	
	for zbin,proj_histo in proj_histos.iteritems() :
		proj_histo.Rebin2D(xrbin,yrbin)
		rebinned_proj_histos[zbin]=proj_histo
	#build the new, coarsely-binned total template back up from its projections
	temp_coarse = Template(name_coarse,name_coarse,None,binning='coarse')
	histo3D_coarse = temp_coarse.getHisto3D()
	nxbins = histo3D_coarse.GetXaxis().GetNbins()
	nybins = histo3D_coarse.GetYaxis().GetNbins()
	for k in range(1,nzbins+1) :
		thisprojhisto=rebinned_proj_histos[k]
		for i in range(1,nxbins+1) :
			for j in range(1,nybins+1) :
				histo3D_coarse.SetBinContent(histo3D_coarse.GetBin(i,j,k),thisprojhisto.GetBinContent(thisprojhisto.GetBin(i,j,k)))
	temp_coarse.setHistos([histo3D_coarse,histo3D_coarse.ProjectionX(),histo3D_coarse.ProjectionY(),histo3D_coarse.ProjectionZ()])
	returntemplate,garbage = temp_coarse.convertTo1D()
	return returntemplate

#given a finely-binned template, rebins and returns the resulting 1D histogram
def rebin2DSlicesFromTemplate(temp_fine,xrbin=5,yrbin=5) :
	#get this template in three dimensions
	histo3D_fine = temp_fine.getHisto3D()
	#make a list of x-y projected 2D histograms, one for each Z bin
	proj_histos = {}
	nzbins = histo3D_fine.GetZaxis().GetNbins()
	for i in range(1,nzbins+1) :
		histo3D_fine.GetZaxis().SetRange(i,i)
		#projstring='yxe' if i==1 else 'yx'
		projstring='yx'
		newprojhisto = histo3D_fine.Project3D(projstring)
		proj_histos[i] = newprojhisto.Clone(newprojhisto.GetName()+'_'+str(i))
	histo3D_fine.GetZaxis().SetRange(-1,nzbins+5)
	#rebin each projection
	rebinned_proj_histos = {}	
	for zbin,proj_histo in proj_histos.iteritems() :
		proj_histo.Rebin2D(xrbin,yrbin)
		rebinned_proj_histos[zbin]=proj_histo
	#build the new, coarsely-binned total template back up from its projections
	name_coarse = temp_fine.getName()+'_coarse'
	temp_coarse = Template(name_coarse,name_coarse,None,binning='coarse')
	histo3D_coarse = temp_coarse.getHisto3D()
	nxbins = histo3D_coarse.GetXaxis().GetNbins()
	nybins = histo3D_coarse.GetYaxis().GetNbins()
	for k in range(1,nzbins+1) :
		thisprojhisto=rebinned_proj_histos[k]
		for i in range(1,nxbins+1) :
			for j in range(1,nybins+1) :
				histo3D_coarse.SetBinContent(histo3D_coarse.GetBin(i,j,k),thisprojhisto.GetBinContent(thisprojhisto.GetBin(i,j,k)))
	temp_coarse.setHistos([histo3D_coarse,histo3D_coarse.ProjectionX(),histo3D_coarse.ProjectionY(),histo3D_coarse.ProjectionZ()])
	return temp_coarse

##########								Main Script								##########

#open old file
if not os.path.isfile(options.infilepath) :
	print 'ERROR: filepath %s cannot be found'%options.infilepath
	os._exit()
input_file = TFile.Open(options.infilepath,'r')

#open new file to cd to its currently empty directory
output_tfile_name = os.path.split(options.infilepath)[1].split('.root')[0]+'_'+options.name+'.root'
garbage_file_name = os.path.split(options.infilepath)[1].split('.root')[0]+'_'+options.name+'_smoothing_garbage.root'
output_aux_file_name = os.path.split(options.infilepath)[1].split('.root')[0]+'_'+options.name+'_smoothing_aux.root'
print 'new templates will be in file %s'%(output_tfile_name)
output_tfile = TFile.Open(output_tfile_name,'recreate')
garbage_file = TFile.Open(garbage_file_name,'recreate')
output_auxfile = TFile.Open(output_aux_file_name,'recreate')
garbage_file.cd()

#iterate through all old histogram objects in the input file 
rebinned_nominal_histos = {}
all_old_histos={}; all_new_histos={}; all_new_templates={}; all_canvs = {}
rebinned_template_histos = {}
keylist = input_file.GetListOfKeys()
for k in keylist :
	#get the histogram name
	n = k.GetName()
	#find the channel, process, and nuisance of the histogram
	nsplit = n.split('__')
	thischan = nsplit[0]; thisproc=''; thisnuis=''; thisud=''
	if len(nsplit) > 1 :
		thisproc=nsplit[1]
		if len(nsplit) > 2 :
			thisnuis=nsplit[2].rstrip('Up').rstrip('Down')
			if nsplit[2].endswith('Up') :
				thisud='Up'
			elif nsplit[2].endswith('Down') :
				thisud='Down'
			else :
				print 'WARNING: Cannot identify whether the %s nuisance is wiggled up or down for this template; nsplit[2]=%s!!!'%(thisnuis,nsplit[2])
				continue
	else :
		print 'WARNING: Cannot identify channel/process/nuisance from key %s, skipping this object entirely!!!'%(n)
		continue
	#get the histogram itself
	all_old_histos[n]=k.ReadObj()
	real_template=None
	#check if this histogram is one that should be smoothed based on the masks
	if (channels=='*' or thischan.lower() in channels) and (processes=='*' or thisproc.lower() in processes) and (nuisances=='*' or (thisnuis!='' and thisnuis.lower() in nuisances)) :
		print 'Smoothing template %s...'%(n)
		#get back the new smoothed template object and its 1D histogram
		all_new_templates[n] = smoothTemplate(all_old_histos[n],n,options.deg,options.kernel,options.robust)
		thisnewhisto,garbage = all_new_templates[n].convertTo1D()
		all_new_histos[n] = thisnewhisto
		#make the comparison plot and save it to the auxiliary file
		nominal_histo = (input_file.Get(thischan+'__'+thisproc)).Clone(thischan+'__'+thisproc+'__'+thisnuis+'NominalFor'+thisud)
		rebinned_nominal_histo = rebin2DSlicesFrom1D(nominal_histo)
		all_canvs[n+'_fine'] = makeComparisonPlots(n+'_fine',all_old_histos[n],all_new_templates[n],nominal_histo,output_auxfile)
		all_canvs[n+'_coarse'] = makeComparisonPlots(n+'_coarse',rebin2DSlicesFrom1D(all_old_histos[n]),rebin2DSlicesFromTemplate(all_new_templates[n]),rebinned_nominal_histo,output_auxfile)
		garbage_file.cd()
		print '	Done.'	
		#add renamed and rebinned template to the list of final histograms
		real_template = rebin2DSlicesFrom1D(all_new_histos[n])
		real_template.SetName(n)
	else :
		real_template = rebin2DSlicesFrom1D(all_old_histos[n])
		real_template.SetName(n)
	rebinned_template_histos[n]=real_template
	#if it was a nominal histogram, add it to the nominal histogram dictionary to fix the zeroed bins later
	if thisnuis=='' :# and thisproc!='data_obs' :
		rebinned_nominal_histos[thischan+'__'+thisproc] = real_template

#next we have to handle the bins that are zeroed
final_real_template_histos = []
for n in rebinned_template_histos :
	nsplit = n.split('__')
	thischan = nsplit[0]; thisproc=''; thisnuis=''
	if len(nsplit) > 1 :
		thisproc=nsplit[1]
		if len(nsplit) > 2 :
			thisnuis=nsplit[2].rstrip('Up').rstrip('Down')
	if thisnuis=='' : #if it's the nominal template, just correct the zero bins
		final_real_template_histos.append(removeZeroBins(rebinned_template_histos[n]))
	else : #otherwise we also need to send the nominal template for this process/channel
		final_real_template_histos.append(removeZeroBins(rebinned_template_histos[n],rebinned_nominal_histos[thischan+'__'+thisproc]))		

#Set the final template histogram locations
for h in final_real_template_histos :
	h.SetDirectory(output_tfile)
#write and close the new files
output_tfile.Write()
output_tfile.Close()
garbage_file.Write()
garbage_file.Close()
output_auxfile.Write()
output_auxfile.Close()
#close the input file
input_file.Close()
#delete the garbage file
os.system('rm -rf '+garbage_file_name)