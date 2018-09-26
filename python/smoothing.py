#Imports
from math import sqrt
from optparse import OptionParser
import os
import numpy as np
from ROOT import *
from template import Template
from lowess import lowess

gROOT.SetBatch()

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

##########						Template Smoothing Function						##########

def smoothTemplate(old_histo1D,name,deg,kernel,robust) :
	#first make a new Template object built from the old 1D histogram
	temp = Template(name,name,None)
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
		##find the dimensionality of the problem
		#nglobalbins = proj_histo.GetSize()
		#ntotalpoints = 0
		#for i in range(nglobalbins) :
		#	if not proj_histo.IsBinOverflow(i) and not proj_histo.IsBinUnderflow(i) :
		#		ntotalpoints+=1
		##make the numpy objects needed for the LOWESS algorithm
		#x=np.zeros((2,len(bincontents)))
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

def makeComparisonPlots(old_histo1D,new_template,nom_histo1D,afp) :
	#first we need to make template objects from the old (unsmoothed) and nominal 1D histograms
	old_template = Template(old_histo1D.GetName(),old_histo1D.GetName(),None)
	old_template.make_from_1D_histo(old_histo1D)
	nom_template = Template(nom_histo1D.GetName(),nom_histo1D.GetName(),None)
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
	new_histo1D = new_template.convertTo1D()
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
	canv = TCanvas(n,n+' smoothing',1200,1200)
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

##########								Main Script								##########

#open old file
if not os.path.isfile(options.infilepath) :
	print 'ERROR: filepath %s cannot be found'%options.infilepath
	os._exit()
input_file = TFile.Open(options.infilepath,'r')

#open new file to cd to its currently empty directory
output_tfile_name = os.path.split(options.infilepath)[1].split('.root')[0]+'_'+options.name+'.root'
output_aux_file_name = os.path.split(options.infilepath)[1].split('.root')[0]+'_'+options.name+'_smoothing_aux.root'
print 'new templates will be in file %s'%(output_tfile_name)
output_tfile = TFile.Open(output_tfile_name,'recreate')
output_auxfile = TFile.Open(output_aux_file_name,'recreate')
output_tfile.cd()

#iterate through all old histogram objects in the input file
all_old_histos = {}; all_new_histos={}; all_new_templates={}; all_canvs = {}
keylist = input_file.GetListOfKeys()
for k in keylist :
	#get the histogram name
	n = k.GetName()
	#find the channel, process, and nuisance of the histogram
	nsplit = n.split('__')
	thischan = nsplit[0]; thisproc=''; thisnuis=''
	if len(nsplit) > 1 :
		thisproc=nsplit[1]
		if len(nsplit) > 2 :
			thisnuis=nsplit[2].rstrip('Up').rstrip('Down')
	else :
		print 'WARNING: Cannot identify channel/process/nuisance from key %s, skipping this object entirely!!!'%(n)
		continue
	#get the histogram itself
	all_old_histos[n]=k.ReadObj()
	#check if this histogram is one that should be smoothed based on the masks
	if (channels=='*' or thischan.lower() in channels) and (processes=='*' or thisproc.lower() in processes) and (nuisances=='*' or (thisnuis!='' and thisnuis.lower() in nuisances)) :
		print 'Smoothing template %s...'%(n)
		#get back the new smoothed template object and its 1D histogram
		all_new_templates[n] = smoothTemplate(all_old_histos[n],n,options.deg,options.kernel,options.robust)
		all_new_histos[n] = all_new_templates[n].convertTo1D()
		#make the comparison plot and save it to the auxiliary file
		all_canvs[n] = makeComparisonPlots(all_old_histos[n],all_new_templates[n],input_file.Get(thischan+'__'+thisproc),output_auxfile)
		output_tfile.cd()
		print '	Done.'	
	#if this histogram is one of the new ones, put its new histogram in the output file, otherwise just put the old one	
	all_new_histos[n].SetDirectory(output_tfile) if n in all_new_histos else all_old_histos[n].SetDirectory(output_tfile)
#write and close the new files
output_tfile.Write()
output_tfile.Close()
output_auxfile.Write()
output_auxfile.Close()
#close the input file
input_file.Close()