#Imports
from math import sqrt
import os
import random
from ROOT import *
from template import Template

gROOT.SetBatch()

#input/output total template file names
infilepath = 'templates_powheg_iter_5_smoothed_all.root'
outfilepath = 'templates_powheg_iter_5_regularized_all.root'

#process and systematic masks: which histograms are we allowed to modify? 
process_masks = ['fqp0','fqm0','fq0','fq1','fq2','fg0','fg1','fg2','fg3','fg4','fbck','fwjets'] #note no QCD
sys_masks = ['JES','JER'] #right now only worried about jet energy scale/resolution

#open the input file
infilep = TFile.Open(infilepath,'r')

#make the dictionary of histogram sets by looking at every nuisance-wiggled histogram in the file
histo_sets = {}
keylist = infilep.GetListOfKeys()
for k in keylist :
	#get the histogram name
	n = k.GetName()
	#find the channel, process, and nuisance of the histogram
	nsplit = n.split('__')
	thischan = nsplit[0]; thisproc=''; thisnuis=''; thisud=''
	thisproc=nsplit[1]
	if len(nsplit)>2 :
		thisnuis=nsplit[2].rstrip('Up').rstrip('Down')
		if nsplit[2].endswith('Up') :
			thisud='Up'
		elif nsplit[2].endswith('Down') :
			thisud='Down'
	#dictionary key is the channel__process
	dictkey = thischan+'__'+thisproc
	if not dictkey in histo_sets :
		histo_sets[dictkey]={}
	#put the histogram itself in the dictionary
	#if this is a nominal histogram just put it in under 'nominal'
	if thisnuis=='' :
		histo_sets[dictkey]['nominal']=k.ReadObj()
	else :
		#otherwise second key is thisnuisance
		if not thisnuis in histo_sets[dictkey] :
			histo_sets[dictkey][thisnuis]={}
		#and third key is 'Up'/'Down'
		histo_sets[dictkey][thisnuis][thisud]=k.ReadObj()

#open new file
print 'new templates will be in file %s'%(outfilepath)
outfilep = TFile.Open(outfilepath,'recreate')

#get the average/std. dev. of the total up/down fluctuations in each bin for each channel__process and systematic
flucts={}
for dk in histo_sets :
	if not dk in flucts :
		flucts[dk]={}
	#the nominal histogram
	nomhisto = histo_sets[dk]['nominal']
	nbins = nomhisto.GetNbinsX()
	#for each systematic:
	for s in histo_sets[dk] :
		if s=='nominal' : 
			continue
		#look in each bin to add to the average and std. dev.
		thisavg=0.; thisavg2=0.
		for bin in range(1,nbins+1) :
			fluct = (histo_sets[dk][s]['Up'].GetBinContent(bin)-histo_sets[dk][s]['Down'].GetBinContent(bin))/nomhisto.GetBinContent(bin)
			thisavg+=fluct
			thisavg2+=fluct**2
		thisavg/=nbins; thisavg2/=nbins
		flucts[dk][s]={'avg':thisavg,'stddev':sqrt(abs(thisavg2-thisavg**2))}

#loop back through again to make the modifications
for dk in histo_sets :
	print 'checking '+dk+'...'
	#see if we want to skip this process
	if dk.split('__')[-1] not in process_masks :
		#print '	this process is skipped.'
		#just get ready to copy over the original histograms 
		for s in histo_sets[dk] :
			if s=='nominal' :
				histo_sets[dk][s].SetDirectory(outfilep)
			else :
				histo_sets[dk][s]['Up'].SetDirectory(outfilep)
				histo_sets[dk][s]['Down'].SetDirectory(outfilep)
		#skip to the next channel/process
		continue
	#otherwise start looking through the systematics
	nomhisto = histo_sets[dk]['nominal']
	nomhisto.SetDirectory(outfilep) #get ready to copy over the nominal histogram
	nbins = nomhisto.GetNbinsX()
	for s in histo_sets[dk] :
		if s=='nominal' :
			continue
		print '	systematic %s (avg fluctuation = %.6f +/- %.6f)...'%(s,flucts[dk][s]['avg'],flucts[dk][s]['stddev'])
		#see if we want to modify histograms for this systematic
		if s in sys_masks :
			#loop through the bins
			for bin in range(1,nbins+1) :
				#get this bin's fluctuation 
				oldupcont=histo_sets[dk][s]['Up'].GetBinContent(bin)
				olddowncont=histo_sets[dk][s]['Down'].GetBinContent(bin)
				thisfluct = (oldupcont-olddowncont)/nomhisto.GetBinContent(bin)
				#if the abs(fluctuation)>20. it's an outlier bin; reset the values in the shifted histograms
				if abs(thisfluct)>1. :
				#if abs(thisfluct-flucts[dk][s]['avg'])>20.*flucts[dk][s]['stddev'] :
					nomcont=nomhisto.GetBinContent(bin)
					newupcont=nomcont*(1.+0.05*random.random()); newdncont=nomcont*(1.-0.05*random.random())
					print '		bin %d has a fluctuation %.6f: setting Up/Down content to %.6f / %.6f from %.6f / %.6f '%(bin,thisfluct,newupcont,newdncont,oldupcont,olddowncont)
					histo_sets[dk][s]['Up'].SetBinContent(bin,newupcont)
					histo_sets[dk][s]['Down'].SetBinContent(bin,newdncont)
			print '	Done.'
		#else :
			#print '		this systematic is skipped.'
		#set the histogram directories
		histo_sets[dk][s]['Up'].SetDirectory(outfilep)
		histo_sets[dk][s]['Down'].SetDirectory(outfilep)
	print 'Done.'

#write the new file and close it
outfilep.Write()
outfilep.Close()

#close the input file
infilep.Close()