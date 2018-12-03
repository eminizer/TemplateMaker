#makes a "total templates" file from the JEC systematic-split files

#imports
from ROOT import TFile
import os

#filename for total output
outfilename = 'templates_powheg_NEW_JEC_all.root'

#channel names
channel_names = []
for t in ['t1','t2','t3'] :
	for leptype in ['muplus','muminus','elplus','elminus'] :
		for r in ['SR','WJets_CR'] :
			channel_names.append(t+'_'+leptype+'_'+r)
#process names for which JEC systematics matter
process_names = ['fqp0','fqm0','fq0','fq1','fq2','fg0','fg1','fg2','fg3','fg4','fbck','fwjets','fqcd']

#dict of filenames keyed by included JEC wiggles
pre = '../template_dirs/powheg_'
post = '/templates.root'
all_files = {pre+'JECnominal'+post: 	 'nominal',
			 pre+'AK4JESPU_Eta'+post: 	 ['AK4JESPU','AK4JESEta'],
			 pre+'AK4JESPt_Scale'+post:  ['AK4JESPt','AK4JESScale'],
			 pre+'AK4JESTime_Flav'+post: ['AK4JESTime','AK4JESFlav'],
			 pre+'AK4JER'+post: 		 ['AK4JERStat','AK4JERSys'],
			 pre+'AK8JESPU_Eta'+post: 	 ['AK8JESPU','AK8JESEta'],
			 pre+'AK8JESPt_Scale'+post:  ['AK8JESPt','AK8JESScale'],
			 pre+'AK8JESTime_Flav'+post: ['AK8JESTime','AK8JESFlav'],
			 pre+'AK8JER'+post: 		 ['AK8JERStat','AK8JERSys'],
			}
#dictionary of all histograms keyed by name only
all_template_histograms = {}

#open the output file
output_tfile = TFile.Open(outfilename,'recreate')
garbage_file = TFile.Open(outfilename+'_garbage','recreate')
garbage_file.cd()

#loop through files and add templates
fps = []
for fn,jecsys in all_files.iteritems() :
	#open the file and get its list of keys
	print 'doing file '+fn
	fp = TFile(fn,'r')
	fps.append(fp)
	lok = fp.GetListOfKeys()
	#if this the the nominal file, literally just get all of the histograms
	if jecsys=='nominal' :
		for k in lok :
			print 'saving non-JEC template '+str(k.GetName())
			all_template_histograms[k.GetName()] = k.ReadObj()
	#otherwise find only the necessary JEC-wiggled histograms
	else :
		#list of histogram names to look for
		histonames = []
		for cn in channel_names :
			for pn in process_names :
				for js in jecsys :
					histonames.append(cn+'__'+pn+'__'+js+'Up')
					histonames.append(cn+'__'+pn+'__'+js+'Down')
		for k in lok :
			n = k.GetName()
			if n in histonames :
				print 'saving JEC template '+n
				all_template_histograms[n]=k.ReadObj()

#save all the template histograms to the output file
for th in all_template_histograms.values() :
	th.SetDirectory(output_tfile)
output_tfile.Write()
output_tfile.Close()
for fp in fps :
	fp.Close()
garbage_file.Write()
garbage_file.Close()
os.system('rm -rf '+outfilename+'_garbage')

