#Imports
from ROOT import *

#filenames
majority_fn = '../new_binning_v4/templates.root'
t3_fn = '../new_binning_v2/templates.root'
t1_el_CR_fn = '../new_binning_v1/templates.root'
new_fn = 'templates.root'

#open new file 
new_file = TFile.Open(new_fn,'recreate')

#all histos so you don't lose 'em by accident
all_histos = {}

#add most of the templates
maj_f = TFile.Open(majority_fn,'r')
keylist = maj_f.GetListOfKeys()
for k in  keylist :
	n = k.GetName()
	if n.startswith('t3_') or n.startswith('t1_elplus_WJets_CR__') or n.startswith('t1_elminus_WJets_CR__') :
		continue
	else :
		all_histos[n]=k.ReadObj()
		all_histos[n].SetDirectory(new_file)

#add the type 3 templates
t3_f = TFile.Open(t3_fn,'r')
keylist = t3_f.GetListOfKeys()
for k in  keylist :
	n = k.GetName()
	if n.startswith('t3_') :
		all_histos[n]=k.ReadObj()
		all_histos[n].SetDirectory(new_file)
	else :
		continue

#add the type-1 electrons CR templates
t1_el_CR_f = TFile.Open(t1_el_CR_fn,'r')
keylist = t1_el_CR_f.GetListOfKeys()
for k in  keylist :
	n = k.GetName()
	if n.startswith('t1_elplus_WJets_CR__') or n.startswith('t1_elminus_WJets_CR__') :
		all_histos[n]=k.ReadObj()
		all_histos[n].SetDirectory(new_file)
	else :
		continue

#write and close the new file
new_file.Write()
new_file.Close()
#close the old files
maj_f.Close()
t3_f.Close()
t1_el_CR_f.Close()