#Imports
from ROOT import *

#filenames
old_file_name = 'templates_powheg_all.root'
new_file_name = 'templates_powheg_new_all.root'

#open old file
old_file = TFile.Open(old_file_name,'r')

#open new file to cd to its currently empty directory
new_file = TFile.Open(new_file_name,'recreate')

#make dictionary of all old histogram objects and put them in the new file
all_histos = {}
keylist = old_file.GetListOfKeys()
for k in keylist :
	#get the histogram name
	n = k.GetName()
	#get the histogram itself
	all_histos[n]=k.ReadObj()
	#check if it needs to be renamed and do so if that's the case
	JECstems = [('JES','Up'),('JES','Down'),('JER','Up'),('JER','Down')]
	for stemtuple in JECstems :
		ostem = stemtuple[0]+stemtuple[1]
		if n.find(ostem)!=-1 :
			print 'found name %s corresponding to stem %s'%(n,ostem)
			t = n.split('__')[0].split('_')[0]
			newname = n[:n.find(ostem)]+stemtuple[0]+'_'+t+stemtuple[1]
			print '	new name will be %s'%(newname)
			all_histos[n].SetName(newname)
	#set the histogram's disrectory to the new files
	all_histos[n].SetDirectory(new_file)
#write and close the new file
new_file.Write()
new_file.Close()
#close the old file
old_file.Close()