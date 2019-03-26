#Imports
from ROOT import *
from template import Template

ZEROED_BIN_CONTENT=0.000001

gROOT.SetBatch()

#name of file with templates to skim from
#ifn = 'templates_powheg_dynamic_binning_corrected_all.root'
#ifn = 'templates_powheg_dynamic_binning_vetoed_v7_alt_in_halves.root'
#ifn = 'templates_powheg_dynamic_binning_aggregated_v5_all.root'
ifn = 'templates_powheg_dynamic_binning_vetoed_v10_in_halves.root'
#name of file to put new templates in
#ofn = 'templates_powheg_dynamic_binning_vetoed_v10_in_halves.root'
ofn = 'templates_powheg_dynamic_binning_vetoed_v10.root'

#dict of channels with single/ranges of bins to remove from the fit
veto = {'t1_mu':{'singles':[],
				 'ranges':[]},
		't1_el':{'singles':[],
				 'ranges':[]},
		't2_mu':{'singles':[],
				 'ranges':[]},
		't2_el':{'singles':[],
				 'ranges':[]},
		't3_mu':{'singles':[],
				 'ranges':[]},
		't3_el':{'singles':[],
				 'ranges':[]},
}
#dicts of channels with bins to sum together
cstarpmbins = {'t1_mu':[(1,4),(4,7),
						(7,8),(8,9),(9,10),(10,11),
						(11,12),(12,13),
						],
			   't1_el':[(1,3),(3,5),
						(5,6),(6,7),
						(7,8),(8,9),
						],
			   't2_mu':[(1,4),(4,7),(7,10),(10,13),(13,16),(16,19),
						(19,22),(22,25),(25,28),(28,31),
						(31,34),(34,37),
						(37,38),(38,39),(39,40),(40,41),(41,42),(42,43),
						(43,44),(44,45),(45,46),(46,47),
						],
			   't2_el':[(1,3),(3,5),
						(5,6),(6,7),
						],
			   't3_mu':[(1,11),(11,21),(21,31),(31,41),(41,51),(51,61),
			   			(61,66),(66,71),(71,76),(76,81),(81,86),(86,91),(91,96),(96,101),(101,106),(106,111),
			   			(111,119),(119,127),(127,135),(135,143),
			   			(143,151),(151,159),(159,167),(167,175),
			   			(175,179),(179,183),(183,187),(187,191),
			   			(191,193),(193,195),(195,197),(197,199),
			   			(199,201),(201,203),(203,205),(205,207),
			   			(207,210),(210,213),
			   			(213,216),(216,219),
			   			(219,221),(221,223),
			   			],
			   't3_el':[(1,6),(6,11),(11,16),(16,21),(21,26),(26,31),(31,36),(36,41),(41,46),(46,51),
			   			(51,61),(61,71),(71,81),(81,91),
			   			(91,99),(99,107),(107,115),(115,123),
			   			(123,126),(126,129),(129,132),(132,135),(135,138),(138,141),(141,144),(144,147),
			   			(147,149),(149,151),(151,153),(153,155),
			   			(155,158),(158,161),
			   			(161,163),(163,165),
						],
}
tosuminhalves = {'t1_mu':[(1,4),  (7,8),(9,10),  (11,12),
						  (4,7),  (8,9),(10,11),  (12,13),
						  ],
				 't1_el':[(1,3),  (5,6),  (7,8),
				 		  (3,5),  (6,7),  (8,9),
						  ],
				 't2_mu':[(1,4),(7,10),(13,16),  (19,22),(25,28),  (31,34),  (37,38),(39,40),(41,42),  (43,44),(45,46),
				 		  (4,7),(10,13),(16,19),  (22,25),(28,31),  (34,37),  (38,39),(40,41),(42,43),  (44,45),(46,47),
						  ],
				 't2_el':[(1,3),  (5,6),
				 		  (3,5),  (6,7),
						  ],
				 't3_mu':[(1,11),(21,31),(41,51),  (61,66),(71,76),(81,86),(91,96),(101,106),  (111,119),(127,135),  (143,151),(159,167),
						  (175,179),(183,187),  (191,193),(195,197),  (199,201),(203,205),  (207,210),  (213,216),  (219,221),
				 		  (11,21),(31,41),(51,61),  (66,71),(76,81),(86,91),(96,101),(106,111),  (119,127),(135,143),  (151,159),(167,175),
				 		  (179,183),(187,191),  (193,195),(197,199),  (201,203),(205,207),  (210,213),  (216,219),  (221,223),
						  ],
				 't3_el':[(1,6),(11,16),(21,26),(31,36),(41,46),  (51,61),(71,81),  (91,99),(107,115),  (123,126),(129,132),(135,138),(141,144),
						  (147,149),(151,153),  (155,158),  (161,163),
				 		  (6,11),(16,21),(26,31),(36,41),(46,51),  (61,71),(81,91),  (99,107),(115,123),  (126,129),(132,135),(138,141),(144,147),
				 		  (149,151),(153,155),  (158,161),  (163,165),
						  ],
}
toremovemassbins = {'t1_mu':[(1,2),(2,4),(4,5),
							 (5,6),(6,8),(8,9),
							 ],
					't1_el':[], #already no mass binning
					't2_mu':[(1,4),(4,6),(6,7),(7,10),(10,12),
							 (12,15),(15,17),(17,18),(18,21),(21,23)],
					't2_el':[], #again already no mass binning
					't3_mu':[(1,4),(4,9),(9,11),(11,13),(13,15),(15,17),(17,19),(19,20),(20,21),(21,22),
							 (22,25),(25,30),(30,32),(32,34),(34,36),(36,38),(38,40),(40,41),(41,42),(42,43),],
					't3_el':[(1,6),(6,8),(8,10),(10,14),(14,16),(16,17),(17,18),
							 (18,23),(23,25),(25,27),(27,31),(31,33),(33,34),(34,35),],
}
cstarpmonly = {'t1_mu':[(1,5),(5,9)],
			   't1_el':[(1,4),(4,7)],
			   't2_mu':[(1,12),(12,23)],
			   't2_el':[(1,3),(3,5)],
			   't3_mu':[(1,22),(22,43)],
			   't3_el':[(1,18),(18,35)],
}


#set the one to use
#consolidate = tosuminhalves
consolidate = cstarpmonly

#check that you weren't very silly with the veto/consolidate options
goodopts = True
for c in veto.keys() :
	for cr in consolidate[c] :
		for s in veto[c]['singles'] :
			if s>=cr[0] and s<cr[1] :
				print('ERROR: requested bin {} to veto but also consolidate in channel {}'.format(s,c))
				goodopts=False
		for rt in veto[c]['ranges'] :
			for i in range(rt[0],rt[1]) :
				if i>=cr[0] and i<cr[1] :
					print('ERROR: requested bin {} to veto but also consolidate in channel {}'.format(s,c))
					goodopts=False
if not goodopts :
	print('Errors in binning options; can\'t run : (')
	exit()

#open the old file
ifp = TFile.Open(ifn,'r')

#open the new file
ofp = TFile.Open(ofn,'recreate')

#for every template in the file
all_new_histos = []
for k in ifp.GetListOfKeys() :
	n = k.GetName()
	#otherwise make the string to index by and get the single bins and bin ranges to veto
	csplits = k.GetName().split('__')[0].split('_')
	s = csplits[0]+'_'+csplits[1][:2]
	singles = veto[s]['singles']
	ranges = veto[s]['ranges']
	#make the list of all the bin numbers to veto and consolidate
	all_bins_to_veto = singles
	for r in ranges :
		for i in range(r[0],r[1]) :
			if not i in all_bins_to_veto :
				all_bins_to_veto.append(i)
	#get the old 1D template
	old_histo = k.ReadObj()
	if n.find('WJets_CR__')!=-1 or (len(all_bins_to_veto)==0 and len(consolidate[s])==0) :
		old_histo.SetDirectory(ofp)
		continue
	old_nbins = old_histo.GetNbinsX()
	#make the new 1D template
	new_nbins = old_nbins-len(all_bins_to_veto)
	for cr in consolidate[s] :
		new_nbins-=(cr[1]-cr[0])
	first_consolidated_bin = new_nbins+1
	new_nbins+=len(consolidate[s])
	new_histo = TH1D(n,old_histo.GetTitle(),new_nbins,0.,new_nbins-1.)
	all_new_histos.append(new_histo)
	#copy over bin contents skipping the bins to veto and summing over bins to consolidate
	new_bincounter = 1
	for i in range(1,old_nbins+1) :
		#only touch it if we don't want to veto it completely
		if not i in all_bins_to_veto :
			#first check if it's going in one of the consolidated bins
			consolidated=False
			for j in range(len(consolidate[s])) :
				if i>=consolidate[s][j][0] and i<consolidate[s][j][1] :
					tb = first_consolidated_bin+j
					#print('consolidating bin {} into bin {}'.format(i,tb))
					if old_histo.GetBinContent(i)!=ZEROED_BIN_CONTENT :
						new_histo.SetBinContent(tb,new_histo.GetBinContent(tb)+old_histo.GetBinContent(i))
					consolidated=True
			#if it wasn't added to one of the bins we're consolidating copy it over
			if not consolidated :
				new_histo.SetBinContent(new_bincounter,old_histo.GetBinContent(i))
				new_histo.SetBinError(new_bincounter,old_histo.GetBinError(i))
				new_bincounter+=1
	#set the new histogram's directory to the output file
	new_histo.SetDirectory(ofp)

#write the new file
print('writing new file')
ofp.Write()
ofp.Close()

ifp.Close()