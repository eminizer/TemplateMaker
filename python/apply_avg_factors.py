#imports
from ROOT import TFile

#name of finer-binned template file
old_fn='templates_powheg_dynamic_binning_aggregated_v5_all.root'

#name of coarsely binned template file (used to calculate average factors)
avg_fac_fn='templates_powheg_dynamic_binning_vetoed_v10.root'

#the final template filename
new_fn='templates_powheg_dynamic_binning_avg_facs_v1_all.root'

#dict of bin ranges that will all have the same factor for the separate sample systematics (i.e. those that were consolidated into the coarsely-binned file)
bins = {
  't1_mu':{1:[(1,4),  (7,8),(9,10),  (11,12)],
           2:[(4,7),  (8,9),(10,11),  (12,13)]},
  't1_el':{1:[(1,3),  (5,6),  (7,8)],
           2:[(3,5),  (6,7),  (8,9)]},
  't2_mu':{1:[(1,4),(7,10),(13,16),  (19,22),(25,28),  (31,34),  (37,38),(39,40),(41,42),  (43,44),(45,46)],
           2:[(4,7),(10,13),(16,19),  (22,25),(28,31),  (34,37),  (38,39),(40,41),(42,43),  (44,45),(46,47)]},
  't2_el':{1:[(1,3),  (5,6)],
           2:[(3,5),  (6,7)]},
  't3_mu':{1:[(1,11),(21,31),(41,51),  (61,66),(71,76),(81,86),(91,96),(101,106),  (111,119),(127,135),  (143,151),(159,167),
              (175,179),(183,187),  (191,193),(195,197),  (199,201),(203,205),  (207,210),  (213,216),  (219,221)],
           2:[(11,21),(31,41),(51,61),  (66,71),(76,81),(86,91),(96,101),(106,111),  (119,127),(135,143),  (151,159),(167,175),
              (179,183),(187,191),  (193,195),(197,199),  (201,203),(205,207),  (210,213),  (216,219),  (221,223)]},
  't3_el':{1:[(1,6),(11,16),(21,26),(31,36),(41,46),  (51,61),(71,81),  (91,99),(107,115),  (123,126),(129,132),(135,138),(141,144),
              (147,149),(151,153),  (155,158),  (161,163)],
           2:[(6,11),(16,21),(26,31),(36,41),(46,51),  (61,71),(81,91),  (99,107),(115,123),  (126,129),(132,135),(138,141),(144,147),
              (149,151),(153,155),  (158,161),  (163,165)]}
}

#which systematics should we check?
sys_to_correct = ['isr','fsr','hdamp','tune','cr']

#channel names
cnames = []
for t in ['t1','t2','t3'] :
  for leptype in ['muplus','muminus','elplus','elminus'] :
    for r in ['SR'] : #,'WJets_CR'] :
      if t=='t3' and r=='WJets_CR' :
        continue
      cnames.append(t+'_'+leptype+'_'+r)

#process names
process_names = ['fqq','fgg','fqp0','fqm0','fq0','fq1','fq2','fg0','fg1','fg2','fg3','fg4','fbck']

#open the old files
fine_f = TFile.Open(old_fn,'r')
coarse_f = TFile.Open(avg_fac_fn,'r')

#open the new file
new_f = TFile.Open(new_fn,'recreate')

#for each of the systematics
altered_keys = []
for sys in sys_to_correct :
  pnames = process_names+['fwjets'] if sys in ['JES','JER'] else process_names
  #for each channel and process
  for c in cnames :
    for p in  pnames :
      #get the old, finely-binned templates
      nom=fine_f.Get(c+'__'+p)
      up=fine_f.Get(c+'__'+p+'__'+sys+'Up'); altered_keys.append(c+'__'+p+'__'+sys+'Up')
      dn=fine_f.Get(c+'__'+p+'__'+sys+'Down'); altered_keys.append(c+'__'+p+'__'+sys+'Down')
      #get the coarsely-binned templates
      nom_coarse=coarse_f.Get(c+'__'+p)
      up_coarse=coarse_f.Get(c+'__'+p+'__'+sys+'Up')
      dn_coarse=coarse_f.Get(c+'__'+p+'__'+sys+'Down')
      #find which bin dict to use
      thisbindict = bins[c[:5]]
      #for each coarse bin in the dictionary
      for coarsebin,finebinranges in thisbindict.items() :
        #find the up/down SFs
        upsf = up_coarse.GetBinContent(coarsebin)/nom_coarse.GetBinContent(coarsebin)
        dnsf = dn_coarse.GetBinContent(coarsebin)/nom_coarse.GetBinContent(coarsebin)
        print('channel: {}, proc: {}, sys: {}, coarsebin: {}, upsf: {}, dnsf: {}'.format(c,p,sys,coarsebin,upsf,dnsf)) #DEBUG
        #loop over the fine bin ranges
        for br in finebinranges :
          for i in range(br[0],br[1]) :
            nomcont=nom.GetBinContent(i)
            up.SetBinContent(i,upsf*nomcont)
            dn.SetBinContent(i,dnsf*nomcont)
      #set them to write to the new file
      up.SetDirectory(new_f)
      dn.SetDirectory(new_f)

#copy over all the other stuff
for k in fine_f.GetListOfKeys() :
  n=k.GetName()
  if n in altered_keys :
    continue
  else :
    obj=k.ReadObj()
    obj.SetDirectory(new_f)

#write the new file
new_f.Write()
new_f.Close()

fine_f.Close()
coarse_f.Close()

