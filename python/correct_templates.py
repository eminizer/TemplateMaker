#imports
import random
from ROOT import TFile

ZEROED_BIN_CONTENT=0.000001

#correction function
def correct_template(nom,up,dn) :
  #loop over the bins
  for i in range(1,nom.GetNbinsX()+1) :
    #if this bin is zeroed in the nominal template, zero it in the shifted templates
    if nom.GetBinContent(i)==ZEROED_BIN_CONTENT :
      up.SetBinContent(i,ZEROED_BIN_CONTENT)
      dn.SetBinContent(i,ZEROED_BIN_CONTENT)
    #otherwise if it's zeroed in a shifted template, set that bin's content to the nominal content
    elif abs(up.GetBinContent(i)-ZEROED_BIN_CONTENT)<0.000005 :
      up.SetBinContent(i,nom.GetBinContent(i))
    elif abs(dn.GetBinContent(i)-ZEROED_BIN_CONTENT)<0.000005 :
      dn.SetBinContent(i,nom.GetBinContent(i))
  #return the shifted templates
  return up,dn

#the aggregated template filename
old_fn='templates_powheg_dynamic_binning_aggregated_v4_all.root'

#the corrected template filename
new_fn='templates_powheg_dynamic_binning_corrected_v4_all.root'

#which systematics should we check?
sys_to_correct = ['JES','JER','isr','fsr','hdamp','tune','cr']

#channel names
cnames = []
for t in ['t1','t2','t3'] :
  for leptype in ['muplus','muminus','elplus','elminus'] :
    for r in ['SR','WJets_CR'] :
      if t=='t3' and r=='WJets_CR' :
        continue
      cnames.append(t+'_'+leptype+'_'+r)

#process names
process_names = ['fqq','fgg','fqp0','fqm0','fq0','fq1','fq2','fg0','fg1','fg2','fg3','fg4','fbck']

#open the old file
old_f = TFile.Open(old_fn,'r')

#open the new file
new_f = TFile.Open(new_fn,'recreate')

#for each of the systematics
altered_keys = []
for sys in sys_to_correct :
  pnames = process_names+['fwjets'] if sys in ['JES','JER'] else process_names
  #for each channel and process
  for c in cnames :
    for p in  pnames :
      #get the old templates
      nom=old_f.Get(c+'__'+p); #altered_keys.append(c+'__'+p)
      up=old_f.Get(c+'__'+p+'__'+sys+'Up'); altered_keys.append(c+'__'+p+'__'+sys+'Up')
      dn=old_f.Get(c+'__'+p+'__'+sys+'Down'); altered_keys.append(c+'__'+p+'__'+sys+'Down')
      #get the corrected templates
      newup, newdown = correct_template(nom,up,dn)
      #set them to write to the new file
      #nom.SetDirectory(new_f)
      newup.SetDirectory(new_f)
      newdown.SetDirectory(new_f)

#copy over all the other stuff
for k in old_f.GetListOfKeys() :
  n=k.GetName()
  if n in altered_keys :
    continue
  else :
    obj=k.ReadObj()
    obj.SetDirectory(new_f)

#write the new file
new_f.Write()
new_f.Close()

old_f.Close()

