#imports
from ROOT import *
from math import sqrt

gROOT.SetBatch()

#name of finer-binned template file
old_fn='templates_powheg_dynamic_binning_aggregated_sep_sample_symmetrized_all.root'

#name of coarsely binned template file (used to calculate average factors)
avg_fac_fn='templates_powheg_dynamic_binning_consolidated_for_smoothing.root'

#the final template filename
new_fn='templates_powheg_dynamic_binning_smoothed_with_JEC.root'

#the comparison plots filename
plot_fn='avg_fac_comparison_plots_with_JEC.root'

#dict of bin ranges that will all have the same factor for the separate sample systematics (i.e. those that were consolidated into the coarsely-binned file)
##for coarse templates with two bins per channel (just c* >/< 0)
#bins = {
#  't1_mu':{1:[(1,4),  (7,8),(9,10),  (11,12)],
#           2:[(4,7),  (8,9),(10,11),  (12,13)]},
#  't1_el':{1:[(1,3),  (5,6),  (7,8)],
#           2:[(3,5),  (6,7),  (8,9)]},
#  't2_mu':{1:[(1,4),(7,10),(13,16),  (19,22),(25,28),  (31,34),  (37,38),(39,40),(41,42),  (43,44),(45,46)],
#           2:[(4,7),(10,13),(16,19),  (22,25),(28,31),  (34,37),  (38,39),(40,41),(42,43),  (44,45),(46,47)]},
#  't2_el':{1:[(1,3),  (5,6)],
#           2:[(3,5),  (6,7)]},
#  't3_mu':{1:[(1,11),(21,31),(41,51),  (61,66),(71,76),(81,86),(91,96),(101,106),  (111,119),(127,135),  (143,151),(159,167),
#              (175,179),(183,187),  (191,193),(195,197),  (199,201),(203,205),  (207,210),  (213,216),  (219,221)],
#           2:[(11,21),(31,41),(51,61),  (66,71),(76,81),(86,91),(96,101),(106,111),  (119,127),(135,143),  (151,159),(167,175),
#              (179,183),(187,191),  (193,195),(197,199),  (201,203),(205,207),  (210,213),  (216,219),  (221,223)]},
#  't3_el':{1:[(1,6),(11,16),(21,26),(31,36),(41,46),  (51,61),(71,81),  (91,99),(107,115),  (123,126),(129,132),(135,138),(141,144),
#              (147,149),(151,153),  (155,158),  (161,163)],
#           2:[(6,11),(16,21),(26,31),(36,41),(46,51),  (61,71),(81,91),  (99,107),(115,123),  (126,129),(132,135),(138,141),(144,147),
#              (149,151),(153,155),  (158,161),  (163,165)]}
#}
##for coarse templates that are just summed over all +/- cstar bins (i.e. intact x_F and M binning)
#bins = {
#  't1_mu':{1:[(1,4)],2:[(4,7)],3:[(7,8)],4:[(8,9)],5:[(9,10)],6:[(10,11)],7:[(11,12)],8:[(12,13)],},
#  't1_el':{1:[(1,3)],2:[(3,5)],3:[(5,6)],4:[(6,7)],5:[(7,8)],6:[(8,9)],},
#  't2_mu':{1:[(1,4)],2:[(4,7)],3:[(7,10)],4:[(10,13)],5:[(13,16)],6:[(16,19)],7:[(19,22)],8:[(22,25)],9:[(25,28)],10:[(28,31)],
#           11:[(31,34)],12:[(34,37)],13:[(37,38)],14:[(38,39)],15:[(39,40)],16:[(40,41)],17:[(41,42)],18:[(42,43)],19:[(43,44)],20:[(44,45)],
#           21:[(45,46)],22:[(46,47)],},
#  't2_el':{1:[(1,3)],2:[(3,5)],3:[(5,6)],4:[(6,7)],},
#  't3_mu':{1:[(1,11)],2:[(11,21)],3:[(21,31)],4:[(31,41)],5:[(41,51)],6:[(51,61)],7:[(61,66)],8:[(66,71)],9:[(71,76)],10:[(76,81)],
#           11:[(81,86)],12:[(86,91)],13:[(91,96)],14:[(96,101)],15:[(101,106)],16:[(106,111)],17:[(111,119)],18:[(119,127)],19:[(127,135)],20:[(135,143)],
#           21:[(143,151)],22:[(151,159)],23:[(159,167)],24:[(167,175)],25:[(175,179)],26:[(179,183)],27:[(183,187)],28:[(187,191)],29:[(191,193)],30:[(193,195)],
#           31:[(195,197)],32:[(197,199)],33:[(199,201)],34:[(201,203)],35:[(203,205)],36:[(205,207)],37:[(207,210)],38:[(210,213)],39:[(213,216)],40:[(216,219)],
#           41:[(219,221)],42:[(221,223)],},
#  't3_el':{1:[(1,6)],2:[(6,11)],3:[(11,16)],4:[(16,21)],5:[(21,26)],6:[(26,31)],7:[(31,36)],8:[(36,41)],9:[(41,46)],10:[(46,51)],
#           11:[(51,61)],12:[(61,71)],13:[(71,81)],14:[(81,91)],15:[(91,99)],16:[(99,107)],17:[(107,115)],18:[(115,123)],19:[(123,126)],20:[(126,129)],
#           21:[(129,132)],22:[(132,135)],23:[(135,138)],24:[(138,141)],25:[(141,144)],26:[(144,147)],27:[(147,149)],28:[(149,151)],29:[(151,153)],30:[(153,155)],
#           31:[(155,158)],32:[(158,161)],33:[(161,163)],34:[(163,165)],}
#}
#for coarse templates that are have no mass bins and only two bins in cstar (the bin order is kinda weird also)
bins = {
  't1_mu':{1:[(1,4),],2:[(7,8),(9,10),],3:[(11,12),],4:[(4,7),],5:[(8,9),(10,11),],6:[(12,13),],},
  't1_el':{1:[(1,3),],2:[(5,6),],3:[(7,8),],4:[(3,5),],5:[(6,7),],6:[(8,9),],},
  't2_mu':{1:[(1,4),(7,10),(13,16),],2:[(19,22),(25,28),],3:[(31,34),],4:[(37,38),(39,40),(41,42),],5:[(43,44),(45,46),],6:[(4,7),(10,13),(16,19),],7:[(22,25),(28,31),],8:[(34,37),],9:[(38,39),(40,41),(42,43),],10:[(44,45),(46,47),],},
  't2_el':{1:[(1,3)],2:[(5,6)],3:[(3,5)],4:[(6,7)],},
  't3_mu':{1:[(1,11),(21,31),(41,51)],2:[(61,66),(71,76),(81,86),(91,96),(101,106)],3:[(111,119),(127,135)],4:[(143,151),(159,167)],5:[(175,179),(183,187)],
           6:[(191,193),(195,197)],7:[(199,201),(203,205)],8:[(207,210)],9:[(213,216)],10:[(219,221)],
           11:[(11,21),(31,41),(51,61)],12:[(66,71),(76,81),(86,91),(96,101),(106,111)],13:[(119,127),(135,143)],14:[(151,159),(167,175)],15:[(179,183),(187,191)],
           16:[(193,195),(197,199)],17:[(201,203),(205,207)],18:[(210,213)],19:[(216,219)],20:[(221,223),],},
  't3_el':{1:[(1,6),(11,16),(21,26),(31,36),(41,46)],2:[(51,61),(71,81)],3:[(91,99),(107,115)],4:[(123,126),(129,132),(135,138),(141,144)],5:[(147,149),(151,153)],
           6:[(155,158)],7:[(161,163)],8:[(6,11),(16,21),(26,31),(36,41),(46,51)],9:[(61,71),(81,91)],10:[(99,107),(115,123)],
           11:[(126,129),(132,135),(138,141),(144,147)],12:[(149,151),(153,155)],13:[(158,161)],14:[(163,165)],}
}

#which systematics should we check?
sys_to_correct = ['JES','JER','isr','fsr','hdamp','tune','cr']

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

#dictionary of the colors per channel
ccolors = {'t1_muplus_SR':kOrange+9,
       't1_muminus_SR':kOrange+9,
       't1_elplus_SR':kViolet-4,
       't1_elminus_SR':kViolet-4,
       't1_muplus_WJets_CR':kOrange-3,
       't1_muminus_WJets_CR':kOrange-3,
       't1_elplus_WJets_CR':kViolet-6,
       't1_elminus_WJets_CR':kViolet-6,
       't2_muplus_SR':kGreen+3,
       't2_muminus_SR':kGreen+3,
       't2_elplus_SR':kMagenta+3,
       't2_elminus_SR':kMagenta+3,
       't2_muplus_WJets_CR':kGreen-7,
       't2_muminus_WJets_CR':kGreen-7,
       't2_elplus_WJets_CR':kMagenta-7,
       't2_elminus_WJets_CR':kMagenta-7,
       't3_muplus_SR':kRed+2,
       't3_muminus_SR':kRed+2,
       't3_elplus_SR':kBlue+2,
       't3_elminus_SR':kBlue+2,
}

#open the old files
fine_f = TFile.Open(old_fn,'r')
coarse_f = TFile.Open(avg_fac_fn,'r')

#open the new files
new_f = TFile.Open(new_fn,'recreate')
plot_f = TFile.Open(plot_fn,'recreate')

#for each of the systematics
altered_keys = []
all_canvs = []
all_histos = []
for sys in sys_to_correct :
  pnames = process_names+['fwjets'] if sys in ['JES','JER'] else process_names
  totaldevcanv = TCanvas(sys+'_all_devplots',sys+'_all_devplots',1200,900); all_canvs.append(totaldevcanv)
  updev_stack = THStack(sys+'_updev_stack',sys+' up smoothing changes; raw/smoothed #Delta/#sigma; number of bins')
  dndev_stack = THStack(sys+'_dndev_stack',sys+' down smoothing changes; raw/smoothed #Delta/#sigma; number of bins')
  all_histos.append(updev_stack); all_histos.append(dndev_stack)
  #for each channel and process
  for c in cnames :
    #make the shift deviation plot objects
    plot_f.cd()
    devcanv = TCanvas(sys+'_'+c+'_devplots',sys+'_'+c+'_devplots',1200,900); all_canvs.append(devcanv)
    updev_h = TH1D(sys+'_'+c+'_updevs',sys+' up '+c+' smoothing changes; raw/smoothed #Delta/#sigma; number of bins',100,-5.,5.)
    dndev_h = TH1D(sys+'_'+c+'_dndevs',sys+' down '+c+' smoothing changes; raw/smoothed #Delta/#sigma; number of bins',100,-5.,5.)
    updev_h.SetMarkerStyle(21); dndev_h.SetMarkerStyle(21)
    updev_h.SetMarkerColor(ccolors[c]); dndev_h.SetMarkerColor(ccolors[c])
    updev_h.SetLineColor(ccolors[c]); dndev_h.SetLineColor(ccolors[c])
    updev_h.SetFillColor(ccolors[c]); dndev_h.SetFillColor(ccolors[c])
    all_histos.append(updev_h); all_histos.append(dndev_h)
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
      #make the shift plot objects
      plot_f.cd()
      shiftcanv = TCanvas(sys+'_'+c+'_'+p+'_shiftplots',sys+'_'+c+'_'+p+'_shiftplots',1500,900); all_canvs.append(shiftcanv)
      fine_upshifts_h = TH1D(sys+'_'+c+'_'+p+'_fine_upshifts',sys+'_'+c+'_'+p+'_fine_upshifts',nom.GetNbinsX(),0,nom.GetXaxis().GetBinUpEdge(nom.GetNbinsX()))
      fine_dnshifts_h = TH1D(sys+'_'+c+'_'+p+'_fine_dnshifts',sys+'_'+c+'_'+p+'_fine_dnshifts',nom.GetNbinsX(),0,nom.GetXaxis().GetBinUpEdge(nom.GetNbinsX()))
      coarse_upshifts_h = TH1D(sys+'_'+c+'_'+p+'_coarse_upshifts',sys+' '+c+' '+p+' shift comparison; template bin; shifted/nominal',nom.GetNbinsX(),0,nom.GetXaxis().GetBinUpEdge(nom.GetNbinsX()))
      coarse_dnshifts_h = TH1D(sys+'_'+c+'_'+p+'_coarse_dnshifts',sys+'_'+c+'_'+p+'_coarse_dnshifts',nom.GetNbinsX(),0,nom.GetXaxis().GetBinUpEdge(nom.GetNbinsX()))
      fine_upshifts_h.SetMarkerStyle(22); fine_upshifts_h.SetMarkerColor(kRed+2); fine_upshifts_h.SetLineColor(kRed+2)
      fine_dnshifts_h.SetMarkerStyle(23); fine_dnshifts_h.SetMarkerColor(kBlue+2); fine_dnshifts_h.SetLineColor(kBlue+2)
      coarse_upshifts_h.SetLineColor(kRed+2); coarse_upshifts_h.SetLineWidth(2)
      coarse_dnshifts_h.SetLineColor(kBlue+2); coarse_dnshifts_h.SetLineWidth(2)
      all_histos.append(fine_upshifts_h)
      all_histos.append(fine_dnshifts_h)
      all_histos.append(coarse_upshifts_h)
      all_histos.append(coarse_dnshifts_h)
      #for each coarse bin in the dictionary
      for coarsebin,finebinranges in thisbindict.items() :
        #find the up/down SFs
        upsf = up_coarse.GetBinContent(coarsebin)/nom_coarse.GetBinContent(coarsebin)
        dnsf = dn_coarse.GetBinContent(coarsebin)/nom_coarse.GetBinContent(coarsebin)
        print('channel: {}, proc: {}, sys: {}, coarsebin: {}, upsf: {}, dnsf: {}'.format(c,p,sys,coarsebin,upsf,dnsf)) #DEBUG
        #loop over the fine bin ranges
        for br in finebinranges :
          #first add to the plot histograms based on the old, raw templates
          for i in range(br[0],br[1]) :
            finenomcont = nom.GetBinContent(i); finenomerr = nom.GetBinError(i)
            fineupcont = up.GetBinContent(i); fineuperr = up.GetBinError(i)
            finedncont = dn.GetBinContent(i); finednerr = dn.GetBinError(i)
            thisfineupshift = fineupcont/finenomcont
            thisfineupshifterr = thisfineupshift*sqrt((fineuperr/fineupcont)**2+(finenomerr/finenomcont)**2)
            thisfinednshift = finedncont/finenomcont
            thisfinednshifterr = thisfinednshift*sqrt((finednerr/finedncont)**2+(finenomerr/finenomcont)**2)
            fine_upshifts_h.SetBinContent(i,thisfineupshift); fine_upshifts_h.SetBinError(i,thisfineupshifterr)
            fine_dnshifts_h.SetBinContent(i,thisfinednshift); fine_dnshifts_h.SetBinError(i,thisfinednshifterr)
            coarse_upshifts_h.SetBinContent(i,upsf); coarse_dnshifts_h.SetBinContent(i,dnsf)
            updev_h.Fill((thisfineupshift-upsf)/thisfineupshifterr)
            dndev_h.Fill((thisfinednshift-dnsf)/thisfinednshifterr)
          #then alter the template bin contents
          for i in range(br[0],br[1]) :
            nomcont=nom.GetBinContent(i)
            up.SetBinContent(i,upsf*nomcont)
            dn.SetBinContent(i,dnsf*nomcont)
      #make and write the plots
      shiftcanv.cd()
      coarse_upshifts_h.SetMinimum(min([coarse_upshifts_h.GetMinimum(),coarse_dnshifts_h.GetMinimum(),fine_upshifts_h.GetMinimum(),fine_dnshifts_h.GetMinimum()]))
      coarse_upshifts_h.SetMaximum(max([coarse_upshifts_h.GetMaximum(),coarse_dnshifts_h.GetMaximum(),fine_upshifts_h.GetMaximum(),fine_dnshifts_h.GetMaximum()]))
      coarse_upshifts_h.Draw("HIST")
      coarse_dnshifts_h.Draw("HIST SAME")
      fine_upshifts_h.Draw("PE1X0 SAME")
      fine_dnshifts_h.Draw("PE1X0 SAME")
      shiftcanv.Write()
      #set the templates to write to the new file
      up.SetDirectory(new_f)
      dn.SetDirectory(new_f)
    #plot the delta/sigma deviations due to the smoothing in  this channel overall
    devcanv.cd()
    devcanv.Divide(1,2)
    devcanv.cd(1)
    updev_h.Draw()
    devcanv.cd(2)
    dndev_h.Draw()
    devcanv.Write()
    #add to the stakcs over all channels
    updev_stack.Add(updev_h)
    dndev_stack.Add(dndev_h)
  #plot the stacked deviation plots
  totaldevcanv.cd()
  totaldevcanv.Divide(1,2)
  totaldevcanv.cd(1)
  updev_stack.Draw()
  totaldevcanv.cd(2)
  dndev_stack.Draw()
  totaldevcanv.Write()

#copy over all the other stuff
for k in fine_f.GetListOfKeys() :
  n=k.GetName()
  if n in altered_keys :
    continue
  else :
    obj=k.ReadObj()
    obj.SetDirectory(new_f)

#write the new files
new_f.Write()
new_f.Close()
plot_f.Write()
plot_f.Close()

#close the old files
fine_f.Close()
coarse_f.Close()

