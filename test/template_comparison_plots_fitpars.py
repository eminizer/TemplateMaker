from ROOT import *
import CMS_lumi, tdrstyle
from datetime import date
from optparse import OptionParser

# COMMAND LINE OPTIONS
parser = OptionParser()
parser.add_option('--template_file',  metavar='F', type='string', action='store', 
                              dest='templatefile',  help='which template auxiliary file to pull from?')
parser.add_option('--outtag', metavar='F', type='string', action='store', 
                              default='', 
                              dest='outtag') ## name for output file
parser.add_option('-M','--mode', type='choice', action='store', dest='mode', choices=['points','lines'], 
							  default='lines',
							  help='Plot as "lines" or as "points" with error bars?')
(options, args) = parser.parse_args()

#TDR plot style stuff
gROOT.SetBatch()
tdrstyle.setTDRStyle()
iPeriod = 4 #13TeV iPeriod = 1*(0/1 7 TeV) + 2*(0/1 8 TeV)  + 4*(0/1 13 TeV)
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Preliminary"

#function to get a histogram from any subdirectory of the file
def getHistFromFile(filep,histname) :
	print 'looking for histogram %s...'%(histname)
	for k in filep.GetListOfKeys() :
		if k.GetName()==histname :
			#print 'returning based on key %s'%(k)
			newhist = filep.Get(histname).Clone()
			newhist.SetTitle('')
			newhist.GetYaxis().SetTitle('Events')
			newhist.SetStats(0)
			return newhist
	print 'failed to find histogram %s'%(histname)
	return None


#open input file
infilep = TFile(options.templatefile)

#lepton types to sum over
leptypes = ['muplus','muminus','elplus','elminus']

#histogram name stems
noms = 		['__fqq', 			   '__fgg', 			'__fbck', '__fwjets', '__fqcd']
afb_ups = 	['__fqq__par_AfbUp',   None, 				None, 	  None, 	  '__fqcd__par_AfbUp']
afb_downs = ['__fqq__par_AfbDown', None, 				None, 	  None, 	  '__fqcd__par_AfbDown']
d_ups = 	['__fqq__par_dUp', 	   '__fgg__par_dUp', 	None, 	  None, 	  '__fqcd__par_dUp']
d_downs = 	['__fqq__par_dDown',   '__fgg__par_dDown', 	None, 	  None, 	  '__fqcd__par_dDown']
mu_ups = 	['__fqq__par_muUp',    '__fgg__par_muUp', 	None, 	  None, 	  '__fqcd__par_muUp']
mu_downs = 	['__fqq__par_muDown',  '__fgg__par_muDown', None, 	  None, 	  '__fqcd__par_muDown']
allnames = [noms,afb_ups,afb_downs,d_ups,d_downs,mu_ups,mu_downs]
colors_dists = [kRed,kBlue,kMagenta,kGreen,kYellow]
colors_types = [0,3,3,-3,-3,-9,-9]
leg_names_dists = ['q#bar{q} #rightarrow t#bar{t}','gg/qg etc. #rightarrow t#bar{t}','other top background','W+Jets background','QCD background']
leg_names_types = ['nominal','Afb up','Afb down','d up','d down','#mu up','#mu down']

#dictionary that will eventually hold all of the histos organized by topology and region
hists = {'t1':{'SR':[],'WJets_CR':[]},'t2':{'SR':[],'WJets_CR':[]},'t3':{'SR':[]}}

#open the output file
outname = 'template_comparison_plots_'+options.mode+'_'+str(date.today())
if options.outtag!='' :
	outname+='_'+options.outtag
outname+='.root'
outfilep = TFile(outname,'recreate')
#outfilep = TFile('template_comparison_plots_error_bars.root','recreate')

#make canvases, legends, and CMS lumi objects for each plot
canvs 	  = {'t1':{'SR':[],'WJets_CR':[]},'t2':{'SR':[],'WJets_CR':[]},'t3':{'SR':[]}}
legs 	  = {'t1':{'SR':[],'WJets_CR':[]},'t2':{'SR':[],'WJets_CR':[]},'t3':{'SR':[]}}
lumi_objs = {'t1':{'SR':[],'WJets_CR':[]},'t2':{'SR':[],'WJets_CR':[]},'t3':{'SR':[]}}
for top in hists :
	for reg in hists[top] :
		for nom in noms :
			if (top=='t1' or (top=='t2' and reg=='SR')) and nom=='__fqcd' :
				continue
			name = top+'_'+reg+'_'+nom.lstrip('__')+'_canv'
			canvs[top][reg].append(TCanvas(name,name,1100,900))
			legs[top][reg].append(TLegend(0.62,0.67,0.9,0.9))

#for each of the topologies
for top in hists :
	#and each of the regions
	for reg in hists[top] :
		#for each of the plot types
		for i in range(len(noms)) :
			if (top=='t1' or (top=='t2' and reg=='SR')) and noms[i]=='__fqcd' :
				continue
			thishistlist = hists[top][reg]
			thishistlist.append([])
			#get the nominal histogram in the first channel
			thishistlist[i].append(getHistFromFile(infilep,top+'_'+leptypes[0]+'_'+reg+noms[i]+'_x'))
			#sum over the other lepton types
			for k in range(1,len(leptypes)) :
				newname = top+'_'+leptypes[k]+'_'+reg+noms[i]+'_x'
				thishistlist[i][0].Add(getHistFromFile(infilep,newname).Clone())
			if reg!='WJets_CR' :
				#also get the other histograms in the first channel
				for j in range(1,len(allnames)) :
					if allnames[j][i]!=None :
						thishistlist[i].append(getHistFromFile(infilep,top+'_'+leptypes[0]+'_'+reg+allnames[j][i]+'_x').Clone())
						#and for each of those sum up from all of the channels
						for k in range(1,len(leptypes)) :
							newname = top+'_'+leptypes[k]+'_'+reg+allnames[j][i]+'_x'
							thishistlist[i][-1].Add(getHistFromFile(infilep,newname).Clone())

#Set plot attributes and add to legend
for top in hists :
	#and each of the regions
	for reg in hists[top] :
		#for each type of histogram
		for i in range(len(leg_names_dists)) :
			if (top=='t1' or (top=='t2' and reg=='SR')) and noms[i]=='__fqcd' :
				continue
			thishistlist = thishistlist = hists[top][reg]
			histj = 0
			for j in range(len(leg_names_types)) :
				if allnames[j][i]==None or (reg=='WJets_CR' and j!=0) :
					continue
				thishistlist[i][histj].SetLineWidth(4)
				thishistlist[i][histj].SetMarkerStyle(21)
				thishistlist[i][histj].SetFillStyle(0)
				if (j-1)%2==0 :
					thishistlist[i][histj].SetLineStyle(7)
				elif j>0 :
					thishistlist[i][histj].SetLineStyle(3)
				thishistlist[i][histj].SetLineColor(colors_dists[i]+colors_types[j])
				legs[top][reg][i].AddEntry(thishistlist[i][histj],leg_names_dists[i]+', '+leg_names_types[j],'L')
				histj+=1

#Plot plots
for top in hists :
	for reg in hists[top] :
		for i in range(len(hists[top][reg])) :
			thismax = hists[top][reg][i][0].GetMaximum()
			if reg!='WJets_CR' :
				for j in range(1,len(hists[top][reg][i])) :
					newmax = hists[top][reg][i][j].GetMaximum()
					if newmax>thismax :
						thismax = newmax
			hists[top][reg][i][0].GetYaxis().SetRangeUser(0.,1.1*thismax)
			canvs[top][reg][i].cd()
			if options.mode=='lines' :
				hists[top][reg][i][0].Draw('HIST')
			elif options.mode=='points' :
				hists[top][reg][i][0].Draw('E')
			if reg!='WJets_CR' :
				for j in range(1,len(hists[top][reg][i])) :
					if options.mode=='lines' :
						hists[top][reg][i][j].Draw('SAME HIST')
					elif options.mode=='points' :
						hists[top][reg][i][j].Draw('SAME E')
			legs[top][reg][i].Draw()
			lumi_objs[top][reg].append(CMS_lumi.CMS_lumi(canvs[top][reg][i], iPeriod, 0))

#Save plots
outfilep.cd()
for top in hists :
	for reg in hists[top] :
		for canv in canvs[top][reg] :
			canv.Write()