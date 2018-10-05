from ROOT import *
import tdrstyle
from array import array
from template import Template

#TDR plot style stuff
gROOT.SetBatch()
tdrstyle.setTDRStyle()

#output file name
outfilen = 'binning_investigation.root'

#open input file
infilep = TFile('../total_template_files/templates_powheg_all.root')

#lepton types with fill styles
leptypes = ['elplus',     'elminus',     'muplus',     'muminus']
lepfills = {'elplus':3395,'elminus':3305,'muplus':3354,'muminus':3345}
gStyle.SetHatchesLineWidth(4)
gStyle.SetHatchesSpacing(0.5)
#processes and fill colors
procs = 	['fqp0',       'fqm0',          'fg0',        'fbck',          'fwjets']#,         'fqcd']
procfills = {'fqp0':kRed+2,'fqm0':kMagenta+2,'fg0':kBlue+2,'fbck':kYellow+2,'fwjets':kGreen+2}#,'fqcd':kCyan+2}
#analysis channels
templates = {'t1':{'SR':{},'WJets_CR':{}},'t2':{'SR':{},'WJets_CR':{}},'t3':{'SR':{}}}
for t in templates :
	for r in templates[t] :
		templates[t][r] = {}
		for lt in leptypes :
			templates[t][r][lt] = {}
			for p in procs :
				templates[t][r][lt][p] = {}
				n = t+'_'+lt+'_'+r+'__'+p
				templates[t][r][lt][p]['name']=n
				thishist = infilep.Get(n)
				templates[t][r][lt][p]['hist']=thishist
				newtemp = Template(n+'_temp',n+'_temp',None)
				newtemp.make_from_1D_histo(thishist)
				templates[t][r][lt][p]['temp']=newtemp

#hold on to all the plots
allplots = []
#declare the individual channel plots
ind_plots = {}
for t in templates :
	ind_plots[t]={}
	for r in templates[t] :
		ind_plots[t][r] = {}
		for lt in templates[t][r] :
			ind_plots[t][r][lt] = {}
			for p in templates[t][r][lt] :
				ind_plots[t][r][lt][p] = {}
				n = templates[t][r][lt][p]['name']+'_bad_bins_'
				this1Dhist = templates[t][r][lt][p]['temp'].convertTo1D()
				this3Dhist = templates[t][r][lt][p]['temp'].getHisto3D()
				ind_plots[t][r][lt][p]['bad_bins_x'] = TH1D(n+'x',n+'x',this3Dhist.GetNbinsX(),0.5,this3Dhist.GetNbinsX()+0.5)
				ind_plots[t][r][lt][p]['bad_bins_y'] = TH1D(n+'y',n+'y',this3Dhist.GetNbinsY(),0.5,this3Dhist.GetNbinsY()+0.5)
				ind_plots[t][r][lt][p]['bad_bins_z'] = TH1D(n+'z',n+'z',this3Dhist.GetNbinsZ(),0.5,this3Dhist.GetNbinsZ()+0.5)
				ind_plots[t][r][lt][p]['bad_bins_b'] = TH1D(n+'b',n+'b',this1Dhist.GetNbinsX(),0.5,this1Dhist.GetNbinsX()+0.5)
				ind_plots[t][r][lt][p]['bad_bins_x'].SetFillColor(procfills[p]); ind_plots[t][r][lt][p]['bad_bins_x'].SetFillStyle(lepfills[lt])
				ind_plots[t][r][lt][p]['bad_bins_y'].SetFillColor(procfills[p]); ind_plots[t][r][lt][p]['bad_bins_y'].SetFillStyle(lepfills[lt])
				ind_plots[t][r][lt][p]['bad_bins_z'].SetFillColor(procfills[p]); ind_plots[t][r][lt][p]['bad_bins_z'].SetFillStyle(lepfills[lt])
				ind_plots[t][r][lt][p]['bad_bins_b'].SetFillColor(procfills[p]); ind_plots[t][r][lt][p]['bad_bins_b'].SetFillStyle(lepfills[lt])
				allplots.append(ind_plots[t][r][lt][p]['bad_bins_x'])
				allplots.append(ind_plots[t][r][lt][p]['bad_bins_y'])
				allplots.append(ind_plots[t][r][lt][p]['bad_bins_z'])
				allplots.append(ind_plots[t][r][lt][p]['bad_bins_b'])
				#search template for zeroed bins
				content_zeroed = this1Dhist.GetMinimum()
				for gbin in range(1,this1Dhist.GetNbinsX()+1) :
					if this1Dhist.GetBinContent(gbin)==content_zeroed :
						ind_plots[t][r][lt][p]['bad_bins_b'].Fill(gbin)
				for gbin in range(this3Dhist.GetSize()) :
					if not this3Dhist.IsBinOverflow(gbin) and not this3Dhist.IsBinUnderflow(gbin) :
						#if this bin was zeroed
						if this3Dhist.GetBinContent(gbin)==content_zeroed :
							#increment the bad bin histograms in this projection bin
							binx = array('i',[0]); biny = array('i',[0]); binz = array('i',[0])
							this3Dhist.GetBinXYZ(gbin,binx,biny,binz)
							ind_plots[t][r][lt][p]['bad_bins_x'].Fill(float(binx[0]))
							ind_plots[t][r][lt][p]['bad_bins_y'].Fill(float(biny[0]))
							ind_plots[t][r][lt][p]['bad_bins_z'].Fill(float(binz[0]))

#stack the plots over lepton type for each process
ltype_summed_plots = {}
for t in ind_plots :
	ltype_summed_plots[t]={}
	for r in ind_plots[t] :
		ltype_summed_plots[t][r] = {}
		#get the first lepton type in the list
		ltypekeys = ind_plots[t][r].keys()
		#for each process
		for p in ind_plots[t][r][ltypekeys[0]] :
			n = t+'_'+r+'__'+p+'_bad_bins_'
			#create the new plots 
			this3Dhist = templates[t][r][ltypekeys[0]][p]['temp'].getHisto3D()
			ltype_summed_plots[t][r][p] = {}
			ltype_summed_plots[t][r][p]['bad_bins_x'] = THStack(n+'x',n+'x')
			ltype_summed_plots[t][r][p]['bad_bins_y'] = THStack(n+'y',n+'y')
			ltype_summed_plots[t][r][p]['bad_bins_z'] = THStack(n+'z',n+'z')
			ltype_summed_plots[t][r][p]['bad_bins_b'] = THStack(n+'b',n+'b')
			#add to them from the rest of the lepton types
			for lt in ltypekeys :
				ltype_summed_plots[t][r][p]['bad_bins_x'].Add(ind_plots[t][r][lt][p]['bad_bins_x'],'hist')
				ltype_summed_plots[t][r][p]['bad_bins_y'].Add(ind_plots[t][r][lt][p]['bad_bins_y'],'hist')
				ltype_summed_plots[t][r][p]['bad_bins_z'].Add(ind_plots[t][r][lt][p]['bad_bins_z'],'hist')
				ltype_summed_plots[t][r][p]['bad_bins_b'].Add(ind_plots[t][r][lt][p]['bad_bins_b'],'hist')

#stack the plots over process for each lepton type
proc_summed_plots = {}
for t in ind_plots :
	proc_summed_plots[t]={}
	for r in ind_plots[t] :
		proc_summed_plots[t][r] = {}
		#for each lepton type
		for lt in ind_plots[t][r] :
			#get the list of process keys
			pkeys = ind_plots[t][r][lt].keys()
			n = t+'_'+lt+'_'+r+'_bad_bins_'
			#create the new plots 
			this3Dhist = templates[t][r][lt][pkeys[0]]['temp'].getHisto3D()
			proc_summed_plots[t][r][lt] = {}
			proc_summed_plots[t][r][lt]['bad_bins_x'] = THStack(n+'x',n+'x')
			proc_summed_plots[t][r][lt]['bad_bins_y'] = THStack(n+'y',n+'y')
			proc_summed_plots[t][r][lt]['bad_bins_z'] = THStack(n+'z',n+'z')
			proc_summed_plots[t][r][lt]['bad_bins_b'] = THStack(n+'b',n+'b')
			#add to them from the rest of the lepton types
			for p in pkeys :
				proc_summed_plots[t][r][lt]['bad_bins_x'].Add(ind_plots[t][r][lt][p]['bad_bins_x'],'hist')
				proc_summed_plots[t][r][lt]['bad_bins_y'].Add(ind_plots[t][r][lt][p]['bad_bins_y'],'hist')
				proc_summed_plots[t][r][lt]['bad_bins_z'].Add(ind_plots[t][r][lt][p]['bad_bins_z'],'hist')
				proc_summed_plots[t][r][lt]['bad_bins_b'].Add(ind_plots[t][r][lt][p]['bad_bins_b'],'hist')

#stack the plots over lepton type AND process
ltype_proc_summed_plots = {}
for t in ind_plots :
	ltype_proc_summed_plots[t]={}
	for r in ind_plots[t] :
		#get the lists of lepton type and process keys
		ltkeys=ind_plots[t][r].keys()
		pkeys =ind_plots[t][r][ltkeys[0]].keys()
		n = t+'_'+r+'_bad_bins_'
		#create the new plots 
		this3Dhist = templates[t][r][ltkeys[0]][pkeys[0]]['temp'].getHisto3D()
		ltype_proc_summed_plots[t][r] = {}
		ltype_proc_summed_plots[t][r]['bad_bins_x'] = THStack(n+'x',n+'x')
		ltype_proc_summed_plots[t][r]['bad_bins_y'] = THStack(n+'y',n+'y')
		ltype_proc_summed_plots[t][r]['bad_bins_z'] = THStack(n+'z',n+'z')
		ltype_proc_summed_plots[t][r]['bad_bins_b'] = THStack(n+'b',n+'b')
		#add to them from the rest of the lepton and process types
		for lt in ltkeys :
			for p in pkeys :
				ltype_proc_summed_plots[t][r]['bad_bins_x'].Add(ind_plots[t][r][lt][p]['bad_bins_x'],'hist')
				ltype_proc_summed_plots[t][r]['bad_bins_y'].Add(ind_plots[t][r][lt][p]['bad_bins_y'],'hist')
				ltype_proc_summed_plots[t][r]['bad_bins_z'].Add(ind_plots[t][r][lt][p]['bad_bins_z'],'hist')
				ltype_proc_summed_plots[t][r]['bad_bins_b'].Add(ind_plots[t][r][lt][p]['bad_bins_b'],'hist')

#set the minima all at zero
for plot in allplots :
	plot.SetMinimum(0)

#open the output file
outfilep = TFile(outfilen,'recreate')
#make its directory structure
ltype_proc_summed_dir = outfilep.mkdir('ltype_proc_summed','ltype_proc_summed')
proc_summed_dir = outfilep.mkdir('proc_summed','proc_summed')
proc_summed_subdirs = {}
for ltype in leptypes :
	proc_summed_subdirs[ltype] = proc_summed_dir.mkdir(ltype,ltype)
ltype_summed_dir = outfilep.mkdir('ltype_summed','ltype_summed')
ltype_summed_subdirs = {}
for p in procs :
	ltype_summed_subdirs[p] = ltype_summed_dir.mkdir(p,p)
ind_dir = outfilep.mkdir('individual','individual')
ind_subdirs = {}
for t in templates :
	for r in templates[t] :
		for lt in templates[t][r] :
			ind_subdirs[t+'_'+lt+'_'+r] = ind_dir.mkdir(t+'_'+lt+'_'+r,t+'_'+lt+'_'+r)
#print 'ind_subdirs = %s'%(ind_subdirs)
outfilep.cd()

allcanvs = []
#set plots directories
for t in templates :
	for r in templates[t] :
		#lepton type/process-summed plots
		ltype_proc_summed_dir.cd()
		newcanv = TCanvas(t+'_'+r,t+'_'+r,1000,1000)
		newcanv.Divide(2,2)
		allcanvs.append(newcanv)
		i=1
		for plot in ltype_proc_summed_plots[t][r].values() :
			newcanv.cd(i); plot.Draw(); i+=1
		newcanv.Write()
		#lepton-type summed plots
		ltkeys = templates[t][r].keys()
		for p in templates[t][r][ltkeys[0]] :
			ltype_summed_subdirs[p].cd()
			newcanv = TCanvas(t+'_'+r+'_'+p,t+'_'+r+'_'+p,1000,1000)
			newcanv.Divide(2,2)
			allcanvs.append(newcanv)
			i=1
			for plot in ltype_summed_plots[t][r][p].values() :
				newcanv.cd(i); plot.Draw(); i+=1
			newcanv.Write()
		#down into lepton type
		for lt in templates[t][r] :
			#process-summed plots
			proc_summed_subdirs[lt].cd()
			newcanv = TCanvas(t+'_'+lt+'_'+r,t+'_'+lt+'_'+r,1000,1000)
			newcanv.Divide(2,2)
			allcanvs.append(newcanv)
			i=1
			for plot in proc_summed_plots[t][r][lt].values() :
				newcanv.cd(i); plot.Draw(); i+=1
			newcanv.Write()
			#down into process
			for p in templates[t][r][lt] :
				#individual plots
				for plot in ind_plots[t][r][lt][p].values() :
					plot.SetDirectory(ind_subdirs[t+'_'+lt+'_'+r])

#write and close output file
outfilep.Write()
outfilep.Close()

#close the input file
infilep.Close()


