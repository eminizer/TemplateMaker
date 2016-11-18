from ROOT import *
from math import *
from optparse import OptionParser
from template import Template

parser = OptionParser()
parser.add_option('--tfilename',   type='string', action='store', dest='tfilename')
parser.add_option('--pfhfilename', type='string', action='store', dest='pfhfilename')
parser.add_option('--outfilename', type='string', action='store', dest='outfilename')
(options, args) = parser.parse_args()

#Make a list of all the channels in the file
pfhfile = TFile(options.pfhfilename)
keylist = pfhfile.GetListOfKeys()
channel_names = ['allchannels']
disttypes = []
disttypecolors = []
for i in range(len(keylist)) :
	hname = keylist.At(i).GetName()
	chan = hname.split('__')[0]
	if not chan in channel_names :
		channel_names.append(chan)
	disttype = hname.split('__')[1]
	if not disttype in disttypes :
		disttypes.append(disttype)
		if disttype == 'fqq' :
			disttypecolors.append(kRed+2)
		elif disttype == 'fgg' :
			disttypecolors.append(kBlue)
		elif disttype == 'fbck' :
			disttypecolors.append(kYellow)
		elif disttype == 'fntmj' :
			disttypecolors.append(kGreen)
#output file
outfile = TFile(options.outfilename,'recreate')
#Make lists of histogram stacks, MC uncertainty graphs, residual plots, canvases, and pads
x_stacks = []; y_stacks = []; z_stacks = []
x_err_hs = []; y_err_hs = []; z_err_hs = []
x_resids = []; y_resids = []; z_resids = []
x_canvs  = []; y_canvs  = []; z_canvs  = []
x_histo_pads = []; y_histo_pads = []; z_histo_pads = []
x_resid_pads = []; y_resid_pads = []; z_resid_pads = []
#make a new template to clone
dummy_template = Template('dummy','dummy',None)
for channame in channel_names :
	#Projection histo stacks
	if channame == 'allchannels' :
		x_stacks.append(THStack(channame+'_x_stack','all channels comparison plot, c* projection;;Events/0.2'))
		y_stacks.append(THStack(channame+'_y_stack','all channels comparison plot, x_{F} projection;;Events/bin'))
		z_stacks.append(THStack(channame+'_z_stack','all channels comparison plot, M projection;;Events/bin'))
	else :
		x_stacks.append(THStack(channame+'_x_stack',channame+' channel comparison plot, c* projection;;Events/0.2'))
		y_stacks.append(THStack(channame+'_y_stack',channame+' channel comparison plot, x_{F} projection;;Events/bin'))
		z_stacks.append(THStack(channame+'_z_stack',channame+' channel comparison plot, M projection;;Events/bin'))
	#MC Uncertainty graphs
	x_err_hs.append(dummy_template.getHistoX().Clone()); x_err_hs[len(x_err_hs)-1].SetFillColor(kBlack)
	y_err_hs.append(dummy_template.getHistoY().Clone()); y_err_hs[len(y_err_hs)-1].SetFillColor(kBlack)
	z_err_hs.append(dummy_template.getHistoZ().Clone()); z_err_hs[len(z_err_hs)-1].SetFillColor(kBlack)
	#Residual plots
	x_resids.append(dummy_template.getHistoX().Clone())
	x_resids[len(x_resids)-1].SetTitle(';c*;(Data-MC)/#sigma'); x_resids[len(x_resids)-1].SetFillColor(kMagenta+3)
	y_resids.append(dummy_template.getHistoY().Clone())
	y_resids[len(x_resids)-1].SetTitle(';x_{F};(Data-MC)/#sigma'); y_resids[len(x_resids)-1].SetFillColor(kMagenta+3)
	z_resids.append(dummy_template.getHistoZ().Clone())
	z_resids[len(x_resids)-1].SetTitle(';M (GeV);(Data-MC)/#sigma'); z_resids[len(x_resids)-1].SetFillColor(kMagenta+3)
	#Canvases
	x_canvs.append(TCanvas(channame+'_x_canvas',channame+'_x_canvas',900,900))
	y_canvs.append(TCanvas(channame+'_y_canvas',channame+'_y_canvas',900,900))
	z_canvs.append(TCanvas(channame+'_z_canvas',channame+'_z_canvas',900,900))
#Set plot directories
for i in range(len(channel_names)) :
	x_err_hs[i].SetDirectory(0); y_err_hs[i].SetDirectory(0); z_err_hs[i].SetDirectory(0)
	x_resids[i].SetDirectory(0); y_resids[i].SetDirectory(0); z_resids[i].SetDirectory(0)
nxbins = x_resids[0].GetNbinsX(); nybins = y_resids[0].GetNbinsX(); nzbins = z_resids[0].GetNbinsX()
#build 3D histograms out of 1D post-fit histograms from theta; add to histogram stacks and MC uncertainty graphs
disthistos = []
for k in range(len(disttypes)) :
	for i in range(len(channel_names)) :
		channame = channel_names[i]
		if channame == 'allchannels' :
			continue
		#Get histogram from post-fit file and make a new template
		newname = channame+'__'+disttypes[k]
		newtemp = Template(newname+'__POSTFIT',newname+'__POSTFIT',None)
		new1Dhisto = pfhfile.Get(newname)
		if i==len(channel_names)-1 :
			disthistos.append(new1Dhisto)
		newtemp.make_from_1D_histo(new1Dhisto)
		#Set attributes and add to histogram stack
		x_histo = newtemp.getHistoX()
		x_histo.SetFillColor(disttypecolors[k]); x_histo.SetLineColor(disttypecolors[k]); x_histo.SetMarkerStyle(21); x_histo.SetMarkerColor(disttypecolors[k])
		x_stacks[i].Add(x_histo,'hist')
		x_stacks[0].Add(x_histo,'hist')
		#increment error values
		for j in range(x_histo.GetSize()) :
			if not x_histo.IsBinOverflow(j) and not x_histo.IsBinUnderflow(j) :
				x_err_hs[i].SetBinContent(j,x_err_hs[i].GetBinContent(j)+x_histo.GetBinContent(j))
				x_err_hs[i].SetBinError(j,x_err_hs[i].GetBinError(j)+x_histo.GetBinError(j)**2)
				x_err_hs[0].SetBinContent(j,x_err_hs[0].GetBinContent(j)+x_histo.GetBinContent(j))
				x_err_hs[0].SetBinError(j,x_err_hs[0].GetBinError(j)+x_histo.GetBinError(j)**2)
		#Repeat all of the above for y
		y_histo = newtemp.getHistoY()
		y_histo.SetFillColor(disttypecolors[k]); y_histo.SetLineColor(disttypecolors[k]); y_histo.SetMarkerStyle(21); y_histo.SetMarkerColor(disttypecolors[k])
		y_stacks[i].Add(y_histo,'hist')
		y_stacks[0].Add(y_histo,'hist')
		#increment error values
		for j in range(y_histo.GetSize()) :
			if not y_histo.IsBinOverflow(j) and not y_histo.IsBinUnderflow(j) :
				y_err_hs[i].SetBinContent(j,y_err_hs[i].GetBinContent(j)+y_histo.GetBinContent(j))
				y_err_hs[i].SetBinError(j,y_err_hs[i].GetBinError(j)+y_histo.GetBinError(j)**2)
				y_err_hs[0].SetBinContent(j,y_err_hs[0].GetBinContent(j)+y_histo.GetBinContent(j))
				y_err_hs[0].SetBinError(j,y_err_hs[0].GetBinError(j)+y_histo.GetBinError(j)**2)
		#Repeat all of the above for z
		z_histo = newtemp.getHistoZ()
		z_histo.SetFillColor(disttypecolors[k]); z_histo.SetLineColor(disttypecolors[k]); z_histo.SetMarkerStyle(21); z_histo.SetMarkerColor(disttypecolors[k])
		z_stacks[i].Add(z_histo,'hist')
		z_stacks[0].Add(z_histo,'hist')
		#increment error values
		for j in range(z_histo.GetSize()) :
			if not z_histo.IsBinOverflow(j) and not z_histo.IsBinUnderflow(j) :
				z_err_hs[i].SetBinContent(j,z_err_hs[i].GetBinContent(j)+z_histo.GetBinContent(j))
				z_err_hs[i].SetBinError(j,z_err_hs[i].GetBinError(j)+z_histo.GetBinError(j)**2)
				z_err_hs[0].SetBinContent(j,z_err_hs[0].GetBinContent(j)+z_histo.GetBinContent(j))
				z_err_hs[0].SetBinError(j,z_err_hs[0].GetBinError(j)+z_histo.GetBinError(j)**2)
#Take the root of all of the final MC uncertainties
for i in range(len(channel_names)) :
	for j in range(nxbins) :
		x_err_hs[i].SetBinError(j,sqrt(x_err_hs[i].GetBinError(j)))
	for j in range(nybins) :
		y_err_hs[i].SetBinError(j,sqrt(y_err_hs[i].GetBinError(j)))
	for j in range(nzbins) :
		z_err_hs[i].SetBinError(j,sqrt(z_err_hs[i].GetBinError(j)))
#build residuals plots
maxxdeviations = []; maxydeviations = []; maxzdeviations = []
#get the very first template file
newdatatemplates = []
initial_templates_file = TFile(options.tfilename)
for i in range(len(channel_names)) :
	channame = channel_names[i]
	maxxdeviations.append(0.0); maxydeviations.append(0.0); maxzdeviations.append(0.0)
	newname = channame+'__DATA'
	newdatatemp = Template(newname+'__POSTFIT',newname+'__POSTFIT',None)
	if channame=='allchannels' :
		data1Dhisto = initial_templates_file.Get(channel_names[1]+'__DATA').Clone()
		for j in range(2,len(channel_names)) :
			newname = channel_names[j]+'__DATA'
			data1Dhisto.Add(initial_templates_file.Get(newname).Clone())
	else :
		data1Dhisto = initial_templates_file.Get(newname).Clone()
	newdatatemp.make_from_1D_histo(data1Dhisto)
	newdatatemplates.append(newdatatemp)
	for j in range(newdatatemp.getHistoX().GetSize()) :
		if newdatatemp.getHistoX().IsBinUnderflow(j) or newdatatemp.getHistoX().IsBinOverflow(j) :
			continue
		if x_err_hs[i].GetBinContent(j)!=0. and newdatatemp.getHistoX().GetBinContent(j)!=0. :
			delta = newdatatemp.getHistoX().GetBinContent(j) - x_err_hs[i].GetBinContent(j)
			sigma = sqrt(newdatatemp.getHistoX().GetBinError(j)**2 + x_err_hs[i].GetBinError(j)**2)
			content = delta/sigma
			#print '		%s XHISTO BIN %d = (%.2f - %.2f)/sqrt(%.2f^2 + %.2f^2) = %.2f/%.2f = %.2f'%(channame,j,newdatatemp.getHistoX().GetBinContent(j),x_err_hs[i].GetBinContent(j),newdatatemp.getHistoX().GetBinError(j),x_err_hs[i].GetBinError(j),delta,sigma,content) #DEBUG
			x_resids[i].SetBinContent(j,content)
			maxxdeviations[i] = max(maxxdeviations[i],abs(content))
	for j in range(newdatatemp.getHistoY().GetSize()) :
		if newdatatemp.getHistoY().IsBinUnderflow(j) or newdatatemp.getHistoY().IsBinOverflow(j) :
			continue
		if y_err_hs[i].GetBinContent(j)!=0. and newdatatemp.getHistoY().GetBinContent(j)!=0. :
			delta = newdatatemp.getHistoY().GetBinContent(j) - y_err_hs[i].GetBinContent(j)
			sigma = sqrt(newdatatemp.getHistoY().GetBinError(j)**2 + y_err_hs[i].GetBinError(j)**2)
			content = delta/sigma
			#print '		%s YHISTO BIN %d = (%.2f - %.2f)/sqrt(%.2f^2 + %.2f^2) = %.2f/%.2f = %.2f'%(channame,j,newdatatemp.getHistoY().GetBinContent(j),y_err_hs[i].GetBinContent(j),newdatatemp.getHistoY().GetBinError(j),y_err_hs[i].GetBinError(j),delta,sigma,content) #DEBUG
			y_resids[i].SetBinContent(j,content)
			maxydeviations[i] = max(maxydeviations[i],abs(content))
	for j in range(newdatatemp.getHistoZ().GetSize()) :
		if newdatatemp.getHistoZ().IsBinUnderflow(j) or newdatatemp.getHistoZ().IsBinOverflow(j) :
			continue
		if z_err_hs[i].GetBinContent(j)!=0. and newdatatemp.getHistoZ().GetBinContent(j)!=0. :
			delta = newdatatemp.getHistoZ().GetBinContent(j) - z_err_hs[i].GetBinContent(j)
			sigma = sqrt(newdatatemp.getHistoZ().GetBinError(j)**2 + z_err_hs[i].GetBinError(j)**2)
			content = delta/sigma
			#print '		%s ZHISTO BIN %d = (%.2f - %.2f)/sqrt(%.2f^2 + %.2f^2) = %.2f/%.2f = %.2f'%(channame,j,newdatatemp.getHistoZ().GetBinContent(j),z_err_hs[i].GetBinContent(j),newdatatemp.getHistoZ().GetBinError(j),z_err_hs[i].GetBinError(j),delta,sigma,content) #DEBUG
			z_resids[i].SetBinContent(j,content)
			maxzdeviations[i] = max(maxzdeviations[i],abs(content))
	#reset stack maxima
	xmaxdata = newdatatemp.getHistoX().GetMaximum() 
	ymaxdata = newdatatemp.getHistoY().GetMaximum() 
	zmaxdata = newdatatemp.getHistoZ().GetMaximum()
	x_stacks[i].SetMaximum(1.02*max(x_stacks[i].GetMaximum(),xmaxdata+sqrt(xmaxdata)))
	y_stacks[i].SetMaximum(1.02*max(y_stacks[i].GetMaximum(),ymaxdata+sqrt(ymaxdata)))
	z_stacks[i].SetMaximum(1.02*max(z_stacks[i].GetMaximum(),zmaxdata+sqrt(zmaxdata)))
#Set histogram, MC error graph, and residual plot properties
for i in range(len(channel_names)) :
	x_resids[i].SetStats(0)
	x_resids[i].GetXaxis().SetLabelSize((0.05*0.72)/0.28); x_resids[i].GetXaxis().SetTitleOffset(1.0)
	x_resids[i].GetYaxis().SetLabelSize((0.05*0.72)/0.28); x_resids[i].GetYaxis().SetTitleOffset(0.4)
	x_resids[i].GetXaxis().SetTitleSize((0.72/0.28)*x_resids[i].GetXaxis().GetTitleSize())
	x_resids[i].GetYaxis().SetTitleSize((0.72/0.28)*x_resids[i].GetYaxis().GetTitleSize())
	maxx = 0.1+ceil(maxxdeviations[i])
	minx = -0.1-ceil(maxxdeviations[i])
	x_resids[i].GetYaxis().SetRangeUser(minx,maxx)
	x_resids[i].GetYaxis().SetNdivisions(503)
	x_resids[i].SetMarkerStyle(20)
	x_err_hs[i].SetFillStyle(3005)
	y_resids[i].SetStats(0)
	y_resids[i].GetXaxis().SetLabelSize((0.05*0.72)/0.28); y_resids[i].GetXaxis().SetTitleOffset(1.0)
	y_resids[i].GetYaxis().SetLabelSize((0.05*0.72)/0.28); y_resids[i].GetYaxis().SetTitleOffset(0.4)
	y_resids[i].GetXaxis().SetTitleSize((0.72/0.28)*y_resids[i].GetXaxis().GetTitleSize())
	y_resids[i].GetYaxis().SetTitleSize((0.72/0.28)*y_resids[i].GetYaxis().GetTitleSize())
	maxy = 0.1+ceil(maxydeviations[i])
	miny = -0.1-ceil(maxydeviations[i])
	y_resids[i].GetYaxis().SetRangeUser(miny,maxy)
	y_resids[i].GetYaxis().SetNdivisions(503)
	y_resids[i].SetMarkerStyle(20)
	y_err_hs[i].SetFillStyle(3005)
	z_resids[i].SetStats(0)
	z_resids[i].GetXaxis().SetLabelSize((0.05*0.72)/0.28); z_resids[i].GetXaxis().SetTitleOffset(1.0)
	z_resids[i].GetYaxis().SetLabelSize((0.05*0.72)/0.28); z_resids[i].GetYaxis().SetTitleOffset(0.4)
	z_resids[i].GetXaxis().SetTitleSize((0.72/0.28)*z_resids[i].GetXaxis().GetTitleSize())
	z_resids[i].GetYaxis().SetTitleSize((0.72/0.28)*z_resids[i].GetYaxis().GetTitleSize())
	maxz = 0.1+ceil(maxzdeviations[i])
	minz = -0.1-ceil(maxzdeviations[i])
	z_resids[i].GetYaxis().SetRangeUser(minz,maxz)
	z_resids[i].GetYaxis().SetNdivisions(503)
	z_resids[i].SetMarkerStyle(20)
	z_err_hs[i].SetFillStyle(3005)
#Build a legend
leg = TLegend(0.62,0.67,0.9,0.9)
for i in range(len(disttypes)) :
	disthistos[i].SetFillColor(disttypecolors[i])
	leg.AddEntry(disthistos[i],disttypes[i],'F')
newdatatemplates[0].getHistoX().SetMarkerStyle(20)
leg.AddEntry(newdatatemplates[0].getHistoX(),'DATA','PE')
leg.AddEntry(x_err_hs[0],'MC uncertainty','F')
#plot stacks with data overlaid and residuals
for i in range(len(channel_names)) :
	channame = channel_names[i]
	x_canvs[i].cd() 
	x_histo_pad=TPad(channame+'_x_histo_pad',channame+'_x_histo_pad',0,0.25,1,1)
	x_resid_pad=TPad(channame+'_x_residuals_pad',channame+'_x_residuals_pad',0,0,1.,0.25)
	x_histo_pad.SetCanvas(x_canvs[i]); x_resid_pad.SetCanvas(x_canvs[i])
	x_histo_pad.SetLeftMargin(0.16); x_histo_pad.SetRightMargin(0.05) 
	x_histo_pad.SetTopMargin(0.11);	 x_histo_pad.SetBottomMargin(0.02)
	x_histo_pad.SetBorderMode(0)
	x_resid_pad.SetLeftMargin(0.16); x_resid_pad.SetRightMargin(0.05)
	x_resid_pad.SetTopMargin(0.0);   x_resid_pad.SetBottomMargin(0.3)
	x_resid_pad.SetBorderMode(0)
	x_resid_pad.Draw(); x_histo_pad.Draw()
	x_histo_pad.cd(); 
	newdatatemplates[i].getHistoX().SetMarkerStyle(20)
	x_stacks[i].Draw(); newdatatemplates[i].getHistoX().Draw('SAME PE1X0'); x_err_hs[i].Draw('SAME E2'); x_stacks[i].GetXaxis().SetLabelOffset(999)
	leg.Draw()
	x_resid_pad.cd(); 
	x_resids[i].Draw('B')
	x_canvs[i].Update()
	outfile.cd()
	x_canvs[i].Write()
	
	y_canvs[i].cd() 
	y_histo_pad=TPad(channame+'_y_histo_pad',channame+'_y_histo_pad',0,0.3,1,1)
	y_resid_pad=TPad(channame+'_y_residuals_pad',channame+'_y_residuals_pad',0,0,1.,0.3)
	y_histo_pad.SetCanvas(y_canvs[i]); y_resid_pad.SetCanvas(y_canvs[i])
	y_histo_pad.SetLeftMargin(0.16); y_histo_pad.SetRightMargin(0.05) 
	y_histo_pad.SetTopMargin(0.11);	 y_histo_pad.SetBottomMargin(0.02)
	y_histo_pad.SetBorderMode(0)
	y_resid_pad.SetLeftMargin(0.16); y_resid_pad.SetRightMargin(0.05)
	y_resid_pad.SetTopMargin(0.0);   y_resid_pad.SetBottomMargin(0.3)
	y_resid_pad.SetBorderMode(0)
	y_resid_pad.Draw(); y_histo_pad.Draw()
	y_histo_pad.cd(); 
	newdatatemplates[i].getHistoY().SetMarkerStyle(20)
	y_stacks[i].Draw(); newdatatemplates[i].getHistoY().Draw('SAME PE1X0'); y_err_hs[i].Draw('SAME E2'); y_stacks[i].GetXaxis().SetLabelOffset(999)
	leg.Draw()
	y_resid_pad.cd(); 
	y_resids[i].Draw('B')
	y_canvs[i].Update()
	outfile.cd()
	y_canvs[i].Write()
	
	z_canvs[i].cd() 
	z_histo_pad=TPad(channame+'_z_histo_pad',channame+'_z_histo_pad',0,0.3,1,1)
	z_resid_pad=TPad(channame+'_z_residuals_pad',channame+'_z_residuals_pad',0,0,1.,0.3)
	z_histo_pad.SetCanvas(z_canvs[i]); z_resid_pad.SetCanvas(z_canvs[i])
	z_histo_pad.SetLeftMargin(0.16); z_histo_pad.SetRightMargin(0.05) 
	z_histo_pad.SetTopMargin(0.11);	 z_histo_pad.SetBottomMargin(0.02)
	z_histo_pad.SetBorderMode(0)
	z_resid_pad.SetLeftMargin(0.16); z_resid_pad.SetRightMargin(0.05)
	z_resid_pad.SetTopMargin(0.0);   z_resid_pad.SetBottomMargin(0.3)
	z_resid_pad.SetBorderMode(0)
	z_resid_pad.Draw(); z_histo_pad.Draw()
	z_histo_pad.cd(); 
	newdatatemplates[i].getHistoZ().SetMarkerStyle(20)
	z_stacks[i].Draw(); newdatatemplates[i].getHistoZ().Draw('SAME PE1X0'); z_err_hs[i].Draw('SAME E2'); z_stacks[i].GetXaxis().SetLabelOffset(999)
	leg.Draw()
	z_resid_pad.cd(); 
	z_resids[i].Draw('B')
	z_canvs[i].Update()
	outfile.cd()
	z_canvs[i].Write()