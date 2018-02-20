from ROOT import TFile, TCanvas, kRed, kBlue
from template import Template

#path to files holding the templates to test on
tfn = '../templates_new/templates_all.root'
tfn_aux = '../templates_new/templates_all_aux.root'

#channels of templates to compare pre/post shapes
channels = ['t1_muplus_SR',
			't1_muminus_SR',
			't1_elplus_SR',
			't1_elminus_SR',
			't2_muplus_SR',
			't2_muminus_SR',
			't2_elplus_SR',
			't2_elminus_SR',
			't3_muplus_SR',
			't3_muminus_SR',
			't3_elplus_SR',
			't3_elminus_SR',]

#names of templates to compare pre/post shapes
names = ['fgg',
		 'fqq',
		 'fbck',
		 'fwjets',
		 'fqp0',
		 'fqm0',
		 'fg0']

#list of all the complete template names to find
all_template_names = []
for channel in channels :
	for name in names :
		all_template_names.append(channel+'__'+name)

#get 1D templates from template file
all_orig_1D_templates = {}
tf = TFile.Open(tfn)
for tn in all_template_names :
	newhisto = tf.Get(tn)
	newhisto.SetDirectory(0)
	all_orig_1D_templates[tn] = newhisto
tf.Close()
#print 'all_orig_1D_templates = %s'%(all_orig_1D_templates) #DEBUG

#get original 3 template projections from the auxiliary file
all_orig_3D_templates = {}
tfa = TFile.Open(tfn_aux)
for tn in all_template_names :
	all_orig_3D_templates[tn]={}
	new_histo_x = tfa.Get(tn+'_x')
	new_histo_x.SetDirectory(0)
	all_orig_3D_templates[tn]['x'] = new_histo_x
	new_histo_y = tfa.Get(tn+'_y')
	new_histo_y.SetDirectory(0)
	all_orig_3D_templates[tn]['y'] = new_histo_y
	new_histo_z = tfa.Get(tn+'_z')
	new_histo_z.SetDirectory(0)
	all_orig_3D_templates[tn]['z'] = new_histo_z
tfa.Close()
#print 'all_orig_3D_templates = %s'%(all_orig_3D_templates) #DEBUG

#make new templates from the 1D histograms
all_new_templates = {}
for tn in all_orig_1D_templates :
	new_template = Template(tn+'__NEW',tn+'__NEW',None)
	new_template.make_from_1D_histo(all_orig_1D_templates[tn])
	all_new_templates[tn]=new_template

#open the output file
ofn = 'conversion_check_plots.root'
of = TFile.Open(ofn,'recreate')

#create canvases
all_canvs = {}
for tn in all_template_names :
	all_canvs[tn] = {}
	all_canvs[tn]['x'] = TCanvas(tn+'_x',tn+'_x',1100,900)
	all_canvs[tn]['y'] = TCanvas(tn+'_y',tn+'_y',1100,900)
	all_canvs[tn]['z'] = TCanvas(tn+'_z',tn+'_z',1100,900)

#plot on canvases
for tn in all_template_names :
	print 'plotting for template name %s...'%(tn)
	#x plots
	old_x_plot = all_orig_3D_templates[tn]['x']
	new_x_plot = all_new_templates[tn].getHistoX()
	old_x_plot.SetLineWidth(3); old_x_plot.SetLineColor(kRed); old_x_plot.SetLineStyle(9)
	new_x_plot.SetLineWidth(3); new_x_plot.SetLineColor(kBlue); new_x_plot.SetLineStyle(7)
	all_canvs[tn]['x'].cd()
	old_x_plot.Draw('HIST'); new_x_plot.Draw('SAME')
	of.cd()
	all_canvs[tn]['x'].Write()
	#y plots
	old_y_plot = all_orig_3D_templates[tn]['y']
	new_y_plot = all_new_templates[tn].getHistoY()
	old_y_plot.SetLineWidth(3); old_y_plot.SetLineColor(kRed); old_y_plot.SetLineStyle(9)
	new_y_plot.SetLineWidth(3); new_y_plot.SetLineColor(kBlue); new_y_plot.SetLineStyle(7)
	all_canvs[tn]['y'].cd()
	old_y_plot.Draw('HIST'); new_y_plot.Draw('SAME')
	of.cd()
	all_canvs[tn]['y'].Write()
	#z plots
	old_z_plot = all_orig_3D_templates[tn]['z']
	new_z_plot = all_new_templates[tn].getHistoZ()
	old_z_plot.SetLineWidth(3); old_z_plot.SetLineColor(kRed); old_z_plot.SetLineStyle(9)
	new_z_plot.SetLineWidth(3); new_z_plot.SetLineColor(kBlue); new_z_plot.SetLineStyle(7)
	all_canvs[tn]['z'].cd()
	old_z_plot.Draw('HIST'); new_z_plot.Draw('SAME')
	of.cd()
	all_canvs[tn]['z'].Write()
print 'Done.'

#close the output file
of.Close()

