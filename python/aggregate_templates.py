#makes a "total templates" file from the powheg_TT nominal and systematic-split files

#imports
from ROOT import TFile
import os
from template import Template

#filename for total output
outfilenames = {'temp':'templates_powheg_aggregated_all.root','aux':'templates_powheg_aggregated_all_aux.root'}

#dict of filenames holding otherwise nominal templates and all other systematic variations
base_templates_filenames = {'temp':'../template_dirs/powheg_TT_ss_and_JEC/templates.root','aux':'../template_dirs/powheg_TT_ss_and_JEC/templates_aux.root'}

#channel names
channel_names = []
for t in ['t1','t2','t3'] :
	for leptype in ['muplus','muminus','elplus','elminus'] :
		for r in ['SR','WJets_CR'] :
			if t=='t3' and r=='WJets_CR' :
				continue
			channel_names.append(t+'_'+leptype+'_'+r)
#process names for which generation systematics matter
process_names = ['fqq','fgg','fqp0','fqm0','fq0','fq1','fq2','fg0','fg1','fg2','fg3','fg4','fbck']

#dict of filenames keyed by included systematic wiggle
pre = '../template_dirs/powheg_TT_'
post = '/templates'
#for those whose templates just need to be added to the final file
simple_add_files_dict = {'isrUp':{'temp':pre+'isrup'+post+'.root','aux':pre+'isrup'+post+'_aux.root'},
						 'isrDown':{'temp':pre+'isrdown'+post+'.root','aux':pre+'isrdown'+post+'_aux.root'},
						 'fsrUp':{'temp':pre+'fsrup'+post+'.root','aux':pre+'fsrup'+post+'_aux.root'},
						 'fsrDown':{'temp':pre+'fsrdown'+post+'.root','aux':pre+'fsrdown'+post+'_aux.root'},
						 'hdampUp':{'temp':pre+'hdampup'+post+'.root','aux':pre+'hdampup'+post+'_aux.root'},
						 'hdampDown':{'temp':pre+'hdampdown'+post+'.root','aux':pre+'hdampdown'+post+'_aux.root'},
						 'tuneUp':{'temp':pre+'tuneup'+post+'.root','aux':pre+'tuneup'+post+'_aux.root'},
						 'tuneDown':{'temp':pre+'tunedown'+post+'.root','aux':pre+'tunedown'+post+'_aux.root'},
						}
#for those whose templates all have to be compared with the nominal to decide on the envelope : /
CR_wiggle_files_dict = {'mpicr':{'temp':pre+'mpicr'+post+'.root','aux':pre+'mpicr'+post+'.root'},
						'qcdcr':{'temp':pre+'qcdcr'+post+'.root','aux':pre+'qcdcr'+post+'.root'},
						'gmcr':{'temp':pre+'gmcr'+post+'.root','aux':pre+'gmcr'+post+'.root'},
					   }

#all the template and auxiliary objects to add to the files
all_objs = {'temp':[],'aux':[]}

#all the filepointers to close at the very end
all_fps = {}

#open the output file
output_tempfile = TFile.Open(outfilenames['temp'],'recreate')
output_auxfile = TFile.Open(outfilenames['aux'],'recreate')
garbage_file = TFile.Open('template_aggregator_garbage.root','recreate')
garbage_file.cd()

#get objects from the nominal files
filetypes = ['temp','aux']
for ft in filetypes :
	print 'getting objects from %s...'%(base_templates_filenames[ft])
	all_fps['nominal_'+ft] = TFile(base_templates_filenames[ft],'r')
	lok = all_fps['nominal_'+ft].GetListOfKeys()
	for k in lok :
		#print '	reading key %s'%(k) #DEBUG
		newobj = k.ReadObj()
		all_objs[ft].append(newobj)
		#check if this is one of the templates in the nominal file that needs to be renormalized
		n=k.GetName()
		if n.find('__JES')!=-1 or n.find('__JER')!=-1 or n.find('__B_br_weight')!=-1 or n.find('__pdfas_weight')!=-1 :
			norm = all_fps['nominal_'+ft].Get(n.split('__')[0]+'__'+n.split('__')[1]).Integral()
			all_objs[ft][-1].Scale(norm/all_objs[ft][-1].Integral())
	print 'Done.'

#get objects from the files that just need to be renamed and added
for syswig in simple_add_files_dict.keys() :
	for ft in filetypes :
		print 'getting objects from %s...'%(simple_add_files_dict[syswig][ft])
		all_fps[syswig+'_'+ft] = TFile(simple_add_files_dict[syswig][ft],'r')
		#loop over channels and processes
		for cn in channel_names :
			for pn in process_names :
				#build the old and new template name
				old_template_name = cn+'__'+pn
				new_template_name = cn+'__'+pn+'__'+syswig
				#print old_template_name #DEBUG
				#get the nominal template to renormalize these to the same expected yield
				nominal_template = all_fps['nominal_temp'].Get(old_template_name)
				#if we're doing templates just get the old, rename, renormalize, and add it to the list of new
				if ft=='temp' :
					histo = all_fps[syswig+'_'+ft].Get(old_template_name)
					histo.SetName(new_template_name)
					histo.Scale(nominal_template.Integral()/histo.Integral())
					all_objs[ft].append(histo)
				#if we're doing the auxiliary objects we need the 3D template and its projections
				elif ft=='aux' :
					histo_3D = all_fps[syswig+'_'+ft].Get(old_template_name)
					histo_x  = all_fps[syswig+'_'+ft].Get(old_template_name+'_x')
					histo_y  = all_fps[syswig+'_'+ft].Get(old_template_name+'_y')
					histo_z  = all_fps[syswig+'_'+ft].Get(old_template_name+'_z')
					histo_3D.SetName(new_template_name)
					histo_x.SetName(new_template_name+'_x')
					histo_y.SetName(new_template_name+'_y')
					histo_z.SetName(new_template_name+'_z')
					histo_3D.Scale(nominal_template.Integral()/histo_3D.Integral())
					histo_x.Scale(nominal_template.Integral()/histo_x.Integral())
					histo_y.Scale(nominal_template.Integral()/histo_y.Integral())
					histo_z.Scale(nominal_template.Integral()/histo_z.Integral())
					all_objs[ft].append(histo_3D)
					all_objs[ft].append(histo_x)
					all_objs[ft].append(histo_y)
					all_objs[ft].append(histo_z)
		print 'Done.'

#handle the CR systematics templates, for which an envelope must be chosen and new templates created
#start with loop over channels/processes
print 'Building color reconnection templates...'
for cn in channel_names :
	for pn in process_names :
		#build the template names
		old_template_name = cn+'__'+pn
		new_template_name_up = cn+'__'+pn+'__crUp'
		new_template_name_dn = cn+'__'+pn+'__crDown'
		print '	doing templates for '+old_template_name
		#get the nominal template for this channel/process
		nominal_template = all_fps['nominal_temp'].Get(old_template_name)
		#get unrolled template histograms from all the different CR scenarios
		cr_sys_templates = []
		for crscenario in CR_wiggle_files_dict :
			if not crscenario+'_temp' in all_fps.keys() :
				all_fps[crscenario+'_temp'] = TFile(CR_wiggle_files_dict[crscenario]['temp'],'r')
			cr_sys_templates.append(all_fps[crscenario+'_temp'].Get(old_template_name))
		#renormalize them to the same event yield as the nominal template
		for cr_sys_template in cr_sys_templates :
			cr_sys_template.Scale(nominal_template.Integral()/cr_sys_template.Integral())
		#create new Template objects for color reconnection up/down scenarios and get their unfilled 1D templates
		new_up_template_obj = Template(new_template_name_up,old_template_name+' color reconnection up template',None)
		new_dn_template_obj = Template(new_template_name_dn,old_template_name+' color reconnection down template',None)
		new_up_template_histo, trash1 = new_up_template_obj.convertTo1D()
		new_dn_template_histo, trash2 = new_dn_template_obj.convertTo1D()
		#loop over the bins 
		for i in range(1,nominal_template.GetNbinsX()+1) :
			#set the up/down templates' bin contents equal to +/- the average deviation from the nominal
			nomcont=nominal_template.GetBinContent(i)
			devs = []
			for cr_sys_temp in cr_sys_templates :
				devs.append(cr_sys_temp.GetBinContent(i)-nomcont)
			avgdev=sum(devs)/float(len(devs))
			new_up_template_histo.SetBinContent(i,nomcont+avgdev)
			new_dn_template_histo.SetBinContent(i,nomcont-avgdev)
		#add the 1D template to the list of templates
		all_objs['temp'].append(new_up_template_histo)
		all_objs['temp'].append(new_dn_template_histo)
		#rebuild the templates from the modified 1D histograms and add the auxiliary objects to the lists
		new_up_template_obj.make_from_1D_histo(new_up_template_histo)
		new_dn_template_obj.make_from_1D_histo(new_dn_template_histo)
		all_objs['aux']+=new_up_template_obj.getHistos()
		all_objs['aux']+=new_dn_template_obj.getHistos()
print 'Done.'

#save all the objects to the appropriate output file
print 'Writing templates and objects to files...'
for temphisto in all_objs['temp'] :
	temphisto.SetDirectory(output_tempfile)
for auxobj in all_objs['aux'] :
	auxobj.SetDirectory(output_auxfile)
print '	writing template file...'
#write and close the output files
output_tempfile.Write(); output_tempfile.Close()
print '	done. writing aux file...'
output_auxfile.Write(); output_auxfile.Close()
print 'Done.'
#close all the files that were opened
for fp in all_fps.values() :
	fp.Close()
#write, close, and delete the garbage file
print 'Getting rid of garbage file...'
garbage_file.Write(); garbage_file.Close()
os.system('rm -rf template_aggregator_garbage.root')
print 'Done!'
