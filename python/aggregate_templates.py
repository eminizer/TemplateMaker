#makes a "total templates" file from the powheg_TT nominal and systematic-split files

#imports
from ROOT import TFile
import os
from template import Template

#filename for total output
outfilenames = {'temp':'templates_powheg_dynamic_binning_aggregated_sep_sample_symmetrized_all.root','aux':'templates_powheg_dynamic_binning_aggregated_sep_sample_symmetrized_all_aux.root'}

#dict of filenames holding otherwise nominal templates and all other systematic variations
base_templates_filenames = {'temp':'templates_powheg_dynamic_binning_JEC_and_simple_sys_all.root','aux':'templates_powheg_dynamic_binning_JEC_and_simple_sys_all_aux.root'}

ZEROED_BIN_CONTENT=0.000001

#correction function
def correct_template(nom,corr) :
  #loop over the bins
  for i in range(1,nom.GetNbinsX()+1) :
    #if this bin is zeroed in the nominal template, zero it in the shifted template
    if nom.GetBinContent(i)==ZEROED_BIN_CONTENT :
      corr.SetBinContent(i,ZEROED_BIN_CONTENT)
    #otherwise if it's zeroed in a shifted template, set that bin's content to the nominal content
    elif abs(corr.GetBinContent(i)-ZEROED_BIN_CONTENT)<0.000005 :
      corr.SetBinContent(i,nom.GetBinContent(i))
  #return the shifted template
  return corr

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
simple_add_files_dict = {'isr':{'Up':{'temp':pre+'isrup'+post+'.root','aux':pre+'isrup'+post+'_aux.root'},
							    'Down':{'temp':pre+'isrdown'+post+'.root','aux':pre+'isrdown'+post+'_aux.root'},},
						 'fsr':{'Up':{'temp':pre+'fsrup'+post+'.root','aux':pre+'fsrup'+post+'_aux.root'},
							    'Down':{'temp':pre+'fsrdown'+post+'.root','aux':pre+'fsrdown'+post+'_aux.root'},},
						 'hdamp':{'Up':{'temp':pre+'hdampup'+post+'.root','aux':pre+'hdampup'+post+'_aux.root'},
								  'Down':{'temp':pre+'hdampdown'+post+'.root','aux':pre+'hdampdown'+post+'_aux.root'},},
						 'tune':{'Up':{'temp':pre+'tuneup'+post+'.root','aux':pre+'tuneup'+post+'_aux.root'},
								 'Down':{'temp':pre+'tunedown'+post+'.root','aux':pre+'tunedown'+post+'_aux.root'},},
						}
#simple_add_files_dict = {'isrUp':{'temp':pre+'isrup'+post+'.root','aux':pre+'isrup'+post+'_aux.root'},
#						 'isrDown':{'temp':pre+'isrdown'+post+'.root','aux':pre+'isrdown'+post+'_aux.root'},
#						 'fsrUp':{'temp':pre+'fsrup'+post+'.root','aux':pre+'fsrup'+post+'_aux.root'},
#						 'fsrDown':{'temp':pre+'fsrdown'+post+'.root','aux':pre+'fsrdown'+post+'_aux.root'},
#						 'hdampUp':{'temp':pre+'hdampup'+post+'.root','aux':pre+'hdampup'+post+'_aux.root'},
#						 'hdampDown':{'temp':pre+'hdampdown'+post+'.root','aux':pre+'hdampdown'+post+'_aux.root'},
#						 'tuneUp':{'temp':pre+'tuneup'+post+'.root','aux':pre+'tuneup'+post+'_aux.root'},
#						 'tuneDown':{'temp':pre+'tunedown'+post+'.root','aux':pre+'tunedown'+post+'_aux.root'},
#						}
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
filetypes = ['temp']#,'aux']
for ft in filetypes :
	print 'getting objects from %s...'%(base_templates_filenames[ft])
	all_fps['nominal_'+ft] = TFile(base_templates_filenames[ft],'r')
	lok = all_fps['nominal_'+ft].GetListOfKeys()
	for k in lok :
		n=k.GetName()
		#check if this is one of the templates in the file that needs to be corrected
		if n.find('__JES')!=-1 or n.find('__JER')!=-1 :
			nsplit=n.split('__')
			nom=all_fps['nominal_'+ft].Get(nsplit[0]+'__'+nsplit[1])
			all_objs[ft].append(correct_template(nom,k.ReadObj()))
		else :
			#print '	reading key %s'%(k) #DEBUG
			newobj = k.ReadObj()
			all_objs[ft].append(newobj)
		#	#check if this is one of the templates in the nominal file that needs to be renormalized
		#	if n.find('__JES')!=-1 or n.find('__JER')!=-1 or n.find('__B_br_weight')!=-1 or n.find('__pdfas_weight')!=-1 :
		#		norm = all_fps['nominal_'+ft].Get(n.split('__')[0]+'__'+n.split('__')[1]).Integral()
		#		all_objs[ft][-1].Scale(norm/all_objs[ft][-1].Integral())
	print 'Done.'

#handle the other separate sample systematics, for which we want an envelope around the up/down scenarios
print 'Building separate sample systematic templates...'
for sn in simple_add_files_dict.keys() :
	for cn in channel_names :
		for pn in process_names :
			#build the template names
			old_template_name = cn+'__'+pn
			new_template_name_up = cn+'__'+pn+'__'+sn+'Up'
			new_template_name_dn = cn+'__'+pn+'__'+sn+'Down'
			print '	doing templates for '+old_template_name
			#get the nominal template for this channel/process
			nominal_template = all_fps['nominal_temp'].Get(old_template_name)
			#get template histograms from the up/down scenarios
			if not sn+'_up' in all_fps.keys() :
				all_fps[sn+'_up'] = TFile(simple_add_files_dict[sn]['Up']['temp'],'r')
			orig_up_template = all_fps[sn+'_up'].Get(old_template_name)
			if not sn+'_dn' in all_fps.keys() :
				all_fps[sn+'_dn'] = TFile(simple_add_files_dict[sn]['Down']['temp'],'r')
			orig_dn_template = all_fps[sn+'_dn'].Get(old_template_name)
			#renormalize to the same event yield as the nominal template
			orig_up_template.Scale(nominal_template.Integral()/orig_up_template.Integral())
			orig_dn_template.Scale(nominal_template.Integral()/orig_dn_template.Integral())
			#create new Template objects for up/down scenarios and get their unfilled 1D templates
			new_up_template_obj = Template(new_template_name_up,old_template_name+' '+sn+' up template',None)
			new_dn_template_obj = Template(new_template_name_dn,old_template_name+' '+sn+' down template',None)
			new_up_template_histo = new_up_template_obj.convertTo1D()
			new_dn_template_histo = new_dn_template_obj.convertTo1D()
			#loop over the bins 
			for i in range(1,nominal_template.GetNbinsX()+1) :
				nomcont=nominal_template.GetBinContent(i)
				upcont=orig_up_template.GetBinContent(i)
				dncont=orig_dn_template.GetBinContent(i)
				#set the up/down templates' bin contents equal to nominal +/- 1/2 up-down 
				shift = 0.5*(upcont-dncont)
				new_up_template_histo.SetBinContent(i,nomcont+shift)
				new_dn_template_histo.SetBinContent(i,nomcont-shift)
				##set the up/down templates' bin contents equal to their original values 
				#new_up_template_histo.SetBinContent(i,upcont)
				#new_dn_template_histo.SetBinContent(i,dncont)
			#correct the templates
			new_up_template_histo = correct_template(nominal_template,new_up_template_histo)
			new_dn_template_histo = correct_template(nominal_template,new_dn_template_histo)
			#add the 1D template to the list of templates
			all_objs['temp'].append(new_up_template_histo)
			all_objs['temp'].append(new_dn_template_histo)
			#rebuild the templates from the modified 1D histograms and add the auxiliary objects to the lists
			new_up_template_obj.make_from_1D_histo(new_up_template_histo)
			new_dn_template_obj.make_from_1D_histo(new_dn_template_histo)
			all_objs['aux']+=new_up_template_obj.getHistos()
			all_objs['aux']+=new_dn_template_obj.getHistos()
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
		new_up_template_histo = new_up_template_obj.convertTo1D()
		new_dn_template_histo = new_dn_template_obj.convertTo1D()
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
		#correct the templates
		new_up_template_histo = correct_template(nominal_template,new_up_template_histo)
		new_dn_template_histo = correct_template(nominal_template,new_dn_template_histo)
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
#for auxobj in all_objs['aux'] :
#	auxobj.SetDirectory(output_auxfile)
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
