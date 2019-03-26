from ROOT import TFile
from optparse import OptionParser
import os
from template_maker import Template_Group

##########								Parser Options								##########

parser = OptionParser()
#run options
parser.add_option('--input', 	type='string', action='store', default='input', 	dest='input', 	 help='Path to input file holding list of files to run on')
parser.add_option('--on_grid', 	type='string', action='store', default='no', 		dest='on_grid',  help='Path to on_grid file holding list of files to run on')
parser.add_option('--out_name', type='string', action='store', default='templates', dest='out_name', help='Name of output file that will have all the templates in it')
parser.add_option('--force_build_proctrees', type='string', action='store', default='no', dest='force_rb_ptrees', help='Set to "yes" to re build process trees even if they already exist')
parser.add_option('--n_procs', 	 type='int', action='store', default=10, dest='n_procs', help='how many processes can be allowed to run in parallel?')
#time-saving measures for later
parser.add_option('--sum_charges', type='string', action='store', default='no',  dest='sum_charges', help='Whether or not to integrate over the lepton charge in building templates')
parser.add_option('--include_mu',  type='string', action='store', default='yes', dest='include_mu',  help='Whether or not to include muons in the templates')
parser.add_option('--include_el',  type='string', action='store', default='yes', dest='include_el',  help='Whether or not to include electrons in the templates')
parser.add_option('--include_JEC', type='string', action='store', default='yes', dest='include_JEC', help='Whether or not to include JEC systematics in building templates')
parser.add_option('--include_sss', type='string', action='store', default='yes', dest='include_sss', help='Whether or not to include simple systematics in building templates')
parser.add_option('--topologies',  type='string', action='store', default='t1_t2_t3', dest='topologies', help='List of topologies to make templates for separated by underscores (default: t1_t2_t3)')
#fit parameters
parser.add_option('--Rqqbar', type='float', action='store', default=1.0, dest='Rqqbar', help='Rqqbar scale factor value')
parser.add_option('--Rbck',   type='float', action='store', default=1.0, dest='Rbck', 	help='Rbck scale factor value')
parser.add_option('--Rwjets', type='float', action='store', default=1.0, dest='Rwjets', help='Rwjets scale factor value')
parser.add_option('--Rqcd',   type='float', action='store', default=1.0, dest='Rqcd', 	help='Rqcd scale factor value')
parser.add_option('--Afb', 	  type='float', action='store', default=0.0, dest='Afb', 	help='Afb value')
parser.add_option('--mu', 	  type='float', action='store', default=0.0, dest='mu', 	help='mu value')
parser.add_option('--d', 	  type='float', action='store', default=0.0, dest='d', 		help='d value')
#fit parameters uncertainties
parser.add_option('--Rqqbar_sigma', type='float',  action='store', default=0.01,		 dest='Rqqbar_sigma', help='Rqqbar scale factor uncertainty')
parser.add_option('--Rbck_sigma', 	type='float',  action='store', default=0.1,		 dest='Rbck_sigma',   help='Rbck scale factor uncertainty')
parser.add_option('--Rwjets_sigma', type='float',  action='store', default=0.1,		 dest='Rwjets_sigma',   help='Rwjets scale factor uncertainty')
parser.add_option('--Rqcd_sigma', 	type='float',  action='store', default=0.3,		 dest='Rqcd_sigma',   help='Rqcd scale factor uncertainty')
parser.add_option('--Afb_sigma', 	type='float',  action='store', default=1.0,		 dest='Afb_sigma', 	  help='Afb uncertainty')
parser.add_option('--mu_sigma', 	type='float',  action='store', default=1.0,		 dest='mu_sigma', 	  help='mu uncertainty')
parser.add_option('--d_sigma', 		type='float',  action='store', default=1.0,		 dest='d_sigma', 	  help='d uncertainty')
#fix parameters?
parser.add_option('--fix_Rqqbar', type='string', action='store', default='yes', dest='fix_Rqqbar')
parser.add_option('--fix_Rbck',   type='string', action='store', default='yes', dest='fix_Rbck')
parser.add_option('--fix_Rwjets', type='string', action='store', default='yes', dest='fix_Rwjets')
parser.add_option('--fix_Rqcd',   type='string', action='store', default='yes', dest='fix_Rqcd')
parser.add_option('--fix_Afb', 	  type='string', action='store', default='no', dest='fix_Afb')
parser.add_option('--fix_mu', 	  type='string', action='store', default='no', dest='fix_mu')
parser.add_option('--fix_d', 	  type='string', action='store', default='no', dest='fix_d')
(options, args) = parser.parse_args()
topology_list = options.topologies.split('_')

##########								Script								##########
#Get the output file name and the auxiliary output file name
output_name = options.out_name
if not output_name.endswith('.root') :
	output_name+='.root'
output_name_aux = output_name.split('.root')[0]+'_aux.root'
#Create the template and auxiliary output files
outputfile = TFile(output_name,'recreate')
outputfile_aux = TFile(output_name_aux,'recreate')
#make the list of fit parameters based on input options
rqqbar = ('Rqqbar',options.Rqqbar,options.Rqqbar_sigma,options.fix_Rqqbar.lower()=='yes')
rbck   = ('Rbck',  options.Rbck,  options.Rbck_sigma,  options.fix_Rbck.lower()=='yes')
rwjets = ('Rwjets',options.Rwjets,options.Rwjets_sigma,options.fix_Rwjets.lower()=='yes')
rqcd   = ('Rqcd',  options.Rqcd,  options.Rqcd_sigma,  options.fix_Rqcd.lower()=='yes')
afb    = ('Afb',   options.Afb,   options.Afb_sigma,   options.fix_Afb.lower()=='yes')
mu 	   = ('mu',    options.mu,    options.mu_sigma,    options.fix_mu.lower()=='yes')
d 	   = ('d', 	   options.d, 	  options.d_sigma, 	   options.fix_d.lower()=='yes')
fit_parameter_tuple = (rqqbar,rbck,rwjets,rqcd,afb,mu,d)
#Start up the group of templates
print 'Creating template group'
templates = Template_Group(fit_parameter_tuple,options.out_name,options.sum_charges.lower()=='yes',options.include_mu.lower()=='yes',options.include_el.lower()=='yes',
							topology_list,options.include_JEC.lower()=='yes',options.include_sss.lower()=='yes',options.n_procs)
print 'Done'
#Build process trees from reconstructor trees if necessary
proc_tree_file_path = ''
if options.on_grid.lower() == 'yes' :
	proc_tree_file_path+='./tardir/'
else:
	proc_tree_file_path+='./'
proc_tree_file_path+=options.out_name+'_process_trees.root'
if (not os.path.isfile(proc_tree_file_path)) or options.force_rb_ptrees.lower()=='yes' :
	#Open the input file
	input_file_path =''
	if options.on_grid.lower() == 'yes' :
		input_file_path+='./tardir/'
	else:
		input_file_path+='./'
	input_file_path+=options.input
	if not '.txt' in options.input :
		input_file_path+='.txt'
	input_file = open(input_file_path,'r')
	#Read the files at the path
	print 'Reading from files'
	cmd = 'hadd -f '+options.out_name+'_process_trees.root '
	for line in input_file :
		if not line.startswith('#') :
			ttree_file_path = line.rstrip()
			print '	Adding files from '+ttree_file_path
			templates.add_file_to_processes(ttree_file_path)
			if os.path.isfile(ttree_file_path.split('/')[-1].rstrip('skim_all.root')+'_'+options.out_name+'_process_trees_all.root') :
				cmd+=ttree_file_path.split('/')[-1].rstrip('skim_all.root')+'_'+options.out_name+'_process_trees_all.root '
			print '	Done'
	#aggregate all the different process tree files
	print 'Aggregating all process tree files'
	os.system(cmd)
	os.system('rm -rf *_process_trees_all.root')
	print 'Done'
#Build the templates
print 'Building data-driven QCD templates'
templates.build_QCD_templates(proc_tree_file_path)
print 'Building DATA and MC-based templates'
templates.build_templates()
print 'Done'
#Save all the templates for this group
print 'Writing templates to file'
aux_obj_list = templates.get_list_of_auxiliary_objects()
outputfile_aux.cd()
for obj in aux_obj_list :
	obj.Write()
histo_list = templates.get_list_of_1D_histos()
outputfile.cd()
for histo in histo_list :
	histo.Write()
print 'Done'
