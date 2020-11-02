#!/usr/bin/python

#Python port of mcproduzh used for running evtgen with events scaled for experiment luminosities
#Uses analysis config file to set locations of Evtgen and GSIM dirs
#Tristan Bloomfield (t.bloomfield@student.unimelb.edu.au)
#26/02/2015

import ConfigParser
from optparse import OptionParser, Option, IndentedHelpFormatter
import os
import os.path as path
from os.path import exists
from os.path import join
from os.path import expanduser
import string
import sys
import time

class PosOptionParser(OptionParser):
    def format_help(self, formatter=None):
        class Positional(object):
            def __init__(self, args):
                self.option_groups = []
                self.option_list = args

        positional = Positional(self.positional)
        formatter = IndentedHelpFormatter()
        formatter.store_option_strings(positional)
        output = ['\n', formatter.format_heading("  Positional Arguments")]
        formatter.indent()
        pos_help = [formatter.format_option(opt) for opt in self.positional]
        pos_help = [line.replace('--','  ') for line in pos_help]
        output += pos_help
        return OptionParser.format_help(self, formatter) + ''.join(output)

    def add_positional_argument(self, option):
        try:
            args = self.positional
        except AttributeError:
            args = []
        args.append(option)
        self.positional = args

    def set_out(self, out):
        self.out = out

#**********************Parse Arguments and config file********************

if not exists('analysis.cfg'):
        print 'analysis.cfg not found. Run Setup.py to generate it.'
        sys.exit()
config = ConfigParser.ConfigParser()
config.read('analysis.cfg')

#Setup argument parser
usage = 'usage: %prog [options] DECAY_FILE TOT_EVTS JOB_EVTS'
parser = PosOptionParser(usage)
parser.add_option('-i','--interactive', action='store_true', help='Run the script in interactive mode (by default scripts are submitted to KEK batch queues).')
parser.add_option('-d','--debug', action='store_true', help='Debug test. Print script command line calls only, does not actually run anything.')
parser.add_positional_argument(Option('--DECAY_FILE', action='store_true', help='    Name of the decay file you wish to generate MC with.'))
parser.add_positional_argument(Option('--TOT_EVTS', action='store_true', help='    Total number of events you wish to be generated.'))
parser.add_positional_argument(Option('--JOB_EVTS', action='store_true', help='    Number of events to generate per job.'))


(options, args) = parser.parse_args()
if len(args) != 3:
        parser.error("Incorrect number of arguments. Run Signal_Evtgen.py --help for more information")

#Decay File
decay_file = args[0]
decay_folder = config.get('Evtgen','DEC_DIR')
if(decay_folder != ''):
	decay_folder = expanduser(decay_folder)
dec = join(decay_folder,decay_file)

if not path.isfile(dec):
	print "Decay file: \'"+dec+"\' does not exist."
	sys.exit()

#Job Numbers
tot_num = int(args[1])
#if not isinstance( tot_num, int ):
#	print "Total number of events is not an integer"
#	sys.exit()
job_num = int(args[2])
#if not isinstance( job_num, int ):
#	print 'Max events per job is not an integer'
#	sys.exit()

#Set MC type file and read in luminosity file
Lum_loc = config.get('Evtgen', 'nBB', 0)
if(Lum_loc!=''):
	Lum_loc = expanduser(Lum_loc)
if not path.isfile(Lum_loc):
	print "Luminosity file: \'"+Lum_loc+"\' does not exist."
	sys.exit()
Lum_file = open(Lum_loc)
Lum = Lum_file.readlines()
Lum_file.close()
MC = config.get('Evtgen', 'MC_type', 0)
if(MC!=''):
	MC = expanduser(MC)
if not path.isfile(MC):
	print "MC .conf file: \'"+MC+"\' does not exist."
#Output dir, ask for create dir if it doesn't exist
gen_dir = config.get('Output_Locations', 'evt_dir', 0)
if(gen_dir!=''):
		gen_dir = expanduser(gen_dir)
if not exists(gen_dir):
	reply =  raw_input('Output file location in config file '+gen_dir+' does not exist, create it now? y/[n] ').lower().strip()
	if reply=='yes' or reply=='y':
		print 'Creating directory '+gen_dir
		if(options.debug):
			print 'mkdir -p '+gen_dir
		else:
			os.system('mkdir -p '+gen_dir)
	else:
		print 'Aborting...'
		sys.exit()

#***********Check if output dir exists, if so, add version number***********
out_dir = path.split(path.splitext(decay_file)[0])[1]
out_dir_abs = path.normpath(join(gen_dir,out_dir))
out_dir_ver = out_dir_abs
ver = 1
while exists(out_dir_ver):
	ver = ver+1
	out_dir_ver = out_dir_abs+'_'+str(ver)
print 'Creating directory '+out_dir_ver
if(options.debug):
	print "mkdir -p "+out_dir_ver
else:
	os.system("mkdir -p "+out_dir_ver)

#***************Sum the luminosities for each experiment********************
sum = float(0)
for l in Lum:
	if not string.split(l)[0].startswith('#'):		#if line not a comment
		sum = float(sum)+float(string.split(l)[1])
print 'Total Luminosity = '+str(sum)
print 'Will scale number of events for each experiment accordingly. Total sum of generated events will be '+str(tot_num)

#**************************Submit scaled scripts****************************

for l in Lum:
	if l[0] != '#':
		#Get the experiment number and luminosity
		exp = string.split(l)[0]
		explum = float(string.split(l)[1])
		
		#Calculate number of events to generate
		evts = int(round(explum/sum*float(tot_num)))
		print '**** Exp = '+exp+': generating '+str(evts)+' events'

		#Set string of experiment number only, subdivions only and value with 2 digits
		expno = ''.join([i for i in exp if i.isdigit()])
		expdiv = ''.join([i for i in exp if i.isalpha()])
		basfexp = expno.zfill(2)+expdiv

		#Set counter for number of events submitted and number of jobs submitted
		startedEvts = int(0)
		job = int(0)
		
		#Loop for jobs with full number of events
		while (startedEvts <= (evts-job_num)):
			print 'Starting job '+str(job)+' with '+str(job_num)+' events'

			#Set Output file name and submit job
			File_Name = 'evtgen_exp_'+basfexp+'_'+path.splitext(decay_file)[0]+'-'+str(job)
			Output = join(out_dir_ver,File_Name)
                        if (options.debug):
                                if (options.interactive):
                                        print './core_scripts/mcproduzh/evtgen/EvtgenBASFScript.sh '+MC+' '+dec+' '+str(job_num)+' '+Output
                                else:
                                        print 'bsub -q b_a ./core_scripts/mcproduzh/evtgen/EvtgenBASFScript.sh '+MC+' '+dec+' '+str(job_num)+' '+Output
                        else:
                                if (options.interactive):
                                        os.system('./core_scripts/mcproduzh/evtgen/EvtgenBASFScript.sh '+MC+' '+dec+' '+str(job_num)+' '+Output)
                                else:
                                        os.system('bsub -q b_a ./core_scripts/mcproduzh/evtgen/EvtgenBASFScript.sh '+MC+' '+dec+' '+str(job_num)+' '+Output)	
			#Iterate submitted jobs and events
			startedEvts = startedEvts+job_num
			job = str(int(job)+1)

		#Submitted final job if there are left over events
		remEvts = evts-startedEvts
		if (remEvts>0):
                        print 'Starting job '+str(job)+' with '+str(remEvts)+' events'
			
			#Set Output file name and submit job
                        File_Name = 'evtgen_exp_'+basfexp+'_'+path.splitext(decay_file)[0]+'-'+str(job)
                        Output = join(out_dir_ver,File_Name)
			if (options.debug):
				if (options.interactive):
					print './core_scripts/mcproduzh/evtgen/EvtgenBASFScript.sh '+MC+' '+dec+' '+str(remEvts)+' '+Output
				else:
					print 'bsub -q b_a ./core_scripts/mcproduzh/evtgen/EvtgenBASFScript.sh '+MC+' '+dec+' '+str(remEvts)+' '+Output
                        else:
				if (options.interactive):
                                        os.system('./core_scripts/mcproduzh/evtgen/EvtgenBASFScript.sh '+MC+' '+dec+' '+str(remEvts)+' '+Output)
                                else:
                                        os.system('bsub -q b_a ./core_scripts/mcproduzh/evtgen/EvtgenBASFScript.sh '+MC+' '+dec+' '+str(remEvts)+' '+Output)
		#wait for 2 seconds between exp submits (needed?)
		#time.sleep(1)
	
