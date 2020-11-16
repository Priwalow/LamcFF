#!/usr/bin/python

#Python wrapper for running GSIM using mcproduzh backend scripts
#Uses analysis config file to set location of Evtgen and GSIM dirs

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
from random import randrange

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

#Function to ask to overwrite a folder
def delete_dir(dir_name,skim):
    if os.path.exists(dir_name):
        reply = raw_input("GSIM mdst already exists at "+skim+" do you want to overwrite it? [y/[n]] ").lower().strip()
        print reply
        if reply=='yes' or reply=='y':
		if(options.debug):
			print 'rm -r '+dir_name
		else:
        		os.system('rm -r '+dir_name)
        else:
            print "Aborting..."
            sys.exit()

#Set config file
if not exists('analysis.cfg'):
	print 'analysis.cfg not found. Run Setup.py to generate it.'
	sys.exit()
config = ConfigParser.ConfigParser()
config.read('analysis.cfg')

#Setup argument parser
usage = 'usage: %prog [options] EVTGEN_DIR'
parser = PosOptionParser(usage)
parser.add_option('-i','--interactive', action='store_true', help='run the script in interactive mode (by default scripts are submitted to KEK batch queues)')
parser.add_option('-d','--debug', action='store_true', help='Debug test. Print script command line calls only, does not actually run anything.')
parser.add_positional_argument(Option('--EVTGEN_DIR', action='store_true',help='  Directory containing the data to run GSIM over (relative to the Evtgen output directory set in the config). This will be the name of the decay file with a version number appended'))
(options, args) = parser.parse_args()

#Check number of arguments
if len(args) != 1:
	parser.error("Incorrect number of arguments. Run 'Signal_GSIM.py --help' for more information")

#Get directory for Evtgen files from the config file
gen_dir = config.get('Output_Locations', 'evt_dir', 0)
if(gen_dir!=''):
	gen_dir = expanduser(gen_dir)
if not exists(gen_dir):					#check valid dir
        print 'Evtgen input file location in config file is not a valid folder'
        sys.exit()

#Get directory for GSIM output from the config file
gsim_dir = config.get('Output_Locations', 'gsim_dir', 0)
if(gsim_dir!=''):
	gsim_dir = expanduser(gsim_dir)
if not exists(gsim_dir):				#check valid dir
	reply =  raw_input('GSIM Output file location in config file '+gsim_dir+' does not exist, create it now? y/[n] ').lower().strip()
        if reply=='yes' or reply=='y':
                print 'Creating directory '+gsim_dir
                if(options.debug):
                        print 'mkdir -p '+gsim_dir
                else:
                        os.system('mkdir -p '+gsim_dir)
        else:
                print 'Aborting...'
                sys.exit()

#Get input folder from user argument and add to Evtgen input dir
shortgenpath = args[0]
evtgenpath = join(gen_dir,shortgenpath)
if not exists(evtgenpath):                                 #check valid dir
        print 'Evtgen folder for analysis is not a valid folder'
        sys.exit()

#Define output dir and add version number if folder already exists
GSIMout = path.normpath(join(gsim_dir,path.basename(evtgenpath))) #attach input folder name to GSIM out
delete_dir(GSIMout,shortgenpath)
if(options.debug):
	print 'mkdir '+GSIMout
else:
	os.system('mkdir '+GSIMout)

#list of files in input path
mdst = [ f for f in os.listdir(evtgenpath) if path.isfile(join(evtgenpath,f)) ]

#os.chdir(join(os.getcwd(),'core_scripts/mcproduzh/gsim'))
os.chdir('core_scripts/mcproduzh/gsim')

#submit scripts to queue for all .mdsts in input dir
for m in mdst:
	if string.lower(path.splitext(m)[1]) == '.mdst':	#check file is mdst
		#get exp from input file name
		exp = string.split(m,'_')[2]

		#Manipulate exp Variables for basf
		exp = exp.strip("0")
                expno = ''.join([i for i in exp if i.isdigit()])
                expdiv = ''.join([i for i in exp if i.isalpha()])
                exp_2 = expno.zfill(2)+expdiv
		exp_6 = expno.zfill(6)

		#location of basf and dat files to read
		basf = "basf/mcprod5-e"+exp_2+".basf"
		gsim = "gsim/gsim."+expno+".dat"
		if not exists(basf):
			print 'basf setting \"'+basf+'\" not found'
			sys.exit()
		if not exists(gsim):
			print 'gsim file \"'+gsim+'\" not found'
			sys.exit()

		#Generate random run number for BG
		run = str(randrange(110))

		#Set belle level based on exp no. 
		if int(expno) <= 27:
			BL = "b20030807_1600"
		elif int(expno) <= 71:
			BL = "b20090127_0910"
		else:
			print "Exp $expno is not implemented in this package"
			sys.exit()
		if(options.debug):
			if(options.interactive):
                		print './GSIM_test.sh '+BL+' '+exp+' '+basf+' '+exp_6+' '+run+' '+join(evtgenpath,m)+' '+join(GSIMout,path.splitext(m)[0])
			else:
				print 'bsub -q l ./GSIM_test.sh '+BL+' '+exp+' '+basf+' '+exp_6+' '+run+' '+join(evtgenpath,m)+' '+join(GSIMout,path.splitext(m)[0])
		else:
			if(options.interactive):
				os.system('./GSIM.sh '+BL+' '+exp+' '+basf+' '+exp_6+' '+run+' '+join(evtgenpath,m)+' '+join(GSIMout,path.splitext(m)[0]))
			else:
				os.system('bsub -q l ./GSIM.sh '+BL+' '+exp+' '+basf+' '+exp_6+' '+run+' '+join(evtgenpath,m)+' '+join(GSIMout,path.splitext(m)[0]))
