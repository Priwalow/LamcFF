#!/usr/bin/python

#Python script to setup config file for MC package.

#Tristan Bloomfield (t.bloomfield@student.unimelb.edu.au)
#30/01/2015

import ConfigParser
import os
import os.path as path
from os.path import exists
from os.path import join
from os.path import expanduser
import string
import sys

#Set config file
config = ConfigParser.ConfigParser()
config.read('analysis.cfg')

print "This script will set up the configuration file for use with the Evtgen and GSIM scripts."
#Output Existing values
if exists('analysis.cfg'):
	print 'Configuration files already exists. Current values for are:'
	print 'Path to evtgen files: '+config.get('Output_Locations','evt_dir',0)
	print 'Path to GSIM files: '+config.get('Output_Locations','gsim_dir',0)
	print 'Decay file directory: '+config.get('Evtgen','DEC_DIR',0)
	print 'Luminosity file: '+config.get('Evtgen','nBB',0)
	print 'MC Config file: '+config.get('Evtgen','MC_type',0)
	change = raw_input('Would you like to change these values? [y/[n]]: ').lower().strip()
	if change=='yes' or change=='y':
		print 'Overwriting existing values'
	else:
		print 'Exiting'
		sys.exit()
#Add Sections
else:
	print 'Setting values for the first time'
	config.add_section('Output_Locations')
	config.add_section('Evtgen')

#______Set Output dirs_________

#Set evtgen directory
evt = expanduser(raw_input("Path to the directory you would like evtgen files to be saved to: "))
while not exists(evt):
	reply = raw_input('Directory \"'+evt+'\" does not exists, create it now? [y/[n]]: ').lower().strip()
	if reply=='yes' or reply=='y':
		print "Creating directory "+evt
		os.system('mkdir '+evt)
	else:
		evt = expanduser(raw_input('Specify a different directory for output of Evtgen files: '))
config.set('Output_Locations','evt_dir',evt)

#Set GSIM directory
gsim = expanduser(raw_input("Path of the directory you would like GSIM files to be saved to: "))
while not exists(gsim):
        reply2 = raw_input('Directory \"'+gsim+'\" does not exists, create it now? [y/[n]]: ').lower().strip()
        if reply2=='yes' or reply2=='y':
                print "Creating directory "+gsim
                os.system('mkdir '+gsim)
        else:
                gsim = expanduser(raw_input('Specify a different directory for output of Evtgen files: '))
config.set('Output_Locations','gsim_dir',gsim)

#Set decay file directory
decay = expanduser(raw_input('Path of the directory where your decay files for evtgen are located: '))
while not exists(decay):
        reply3 = raw_input('Directory \"'+decay+'\" does not exists, create it now? [y/[n]]: ').lower().strip()
        if reply3=='yes' or reply3=='y':
                print "Creating directory "+evt
                os.system('mkdir '+evt)
        else:
                decay = expanduser(raw_input('Specify a different directory for output of Evtgen files: '))
config.set('Evtgen','DEC_DIR',decay)

#_______Set Evtgen configuration files_______
detailed_conf = raw_input('Would you like to change the default luminosity and MC config (Y4S)? [y/n]: ').lower().strip()
while (detailed_conf!='yes') and (detailed_conf!='y') and (detailed_conf!='no') and (detailed_conf!='n'):
	detailed_conf = raw_input('Please answer yes/no: ').lower().strip()
if detailed_conf=='yes' or detailed_conf=='y':
	lum = expanduser(raw_input('Location of luminosity file: '))
	while not exists(lum):
		lum = raw_input('\"'+lum+'\" does not exist, enter a valid file: ')
	config.set('Evtgen','nBB',lum)
	mct = expanduser(raw_input('Location of MC .conf file: '))
	while not exists(mct):
		mct = raw_input('\"'+mct+'\" does not exist, enter a valid file: ')
	config.set('Evtgen','MC_type',mct)
if detailed_conf=='no' or detailed_conf=='n':
	config.set('Evtgen','nBB','core_scripts/mcproduzh/evtgen/nBB-Y4S.txt')
	config.set('Evtgen','MC_type','core_scripts/mcproduzh/evtgen/Y4S.conf')
	print "Using default setup for Y4S evtgen"

#_______Write Config File_______
configfile =  open('analysis.cfg', 'w+')
try:
	config.write(configfile)
finally:
	configfile.close()
print "MC Production config file is now set up"
sys.exit()
