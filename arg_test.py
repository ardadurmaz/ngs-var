#!/usr/bin/python

######################################################
import sys, os, argparse
import os.path, re
import glob

dir_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(dir_path)
import ngs_classes
######################################################

def check_inputs(inputs):
	
	if inputs.verbose: print "[INFO] Checking provided arguments"
	
	if inputs.trimmer == 'CUTADAPT':
		if(os.path.exists(inputs.cutadapt)):
			if inputs.verbose: print "[INFO:CUTADAPT] Available"
		else:
			if inputs.verbose: print "[ERROR:CUTADAPT] Not Available"
			return(None)
	elif inputs.trimmer == 'AFTERQC':
		if(os.path.exists(inputs.afterqc)):
			if inputs.verbose: print "[INFO:AFTERQC] Available"
		else:
			if inputs.verbose: print "[ERROR:AFTERQC] Not Available"
			return(None)
	elif inputs.trimmer == 'SKEWER':
		if(os.path.exists(inputs.skewer)):
			if inputs.verbose: print "[INFO:SKEWER] Available"
		else:
			if inputs.verbose: print "[ERROR:SKEWER] Not Available"
			return(None)
	
	if not (inputs.aligner == 'BWA' or inputs.aligner == 'BOWTIE2'):
		if inputs.verbose: print "[ERROR] Name of the aligner does not match available tools"
		return(None)
	else:
		if inputs.aligner == 'BWA':
			if(os.path.exists(inputs.aligner)):
				if inputs.verbose: print "[INFO:BWA] Available"
			else:
				if inputs.verbose: print "[ERROR:BWA] Not Available"
				return(None)
		elif inputs.aligner == 'BOWTIE2':
			if(os.path.exists(inputs.bowtie2)):
				if inputs.verbose: print "[INFO:BOWTIE2] Available"
			else:
				if inputs.verbose: print "[ERROR:BOWTIE2] Not Available"
				return(None)

	if inputs.workflow == 'GERMLINESNVINDEL':
		if not (inputs.tool == 'HAPLOTYPECALLER' or inputs.tool == 'STRELKA2'):
			if inputs.verbose: print "[ERROR] Name of the caller does not match available tools"
			return(None)
	elif inputs.workflow == 'SOMATICSNVINDEL':
		if not (inputs.tool == 'STRELKA2'):
			if inputs.verbose: print "[ERROR] Name of the caller does not match available tools"
			return(None)
	else:
		if inputs.verbose: print "[ERROR] Name of the workflow does not match available workflows"
	
	for tool in [inputs.gatk,
				 inputs.picard,
				 inputs.sambamba]:
		if(os.path.exists(os.path.abspath(tool)) and os.access(tool, os.R_OK)):
			if inputs.verbose: print "[INFO] Available: %s" % (tool)
		else:
			if inputs.verbose: print "[ERROR] Not Available: %s, exiting" % (tool)
			return(None)
	
	for knownSite in inputs.knownsites:
		if(os.path.exists(os.path.abspath(knownSite)) and os.access(knownSite, os.R_OK)):
			if inputs.verbose: print "[INFO] Available %s" % (knownSite)
		else:
			if inputs.verbose: print "[ERROR] Not Available: %s, exiting" % (knownSite)
			return(None)
		
	if inputs.verbose: print "[INFO] Done."
		
	return(True)

def parse_config(config_file, inputs):
	
	config_count = 0
	if not(os.path.exists(os.path.abspath(config_file))):
		print "\n!Error: Configuration file does not exists\n"
		return(True)
	fh = open(config_file, "rU")
	try:
		for l in fh:
			ss = l.rstrip("\n")
			match = re.search('^(\w+)\:\t(\S+)$', ss)
			if(match):
				key = str.upper(match.group(1))
				val = match.group(2)
				#print "KEY: %s, VAL: %s" % (key, val)
				if key == 'REFERENCE':
					inputs.reference = val
					if inputs.verbose: print "[INFO] Parsed reference resource"
					config_count = config_count + 1
				elif key == 'GATK':
					inputs.gatk = val
					if inputs.verbose: print "[INFO] Parsed gatk resource"
					config_count = config_count + 1
				elif key == 'PICARD':
					inputs.picard = val
					if inputs.verbose: print "[INFO] Parsed picard resource"
					config_count = config_count + 1
				elif key == 'SAMTOOLS':
					inputs.samtools = val
					if inputs.verbose: print "[INFO] Parsed samtools resource"
					config_count = config_count + 1
				elif key ==	'BCFTOOLS':
					inputs.bcftools = val
					if inputs.verbose: print "[INFO] Parsed bcftools resource"
					config_count = config_count + 1
				elif key == 'SAMBAMBA':
					inputs.sambamba = val
					if inputs.verbose: print "[INFO] Parsed sambamba resource"
					config_count = config_count + 1
				elif key == 'BOWTIE2':
					inputs.bowtie2 = val
					if inputs.verbose: print "[INFO] Parsed bowtie2 resource"
					config_count = config_count + 1
				elif key == 'BOWTIE2BUILD':
					inputs.bowtie2build = val
					if inputs.verbose: print "[INFO] Parsed bowtie2build resource"
					config_count = config_count + 1
				elif key == 'BWA':
					inputs.bwa = val
					if inputs.verbose: print "[INFO] Parsed bwa resource"
					config_count = config_count + 1
				elif key == 'CUTADAPT':
					inputs.cutadapt = val
					if inputs.verbose: print "[INFO] Parsed cutadapt resource"
					config_count = config_count + 1
				elif key == 'SKEWER':
					inputs.skewer = val
					if inputs.verbose: print "[INFO] Parsed skewer resource"
					config_count = config_count + 1
				elif key == 'AFTERQC':
					inputs.afterqc = val
					if inputs.verbose: print "[INFO] Parsed afterqc resource"
					config_count = config_count + 1
				elif key == 'STRELKABIN':
					inputs.strelkabin = val
					if inputs.verbose: print "[INFO] Parsed strelkabin resource"
					config_count = config_count + 1
				elif key == 'MANTABIN':
					inputs.mantabin = val
					if inputs.verbose: print "[INFO] Parsed mantabin resource"
					config_count = config_count + 1
				elif key == 'KNOWNSITES':
					inputs.knownsites = val.split(",")
					if inputs.verbose: print "[INFO] Parsed knownsites resource"
					config_count = config_count + 1
				else:
					if inputs.verbose: print "[WARNING] Unrecognized resource in configuration file"
	finally:
		fh.close()
	
	if(config_count > 0):
		if inputs.verbose: print "[INFO] Processed %d value pairs in configuration file" % (config_count)
		return(False)
	else:
		return(None)


def get_inputs():
	
	## Get Arguments ##
	parser = argparse.ArgumentParser(prog='NGS',description='*** Epigenetiks WXS Pipeline ***')
	
	parser.add_argument('--dir',default='NGS_WD',help='Path to working directory (Will be created if doesn\'t exists)')
	parser.add_argument('--config',default='ngs.config',help="Configuration file")
	parser.add_argument('--threads',default=4,type=int,help="Number of threads to use")
	parser.add_argument('--title',default='NGS_Analysis',help='Title of the run')

	parser.add_argument('--verbose',action="store_true",help='Verbosity')
	parser.add_argument('--dry',action="store_true",help='Dry-Run')
	parser.add_argument('--clear',action='store_true',default=False,help='Clear workspace')
	parser.add_argument('--exome',action='store_true',default=False,help='Whole Exome Data')

	parser.add_argument('--in_file',default=None,required=True,help='Targets file containing sample ids and associated fastq files')
	parser.add_argument('--trimmer',default=None,required=False,help='Name of the trimmer <cutadapt|afterqc|skewer>')
	parser.add_argument('--aligner',default='bwa',help='Name of the aligner <bowtie2|bwa>')
	parser.add_argument('--workflow',default='germlinesnvindel',help='Type of analysis to run <GermlineSNVIndel|SomaticSNVIndel>')
	parser.add_argument('--tool',default='Strelka2',help='Name of the caller <HaplotypeCaller|Strelka2')
		
	parser.add_argument('--bed',default=None,help='Bed file for regions')
	args = parser.parse_args()
	
	
	## Collect Inputs ##
	inputs = ngs_classes.ngsInputs()
	inputs.dir = args.dir
	inputs.threads = args.threads
	inputs.title = args.title
	
	inputs.verbose = args.verbose
	inputs.dry = args.dry
	inputs.clear = args.clear
	inputs.exome = args.exome
	
	inputs.in_file = args.in_file
	inputs.trimmer = str.upper(args.trimmer) if not args.trimmer is None else None
	inputs.aligner = str.upper(args.aligner)
	inputs.workflow = str.upper(args.workflow)
	inputs.tool = str.upper(args.tool)

	## Get Config Data ##
	if parse_config(args.config, inputs) is None:
		print "!Error in parsing configuration file\n"
		return(None)
	
	
	if inputs.verbose:
		print "\n"
		print "\/" * 40
		
		for key, val in vars(inputs).iteritems():
			if(key == 'knownsites'):
				for site in val:
					print "\/ %s: %s" % ("KnownSite", site)
			else:
				print "\/ %s: %s" % (key, val)
				
		print "\/" * 40
		print "\n"
		
	return inputs





if __name__ == '__main__':
	
	print "\n\n *** NGS Analysis ***\n\n"

	inputs = get_inputs();
	if inputs is None:
		print "!! Error in reading inputs"
		sys.exit(1)

	if check_inputs(inputs) is None:
		print "[ERROR] Error in arguments"
		sys.exit(1)
		
		
		
		
	sys.exit(0)