#!/usr/bin/python

######################################################
import sys, os, argparse
import os.path, re
import glob

dir_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(dir_path)
import ngs_classes
from ngs_functions import print_log, run_command
######################################################

def check_inputs(inputs):
	
	if inputs.verbose: print_log(inputs, "[INFO] Checking provided arguments")
	if inputs.trimmer == 'FASTP':
		if(os.path.exists(inputs.cutadapt)):
			if inputs.verbose: print_log(inputs, "[INFO:FASTP] Available")
		else:
			if inputs.verbose: print_log(inputs, "[ERROR:FASTP] Not Available")
			return(None)
	elif inputs.trimmer == 'CUTADAPT':
		if(os.path.exists(inputs.cutadapt)):
			if inputs.verbose: print_log(inputs, "[INFO:CUTADAPT] Available")
		else:
			if inputs.verbose: print_log(inputs, "[ERROR:CUTADAPT] Not Available")
			return(None)
	elif inputs.trimmer == 'AFTERQC':
		if(os.path.exists(inputs.afterqc)):
			if inputs.verbose: print_log(inputs, "[INFO:AFTERQC] Available")
		else:
			if inputs.verbose: print_log(inputs, "[ERROR:AFTERQC] Not Available")
			return(None)
	elif inputs.trimmer == 'SKEWER':
		if(os.path.exists(inputs.skewer)):
			if inputs.verbose: print_log(inputs, "[INFO:SKEWER] Available")
		else:
			if inputs.verbose: print_log(inputs, "[ERROR:SKEWER] Not Available")
			return(None)
	
	if not (inputs.aligner == 'BWA' or inputs.aligner == 'BOWTIE2'):
		if inputs.verbose: print_log(inputs, "[ERROR] Name of the aligner does not match available tools")
		return(None)
	else:
		if inputs.aligner == 'BWA':
			if(os.path.exists(inputs.aligner)):
				if inputs.verbose: print_log(inputs, "[INFO:BWA] Available")
			else:
				if inputs.verbose: print_log(inputs, "[ERROR:BWA] Not Available")
				return(None)
		elif inputs.aligner == 'BOWTIE2':
			if(os.path.exists(inputs.bowtie2)):
				if inputs.verbose: print_log(inputs, "[INFO:BOWTIE2] Available")
			else:
				if inputs.verbose: print_log(inputs, "[ERROR:BOWTIE2] Not Available")
				return(None)

	if inputs.workflow == 'GERMLINESNVINDEL':
		if not (inputs.tool == 'HAPLOTYPECALLER' or inputs.tool == 'STRELKA2'):
			if inputs.verbose: print_log(inputs, "[ERROR] Name of the caller does not match available tools")
			return(None)
	elif inputs.workflow == 'SOMATICSNVINDEL':
		if not (inputs.tool == 'STRELKA2'):
			if inputs.verbose: print_log(inputs, "[ERROR] Name of the caller does not match available tools")
			return(None)
	else:
		if inputs.verbose: print_log(inputs, "[ERROR] Name of the workflow does not match available workflows")
	
	for tool in [inputs.gatk,
				 inputs.sambamba]:
		if(os.path.exists(os.path.abspath(tool)) and os.access(tool, os.R_OK)):
			if inputs.verbose: print_log(inputs, f"[INFO] Available: {tool}")
		else:
			if inputs.verbose: print_log(inputs, f"[ERROR] Not Available: {tool}, exiting")
			return(None)
	
	for knownSite in inputs.knownsites:
		if(os.path.exists(os.path.abspath(knownSite)) and os.access(knownSite, os.R_OK)):
			if inputs.verbose: print_log(inputs, f"[INFO] Available {knownSite}")
		else:
			if inputs.verbose: print_log(inputs, f"[ERROR] Not Available: {knownSite}, exiting")
			return(None)
		
	if inputs.verbose: print_log(inputs, "[INFO] Done.")
		
	return(True)

def parse_config(config_file, inputs):
	
	config_count = 0
	if not(os.path.exists(os.path.abspath(config_file))):
		print_log(inputs, "\n!Error: Configuration file does not exists\n")
		return(True)
	fh = open(config_file, "rU")
	try:
		for l in fh:
			ss = l.rstrip("\n")
			match = re.search('^(\w+)\:\t(\S+)$', ss)
			if(match):
				key = str.upper(match.group(1))
				val = match.group(2)
				if key == 'REFERENCE':
					inputs.reference = val
					if inputs.verbose: print_log(inputs, "[INFO] Parsed reference resource")
					config_count = config_count + 1
				elif key == 'GATK':
					inputs.gatk = val
					if inputs.verbose: print_log(inputs, "[INFO] Parsed gatk resource")
					config_count = config_count + 1
				elif key == 'SAMTOOLS':
					inputs.samtools = val
					if inputs.verbose: print_log(inputs, "[INFO] Parsed samtools resource")
					config_count = config_count + 1
				elif key ==	'BCFTOOLS':
					inputs.bcftools = val
					if inputs.verbose: print_log(inputs, "[INFO] Parsed bcftools resource")
					config_count = config_count + 1
				elif key == 'SAMBAMBA':
					inputs.sambamba = val
					if inputs.verbose: print_log(inputs, "[INFO] Parsed sambamba resource")
					config_count = config_count + 1
				elif key == 'BOWTIE2':
					inputs.bowtie2 = val
					if inputs.verbose: print_log(inputs, "[INFO] Parsed bowtie2 resource")
					config_count = config_count + 1
				elif key == 'BOWTIE2BUILD':
					inputs.bowtie2build = val
					if inputs.verbose: print_log(inputs, "[INFO] Parsed bowtie2build resource")
					config_count = config_count + 1
				elif key == 'BWA':
					inputs.bwa = val
					if inputs.verbose: print_log(inputs, "[INFO] Parsed bwa resource")
					config_count = config_count + 1
				elif key == 'FASTP':
					inputs.fastp = val
					if inputs.verbose: print_log(inputs, "[INFO] Parsed fastp resource")
					config_count = config_count + 1
				elif key == 'CUTADAPT':
					inputs.cutadapt = val
					if inputs.verbose: print_log(inputs, "[INFO] Parsed cutadapt resource")
					config_count = config_count + 1
				elif key == 'SKEWER':
					inputs.skewer = val
					if inputs.verbose: print_log(inputs, "[INFO] Parsed skewer resource")
					config_count = config_count + 1
				elif key == 'AFTERQC':
					inputs.afterqc = val
					if inputs.verbose: print_log(inputs, "[INFO] Parsed afterqc resource")
					config_count = config_count + 1
				elif key == 'STRELKABIN':
					inputs.strelkabin = val
					if inputs.verbose: print_log(inputs, "[INFO] Parsed strelkabin resource")
					config_count = config_count + 1
				elif key == 'MANTABIN':
					inputs.mantabin = val
					if inputs.verbose: print_log(inputs, "[INFO] Parsed mantabin resource")
					config_count = config_count + 1
				elif key == 'KNOWNSITES':
					inputs.knownsites = val.split(",")
					if inputs.verbose: print_log(inputs, "[INFO] Parsed knownsites resource")
					config_count = config_count + 1
				elif key == 'CNVKIT':
					inputs.cnvkit = val
					if inputs.verbose: print_log(inputs, "[INFO] Parsed Cnvkit resource")
					config_count = config_count + 1
				else:
					if inputs.verbose: print_log(inputs, f"[WARNING] Unrecognized resource in configuration file: {key}")
	finally:
		fh.close()
	
	if(config_count > 0):
		if inputs.verbose: print_log(inputs, f"[INFO] Processed {config_count} value pairs in configuration file")
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
		print_log(inputs, "!Error in parsing configuration file\n")
		return(None)
	
	
	if inputs.verbose:
		print_log(inputs, "\n")
		print_log(inputs, "\/" * 40)
		
		for key, val in vars(inputs).items():
			if(key == 'knownsites'):
				for site in val:
					print_log(inputs, f"\/ KnownSite: {site}")
			else:
				print_log(inputs, f"\/ {key}: {val}")
				
		print_log(inputs, "\/" * 40)
		print_log(inputs, "\n")
		
	return inputs





if __name__ == '__main__':
	
	print_log(inputs, "\n\n *** NGS Analysis ***\n\n")

	inputs = get_inputs();
	if inputs is None:
		print_log(inputs, "!! Error in reading inputs")
		sys.exit(1)

	if check_inputs(inputs) is None:
		print(inputs, "[ERROR] Error in arguments")
		sys.exit(1)
		
		
		
		
	sys.exit(0)
