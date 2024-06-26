#!/usr/bin/python

######################################################
import sys, os, argparse
import os.path, re
import glob

dir_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(dir_path)
import ngs_classes, ngs_germline, ngs_somatic, ngs_preprocess, ngs_rna
from ngs_functions import print_log, run_command
######################################################
	
def ngs_clean(inputs):
	
	if not inputs.dry:
		if inputs.upload:
			if inputs.verbose: print_log(inputs, f"[INFO:UPLOAD] sudo aws s3 sync {inputs.dir} s3://rnaseq-bucket/{inputs.title}")
			os.system(f"aws s3 sync {inputs.dir} s3://rnaseq-bucket/{inputs.title}")
	return 0

def read_targets(inputs):
	
	in_data = []
	fh = open(inputs.in_file)
	try:
		indx = {}    
		_header = True
		_mixed = False
		for l in fh:
			if not(re.match('^#', l) or re.match('^\\s+$', l)):
				local = l.rstrip("\n").split('\t')
				if _header:
					for i in range(len(local)):
						if re.search('read1', str.lower(local[i])):
							indx["Read1"] = i
						elif re.search('read2', str.lower(local[i])):
							indx["Read2"] = i
						elif re.search('lane', str.lower(local[i])):
							indx["Lane"] = i
						elif re.search('id', str.lower(local[i])):
							indx["ID"] = i
						elif re.search('plat', str.lower(local[i])):
							indx["Platform"] = i
						elif re.search('lib', str.lower(local[i])):
							indx["Library"] = i
						if re.search('type', str.lower(local[i])):
							if(inputs.workflow == "GERMLINESNVINDEL" and inputs.verbose): print_log(inputs, "[WARNING] Workflow type "\
							"is GERMLINESNVINDEL but tumor/normal type column found (Using only Normal tags)")
							indx["Type"] = i
							_mixed = True
					_header = False
				else:
					joined = "\t".join(local)
					if len(list([x for x in local if re.search('^[^\\s]+$', x)])) != len(indx):
						print_log(inputs, f'[WARNING:FORMAT] {joined}')
					else:
						if(str.upper(inputs.workflow) == "GERMLINESNVINDEL"):
							if not(_mixed):
								in_data.append(ngs_classes.WxsGermline(_r1 = local[indx["Read1"]],
																	   _r2 = local[indx["Read2"]],
																	   _lane = local[indx["Lane"]],
																	   _id = local[indx["ID"]],
																	   _plat = local[indx["Platform"]],
																	   _lib = local[indx["Library"]]))
							else:
								if(str.lower(local[indx["Type"]]) == 'normal'):
									in_data.append(ngs_classes.WxsGermline(_r1 = local[indx["Read1"]],
																		   _r2 = local[indx["Read2"]],
																		   _lane = local[indx["Lane"]],
																		   _id = local[indx["ID"]],
																		   _plat = local[indx["Platform"]],
																		   _lib = local[indx["Library"]]))
									
						elif(str.upper(inputs.workflow) == "SOMATICSNVINDEL"):
							in_data.append(ngs_classes.WxsSomatic(_r1 = local[indx["Read1"]],
																  _r2 = local[indx["Read2"]],
																  _lane = local[indx["Lane"]],
																  _id = local[indx["ID"]],
																  _plat = local[indx["Platform"]],
																  _lib = local[indx["Library"]],
																  _type = local[indx["Type"]]))
							
						elif(str.upper(inputs.workflow) == "RNASEQ"):
							in_data.append(ngs_classes.WxsRNA(_r1 = local[indx["Read1"]],
															  _r2 = local[indx["Read2"]],
															  _lane = local[indx["Lane"]],
															  _id = local[indx["ID"]],
															  _plat = local[indx["Platform"]],
															  _lib = local[indx["Library"]]))
	finally:
		fh.close()

	if inputs.verbose:
		for r in in_data:
			print_log(inputs, "\n    *** Run Info ***")
			for arg, val in vars(r).items():
				if val:
					print_log(inputs, f"  [{arg} --> {val}]")   
		print_log(inputs, "\n") 

	targets_count = 0
	if str.upper(inputs.workflow) == "GERMLINESNVINDEL":
		for entry in in_data:
			if hasattr(entry, 'type'):
				if str.lower(entry.type) == 'normal':
					targets_count += 1
			else:
				targets_count += 1
	else:
		targets_count = len(in_data)

	if targets_count == 0:
		print_log(inputs, "[INFO] Read 0 targets... exiting")
		sys.exit(1)
	else:
		print_log(inputs, f"[INFO] Targets read: {targets_count}")
	
	print_log(inputs, "\n") 
	return in_data

def check_inputs(inputs):
    if inputs.verbose: print_log(inputs, "[INFO] Checking provided arguments")
	
	## Check Trimmer ##
    if inputs.trimmer == 'CUTADAPT':
        if not(os.path.exists(inputs.cutadapt)): raise ngs_classes.ngsExcept("[ERROR:CUTADAPT] Not Available")
    elif inputs.trimmer == 'AFTERQC':
        if not(os.path.exists(inputs.afterqc)): raise ngs_classes.ngsExcept("[ERROR:AFTERQC] Not Available")
    elif inputs.trimmer == 'SKEWER':
        if not(os.path.exists(inputs.skewer)): raise ngs_classes.ngsExcept("[ERROR:SKEWER] Not Available")
    elif inputs.trimmer == 'FASTP':
        if not(os.path.exists(inputs.fastp)): raise ngs_classes.ngsExcept("[ERROR:FASTP] Not Available")
    else:
        raise ngs_classes.ngsExcept("[ERROR:TRIMMER] Name of the trimmer does not match available tools")
	
	## Check Aligner ##
    if inputs.aligner == 'BWA':
        if not(os.path.exists(inputs.bwa)): raise ngs_classes.ngsExcept("[ERROR:BWA] Not Available")
    elif inputs.aligner == 'STAR':
        if not(os.path.exists(inputs.star)): raise ngs_classes.ngsExcept("[ERROR:STAR] Not Available")
    else:
        raise ngs_classes.ngsExcept("[ERROR:ALIGNER] Name of the aligner does not match available tools")

    if inputs.workflow == 'GERMLINESNVINDEL':
        if(inputs.tool == 'HAPLOTYPECALLER'):
            if not(os.path.exists(inputs.gatk) or os.access(inputs.gatk,os.R_OK)): raise ngs_classes.ngsExcept("[ERROR] Failed to access GATK tool for HaplotypeCaller")
        elif(inputs.tool == 'STRELKA2'):
            if not(os.path.exists(inputs.strelkabin) or os.access(inputs.strelkabin,os.R_OK)): raise ngs_classes.ngsExcept("[ERROR] Failed to access Strelka2 directory")
            elif not(os.path.exists(inputs.mantabin) or os.access(inputs.mantabin,os.R_OK)): raise ngs_classes.ngsExcept("[ERROR] Failed to access Manta directory")
        else:
            raise ngs_classes.ngsExcept("[ERROR:TOOL] Name of the caller tool does not match available tools")
    elif inputs.workflow == 'SOMATICSNVINDEL':
        if(inputs.tool == 'MUTECT2'):
            if not(os.path.exists(inputs.gatk) or os.access(inputs.gatk,os.R_OK)): raise ngs_classes.ngsExcept("[ERROR] Failed to access GATK tool for MuTect2")
        elif(inputs.tool == 'STRELKA2'):
            if not(os.path.exists(inputs.strelkabin) or os.access(inputs.strelkabin,os.R_OK)): raise ngs_classes.ngsExcept("[ERROR] Failed to access Strelka2 directory")
        else:
            raise ngs_classes.ngsExcept("[ERROR:TOOL] Name of the caller tool does not match available tools")
    elif inputs.workflow == 'RNASEQ':
        pass
    else:
        raise ngs_classes.ngsExcept("[ERROR:WORKFLOW] Provided workflow is not available")
	
	## Check SAMBAMBA ##
    if not(os.path.exists(inputs.sambamba) or os.access(inputs.sambamba,os.R_OK)): raise ngs_classes.ngsExcept("[ERROR] Failed to access Sambamba tool")
	
	## Check Known Sites ##
    if(inputs.bqsr):
        for knownSite in inputs.knownsites:
            if not(os.path.exists(os.path.abspath(knownSite)) or os.access(knownSite, os.R_OK)): raise ngs_classes.ngsExcept(f"[ERROR] Failed to access resource {knownSite}")

	## Create Working Dir ##
    try:
        prep_wd(inputs)
    except ngs_classes.ngsExcept as err:
        print_log(inputs, err.msg)
        raise ngs_classes.ngsExcept("[ERROR] Input Check Failed")
        
    if inputs.verbose: print_log(inputs, "[INFO] Done.")
    
    return 0

def parse_config(config_file, inputs):
	
	config_count = 0
	if not(os.path.exists(os.path.abspath(config_file)) or os.access(os.path.abspath(config_file),os.R_OK)):
		raise ngs_classes.ngsExcept("[ERROR] Failed to access .config file")
	
	fh = open(config_file)
	try:
		for l in fh:
			ss = l.rstrip("\n")
			match = re.search('^(\\w+)\\:\t(\\S+)$', ss)
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
				elif key == 'SAMBAMBA':
					inputs.sambamba = val
					if inputs.verbose: print_log(inputs, "[INFO] Parsed sambamba resource")
					config_count = config_count + 1
				elif key == 'BWA':
					inputs.bwa = val
					if inputs.verbose: print_log(inputs, "[INFO] Parsed bwa resource")
					config_count = config_count + 1
				elif key == 'STAR':
					inputs.star = val
					if inputs.verbose: print_log(inputs, "[INFO] Parsed STAR resource")
					config_count = config_count + 1
				elif key == 'SAMTOOLS':
					inputs.samtools = val
					if inputs.verbose: print_log(inputs, "[INFO] Parsed samtools resource")
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
				elif key == 'RSCRIPT':
					inputs.rscript = val
					if inputs.verbose: print_log(inputs, "[INFO] Parsed rscript resource")
					config_count = config_count + 1
				elif key == 'SALMON':
					inputs.salmon = val
					if inputs.verbose: print_log(inputs, "[INFO] Parsed salmon resource")
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
				elif key == "ANNOTATION":
					inputs.annotation = val
					if inputs.verbose: print_log(inputs, "[INFO] Parsed annotation resource")
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
		return 0
	else:
		raise ngs_classes.ngsExcept("[ERROR] Failed to parse .config file")

def get_inputs(inputs):

	## Get Arguments ##
	parser = argparse.ArgumentParser(prog='NGS',description="NGS Pipeline")
	
	parser.add_argument('--dir',default='NGS_WD',help='Name of working directory')
	parser.add_argument('--config',default='ngs.config',help="Configuration file")
	parser.add_argument('--threads',default=4,type=int,help="Number of threads to use")
	parser.add_argument('--title',default='NGS_Analysis',help='Title of the run')
	parser.add_argument('--log',default="log.txt",help='Name of log file')
	
	parser.add_argument('--verbose',action="store_true",help='Verbosity')
	parser.add_argument('--dry',action="store_true",help='Dry-Run')
	parser.add_argument('--clear',action='store_true',default=False,help='Clear workspace')
	parser.add_argument('--upload',action='store_true',default=False,help='Upload intermediate results to S3 storage (aws must be set for non-root user)')
	parser.add_argument('--exome',action='store_true',default=False,help='Whole Exome Data')
	parser.add_argument('--bqsr',action='store_true',default=False,help='Do BQSR')
	
	parser.add_argument('--trimmer',default='cutadapt',required=False,help='Name of the trimmer <Cutadapt|Fastp|Skewer|AfterQC>')
	parser.add_argument('--aligner',default='bwa',help='Name of the aligner <BWA|STAR>')
	parser.add_argument('--tool',default='HaplotypeCaller',help='Name of the caller <HaplotypeCaller|Strelka2|MuTect2>')
	parser.add_argument('--workflow',default='GermlineSNVIndel',help='Type of analysis to run <GermlineSNVIndel|SomaticSNVIndel|RNASeq>')
	parser.add_argument('--in_file',default='None',required=True,help='Targets file containing sample ids and associated fastq files')
	
	parser.add_argument('--bed',default=None,help='Bed file for regions')
	parser.add_argument('--cbed',default=None,help='Bgzip compressed and tabix indexed bed file for regions')

	args = parser.parse_args()
	
	## Collect Inputs ##
	inputs.dir = os.path.abspath(args.dir)
	inputs.threads = args.threads
	inputs.title = args.title
	inputs.log = os.path.abspath(args.log)
	inputs.bed = args.bed
	inputs.cbed = args.cbed	
	
	inputs.verbose = args.verbose
	inputs.dry = args.dry
	inputs.clear = args.clear
	inputs.exome = args.exome
	inputs.bqsr = args.bqsr
	inputs.upload = args.upload
	
	inputs.in_file = args.in_file
	inputs.trimmer = str.upper(args.trimmer) if args.trimmer else None
	inputs.aligner = str.upper(args.aligner)
	inputs.workflow = str.upper(args.workflow)
	inputs.tool = str.upper(args.tool)

	prep_log(inputs)
	## Get Config Data ##
	try:
		parse_config(args.config, inputs)
	except ngs_classes.ngsExcept as err:
		print_log(inputs, err.msg)
		raise ngs_classes.ngsExcept("[ERROR] Failed to parse inputs")

	if inputs.verbose:
		print_log(inputs, "\\n")
		print_log(inputs, "\\/" * 40)

		for key, val in vars(inputs).items():
			if(key == 'knownsites'):
				for site in val:
					print_log(inputs, f"\\/ KnownSite: {site}")
			else:
				print_log(inputs, f"\\/ {key}: {val}")

		print_log(inputs, "\\/" * 40)
		print_log(inputs, "\n")

	return inputs


def prep_wd(inputs):
	if inputs.clear:
		for l_dir in ['readsTrimmed', 'readsAligned' ,'bamProcessed', 'bamMerged', 'recalib', 'bamCalib']:
			ll_dir = inputs.dir + '/' + l_dir
			if os.path.exists(ll_dir): os.system(f"rm -ir {ll_dir}")

	for l_dir in ['readsTrimmed', 'readsAligned' ,'bamProcessed', 'bamMerged', 'recalib', 'bamCalib']:
		ll_dir = inputs.dir + '/' + l_dir
		if not(os.path.exists(ll_dir)):
			try:
				os.makedirs(ll_dir)
			except OSError:
				raise ngs_classes.ngsExcept(f"[ERROR] Failed to create directory {ll_dir}")
			
	return 0

def prep_log(inputs):
	with open(inputs.log, 'w') as f:
		f.write("")

if __name__ == '__main__':
	
	inputs = ngs_classes.ngsInputs()
	try:
		get_inputs(inputs)
	except ngs_classes.ngsExcept as err:
		print_log(inputs, err.msg)
		sys.exit(1)

	print_log(inputs, "\n\n *** NGS Analysis ***\n\n")
	print_log(inputs, "[INFO:MAN] Create directories: readsTrimmed, readsAligned, bamProcessed, bamMerged, recalib, bamCalib\n\n")

	try:
		check_inputs(inputs)
	except ngs_classes.ngsExcept as err:
		print_log(inputs, err.msg)
		sys.exit(1)
	else:
		print_log(inputs, "[INFO] Input check complete.")

	targets_data = read_targets(inputs)

	## Index ##
	try:
		ngs_preprocess.ngs_index(inputs)
	except ngs_classes.ngsExcept as err:
		print_log(inputs, err.msg)
		sys.exit(1)
	else:
		print_log(inputs, "[INFO:PASS] Reference indexing complete.")
	finally:
		ngs_clean(inputs)
	
	
	## Trim ##
	try:
		ngs_preprocess.ngs_trim(inputs,targets_data)
	except ngs_classes.ngsExcept as err:
		print_log(inputs, err.msg)
		sys.exit(1)
	else:
		print_log(inputs, "[INFO:PASS] Read trimming complete.")
	finally:
		ngs_clean(inputs)

	## Align ##
	try:
		ngs_preprocess.ngs_align(inputs,targets_data)
	except ngs_classes.ngsExcept as err:
		print_log(inputs, err.msg)
		sys.exit(1)
	else:
		print_log(inputs, "[INFO:PASS] Read alignment complete.")
	finally:
		ngs_clean(inputs)

	## Sort & Mark ##    
	try:
		ngs_preprocess.ngs_mark_sort(inputs,targets_data)
	except ngs_classes.ngsExcept as err:
		print_log(inputs, err.msg)
		sys.exit(1)
	else:
		print_log(inputs, "[INFO:PASS] Duplicate Marking & Sorting complete.")
	finally:
		ngs_clean(inputs)

	## Merge Bam Files ##
	targets_data_processed = []
	try:
		ngs_preprocess.ngs_merge(inputs,targets_data,targets_data_processed)
	except ngs_classes.ngsExcept as err:
		print_log(inputs, err.msg)
		sys.exit(1)
	else:
		print_log(inputs, "[INFO:PASS] Bam merging complete.")
	finally:
		ngs_clean(inputs)

	## BQSR ##
	if(inputs.bqsr):
		try:
			ngs_preprocess.bqsr(inputs,targets_data_processed)
		except ngs_classes.ngsExcept as err:
			print_log(inputs, err.msg)
			sys.exit(1)
		else:
			print_log(inputs, "[INFO:PASS] BQSR complete.")
		finally:
			ngs_clean(inputs)
	else:
		for s in targets_data_processed:
			if(inputs.workflow == "GERMLINESNVINDEL" or inputs.workflow == "RNASEQ"):
				s._bamCalib = s._bam
			elif(inputs.workflow == "SOMATICSNVINDEL"):
				s._bamCalibTumor = s._bamTumor
				s._bamCalibNormal = s._bamNormal

	if inputs.verbose: print_log(inputs, "[INFO] Preprocessing done.")
	
	## Workflow ##
	if(inputs.workflow == "GERMLINESNVINDEL"):
		
		if inputs.verbose: print_log(inputs, "[INFO] GERMLINESNVINDEL Analysis")
		
		## SNV Indel ##
		if(re.search("HAPLOTYPECALLER",inputs.tool)):
			try:
				ngs_germline.ngs_haplotypecaller(inputs,targets_data_processed)
			except ngs_classes.ngsExcept as err:
				print_log(inputs, err.msg)
				sys.exit(1)
			finally:
				ngs_clean(inputs)

			try:
				ngs_germline.ngs_haplotypecaller_combine(inputs,targets_data_processed)
			except ngs_classes.ngsExcept as err:
				print_log(inputs, err.msg)
				sys.exit(1)
			finally:
				ngs_clean(inputs)

			try:
				ngs_germline.ngs_haplotypecaller_genotype(inputs)
			except ngs_classes.ngsExcept as err:
				print_log(inputs, err.msg)
				sys.exit(1)
			finally:
				ngs_clean(inputs)

			try:
				ngs_germline.ngs_haplotypecaller_filter(inputs)
			except ngs_classes.ngsExcept as err:
				print_log(inputs, err.msg)
				sys.exit(1)
			finally:
				ngs_clean(inputs)

		elif(re.search("STRELKA2", inputs.tool)):
			try:
				ngs_germline.ngs_strelka_germline(inputs,targets_data_processed)
			except ngs_classes.ngsExcept as err:
				print_log(inputs, err.msg)
				sys.exit(1)
			finally:
				ngs_clean(inputs)

		## CNV ##
		if inputs.verbose: 
			try:
				ngs_germline.ngs_cnvkit_germline(inputs,targets_data_processed)
			except ngs_classes.ngsExcept as err:
				print_log(inputs, err.msg)
				sys.exit(1)
			finally:
				ngs_clean(inputs)

	elif(inputs.workflow == "SOMATICSNVINDEL"):
	
		if inputs.verbose: print_log(inputs, "[INFO] SOMATICSNVINDEL Analysis")

		## SNV Indel ##
		if(re.search("STRELKA2", inputs.tool)):
			try:
				ngs_somatic.ngs_strelka_somatic(inputs,targets_data_processed)
			except ngs_classes.ngsExcept as err:
				print_log(inputs, err.msg)
				sys.exit(1)
			finally:
				ngs_clean(inputs)

		elif(re.search("MUTECT2", inputs.tool)):
			try:
				ngs_somatic.ngs_mutect(inputs,targets_data_processed)
			except ngs_classes.ngsExcept as err:
				print_log(inputs, err.msg)
				sys.exit(1)
			finally:
				ngs_clean(inputs)
		
		## CNV ##
		try:
			ngs_somatic.ngs_cnvkit_somatic(inputs,targets_data_processed)
		except ngs_classes.ngsExcept as err:
			print_log(inputs, err.msg)
			sys.exit(1)
		finally:
			ngs_clean(inputs)
	
	elif inputs.workflow == "RNASEQ":
		
		if inputs.verbose: print_log(inputs, "[INFO] RNASEQ Analysis")
		try:
			ngs_rna.ngs_quantification(inputs)
		except ngs_classes.ngsExcept as err:
			print_log(inputs, err.msg)
			sys.exit(1)

	sys.exit(0)
