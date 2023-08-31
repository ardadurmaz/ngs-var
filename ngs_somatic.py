######################################################
import sys, os, argparse
import os.path, re
import glob
import ngs_classes
from ngs_functions import print_log, run_command
######################################################


def ngs_mutect(inputs, targets_data_processed):

	print_log(inputs, "\n<---> MuTect 2 <--->\n")
	
	runDir = inputs.dir + '/MUTECT_SOMATIC'
	if(os.path.exists(runDir)):
		print_log(inputs, "[INFO] MuTect2 run directory exists, skipping.")
	else:
		try:
			os.mkdir(runDir)
		except OSError:
			raise ngs_classes.ngsExcept(f"[ERROR] Failed to create directory {runDir}")

	for s in targets_data_processed:
		outFile = runDir + '/Sample_' + s._id + '.vcf'

		if(os.path.exists(outFile)):
			if inputs.verbose: print_log(inputs, f"[INFO] MuTect2 output file found for sample {s._id}, skipping.")
		else:
			if inputs.bed:
				run_command(inputs, f"java -jar {inputs.gatk }-T MuTect2 -R {inputs.reference} -nct {inputs.threads} -I:tumor {s._bamCalibTumor} -I:normal {s._bamCalibNormal} -L {inputs.bed} -o {outFile}", f"[ERROR] Failed to run MuTect2 on files {s._bamCalibTumor} {s._bamCalibNormal}")
			else:
				run_command(inputs, f"java -jar {inputs.gatk} -T MuTect2 -R {inputs.reference} -nct {inputs.threads} -I:tumor {s._bamCalibTumor} -I:normal {s._bamCalibNormal} -o {outFile}", f"[ERROR] Failed to run MuTect2 on files {s._bamCalibTumor} {s._bamCalibNormal}")
				
	print_log(inputs, "\n<---> Done <--->\n")
	return 0


def ngs_strelka_somatic(inputs, targets_data_processed):
	
	print_log(inputs, "\n<---> Strelka2 Somatic Workflow <--->\n")
	
	runDir = inputs.dir + '/STRELKA_SOMATIC'
	mantaRun = os.path.abspath(inputs.mantabin + '/configManta.py')
	strelkaRun = os.path.abspath(inputs.strelkabin + '/configureStrelkaSomaticWorkflow.py')
	
	## Index ##
	for s in targets_data_processed:
		if(len(glob.glob(os.path.abspath(s._bamCalibTumor + '.bai'))) > 0):
			print_log(inputs, "[INFO] Indices found, skipping.")
		else:
			run_command(inputs, f'{inputs.sambamba} index -t {inputs.threads} {s._bamCalibTumor}')
		if(len(glob.glob(os.path.abspath(s._bamCalibNormal + '.bai'))) > 0):
			print_log(inputs, "[INFO] Indices found, skipping.")
		else:
			run_command(inputs, f'{inputs.sambamba} index -t {inputs.threads} {s._bamCalibNormal}')
			
	if(os.path.exists(runDir)):
		if inputs.verbose: print_log(inputs, "[INFO] Strelka run directory exists")
	else:
		try:
			os.mkdir(runDir)
		except OSError:
			raise ngs_classes.ngsExcept(f"[ERROR] Failed to create directory {runDir}")
	## Manta ##
	for s in targets_data_processed:
		
		## Config ##
		s._mantaWD = runDir + '/MANTA_WD_Sample_' + s._id
		mantaIndel = s._mantaWD + '/results/variants/candidateSmallIndels.vcf.gz'

		if(os.path.exists(s._mantaWD)):
			print_log(inputs, f"[INFO] Manta run directory found {s._mantaWD}")
			if(os.path.exists(mantaIndel)):
				print_log(inputs, f"[INFO] Manta results found for sample {s._id} skipping")
				continue
		else:
			try:
				os.mkdir(s._mantaWD)
			except OSError:
				raise ngs_classes.ngsExcept(f"[ERROR] Failed to create directory {s._mantaWD}")

		base_command = f"python2.7 {mantaRun} --normalBam {s._bamCalibNormal} --tumorBam {s._bamCalibTumor} --referenceFasta {inputs.reference} --runDir {s._mantaWD}"
		if inputs.cbed is not None:
			base_command += f" --callRegions {inputs.cbed}"

		if inputs.exome:
			base_command += " --exome"

		error_msg = f"[ERROR] Failed to run manta workflow for files {s._bamCalibNormal} {s._bamCalibTumor}"
		run_command(inputs, base_command, error_msg)

		mantaPath = os.path.abspath(s._mantaWD + '/runWorkflow.py')

		run_workflow_command = f"python2.7 {mantaPath} -m local -j {inputs.threads}"
		error_msg_workflow = "[ERROR] Failed to execute Manta workflow script"
		run_command(inputs, run_workflow_command, error_msg_workflow)
			
	## Strelka ##
	for s in targets_data_processed:
		
		## Config ##
		s._strelkaWD = runDir + '/STRELKA_WD_Sample_' + s._id
		mantaIndel = s._mantaWD + '/results/variants/candidateSmallIndels.vcf.gz'
		
		if (os.path.exists(s._strelkaWD)):
			print_log(inputs, f"[INFO] Strelka run directory found {s._strelkaWD}")
			if(os.path.exists(s._strelkaWD + '/results/variants/somatic.snvs.vcf.gz')):
				print_log(inputs, f"[INFO] Strelka results found for sample {s._id} skipping.")
				continue
		else:
			try:
				os.mkdir(s._strelkaWD)
			except OSError:
				raise ngs_classes.ngsExcept(f"[ERROR] Failed to create directory {s._strelkaWD}")
		
		base_command = f"python2.7 {strelkaRun} --normalBam {s._bamCalibNormal} --tumorBam {s._bamCalibTumor} --referenceFasta {inputs.reference} --indelCandidates {mantaIndel} --runDir {s._strelkaWD}"

		if inputs.exome:
			base_command += " --exome"
		
		if inputs.cbed is not None:
			base_command += f" --callRegions {inputs.cbed}"

		error_msg = f"[ERROR] Failed to run strelka workflow for files {s._bamCalibNormal} {s._bamCalibTumor}"
		run_command(inputs, base_command, error_msg)

		strelkaPath = os.path.abspath(s._strelkaWD + '/runWorkflow.py')
		run_workflow_command = f"python2.7 {strelkaPath} -m local -j {inputs.threads}"
		error_msg_workflow = "[ERROR] Failed to execute Strelka workflow script"
		run_command(inputs, run_workflow_command, error_msg_workflow)
			
			
	print_log(inputs, "\n<---> Done <--->\n")
	return 0


def ngs_cnvkit_somatic(inputs, targets_data_processed):
	
	print_log(inputs, "\n<---> Cnvkit Somatic Workflow <--->\n")
	runDir = inputs.dir + '/CNVKIT_SOMATIC/'
	if(os.path.exists(runDir)):
		if inputs.verbose: print_log(inputs, "[INFO] Cnvkit somatic directory found")
	else:
		try:
			os.mkdir(runDir)
		except OSError:
			raise ngs_classes.ngsExcept(f"[ERROR] Failed to create directory {runDir}")
	
	## Index ##
	for s in targets_data_processed:
		if(len(glob.glob(os.path.abspath(s._bamCalibTumor + '.bai'))) > 0):
			print_log(inputs, "[INFO] Indices found, skipping.")
		else:
			run_command(inputs, f'{inputs.sambamba} index -t {inputs.threads} {s._bamCalibTumor}')
		if(len(glob.glob(os.path.abspath(s._bamCalibNormal + '.bai'))) > 0):
			print_log(inputs, "[INFO] Indices found, skipping.")
		else:
			run_command(inputs, f'{inputs.sambamba} index -t {inputs.threads} {s._bamCalibNormal}')
	
	outBed = runDir + 'access.bed'
	bedPathParts = inputs.bed.rsplit('/', 1)
	dir, fileName = bedPathParts
	localBed = dir + '/' + 'local.' + fileName
	targets = localBed.split('.bed')[0] + '.target.bed'
	antitargets = localBed.split('.bed')[0] + '.antitarget.bed'
	outRef = runDir + 'reference.cnn'
	error_msg = "[ERROR:CNVKIT] Failed to run CNVkit configuration"
	files = []
	for s in targets_data_processed:
		files.append(s._bamCalibNormal)
		files.append(s._bamCalibTumor)

	run_command(inputs, f"{inputs.cnvkit} target {inputs.bed} --split -o {localBed}", error_msg)
	run_command(inputs, f"{inputs.cnvkit} access {inputs.reference} -o {outBed}", error_msg)
	for s in targets_data_processed:

		if inputs.exome:
			run_command(inputs, f"{inputs.cnvkit} autobin {' '.join(files)} -t {inputs.bed} -g {outBed}", error_msg)
		else:
			run_command(inputs, f"{inputs.cnvkit} autobin {' '.join(files)} -t {inputs.bed} -g {outBed} -m wgs", error_msg)

	for s in targets_data_processed:
		targetCoverageNormal = runDir + s._id + ".Normal" + ".targetcoverage.cnn"
		antitargetCoverageNormal = runDir + s._id + ".Normal" + ".antitargetcoverage.cnn"
		targetCoverageTumor = runDir + s._id + ".Tumor" + ".targetcoverage.cnn"
		antitargetCoverageTumor = runDir + s._id + ".Tumor" + ".antitargetcoverage.cnn"
		run_command(inputs, f"{inputs.cnvkit} coverage {s._bamCalibNormal} {targets} -o {targetCoverageNormal} -p {inputs.threads}", error_msg)
		run_command(inputs, f"{inputs.cnvkit} coverage {s._bamCalibNormal} {antitargets} -o {antitargetCoverageNormal} -p {inputs.threads}", error_msg)
		run_command(inputs, f"{inputs.cnvkit} coverage {s._bamCalibTumor} {targets} -o {targetCoverageTumor} -p {inputs.threads}", error_msg)
		run_command(inputs, f"{inputs.cnvkit} coverage {s._bamCalibTumor} {antitargets} -o {antitargetCoverageTumor} -p {inputs.threads}", error_msg)
		run_command(inputs, f"{inputs.cnvkit} reference {targetCoverageNormal + ' ' + antitargetCoverageNormal} --fasta {inputs.reference} -o {outRef}", error_msg)
	
	for s in targets_data_processed:
		cnr = runDir + s._id + ".Tumor" + ".cnr"
		cns = runDir + s._id + ".Tumor" + ".cns"
		targetCoverageTumor = runDir + s._id + ".Tumor" + ".targetcoverage.cnn"
		antitargetCoverageTumor = runDir + s._id + ".Tumor" + ".antitargetcoverage.cnn"
		run_command(inputs, f"{inputs.cnvkit} fix {targetCoverageTumor} {antitargetCoverageTumor} {outRef} -o {cnr}", error_msg)
		run_command(inputs, f"{inputs.cnvkit} segment {cnr} -o {cns} -m cbs --rscript-path {inputs.rscript} --smooth-cbs --drop-low-coverage --drop-outliers 10", error_msg)

	print_log(inputs, "\n<---> Done. <--->\n")
	return 0
