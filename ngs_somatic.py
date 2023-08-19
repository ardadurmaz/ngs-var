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
			raise ngs_classes.ngsExcept("[ERROR] Failed to create directory" % (runDir))

	for s in targets_data_processed:
		outFile = runDir + '/Sample_' + s._id + '.vcf'

		if(os.path.exists(outFile)):
			if inputs.verbose: print_log(inputs, "[INFO] MuTect2 output file found for sample %s, skipping." % (s._id))
		else:
			if not(inputs.bed is None):
				run_command(inputs, "java -jar %s -T MuTect2 -R %s -nct %d -I:tumor %s -I:normal %s -L %s -o %s" % (inputs.gatk,inputs.reference,inputs.threads,
																														s._bamCalibTumor,s._bamCalibNormal,
																														inputs.bed,outFile), "[ERROR] Failed to run MuTect2 on files %s %s" % (s._bamCalibTumor, s._bamCalibNormal))
			else:
				run_command(inputs, "java -jar %s -T MuTect2 -R %s -nct %d -I:tumor %s -I:normal %s -o %s" % (inputs.gatk,inputs.reference,inputs.threads,
																												  s._bamCalibTumor,s._bamCalibNormal,outFile), "[ERROR] Failed to run MuTect2 on files %s %s" % (s._bamCalibTumor,s._bamCalibNormal))
				
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
			run_command(inputs, '%s index -t %d %s' % (inputs.sambamba, inputs.threads, s._bamCalibTumor))
		if(len(glob.glob(os.path.abspath(s._bamCalibNormal + '.bai'))) > 0):
			print_log(inputs, "[INFO] Indices found, skipping.")
		else:
			run_command(inputs, '%s index -t %d %s' % (inputs.sambamba, inputs.threads, s._bamCalibNormal))
			
	if(os.path.exists(runDir)):
		if inputs.verbose: print_log(inputs, "[INFO] Strelka run directory exists")
	else:
		try:
			os.mkdir(runDir)
		except OSError:
			raise ngs_classes.ngsExcept("[ERROR] Failed to create directory %s" % (runDir))
	## Manta ##
	for s in targets_data_processed:
		
		## Config ##
		s._mantaWD = runDir + '/MANTA_WD_Sample_' + s._id
		mantaIndel = s._mantaWD + '/results/variants/candidateSmallIndels.vcf.gz'

		if(os.path.exists(s._mantaWD)):
			print_log(inputs, "[INFO] Manta run directory found %s" %(s._mantaWD))
			if(os.path.exists(mantaIndel)):
				print_log(inputs, "[INFO] Manta results found for sample %s skipping" % (s._id))
				continue
		else:
			try:
				os.mkdir(s._mantaWD)
			except OSError:
				raise ngs_classes.ngsExcept("[ERROR] Failed to create directory %s" % (s._mantaWD))

		base_command = "python2.7 %s --normalBam %s --tumorBam %s --referenceFasta %s --runDir %s" % (mantaRun,
                                                                                                  s._bamCalibNormal,
                                                                                                  s._bamCalibTumor,
                                                                                                  inputs.reference,
                                                                                                  s._mantaWD)
		if inputs.cbed is not None:
			base_command += " --callRegions %s" % inputs.cbed

		if inputs.exome:
			base_command += " --exome"

		error_msg = "[ERROR] Failed to run manta workflow for files %s %s" % (s._bamCalibNormal, s._bamCalibTumor)
		run_command(inputs, base_command, error_msg)

		run_workflow_command = "python2.7 %s -m local -j %d" % (os.path.abspath(s._mantaWD + '/runWorkflow.py'), inputs.threads)
		error_msg_workflow = "[ERROR] Failed to execute Manta workflow script"
		run_command(inputs, run_workflow_command, error_msg_workflow)
			
	## Strelka ##
	for s in targets_data_processed:
		
		## Config ##
		s._strelkaWD = runDir + '/STRELKA_WD_Sample_' + s._id
		mantaIndel = s._mantaWD + '/results/variants/candidateSmallIndels.vcf.gz'
		
		if (os.path.exists(s._strelkaWD)):
			print_log(inputs, "[INFO] Strelka run directory found %s" % (s._strelkaWD))
			if(os.path.exists(s._strelkaWD + '/results/variants/somatic.snvs.vcf.gz')):
				print_log(inputs, "[INFO] Strelka results found for sample %s skipping." % (s._id))
				continue
		else:
			try:
				os.mkdir(s._strelkaWD)
			except OSError:
				raise ngs_classes.ngsExcept("[ERROR] Failed to create directory %s" % (s._strelkaWD))
		
		base_command = "python2.7 %s --normalBam %s --tumorBam %s --referenceFasta %s --indelCandidates %s --runDir %s" % (strelkaRun,
																														s._bamCalibNormal,
																														s._bamCalibTumor,
																														inputs.reference,
																														mantaIndel,
																														s._strelkaWD)

		if inputs.exome:
			base_command += " --exome"
		
		if inputs.cbed is not None:
			base_command += " --callRegions %s" % inputs.cbed

		error_msg = "[ERROR] Failed to run strelka workflow for files %s %s" % (s._bamCalibNormal, s._bamCalibTumor)
		run_command(inputs, base_command, error_msg)

		run_workflow_command = "python2.7 %s -m local -j %d" % (os.path.abspath(s._strelkaWD + '/runWorkflow.py'), inputs.threads)
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
			raise ngs_classes.ngsExcept("[ERROR] Failed to create directory %s" % (runDir))
	
	## Index ##
	for s in targets_data_processed:
		if(len(glob.glob(os.path.abspath(s._bamCalibTumor + '.bai'))) > 0):
			print_log(inputs, "[INFO] Indices found, skipping.")
		else:
			run_command(inputs, '%s index -t %d %s' % (inputs.sambamba, inputs.threads, s._bamCalibTumor))
		if(len(glob.glob(os.path.abspath(s._bamCalibNormal + '.bai'))) > 0):
			print_log(inputs, "[INFO] Indices found, skipping.")
		else:
			run_command(inputs, '%s index -t %d %s' % (inputs.sambamba, inputs.threads, s._bamCalibNormal))
	
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

	run_command(inputs, "%s target %s --split -o %s" % (inputs.cnvkit, inputs.bed, localBed), error_msg)
	run_command(inputs, "%s access %s -o %s" % (inputs.cnvkit, inputs.reference, outBed), error_msg)
	for s in targets_data_processed:

		if inputs.exome:
			run_command(inputs, "%s autobin %s -t %s -g %s" % (inputs.cnvkit, " ".join(files), inputs.bed, outBed), error_msg)
		else:
			run_command(inputs, "%s autobin %s -t %s -g %s -m wgs" % (inputs.cnvkit, " ".join(files), inputs.bed, outBed), error_msg)

	for s in targets_data_processed:
		targetCoverageNormal = runDir + s._id + ".Normal" + ".targetcoverage.cnn"
		antitargetCoverageNormal = runDir + s._id + ".Normal" + ".antitargetcoverage.cnn"
		targetCoverageTumor = runDir + s._id + ".Tumor" + ".targetcoverage.cnn"
		antitargetCoverageTumor = runDir + s._id + ".Tumor" + ".antitargetcoverage.cnn"
		run_command(inputs, "%s coverage %s %s -o %s -p %d" % (inputs.cnvkit, s._bamCalibNormal, targets, targetCoverageNormal, inputs.threads), error_msg)
		run_command(inputs, "%s coverage %s %s -o %s -p %d" % (inputs.cnvkit, s._bamCalibNormal, antitargets, antitargetCoverageNormal, inputs.threads), error_msg)
		run_command(inputs, "%s coverage %s %s -o %s -p %d" % (inputs.cnvkit, s._bamCalibTumor, targets, targetCoverageTumor, inputs.threads), error_msg)
		run_command(inputs, "%s coverage %s %s -o %s -p %d" % (inputs.cnvkit, s._bamCalibTumor, antitargets, antitargetCoverageTumor, inputs.threads), error_msg)

		run_command(inputs, "%s reference %s --fasta %s -o %s" % (inputs.cnvkit, targetCoverageNormal + " " + antitargetCoverageNormal, inputs.reference, outRef), error_msg)
	
	for s in targets_data_processed:
		cnr = runDir + s._id + ".Tumor" + ".cnr"
		cns = runDir + s._id + ".Tumor" + ".cns"
		targetCoverageTumor = runDir + s._id + ".Tumor" + ".targetcoverage.cnn"
		antitargetCoverageTumor = runDir + s._id + ".Tumor" + ".antitargetcoverage.cnn"
		run_command(inputs, "%s fix %s %s %s -o %s" % (inputs.cnvkit, targetCoverageTumor, antitargetCoverageTumor, outRef, cnr), error_msg)
		run_command(inputs, "%s segment %s -o %s -m cbs --rscript-path %s --smooth-cbs --drop-low-coverage --drop-outliers 10" % (inputs.cnvkit, cnr, cns, inputs.rscript), error_msg)

	print_log(inputs, "\n<---> Done. <--->\n")
	return 0
