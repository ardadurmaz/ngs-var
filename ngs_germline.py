######################################################
import sys, os, argparse
import os.path, re
import glob
import ngs_classes
######################################################


def ngs_haplotypecaller(inputs, targets_data_processed):
	
	print "\n<---> HaplotypeCaller <--->\n" ## 0 return code for success
	runDir = os.path.abspath(inputs.dir + '/HCaller/')

	if(os.path.exists(runDir)):
		if inputs.verbose: print "[INFO] HaplotypeCaller directory found"
	else:
		if inputs.verbose: print "[INFO] Creating HaplotypeCaller directory: %s" % (runDir)
		try:
			os.mkdir(runDir)
		except OSError:
			raise ngs_classes.ngsExcept("[ERROR] Failed to create HaplotypeCaller directory %s" % (runDir))
	
	## Index ##
	for s in targets_data_processed:
		if inputs.verbose: print '[INFO:COMMAND] %s index -t %d %s' % (inputs.sambamba, inputs.threads, s._bamCalib)
		if not inputs.dry: os.system('%s index -t %d %s' % (inputs.sambamba, inputs.threads, s._bamCalib))
	
	files = []
	for s in targets_data_processed:
		files.append(s._bamCalib)
	
	outFile = runDir + '/HCaller.raw.snps.indels.vcf'
	inFile = " -I ".join(files)
	
	if(os.path.exists(outFile)):
		if inputs.verbose: print "[INFO] HaplotypeCaller raw file found, skipping."
		return 0
	else:		
		if inputs.verbose:
			print "[INFO:COMMAND] java -jar %s -T HaplotypeCaller -R %s -I %s "\
			"--num_cpu_threads_per_data_thread %d -stand_call_conf 10 -o %s "\
			"-A BaseQualityRankSumTest -A Coverage -A DepthPerAlleleBySample "\
			"-A FisherStrand -A QualByDepth -A ReadPosRankSumTest" % (inputs.gatk,inputs.reference,inFile,inputs.threads,outFile)
		if not inputs.dry:
			ret_val = os.system("java -jar %s -T HaplotypeCaller -R %s -I %s "\
					  "--num_cpu_threads_per_data_thread %d -stand_call_conf 10 -o %s "\
					  "-A BaseQualityRankSumTest -A Coverage -A DepthPerAlleleBySample "\
					  "-A FisherStrand -A QualByDepth -A ReadPosRankSumTest" % (inputs.gatk,inputs.reference,inFile,inputs.threads,outFile))
			if(ret_val >> 8 != 0): raise ngs_classes.ngsExcept("[ERROR:HaplotypeCaller] Failed to call germline variants")
		
	print "\n<---> Done <--->\n"
	return 0


def ngs_haplotypecaller_filter(inputs):
	
	print "\n<---> Variant Filtration <--->\n"
	runDir = os.path.abspath(inputs.dir + '/HCaller/')
	inFile = runDir + '/HCaller.raw.snps.indels.vcf'
	outFile = runDir + '/HCaller.pre.filtered.snps.indels.vcf'
	mainOutFile = runDir + '/HCaller.filtered.snps.indels.vcf'
	
	if(os.path.exists(mainOutFile)):
		if inputs.verbose: print "[INFO] HaplotypeCaller filtered vcf file found, skipping."
	else:	
		if inputs.verbose:
			print "[INFO:COMMAND] java -jar %s -T VariantFiltration -R %s -o %s -V %s "\
			"--filterExpression \'QD < 2.0\' --filterName \"LowQD\" "\
			"--filterExpression \'FS > 60.0 && vc.getType() == \"SNP\"\' --filterName \"HighFSSNP\" "\
			"--filterExpression \'FS > 200.0 && vc.getType() == \"INDEL\"\' --filterName \"HighFSINDEL\" "\
			"--filterExpression \'MQ < 40.0 && vc.getType() == \"SNP\"\' --filterName \"LowMQSNP\" "\
			"--filterExpression \'ReadPosRankSum < -8.0 && vc.getType() == \"SNP\"\' --filterName \"LowRPRSSNP\" "\
			"--filterExpression \'ReadPosRankSum < -20.0 && vc.getType() == \"INDEL\"\' --filterName \"LowRPRSINDEL\" "\
			"--genotypeFilterExpression \"DP < 10\" --genotypeFilterName \"LowDPGT\" --setFilteredGtToNocall" % (inputs.gatk,inputs.reference,outFile,inFile)
		
		if not inputs.dry:
			ret_val = os.system("java -jar %s -T VariantFiltration -R %s -o %s -V %s "\
					  "--filterExpression \'QD < 2.0\' --filterName \"LowQD\" "\
					  "--filterExpression \'FS > 60.0 && vc.getType() == \"SNP\"\' --filterName \"HighFSSNP\" "\
					  "--filterExpression \'FS > 200.0 && vc.getType() == \"INDEL\"\' --filterName \"HighFSINDEL\" "\
					  "--filterExpression \'MQ < 40.0 && vc.getType() == \"SNP\"\' --filterName \"LowMQSNP\" "\
					  "--filterExpression \'ReadPosRankSum < -8.0 && vc.getType() == \"SNP\"\' --filterName \"LowRPRSSNP\" "\
					  "--filterExpression \'ReadPosRankSum < -20.0 && vc.getType() == \"INDEL\"\' --filterName \"LowRPRSINDEL\" "\
					  "--genotypeFilterExpression \"DP < 10\" --genotypeFilterName \"LowDPGT\" --setFilteredGtToNocall" % (inputs.gatk,inputs.reference,outFile,inFile))
			if(ret_val >> 8 != 0): raise ngs_classes.ngsExcept("[ERROR:VariantFiltration] Failed to process raw vcf file")
		
		if inputs.verbose:
			print "[INFO:COMMAND] java -jar %s -T SelectVariants -R %s "\
			"--excludeNonVariants --excludeFiltered -selectType SNP -selectType INDEL "\
			"-o %s -V %s" % (inputs.gatk,inputs.reference,mainOutFile,outFile)
		
		if not inputs.dry:
			ret_val = os.system("java -jar %s -T SelectVariants -R %s "\
								"--excludeNonVariants --excludeFiltered -selectType SNP -selectType INDEL "\
								"-o %s -V %s" % (inputs.gatk,inputs.reference,mainOutFile,outFile))
			if(ret_val >> 8 != 0): raise ngs_classes.ngsExcept("[ERROR:VariantFiltration] Failed to process pre filtered vcf file")
		
	print "\n<---> Done <--->\n"
	return 0


def ngs_strelka_germline(inputs, targets_data_processed):
	
	print "\n<---> Strelka2 Germline Workflow <--->\n"
	
	runDir = os.path.abspath(inputs.dir + '/STRELKA_GERMLINE/')
	runWork = inputs.strelkabin + '/configureStrelkaGermlineWorkflow.py'
	
	if(os.path.exists(runDir)):
		if inputs.verbose: print "[INFO] Strelka germline run directory found, skipping."
	else:
		try:
			os.mkdir(runDir)
		except OSError:
			raise ngs_classes.ngsExcept("[ERROR] Failed to create HaplotypeCaller directory %s" % (runDir))

	## Index Bam Files ##
	for s in targets_data_processed:
		if inputs.verbose: print '[INFO:COMMAND] %s index -t %d %s' % (inputs.sambamba, inputs.threads, s._bamCalib)
		if not inputs.dry: os.system('%s index -t %d %s' % (inputs.sambamba, inputs.threads, s._bamCalib))
		
	
	
	
	files = []
	for s in targets_data_processed:
		files.append(s._bamCalib)
	inFile = " --bam ".join(files)
	
	## Configure ##
	if(os.path.exists(runDir + '/results')):
		if inputs.verbose: "[INFO] Strelka germline results found, skipping."
		return 0
	else:
		if inputs.exome:
			if inputs.verbose:
				if not(inputs.cbed is None):
					print "[INFO:COMMAND] python %s --referenceFasta %s --exome --runDir %s --callRegions %s --bam %s" % (runWork,inputs.reference,runDir,inputs.cbed,inFile)
				else:
					print "[INFO:COMMAND] python %s --referenceFasta %s --exome --runDir %s --bam %s" % (runWork,inputs.reference,runDir,inFile)
			if not inputs.dry:
				if not(inputs.cbed is None):
					ret_val = os.system("python %s --referenceFasta %s --exome --runDir %s --callRegions %s --bam %s" % (runWork,inputs.reference,runDir,inputs.cbed,inFile))
					if(ret_val >> 8 != 0): raise ngs_classes.ngsExcept("[ERROR:STRELKA] Failed to configure strelka2 run")
				else:
					ret_val = os.system("python %s --referenceFasta %s --exome --runDir %s --bam %s" % (runWork,inputs.reference,runDir,inFile))
					if(ret_val >> 8 != 0): raise ngs_classes.ngsExcept("[ERROR:STRELKA] Failed to configure strelka2 run")
				
		else:
			if inputs.verbose:
				if not(inputs.cbed is None):
					print "[INFO:COMMAND] python %s --referenceFasta %s --runDir %s --callRegions %s --bam %s" % (runWork,inputs.reference,runDir,inputs.cbed,inFile)
				else:
					print "[INFO:COMMAND] python %s --referenceFasta %s --runDir %s --bam %s" % (runWork,inputs.reference,runDir,inFile)
				
			if not inputs.dry:
				if not(inputs.cbed is None):
					ret_val = os.system("python %s --referenceFasta %s --runDir %s --callRegions %s --bam %s" % (runWork,inputs.reference,runDir,inputs.cbed,inFile))
					if(ret_val >> 8 != 0): raise ngs_classes.ngsExcept("[ERROR:STRELKA] Failed to configure strelka2 run")
				else:
					ret_val = os.system("python %s --referenceFasta %s --runDir %s --bam %s" % (runWork,inputs.reference,runDir,inFile))
					if(ret_val >> 8 != 0): raise ngs_classes.ngsExcept("[ERROR:STRELKA] Failed to configure strelka2 run")
	
	
		## Run ##
		if inputs.verbose:
			print "[INFO:COMMAND] python %s -m local -j %d" % (runDir + '/runWorkflow.py',inputs.threads)
		if not inputs.dry:
			ret_val = os.system("python %s -m local -j %d" % (runDir + '/runWorkflow.py', inputs.threads))
			if(ret_val >> 8 != 0): raise ngs_classes.ngsExcept("[ERROR:STRELKA] Failed to run strelka2 configuration")
		
	print "\n<---> Done <--->\n"

	return 0


def ngs_cnvkit_germline(inputs, targets_data_processed):
	
	print "\n<---> Cnvkit Germline Workflow <--->\n"
	runDir = inputs.dir + '/CNVKIT_GERMLINE/'
	if(os.path.exists(runDir)):
		print "[INFO] Cnvkit germline directory found, skipping."
	else:
		try:
			os.mkdir(runDir)
		except OSError:
			raise ngs_classes.ngsExcept("[ERROR] Failed to create directory %s" % (runDir))

	## Index Bam Files ##
	for s in targets_data_processed:
		if inputs.verbose: print '[INFO:COMMAND] %s index -t %d %s' % (inputs.sambamba, inputs.threads, s._bamCalib)
		if not inputs.dry: os.system('%s index -t %d %s' % (inputs.sambamba, inputs.threads, s._bamCalib))
	
	tumor_files = []
	for s in targets_data_processed:
		tumor_files.append(s._bamCalib)
	
	if inputs.exome:
		if inputs.verbose:
			if not(inputs.bed is None):
				print "[INFO:COMMAND] %s batch --drop-low-coverage -p %d --normal -f %s -t %s -d %s %s" % (inputs.cnvkit,inputs.threads,inputs.reference,inputs.bed,runDir," ".join(tumor_files))
			else:
				print "[ERROR:COMMAND] Cannot run copy number analysis without target file for exome data"
				return False
		if not inputs.dry:
			if not(inputs.bed is None):
				os.system("%s batch --drop-low-coverage -p %d --normal  -f %s -t %s -d %s %s" % (inputs.cnvkit,inputs.threads,inputs.reference,inputs.bed,runDir," ".join(tumor_files)))
			else:
				print "[ERROR:COMMAND] Cannot run copy number analysis without target file for exome data"
				return False
	else:
		if inputs.verbose:
			if not(inputs.bed is None):
				print "[INFO:COMMAND] %s batch -m wgs --drop-low-coverage -p %d --normal -f %s -t %s -d %s %s" % (inputs.cnvkit,inputs.threads,inputs.reference,inputs.bed,runDir," ".join(tumor_files))
			else:
				print "[INFO:COMMAND] %s batch -m wgs --drop-low-coverage -p %d --normal -f %s -d %s %s" % (inputs.cnvkit,inputs.threads,inputs.reference,runDir," ".join(tumor_files))
		if not inputs.dry:
			if not(inputs.bed is None):
				os.system("%s batch -m wgs --drop-low-coverage -p %d --normal -f %s -t %s -d %s %s" % (inputs.cnvkit,inputs.threads,inputs.reference,inputs.bed,runDir," ".join(tumor_files)))
			else:
				os.system("%s batch -m wgs --drop-low-coverage -p %d --normal -f %s -d %s %s" % (inputs.cnvkit,inputs.threads,inputs.reference,runDir," ".join(tumor_files)))
				
	print "\n<---> Done. <--->\n"
	return 0
