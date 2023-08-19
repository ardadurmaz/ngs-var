######################################################
import sys, os, argparse
import os.path, re
import glob
import ngs_classes
from ngs_functions import print_log, run_command
######################################################

def ngs_haplotypecaller(inputs, targets_data_processed):
	
	print_log(inputs, "\n<---> HaplotypeCaller <--->\n") ## 0 return code for success
	runDir = os.path.abspath(inputs.dir + '/HCaller/')

	if(os.path.exists(runDir)):
		if inputs.verbose: print_log(inputs, "[INFO] HaplotypeCaller directory found")
	else:
		if inputs.verbose: print_log(inputs, "[INFO] Creating HaplotypeCaller directory: %s" % (runDir))
		try:
			os.mkdir(runDir)
		except OSError:
			raise ngs_classes.ngsExcept("[ERROR] Failed to create HaplotypeCaller directory %s" % (runDir))
	
	## Index ##
	for s in targets_data_processed:
		if(len(glob.glob(os.path.abspath(s._bamCalib + '.bai'))) > 0):
			print_log(inputs, "[INFO] Indices found, skipping.")
		else:
			run_command(inputs, '%s index -t %d %s' % (inputs.sambamba, inputs.threads, s._bamCalib))

	for s in targets_data_processed:
		outFile = runDir + f'/{s._id}.raw.snps.indels.vcf'
		inFile = s._bamCalib
		sampleName = s._id
		
		if(os.path.exists(outFile)):
			if inputs.verbose: print_log(inputs, f"[INFO] HaplotypeCaller raw file found for sample {s._id}, skipping.")
			return 0
		else:		
			run_command(inputs, "java -Xmx8g -jar %s HaplotypeCaller --native-pair-hmm-threads %d -ERC GVCF -R %s -I %s "\
				"--sample-name %s "\
				"-stand-call-conf 10 -O %s "\
				"-A BaseQualityRankSumTest -A Coverage -A DepthPerAlleleBySample "\
				"-A FisherStrand -A QualByDepth -A ReadPosRankSumTest" % (inputs.gatk, inputs.threads, inputs.reference, inFile, sampleName, outFile), f"[ERROR:HaplotypeCaller] Failed to call germline variants for sample {s._id}")
	
	print_log(inputs, "\n<---> Done <--->\n")
	return 0

def ngs_haplotypecaller_combine(inputs, targets_data_processed):
	print_log(inputs, "\n<---> CombineGVCFs <--->\n")
	runDir = os.path.abspath(inputs.dir + '/HCaller/')
	
	files = [runDir + f'/{s._id}.raw.snps.indels.vcf' for s in targets_data_processed]
	inFile = " --variant ".join(files)
	outFile = runDir + '/cohort.g.vcf'

	if(os.path.exists(outFile)):
		if inputs.verbose: print_log(inputs, "[INFO] Skipping")
		return 0
	else:
		run_command(inputs, "java -Xmx8g -jar %s CombineGVCFs -R %s --variant %s -O %s" % \
		       (inputs.gatk, inputs.reference, inFile, outFile), "[ERROR:CombineGVCFs] Failed to combine GVCFs")

	print_log(inputs, "\n<---> Done <--->\n")
	return 0

def ngs_haplotypecaller_genotype(inputs):
    print_log(inputs, "\n<---> GenotypeGVCFs <--->\n")
    runDir = os.path.abspath(inputs.dir + '/HCaller/')

    inFile = runDir + '/cohort.g.vcf'
    outFile = runDir + '/cohort_jointcall.vcf'

    if(os.path.exists(outFile)):
        if inputs.verbose: print_log(inputs, "[INFO] Skipping")
        return 0
    else:
        run_command(inputs, "java -Xmx8g -jar %s GenotypeGVCFs -R %s -V %s -O %s -A BaseQualityRankSumTest -A MappingQualityRankSumTest -A FisherStrand -A QualByDepth -A ReadPosRankSumTest -A RMSMappingQuality -A StrandOddsRatio -A DepthPerAlleleBySample -A Coverage -ip 100" % \
                (inputs.gatk, inputs.reference, inFile, outFile), "[ERROR:GenotypeGVCFs] Failed to perform joint genotyping")
    print_log(inputs, "\n<---> Done <--->\n")
    return 0

def ngs_haplotypecaller_filter(inputs):
    print_log(inputs, "\n<---> Variant Filtration <--->\n")
    runDir = os.path.abspath(inputs.dir + '/HCaller/')
    inFile = runDir + '/cohort_jointcall.vcf'
    normFile = runDir + '/norm.cohort.jointcall.vcf'
    snpFile = runDir + '/snp.norm.cohort.jointcall.vcf'
    indelFile = runDir + '/indel.norm.cohort.jointcall.vcf'
    filtSnpFile = runDir + '/snp.filt.norm.cohort.jointcall.vcf'
    filtIndelFile = runDir + '/indel.filt.norm.cohort.jointcall.vcf'
    mergedFile = runDir + '/filt.norm.cohort.jointcall.vcf'
    mainOutFile = runDir + '/filt.gatk.calls.vcf'
    
    commands = [
        "java -Xmx8g -jar %s LeftAlignAndTrimVariants -R %s -V %s -O %s --split-multi-allelics" % (inputs.gatk, inputs.reference, inFile, normFile),
        "java -Xmx8g -jar %s SelectVariants -R %s -V %s -O %s --exclude-non-variants --select-type-to-include \"SNP\"" % (inputs.gatk, inputs.reference, normFile, snpFile),
        "java -Xmx8g -jar %s SelectVariants -R %s -V %s -O %s --exclude-non-variants --select-type-to-include \"INDEL\"" % (inputs.gatk, inputs.reference, normFile, indelFile),
        "java -Xmx8g -jar %s VariantFiltration -R %s -V %s -O %s --filter-expression \"QD < 2.0\" --filter-name \"LowQD\" --filter-expression \"QUAL < 30.0\" --filter-name \"QUAL30\" --filter-expression 'FS > 60.0 && SOR > 9.0' --filter-name 'HighSBSNP' --filter-expression 'MQ < 40.0' --filter-name \"LowMQSNP\" --filter-expression 'ReadPosRankSum < -8.0' --filter-name \"LowRPRSSNP\" --genotype-filter-expression \"GQ < 30.0\" --genotype-filter-name \"GQ\" --genotype-filter-expression \"DP < 10.0\" --genotype-filter-name \"LowDPGT\"" % (inputs.gatk, inputs.reference, snpFile, filtSnpFile),
        "java -Xmx8g -jar %s VariantFiltration -R %s -V %s -O %s --filter-expression \"QD < 2.0\" --filter-name \"LowQD\" --filter-expression \"QUAL < 30.0\" --filter-name \"QUAL30\" --filter-expression 'FS > 200.0 && SOR > 9.0' --filter-name \"HighSBINDEL\" --filter-expression 'ReadPosRankSum < -20.0' --filter-name \"LowRPRSINDEL\" --genotype-filter-expression \"GQ < 30.0\" --genotype-filter-name \"GQ\" --genotype-filter-expression \"DP < 10.0\" --genotype-filter-name \"LowDPGT\"" % (inputs.gatk, inputs.reference, indelFile, filtIndelFile),
        "java -Xmx8g -jar %s MergeVcfs -I %s -I %s -O %s" % (inputs.gatk, filtSnpFile, filtIndelFile, mergedFile),
        "java -Xmx8g -jar %s SelectVariants -R %s -V %s -O %s --exclude-non-variants --exclude-filtered --select-type-to-include \"SNP\" --select-type-to-include \"INDEL\"" % (inputs.gatk, inputs.reference, mergedFile, mainOutFile)
    ]

    for command in commands:
        run_command(inputs, command, "[ERROR:VariantFiltration] Failed to process command: %s" % command)
        
    print_log(inputs, "\n<---> Done <--->\n")
    return 0

def ngs_strelka_germline(inputs, targets_data_processed):
	
	print_log(inputs, "\n<---> Strelka2 Germline Workflow <--->\n")
	
	runDir = os.path.abspath(inputs.dir + '/STRELKA_GERMLINE/')
	runWork = inputs.strelkabin + '/configureStrelkaGermlineWorkflow.py'
	
	if(os.path.exists(runDir)):
		if inputs.verbose: print_log(inputs, "[INFO] Strelka germline run directory found, skipping.")
	else:
		try:
			os.mkdir(runDir)
		except OSError:
			raise ngs_classes.ngsExcept("[ERROR] Failed to create HaplotypeCaller directory %s" % (runDir))

	## Index Bam Files ##
	for s in targets_data_processed:
		if(len(glob.glob(os.path.abspath(s._bamCalib + '.bai'))) > 0):
			print_log(inputs, "[INFO] Indices found, skipping.")
		else:
			run_command(inputs, '%s index -t %d %s' % (inputs.sambamba, inputs.threads, s._bamCalib))
	
	files = []
	for s in targets_data_processed:
		files.append(s._bamCalib)
	inFile = " --bam ".join(files)
	
	## Configure ##
	if os.path.exists(runDir + '/results'):
		if inputs.verbose:
			print("[INFO] Strelka germline results found, skipping.")
		return 0
	
	base_command = f"python2.7 {runWork} --referenceFasta {inputs.reference} --runDir {runDir} --bam {inFile}"
	if inputs.exome:
		base_command += " --exome"
	if inputs.cbed is not None:
		base_command += f" --callRegions {inputs.cbed}"

	error_msg = "[ERROR:STRELKA] Failed to configure strelka2 run"
	run_command(inputs, base_command, error_msg)

	## Run ##
	run_workflow_command = f"python2.7 {runDir}/runWorkflow.py -m local -j {inputs.threads}"
	error_msg_run = "[ERROR:STRELKA] Failed to run strelka2 configuration"
	run_command(inputs, run_workflow_command, error_msg_run)

		
	print_log(inputs, "\n<---> Done <--->\n")

	return 0


def ngs_cnvkit_germline(inputs, targets_data_processed):
	
	print_log(inputs, "\n<---> Cnvkit Germline Workflow <--->\n")
	runDir = inputs.dir + '/CNVKIT_GERMLINE/'
	if(os.path.exists(runDir)):
		print_log(inputs, "[INFO] Cnvkit germline directory found, skipping.")
	else:
		try:
			os.mkdir(runDir)
		except OSError:
			raise ngs_classes.ngsExcept("[ERROR] Failed to create directory %s" % (runDir))

	## Index Bam Files ##
	for s in targets_data_processed:
		if(len(glob.glob(os.path.abspath(s._bamCalib + '.bai'))) > 0):
			print_log(inputs, "[INFO] Indices found, skipping.")
		else:
			run_command(inputs, '%s index -t %d %s' % (inputs.sambamba, inputs.threads, s._bamCalib))
	
	files = []
	for s in targets_data_processed:
		files.append(s._bamCalib)
	
	outBed = runDir + 'access.bed'
	bedPathParts = inputs.bed.rsplit('/', 1)
	dir, fileName = bedPathParts
	localBed = dir + '/' + 'local.' + fileName
	targets = localBed.split('.bed')[0] + '.target.bed'
	antitargets = localBed.split('.bed')[0] + '.antitarget.bed'
	targetCoverage = runDir + s._id + ".Normal" + ".targetcoverage.cnn"
	antitargetCoverage = runDir + s._id + ".Normal" + ".antitargetcoverage.cnn"
	outRef = runDir + 'flatReference.cnn'
	cnr = runDir + s._id + ".Normal" + ".cnr"
	cns = runDir + s._id + ".Normal" + ".cns"
	error_msg = "[ERROR:CNVKIT] Failed to run CNVkit configuration"

	run_command(inputs, "%s target %s --split -o %s" % (inputs.cnvkit, inputs.bed, localBed), error_msg)
	run_command(inputs, "%s access %s -o %s" % (inputs.cnvkit, inputs.reference, outBed), error_msg)

	if inputs.exome:
		run_command(inputs, "%s autobin %s -t %s -g %s" % (inputs.cnvkit, " ".join(files), localBed, outBed), error_msg)
	else:
		run_command(inputs, "%s autobin %s -t %s -g %s -m wgs" % (inputs.cnvkit, " ".join(files), localBed, outBed), error_msg)

	for s in targets_data_processed:
		run_command(inputs, "%s coverage %s %s -o %s -p %d" % (inputs.cnvkit, s._bamCalib, targets, targetCoverage, inputs.threads), error_msg)
		run_command(inputs, "%s coverage %s %s -o %s -p %d" % (inputs.cnvkit, s._bamCalib, antitargets, antitargetCoverage, inputs.threads), error_msg)
		run_command(inputs, "%s reference -o %s -f %s -t %s -a %s" % (inputs.cnvkit, outRef, inputs.reference, targets, antitargets), error_msg)
		run_command(inputs, "%s fix %s %s %s -o %s" % (inputs.cnvkit, targetCoverage, antitargetCoverage, outRef, cnr), error_msg)
		run_command(inputs, "%s segment %s -o %s -m cbs --rscript-path %s --smooth-cbs --drop-low-coverage --drop-outliers 10" % (inputs.cnvkit, cnr, cns, inputs.rscript), error_msg)

	print_log(inputs, "\n<---> Done. <--->\n")
	return 0
