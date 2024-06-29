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
		if inputs.verbose: print_log(inputs, f"[INFO] Creating HaplotypeCaller directory: {runDir}")
		try:
			os.mkdir(runDir)
		except OSError:
			raise ngs_classes.ngsExcept(f"[ERROR] Failed to create HaplotypeCaller directory {runDir}")
	
	## Index ##
	for s in targets_data_processed:
		if(len(glob.glob(os.path.abspath(s._bamCalib + '.bai'))) > 0):
			print_log(inputs, "[INFO] Indices found, skipping.")
		else:
			run_command(inputs, f"{inputs.sambamba} index -t {inputs.threads} {s._bamCalib}")

	for s in targets_data_processed:
		outFile = runDir + f'/{s._id}.raw.snps.indels.vcf'
		inFile = s._bamCalib
		sampleName = s._id
		
		if(os.path.exists(outFile)):
			if inputs.verbose: print_log(inputs, f"[INFO] HaplotypeCaller raw file found for sample {s._id}, skipping.")
			return 0
		else:		
			run_command(inputs, f"{inputs.gatk} HaplotypeCaller --native-pair-hmm-threads {inputs.threads} -ERC GVCF -R {inputs.reference} -I {inFile} "\
				f"--sample-name {sampleName} "\
				f"-stand-call-conf 10 -O {outFile} "\
				"-A BaseQualityRankSumTest -A Coverage -A DepthPerAlleleBySample "\
				"-A FisherStrand -A QualByDepth -A ReadPosRankSumTest", f"[ERROR:HaplotypeCaller] Failed to call germline variants for sample {s._id}")
	
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
		run_command(inputs, f"{inputs.gatk} CombineGVCFs -R {inputs.reference} --variant {inFile} -O {outFile}", "[ERROR:CombineGVCFs] Failed to combine GVCFs")

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
        run_command(inputs, f"{inputs.gatk} GenotypeGVCFs -R {inputs.reference} -V {inFile} -O {outFile} -A BaseQualityRankSumTest -A MappingQualityRankSumTest -A FisherStrand -A QualByDepth -A ReadPosRankSumTest -A RMSMappingQuality -A StrandOddsRatio -A DepthPerAlleleBySample -A Coverage -ip 100", "[ERROR:GenotypeGVCFs] Failed to perform joint genotyping")
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
        f"{inputs.gatk} LeftAlignAndTrimVariants -R {inputs.reference} -V {inFile} -O {normFile} --split-multi-allelics",
        f"{inputs.gatk} SelectVariants -R {inputs.reference} -V {normFile} -O {snpFile} --exclude-non-variants --select-type-to-include \"SNP\"",
        f"{inputs.gatk} SelectVariants -R {inputs.reference} -V {normFile} -O {indelFile} --exclude-non-variants --select-type-to-include \"INDEL\"",
        f"{inputs.gatk} VariantFiltration -R {inputs.reference} -V {snpFile} -O {filtSnpFile} --filter-expression \"QD < 2.0\" --filter-name \"LowQD\" --filter-expression \"QUAL < 30.0\" --filter-name \"QUAL30\" --filter-expression 'FS > 60.0 && SOR > 9.0' --filter-name 'HighSBSNP' --filter-expression 'MQ < 40.0' --filter-name \"LowMQSNP\" --filter-expression 'ReadPosRankSum < -8.0' --filter-name \"LowRPRSSNP\" --genotype-filter-expression \"GQ < 30.0\" --genotype-filter-name \"GQ\" --genotype-filter-expression \"DP < 10.0\" --genotype-filter-name \"LowDPGT\"",
        f"{inputs.gatk} VariantFiltration -R {inputs.reference} -V {indelFile} -O {filtIndelFile} --filter-expression \"QD < 2.0\" --filter-name \"LowQD\" --filter-expression \"QUAL < 30.0\" --filter-name \"QUAL30\" --filter-expression 'FS > 200.0 && SOR > 9.0' --filter-name \"HighSBINDEL\" --filter-expression 'ReadPosRankSum < -20.0' --filter-name \"LowRPRSINDEL\" --genotype-filter-expression \"GQ < 30.0\" --genotype-filter-name \"GQ\" --genotype-filter-expression \"DP < 10.0\" --genotype-filter-name \"LowDPGT\"",
        f"{inputs.gatk} MergeVcfs -I {filtSnpFile} -I {filtIndelFile} -O {mergedFile}",
        f"{inputs.gatk} SelectVariants -R {inputs.reference} -V {mergedFile} -O {mainOutFile} --exclude-non-variants --exclude-filtered --select-type-to-include \"SNP\" --select-type-to-include \"INDEL\""
    ]

    for command in commands:
        run_command(inputs, command, f"[ERROR:VariantFiltration] Failed to process command: {command}")
        
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
			raise ngs_classes.ngsExcept(f"[ERROR] Failed to create HaplotypeCaller directory {runDir}")

	## Index Bam Files ##
	for s in targets_data_processed:
		if(len(glob.glob(os.path.abspath(s._bamCalib + '.bai'))) > 0):
			print_log(inputs, "[INFO] Indices found, skipping.")
		else:
			run_command(inputs, f'{inputs.sambamba} index -t {inputs.threads} {s._bamCalib}')
	
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
			raise ngs_classes.ngsExcept(f"[ERROR] Failed to create directory {runDir}")

	## Index Bam Files ##
	for s in targets_data_processed:
		if(len(glob.glob(os.path.abspath(s._bamCalib + '.bai'))) > 0):
			print_log(inputs, "[INFO] Indices found, skipping.")
		else:
			run_command(inputs, f'{inputs.sambamba} index -t {inputs.threads} {s._bamCalib}')
	
	files = []
	for s in targets_data_processed:
		files.append(s._bamCalib)
	
	outBed = runDir + 'access.bed'
	bedPathParts = inputs.bed.rsplit('/', 1)
	dir, fileName = bedPathParts
	localBed = dir + '/' + 'local.' + fileName
	targets = localBed.split('.bed')[0] + '.target.bed'
	antitargets = localBed.split('.bed')[0] + '.antitarget.bed'
	outRef = runDir + 'flatReference.cnn'
	error_msg = "[ERROR:CNVKIT] Failed to run CNVkit configuration"

	run_command(inputs, f"{inputs.cnvkit} target {inputs.bed} --split -o {localBed}", error_msg)
	run_command(inputs, f"{inputs.cnvkit} access {inputs.reference} -o {outBed}", error_msg)

	working_dir = os.getcwd()
	output_dir = os.path.dirname(localBed)

	os.chdir(output_dir)
	if inputs.exome:
		run_command(inputs, f"{inputs.cnvkit} autobin {' '.join(files)} -t {localBed} -g {outBed}", error_msg)
	else:
		run_command(inputs, f"{inputs.cnvkit} autobin {' '.join(files)} -t {localBed} -g {outBed} -m wgs", error_msg)
	os.chdir(working_dir)

	for s in targets_data_processed:
		targetCoverage = runDir + s._id + ".Normal" + ".targetcoverage.cnn"
		antitargetCoverage = runDir + s._id + ".Normal" + ".antitargetcoverage.cnn"
		cnr = runDir + s._id + ".Normal" + ".cnr"
		cns = runDir + s._id + ".Normal" + ".cns"
		run_command(inputs, f"{inputs.cnvkit} coverage {s._bamCalib} {targets} -o {targetCoverage} -p {inputs.threads}", error_msg)
		run_command(inputs, f"{inputs.cnvkit} coverage {s._bamCalib} {antitargets} -o {antitargetCoverage} -p {inputs.threads}", error_msg)
		run_command(inputs, f"{inputs.cnvkit} reference -o {outRef} -f {inputs.reference} -t {targets} -a {antitargets}", error_msg)
		run_command(inputs, f"{inputs.cnvkit} fix {targetCoverage} {antitargetCoverage} {outRef} -o {cnr}", error_msg)
		run_command(inputs, f"{inputs.cnvkit} segment {cnr} -o {cns} -m cbs --rscript-path {inputs.rscript} --smooth-cbs --drop-low-coverage --drop-outliers 10", error_msg)

	print_log(inputs, "\n<---> Done. <--->\n")
	return 0
