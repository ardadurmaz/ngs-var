######################################################
import sys, os, argparse
import os.path, re
import glob
import ngs_classes
######################################################

def ngs_haplotypecaller(inputs, targets_data_processed):
	
	print("\n<---> HaplotypeCaller <--->\n") ## 0 return code for success
	runDir = os.path.abspath(inputs.dir + '/HCaller/')

	if(os.path.exists(runDir)):
		if inputs.verbose: print("[INFO] HaplotypeCaller directory found")
	else:
		if inputs.verbose: print("[INFO] Creating HaplotypeCaller directory: %s" % (runDir))
		try:
			os.mkdir(runDir)
		except OSError:
			raise ngs_classes.ngsExcept("[ERROR] Failed to create HaplotypeCaller directory %s" % (runDir))
	
	## Index ##
	for s in targets_data_processed:
		if inputs.verbose: print('[INFO:COMMAND] %s index -t %d %s' % (inputs.sambamba, inputs.threads, s._bamCalib))
		if not inputs.dry: os.system('%s index -t %d %s' % (inputs.sambamba, inputs.threads, s._bamCalib))

	for s in targets_data_processed:
		outFile = runDir + f'/{s._id}.raw.snps.indels.vcf'
		inFile = s._bamCalib
		sampleName = s._id
		
		if(os.path.exists(outFile)):
			if inputs.verbose: print(f"[INFO] HaplotypeCaller raw file found for sample {s._id}, skipping.")
			return 0
		else:		
			if inputs.verbose:
				print("[INFO:COMMAND] java -Xmx8g -jar %s HaplotypeCaller --native-pair-hmm-threads %d -ERC GVCF -R %s -I %s "\
				"--sample-name %s "\
				"-stand-call-conf 10 -O %s "\
				"-A BaseQualityRankSumTest -A Coverage -A DepthPerAlleleBySample "\
				"-A FisherStrand -A QualByDepth -A ReadPosRankSumTest" % (inputs.gatk, inputs.threads, inputs.reference, inFile, sampleName, outFile))
			if not inputs.dry:
				ret_val = os.system("java -Xmx8g -jar %s HaplotypeCaller --native-pair-hmm-threads %d -ERC GVCF -R %s -I %s "\
				"--sample-name %s "\
				"-stand-call-conf 10 -O %s "\
				"-A BaseQualityRankSumTest -A Coverage -A DepthPerAlleleBySample "\
				"-A FisherStrand -A QualByDepth -A ReadPosRankSumTest" % (inputs.gatk, inputs.threads, inputs.reference, inFile, sampleName, outFile))
				if(ret_val >> 8 != 0): raise ngs_classes.ngsExcept(f"[ERROR:HaplotypeCaller] Failed to call germline variants for sample {s._id}")
	
	print("\n<---> Done <--->\n")
	return 0

def ngs_haplotypecaller_combine(inputs, targets_data_processed):
	print("\n<---> CombineGVCFs <--->\n")
	runDir = os.path.abspath(inputs.dir + '/HCaller/')
	
	files = [runDir + f'/{s._id}.raw.snps.indels.vcf' for s in targets_data_processed]
	inFile = " --variant ".join(files)
	outFile = runDir + '/cohort.g.vcf'

	if(os.path.exists(outFile)):
		if inputs.verbose: print("[INFO] Skipping")
		return 0
	else:
		if inputs.verbose:
			print("[INFO:COMMAND] java -Xmx8g -jar %s CombineGVCFs -R %s --variant %s -O %s" % \
		       (inputs.gatk, inputs.reference, inFile, outFile))
		if not inputs.dry:
			ret_val = os.system("java -Xmx8g -jar %s CombineGVCFs -R %s --variant %s -O %s" % \
		       (inputs.gatk, inputs.reference, inFile, outFile))
			if(ret_val >> 8 != 0): raise ngs_classes.ngsExcept("[ERROR:CombineGVCFs] Failed")

	print("\n<---> Done <--->\n")
	return 0

def ngs_haplotypecaller_genotype(inputs):
    print("\n<---> GenotypeGVCFs <--->\n")
    runDir = os.path.abspath(inputs.dir + '/HCaller/')

    inFile = runDir + '/cohort.g.vcf'
    outFile = runDir + '/cohort_jointcall.vcf'

    if(os.path.exists(outFile)):
        if inputs.verbose: print("[INFO] Skipping")
        return 0
    else:
        if inputs.verbose:
            print("java -Xmx8g -jar %s GenotypeGVCFs -R %s -V %s -O %s -A BaseQualityRankSumTest -A MappingQualityRankSumTest -A FisherStrand -A QualByDepth -A ReadPosRankSumTest -A RMSMappingQuality -A StrandOddsRatio -A DepthPerAlleleBySample -A Coverage -ip 100" % \
                (inputs.gatk, inputs.reference, inFile, outFile))
        if not inputs.dry:
            ret_val = os.system("java -Xmx8g -jar %s GenotypeGVCFs -R %s -V %s -O %s -A BaseQualityRankSumTest -A MappingQualityRankSumTest -A FisherStrand -A QualByDepth -A ReadPosRankSumTest -A RMSMappingQuality -A StrandOddsRatio -A DepthPerAlleleBySample -A Coverage -ip 100" % \
                (inputs.gatk, inputs.reference, inFile, outFile))
            if(ret_val >> 8 != 0): raise ngs_classes.ngsExcept("[ERROR:GenotypeGVCFs] Failed")
    print("\n<---> Done <--->\n")
    return 0

def ngs_haplotypecaller_filter(inputs):
    print("\n<---> Variant Filtration <--->\n")
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
        if inputs.verbose:
            print("[INFO:COMMAND] %s" % command) 
        if not inputs.dry:
            ret_val = os.system(command)
            if(ret_val >> 8 != 0): 
                raise ngs_classes.ngsExcept("[ERROR:VariantFiltration] Failed to process command: %s" % command)
        
    print("\n<---> Done <--->\n")
    return 0

def ngs_strelka_germline(inputs, targets_data_processed):
	
	print("\n<---> Strelka2 Germline Workflow <--->\n")
	
	runDir = os.path.abspath(inputs.dir + '/STRELKA_GERMLINE/')
	runWork = inputs.strelkabin + '/configureStrelkaGermlineWorkflow.py'
	
	if(os.path.exists(runDir)):
		if inputs.verbose: print("[INFO] Strelka germline run directory found, skipping.")
	else:
		try:
			os.mkdir(runDir)
		except OSError:
			raise ngs_classes.ngsExcept("[ERROR] Failed to create HaplotypeCaller directory %s" % (runDir))

	## Index Bam Files ##
	for s in targets_data_processed:
		if inputs.verbose: print('[INFO:COMMAND] %s index -t %d %s' % (inputs.sambamba, inputs.threads, s._bamCalib))
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
					print("[INFO:COMMAND] python %s --referenceFasta %s --exome --runDir %s --callRegions %s --bam %s" % (runWork,inputs.reference,runDir,inputs.cbed,inFile))
				else:
					print("[INFO:COMMAND] python %s --referenceFasta %s --exome --runDir %s --bam %s" % (runWork,inputs.reference,runDir,inFile))
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
					print("[INFO:COMMAND] python %s --referenceFasta %s --runDir %s --callRegions %s --bam %s" % (runWork,inputs.reference,runDir,inputs.cbed,inFile))
				else:
					print("[INFO:COMMAND] python %s --referenceFasta %s --runDir %s --bam %s" % (runWork,inputs.reference,runDir,inFile))
				
			if not inputs.dry:
				if not(inputs.cbed is None):
					ret_val = os.system("python %s --referenceFasta %s --runDir %s --callRegions %s --bam %s" % (runWork,inputs.reference,runDir,inputs.cbed,inFile))
					if(ret_val >> 8 != 0): raise ngs_classes.ngsExcept("[ERROR:STRELKA] Failed to configure strelka2 run")
				else:
					ret_val = os.system("python %s --referenceFasta %s --runDir %s --bam %s" % (runWork,inputs.reference,runDir,inFile))
					if(ret_val >> 8 != 0): raise ngs_classes.ngsExcept("[ERROR:STRELKA] Failed to configure strelka2 run")
	
	
		## Run ##
		if inputs.verbose:
			print("[INFO:COMMAND] python %s -m local -j %d" % (runDir + '/runWorkflow.py',inputs.threads))
		if not inputs.dry:
			ret_val = os.system("python %s -m local -j %d" % (runDir + '/runWorkflow.py', inputs.threads))
			if(ret_val >> 8 != 0): raise ngs_classes.ngsExcept("[ERROR:STRELKA] Failed to run strelka2 configuration")
		
	print("\n<---> Done <--->\n")

	return 0


def ngs_cnvkit_germline(inputs, targets_data_processed):
	
	print("\n<---> Cnvkit Germline Workflow <--->\n")
	runDir = inputs.dir + '/CNVKIT_GERMLINE/'
	if(os.path.exists(runDir)):
		print("[INFO] Cnvkit germline directory found, skipping.")
	else:
		try:
			os.mkdir(runDir)
		except OSError:
			raise ngs_classes.ngsExcept("[ERROR] Failed to create directory %s" % (runDir))

	## Index Bam Files ##
	for s in targets_data_processed:
		if inputs.verbose: print('[INFO:COMMAND] %s index -t %d %s' % (inputs.sambamba, inputs.threads, s._bamCalib))
		if not inputs.dry: os.system('%s index -t %d %s' % (inputs.sambamba, inputs.threads, s._bamCalib))
	
	tumor_files = []
	for s in targets_data_processed:
		tumor_files.append(s._bamCalib)
	
	if inputs.exome:
		if inputs.verbose:
			if not(inputs.bed is None):
				print("[INFO:COMMAND] %s batch --drop-low-coverage -p %d --normal -f %s -t %s -d %s %s" % (inputs.cnvkit,inputs.threads,inputs.reference,inputs.bed,runDir," ".join(tumor_files)))
			else:
				print("[ERROR:COMMAND] Cannot run copy number analysis without target file for exome data")
				return False
		if not inputs.dry:
			if not(inputs.bed is None):
				os.system("%s batch --drop-low-coverage -p %d --normal  -f %s -t %s -d %s %s" % (inputs.cnvkit,inputs.threads,inputs.reference,inputs.bed,runDir," ".join(tumor_files)))
			else:
				print("[ERROR:COMMAND] Cannot run copy number analysis without target file for exome data")
				return False
	else:
		if inputs.verbose:
			if not(inputs.bed is None):
				print("[INFO:COMMAND] %s batch -m wgs --drop-low-coverage -p %d --normal -f %s -t %s -d %s %s" % (inputs.cnvkit,inputs.threads,inputs.reference,inputs.bed,runDir," ".join(tumor_files)))
			else:
				print("[INFO:COMMAND] %s batch -m wgs --drop-low-coverage -p %d --normal -f %s -d %s %s" % (inputs.cnvkit,inputs.threads,inputs.reference,runDir," ".join(tumor_files)))
		if not inputs.dry:
			if not(inputs.bed is None):
				os.system("%s batch -m wgs --drop-low-coverage -p %d --normal -f %s -t %s -d %s %s" % (inputs.cnvkit,inputs.threads,inputs.reference,inputs.bed,runDir," ".join(tumor_files)))
			else:
				os.system("%s batch -m wgs --drop-low-coverage -p %d --normal -f %s -d %s %s" % (inputs.cnvkit,inputs.threads,inputs.reference,runDir," ".join(tumor_files)))
				
	print("\n<---> Done. <--->\n")
	return 0
