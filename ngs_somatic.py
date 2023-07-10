######################################################
import sys, os, argparse
import os.path, re
import glob
import ngs_classes
######################################################


def ngs_mutect(inputs, targets_data_processed):

	print("\n<---> MuTect 2 <--->\n")
	
	runDir = inputs.dir + '/MUTECT_SOMATIC'
	if(os.path.exists(runDir)):
		print("[INFO] MuTect2 run directory exists, skipping.")
	else:
		try:
			os.mkdir(runDir)
		except OSError:
			raise ngs_classes.ngsExcept("[ERROR] Failed to create directory" % (runDir))

	for s in targets_data_processed:
		outFile = runDir + '/Sample_' + s._id + '.vcf'
		tumorFile = runDir + '/TempTumor_' + s._id + '.bam'
		normalFile = runDir + '/TempNormal_' + s._id + '.bam'
		
		if(os.path.exists(tumorFile) and os.path.exists(normalFile)):
			if inputs.verbose: print("[INFO] Temporary tumor and normal bam files found for sample %s, skipping." % (s._id))
		else:
			if inputs.verbose:
				print("[INFO:COMMAND] java -jar %s AddOrReplaceReadGroups I=%s O=%s RGID=%s RGLB=%s RGPL=%s RGPU=%s RGSM=%s" % (inputs.gatk,s._bamCalibTumor,tumorFile,s._id + '_Tumor_Lib1','Lib1','ILLUMINA',s._id + '_Tumor_Lib1',s._id + '_Tumor'))
			if not inputs.dry:
				ret_val = os.system("java -jar %s AddOrReplaceReadGroups I=%s O=%s RGID=%s RGLB=%s RGPL=%s RGPU=%s RGSM=%s" % (inputs.gatk,s._bamCalibTumor,tumorFile,s._id + '_Tumor_Lib1','Lib1','ILLUMINA',s._id + '_Tumor_Lib1',s._id + '_Tumor'))
			if inputs.verbose:
				print("[INFO:COMMAND] java -jar %s AddOrReplaceReadGroups I=%s O=%s RGID=%s RGLB=%s RGPL=%s RGPU=%s RGSM=%s" % (inputs.gatk,s._bamCalibNormal,normalFile,s._id + '_Normal_Lib2','Lib2','ILLUMINA',s._id + '_Normal_Lib2',s._id + '_Normal'))
			if not inputs.dry:
				ret_val = os.system("java -jar %s AddOrReplaceReadGroups I=%s O=%s RGID=%s RGLB=%s RGPL=%s RGPU=%s RGSM=%s" % (inputs.gatk,s._bamCalibNormal,normalFile,s._id + '_Normal_Lib2','Lib2','ILLUMINA',s._id + '_Normal_Lib2',s._id + '_Normal'))
	
			if inputs.verbose: print('[INFO:COMMAND] %s index -t %d %s' % (inputs.sambamba, inputs.threads, tumorFile))
			if not inputs.dry: os.system('%s index -t %d %s' % (inputs.sambamba, inputs.threads, tumorFile))
			if inputs.verbose: print('[INFO:COMMAND] %s index -t %d %s' % (inputs.sambamba, inputs.threads, normalFile))
			if not inputs.dry: os.system('%s index -t %d %s' % (inputs.sambamba, inputs.threads, normalFile))
		
		if(os.path.exists(outFile)):
			if inputs.verbose: print("[INFO] MuTect2 output file found for sample %s, skipping." % (s._id))
		else:
			if inputs.verbose:
				if not(inputs.bed is None):
					print("[INFO:COMMAND] java -jar %s -T MuTect2 -R %s -nct %d -I:tumor %s -I:normal %s -L %s -o %s" % (inputs.gatk,inputs.reference,inputs.threads,
																														 tumorFile,normalFile,
																														 inputs.bed,outFile))
				else:
					print("[INFO:COMMAND] java -jar %s -T MuTect2 -R %s -nct %d -I:tumor %s -I:normal %s -o %s" % (inputs.gatk,inputs.reference,inputs.threads,
																												   tumorFile,normalFile,outFile))
			if not(inputs.dry):
				if not(inputs.bed is None):
					ret_val = os.system("java -jar %s -T MuTect2 -R %s -nct %d -I:tumor %s -I:normal %s -L %s -o %s" % (inputs.gatk,inputs.reference,inputs.threads,
																														tumorFile,normalFile,
																														inputs.bed,outFile))
					if(ret_val >> 8 != 0): raise ngs_classes.ngsExcept("[ERROR] Failed to run MuTect2 on files %s %s" % (tumorFile, normalFile))
				else:
					ret_val = os.system("java -jar %s -T MuTect2 -R %s -nct %d -I:tumor %s -I:normal %s -o %s" % (inputs.gatk,inputs.reference,inputs.threads,
																												  tumorFile,normalFile,outFile))
					if(ret_val >> 8 != 0): raise ngs_classes.ngsExcept("[ERROR] Failed to run MuTect2 on files %s %s" % (tumorFile,normalFile))
				
	print("\n<---> Done <--->\n")
	return 0




def ngs_strelka_somatic(inputs, targets_data_processed):
	
	print("\n<---> Strelka2 Somatic Workflow <--->\n")
	
	runDir = inputs.dir + '/STRELKA_SOMATIC'
	mantaRun = os.path.abspath(inputs.mantabin + '/configManta.py')
	strelkaRun = os.path.abspath(inputs.strelkabin + '/configureStrelkaSomaticWorkflow.py')
	
	## Index ##
	for s in targets_data_processed:
		if inputs.verbose: print('[INFO:COMMAND] %s index -t %d %s' % (inputs.sambamba, inputs.threads, s._bamCalibTumor))
		if not inputs.dry: os.system('%s index -t %d %s' % (inputs.sambamba, inputs.threads, s._bamCalibTumor))
		if inputs.verbose: print('[INFO:COMMAND] %s index -t %d %s' % (inputs.sambamba, inputs.threads, s._bamCalibNormal))
		if not inputs.dry: os.system('%s index -t %d %s' % (inputs.sambamba, inputs.threads, s._bamCalibNormal))
			
			
	if(os.path.exists(runDir)):
		if inputs.verbose: print("[INFO] Strelka run directory exists")
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
		tumorFile = s._mantaWD + '/TempTumor_' + s._id + '.bam'
		normalFile = s._mantaWD + '/TempNormal_' + s._id + '.bam'

		if(os.path.exists(s._mantaWD)):
			print("[INFO] Manta run directory found %s" %(s._mantaWD))
			if(os.path.exists(mantaIndel)):
				print("[INFO] Manta results found for sample %s skipping" % (s._id))
				continue
		else:
			try:
				os.mkdir(s._mantaWD)
			except OSError:
				raise ngs_classes.ngsExcept("[ERROR] Failed to create directory %s" % (s._mantaWD))

		if(os.path.exists(tumorFile) and os.path.exists(normalFile)):
			if inputs.verbose: print("[INFO] Temporary tumor and normal bam files found for sample %s, skipping." % (s._id))
		else:
			if inputs.verbose:
				print("[INFO:COMMAND] java -jar %s AddOrReplaceReadGroups I=%s O=%s RGID=%s RGLB=%s RGPL=%s RGPU=%s RGSM=%s" % (inputs.gatk,s._bamCalibTumor,tumorFile,s._id + '_Tumor_Lib1','Lib1','ILLUMINA',s._id + '_Tumor_Lib1',s._id + '_Tumor'))
			if not inputs.dry:
				ret_val = os.system("java -jar %s AddOrReplaceReadGroups I=%s O=%s RGID=%s RGLB=%s RGPL=%s RGPU=%s RGSM=%s" % (inputs.gatk,s._bamCalibTumor,tumorFile,s._id + '_Tumor_Lib1','Lib1','ILLUMINA',s._id + '_Tumor_Lib1',s._id + '_Tumor'))
			if inputs.verbose:
				print("[INFO:COMMAND] java -jar %s AddOrReplaceReadGroups I=%s O=%s RGID=%s RGLB=%s RGPL=%s RGPU=%s RGSM=%s" % (inputs.gatk,s._bamCalibNormal,normalFile,s._id + '_Normal_Lib2','Lib2','ILLUMINA',s._id + '_Normal_Lib2',s._id + '_Normal'))
			if not inputs.dry:
				ret_val = os.system("java -jar %s AddOrReplaceReadGroups I=%s O=%s RGID=%s RGLB=%s RGPL=%s RGPU=%s RGSM=%s" % (inputs.gatk,s._bamCalibNormal,normalFile,s._id + '_Normal_Lib2','Lib2','ILLUMINA',s._id + '_Normal_Lib2',s._id + '_Normal'))
	
			if inputs.verbose: print('[INFO:COMMAND] %s index -t %d %s' % (inputs.sambamba, inputs.threads, tumorFile))
			if not inputs.dry: os.system('%s index -t %d %s' % (inputs.sambamba, inputs.threads, tumorFile))
			if inputs.verbose: print('[INFO:COMMAND] %s index -t %d %s' % (inputs.sambamba, inputs.threads, normalFile))
			if not inputs.dry: os.system('%s index -t %d %s' % (inputs.sambamba, inputs.threads, normalFile))
		
		if inputs.exome:
			if inputs.verbose:
				if not(inputs.cbed is None):
					print("[INFO:COMMAND] python %s --normalBam %s --tumorBam %s --exome --referenceFasta %s --runDir %s --callRegions %s" % (mantaRun,
																																			  normalFile,
																																			  tumorFile,
																																			  inputs.reference,
																																			  s._mantaWD,
																																			  inputs.cbed))
				else:
					print("[INFO:COMMAND] python %s --normalBam %s --tumorBam %s --exome --referenceFasta %s --runDir %s" % (mantaRun,
																															 normalFile,
																															 tumorFile,
																															 inputs.reference,
																															 s._mantaWD))

			if not inputs.dry:
				if not(inputs.cbed is None):
					ret_val = os.system("python %s --normalBam %s --tumorBam %s --exome --referenceFasta %s --runDir %s --callRegions %s" % (mantaRun,
																																			 normalFile,
																																			 tumorFile,
																																			 inputs.reference,
																																			 s._mantaWD,
																																			 inputs.cbed))
					if(ret_val >> 8 != 0): raise ngs_classes.ngsExcept("[ERROR] Failed to run manta workflow for files %s %s" % (normalFile, tumorFile))
				else:
					ret_val = os.system("python %s --normalBam %s --tumorBam %s --exome --referenceFasta %s --runDir %s" % (mantaRun,
																															normalFile,
																															tumorFile,
																															inputs.reference,
																															s._mantaWD))
					if(ret_val >> 8 != 0): raise ngs_classes.ngsExcept("[ERROR] Failed to run manta workflow for files %s %s" % (normalFile, tumorFile))
		else:
			if inputs.verbose:
				if not(inputs.cbed is None):
					print("[INFO:COMMAND] python %s --normalBam %s --tumorBam %s --referenceFasta %s --runDir %s --callRegions %s" % (mantaRun,
																																	  normalFile,
																																	  tumorFile,
																																	  inputs.reference,
																																	  s._mantaWD,
																																	  inputs.cbed))
				else:
					print("[INFO:COMMAND] python %s --normalBam %s --tumorBam %s --referenceFasta %s --runDir %s" % (mantaRun,
																													 normalFile,
																													 tumorFile,
																													 inputs.reference,
																													 s._mantaWD))

			if not inputs.dry:
				if not(inputs.cbed is None):
					ret_val = os.system("python %s --normalBam %s --tumorBam %s --referenceFasta %s --runDir %s --callRegions %s" % (mantaRun,
																																	 normalFile,
																																	 tumorFile,
																																	 inputs.reference,
																																	 s._mantaWD,
																																	 inputs.cbed))
					if(ret_val >> 8 != 0): raise ngs_classes.ngsExcept("[ERROR] Failed to run manta workflow for files %s %s" % (normalFile,tumorFile))
				else:
					ret_val = os.system("python %s --normalBam %s --tumorBam %s --referenceFasta %s --runDir %s" % (mantaRun,
																													normalFile,
																													tumorFile,
																													inputs.reference,
																													s._mantaWD))
					if(ret_val >> 8 != 0): raise ngs_classes.ngsExcept("[ERROR] Failed to run manta workflow for files %s %s" % (normalFile,tumorFile))
		## Run ##			
		if inputs.verbose:
			print("[INFO:COMMAND] python %s -m local -j %d" % (os.path.abspath(s._mantaWD + '/runWorkflow.py'),inputs.threads))
			
		if not inputs.dry:
			ret_val = os.system("python %s -m local -j %d" % (os.path.abspath(s._mantaWD + '/runWorkflow.py'),inputs.threads))
			if(ret_val >> 8 != 0): raise ngs_classes.ngsExcept("[ERROR] Failed to execute Manta workflow script")
			
	
	## Strelka ##
	for s in targets_data_processed:
		
		## Config ##
		s._strelkaWD = runDir + '/STRELKA_WD_Sample_' + s._id
		mantaIndel = s._mantaWD + '/results/variants/candidateSmallIndels.vcf.gz'
		tumorFile = s._mantaWD + '/TempTumor_' + s._id + '.bam'
		normalFile = s._mantaWD + '/TempNormal_' + s._id + '.bam'
		
		if (os.path.exists(s._strelkaWD)):
			print("[INFO] Strelka run directory found %s" % (s._strelkaWD))
			if(os.path.exists(s._strelkaWD + '/results/variants/somatic.snvs.vcf.gz')):
				print("[INFO] Strelka results found for sample %s skipping." % (s._id))
				continue
		else:
			try:
				os.mkdir(s._strelkaWD)
			except OSError:
				raise ngs_classes.ngsExcept("[ERROR] Failed to create directory %s" % (s._strelkaWD))
		
		if inputs.exome:
			if inputs.verbose:
				if not(inputs.cbed is None):
					print("[INFO:COMMAND] python %s --normalBam %s --tumorBam %s --referenceFasta %s --indelCandidates %s --exome --callRegions %s --runDir %s" % (strelkaRun,
																																								   normalFile,
																																								   tumorFile,
																																								   inputs.reference,
																																								   mantaIndel,
																																								   inputs.cbed,
																																								   s._strelkaWD))
				else:
					print("[INFO:COMMAND] python %s --normalBam %s --tumorBam %s --referenceFasta %s --indelCandidates %s --exome --runDir %s" % (strelkaRun,
																																				  normalFile,
																																				  tumorFile,
																																				  inputs.reference,
																																				  mantaIndel,
																																				  s._strelkaWD))
			if not inputs.dry:
				if not(inputs.cbed is None):
					ret_val = os.system("python %s --normalBam %s --tumorBam %s --referenceFasta %s --indelCandidates %s --exome --callRegions %s --runDir %s" % (strelkaRun,
																																								  normalFile,
																																								  tumorFile,
																																								  inputs.reference,
																																								  mantaIndel,
																																								  inputs.cbed,
																																								  s._strelkaWD))
					if(ret_val >> 8 != 0): raise ngs_classes.ngsExcept("[ERROR] Failed to run strelka workflow for files %s %s" % (normalFile,tumorFile))
				else:
					ret_val = os.system("python %s --normalBam %s --tumorBam %s --referenceFasta %s --indelCandidates %s --exome --runDir %s" % (strelkaRun,
																																				 normalFile,
																																				 tumorFile,
																																				 inputs.reference,
																																				 mantaIndel,
																																				 s._strelkaWD))
					if(ret_val >> 8 != 0): raise ngs_classes.ngsExcept("[ERROR] Failed to run strelka workflow for files %s %s" % (normalFile,tumorFile))
		else:
			if inputs.verbose:
				if not(inputs.cbed is None):
					print("[INFO:COMMAND] python %s --normalBam %s --tumorBam %s --referenceFasta %s --indelCandidates %s --callRegions %s --runDir %s" % (strelkaRun,
																																						   normalFile,
																																						   tumorFile,
																																						   inputs.reference,
																																						   mantaIndel,
																																						   inputs.cbed,
																																						   s._strelkaWD))
				else:
					print("[INFO:COMMAND] python %s --normalBam %s --tumorBam %s --referenceFasta %s --indelCandidates %s --runDir %s" % (strelkaRun,
																																		  normalFile,
																																		  tumorFile,
																																		  inputs.reference,
																																		  mantaIndel,
																																		  s._strelkaWD))
			if not inputs.dry:
				if not(inputs.cbed is None):
					ret_val = os.system("python %s --normalBam %s --tumorBam %s --referenceFasta %s --indelCandidates %s --callRegions %s --runDir %s" % (strelkaRun,
																																						  normalFile,
																																						  tumorFile,
																																						  inputs.reference,
																																						  mantaIndel,
																																						  inputs.cbed,
																																						  s._strelkaWD))
					if(ret_val >> 8 != 0): raise ngs_classes.ngsExcept("[ERROR] Failed to run strelka workflow for files %s %s" % (normalFile,tumorFile))
				else:
					ret_val = os.system("python %s --normalBam %s --tumorBam %s --referenceFasta %s --indelCandidates %s --runDir %s" % (strelkaRun,
																																		 normalFile,
																																		 tumorFile,
																																		 inputs.reference,
																																		 mantaIndel,
																																		 s._strelkaWD))
					if(ret_val >> 8 != 0): raise ngs_classes.ngsExcept("[ERROR] Failed to run strelka workflow for files %s %s" % (normalFile,tumorFile))
		## Run ##
		if inputs.verbose:
			print("[INFO:COMMAND] python %s -m local -j %d" % (os.path.abspath(s._strelkaWD + '/runWorkflow.py'),inputs.threads))
			
		if not inputs.dry:
			ret_val = os.system("python %s -m local -j %d" % (os.path.abspath(s._strelkaWD + '/runWorkflow.py'),inputs.threads))
			if(ret_val >> 8 != 0): raise ngs_classes.ngsExcept("[ERROR] Failed to execute Strelka workflow script")
			
			
	print("\n<---> Done <--->\n")
	return 0


def ngs_cnvkit_somatic(inputs, targets_data_processed):
	
	print("\n<---> Cnvkit Somatic Workflow <--->\n")
	runDir = inputs.dir + '/CNVKIT_SOMATIC/'
	if(os.path.exists(runDir)):
		if inputs.verbose: print("[INFO] Cnvkit somatic directory found")
	else:
		try:
			os.mkdir(runDir)
		except OSError:
			raise ngs_classes.ngsExcept("[ERROR] Failed to create directory %s" % (runDir))

	
	## Index ##
	for s in targets_data_processed:
		if inputs.verbose: print('[INFO:COMMAND] %s index -t %d %s' % (inputs.sambamba, inputs.threads, s._bamCalibTumor))
		if not inputs.dry: os.system('%s index -t %d %s' % (inputs.sambamba, inputs.threads, s._bamCalibTumor))
		if inputs.verbose: print('[INFO:COMMAND] %s index -t %d %s' % (inputs.sambamba, inputs.threads, s._bamCalibNormal))
		if not inputs.dry: os.system('%s index -t %d %s' % (inputs.sambamba, inputs.threads, s._bamCalibNormal))
	
	normal_files = []
	tumor_files = []
	for s in targets_data_processed:
		normal_files.append(s._bamCalibNormal)
		tumor_files.append(s._bamCalibTumor)
	
	if inputs.exome:
		if inputs.verbose:
			if not(inputs.bed is None):
				print("[INFO:COMMAND] %s batch --drop-low-coverage -p %d --normal %s -f %s -t %s -d %s %s" % (inputs.cnvkit,inputs.threads," ".join(normal_files),inputs.reference,inputs.bed,runDir," ".join(tumor_files)))
			else:
				print("[ERROR:COMMAND] Cannot run copy number analysis without target file for exome data")
				return False
		if not inputs.dry:
			if not(inputs.bed is None):
				os.system("%s batch --drop-low-coverage -p %d --normal %s -f %s -t %s -d %s %s" % (inputs.cnvkit,inputs.threads," ".join(normal_files),inputs.reference,inputs.bed,runDir," ".join(tumor_files)))
			else:
				print("[ERROR:COMMAND] Cannot run copy number analysis without target file for exome data")
				return False
	else:
		if inputs.verbose:
			if not(inputs.bed is None):
				print("[INFO:COMMAND] %s batch -m wgs --drop-low-coverage -p %d --normal %s -f %s -t %s -d %s %s" % (inputs.cnvkit,inputs.threads," ".join(normal_files),inputs.reference,inputs.bed,runDir," ".join(tumor_files)))
			else:
				print("[INFO:COMMAND] %s batch -m wgs --drop-low-coverage -p %d --normal %s -f %s -d %s %s" % (inputs.cnvkit,inputs.threads," ".join(normal_files),inputs.reference,runDir," ".join(tumor_files)))
		if not inputs.dry:
			if not(inputs.bed is None):
				os.system("%s batch -m wgs --drop-low-coverage -p %d --normal %s -f %s -t %s -d %s %s" % (inputs.cnvkit,inputs.threads," ".join(normal_files),inputs.reference,inputs.bed,runDir," ".join(tumor_files)))
			else:
				os.system("%s batch -m wgs --drop-low-coverage -p %d --normal %s -f %s -d %s %s" % (inputs.cnvkit,inputs.threads," ".join(normal_files),inputs.reference,runDir," ".join(tumor_files)))
				
	print("\n<---> Done. <--->\n")
	return 0
