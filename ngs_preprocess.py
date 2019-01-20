######################################################
import sys, os, argparse
import os.path, re
import glob
import ngs_classes
######################################################


def ngs_index(inputs):
	if re.search('bowtie2', str.lower(inputs.aligner)):
		if(len(glob.glob(os.path.abspath(inputs.dir + '/indices/ref_indices_bowtie2*'))) > 0):
			print "[INFO] Indices found, skipping."
		else:
			idx_Name = os.path.abspath(inputs.dir + '/indices/ref_indices_bowtie2')
			if inputs.verbose:
				print "[INFO:COMMAND] %s --threads %d %s %s" % (inputs.bowtie2build,inputs.threads,inputs.reference,idx_Name)
			if not inputs.dry:
				ret_val = os.system("%s --threads %d %s %s" % (inputs.bowtie2build,inputs.threads,inputs.reference,idx_Name)) ## 0 return code for success
				if(ret_val >> 8 != 0): raise ngs_classes.ngsExcept("[ERROR:BOWTIE2] Failed to create fasta index")
	elif re.search("bwa", str.lower(inputs.aligner)):
		if(len(glob.glob(os.path.abspath(inputs.reference + '.bwt'))) > 0):
			print "[INFO] Indices found, skipping."
		else:
			if inputs.verbose:
				print "[INFO:COMMAND] %s index %s" % (inputs.bwa,inputs.reference)
			if not inputs.dry:
				ret_val = os.system("%s index %s" % (inputs.bwa,inputs.reference))
				if(ret_val >> 8 != 0): raise ngs_classes.ngsExcept("[ERROR:BWA] Failed to create fasta index")	
	return True

def ngs_trim(inputs, targets_data):

	if re.search('cutadapt', str.lower(inputs.trimmer)):
		for r in targets_data:
			trimmedR1_Name = os.path.abspath(inputs.dir + '/readsTrimmed/Trimmed_Read1_' + r._id + '_' + r._lib + '_' + r._lane + '.fastq')
			trimmedR2_Name = os.path.abspath(inputs.dir + '/readsTrimmed/Trimmed_Read2_' + r._id + '_' + r._lib + '_' + r._lane + '.fastq')
			
			if(os.path.exists(trimmedR1_Name) and os.path.exists(trimmedR2_Name)):
				print "[INFO] Trimmed reads found, skipping"	
			else:
				if inputs.verbose:
					print "[INFO:COMMAND] %s -q 15,10 -A XXX -m 30 -o %s -p %s %s %s" % (inputs.cutadapt,trimmedR1_Name,trimmedR2_Name,r._r1,r._r2)
				if not inputs.dry :
					ret_val = os.system("%s -q 15,10 -A XXX -m 30 -o %s -p %s %s %s" % (inputs.cutadapt,trimmedR1_Name,trimmedR2_Name,r._r1,r._r2))
					if(ret_val >> 8 != 0): raise ngs_classes.ngsExcept("[ERROR:CUTADAPT] Failed to trim reads")
			r._trimmedR1 = trimmedR1_Name
			r._trimmedR2 = trimmedR2_Name
			
	elif(re.search('afterqc', str.lower(inputs.trimmer))):
		for r in targets_data:
			stagr1 = os.path.basename(r_r1).split(".")[0]
			stagr1 = os.path.basename(r_r2).split(".")[0]
			if(os.path.exists(os.path.abspath(inputs.dir + '/readsTrimmed/' + stagr1 + '.good.fq')) and
			os.path.exists(os.path.abspath(inputs.dir + '/readsTrimmed/' + stagr2 + '.good.fq'))):
				print "    *** trimmed reads found, skipping..."
			else:
				if inputs.verbose:
					print "Running: python %s -1 %s -2 %s --debubble False -g %s --no_correction" % (inputs.afterqc,r._r1, r._r2,os.path.abspath(inputs.dir + '/readsTrimmed/'))
				if not inputs.dry:
					ret_val = os.system("python %s -1 %s -2 %s --debubble False -g %s --no_correction" % (inputs.afterqc,r._r1, r._r2,os.path.abspath(inputs.dir + '/readsTrimmed/')))
					if(ret_val >> 8 != 0): raise ngs_classes.ngsExcept("[ERROR:AFTERQC] Failed to trim reads")
			r._trimmedR1 = os.path.abspath(inputs.dir + '/readsTrimmed/' + stagr1 + '.good.fq')
			r._trimmedR2 = os.path.abspath(inputs.dir + '/readsTrimmed/' + stagr2 + '.good.fq')
	
	elif(re.search('skewer', str.lower(inputs.trimmer))):
		for r in targets_data:
			trimmedR1_Name = os.path.abspath(inputs.dir + '/readsTrimmed/' + r._id + '_' + r._lib + '_' + r._lane + '-trimmed-pair1.fastq')
			trimmedR2_Name = os.path.abspath(inputs.dir + '/readsTrimmed/' + r._id + '_' + r._lib + '_' + r._lane + '-trimmed-pair2.fastq')
			trimmed_Name = os.path.abspath(inputs.dir + '/readsTrimmed/' + r._id + '_' + r._lib + '_' + r._lane)
			
			if(os.path.exists(trimmedR1_Name) and os.path.exists(trimmedR2_Name)):
				print "[INFO] Trimmed reads found, skipping."
			else:
				if inputs.verbose:					
					print "[INFO:COMMAND] %s -m pe -q 10 -Q 20 -l 35 -n -t %d -o %s %s %s" % (inputs.skewer,inputs.threads,trimmed_Name,r._r1,r._r2)
				if not inputs.dry:
					ret_val = os.system("%s -m pe -q 10 -Q 20 -l 35 -n -t %d -o %s %s %s" % (inputs.skewer,inputs.threads,trimmed_Name,r._r1,r._r2))
					if(ret_val >> 8 != 0): raise ngs_classes.ngsExcept("[ERROR:SKEWER] Failed to trim reads %s %s" % (r._r1,r._r2))
			r._trimmedR1 = trimmedR1_Name
			r._trimmedR2 = trimmedR2_Name
			
	return True

def ngs_align(inputs, targets_data):
	
	if re.search('bowtie2', str.lower(inputs.aligner)):
		for r in targets_data:
			aligned_Name = os.path.abspath(inputs.dir + '/readsAligned/Aligned_' + r._id + '_' + r._lib + '_' + r._lane + '.sam')
			if(os.path.exists(aligned_Name)):
				print "[INFO] Aligned reads found, skipping."
			else:
				if inputs.verbose:
					print "[INFO:COMMAND] %s --threads %d -x %s --quiet --phred33 --no-mixed --no-unal "\
					"--very-sensitive-local -1 %s -2 %s -S %s --rg-id %s --rg %s --rg %s --rg %s" % (inputs.bowtie2,inputs.threads,os.path.abspath(inputs.dir + '/indices/ref_indices_bowtie2'),
					r._trimmedR1,r._trimmedR2,aligned_Name,
					r._id + '_' + r._lib + '_' + r._lane,
					'SM:' + r._id,'PL:' + str.upper(r._plat),'LB:' + r._lib)
				if not inputs.dry:
					ret_val = os.system("%s --threads %d -x %s --quiet --phred33 --no-mixed --no-unal "\
					"--very-sensitive-local -1 %s -2 %s -S %s --rg-id %s --rg %s --rg %s --rg %s" % (inputs.bowtie2,inputs.threads,os.path.abspath(inputs.dir + '/indices/ref_indices_bowtie2'),
					r._trimmedR1,r._trimmedR2,aligned_Name,
					r._id + '_' + r._lib + '_' + r._lane,'SM:' + r._id,'PL:' + str.upper(r._plat),'LB:' + r._lib)) ## 0 return code for success
					if(ret_val >> 8 != 0): raise ngs_classes.ngsExcept("[ERROR:BOWTIE2] Failed to align reads %s %s" % (r._trimmedR1, r._trimmedR2))
			r._aligned = aligned_Name
	
	elif re.search('bwa', str.lower(inputs.aligner)):
		for r in targets_data:
			aligned_Name = os.path.abspath(inputs.dir + '/readsAligned/Aligned_' + r._id + '_' + r._lib + '_' + r._lane + '.sam')
			if(os.path.exists(aligned_Name)):
				print "[INFO] Aligned reads found, skipping"
			else:
				if inputs.verbose:
					print "[INFO:COMMAND] %s mem -t %d -M -R \'@RG\\tID:%s\\tSM:%s\\tPL:%s\\tLB:%s\' %s %s %s > %s" % (inputs.bwa,inputs.threads,
					r._id + '_' + r._lib + '_' + r._lane,r._id,str.upper(r._plat),r._lib,
					inputs.reference,r._trimmedR1,r._trimmedR2,aligned_Name)
				if not inputs.dry:
					ret_val = os.system("%s mem -t %d -M -R \'@RG\\tID:%s\\tSM:%s\\tPL:%s\\tLB:%s\' %s %s %s > %s" % (inputs.bwa,inputs.threads,
					r._id + '_' + r._lib + '_' + r._lane,r._id,str.upper(r._plat),r._lib,
					inputs.reference,r._trimmedR1,r._trimmedR2,aligned_Name))
					if(ret_val >> 8 != 0): raise ngs_classes.ngsExcept("[ERROR:BWA] Failed to align reads %s %S" % (r._trimmedR1, r._trimmedR2))
			r._aligned = aligned_Name
			
	return True

def ngs_mark_sort(inputs, targets_data):
	for r in targets_data:
		bam_Name = os.path.abspath(inputs.dir + '/bamProcessed/' + r._id + '_' + r._lib + '_' + r._lane + '.bam')
		rmdup_Name = os.path.abspath(inputs.dir + '/bamProcessed/Rmdup_' + r._id + '_' + r._lib + '_' + r._lane + '.bam')
		processed_Name = os.path.abspath(inputs.dir + '/bamProcessed/Processed_' + r._id + '_' + r._lib + '_' + r._lane + '.bam')
		
		if(os.path.exists(processed_Name)):
			print "[INFO] Processed bam files found, skipping."
		else:
			
			## Convert to Bam ##
			if inputs.verbose: print "[INFO:COMMAND] %s view -h -S -t %d -f bam -o %s %s" % (inputs.sambamba,inputs.threads,bam_Name,r._aligned)
			if not inputs.dry:
				ret_val = os.system("%s view -h -S -t %d -f bam -o %s %s" % (inputs.sambamba,inputs.threads,bam_Name,r._aligned)) ## 0 return code for success
				if(ret_val >> 8 != 0): raise ngs_classes.ngsExcept("[ERROR:SAMBAMBA] Failed to convert .sam format %s" % (r._aligned))
			
			## Remove Duplicates ##
			if inputs.verbose: print "[INFO:COMMAND] %s markdup -r -t %d %s %s" % (inputs.sambamba,inputs.threads,bam_Name,rmdup_Name)
			if not inputs.dry:
				ret_val = os.system("%s markdup -r -t %d %s %s" % (inputs.sambamba,inputs.threads,bam_Name,rmdup_Name)) ## 0 return code for success
				if(ret_val >> 8 != 0): raise ngs_classes.ngsExcept("[ERROR:SAMBAMBA] Failed to mark&remove duplicates in .bam file %s" %(bam_Name))
			
			## Sort Bam File ##
			if inputs.verbose: print "[INFO:COMMAND] %s sort -t %d --tmpdir %s -o %s %s" % (inputs.sambamba,inputs.threads,inputs.dir,processed_Name,rmdup_Name)
			if not inputs.dry: 
				ret_val = os.system("%s sort -t %d --tmpdir %s -o %s %s" % (inputs.sambamba,inputs.threads,inputs.dir,processed_Name,rmdup_Name)) ## 0 return code for success
				if(ret_val >> 8 != 0): raise ngs_classes.ngsExcept("[ERROR:SAMBAMBA] Failed to sort .bam file %s" %(rmdup_Name))
			
		r._processed = processed_Name
	
	return True

def ngs_merge(inputs,targets_data,targets_data_processed):
	if str.upper(inputs.workflow) == 'SOMATICSNVINDEL':
		sample_ids = {}
		for r in targets_data:
			if r._id in sample_ids.keys():
				if str.upper(r._type) == "TUMOR":
					sample_ids[r._id]["TUMOR"].append(r._processed)
				elif str.upper(r._type) == "NORMAL":
					sample_ids[r._id]["NORMAL"].append(r._processed)
			else:
				sample_ids[r._id] = {"TUMOR" : [], "NORMAL" : []}
				if str.upper(r._type) == "TUMOR":
					sample_ids[r._id]["TUMOR"].append(r._processed)
				elif str.upper(r._type) == "NORMAL":
					sample_ids[r._id]["NORMAL"].append(r._processed)
		
		for s in sample_ids:
			## TUMOR ##
			if(len(sample_ids[s]["TUMOR"]) == 1):
				print "[INFO] Single bam file for sample %s skipping." % (s)
				merged_tumor_bam_name = sample_ids[s]["TUMOR"][0]
			else:
				merged_tumor_bam_name = os.path.abspath(inputs.dir + '/bamMerged/Merged_Tumor_' + s + '.bam')
				if(os.path.exists(merged_tumor_bam_name)):
					print "[INFO] Merged .bam file found for tumor sample, skipping."
				else:
					if inputs.verbose: print "[INFO:COMMAND] %s merge -t %d %s %s" % (inputs.sambamba,inputs.threads,merged_tumor_bam_name,' '.join(sample_ids[s]["TUMOR"]))
					if not inputs.dry:
						ret_val = os.system("%s merge -t %d %s %s" % (inputs.sambamba,inputs.threads,merged_tumor_bam_name,' '.join(sample_ids[s]["TUMOR"])))
						if(ret_val >> 8 != 0): raise ngs_classes.ngsExcept("[ERROR:SAMBAMBA] Failed to merge files %s" % (' '.join(sample_ids[s]["TUMOR"])))

			## NORMAL ##
			if(len(sample_ids[s]["NORMAL"]) == 1):
				print "[INFO] Single bam file for sample %s skipping." % (s)
				merged_normal_bam_name = sample_ids[s]["NORMAL"][0]
			else:
				merged_normal_bam_name = os.path.abspath(inputs.dir + '/bamMerged/Merged_Normal_' + s + '.bam')
				if(os.path.exists(merged_normal_bam_name)):
					print "[INFO] Merged .bam file found for normal sample, skipping"
				else:
					if inputs.verbose: print "[INFO:COMMAND] %s merge -t %d %s %s" % (inputs.sambamba,inputs.threads,merged_normal_bam_name,' '.join(sample_ids[s]["NORMAL"]))
					if not inputs.dry:
						ret_val = os.system("%s merge -t %d %s %s" % (inputs.sambamba,inputs.threads,merged_normal_bam_name,' '.join(sample_ids[s]["NORMAL"])))
						if(ret_val >> 8 != 0): raise ngs_classes.ngsExcept("[ERROR:SAMBAMA] Failed to merge files %s" % (' '.join(sample_ids[s]["TUMOR"])))
			targets_data_processed.append(ngs_classes.WxsSample(_id=s,_bamTumor = merged_tumor_bam_name,_bamNormal=merged_normal_bam_name))            
	elif str.upper(inputs.workflow) == "GERMLINESNVINDEL":
		sample_ids = {}
		for r in targets_data:
			if r._id in sample_ids.keys():
				sample_ids[r._id].append(r._processed)
			else:
				sample_ids[r._id] = []
				sample_ids[r._id].append(r._processed)
		for s in sample_ids:
			if(len(sample_ids[s]) == 1):
				print "[INFO] Single bam file for sample %s skipping." % (s)
				merged_bam_name = sample_ids[s][0]
			else:
				merged_bam_name = os.path.abspath(inputs.dir + '/bamMerged/Merged_Normal_' + s + '.bam')
				if(os.path.exists(merged_bam_name)):
					print "[INFO] Merged bam files found, skipping."
				else:
					if inputs.verbose: print "[INFO:COMMAND] %s merge -t %d %s %s" % (inputs.sambamba,inputs.threads,merged_bam_name,' '.join(sample_ids[s]))
					if not inputs.dry: 
						ret_val = os.system("%s merge -t %d %s %s" % (inputs.sambamba,inputs.threads,merged_bam_name,' '.join(sample_ids[s])))
						if(ret_val >> 8 != 0): raise ngs_classes.ngsExcept("[ERROR:SAMBAMBA] Failed to merge files %s" % (' '.join(sample_ids[s])))
			targets_data_processed.append(ngs_classes.WxsSample(_id = s,_bam = merged_bam_name))
	
	return 0

def bqsr(inputs, targets_data_processed):
	
	for s in targets_data_processed:
		if str.upper(inputs.workflow) == "SOMATICSNVINDEL":
			if inputs.verbose: print "[INFO] %s index -t %d %s" % (inputs.sambamba,inputs.threads,s._bamTumor)
			if not inputs.dry: os.system("%s index -t %d %s" % (inputs.sambamba,inputs.threads,s._bamTumor))
			if inputs.verbose: print "[INFO] %s index -t %d %s" % (inputs.sambamba,inputs.threads,s._bamNormal)
			if not inputs.dry: os.system("%s index -t %d %s" % (inputs.sambamba,inputs.threads,s._bamNormal))
		elif str.upper(inputs.workflow) == "GERMLINESNVINDEL":
			if inputs.verbose: print "[INFO] %s index -t %d %s" % (inputs.sambamba,inputs.threads,s._bam)
			if not inputs.dry: os.system("%s index -t %d %s" % (inputs.sambamba,inputs.threads,s._bam))
			
	## BQSR (Step 1) ##
	for s in targets_data_processed:
		if str.upper(inputs.workflow) == "SOMATICSNVINDEL":
			tumorRecalibTable = os.path.abspath(inputs.dir + '/recalib/RecalibrationTableTumor_' + s._id + '.table')
			normalRecalibTable = os.path.abspath(inputs.dir + '/recalib/RecalibrationTableNormal_' + s._id + '.table')
			if(os.path.exists(tumorRecalibTable) and os.path.exists(normalRecalibTable)):
				print "[INFO] Recalibration tables found, skipping."
			else: ## 0 return code for success
				if inputs.verbose:
					if not(inputs.bed is None):
						print "[INFO:COMMAND] java -jar %s -T BaseRecalibrator -nct %d -R %s -I %s -knownSites %s -o %s -L %s" % (inputs.gatk,inputs.threads,inputs.reference,s._bamTumor,' -knownSites '.join(inputs.knownsites),tumorRecalibTable,inputs.bed)
					else:
						print "[INFO:COMMAND] java -jar %s -T BaseRecalibrator -nct %d -R %s -I %s -knownSites %s -o %s" % (inputs.gatk,inputs.threads,inputs.reference,s._bamTumor,' -knownSites '.join(inputs.knownsites),tumorRecalibTable)
				if not inputs.dry:
					if not(inputs.bed is None):
						os.system("java -jar %s -T BaseRecalibrator -nct %d -R %s -I %s -knownSites %s -o %s -L %s" % (inputs.gatk,inputs.threads,inputs.reference,s._bamTumor,' -knownSites '.join(inputs.knownsites),tumorRecalibTable,inputs.bed))
					else:
						os.system("java -jar %s -T BaseRecalibrator -nct %d -R %s -I %s -knownSites %s -o %s" % (inputs.gatk,inputs.threads,inputs.reference,s._bamTumor,' -knownSites '.join(inputs.knownsites),tumorRecalibTable))
				s._recalibTableTumor = tumorRecalibTable
				
				if inputs.verbose:
					if not(inputs.bed is None):
						print "[INFO:COMMAND] java -jar %s -T BaseRecalibrator -nct %d -R %s -I %s -knownSites %s -o %s -L %s" % (inputs.gatk,inputs.threads,inputs.reference,s._bamNormal,' -knownSites '.join(inputs.knownsites),normalRecalibTable,inputs.bed)
					else:
						print "[INFO:COMMAND] java -jar %s -T BaseRecalibrator -nct %d -R %s -I %s -knownSites %s -o %s" % (inputs.gatk,inputs.threads,inputs.reference,s._bamNormal,' -knownSites '.join(inputs.knownsites),normalRecalibTable)
				if not inputs.dry:
					if not(inputs.bed is None):
						os.system("java -jar %s -T BaseRecalibrator -nct %d -R %s -I %s -knownSites %s -o %s -L %s" % (inputs.gatk,inputs.threads,inputs.reference,s._bamNormal,' -knownSites '.join(inputs.knownsites),normalRecalibTable,inputs.bed))
					else:
						os.system("java -jar %s -T BaseRecalibrator -nct %d -R %s -I %s -knownSites %s -o %s" % (inputs.gatk,inputs.threads,inputs.reference,s._bamNormal,' -knownSites '.join(inputs.knownsites),normalRecalibTable))
				s._recalibTableNormal = normalRecalibTable
		elif str.upper(inputs.workflow) == "GERMLINESNVINDEL":
			RecalibTable = os.path.abspath(inputs.dir + '/recalib/RecalibrationTableNormal_' + s._id + '.table')
			if(os.path.exists(RecalibTable)):
				print "[INFO] Recalibration table found, skipping."
			else:
				if inputs.verbose:
					if not(inputs.bed is None):
						print "[INFO:COMMAND] java -jar %s -T BaseRecalibrator -nct %d -R %s -I %s -knownSites %s -o %s -L %s" % (inputs.gatk,inputs.threads,inputs.reference,s._bam,' -knownSites '.join(inputs.knownsites),RecalibTable,inputs.bed)
					else:
						print "[INFO:COMMAND] java -jar %s -T BaseRecalibrator -nct %d -R %s -I %s -knownSites %s -o %s" % (inputs.gatk,inputs.threads,inputs.reference,s._bam,' -knownSites '.join(inputs.knownsites),RecalibTable)
				if not inputs.dry:
					if not(inputs.bed is None):
						os.system("java -jar %s -T BaseRecalibrator -nct %d -R %s -I %s -knownSites %s -o %s -L %s" % (inputs.gatk,inputs.threads,inputs.reference,s._bam,' -knownSites '.join(inputs.knownsites),RecalibTable,inputs.bed))
					else:
						os.system("java -jar %s -T BaseRecalibrator -nct %d -R %s -I %s -knownSites %s -o %s" % (inputs.gatk,inputs.threads,inputs.reference,s._bam,' -knownSites '.join(inputs.knownsites),RecalibTable))
				s._recalibTable = RecalibTable
				
	## BQSR (Step 2) ##
	for s in targets_data_processed:
		if str.upper(inputs.workflow) == "SOMATICSNVINDEL":
			tumorCalibBam = os.path.abspath(inputs.dir + '/bamCalib/CalibratedTumor_' + s._id + '.bam')
			normalCalibBam = os.path.abspath(inputs.dir + '/bamCalib/CalibratedNormal_' + s._id + '.bam')
			if(os.path.exists(tumorCalibBam) and os.path.exists(normalCalibBam)):
				print "[INFO] BQSR adjusted files found, skipping."
			else:
				if inputs.verbose: ## 0 return code for Success
					if not(inputs.bed is None):
						print "[INFO:COMMAND] java -jar %s -T PrintReads -nct %d -R %s -I %s -BQSR %s -o %s -L %s" % (inputs.gatk,inputs.threads,inputs.reference,s._bamTumor,s._recalibTableTumor,tumorCalibBam,inputs.bed)
					else:
						print "[INFO:COMMAND] java -jar %s -T PrintReads -nct %d -R %s -I %s -BQSR %s -o %s" % (inputs.gatk,inputs.threads,inputs.reference,s._bamTumor,s._recalibTableTumor,tumorCalibBam)
				if not inputs.dry:
					if not(inputs.bed is None):
						os.system("java -jar %s -T PrintReads -nct %d -R %s -I %s -BQSR %s -o %s -L %s" % (inputs.gatk,inputs.threads,inputs.reference,s._bamTumor,s._recalibTableTumor,tumorCalibBam,inputs.bed))
					else:
						os.system("java -jar %s -T PrintReads -nct %d -R %s -I %s -BQSR %s -o %s" % (inputs.gatk,inputs.threads,inputs.reference,s._bamTumor,s._recalibTableTumor,tumorCalibBam))
				
				if inputs.verbose:
					if not(inputs.bed is None):
						print "[INFO:COMMAND] java -jar %s -T PrintReads -nct %d -R %s -I %s -BQSR %s -o %s -L %s" % (inputs.gatk,inputs.threads,inputs.reference,s._bamNormal,s._recalibTableNormal,normalCalibBam,inputs.bed)
					else:
						print "[INFO:COMMAND] java -jar %s -T PrintReads -nct %d -R %s -I %s -BQSR %s -o %s" % (inputs.gatk,inputs.threads,inputs.reference,s._bamNormal,s._recalibTableNormal,normalCalibBam)
				if not inputs.dry:
					if not(inputs.bed is None):
						os.system("java -jar %s -T PrintReads -nct %d -R %s -I %s -BQSR %s -o %s -L %s" % (inputs.gatk,inputs.threads,inputs.reference,s._bamNormal,s._recalibTableNormal,normalCalibBam,inputs.bed))
					else:
						os.system("java -jar %s -T PrintReads -nct %d -R %s -I %s -BQSR %s -o %s" % (inputs.gatk,inputs.threads,inputs.reference,s._bamTumor,s._recalibTableTumor,normalCalibBam))
			s._bamCalibTumor = tumorCalibBam
			s._bamCalibNormal = normalCalibBam
		elif str.upper(inputs.workflow) == "GERMLINESNVINDEL":
			CalibBam = os.path.abspath(inputs.dir + '/bamCalib/CalibratedNormal_' + s._id + '.bam')
			if(os.path.exists(CalibBam)):
				print "[INFO] BQSR adjusted files found, skipping."
			else:
				if inputs.verbose:
					if not(inputs.bed is None):
						print "[INFO:COMMAND] java -jar %s -T PrintReads -nct %d -R %s -I %s -BQSR %s -o %s -L %s" % (inputs.gatk,inputs.threads,inputs.reference,s._bam,s._recalibTable,CalibBam,inputs.bed)
					else:
						print "[INFO:COMMAND] java -jar %s -T PrintReads -nct %d -R %s -I %s -BQSR %s -o %s" % (inputs.gatk,inputs.threads,inputs.reference,s._bam,s._recalibTable,CalibBam)
				if not inputs.dry:
					if not(inputs.bed is None):
						os.system("java -jar %s -T PrintReads -nct %d -R %s -I %s -BQSR %s -o %s -L %s" % (inputs.gatk,inputs.threads,inputs.reference,s._bam,s._recalibTable,CalibBam,inputs.bed))
					else:
						os.system("java -jar %s -T PrintReads -nct %d -R %s -I %s -BQSR %s -o %s" % (inputs.gatk,inputs.threads,inputs.reference,s._bam,s._recalibTable,CalibBam))
			s._bamCalib = CalibBam
				
	return True
