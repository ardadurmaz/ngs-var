######################################################
import sys, os, argparse
import os.path, re
import glob
import ngs_classes
from ngs_functions import print_log, run_command
######################################################


def ngs_index(inputs):
	if re.search("bwa", str.lower(inputs.aligner)):
		if(len(glob.glob(os.path.abspath(inputs.reference + '.bwt'))) > 0):
			print_log(inputs, "[INFO] Indices found, skipping.")
		else:
			run_command(inputs, "%s index %s" % (inputs.bwa,inputs.reference), "[ERROR:BWA] Failed to create fasta index")
	elif re.search("star", str.lower(inputs.aligner)):
		print_log(inputs, "\n<---> RNA Indexing <--->\n")
		runDir = os.path.abspath(inputs.dir + '/RNA_SEQ/')
		error_msg = "[ERROR] Failed to perform RNA indexing."
		if(os.path.exists(runDir + '/Genome')):
			print_log(inputs, "[INFO] Indices found, skipping.")
			return 0
		else:
			run_command(inputs, "%s --runThreadN %d --runMode genomeGenerate --genomeDir %s --genomeFastaFiles %s --sjdbGTFfile %s" % (inputs.star, inputs.threads, runDir, inputs.reference, inputs.annotation), error_msg)
	return True

def ngs_trim(inputs, targets_data):
	if re.search('fastp', str.lower(inputs.trimmer)):
		for r in targets_data:
			trimmedR1_Name = os.path.abspath(inputs.dir + '/readsTrimmed/Trimmed_Read1_' + r._id + '_' + r._lib + '_' + r._lane + '.fastq.gz')
			trimmedR2_Name = os.path.abspath(inputs.dir + '/readsTrimmed/Trimmed_Read2_' + r._id + '_' + r._lib + '_' + r._lane + '.fastq.gz')
			
			if(os.path.exists(trimmedR1_Name) and os.path.exists(trimmedR2_Name)):
				print_log(inputs, "[INFO] Trimmed reads found, skipping")	
			else:
				run_command(inputs, "%s --thread 12 --compression 9 --trim_front1 10 --trim_front2 10 --cut_tail --detect_adapter_for_pe --in1 %s --in2 %s --out1 %s --out2 %s" % (inputs.fastp,r._r1,r._r2, trimmedR1_Name,trimmedR2_Name), "[ERROR:CUTADAPT] Failed to trim reads")
			r._trimmedR1 = trimmedR1_Name
			r._trimmedR2 = trimmedR2_Name
	elif re.search('cutadapt', str.lower(inputs.trimmer)):
		for r in targets_data:
			trimmedR1_Name = os.path.abspath(inputs.dir + '/readsTrimmed/Trimmed_Read1_' + r._id + '_' + r._lib + '_' + r._lane + '.fastq.gz')
			trimmedR2_Name = os.path.abspath(inputs.dir + '/readsTrimmed/Trimmed_Read2_' + r._id + '_' + r._lib + '_' + r._lane + '.fastq.gz')
			
			if(os.path.exists(trimmedR1_Name) and os.path.exists(trimmedR2_Name)):
				print_log(inputs, "[INFO] Trimmed reads found, skipping")	
			else:
				run_command(inputs, "%s -q 15,10 -A XXX -m 30 -o %s -p %s %s %s" % (inputs.cutadapt,trimmedR1_Name,trimmedR2_Name,r._r1,r._r2), "[ERROR:CUTADAPT] Failed to trim reads")
			r._trimmedR1 = trimmedR1_Name
			r._trimmedR2 = trimmedR2_Name
			
	elif(re.search('afterqc', str.lower(inputs.trimmer))):
		for r in targets_data:
			stagr1 = os.path.basename(r._r1).split(".")[0]
			stagr2 = os.path.basename(r._r2).split(".")[0]
			if(os.path.exists(os.path.abspath(inputs.dir + '/readsTrimmed/' + stagr1 + '.good.fq')) and
			os.path.exists(os.path.abspath(inputs.dir + '/readsTrimmed/' + stagr2 + '.good.fq'))):
				print_log(inputs, "    *** trimmed reads found, skipping...")
			else:
				run_command(inputs, "python %s -1 %s -2 %s --debubble False -g %s --no_correction" % (inputs.afterqc,r._r1, r._r2,os.path.abspath(inputs.dir + '/readsTrimmed/')), "[ERROR:AFTERQC] Failed to trim reads")
			r._trimmedR1 = os.path.abspath(inputs.dir + '/readsTrimmed/' + stagr1 + '.good.fq')
			r._trimmedR2 = os.path.abspath(inputs.dir + '/readsTrimmed/' + stagr2 + '.good.fq')
	
	elif(re.search('skewer', str.lower(inputs.trimmer))):
		for r in targets_data:
			trimmedR1_Name = os.path.abspath(inputs.dir + '/readsTrimmed/' + r._id + '_' + r._lib + '_' + r._lane + '-trimmed-pair1.fastq.gz')
			trimmedR2_Name = os.path.abspath(inputs.dir + '/readsTrimmed/' + r._id + '_' + r._lib + '_' + r._lane + '-trimmed-pair2.fastq.gz')
			trimmed_Name = os.path.abspath(inputs.dir + '/readsTrimmed/' + r._id + '_' + r._lib + '_' + r._lane)
			
			if(os.path.exists(trimmedR1_Name) and os.path.exists(trimmedR2_Name)):
				print_log(inputs, "[INFO] Trimmed reads found, skipping.")
			else:
				run_command(inputs, "%s -m pe -q 10 -Q 20 -l 35 -n -t %d -o %s %s %s" % (inputs.skewer,inputs.threads,trimmed_Name,r._r1,r._r2), "[ERROR:SKEWER] Failed to trim reads %s %s" % (r._r1,r._r2))
			r._trimmedR1 = trimmedR1_Name
			r._trimmedR2 = trimmedR2_Name
			
	return True

def ngs_align(inputs, targets_data):
	if re.search('bwa', str.lower(inputs.aligner)):
		for r in targets_data:
			aligned_Name = os.path.abspath(inputs.dir + '/readsAligned/Aligned_' + r._id + '_' + r._lib + '_' + r._lane + '.bam')
			if(os.path.exists(aligned_Name)):
				print_log(inputs, "[INFO] Aligned reads found, skipping")
			else:
				run_command(inputs, "/bin/bash -c '%s mem -t %d -M -R \"@RG\\tID:%s\\tSM:%s\\tPL:%s\\tLB:%s\" %s <(gunzip -c %s) <(gunzip -c %s) | samtools sort -@%d -l9 -o %s'" % (inputs.bwa,inputs.threads,
						r._id + '_' + r._lib + '_' + r._lane,r._id,str.upper(r._plat),r._lib,
						inputs.reference,r._trimmedR1,r._trimmedR2, inputs.threads, aligned_Name), "[ERROR:BWA] Failed to align reads %s %s" % (r._trimmedR1, r._trimmedR2))
			r._aligned = aligned_Name
	elif re.search("star", str.lower(inputs.aligner)):
		print_log(inputs, "\n<---> RNA Alignment <--->\n") 
		runDir = os.path.abspath(inputs.dir + '/RNA_SEQ/')
		tempDir = "~/star_temp/"
		
		sample_reads = {}
		for s in targets_data:
			if s._id not in sample_reads:
				sample_reads[s._id] = {'r1': [], 'r2': []}
			sample_reads[s._id]['r1'].append(s._trimmedR1)
			sample_reads[s._id]['r2'].append(s._trimmedR2)

		for sample_id, reads in sample_reads.items():
			sampleName = f"RNA_Seq_{sample_id}"
			r1 = ",".join(reads['r1'])
			r2 = ",".join(reads['r2'])
			readFilesIn = " ".join([r1, r2])

			error_msg = f"[ERROR] Failed to perform RNA alignment for sample {sample_id}."
			run_command(inputs, "%s --runThreadN %d --genomeDir %s --readFilesIn %s --outTmpDir %s --readFilesCommand zcat --twopassMode Basic --quantMode TranscriptomeSAM --outSAMtype BAM Unsorted --outBAMcompression 10 --outFileNamePrefix %s" % (inputs.star, inputs.threads, runDir, readFilesIn, tempDir, sampleName), error_msg)
	return True

def ngs_mark_sort(inputs, targets_data):
	if str.upper(inputs.workflow) != 'RNASEQ':
		for r in targets_data:
			rmdup_Name = os.path.abspath(inputs.dir + '/bamProcessed/Rmdup_' + r._id + '_' + r._lib + '_' + r._lane + '.bam')
			
			if(os.path.exists(rmdup_Name)):
				print_log(inputs, "[INFO] Processed bam files found, skipping.")
			else:
				## Remove Duplicates ##
				run_command(inputs, "%s markdup --overflow-list-size=800000 --io-buffer-size=512 -l9 -p -t %d %s %s" % (inputs.sambamba,inputs.threads,r._aligned,rmdup_Name), "[ERROR:SAMBAMBA] Failed to mark&remove duplicates in .bam file")
			r._processed = rmdup_Name
	
	return True

def ngs_merge(inputs,targets_data,targets_data_processed):
	if str.upper(inputs.workflow) == 'SOMATICSNVINDEL':
		sample_ids = {}
		for r in targets_data:
			if r._id in list(sample_ids.keys()):
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
				print_log(inputs, "[INFO] Single bam file for sample %s skipping." % (s))
				merged_tumor_bam_name = sample_ids[s]["TUMOR"][0]
			else:
				merged_tumor_bam_name = os.path.abspath(inputs.dir + '/bamMerged/Merged_Tumor_' + s + '.bam')
				if(os.path.exists(merged_tumor_bam_name)):
					print_log(inputs, "[INFO] Merged .bam file found for tumor sample, skipping.")
				else:
					run_command(inputs, "%s merge -l9 -t %d %s %s" % (inputs.sambamba,inputs.threads,merged_tumor_bam_name,' '.join(sample_ids[s]["TUMOR"])), "[ERROR:SAMBAMBA] Failed to merge files %s" % (' '.join(sample_ids[s]["TUMOR"])))

			## NORMAL ##
			if(len(sample_ids[s]["NORMAL"]) == 1):
				print_log(inputs, "[INFO] Single bam file for sample %s skipping." % (s))
				merged_normal_bam_name = sample_ids[s]["NORMAL"][0]
			else:
				merged_normal_bam_name = os.path.abspath(inputs.dir + '/bamMerged/Merged_Normal_' + s + '.bam')
				if(os.path.exists(merged_normal_bam_name)):
					print_log(inputs, "[INFO] Merged .bam file found for normal sample, skipping")
				else:
					run_command(inputs, "%s merge -l9 -t %d %s %s" % (inputs.sambamba,inputs.threads,merged_normal_bam_name,' '.join(sample_ids[s]["NORMAL"])), "[ERROR:SAMBAMA] Failed to merge files %s" % (' '.join(sample_ids[s]["TUMOR"])))
			targets_data_processed.append(ngs_classes.WxsSample(_id=s,_bamTumor = merged_tumor_bam_name,_bamNormal=merged_normal_bam_name))            
	elif str.upper(inputs.workflow) == "GERMLINESNVINDEL":
		sample_ids = {}
		for r in targets_data:
			if r._id in list(sample_ids.keys()):
				sample_ids[r._id].append(r._processed)
			else:
				sample_ids[r._id] = []
				sample_ids[r._id].append(r._processed)
		for s in sample_ids:
			if(len(sample_ids[s]) == 1):
				print_log(inputs, "[INFO] Single bam file for sample %s skipping." % (s))
				merged_bam_name = sample_ids[s][0]
			else:
				merged_bam_name = os.path.abspath(inputs.dir + '/bamMerged/Merged_' + s + '.bam')
				if(os.path.exists(merged_bam_name)):
					print_log(inputs, "[INFO] Merged bam files found, skipping.")
				else:
					run_command(inputs, "%s merge -l9 -t %d %s %s" % (inputs.sambamba,inputs.threads,merged_bam_name,' '.join(sample_ids[s])), "[ERROR:SAMBAMBA] Failed to merge files %s" % (' '.join(sample_ids[s])))
			targets_data_processed.append(ngs_classes.WxsSample(_id = s,_bam = merged_bam_name))
	
	return 0

def bqsr(inputs, targets_data_processed):
	
	for s in targets_data_processed:
		if str.upper(inputs.workflow) == "SOMATICSNVINDEL":
			if(len(glob.glob(os.path.abspath(s._bamTumor + '.bai'))) > 0):
				print_log(inputs, "[INFO] Indices found, skipping.")
			else:
				run_command(inputs, "%s index -t %d %s" % (inputs.sambamba,inputs.threads,s._bamTumor))
			if(len(glob.glob(os.path.abspath(s._bamNormal + '.bai'))) > 0):
				print_log(inputs, "[INFO] Indices found, skipping.")
			else:
				run_command(inputs, "%s index -t %d %s" % (inputs.sambamba,inputs.threads,s._bamNormal))
		elif str.upper(inputs.workflow) == "GERMLINESNVINDEL":
			if(len(glob.glob(os.path.abspath(s._bam + '.bai'))) > 0):
				print_log(inputs, "[INFO] Indices found, skipping.")
			else:
				run_command(inputs, "%s index -t %d %s" % (inputs.sambamba,inputs.threads,s._bam))
			
	## BQSR (Step 1) ##
	for s in targets_data_processed:
		if str.upper(inputs.workflow) == "SOMATICSNVINDEL":
			tumorRecalibTable = os.path.abspath(inputs.dir + '/recalib/RecalibrationTableTumor_' + s._id + '.table')
			normalRecalibTable = os.path.abspath(inputs.dir + '/recalib/RecalibrationTableNormal_' + s._id + '.table')
			if(os.path.exists(tumorRecalibTable) and os.path.exists(normalRecalibTable)):
				print_log(inputs, "[INFO] Recalibration tables found, skipping.")
			else: ## 0 return code for success
				if not(inputs.bed is None):
					run_command(inputs, "java -jar %s -T BaseRecalibrator -nct %d -R %s -I %s -knownSites %s -o %s -L %s" % (inputs.gatk,inputs.threads,inputs.reference,s._bamTumor,' -knownSites '.join(inputs.knownsites),tumorRecalibTable,inputs.bed))
				else: 
					run_command(inputs, "java -jar %s -T BaseRecalibrator -nct %d -R %s -I %s -knownSites %s -o %s" % (inputs.gatk,inputs.threads,inputs.reference,s._bamTumor,' -knownSites '.join(inputs.knownsites),tumorRecalibTable))

				s._recalibTableTumor = tumorRecalibTable
				
				if not(inputs.bed is None):
					run_command(inputs, "java -jar %s -T BaseRecalibrator -nct %d -R %s -I %s -knownSites %s -o %s -L %s" % (inputs.gatk,inputs.threads,inputs.reference,s._bamNormal,' -knownSites '.join(inputs.knownsites),normalRecalibTable,inputs.bed))
				else:
					run_command(inputs, "java -jar %s -T BaseRecalibrator -nct %d -R %s -I %s -knownSites %s -o %s" % (inputs.gatk,inputs.threads,inputs.reference,s._bamNormal,' -knownSites '.join(inputs.knownsites),normalRecalibTable))
				s._recalibTableNormal = normalRecalibTable
		elif str.upper(inputs.workflow) == "GERMLINESNVINDEL":
			RecalibTable = os.path.abspath(inputs.dir + '/recalib/RecalibrationTableNormal_' + s._id + '.table')
			if(os.path.exists(RecalibTable)):
				print_log(inputs, "[INFO] Recalibration table found, skipping.")
			else:
				if not(inputs.bed is None):
					run_command(inputs, "java -jar %s -T BaseRecalibrator -nct %d -R %s -I %s -knownSites %s -o %s -L %s" % (inputs.gatk,inputs.threads,inputs.reference,s._bam,' -knownSites '.join(inputs.knownsites),RecalibTable,inputs.bed))
				else:
					run_command(inputs, "java -jar %s -T BaseRecalibrator -nct %d -R %s -I %s -knownSites %s -o %s" % (inputs.gatk,inputs.threads,inputs.reference,s._bam,' -knownSites '.join(inputs.knownsites),RecalibTable))
				s._recalibTable = RecalibTable
				
	## BQSR (Step 2) ##
	for s in targets_data_processed:
		if str.upper(inputs.workflow) == "SOMATICSNVINDEL":
			tumorCalibBam = os.path.abspath(inputs.dir + '/bamCalib/CalibratedTumor_' + s._id + '.bam')
			normalCalibBam = os.path.abspath(inputs.dir + '/bamCalib/CalibratedNormal_' + s._id + '.bam')
			if(os.path.exists(tumorCalibBam) and os.path.exists(normalCalibBam)):
				print_log(inputs, "[INFO] BQSR adjusted files found, skipping.")
			else:
				if not(inputs.bed is None):
					run_command(inputs, "java -jar %s -T PrintReads -nct %d -R %s -I %s -BQSR %s -o %s -L %s" % (inputs.gatk,inputs.threads,inputs.reference,s._bamTumor,s._recalibTableTumor,tumorCalibBam,inputs.bed))
				else:
					run_command(inputs, "java -jar %s -T PrintReads -nct %d -R %s -I %s -BQSR %s -o %s" % (inputs.gatk,inputs.threads,inputs.reference,s._bamTumor,s._recalibTableTumor,tumorCalibBam))
				
				if not(inputs.bed is None):
					run_command(inputs, "java -jar %s -T PrintReads -nct %d -R %s -I %s -BQSR %s -o %s -L %s" % (inputs.gatk,inputs.threads,inputs.reference,s._bamNormal,s._recalibTableNormal,normalCalibBam,inputs.bed))
				else:
					run_command(inputs, "java -jar %s -T PrintReads -nct %d -R %s -I %s -BQSR %s -o %s" % (inputs.gatk,inputs.threads,inputs.reference,s._bamTumor,s._recalibTableTumor,normalCalibBam))
			s._bamCalibTumor = tumorCalibBam
			s._bamCalibNormal = normalCalibBam
		elif str.upper(inputs.workflow) == "GERMLINESNVINDEL":
			CalibBam = os.path.abspath(inputs.dir + '/bamCalib/CalibratedNormal_' + s._id + '.bam')
			if(os.path.exists(CalibBam)):
				print_log(inputs, "[INFO] BQSR adjusted files found, skipping.")
			else:
				if not(inputs.bed is None):
					run_command(inputs, "java -jar %s -T PrintReads -nct %d -R %s -I %s -BQSR %s -o %s -L %s" % (inputs.gatk,inputs.threads,inputs.reference,s._bam,s._recalibTable,CalibBam,inputs.bed))
				else:
					run_command(inputs, "java -jar %s -T PrintReads -nct %d -R %s -I %s -BQSR %s -o %s" % (inputs.gatk,inputs.threads,inputs.reference,s._bam,s._recalibTable,CalibBam))
			s._bamCalib = CalibBam
				
	return True
