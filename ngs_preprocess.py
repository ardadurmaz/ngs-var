######################################################
import sys, os, argparse
import os.path, re
import glob
import ngs_classes
from ngs_functions import print_log, run_command
######################################################

def ngs_index(inputs):
    if re.search("bwa", str.lower(inputs.aligner)):      
        if len(glob.glob(os.path.abspath(inputs.reference + '.bwt'))) > 0:
            print_log(inputs, "[INFO] BWA indices found, skipping.")
        else:
            run_command(inputs, f"{inputs.bwa} index {inputs.reference}", "[ERROR:BWA] Failed to create fasta index")
        
        if len(glob.glob(os.path.abspath(inputs.reference + '.fai'))) > 0:
            print_log(inputs, "[INFO] FAI index found, skipping.")
        else:
            run_command(inputs, f"{inputs.samtools} faidx {inputs.reference}", "[ERROR:SAMTOOLS] Failed to create fasta index")
    
    elif re.search("star", str.lower(inputs.aligner)):
        print_log(inputs, "\n<---> RNA Indexing <--->\n")
        runDir = os.path.abspath(inputs.dir + '/RNA_SEQ/')
        error_msg = "[ERROR] Failed to perform RNA indexing."
        
        if os.path.exists(runDir + '/Genome'):
            print_log(inputs, "[INFO] Indices found, skipping.")
            return 0
        else:
            run_command(inputs, f"{inputs.star} --runThreadN {inputs.threads} --runMode genomeGenerate --genomeDir {runDir} --genomeFastaFiles {inputs.reference} --sjdbGTFfile {inputs.annotation}", error_msg)
    
    return True


def ngs_trim(inputs, targets_data):
	if re.search('fastp', str.lower(inputs.trimmer)):
		for r in targets_data:
			trimmedR1_Name = os.path.abspath(inputs.dir + '/readsTrimmed/Trimmed_Read1_' + r._id + '_' + r._lib + '_' + r._lane + '.fastq.gz')
			trimmedR2_Name = os.path.abspath(inputs.dir + '/readsTrimmed/Trimmed_Read2_' + r._id + '_' + r._lib + '_' + r._lane + '.fastq.gz')
			
			if(os.path.exists(trimmedR1_Name) and os.path.exists(trimmedR2_Name)):
				print_log(inputs, "[INFO] Trimmed reads found, skipping")	
			else:
				run_command(inputs, f"{inputs.fastp} --thread 12 --compression 9 --trim_front1 10 --trim_front2 10 --cut_tail --detect_adapter_for_pe --in1 {r._r1} --in2 {r._r2} --out1 {trimmedR1_Name} --out2 {trimmedR2_Name}", "[ERROR:CUTADAPT] Failed to trim reads")
			r._trimmedR1 = trimmedR1_Name
			r._trimmedR2 = trimmedR2_Name
	elif re.search('cutadapt', str.lower(inputs.trimmer)):
		for r in targets_data:
			trimmedR1_Name = os.path.abspath(inputs.dir + '/readsTrimmed/Trimmed_Read1_' + r._id + '_' + r._lib + '_' + r._lane + '.fastq.gz')
			trimmedR2_Name = os.path.abspath(inputs.dir + '/readsTrimmed/Trimmed_Read2_' + r._id + '_' + r._lib + '_' + r._lane + '.fastq.gz')
			
			if(os.path.exists(trimmedR1_Name) and os.path.exists(trimmedR2_Name)):
				print_log(inputs, "[INFO] Trimmed reads found, skipping")	
			else:
				run_command(inputs, f"{inputs.cutadapt} -q 15,10 -A XXX -m 30 -o {trimmedR1_Name} -p {trimmedR2_Name} {r._r1} {r._r2}", "[ERROR:CUTADAPT] Failed to trim reads")
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
				run_command(inputs, f"python {inputs.afterqc} -1 {r._r1} -2 {r._r2} --debubble False -g {os.path.abspath(inputs.dir + '/readsTrimmed/')} --no_correction", "[ERROR:AFTERQC] Failed to trim reads")
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
				run_command(inputs, f"{inputs.skewer} -m pe -q 10 -Q 20 -l 35 -n -t {inputs.threads} -o {trimmed_Name} {r._r1} {r._r2}", f"[ERROR:SKEWER] Failed to trim reads {r._r1} {r._r2}")
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
				run_command(inputs, f"/bin/bash -c '{inputs.bwa} mem -t {inputs.threads} -M -R \"@RG\\tID:{r._id}_{r._lib}_{r._lane}\\tSM:{r._id}\\tPL:{str.upper(r._plat)}\\tLB:{r._lib}\" {inputs.reference} <(gunzip -c {r._trimmedR1}) <(gunzip -c {r._trimmedR2}) | {inputs.samtools} sort -@{inputs.threads} -l9 -o {aligned_Name}'", f"[ERROR:BWA] Failed to align reads {r._trimmedR1} {r._trimmedR2}")
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
			run_command(inputs, f"{inputs.star} --runThreadN {inputs.threads} --genomeDir {runDir} --readFilesIn {readFilesIn} --outTmpDir {tempDir} --readFilesCommand zcat --twopassMode Basic --quantMode TranscriptomeSAM --outSAMtype BAM Unsorted --outBAMcompression 10 --outFileNamePrefix {sampleName}", error_msg)
	return True

def ngs_mark_sort(inputs, targets_data):
	if str.upper(inputs.workflow) != 'RNASEQ':
		for r in targets_data:
			rmdup_Name = os.path.abspath(inputs.dir + '/bamProcessed/Rmdup_' + r._id + '_' + r._lib + '_' + r._lane + '.bam')
			
			if(os.path.exists(rmdup_Name)):
				print_log(inputs, "[INFO] Processed bam files found, skipping.")
			else:
				## Remove Duplicates ##
				run_command(inputs, f"{inputs.sambamba} markdup --overflow-list-size=800000 --io-buffer-size=512 -l9 -p -t {inputs.threads} {r._aligned} {rmdup_Name}", "[ERROR:SAMBAMBA] Failed to mark&remove duplicates in .bam file")
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
				print_log(inputs, f"[INFO] Single bam file for sample {s} skipping.")
				merged_tumor_bam_name = sample_ids[s]["TUMOR"][0]
			else:
				merged_tumor_bam_name = os.path.abspath(inputs.dir + '/bamMerged/Merged_Tumor_' + s + '.bam')
				if(os.path.exists(merged_tumor_bam_name)):
					print_log(inputs, "[INFO] Merged .bam file found for tumor sample, skipping.")
				else:
					run_command(inputs, f"{inputs.sambamba} merge -l9 -t {inputs.threads} {merged_tumor_bam_name} {' '.join(sample_ids[s]['TUMOR'])}", f"[ERROR:SAMBAMBA] Failed to merge files {' '.join(sample_ids[s]['TUMOR'])}")

			## NORMAL ##
			if(len(sample_ids[s]["NORMAL"]) == 1):
				print_log(inputs, f"[INFO] Single bam file for sample {s} skipping.")
				merged_normal_bam_name = sample_ids[s]["NORMAL"][0]
			else:
				merged_normal_bam_name = os.path.abspath(inputs.dir + '/bamMerged/Merged_Normal_' + s + '.bam')
				if(os.path.exists(merged_normal_bam_name)):
					print_log(inputs, "[INFO] Merged .bam file found for normal sample, skipping")
				else:
					run_command(inputs, f"{inputs.sambamba} merge -l9 -t {inputs.threads} {merged_normal_bam_name} {' '.join(sample_ids[s]['NORMAL'])}", f"[ERROR:SAMBAMA] Failed to merge files {' '.join(sample_ids[s]['TUMOR'])}")
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
				print_log(inputs, f"[INFO] Single bam file for sample {s} skipping.")
				merged_bam_name = sample_ids[s][0]
			else:
				merged_bam_name = os.path.abspath(inputs.dir + '/bamMerged/Merged_' + s + '.bam')
				if(os.path.exists(merged_bam_name)):
					print_log(inputs, "[INFO] Merged bam files found, skipping.")
				else:
					run_command(inputs, f"{inputs.sambamba} merge -l9 -t {inputs.threads} {merged_bam_name} {' '.join(sample_ids[s])}", f"[ERROR:SAMBAMBA] Failed to merge files {' '.join(sample_ids[s])}")
			targets_data_processed.append(ngs_classes.WxsSample(_id = s,_bam = merged_bam_name))
	
	return 0

def bqsr(inputs, targets_data_processed):
	
	for s in targets_data_processed:
		if str.upper(inputs.workflow) == "SOMATICSNVINDEL":
			if len(glob.glob(os.path.abspath(f"{s._bamTumor}.bai"))) > 0:
				print_log(inputs, "[INFO] Indices found, skipping.")
			else:
				run_command(inputs, f"{inputs.sambamba} index -t {inputs.threads} {s._bamTumor}")
			if len(glob.glob(os.path.abspath(f"{s._bamNormal}.bai"))) > 0:
				print_log(inputs, "[INFO] Indices found, skipping.")
			else:
				run_command(inputs, f"{inputs.sambamba} index -t {inputs.threads} {s._bamNormal}")
		elif str.upper(inputs.workflow) == "GERMLINESNVINDEL":
			if len(glob.glob(os.path.abspath(f"{s._bam}.bai"))) > 0:
				print_log(inputs, "[INFO] Indices found, skipping.")
			else:
				run_command(inputs, f"{inputs.sambamba} index -t {inputs.threads} {s._bam}")
			
	## BQSR (Step 1) ##
	for s in targets_data_processed:
		if str.upper(inputs.workflow) == "SOMATICSNVINDEL":
			tumorRecalibTable = os.path.abspath(f"{inputs.dir}/recalib/RecalibrationTableTumor_{s._id}.table")
			normalRecalibTable = os.path.abspath(f"{inputs.dir}/recalib/RecalibrationTableNormal_{s._id}.table")
			if os.path.exists(tumorRecalibTable) and os.path.exists(normalRecalibTable):
				print_log(inputs, "[INFO] Recalibration tables found, skipping.")
			else:
				recalib_commands = []
				if inputs.bed:
					recalib_commands.append(f"{inputs.gatk} -T BaseRecalibrator -nct {inputs.threads} -R {inputs.reference} -I {s._bamTumor} -knownSites {' -knownSites '.join(inputs.knownsites)} -o {tumorRecalibTable} -L {inputs.bed}")
				else:
					recalib_commands.append(f"{inputs.gatk} -T BaseRecalibrator -nct {inputs.threads} -R {inputs.reference} -I {s._bamTumor} -knownSites {' -knownSites '.join(inputs.knownsites)} -o {tumorRecalibTable}")
				recalib_commands.append(f"{inputs.gatk} -T BaseRecalibrator -nct {inputs.threads} -R {inputs.reference} -I {s._bamNormal} -knownSites {' -knownSites '.join(inputs.knownsites)} -o {normalRecalibTable}")
				
				for cmd in recalib_commands:
					run_command(inputs, cmd)
				
				s._recalibTableTumor = tumorRecalibTable
				s._recalibTableNormal = normalRecalibTable
		elif str.upper(inputs.workflow) == "GERMLINESNVINDEL":
			RecalibTable = os.path.abspath(f"{inputs.dir}/recalib/RecalibrationTableNormal_{s._id}.table")
			if os.path.exists(RecalibTable):
				print_log(inputs, "[INFO] Recalibration table found, skipping.")
			else:
				recalib_command = ""
				if inputs.bed:
					recalib_command = f"{inputs.gatk} -T BaseRecalibrator -nct {inputs.threads} -R {inputs.reference} -I {s._bam} -knownSites {' -knownSites '.join(inputs.knownsites)} -o {RecalibTable} -L {inputs.bed}"
				else:
					recalib_command = f"{inputs.gatk} -T BaseRecalibrator -nct {inputs.threads} -R {inputs.reference} -I {s._bam} -knownSites {' -knownSites '.join(inputs.knownsites)} -o {RecalibTable}"
				
				run_command(inputs, recalib_command)
				
				s._recalibTable = RecalibTable
				
	## BQSR (Step 2) ##
	for s in targets_data_processed:
		if str.upper(inputs.workflow) == "SOMATICSNVINDEL":
			tumorCalibBam = os.path.abspath(f"{inputs.dir}/bamCalib/CalibratedTumor_{s._id}.bam")
			normalCalibBam = os.path.abspath(f"{inputs.dir}/bamCalib/CalibratedNormal_{s._id}.bam")
			
			if os.path.exists(tumorCalibBam) and os.path.exists(normalCalibBam):
				print_log(inputs, "[INFO] BQSR adjusted files found, skipping.")
			else:
				if inputs.bed:
					run_command(inputs, f"{inputs.gatk} -T PrintReads -nct {inputs.threads} -R {inputs.reference} -I {s._bamTumor} -BQSR {s._recalibTableTumor} -o {tumorCalibBam} -L {inputs.bed}")
				else:
					run_command(inputs, f"{inputs.gatk} -T PrintReads -nct {inputs.threads} -R {inputs.reference} -I {s._bamTumor} -BQSR {s._recalibTableTumor} -o {tumorCalibBam}")
				
				if inputs.bed:
					run_command(inputs, f"{inputs.gatk} -T PrintReads -nct {inputs.threads} -R {inputs.reference} -I {s._bamNormal} -BQSR {s._recalibTableNormal} -o {normalCalibBam} -L {inputs.bed}")
				else:
					run_command(inputs, f"{inputs.gatk} -T PrintReads -nct {inputs.threads} -R {inputs.reference} -I {s._bamNormal} -BQSR {s._recalibTableNormal} -o {normalCalibBam}")
			
			s._bamCalibTumor = tumorCalibBam
			s._bamCalibNormal = normalCalibBam
		elif str.upper(inputs.workflow) == "GERMLINESNVINDEL":
			CalibBam = os.path.abspath(f"{inputs.dir}/bamCalib/CalibratedNormal_{s._id}.bam")
			
			if os.path.exists(CalibBam):
				print_log(inputs, "[INFO] BQSR adjusted files found, skipping.")
			else:
				if inputs.bed:
					run_command(inputs, f"{inputs.gatk} -T PrintReads -nct {inputs.threads} -R {inputs.reference} -I {s._bam} -BQSR {s._recalibTable} -o {CalibBam} -L {inputs.bed}")
				else:
					run_command(inputs, f"{inputs.gatk} -T PrintReads -nct {inputs.threads} -R {inputs.reference} -I {s._bam} -BQSR {s._recalibTable} -o {CalibBam}")
			
			s._bamCalib = CalibBam

	return True
