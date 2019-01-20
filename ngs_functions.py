#!/usr/bin/python

######################################################
import sys, os, argparse
import os.path, re
import glob

dir_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(dir_path)
import ngs_classes
######################################################


################## FUNCTIONS #########################
######################################################
def check_inputs(inputs):
	_FLAG = True

	'''
	print "\n---------------------------------------------------"
	print "--------------------- INFO ------------------------"
	for key, val in vars(inputs).iteritems():
		print "%s: %s" % (key, val)

		for arg, val in vars(args).iteritems():
			if isinstance(val, str):
				if not(str.upper(val) == "GERMLINESNVINDEL" or
						str.upper(val) == "SOMATICSNVINDEL" or
						str.upper(val) == "EXPRESSIONQUANTIFICATION"):
				if not (os.path.exists(os.path.abspath(val))):
					print "  Check Argument: [%s --> %s]" % (arg, val)
					_FLAG = False
					'''
					return _FLAG


def read_config(args):
    ## Parse Config ##
    local_data = {}
    fh = open(args.config, "rU")
    try:
        for l in fh:
            ss = l.rstrip("\n")
            if(re.match('^#', ss)):
                try:
                    next(fh)
                except StopIteration:
                    break
            match = re.search('^([^\s]+):\t([^\s]+)', ss)
            if(match):
                local_data[str.upper(match.group(1))] = match.group(2)
    finally:
        fh.close()
    
    if(len(local_data.keys()) == 0):
    	print "!! Error in reading the config file"
    	sys.exit(1)

    if args.verbose:
        print "\n"
        print "\/" * 40
        for key, val in local_data.items():
            print "\/ %s: %s" % (key, val)
            
        print "\/" * 40
        print "\n"
        
    inputs = ngs_classes.ngsInputs()
    
    inputs.workflow = args.workflow
    inputs.dir = args.dir
    inputs.in_file = args.in_file
    inputs.dry = args.dry
    inputs.verbose = args.verbose
    inputs.clear = args.clear
    inputs.threads = args.threads
    inputs.title = args.title
    
    if(str.lower(args.workflow_type) == 'germlinesnvindel'):
        inputs.ped = args.ped
    
    inputs.reference = local_data["REFERENCE"]
    inputs.gatk = local_data["GATK"]
    inputs.picard = local_data["PICARD"]
    inputs.samtools = local_data["SAMTOOLS"]
    inputs.bcftools = local_data["BCFTOOLS"]
    inputs.bowtie2 = local_data["BOWTIE2"]
    inputs.bowtie2build = local_data["BOWTIE2BUILD"]
    inputs.bwa = local_data["BWA"]
    inputs.cutadapt = local_data["CUTADAPT"]
    inputs.afterqc = local_data["AFTERQC"]
    inputs.skewer = local_data["SKEWER"]
    inputs.sambamba = local_data["SAMBAMBA"]
    inputs.strelkabin = local_data["STRELKABIN"]
    inputs.mantabin = local_data["MANTABIN"]
    
    
    if(re.search('bowtie2', str.lower(args.aligner))):
        inputs.aligner = "BOWTIE2"
        inputs.indexer = "BOWTIE2BUILD"
    elif(re.search('bwa', str.lower(args.aligner))):
        inputs.aligner = "BWA"
        inputs.indexer = "BWABUILD"
    
    if not(args.trimmer is None):
        if(re.search('cutadapt', str.lower(args.trimmer))):
            inputs.trimmer = "CUTADAPT"
        elif(re.search('afterqc', str.lower(args.trimmer))):
            inputs.trimmer = "AFTERQC"
        elif(re.search('skewer', str.lower(args.trimmer))):
            inputs.trimmer = "SKEWER"
    else:
        inputs.trimmer = None
    
    if not(args.bed is None):
        inputs.bed = args.bed
    else:
        inputs.bed = None
    
    if(re.search('haplotypecaller', str.lower(args.tool))):
        inputs.tool = "HAPLOTYPECALLER"
    elif(re.search('bcftools', str.lower(args.tool))):
        inputs.tool = "BCFTOOLS"
    elif(re.search('varscan2', str.lower(args.tool))):
        inputs.tool = "VARSCAN2"
    elif(re.search('mutect1', str.lower(args.tool))):
        inputs.tool = "MUTECT1"
    elif(re.search('rsem', str.lower(args.tool))):
        inputs.tool = "RSEM"
    elif(re.search('featurecounts', str.lower(args.tool))):
        inputs.tool = "FEATURECOUNTS"
    

    return inputs

def get_inputs():
    
    parser = argparse.ArgumentParser(prog = 'NGS',
                                     description = '*** Epigenetiks WXS Pipeline ***')   
    subparsers = parser.add_subparsers(description = 'Workflow to process',
                                       dest = 'workflow_type')

    ## Germline ##
    parser_germline = subparsers.add_parser('GermlineSNVIndel',
                                            help = 'Germline SNV & Indel Detection')
    parser_germline.add_argument('--config',
                                 default = None,
                                 required = True,
                                 help = 'Path to config file')  
    parser_germline.add_argument('--tool',
                                 default = 'Strelka2',
                                 required = True,
                                 help = 'Name of the caller <HaplotypeCaller|Strelka2>')
    parser_germline.add_argument('--trimmer',
                                default = None,
                                required = False,
                                help = 'Name of the trimmer <cutadapt|afterqc|skewer>')
    parser_germline.add_argument('--aligner',
                                default = 'bwa',
                                required = True,
                                help = 'Name of the aligner <bowtie2|bwa>')
    parser_germline.add_argument('--dir',
                                 default = 'NGS_WD',
                                 required = True,
                                 help = 'Working Directory')
    parser_germline.add_argument('--in_file',
                                default = None,
                                required = True,
                                help = 'Targets file containing sample ids and associated fastq files')
    parser_germline.add_argument('--verbose',
                                action = "store_true",
                                help = 'Verbosity')
    parser_germline.add_argument('--dry',
                                action = "store_true",
                                help = 'Dry-Run')
    parser_germline.add_argument('--clear',
                                 action = 'store_true',
                                 default = False,
                                 help = 'Clear workspace')
    parser_germline.add_argument('--threads',
                                 default = 2,
                                 type = int,
                                 help = "Number of threads to use")
    parser_germline.add_argument('--title',
                                 default = 'NGS_Analysis',
                                 help = 'Title of the run')
    parser_germline.add_argument('--ped',
                                 default = None,
                                 help = 'Ped file for family analysis')
    parser_germline.add_argument('--bed',
                                 default = None,
                                 help = 'Bed file for regions')
    parser_germline.add_argument('--exome',
                                 action = 'store_true',
                                 default = False,
                                 help = 'Whole Exome Data')
    
    ## Somatic ##
    parser_somatic = subparsers.add_parser('SomaticSNVIndel',
                                            help = 'Somatic SNV & Indel Detection')
    parser_somatic.add_argument('--config',
                                 default = None,
                                 required = True,
                                 help = 'Path to config file')
    parser_somatic.add_argument('--tool',
                                default = 'Strelka2',
                                required = True,
                                help = 'Name of the caller <MuTect1|Strelka2>')
    parser_somatic.add_argument('--trimmer',
                                default = None,
                                required = False,
                                help = 'Name of the trimmer <cutadapt|afterqc|skewer> ')
    parser_somatic.add_argument('--aligner',
                                default = 'bwa',
                                required = True,
                                help = 'Name of the aligner <bowtie2|bwa>')
    parser_somatic.add_argument('--in_file',
                                default = None,
                                required = True,
                                help = 'Targets file containing sample ids and associated fastq files')
    parser_somatic.add_argument('--verbose',
                                action = "store_true",
                                help = 'Verbosity')
    parser_somatic.add_argument('--dir',
                                default = 'NGS_WD',
                                required = True,
                                help = 'Working directory')
    parser_somatic.add_argument('--dry',
                                action = "store_true",
                                help = 'Dry-Run')
    parser_somatic.add_argument('--clear',
                                action = 'store_true',
                                default = False,
                                help = 'Clear workspace')
    parser_somatic.add_argument('--threads',
                                default = 2,
                                type = int,
                                help = 'Number of threads to use')
    parser_somatic.add_argument('--title',
                                default = 'Sample_Analysis',
                                help = 'Title of the run')
    parser_somatic.add_argument('--bed',
                                default = None,
                                help = 'Bed file for regions')
    parser_somatic.add_argument('--exome',
                                action = 'store_true',
                                default = False,
                                help = 'Whole Exome Data')
    
    ## Expression ##
    parser_expr = subparsers.add_parser('ExpressionQuantification',
                                        help = 'Expression Quantification')
    parser_expr.add_argument('--config',
                             default = None,
                             required = True,
                             help = 'Path to config file')
    parser_expr.add_argument('--tool',
                             default = None,
                             required = True,
                             help = 'Name of the quantifier <RSEM|featureCounts>')
    parser_expr.add_argument('--aligner',
                             default = None,
                             required = True,
                             help = 'Name of the aligner <bowtie2|bwa>|STAR')
    parser_expr.add_argument('--in_file',
                            default = None,
                            required = True,
                            help = 'Targets file containing sample ids and associated fastq files')
    parser_expr.add_argument('--dir',
                             default = None,
                             required = True,
                             help = 'Working directory')
    parser_expr.add_argument('--trimmer',
                            default = None,
                            required = False,
                            help = 'Name of the trimmer tool <cutadapt|afterqc|skewer>')
    parser_expr.add_argument('--verbose',
                            action = "store_true",
                            help = 'Verbosity')
    parser_expr.add_argument('--dry',
                            action = "store_true",
                            help = 'Dry-Run')
    parser_expr.add_argument('--clear',
                             action = 'store_true',
                             default = False,
                             help = 'Clear workspace')
    parser_expr.add_argument('--threads',
                            default = 2,
                            type = int,
                            help = 'Number of threads to use')
    parser_expr.add_argument('--title',
                             default = 'Sample_Analysis',
                             help = 'Title of the run')
    
    _res = read_config(parser.parse_args())
    return _res

def read_targets(args):
    
    in_data = []
    fh = open(args.in_file, "rU")
    try:
        indx = {}    
        lCount = 0
        for l in fh:
            lCount += 1
            ss = l.rstrip("\n")
            if(re.match('^#', ss)):
                try:
                    next(fh)
                except StopIteration:
                    break
            local = ss.split('\t')
            if lCount == 1:
                for i in range(len(local)):
                    if re.search('read1', str.lower(local[i])):
                        indx["Read1"] = i
                    elif re.search('read2', str.lower(local[i])):
                        indx["Read2"] = i
                    elif re.search('lane', str.lower(local[i])):
                        indx["Lane"] = i
                    elif re.search('id', str.lower(local[i])):
                        indx["ID"] = i
                    elif re.search('plat', str.lower(local[i])):
                        indx["Platform"] = i
                    elif re.search('lib', str.lower(local[i])):
                        indx["Library"] = i
                    elif re.search('type', str.lower(local[i])):
                        indx["Type"] = i
            else:
                if len(list(filter(lambda x: re.search('^[^\s]+$', x),local))) != len(indx):
                    print '!Warning, format issue. Line %d skipping' % lCount
                    next(fh)
                else:
                    if(str.upper(args.workflow_type) == "GERMLINESNVINDEL"):
                        in_data.append(ngs_classes.WxsGermline(_r1 = local[indx["Read1"]],
                                                               _r2 = local[indx["Read2"]],
                                                               _lane = local[indx["Lane"]],
                                                               _id = local[indx["ID"]],
                                                               _plat = local[indx["Platform"]],
                                                               _lib = local[indx["Library"]]))
                    elif(str.upper(args.workflow_type) == "SOMATICSNVINDEL"):
                        in_data.append(ngs_classes.WxsSomatic(_r1 = local[indx["Read1"]],
                                                              _r2 = local[indx["Read2"]],
                                                              _lane = local[indx["Lane"]],
                                                              _id = local[indx["ID"]],
                                                              _plat = local[indx["Platform"]],
                                                              _lib = local[indx["Library"]],
                                                              _type = local[indx["Type"]]))
                    elif(str.upper(args.workflow_type) == "EXPRESSIONQUANTIFICATION"):
                        in_data.append(ngs_classes.Expression(_r1 = local[indx["Read1"]],
                                                              _r2 = local[indx["Read2"]],
                                                              _id = local[indx["ID"]]))
    finally:
        fh.close()
    
    return in_data


def prep_wd(path_wd, workflow_type, inputs):
	_flag = True
	if(str.upper(workflow_type) == "GERMLINESNVINDEL" or str.upper(workflow_type) == "SOMATICSNVINDEL"):
		for _dir in ["readsTrimmed", "readsAligned", "bamProcessed","bamMerged", "bamCalib", "metrics", "indices", "recalib","rawSNP", "filtSNP"]:
			local_path = os.path.abspath(path_wd + '/' + _dir)
	if not os.path.exists(local_path):
		try:
		if inputs.verbose:
		print "Creating Directory: %s" % (local_path)
		if not inputs.dry:
		os.mkdir(local_path)
		except:
			_flag = False
			break
			else:
			if inputs.clear:
			if inputs.verbose:
			print "Deleting Directory: %s" % (local_path)
			if not inputs.dry:
			os.system("rm -viR %s" % (local_path))
			else:
				if not inputs.dry:
				os.system("rm -R %s" % (local_path))

				if not os.path.exists(local_path):
					try:
					if inputs.verbose:
					print "Creating Directory: %s" % (local_path)
					if not inputs.dry:
					os.mkdir(local_path)
					except:
						_flag = False
						break

						elif(str.upper(workflow_type) == "EXPRESSIONQUANTIFICATION"):
						for _dir in ["readsTrimmed", "readsAligned"]:
						local_path = os.path.abspath(path_wd + '/' + _dir)
						if not os.path.exists(local_path):
							try:
	print "Creating Directory: %s" % (local_path)
	      os.mkdir(local_path)
	      except:
		      _flag = False
		      break

		      return _flag

def ngs_index(inputs):
    if re.search('bowtie2build', str.lower(inputs.indexer)):
        if(len(glob.glob(os.path.abspath(inputs.dir + '/indices/ref_indices_bowtie2*'))) > 0):
            print "    *** indices found, skipping.."
        else:
            if inputs.verbose:
                print "Running: %s --threads %d %s %s" % (inputs.bowtie2build,
                                                          inputs.threads,
                                                          inputs.reference,
                                                          os.path.abspath(inputs.dir + '/indices/ref_indices_bowtie2'))
            if not inputs.dry:
                os.system("%s --threads %d %s %s" % (inputs.bowtie2build,
                                                     inputs.threads,
                                                     inputs.reference,
                                                     os.path.abspath(inputs.dir + '/indices/ref_indices_bowtie2')))
    elif re.search("bwa", str.lower(inputs.indexer)):
        if(len(glob.glob(os.path.abspath(inputs.reference + '.bwt'))) > 0):
            print "    *** indices found, skipping.."
        else:
            if inputs.verbose:
                print "Running: %s index %s" % (inputs.bwa,
                                                inputs.reference)
            if not inputs.dry:
                os.system("%s index %s" % (inputs.bwa,
                                           inputs.reference))
        
    return True

def ngs_trim(inputs, targets_data):
    if inputs.trimmer:        
        if re.search('cutadapt', str.lower(inputs.trimmer)):
            for r in targets_data:
                if(os.path.exists(os.path.abspath(inputs.dir + '/readsTrimmed/Trimmed_Read1_' + r._id + '_' + r._lib + '_' + r._lane + '.fastq')) and
                   os.path.exists(os.path.abspath(inputs.dir + '/readsTrimmed/Trimmed_Read2_' + r._id + '_' + r._lib + '_' + r._lane + '.fastq'))):
                    print "    *** trimmed reads found, skipping..."
                    r._trimmedR1 = os.path.abspath(inputs.dir + '/readsTrimmed/Trimmed_Read1_' + r._id + '_' + r._lib + '_' + r._lane + '.fastq')
                    r._trimmedR2 = os.path.abspath(inputs.dir + '/readsTrimmed/Trimmed_Read2_' + r._id + '_' + r._lib + '_' + r._lane + '.fastq')
                else:
                    if inputs.verbose:
                        print "Running: %s -q 15,10 -A XXX -m 30 -o %s -p %s %s %s" % (inputs.cutadapt,
                                                                                 os.path.abspath(inputs.dir + '/readsTrimmed/Trimmed_Read1_' + r._id + '_' + r._lib + '_' + r._lane + '.fastq'),
                                                                                 os.path.abspath(inputs.dir + '/readsTrimmed/Trimmed_Read2_' + r._id + '_' + r._lib + '_' + r._lane + '.fastq'),
                                                                                 r._r1,
                                                                                 r._r2)
                    if not(inputs.dry):
                        os.system("%s -q 15,10 -A XXX -m 30 -o %s -p %s %s %s" % (inputs.cutadapt,
                                                                            os.path.abspath(inputs.dir + '/readsTrimmed/Trimmed_Read1_' + r._id + '_' + r._lib + '_' + r._lane + '.fastq'),
                                                                            os.path.abspath(inputs.dir + '/readsTrimmed/Trimmed_Read2_' + r._id + '_' + r._lib + '_' + r._lane + '.fastq'),
                                                                            r._r1,
                                                                            r._r2))
                    r._trimmedR1 = os.path.abspath(inputs.dir + '/readsTrimmed/Trimmed_Read1_' + r._id + '_' + r._lib + '_' + r._lane + '.fastq')
                    r._trimmedR2 = os.path.abspath(inputs.dir + '/readsTrimmed/Trimmed_Read2_' + r._id + '_' + r._lib + '_' + r._lane + '.fastq')
        
        elif(re.search('afterqc', str.lower(inputs.trimmer))):
            for r in targets_data:
                stagr1 = os.path.basename(r_r1).split(".")[0]
                stagr1 = os.path.basename(r_r2).split(".")[0]
                if(os.path.exists(os.path.abspath(inputs.dir + '/readsTrimmed/' + stagr1 + '.good.fq.gz')) and
                   os.path.exists(os.path.abspath(inputs.dir + '/readsTrimmed/' + stagr2 + '.good.fq.gz'))):
                    print "    *** trimmed reads found, skipping..."
                    r._trimmedR1 = os.path.abspath(inputs.dir + '/readsTrimmed/' + stagr1 + '.good.fq.gz')
                    r._trimmedR2 = os.path.abspath(inputs.dir + '/readsTrimmed/' + stagr2 + '.good.fq.gz')
                else:
                    if inputs.verbose:
                        print "Running: python %s -1 %s -2 %s --debubble False --gzip --compression 9 -g %s --no_correction" % (inputs.afterqc,
                                                                                                                               r._r1, r._r2,
                                                                                                                               os.path.abspath(inputs.dir + '/readsTrimmed/'))
                    if not inputs.dry:
                        os.system("python %s -1 %s -2 %s --debubble False --gzip --compression 9 -g %s --no_correction" % (inputs.afterqc,
                                                                                                                           r._r1, r._r2,
                                                                                                                           os.path.abspath(inputs.dir + '/readsTrimmed/')))
                    r._trimmedR1 = os.path.abspath(inputs.dir + '/readsTrimmed/' + stagr1 + '.good.fq.gz')
                    r._trimmedR2 = os.path.abspath(inputs.dir + '/readsTrimmed/' + stagr2 + '.good.fq.gz')
        
        elif(re.search('skewer', str.lower(inputs.trimmer))):
            for r in targets_data:
                if(os.path.exists(os.path.abspath(inputs.dir + '/readsTrimmed/' + r._id + '_' + r._lib + '_' + r._lane + '-trimmed-pair1.fastq')) and
                   os.path.exists(os.path.abspath(inputs.dir + '/readsTrimmed/' + r._id + '_' + r._lib + '_' + r._lane + '-trimmed-pair2.fastq'))):
                    print "    *** trimmed reads found, skipping..."
                    r._trimmedR1 = os.path.abspath(inputs.dir + '/readsTrimmed/' + r._id + '_' + r._lib + '_' + r._lane + '-trimmed-pair1.fastq')
                    r._trimmedR2 = os.path.abspath(inputs.dir + '/readsTrimmed/' + r._id + '_' + r._lib + '_' + r._lane + '-trimmed-pair2.fastq')
                else:
                    if inputs.verbose:
                        print "Running: %s -m pe -q 10 -Q 20 -l 35 -n -t %d -o %s %s %s" % (inputs.skewer,
                                                                                            inputs.threads,
                                                                                            os.path.abspath(inputs.dir + '/readsTrimmed/' + r._id + '_' + r._lib + '_' + r._lane),
                                                                                            r._r1, r._r2)
                    if not inputs.dry:
                        os.system("%s -m pe -q 10 -Q 20 -l 35 -n -t %d -o %s %s %s" % (inputs.skewer,
                                                                                       inputs.threads,
                                                                                       os.path.abspath(inputs.dir + '/readsTrimmed/' + r._id + '_' + r._lib + '_' + r._lane),
                                                                                       r._r1, r._r2))
                    r._trimmedR1 = os.path.abspath(inputs.dir + '/readsTrimmed/' + r._id + '_' + r._lib + '_' + r._lane + '-trimmed-pair1.fastq')
                    r._trimmedR2 = os.path.abspath(inputs.dir + '/readsTrimmed/' + r._id + '_' + r._lib + '_' + r._lane + '-trimmed-pair2.fastq')
        
    return True

def ngs_align(inputs, targets_data):
    if re.search('bowtie2', str.lower(inputs.aligner)):
        for r in targets_data:
            if(os.path.exists(os.path.abspath(inputs.dir + '/readsAligned/Aligned_' + r._id + '_' + r._lib + '_' + r._lane + '_bowtie2.sam'))):
                print "    *** aligned reads found, skipping..."
                r._aligned = os.path.abspath(inputs.dir + '/readsAligned/Aligned_' + r._id + '_' + r._lib + '_' + r._lane + '_bowtie2.sam')
            else:
                if inputs.verbose:
                    print "Running: %s --threads %d -x %s --quiet --phred33 --no-mixed --no-unal "\
                    "--very-sensitive-local -1 %s -2 %s -S %s --rg-id %s --rg %s --rg %s --rg %s" % (inputs.bowtie2,
                                                                                                inputs.threads,
                                                                                                os.path.abspath(inputs.dir + '/indices/ref_indices_bowtie2'),
                                                                                                r._trimmedR1,r._trimmedR2,
                                                                                                os.path.abspath(inputs.dir + '/readsAligned/Aligned_' +
                                                                                                                r._id + '_' + r._lib + '_' + r._lane + '_bowtie2.sam'),
                                                                                                r._id + '_' + r._lib + '_' + r._lane,
                                                                                                'SM:' + r._id,
                                                                                                'PL:' + str.upper(r._plat),
                                                                                                'LB:' + r._lib)
                if not(inputs.dry):
                    os.system("%s --threads %d -x %s --quiet --phred33 --no-mixed --no-unal "\
                              "--very-sensitive-local -1 %s -2 %s -S %s --rg-id %s --rg %s --rg %s --rg %s" % (inputs.bowtie2,
                                                                                                          inputs.threads,
                                                                                                          os.path.abspath(inputs.dir + '/indices/ref_indices'),
                                                                                                          r._trimmedR1,r._trimmedR2,
                                                                                                          os.path.abspath(inputs.dir + '/readsAligned/Aligned_' +
                                                                                                                          r._id + '_' + r._lib + '_' + r._lane + '_bowtie2.sam'),
                                                                                                          r._id + '_' + r._lib + '_' + r._lane,
                                                                                                          'SM:' + r._id,
                                                                                                          'PL:' + str.upper(r._plat),
                                                                                                          'LB:' + r._lib))
                r._aligned = os.path.abspath(inputs.dir + '/readsAligned/Aligned_' + r._id + '_' + r._lib + '_' + r._lane + '_bowtie2.sam')
    elif re.search('bwa', str.lower(inputs.aligner)):
        for r in targets_data:
            if(os.path.exists(os.path.abspath(inputs.dir + '/readsAligned/Aligned_' + r._id + '_' + r._lib + '_' + r._lane + '_bwa.sam'))):
                print "    *** aligned reads found, skipping..."
                r._aligned = os.path.abspath(inputs.dir + '/readsAligned/Aligned_' + r._id + '_' + r._lib + '_' + r._lane + '_bwa.sam')
            else:
                if inputs.verbose:
                    print "Running: %s mem -t %d -M -R \'@RG\\tID:%s\\tSM:%s\\tPL:%s\\tLB:%s\' %s %s %s > %s" % (inputs.bwa,
                                                                                                    inputs.threads,
                                                                                                    r._id + '_' + r._lib + '_' + r._lane,
                                                                                                    r._id,
                                                                                                    str.upper(r._plat),
                                                                                                    r._lib,
                                                                                                    inputs.reference,
                                                                                                    r._trimmedR1,
                                                                                                    r._trimmedR2,
                                                                                                    os.path.abspath(inputs.dir + '/readsAligned/Aligned_' +
                                                                                                                    r._id + '_' + r._lib + '_' + r._lane + '_bwa.sam'))
                    
                if not(inputs.dry):
                    os.system("%s mem -t %d -M -R \'@RG\\tID:%s\\tSM:%s\\tPL:%s\\tLB:%s\' %s %s %s > %s" % (inputs.bwa,
                                                                                                    inputs.threads,
                                                                                                    r._id + '_' + r._lib + '_' + r._lane,
                                                                                                    r._id,
                                                                                                    str.upper(r._plat),
                                                                                                    r._lib,
                                                                                                    inputs.reference,
                                                                                                    r._trimmedR1,
                                                                                                    r._trimmedR2,
                                                                                                    os.path.abspath(inputs.dir + '/readsAligned/Aligned_' +
                                                                                                                    r._id + '_' + r._lib + '_' + r._lane + '_bwa.sam')))
                
                r._aligned = os.path.abspath(inputs.dir + '/readsAligned/Aligned_' + r._id + '_' + r._lib + '_' + r._lane + '_bwa.sam')    

    return True


def ngs_mark_sort(inputs, targets_data):
    for r in targets_data:
        if(os.path.exists(os.path.abspath(inputs.dir + '/bamProcessed/Processed_' + r._id + '_' + r._lib + '_' + r._lane + '.bam'))):
            print "    *** marked files found, skipping..."
            r._processed = os.path.abspath(inputs.dir + '/bamProcessed/Processed_' + r._id + '_' + r._lib + '_' + r._lane + '.bam')
        else:
            if inputs.verbose:
                print "Running: %s view -bh -@ %d -o %s %s" % (inputs.samtools,
                                                            inputs.threads,
                                                            os.path.abspath(inputs.dir + '/bamProcessed/' + r._id + '_' + r._lib + '_' + r._lane + '.bam'),
                                                            r._aligned)
            if not(inputs.dry):
                os.system("%s view -bh -@ %d -o %s %s" % (inputs.samtools,
                                                       inputs.threads,
                                                       os.path.abspath(inputs.dir + '/bamProcessed/' + r._id + '_' + r._lib + '_' + r._lane + '.bam'),
                                                       r._aligned))
            if inputs.verbose:
                print "Running: %s markdup -r -p -t %d %s %s" % (inputs.sambamba,
                                                                 inputs.threads,
                                                                 os.path.abspath(inputs.dir + '/bamProcessed/' + r._id + '_' + r._lib + '_' + r._lane + '.bam'),
                                                                 os.path.abspath(inputs.dir + '/bamProcessed/Rmdup_' + r._id + '_' + r._lib + '_' + r._lane + '.bam'))
            if not(inputs.dry):
                os.system("%s markdup -r -p -t %d %s %s" % (inputs.sambamba,
                                                            inputs.threads,
                                                            os.path.abspath(inputs.dir + '/bamProcessed/' + r._id + '_' + r._lib + '_' + r._lane + '.bam'),
                                                            os.path.abspath(inputs.dir + '/bamProcessed/Rmdup_' + r._id + '_' + r._lib + '_' + r._lane + '.bam')))
            if inputs.verbose:
                print "Running: %s sort -@ %d -o %s %s" % (inputs.samtools,
                                                           inputs.threads,
                                                           os.path.abspath(inputs.dir + '/bamProcessed/Processed_' + r._id + '_' + r._lib + '_' + r._lane + '.bam'),
                                                           os.path.abspath(inputs.dir + '/bamProcessed/Rmdup_' + r._id + '_' + r._lib + '_' + r._lane + '.bam'))
            if not(inputs.dry):
                os.system("%s sort -@ %d -o %s %s" % (inputs.samtools,
                                                      inputs.threads,
                                                      os.path.abspath(inputs.dir + '/bamProcessed/Processed_' + r._id + '_' + r._lib + '_' + r._lane + '.bam'),
                                                      os.path.abspath(inputs.dir + '/bamProcessed/Rmdup_' + r._id + '_' + r._lib + '_' + r._lane + '.bam')))
            
            r._processed = os.path.abspath(inputs.dir + '/bamProcessed/Processed_' + r._id + '_' + r._lib + '_' + r._lane + '.bam')
    return True
    


def ngs_merge(inputs, targets_data, targets_data_processed):
    
    if str.upper(inputs.workflow_type) == 'SOMATICSNVINDEL':
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
            if(os.path.exists(os.path.abspath(inputs.dir + '/bamMerged/Merged_Tumor_' + s + '.bam')) and
               os.path.exists(os.path.abspath(inputs.dir + '/bamMerged/Merged_Normal_' + s + '.bam'))):
                print "    *** Merged bam files found, skipping..."
                targets_data_processed.append(ngs_classes.WxsSample(_id = s,
                                                                    _bamTumor = os.path.abspath(inputs.dir + '/bamMerged/Merged_Tumor_' + s + '.bam'),
                                                                    _bamNormal = os.path.abspath(inputs.dir + '/bamMerged/Merged_Normal_' + s + '.bam')))
            else:
                
                if inputs.verbose:
                    print "Running java -jar %s MergeSamFiles I=%s O=%s" % (inputs.picard,
                                                                            ' I='.join(sample_ids[s]["TUMOR"]),
                                                                            os.path.abspath(inputs.dir + '/bamMerged/Merged_Tumor_' + s + '.bam'))
                if not inputs.dry:
                    os.system("java -jar %s MergeSamFiles I=%s O=%s" % (inputs.picard,
                                                                        ' I='.join(sample_ids[s]["TUMOR"]),
                                                                        os.path.abspath(inputs.dir + '/bamMerged/Merged_Tumor_' + s + '.bam')))
                    
                if inputs.verbose:
                    print "Running java -jar %s MergeSamFiles I=%s O=%s" % (inputs.picard,
                                                                            ' I='.join(sample_ids[s]["NORMAL"]),
                                                                            os.path.abspath(inputs.dir + '/bamMerged/Merged_Normal_' + s + '.bam'))
                if not inputs.dry:
                    os.system("java -jar %s MergeSamFiles I=%s O=%s" % (inputs.picard,
                                                                        ' I='.join(sample_ids[s]["NORMAL"]),
                                                                        os.path.abspath(inputs.dir + '/bamMerged/Merged_Normal_' + s + '.bam')))        
                    
                targets_data_processed.append(ngs_classes.WxsSample(_id = s,
                                                                    _bamTumor = os.path.abspath(inputs.dir + '/bamMerged/Merged_Tumor_' + s + '.bam'),
                                                                    _bamNormal = os.path.abspath(inputs.dir + '/bamMerged/Merged_Normal_' + s + '.bam')))
            
    elif str.upper(inputs.workflow_type) == "GERMLINESNVINDEL":
        sample_ids = {}
        for r in targets_data:
            if r._id in sample_ids.keys():
                sample_ids[r._id].append(r._processed)
            else:
                sample_ids[r._id] = []
                sample_ids[r._id].append(r._processed)
    
        for s in sample_ids:
            if(os.path.exists(os.path.abspath(inputs.dir + '/bamMerged/Merged_' + s + '.bam'))):
                print "    *** Merged bam files found, skipping..."
                targets_data_processed.append(ngs_classes.WxsSample(_id = s,
                                                                    _bam = os.path.abspath(inputs.dir + '/bamMerged/Merged_' + s + '.bam')))
            
            else:
                if inputs.verbose:
                    print "Running java -jar %s MergeSamFiles I=%s O=%s" % (inputs.picard,
                                                                          ' I='.join(sample_ids[s]),
                                                                          os.path.abspath(inputs.dir + '/bamMerged/Merged_' + s + '.bam'))
                if not inputs.dry:
                    os.system("java -jar %s MergeSamFiles I=%s O=%s" % (inputs.picard,
                                                                        ' I='.join(sample_ids[s]),
                                                                        os.path.abspath(inputs.dir + '/bamMerged/Merged_' + s + '.bam')))
            
                targets_data_processed.append(ngs_classes.WxsSample(_id = s,
                                                                    _bam = os.path.abspath(inputs.dir + '/bamMerged/Merged_' + s + '.bam')))
    return True


def bqsr(inputs, targets_data_processed):
    
    ## Base Recalibration ##
    for s in targets_data_processed:
        if str.upper(inputs.workflow_type) == "SOMATICSNVINDEL":
            if inputs.verbose:
                print "Running %s index %s" % (inputs.samtools, s._bamTumor)
            
            if not inputs.dry:
                os.system("%s index %s" % (inputs.samtools, s._bamTumor))
            
            if inputs.verbose:
                print "Running %s index %s" % (inputs.samtools, s._bamNormal)
            
            if not inputs.dry:
                os.system("%s index %s" % (inputs.samtools, s._bamNormal))
                
                
        elif str.upper(inputs.workflow_type) == "GERMLINESNVINDEL":
            if inputs.verbose:
                print "Running %s index %s" % (inputs.samtools, s._bam)
            
            if not inputs.dry:
                os.system("%s index %s" % (inputs.samtools, s._bam))    
    
    
    ## BQSR (Step 1) ##
    for s in targets_data_processed:
        if str.upper(inputs.workflow_type) == "SOMATICSNVINDEL":
            if(os.path.exists(os.path.abspath(inputs.dir + 'recalib/RecalibrationTableTumor_' + s._id + '.table')) and
               os.path.exists(os.path.abspath(inputs.dir + 'recalib/RecalibrationTableNormal_' + s._id + '.table'))):
                print "    *** Recalibration tables found, skipping..."
                s._recalibTableTumor = os.path.abspath(inputs.dir + 'recalib/RecalibrationTableTumor_' + s._id + '.table')
                s._recalibTableNormal = os.path.abspath(inputs.dir + 'recalib/RecalibrationTableNormal_' + s._id + '.table')
            
                
            else:
                if inputs.verbose:
                    if not(inputs.bed is None):
                        print "Running java -jar %s -T BaseRecalibrator -nct %d -R %s -I %s -knownSites %s -knownSites %s -o %s -L %s" % (inputs.gatk, inputs.threads,
                                                                                                                  inputs.reference,
                                                                                                                  s._bamTumor,
                                                                                                                  inputs.dbSNP,
                                                                                                                  inputs.mills,
                                                                                                                  os.path.abspath(inputs.dir +
                                                                                                                                  'recalib/RecalibrationTableTumor_' +
                                                                                                                                  s._id + '.table'),
                                                                                                                  inputs.bed)
                    else:
                        print "Running java -jar %s -T BaseRecalibrator -nct %d -R %s -I %s -knownSites %s -knownSites %s -o %s" % (inputs.gatk, inputs.threads,
                                                                                                                  inputs.reference,
                                                                                                                  s._bamTumor,
                                                                                                                  inputs.dbSNP,
                                                                                                                  inputs.mills,
                                                                                                                  os.path.abspath(inputs.dir +
                                                                                                                                  'recalib/RecalibrationTableTumor_' +
                                                                                                                                  s._id + '.table'))
                if not inputs.dry:
                    if not(inputs.bed is None):
                        os.system("java -jar %s -T BaseRecalibrator -nct %d -R %s -I %s -knownSites %s -knownSites %s -o %s -L %s" % (inputs.gatk, inputs.threads,
                                                                                                      inputs.reference,
                                                                                                      s._bamTumor,
                                                                                                      inputs.dbSNP,
                                                                                                      inputs.mills,
                                                                                                      os.path.abspath(inputs.dir + 'recalib/RecalibrationTableTumor_' + s._id + '.table'),
                                                                                                      inputs.bed))
                    else:
                        os.system("java -jar %s -T BaseRecalibrator -nct %d -R %s -I %s -knownSites %s -knownSites %s -o %s -L %s" % (inputs.gatk, inputs.threads,
                                                                                                      inputs.reference,
                                                                                                      s._bamTumor,
                                                                                                      inputs.dbSNP,
                                                                                                      inputs.mills,
                                                                                                      os.path.abspath(inputs.dir + 'recalib/RecalibrationTableTumor_' + s._id + '.table'),
                                                                                                      inputs.bed))
                
                s._recalibTableTumor = os.path.abspath(inputs.dir + 'recalib/RecalibrationTableTumor_' + s._id + '.table')
            
                if inputs.verbose:
                    if not(inputs.bed is None):
                        print "Running java -jar %s -T BaseRecalibrator -nct %d -R %s -I %s -knownSites %s -knownSites %s -o %s -L %s" % (inputs.gatk, inputs.threads,
                                                                                                          inputs.reference,
                                                                                                          s._bamNormal,
                                                                                                          inputs.dbSNP,
                                                                                                          inputs.mills,
                                                                                                          os.path.abspath(inputs.dir + 'recalib/RecalibrationTableNormal_' + s._id + '.table'),
                                                                                                          inputs.bed)
                    else:
                        print "Running java -jar %s -T BaseRecalibrator -nct %d -R %s -I %s -knownSites %s -knownSites %s -o %s" % (inputs.gatk, inputs.threads,
                                                                                                          inputs.reference,
                                                                                                          s._bamNormal,
                                                                                                          inputs.dbSNP,
                                                                                                          inputs.mills,
                                                                                                          os.path.abspath(inputs.dir + 'recalib/RecalibrationTableNormal_' + s._id + '.table'))
                if not inputs.dry:
                    if not(inputs.bed is None):
                        os.system("java -jar %s -T BaseRecalibrator -nct %d -R %s -I %s -knownSites %s -knownSites %s -o %s -L %s" % (inputs.gatk, inputs.threads,
                                                                                                      inputs.reference,
                                                                                                      s._bamNormal,
                                                                                                      inputs.dbSNP,
                                                                                                      inputs.mills,
                                                                                                      os.path.abspath(inputs.dir + 'recalib/RecalibrationTableNormal_' + s._id + '.table'),
                                                                                                      inputs.bed))
                    else:
                        os.system("java -jar %s -T BaseRecalibrator -nct %d -R %s -I %s -knownSites %s -knownSites %s -o %s" % (inputs.gatk, inputs.threads,
                                                                                                      inputs.reference,
                                                                                                      s._bamNormal,
                                                                                                      inputs.dbSNP,
                                                                                                      inputs.mills,
                                                                                                      os.path.abspath(inputs.dir + 'recalib/RecalibrationTableNormal_' + s._id + '.table')))
            
                s._recalibTableNormal = os.path.abspath(inputs.dir + 'recalib/RecalibrationTableNormal_' + s._id + '.table')
            
        elif str.upper(inputs.workflow_type) == "GERMLINESNVINDEL":
            if(os.path.exists(os.path.abspath(inputs.dir + 'recalib/RecalibrationTable_' + s._id + '.table'))):
                print "    *** Recalibration table found, skipping..."
                s._recalibTable = os.path.abspath(inputs.dir + 'recalib/RecalibrationTable_' + s._id + '.table')
        
                
            else:
                if inputs.verbose:
                    print "Running java -jar %s -T BaseRecalibrator -nct %d -R %s -I %s -knownSites %s -knownSites %s -o %s" % (inputs.gatk, inputs.threads,
                                                                                                          inputs.reference,
                                                                                                          s._bam,
                                                                                                          inputs.dbSNP,
                                                                                                          inputs.mills,
                                                                                                          os.path.abspath(inputs.dir + 'recalib/RecalibrationTable_' + s._id + '.table'))
                if not inputs.dry:
                    os.system("java -jar %s -T BaseRecalibrator -nct %d -R %s -I %s -knownSites %s -knownSites %s -o %s" % (inputs.gatk, inputs.threads,
                                                                                                      inputs.reference,
                                                                                                      s._bam,
                                                                                                      inputs.dbSNP,
                                                                                                      inputs.mills,
                                                                                                      os.path.abspath(inputs.dir + 'recalib/RecalibrationTable_' + s._id + '.table')))
                
                s._recalibTable = os.path.abspath(inputs.dir + 'recalib/RecalibrationTable_' + s._id + '.table')

        
    ## BQSR (Step 2) ##
    for s in targets_data_processed:
        if str.upper(inputs.workflow_type) == "SOMATICSNVINDEL":
            if(os.path.exists(os.path.abspath(inputs.dir + 'bamCalib/CalibratedTumor_' + s._id + '.bam')) and
               os.path.exists(os.path.abspath(inputs.dir + 'bamCalib/CalibratedNormal_' + s._id + '.bam'))):
                print "    *** BQSR adjusted files found, skipping..."
                s._bamCalibTumor = os.path.abspath(inputs.dir + 'bamCalib/CalibratedTumor_' + s._id + '.bam')
                s._bamCalibNormal = os.path.abspath(inputs.dir + 'bamCalib/CalibratedNormal_' + s._id + '.bam')
            else:
                if inputs.verbose:
                    print "Running: java -jar %s -T PrintReads -nct %d -R %s -I %s -BQSR %s -o %s" % (inputs.gatk, inputs.threads,
                                                                                              inputs.reference,
                                                                                              s._bamTumor,
                                                                                              s._recalibTableTumor,
                                                                                              os.path.abspath(inputs.dir + 'bamCalib/CalibratedTumor_' + s._id + '.bam'))
                if not inputs.dry:
                    os.system("java -jar %s -T PrintReads -nct %d -R %s -I %s -BQSR %s -o %s" % (inputs.gatk, inputs.threads,
                                                                                         inputs.reference,
                                                                                         s._bamTumor,
                                                                                         s._recalibTableTumor,
                                                                                         os.path.abspath(inputs.dir + 'bamCalib/CalibratedTumor_' + s._id + '.bam')))
                
                s._bamCalibTumor = os.path.abspath(inputs.dir + 'bamCalib/CalibratedTumor_' + s._id + '.bam')
            
                if inputs.verbose:
                    print "Running: java -jar %s -T PrintReads -nct %d -R %s -I %s -BQSR %s -o %s" % (inputs.gatk, inputs.threads,
                                                                                              inputs.reference,
                                                                                              s._bamNormal,
                                                                                              s._recalibTableNormal,
                                                                                              os.path.abspath(inputs.dir + 'bamCalib/CalibratedNormal_' + s._id + '.bam'))
                if not inputs.dry:
                    os.system("java -jar %s -T PrintReads -nct %d -R %s -I %s -BQSR %s -o %s" % (inputs.gatk, inputs.threads,
                                                                                         inputs.reference,
                                                                                         s._bamNormal,
                                                                                         s._recalibTableNormal,
                                                                                         os.path.abspath(inputs.dir + 'bamCalib/CalibratedNormal_' + s._id + '.bam')))
                s._bamCalibNormal = os.path.abspath(inputs.dir + 'bamCalib/CalibratedNormal_' + s._id + '.bam')
            
        elif str.upper(inputs.workflow_type) == "GERMLINESNVINDEL":
            if(os.path.exists(os.path.abspath(inputs.dir + 'bamCalib/Calibrated_' + s._id + '.bam'))):
                print "    *** BQSR adjusted files found, skipping..."
                s._bamCalib = os.path.abspath(inputs.dir + 'bamCalib/Calibrated_' + s._id + '.bam')
            else:
                if inputs.verbose:
                    print "Running: java -jar %s -T PrintReads -nct %d -R %s -I %s -BQSR %s -o %s" % (inputs.gatk, inputs.threads,
                                                                                              inputs.reference,
                                                                                              s._bam,
                                                                                              s._recalibTable,
                                                                                              os.path.abspath(inputs.dir + 'bamCalib/Calibrated_' + s._id + '.bam'))
                if not inputs.dry:
                    os.system("java -jar %s -T PrintReads -nct %d -R %s -I %s -BQSR %s -o %s" % (inputs.gatk, inputs.threads,
                                                                                         inputs.reference,
                                                                                         s._bam,
                                                                                         s._recalibTable,
                                                                                         os.path.abspath(inputs.dir + 'bamCalib/Calibrated_' + s._id + '.bam')))
                s._bamCalib = os.path.abspath(inputs.dir + 'bamCalib/Calibrated_' + s._id + '.bam')
    return True

def ngs_mutect(inputs, targets_data_processed):
    for s in targets_data_processed:
        if(os.path.exists(os.path.abspath(os.path.abspath(inputs.dir + '/rawSNP/' + r._id + '.vcf')))):
            print "    *** MuTect raw call file found, skipping..."
        else:
            if inputs.verbose:
                print "Running %s -jar %s -T MuTect --baq OFF -nct 8 --reference_sequence %s --input_file:normal %s --input_file:tumor %s --out %s --coverage_file %s --vcf %s" % (inputs.java17,
                                                                                                    inputs.mutect1,
                                                                                                    inputs.reference,
                                                                                                    s._bamCalibNormal,
                                                                                                    s._bamCalibTumor,
                                                                                                    os.path.abspath(inputs.dir + '/rawSNP/' + r._id + '.out'),
                                                                                                    os.path.abspath(inputs.dir + '/rawSNP/' + r._id + '.wig.txt'),
                                                                                                    os.path.abspath(inputs.dir + '/rawSNP/' + r._id + '.vcf'))
            if not inputs.dry:
                os.system("%s -jar %s -T MuTect --baq OFF -nct 8 --reference_sequence %s --input_file:normal %s --input_file:tumor %s --out %s --coverage_file %s --vcf %s" % (inputs.java17,
                                                                                                        inputs.mutect1,
                                                                                                        inputs.reference,
                                                                                                        s._bamCalibNormal,
                                                                                                        s._bamCalibTumor,
                                                                                                        os.path.abspath(inputs.dir + '/rawSNP/' + r._id + '.out'),
                                                                                                        os.path.abspath(inputs.dir + '/rawSNP/' + r._id + '.wig.txt'),
                                                                                                        os.path.abspath(inputs.dir + '/rawSNP/' + r._id + '.vcf')))
                
    ## Filter ##
    return True





def ngs_haplotypecaller(inputs, targets_data_processed):
    
    ## HaplotypeCaller ##
    files = []
    for s in targets_data_processed:
        files.append(s._bamCalib)
        
    if(os.path.exists(os.path.abspath(inputs.dir + '/rawSNP/GATK_HCaller_' + inputs.title + '.raw.snps.indels.vcf'))):
        print "    *** HaplotypeCaller raw file found, skipping..."
    else:
        if inputs.verbose:
            print "Running: java -jar %s -T HaplotypeCaller -R %s -I %s --num_cpu_threads_per_data_thread %d -stand_call_conf 20 -o %s "\
            "-A BaseQualityRankSumTest -A Coverage -A DepthPerAlleleBySample -A FisherStrand -A QualByDepth -A ReadPosRankSumTest" % (inputs.gatk,
                                                                                                                                      inputs.reference,
                                                                                                                                      " -I ".join(files),
                                                                                                                                      inputs.threads,
                                                                                                                                      os.path.abspath(inputs.dir +
                                                                                                                                                      '/rawSNP/GATK_HCaller_' +
                                                                                                                                                      inputs.title + '.raw.snps.indels.vcf'))
        if not inputs.dry:
            os.system("java -jar %s -T HaplotypeCaller -R %s -I %s --num_cpu_threads_per_data_thread %d -stand_call_conf 20 -o %s "\
                      "-A BaseQualityRankSumTest -A Coverage -A DepthPerAlleleBySample -A FisherStrand -A QualByDepth -A ReadPosRankSumTest" % (inputs.gatk,
                                                                                                                                                inputs.reference,
                                                                                                                                                " -I ".join(files),
                                                                                                                                                inputs.threads,
                                                                                                                                                os.path.abspath(inputs.dir +
                                                                                                                                                                '/rawSNP/GATK_HCaller_' +
                                                                                                                                                                inputs.title + '.raw.snps.indels.vcf')))
            
    
    
    
    if not(inputs.ped is None):
        if(os.path.exists(os.path.abspath(inputs.dir + '/rawSNP/GATK_HCaller_GTPosteriors_' + inputs.title + '.raw.snps.indels.vcf'))):
            print "    *** HaplotypeCaller GenotypePosteriors file found, skipping..."
        else:
            if inputs.verbose:
                print "Running: java -jar %s -T CalculateGenotypePosteriors -R %s -V %s -o %s --skipPopulationPriors -ped %s" % (inputs.gatk,
                                                                                                                                 inputs.reference,
                                                                                                                                 os.path.abspath(inputs.dir +
                                                                                                                                                 '/rawSNP/GATK_HCaller_' +
                                                                                                                                                 inputs.title + '.raw.snps.indels.vcf'),
                                                                                                                                 os.path.abspath(inputs.dir +
                                                                                                                                                 '/rawSNP/GATK_HCaller_GTPosteriors_' +
                                                                                                                                                 inputs.title + '.raw.snps.indels.vcf'),
                                                                                                                                 inputs.ped)
            if not inputs.dry:
                os.system("java -jar %s -T CalculateGenotypePosteriors -R %s -V %s -o %s --skipPopulationPriors -ped %s" % (inputs.gatk,
                                                                                                                            inputs.reference,
                                                                                                                            os.path.abspath(inputs.dir +
                                                                                                                                            '/rawSNP/GATK_HCaller_' +
                                                                                                                                            inputs.title + '.raw.snps.indels.vcf'),
                                                                                                                            os.path.abspath(inputs.dir +
                                                                                                                                            '/rawSNP/GATK_HCaller_GTPosteriors_' +
                                                                                                                                            inputs.title + '.raw.snps.indels.vcf'),
                                                                                                                            inputs.ped))
    
    
    return True


def ngs_haplotypecaller_filter(inputs):

    if not(inputs.ped is None):
        inFile = os.path.abspath(inputs.dir + '/rawSNP/GATK_HCaller_GTPosteriors_' + inputs.title + '.raw.snps.indels.vcf')
    else:
        inFile = os.path.abspath(inputs.dir + '/rawSNP/GATK_HCaller_' + inputs.title + '.raw.snps.indels.vcf')
    
    if(os.path.exists(outFile)):
        print "    *** HaplotypeCaller filtered vcf file found, skipping..."
    else:
        if inputs.verbose:
            print "Running: java -jar %s -T VariantFiltration -R %s -o %s -V %s "\
            "--filterExpression \'QD < 2.0\' --filterName \"LowQD\" "\
            "--filterExpression \'FS > 60.0 && vc.getType() == \"SNP\"\' --filterName \"HighFSSNP\" "\
            "--filterExpression \'FS > 200.0 && vc.getType() == \"INDEL\"\' --filterName \"HighFSINDEL\" "\
            "--filterExpression \'MQ < 40.0 && vc.getType() == \"SNP\"\' --filterName \"LowMQSNP\" "\
            "--filterExpression \'ReadPosRankSum < -8.0 && vc.getType() == \"SNP\"\' --filterName \"LowRPRSSNP\" "\
            "--filterExpression \'ReadPosRankSum < -20.0 && vc.getType() == \"INDEL\"\' --filterName \"LowRPRSINDEL\" "\
            "--genotypeFilterExpression \"DP < 10\" --genotypeFilterName \"LowDPGT\" --setFilteredGtToNocall" % (inputs.gatk,
                                                                                                                 inputs.reference,
                                                                                                                 os.path.abspath(inputs.dir + '/filtSNP/GATK_HCaller_' + inputs.title + '.pre.filtered.snps.indels.vcf'),
                                                                                                                 inFile)
            
        if not inputs.dry:
            os.system("java -jar %s -T VariantFiltration -R %s -o %s -V %s "\
                      "--filterExpression \'QD < 2.0\' --filterName \"LowQD\" "\
                      "--filterExpression \'FS > 60.0 && vc.getType() == \"SNP\"\' --filterName \"HighFSSNP\" "\
                      "--filterExpression \'FS > 200.0 && vc.getType() == \"INDEL\"\' --filterName \"HighFSINDEL\" "\
                      "--filterExpression \'MQ < 40.0 && vc.getType() == \"SNP\"\' --filterName \"LowMQSNP\" "\
                      "--filterExpression \'ReadPosRankSum < -8.0 && vc.getType() == \"SNP\"\' --filterName \"LowRPRSSNP\" "\
                      "--filterExpression \'ReadPosRankSum < -20.0 && vc.getType() == \"INDEL\"\' --filterName \"LowRPRSINDEL\" "\
                      "--genotypeFilterExpression \"DP < 10\" --genotypeFilterName \"LowDPGT\" --setFilteredGtToNocall" % (inputs.gatk,
                                                                                                                           inputs.reference,
                                                                                                                           os.path.abspath(inputs.dir + '/filtSNP/GATK_HCaller_' + inputs.title + '.pre.filtered.snps.indels.vcf'),
                                                                                                                           inFile))
            
        if inputs.verbose:
            print "Running: java -jar %s -T SelectVariants -R %s "\
            "--maxFractionFilteredGenotypes 0.7 --excludeNonVariants --excludeFiltered -selectType SNP -selectType INDEL "\
            "-o %s -V %s" % (inputs.gatk,
                             inputs.reference,
                             os.path.abspath(inputs.dir + '/filtSNP/GATK_HCaller_' + inputs.title + '.filtered.snps.indels.vcf'),
                             os.path.abspath(inputs.dir + '/filtSNP/GATK_HCaller_' + inputs.title + '.pre.filtered.snps.indels.vcf'))
            
            
        if not inputs.dry:
            os.system("Running: java -jar %s -T SelectVariants -R %s "\
            "--maxFractionFilteredGenotypes 0.7 --excludeNonVariants --excludeFiltered -selectType SNP -selectType INDEL "\
            "-o %s -V %s" % (inputs.gatk,
                             inputs.reference,
                             os.path.abspath(inputs.dir + '/filtSNP/GATK_HCaller_' + inputs.title + '.filtered.snps.indels.vcf'),
                             os.path.abspath(inputs.dir + '/filtSNP/GATK_HCaller_' + inputs.title + '.pre.filtered.snps.indels.vcf')))
            
            
    if not(inputs.ped is None):
        if(os.path.exists(os.path.abspath(inputs.dir + '/filtSNP/GATK_HCaller_' + inputs.title + '.ann.filtered.snps.indels.vcf'))):
            print "    *** Annotated vcf file found, skipping..."
        else:
            if inputs.verbose:
                print "Running: java -jar %s -T VariantAnnotator -R %s -o %s -V %s -A PossibleDeNovo -ped %s" & (inputs.gatk,
                                                                                                                 inputs.reference,
                                                                                                                 os.path.abspath(inputs.dir +
                                                                                                                                 '/filtSNP/GATK_HCaller_' +
                                                                                                                                 inputs.title + '.ann.filtered.snps.indels.vcf'),
                                                                                                                 os.path.abspath(inputs.dir +
                                                                                                                                 '/filtSNP/GATK_HCaller_' +
                                                                                                                                 inputs.title + '.filtered.snps.indels.vcf'),
                                                                                                                 inputs.ped)
                
            if not inputs.dry:
                os.system("java -jar %s -T VariantAnnotator -R %s -o %s -V %s -A PossibleDeNovo -ped %s" & (inputs.gatk,
                                                                                                            inputs.reference,
                                                                                                            os.path.abspath(inputs.dir +
                                                                                                                            '/filtSNP/GATK_HCaller_' +
                                                                                                                            inputs.title + '.ann.filtered.snps.indels.vcf'),
                                                                                                            os.path.abspath(inputs.dir +
                                                                                                                            '/filtSNP/GATK_HCaller_' +
                                                                                                                            inputs.title + '.filtered.snps.indels.vcf'),
                                                                                                            inputs.ped))
            
    return True


def ngs_bcftools(inputs, targets_data_processed):
    
    ## Index Bam Files ##
    for s in targets_data_processed:
        if inputs.verbose:
            print 'Running: %s index %s' % (inputs.samtools, s._bamCalib)
        if not inputs.dry:
            os.system('%s index %s' % (inputs.samtools, s._bamCalib))
    
    ## Samtools MPileup ##
    files = []
    for s in targets_data_processed:
        files.append(s._bamCalib)
        
    if(os.path.exists(os.path.abspath(inputs.dir + '/rawSNP/SAMTools_MPileup_' + inputs.title + '.vcf.gz'))):
        print "    *** Mpileup file found, skipping..."
    else:
        if inputs.verbose:
            print 'Running: %s mpileup -d 100000 -vsBf %s -t AD,ADF,ADR,DP,SP -o %s %s' % (inputs.samtools,
                                                                                           inputs.reference,
                                                                                           os.path.abspath(inputs.dir + '/rawSNP/SAMTools_MPileup_' + inputs.title + '.vcf.gz'),
                                                                                           " ".join(files))
        
        if not inputs.dry:
            os.system('%s mpileup -d 100000 -vsBf %s -t AD,ADF,ADR,DP,SP -o %s %s' % (inputs.samtools,
                                                                                      inputs.reference,
                                                                                      os.path.abspath(inputs.dir + '/rawSNP/SAMTools_MPileup_' + inputs.title + '.vcf.gz'),
                                                                                      " ".join(files)))
        
    ## BCFTools Call ##
    if(os.path.exists(os.path.abspath(inputs.dir + '/rawSNP/BCFTools_Call_' + inputs.title + '.vcf.gz'))):
        print "    *** BCFTools raw call file found, skipping..."
    else:
        if inputs.verbose:
            print "Running: %s call -mv -O z -o %s --threads %d %s" % (inputs.bcftools,
                                                                       os.path.abspath(inputs.dir + '/rawSNP/BCFTools_Call_' + inputs.title + '.vcf.gz'),
                                                                       inputs.threads,
                                                                       os.path.abspath(inputs.dir + '/rawSNP/SAMTools_MPileup_' + inputs.title + '.vcf.gz'))
        if not inputs.dry:
            os.system("%s call -mv -O z -o %s --threads %d %s" % (inputs.bcftools,
                                                                  os.path.abspath(inputs.dir + '/rawSNP/BCFTools_Call_' + inputs.title + '.vcf.gz'),
                                                                  inputs.threads,
                                                                  os.path.abspath(inputs.dir + '/rawSNP/SAMTools_MPileup_' + inputs.title + '.vcf.gz')))
            
            
    ## BCFTools Filter ##
    if(os.path.exists(os.path.abspath(inputs.dir + '/filtSNP/BCFTools_Filt_' + inputs.title + '.vcf'))):
        print "    *** BCFTools filtered call file found, skipping..."
    else:
        if inputs.verbose:
            print "Running: %s filter -O v --threads %d -g20 -G10 -i \'FMT/DP>15 & QUAL>20 & GT!~\"\.\" & MQ>40 & RPB>0.1\' -o %s %s" % (inputs.bcftools,
                                                                                                                       inputs.threads,
                                                                                                                       os.path.abspath(inputs.dir + '/filtSNP/BCFTools_Filt_' + inputs.title + '.vcf'),
                                                                                                                       os.path.abspath(inputs.dir + '/rawSNP/BCFTools_Call_' + inputs.title + '.vcf.gz'))
        if not inputs.dry:
            os.system("%s filter -O v --threads %d -g20 -G10 -i \'FMT/DP>15 & QUAL>20 & GT!~\"\.\" & MQ>40 & RPB>0.1\' -o %s %s" % (inputs.bcftools,
                                                                                                                  inputs.threads,
                                                                                                                  os.path.abspath(inputs.dir + '/filtSNP/BCFTools_Filt_' + inputs.title + '.vcf'),
                                                                                                                  os.path.abspath(inputs.dir + '/rawSNP/BCFTools_Call_' + inputs.title + '.vcf.gz')))
    return True        
        
def ngs_varscan2_germline(inputs, targets_data_processed):
    
    ## Index Bam Files ##
    for s in targets_data_processed:
        if inputs.verbose:
            print 'Running: %s index %s' % (inputs.samtools, s._bamCalib)
        if not inputs.dry:
            os.system('%s index %s' % (inputs.samtools, s._bamCalib))
    
    ## Samtools MPileup ##
    files = []
    for s in targets_data_processed:
        files.append(s._bamCalib)
        
    if(os.path.exists(os.path.abspath(inputs.dir + '/rawSNP/SAMTools_MPileup_' + inputs.title + '.mpileup'))):
        print "    *** Mpileup file found, skipping..."
    else:
        if inputs.verbose:
            print 'Running: %s mpileup -Bf %s -o %s %s' % (inputs.samtools,
                                                           inputs.reference,
                                                           os.path.abspath(inputs.dir + '/rawSNP/SAMTools_MPileup_' + inputs.title + '.mpileup'),
                                                           " ".join(files))
        
        if not inputs.dry:
            os.system('%s mpileup -Bf %s -o %s %s' % (inputs.samtools,
                                                      inputs.reference,
                                                      os.path.abspath(inputs.dir + '/rawSNP/SAMTools_MPileup_' + inputs.title + '.mpileup'),
                                                      " ".join(files)))
    ## VarScan2 ##
    if(os.path.exists(os.path.abspath(inputs.dir + '/rawSNP/VarScan2_Call_' + inputs.title + '.vcf.gz'))):
        print "    *** VarScan2 call file found, skipping..."
    else:
        if inputs.verbose:
            print "Running: java -jar %s mpileup2cns %s --variants --output-vcf --min-coverage 20 --strand-filter 1 > %s" % (inputs.varscan2,
                                                                                                                             os.path.abspath(inputs.dir + '/rawSNP/SAMTools_MPileup_' + inputs.title + '.mpileup'),
                                                                                                                             os.path.abspath(inputs.dir + '/filtSNP/VarScan2_Call_' + inputs.title + '.vcf'))
        if not inputs.dry:
            os.system("java -jar %s mpileup2cns %s --variants --output-vcf --min-coverage 20 --strand-filter 1 > %s" % (inputs.varscan2,
                                                                                                                        os.path.abspath(inputs.dir + '/rawSNP/SAMTools_MPileup_' + inputs.title + '.mpileup'),
                                                                                                                        os.path.abspath(inputs.dir + '/filtSNP/VarScan2_Call_' + inputs.title + '.vcf')))
            
            
    return True



def ngs_varscan2_somatic(inputs, targets_data_processed):
    
    for s in targets_data_processed:
        if(os.path.exists(os.path.abspath(inputs.dir + '/rawSNP/Normal_Tumor_' + s._id + '.mpileup'))):
            print "   *** Pileup file found, skipping..."
            s._pileup = os.path.abspath(inputs.dir + '/rawSNP/Normal_Tumor_' + s._id + '.mpileup')
        else:    
            if inputs.verbose:
                print "Running: %s mpileup -Bf %s -o %s %s %s" % (inputs.samtools,
                                                              inputs.reference,
                                                              os.path.abspath(inputs.dir + '/rawSNP/Normal_Tumor_' + s._id + '.mpileup'),
                                                              s._bamCalibNormal,
                                                              s._bamCalibTumor)
            if not inputs.dry:
                os.system("%s mpileup -d -Bf %s -o %s %s %s" % (inputs.samtools,
                                                              inputs.reference,
                                                              os.path.abspath(inputs.dir + '/rawSNP/Normal_Tumor_' + s._id + '.mpileup'),
                                                              s._bamCalibNormal,
                                                              s._bamCalibTumor)) 
            s._pileup = os.path.abspath(inputs.dir + '/rawSNP/Normal_Tumor_' + s._id + '.mpileup')
                
    for s in targets_data_processed:
        if(os.path.exists(os.path.abspath(inputs.dir + '/rawSNP/' + s._id + '.snp.vcf')) and
           os.path.exists(os.path.abspath(inputs.dir + '/rawSNP/' + s._id + '.indel.vcf'))):
            print "   *** VarScan2 result files found, skipping..."
            s._snpOut = os.path.abspath(inputs.dir + '/rawSNP/' + s._id + '.snp.vcf')
            s._indelOut = os.path.abspath(inputs.dir + '/rawSNP/' + s._id + '.indel.vcf')
        
        else:
            if inputs.verbose:
                print "Running: java -jar %s somatic %s %s --mpileup 1 --min-coverage 10 --min-var-freq 0.01 --p-value 0.1 "\
                "--min-freq-for-hom 0.8 --somatic-p-value 0.1 --strand-filter 1 --normal-purity 0.95 --tumor-purity 0.85 --output-vcf" % (inputs.varscan2,
                                                                                                                              s._pileup,
                                                                                                                              os.path.abspath(inputs.dir + '/rawSNP/' + s._id))
            if not inputs.dry:
                os.system("java -jar %s somatic %s %s --mpileup 1 --min-coverage 10 --min-var-freq 0.01 --p-value 0.1 "\
                          "--min-freq-for-hom 0.8 --somatic-p-value 0.1 --strand-filter 1 --normal-purity 0.95 --tumor-purity 0.85 --output-vcf" %  (inputs.varscan2,
                                                                                                                                        s._pileup,
                                                                                                                                        os.path.abspath(inputs.dir + '/rawSNP/' + s._id)))
            s._snpOut = os.path.abspath(inputs.dir + '/rawSNP/' + s._id + '.snp.vcf')
            s._indelOut = os.path.abspath(inputs.dir + '/rawSNP/' + s._id + '.indel.vcf')
            
    
    ## Filter ##
    for s in targets_data_processed:
        if(os.path.exists(os.path.abspath(inputs.dir + '/filtSNP/' + s._id + '.filtered.snp.vcf'))):
            print "    *** VarScan2 filtered result files found, skipping..."
        else:
            if inputs.verbose:
                print "Running: java -jar %s somaticFilter %s --indel-file %s --output-file %s --p-value 0.1 --min-var-freq 0.01" % (inputs.varscan2,
                                                                                                                                     s._snpOut,
                                                                                                                                     s._indelOut,
                                                                                                                                     os.path.abspath(inputs.dir +
                                                                                                                                                     '/filtSNP/' +
                                                                                                                                                     s._id +
                                                                                                                                                     '.filtered.snp.vcf'))
            if not inputs.dry:
                os.system("java -jar %s somaticFilter %s --indel-file %s --output-file %s --p-value 0.1 --min-var-freq 0.01" % (inputs.varscan2,
                                                                                                                                s._snpOut,
                                                                                                                                s._indelOut,
                                                                                                                                os.path.abspath(inputs.dir +
                                                                                                                                                '/filtSNP/' +
                                                                                                                                                s._id +
                                                                                                                                                '.filtered.snp.vcf')))
                
                
def ngs_strelka_somatic(inputs, targets_data_processed):
    
    ## Index Bam Files ##
    for s in targets_data_processed:
        if inputs.verbose:
            print 'Running: %s index %s' % (inputs.samtools, s._bamCalibNormal)
        if not inputs.dry:
            os.system('%s index %s' % (inputs.samtools, s._bamCalibNormal))
        if inputs.verbose:
            print 'Running: %s index %s' % (inputs.samtools, s._bamCalibTumor)
        if not inputs.dry:
            os.system('%s index %s' % (inputs.samtools, s._bamCalibTumor))
            
    
    os.mkdir(os.path.abspath(inputs.dir + '/STRELKA_RUN/'))
    
    ## Run Manta ##
    for s in targets_data_processed:
        
        s._mantaWD = os.path.abspath(inputs.dir + '/STRELKA_RUN/mantaWD_' + s._id)
        os.mkdir(s._mantaWD)
        
        if inputs.exome:
            if inputs.verbose:
                if not(inputs.cbed is None):
                    print "Running: python %s --normalBam %s --tumorBam %s --exome --referenceFasta %s --runDir %s --callRegions %s" % (os.path.abspath(inputs.mantabin + '/configManta.py'),
                                                                                                                                        s._bamCalibNormal,
                                                                                                                                        s._bamCalibTumor,
                                                                                                                                        inputs.reference,
                                                                                                                                        s._mantaWD,
                                                                                                                                        inputs.cbed)
                else:
                    print "Running: python %s --normalBam %s --tumorBam %s --exome --referenceFasta %s --runDir %s" % (os.path.abspath(inputs.mantabin + '/configManta.py'),
                                                                                                                       s._bamCalibNormal,
                                                                                                                       s._bamCalibTumor,
                                                                                                                       inputs.reference,
                                                                                                                       s._mantaWD)
                    
                    
            if not inputs.dry:
                if not(inputs.cbed is None):
                    os.system("python %s --normalBam %s --tumorBam %s --exome --referenceFasta %s --runDir %s --callRegions %s" % (os.path.abspath(inputs.mantabin + '/configManta.py'),
                                                                                                                                        s._bamCalibNormal,
                                                                                                                                        s._bamCalibTumor,
                                                                                                                                        inputs.reference,
                                                                                                                                        s._mantaWD,
                                                                                                                                        inputs.cbed))
                else:
                    os.system("python %s --normalBam %s --tumorBam %s --exome --referenceFasta %s --runDir %s" % (os.path.abspath(inputs.mantabin + '/configManta.py'),
                                                                                                                       s._bamCalibNormal,
                                                                                                                       s._bamCalibTumor,
                                                                                                                       inputs.reference,
                                                                                                                       s._mantaWD))
        else:
            if inputs.verbose:
                if not(inputs.cbed is None):
                    print "Running: python %s --normalBam %s --tumorBam %s --referenceFasta %s --runDir %s --callRegions %s" % (os.path.abspath(inputs.mantabin + '/configManta.py'),
                                                                                                                                        s._bamCalibNormal,
                                                                                                                                        s._bamCalibTumor,
                                                                                                                                        inputs.reference,
                                                                                                                                        s._mantaWD,
                                                                                                                                        inputs.cbed)
                else:
                    print "Running: python %s --normalBam %s --tumorBam %s --referenceFasta %s --runDir %s" % (os.path.abspath(inputs.mantabin + '/configManta.py'),
                                                                                                                       s._bamCalibNormal,
                                                                                                                       s._bamCalibTumor,
                                                                                                                       inputs.reference,
                                                                                                                       s._mantaWD)
                    
                    
            if not inputs.dry:
                if not(inputs.cbed is None):
                    os.system("python %s --normalBam %s --tumorBam %s --referenceFasta %s --runDir %s --callRegions %s" % (os.path.abspath(inputs.mantabin + '/configManta.py'),
                                                                                                                                        s._bamCalibNormal,
                                                                                                                                        s._bamCalibTumor,
                                                                                                                                        inputs.reference,
                                                                                                                                        s._mantaWD,
                                                                                                                                        inputs.cbed))
                else:
                    os.system("python %s --normalBam %s --tumorBam %s --referenceFasta %s --runDir %s" % (os.path.abspath(inputs.mantabin + '/configManta.py'),
                                                                                                                       s._bamCalibNormal,
                                                                                                                       s._bamCalibTumor,
                                                                                                                       inputs.reference,
                                                                                                                       s._mantaWD))
                    
        if inputs.verbose:
            print "Running: python %s -m local -j %d" % (os.path.abspath(s._mantaWD + '/runWorkflow.py'),
                                                         inputs.threads)
            
        if not inputs.dry:
            os.system("python %s -m local -j %d" % (os.path.abspath(s._mantaWD + '/runWorkflow.py'),
                                                         inputs.threads))
            
            
        ## Run Strelka ##
        for s in targets_data_processed:
            
            s._strelkaWD = os.path.abspath(inputs.dir + '/STRELKA_RUN/strelkaWD_' + s._id)
            os.mkdir(s._strelkaWD)
        
            if inputs.exome:
                if inputs.verbose:
                    if not(inputs.cbed is None):
                        print "Running: python %s --normalBam %s --tumorBam %s --referenceFasta %s --indelCandidates %s --exome --callRegions %s --runDir %s" % (os.path.abspath(inputs.strelkabin + '/configureStrelkaSomaticWorkflow.py'),
                                                                                                                                                              s._bamCalibNormal,
                                                                                                                                                              s._bamCalibTumor,
                                                                                                                                                              inputs.reference,
                                                                                                                                                              os.path.abspath(s._mantaWD + '/results/variants/candidateSmallIndels.vcf.gz'),
                                                                                                                                                              inputs.cbed,
                                                                                                                                                              s._strelkaWD)
                    else:
                        print "Running: python %s --normalBam %s --tumorBam %s --referenceFasta %s --indelCandidates %s --exome --runDir %s" % (os.path.abspath(inputs.strelkabin + '/configureStrelkaSomaticWorkflow.py'),
                                                                                                                                                              s._bamCalibNormal,
                                                                                                                                                              s._bamCalibTumor,
                                                                                                                                                              inputs.reference,
                                                                                                                                                              os.path.abspath(s._mantaWD + '/results/variants/candidateSmallIndels.vcf.gz'),
                                                                                                                                                              s._strelkaWD)
                        
                if not inputs.dry:
                    if not(inputs.cbed is None):
                        os.system("python %s --normalBam %s --tumorBam %s --referenceFasta %s --indelCandidates %s --exome --callRegions %s --runDir %s" % (os.path.abspath(inputs.strelkabin + '/configureStrelkaSomaticWorkflow.py'),
                                                                                                                                                              s._bamCalibNormal,
                                                                                                                                                              s._bamCalibTumor,
                                                                                                                                                              inputs.reference,
                                                                                                                                                              os.path.abspath(s._mantaWD + '/results/variants/candidateSmallIndels.vcf.gz'),
                                                                                                                                                              inputs.cbed,
                                                                                                                                                              s._strelkaWD))
                    else:
                        os.system("python %s --normalBam %s --tumorBam %s --referenceFasta %s --indelCandidates %s --exome --runDir %s" % (os.path.abspath(inputs.strelkabin + '/configureStrelkaSomaticWorkflow.py'),
                                                                                                                                                              s._bamCalibNormal,
                                                                                                                                                              s._bamCalibTumor,
                                                                                                                                                              inputs.reference,
                                                                                                                                                              os.path.abspath(s._mantaWD + '/results/variants/candidateSmallIndels.vcf.gz'),
                                                                                                                                                              s._strelkaWD))
            else:
                if inputs.verbose:
                    if not(inputs.cbed is None):
                        print "Running: python %s --normalBam %s --tumorBam %s --referenceFasta %s --indelCandidates %s --callRegions %s --runDir %s" % (os.path.abspath(inputs.strelkabin + '/configureStrelkaSomaticWorkflow.py'),
                                                                                                                                                              s._bamCalibNormal,
                                                                                                                                                              s._bamCalibTumor,
                                                                                                                                                              inputs.reference,
                                                                                                                                                              os.path.abspath(s._mantaWD + '/results/variants/candidateSmallIndels.vcf.gz'),
                                                                                                                                                              inputs.cbed,
                                                                                                                                                              s._strelkaWD)
                    else:
                        print "Running: python %s --normalBam %s --tumorBam %s --referenceFasta %s --indelCandidates %s --runDir %s" % (os.path.abspath(inputs.strelkabin + '/configureStrelkaSomaticWorkflow.py'),
                                                                                                                                                              s._bamCalibNormal,
                                                                                                                                                              s._bamCalibTumor,
                                                                                                                                                              inputs.reference,
                                                                                                                                                              os.path.abspath(s._mantaWD + '/results/variants/candidateSmallIndels.vcf.gz'),
                                                                                                                                                              s._strelkaWD)
                        
                if not inputs.dry:
                    if not(inputs.cbed is None):
                        os.system("python %s --normalBam %s --tumorBam %s --referenceFasta %s --indelCandidates %s --callRegions %s --runDir %s" % (os.path.abspath(inputs.strelkabin + '/configureStrelkaSomaticWorkflow.py'),
                                                                                                                                                              s._bamCalibNormal,
                                                                                                                                                              s._bamCalibTumor,
                                                                                                                                                              inputs.reference,
                                                                                                                                                              os.path.abspath(s._mantaWD + '/results/variants/candidateSmallIndels.vcf.gz'),
                                                                                                                                                              inputs.cbed,
                                                                                                                                                              s._strelkaWD))
                    else:
                        os.system("python %s --normalBam %s --tumorBam %s --referenceFasta %s --indelCandidates %s --runDir %s" % (os.path.abspath(inputs.strelkabin + '/configureStrelkaSomaticWorkflow.py'),
                                                                                                                                                              s._bamCalibNormal,
                                                                                                                                                              s._bamCalibTumor,
                                                                                                                                                              inputs.reference,
                                                                                                                                                              os.path.abspath(s._mantaWD + '/results/variants/candidateSmallIndels.vcf.gz'),
                                                                                                                                                              s._strelkaWD))
    return True





def ngs_strelka_germline(inputs, targets_data_processed):
    
    if(os.path.exists(os.path.abspath(inputs.dir + 'STRELKA_RUN_GERMLINE'))):
        print "    *** Strelka run directory found, skipping..."
        
    else:
        ## Index Bam Files ##
        for s in targets_data_processed:
            if inputs.verbose:
                print 'Running: %s index %s' % (inputs.samtools, s._bamCalibNormal)
            if not inputs.dry:
                os.system('%s index %s' % (inputs.samtools, s._bamCalibNormal))
            if inputs.verbose:
                print 'Running: %s index %s' % (inputs.samtools, s._bamCalibTumor)
            if not inputs.dry:
                os.system('%s index %s' % (inputs.samtools, s._bamCalibTumor))
            
            
            runDir = os.path.abspath(inputs.dir + '/STRELKA_RUN_GERMLINE/')
            os.mkdir(runDir)
    
            files = []
            for s in targets_data_processed:
                files.append(s._bamCalib)
        
            if inputs.exome:
                if inputs.verbose:
                    if not(inputs.cbed is None):
                        print "Running python %s --referenceFasta %s --exome --runDir %s --callRegions %s --bam %s" % (inputs.strelkabin + '/configureStrelkaGermlineWorkflow.py',
                                                                                                                       inputs.reference,
                                                                                                                       runDir,
                                                                                                                       inputs.cbed,
                                                                                                                       " --bam ".join(files))
                    else:
                        print "Running python %s --referenceFasta %s --exome --runDir %s --bam %s" % (inputs.strelkabin + '/configureStrelkaGermlineWorkflow.py',
                                                                                                      inputs.reference,
                                                                                                      runDir,
                                                                                                      " --bam ".join(files))
                if not(inputs.dry):
                    if not(inputs.cbed is None):
                        os.system("python %s --referenceFasta %s --exome --runDir %s --callRegions %s --bam %s" % (inputs.strelkabin + '/configureStrelkaGermlineWorkflow.py',
                                                                                                                       inputs.reference,
                                                                                                                       runDir,
                                                                                                                       inputs.cbed,
                                                                                                                       " --bam ".join(files)))
                        
                    else:
                        os.system("python %s --referenceFasta %s --exome --runDir %s --bam %s" % (inputs.strelkabin + '/configureStrelkaGermlineWorkflow.py',
                                                                                                  inputs.reference,
                                                                                                  runDir,
                                                                                                  " --bam ".join(files)))
            else:
                if inputs.verbose:
                    if not(inputs.cbed is None):
                        print "Running python %s --referenceFasta %s --runDir %s --callRegions %s --bam %s" % (inputs.strelkabin + '/configureStrelkaGermlineWorkflow.py',
                                                                                                                       inputs.reference,
                                                                                                                       runDir,
                                                                                                                       inputs.cbed,
                                                                                                                       " --bam ".join(files))
                    else:
                        print "Running python %s --referenceFasta %s --runDir %s --bam %s" % (inputs.strelkabin + '/configureStrelkaGermlineWorkflow.py',
                                                                                                      inputs.reference,
                                                                                                      runDir,
                                                                                                      " --bam ".join(files))
                if not(inputs.dry):
                    if not(inputs.cbed is None):
                        os.system("python %s --referenceFasta %s --runDir %s --callRegions %s --bam %s" % (inputs.strelkabin + '/configureStrelkaGermlineWorkflow.py',
                                                                                                                       inputs.reference,
                                                                                                                       runDir,
                                                                                                                       inputs.cbed,
                                                                                                                       " --bam ".join(files)))
                        
                    else:
                        os.system("python %s --referenceFasta %s --runDir %s --bam %s" % (inputs.strelkabin + '/configureStrelkaGermlineWorkflow.py',
                                                                                                  inputs.reference,
                                                                                                  runDir,
                                                                                                  " --bam ".join(files)))
                    
    return True
