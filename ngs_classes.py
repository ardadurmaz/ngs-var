######################################################
import sys, os, argparse
import os.path, re
import glob
######################################################


class ngsExcept(Exception):
	
	def __init__(self, msg):
		self.msg = msg

class ngsInputs:
	
	def __init__(self):
		## Resources ##
		self.bowtie2 = ""
		self.bowtie2build = ""
		self.bwa = ""
		self.cutadapt = ""
		self.afterqc = ""
		self.skewer = ""
		self.sambamba = ""
		self.gatk = ""
		self.reference = ""
		self.knownsites = ""
		self.strelkabin = ""
		self.mantabin = ""
		self.cnvkit = ""
		
		
		## Arguments ##
		self.dir = ""
		self.in_file = ""
		self.verbose = ""
		self.dry = ""
		self.clear = False
		self.title = ""
		self.bed = ""
		self.cbed = ""
		self.bqsr = False
		self.upload = False
		self.aligner = ""
		self.trimmer = ""
		self.workflow = ""
		self.tool = ""
		self.exome = False


class WxsSample:
    def __init__(self, _id, _bam = None, _bamTumor = None, _bamNormal = None):
        self._id = _id
        self._bamTumor = _bamTumor
        self._bamNormal = _bamNormal
        self._bam = _bam
        self._recalibTable = None
        self._recalibTableTumor = None
        self._recalibTableNormal = None
        self._bamCalibTumor = None
        self._bamCalibNormal = None
        self._bamCalib = None
        self._pileup = None
        self._snpOut = None
        self._indelOut = None
        self._hc_raw = None
        self._mantaWD = None
        self._strelkaWD = None

class WxsSomatic:
    def __init__(self, _r1 = None, _r2 = None, _lane = None, _id = None, _plat = None, _lib = None, _type = None):
        self._r1 = _r1
        self._r2 = _r2
        self._lane = _lane
        self._id = _id
        self._plat = _plat
        self._lib = _lib
        self._type = _type
        self._trimmedR1 = None
        self._trimmedR2 = None
        self._aligned = None
        self._sorted = None
        self._marked = None
        self._markedM = None

class WxsGermline:
        def __init__(self, _r1 = None, _r2 = None, _lane = None, _id = None, _plat = None, _lib = None):
            self._r1 = _r1
            self._r2 = _r2
            self._lane = _lane
            self._id = _id
            self._plat = _plat
            self._lib = _lib
            self._trimmedR1 = None
            self._trimmedR2 = None
            self._aligned = None
            self._sorted = None
            self._marked = None
            self._markedM = None
