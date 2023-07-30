######################################################
import sys, os, argparse
import os.path, re
import glob
import ngs_classes
from ngs_functions import print_log, run_command
######################################################

def ngs_rna_index(inputs):
    print_log(inputs, "\n<---> RNA Indexing <--->\n")
    runDir = os.path.abspath(inputs.dir + '/RNA_SEQ/')
    error_msg = "[ERROR] Failed to perform RNA indexing."
    if(os.path.exists(runDir + '/Genome')):
        print_log(inputs, "[INFO] Indices found, skipping.")
        return 0
    else:
        run_command(inputs, "%s --runThreadN %d --runMode genomeGenerate --genomeDir %s --genomeFastaFiles %s --sjdbGTFfile %s" % (inputs.star, inputs.threads, runDir, inputs.reference, inputs.annotation), error_msg)
    return 0

def ngs_rna_align(inputs, targets_data):
    print_log(inputs, "\n<---> RNA Alignment <--->\n") 
    runDir = os.path.abspath(inputs.dir + '/RNA_SEQ/')
    tempDir = "~/star_temp/"
    r1 = ",".join([s._trimmedR1 for s in targets_data])
    r2 = ",".join([s._trimmedR2 for s in targets_data])
    readFilesIn = " ".join([r1, r2])
    error_msg = "[ERROR] Failed to perform RNA alignment."
    run_command(inputs, "%s --runThreadN %d --genomeDir %s --readFilesIn %s --outTmpDir %s" % (inputs.star, inputs.threads, runDir, readFilesIn, tempDir), error_msg)

    return 0
