######################################################
import sys, os, argparse
import os.path, re
import glob
import ngs_classes
from ngs_functions import print_log, run_command
######################################################

def ngs_quantification(inputs):

    print_log(inputs, "\n<---> RNA Quantification <--->\n")
    runDir = os.path.abspath(inputs.dir + '/RNA_SEQ/')
    if os.path.exists(runDir):
        print_log(inputs, "[INFO] RNASeq directory found, skipping.")
    else:
        try:
            os.mkdir(runDir)
        except OSError:
            raise ngs_classes.ngsExcept(f"[ERROR] Failed to create directory {runDir}")

    aligned = [f for f in os.listdir(runDir) if f.endswith('Aligned.toTranscriptome.out.bam')]
    
    for sample in aligned:
        sample_id = sample.split('Aligned.toTranscriptome.out.bam')[0] 
        outFile = os.path.join(runDir, f"Quant_{sample_id}")
        
        error_msg = f"[ERROR] Failed to perform RNA quantification for sample {sample_id}."
        run_command(inputs, "%s quant --posBias --gcBias --seqBias --threads %d --libType ISR -t %s -a %s -o %s --numBootstraps 100" % (inputs.salmon, inputs.threads, inputs.reference, runDir + "/" + sample, outFile), error_msg)
    
    return 0
