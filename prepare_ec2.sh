#!/bin/bash

set -e

PATHTemp=/opt/storage
# Format and mount SSD(s).
if [ ! -d $PATHTemp ]; then
    mkdir $PATHTemp
    mkfs -t ext3 /dev/nvme0n1
    mount -t ext3 /dev/nvme0n1 $PATHTemp
fi


apt-get update -y
apt-get upgrade -y
apt-get install -y build-essential git

## Configure awscli
if [ ! -d /root/.aws ]; then
    mkdir /root/.aws
    echo -e '[default]\naws_access_key_id = ****************\naws_secret_access_key = ****************************\n' > ~/.aws/config
fi

if [ ! -d /mnt/tools ]; then
    mkdir /mnt/tools
fi

echo " *** Installing Anaconda ***\n"
wget https://repo.continuum.io/archive/Anaconda2-5.0.1-Linux-x86_64.sh && bash Anaconda* -b
echo 'export PATH=$HOME/anaconda2/bin:$PATH' >> ~/.bashrc && source ~/.bashrc
conda update -y conda

echo " *** Installing Required Tools ***\n"
conda install -y -c conda-forge -c bioconda -c aroth85\
                                   samtools\
                                   bcftools\
                                   htslib\
                                   bowtie2\
                                   bwa\
								   sambamba\
                                   strelka\
                                   manta\
                                   skewer\
                                   cutadapt\
                                   sambamba\
								   cnvkit\
								   qiime2


git clone https://github.com/atks/vt.git && cd vt && make && cd $HOME


#echo " *** Installing ENSEMBL VEP ***\n"
#conda install -y -c conda-forge -c bioconda ensembl-vep && vep_install -a -cf -s homo_sapiens -y GRCh38 -c $PATHTemp --CONVERT


echo "Done."

