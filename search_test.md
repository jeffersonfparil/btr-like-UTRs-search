# Search for the exclusively-meiocyte-ProphaseI-specifically expressed BTR1-like and BTR2-like genes in GoldenPromise Barley

## Setup working directories
```{sh}
DIR=/data-weedomics-3/BTR-like-UTRs/search_test.md/
DIR_FASTQ=${DIR}/FASTQ ### Assumes the reads are free of adapter sequences
DIR_REF=${DIR}/REF
REF=${DIR_REF}/GoldenPromise.fasta
cd $DIR
```

## Install tools
```{sh}
sudo apt install sra-toolkit
# sudo apt install cutadapt
git clone https://github.com/DaehwanKimLab/hisat2.git
cd hisat2
make
PATH=${PATH}:$(pwd)
cd -
sudo ln -s /usr/bin/python3 /usr/bin/python ### hisat2 needs python and python2 is now fully deprecated and we may have to specify that the default python is python3

```

## Download PRJNA558196 (Barakate et al 2021) anther RNAseq data
```{sh}
prefetch \
    -O ${DIR}/FASTQ \
    PRJNA558196
fastq-dump \
    --split-files ${DIR}/FASTQ/PRJNA558196.sra
```


## Map the RNAseq to the Golden Promise reference genome
```{sh}
### (1) build the reference index: OUTPUT: ${REF%.fasta*}
time \
hisat2-build \
    ${REF} \
    ${REF%.fasta*}
### (2) map
echo '#!/bin/bash
REF=$1
R1=$2
R2=$3
SAMOUT=${R1%.fastq*}.sam
hisat2 \
    -x ${REF%.fasta*} \
    -1 ${R1} \
    -2 ${R2} \
    -S ${SAMOUT}
' > map_RNAseq.sh
chmod +x map_RNAseq.sh
time \
parallel --link \
./map_RNAseq.sh \
    ${REF} \
    {1} \
    {2} \
    ::: $(ls ${DIR_FASTQ}/*_1.fastq.gz) \
    ::: $(ls ${DIR_FASTQ}/*_2.fastq.gz)

```
