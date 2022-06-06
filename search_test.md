# Search for the exclusively-meiocyte-ProphaseI-specifically expressed BTR1-like and BTR2-like genes in GoldenPromise Barley

## Setup working directories
```{sh}
DIR=/data-weedomics-3/BTR-like-UTRs/search_test.md
DIR_SRA=${DIR}/SRA ### Assumes the reads are free of adapter sequences
DIR_FASTQ=${DIR}/FASTQ ### Assumes the reads are free of adapter sequences
DIR_REF=${DIR}/REF
REF=${DIR_REF}/GoldenPromise.fasta
PATH=${PATH}:${DIR}/mmseqs/bin
PATH=${PATH}:${DIR}
cd $DIR
```

## Install tools
```{sh}
sudo apt install -y default-jre
wget https://github.com/soedinglab/MMseqs2/releases/download/13-45111/mmseqs-linux-avx2.tar.gz
tar -xvzf mmseqs-linux-avx2.tar.gz
cd mmseqs/bin/
./mmseqs -h
cd ${DIR}
rm mmseqs-linux-avx2.tar.gz

wget http://www.jstacs.de/download.php?which=GeMoMa
mv 'download.php?which=GeMoMa' GeMoMa.zip
unzip GeMoMa.zip
java -jar GeMoMa-1.8.jar CLI GeMoMaPipeline
rm GeMoMa.zip

sudo apt install -y sra-toolkit ncbi-entrez-direct cmake
# sudo apt install cutadapt
git clone https://github.com/DaehwanKimLab/hisat2.git
cd hisat2
make
PATH=${PATH}:$(pwd)
sudo ln -s /usr/bin/python3 /usr/bin/python ### hisat2 needs python and python2 is now fully deprecated and we may have to specify that the default python is python3
cd ${DIR}

wget https://github.com/ncbi/TPMCalculator/releases/download/0.0.2/TPMCalculator-0.0.2.x86-linux.tar.gz
tar -xvzf TPMCalculator-0.0.2.x86-linux.tar.gz
PATH=${PATH}:$(pwd)
rm TPMCalculator-0.0.2.x86-linux.tar.gz
```

## Map the Morex v3 annotations into the GoldenPromise genome assembly
1. Dowload the genome and annotation files of Morex V3 the representative reference genome of Hordeum vulgare
```{sh}
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/904/849/725/GCF_904849725.1_MorexV3_pseudomolecules_assembly/GCF_904849725.1_MorexV3_pseudomolecules_assembly_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/904/849/725/GCF_904849725.1_MorexV3_pseudomolecules_assembly/GCF_904849725.1_MorexV3_pseudomolecules_assembly_genomic.gff.gz
gunzip -c GCF_904849725.1_MorexV3_pseudomolecules_assembly_genomic.fna.gz > ${DIR_REF}/Morex_V3.fasta
gunzip -c GCF_904849725.1_MorexV3_pseudomolecules_assembly_genomic.gff.gz > ${DIR_REF}/Morex_V3.gff
rm GCF_904849725.1_MorexV3_pseudomolecules_assembly_genomic.fna.gz
rm GCF_904849725.1_MorexV3_pseudomolecules_assembly_genomic.gff.gz
```

2. Map the genome annotations of Morex V3 into the GoldenPromise genome assembly
```{sh}
time \
java -jar -Xmx280G GeMoMa-1.8.jar CLI \
        GeMoMaPipeline \
        threads=16 \
        GeMoMa.Score=ReAlign \
        AnnotationFinalizer.r=NO \
        p=true pc=true pgr=true \
        o=true \
        t=${REF} \
        i=Morex_V3 \
        a=${DIR_REF}/Morex_V3.gff \
        g=${DIR_REF}/Morex_V3.fasta \
        outdir=${DIR_REF}
```

## Download PRJNA558196 (Barakate et al 2021) anther RNAseq data
1. Download SRA sequences (Output: *.sra binary format)
```{sh}
time \
prefetch \
    PRJNA558196 \
    -O ${DIR}/SRA/
```

2. Extract SRA sequence identities (Output: new_names.tmp)
```{sh}
touch new_names.tmp
for f in $(find ${DIR}/SRA -name '*.sra' | sort)
do
    base_name=$(basename $f)
    srr_name=$(echo ${base_name%.sra*})
    new_name=$(esearch -db sra -query $srr_name | \
            esummary | \
            xtract -pattern DocumentSummary -block Library_descriptor -element LIBRARY_NAME | \
            sed 's/ /_/g')
    echo $new_name >> new_names.tmp
done
```

3. Convert SRA into FASTQ sequences, rename, and clean-up (Output: *.fastq.gz)
```{sh}
echo '#!/bin/bash
f=$1
new_name=$2
DIR_FASTQ=$3
fastq-dump \
    $f \
    --accession ${new_name} \
    --split-3 \
    --gzip \
    -O ${DIR_FASTQ}
rm $f
' > extract_RNAseq.sh
chmod +x extract_RNAseq.sh
time \
parallel --link \
./extract_RNAseq.sh \
    {1} \
    {2} \
    ${DIR_FASTQ} \
    ::: $(find ${DIR}/SRA -name '*.sra' | sort) \
    ::: $(cat new_names.tmp)
rm new_names.tmp
rm -R SRA/
```

## Map the RNAseq to the Golden Promise reference genome
1. Build the reference index: OUTPUT: ${REF%.fasta*}
```{sh}
time \
hisat2-build \
    ${REF} \
    ${REF%.fasta*}
```

2. Align in parallel
```{sh}
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

## Assemble, identify transcripts, and estimate transcripts per million (TPM)
```{sh}
time \
TPMCalculator \
    -g GTF_file [-d BAM_files_directory|-b BAM_file]
```