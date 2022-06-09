# Search for the exclusively-meiocyte-ProphaseI-specifically expressed BTR1-like and BTR2-like genes in GoldenPromise Barley

## Setup working directories
```{sh}
DIR=/data-weedomics-3/BTR-like-UTRs/search_test.md
DIR_SRA=${DIR}/SRA ### Assumes the reads are free of adapter sequences
DIR_FASTQ=${DIR}/FASTQ ### Assumes the reads are free of adapter sequences
DIR_REF=${DIR}/REF
DIR_BAM=${DIR}/BAM
DIR_GTF=${DIR}/GTF
REF=${DIR_REF}/GoldenPromise.fasta
PATH=${PATH}:${DIR}/mmseqs/bin
PATH=${PATH}:${DIR}
PATH=${PATH}:${DIR}/stringtie
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

git clone https://github.com/gpertea/stringtie.git
cd stringtie
make release
cd ${DIR}
```

## Map the Morex v3 annotations into the GoldenPromise genome assembly
1. Dowload the genome and annotation files of Morex V3 the representative reference genome of Hordeum vulgare
```{sh}
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/904/849/725/GCF_904849725.1_MorexV3_pseudomolecules_assembly/GCF_904849725.1_MorexV3_pseudomolecules_assembly_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/904/849/725/GCF_904849725.1_MorexV3_pseudomolecules_assembly/GCF_904849725.1_MorexV3_pseudomolecules_assembly_genomic.gtf.gz
gunzip -c GCF_904849725.1_MorexV3_pseudomolecules_assembly_genomic.fna.gz > ${DIR_REF}/Morex_V3.fasta
gunzip -c GCF_904849725.1_MorexV3_pseudomolecules_assembly_genomic.gtf.gz > ${DIR_REF}/Morex_V3.gtf
rm GCF_904849725.1_MorexV3_pseudomolecules_assembly_genomic.fna.gz
rm GCF_904849725.1_MorexV3_pseudomolecules_assembly_genomic.gtf.gz
```

2. Map the genome annotations of Morex V3 into the GoldenPromise genome assembly
```{sh}
time \
java -jar -Xmx280G GeMoMa-1.8.jar CLI \
        GeMoMaPipeline \
        threads=32 \
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

4. Clean-up: remove the m54053_* sequences
```{sh}
rm ${DIR_FASTQ}/m54053_*
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
BAMOUT=${R1%.fastq*}.bam
hisat2 \
    -x ${REF%.fasta*} \
    -1 ${R1} \
    -2 ${R2} \
    -S ${SAMOUT}
samtools \
    sort ${SAMOUT} \
    -O BAM \
    > ${BAMOUT}
' > map_RNAseq.sh
chmod +x map_RNAseq.sh
time \
parallel --link \
./map_RNAseq.sh \
    ${REF} \
    {1} \
    {2} \
    ::: $(ls ${DIR_FASTQ}/*_1.fastq.gz | sort) \
    ::: $(ls ${DIR_FASTQ}/*_2.fastq.gz | sort)
mv ${DIR_FASTQ}/*.bam ${DIR_BAM}
```

## Assemble, identify transcripts, and estimate transcripts per million (TPM)
```{sh}
echo '#!/bin/bash
BAM=$1
GTF=$2
stringtie \
    -G ${GTF} \
    -o ${BAM}.out.gtf \
    ${BAM}
' > stringtie_run.sh
chmod +x stringtie_run.sh
# time \
# parallel \
# ./stringtie_run.sh \
#     {} \
#     ${DIR_REF}/GeMoMa_annotations/final_annotation.gff \
#     ::: $(find $DIR_BAM -name '*.bam')
### Using the mapped Morex V2 annotations into the GoldenPromise V1 genome
time \
parallel \
./stringtie_run.sh \
    {} \
    ${DIR_REF}/Tritex.gff3 \
    ::: $(find $DIR_BAM -name '*.bam')

mv ${DIR_BAM}/*.gtf $DIR_GTF
```

## Mapping the positions of the BTR genes used in Cross et al, 2022
```{sh}
BTR1=HORVU.MOREX.r2.3HG0195510
BTR1_LIKE_a=HORVU.MOREX.r2.3HG0195460
BTR1_LIKE_b1=HORVU.MOREX.r2.3HG0195170
BTR2='KR813335.1(OUH602)'
BTR2_LIKE_a=HORVU.MOREX.r2.3HG0195480
BTR2_LIKE_b1=HORVU.MOREX.r2.3HG0195160
BTR2_LIKE_b2=HORVU.MOREX.r2.3HG0195470
BTR1_LIKE=HORVU.MOREX.r2.1HG0041780
BTR2_LIKE=HORVU.MOREX.r2.2HG0104610
BTR1_2_LIKE_1=HORVU.MOREX.r2.2HG0086890
BTR1_2_LIKE_2=HORVU.MOREX.r2.2HG0109560
BTR1_2_LIKE_3=HORVU.MOREX.r2.6HG0456990


for f in $(find $DIR_GTF -name '*.bam.out.gtf')
do
    head -n2 $f > ${f%.bam*}.BTRs.gtf
    for q in "BTR1" "BTR1_LIKE_a" "BTR1_LIKE_b1" "BTR2" "BTR2_LIKE_a" "BTR2_LIKE_b1" "BTR2_LIKE_b2" "BTR1_LIKE" "BTR2_LIKE" "BTR1_2_LIKE_1" "BTR1_2_LIKE_2" "BTR1_2_LIKE_3"
    do
        # q="BTR2_LIKE_a"
        # f=$(find $DIR_GTF -name '*.gtf' | head -n1)
        grep "${!q}" $f | grep "FPKM" | sed "s/${!q}.1.mrna1/$q/g" >> ${f%.bam*}.BTRs.gtf
    done
done

```

## Gene expression analysis to identify overexpressed transcripts
1. Extract TPM and its associated transcript and gene IDs
```{sh}
echo '# ARGS = ["/data/weedomics/misc/BTR-like_genes_barley/GTF/A-LepZyg1_1.bam.out.gtf"]
fname_gtf = ARGS[1]
sample_name = split(basename(fname_gtf), "_")[1][1:(end-1)]
replicate = split(basename(fname_gtf), "_")[1][end]
fname_out = string(fname_gtf, ".csv")

file_out = open(fname_out, "a")
write(file_out, string("sample,rep,chr,ini,fin,transcript,gene,tpm\n"))

file = open(fname_gtf, "r")
while !eof(file)
    line = readline(file)
    if line[1] == "#"[1]
        continue
    end
    tab_delim = split(line, "\t"[1])
    if tab_delim[3] == "transcript"
        semicol_delim = split(tab_delim[end], ";"[1])
        chr_name = tab_delim[1]
        ini_pos = parse(Int, tab_delim[4])
        fin_pos = parse(Int, tab_delim[5])
        transcript_name = split(semicol_delim[match.(Regex("transcript_id"), semicol_delim) .!= nothing][1], "\""[1])[2]
        gene_name = try
                split(semicol_delim[match.(Regex("reference_id"), semicol_delim) .!= nothing][1], "\""[1])[2]
            catch
                "None"
            end
        gene_name = replace(gene_name, "Morex_V3_rna-"=>"")
        gene_name = replace(gene_name, "_R0"=>"")
        tpm = parse(Float64, split(semicol_delim[match.(Regex("TPM"), semicol_delim) .!= nothing][1], "\""[1])[2])
        write(file_out, string(join([sample_name, replicate, chr_name, ini_pos, fin_pos, transcript_name, gene_name, tpm], ","), "\n"))
    end
end
close(file)
close(file_out)
' > extract_TPM_per_gene.jl
time \
parallel \
julia extract_TPM_per_gene.jl \
    {} ::: $(find $DIR_GTF -name '*_1.BTRs.gtf')
```

2. Merge TPMs
```{sh}
head -n1 $(find ${DIR_GTF} -name '*_1.BTRs.gtf.csv' | head -n1) > ${DIR}/BTR_TPM.csv
for f in $(find ${DIR_GTF} -name '*_1.BTRs.gtf.csv')
do
    tail -n+2 $f >> ${DIR}/BTR_TPM.csv
done
```

3. Analyse
```{R}
dat = read.csv("BTR_TPM.csv")
dat$sample = as.factor(dat$sample)
dat$rep = as.factor(dat$rep)
dat$chr = as.factor(dat$chr)
dat$gene = as.factor(dat$gene)

for (gene in levels(dat$gene)){
    # gene = levels(dat$gene)[6]
    subdf = droplevels(dat[dat$gene==gene, ])
    X = dat$TPM
}

```


## Extract transcript sequences
```{sh}
time \
gffread \
    -g ${REF} \
    ${DIR_GTF}/M-PachDipl1_1.BTRs.gtf \
    -w M-PachDipl1_1.BTRs.gffread.fasta
```