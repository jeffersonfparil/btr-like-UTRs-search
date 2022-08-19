# btr-like-UTRs-search
Search for the exclusively-meiocyte-ProphaseI-expressed BTR1-like and BTR2-like genes in the "Golden Promise" barley cultivar

## Setup working directories
```shell
DIR=/data/weedomics/misc/BTR-like_genes_barley
# DIR=/data-weedomics-1/parilj/BTR_LIKE_GENES_IN_BARLEY
DIR_SRA=${DIR}/SRA ### Assumes the reads are free of adapter sequences
DIR_FASTQ=${DIR}/FASTQ ### Assumes the reads are free of adapter sequences
DIR_TRANSCRIPTOMES=${DIR}/TRANSCRIPTOMES
DIR_BAM=${DIR}/BAM
DIR_GTF=${DIR}/GTF
PATH=${PATH}:${DIR}/sratoolkit.3.0.0-ubuntu64/bin
PATH=${PATH}:${DIR}
MACSE=${DIR}/MACSE/macse_v2.06.jar
PATH=${PATH}:${DIR}/hisat2
PATH=${PATH}:${DIR}/stringtie
cd $DIR
```

## Install tools
```shell
cd $DIR
sudo apt install -y autoconf \
                    default-jre \
                    ncbi-entrez-direct \
                    cmake \
                    bowtie2 \
                    salmon \
                    samtools \
                    bamtools \
                    ncbi-blast+ \
                    r-base-core \
                    seqtk \
                    emboss
wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar -vxzf sratoolkit.tar.gz
PATH=${PATH}:${DIR}/sratoolkit.3.0.0-ubuntu64/bin
vdb-config --interactive ### Then choose ${DIR} as the location of the user-repository
wget https://github.com/trinityrnaseq/trinityrnaseq/releases/download/Trinity-v2.14.0/trinityrnaseq-v2.14.0.FULL.tar.gz
tar -xvzf trinityrnaseq-v2.14.0.FULL.tar.gz
cd trinityrnaseq-v2.14.0
PATH=${PATH}:$(pwd)
Trinity -h
cd -
mkdir MACSE/
cd MACSE/
wget https://bioweb.supagro.inra.fr/macse/releases/macse_v2.06.jar
MACSE=${DIR}/MACSE/macse_v2.06.jar
java -Xmx250G -jar ${MACSE} -help
cd -
git clone https://github.com/DaehwanKimLab/hisat2.git
cd hisat2
make
PATH=${PATH}:${DIR}/hisat2
sudo ln -s /usr/bin/python3 /usr/bin/python ### hisat2 needs python and python2 is now fully deprecated and we may have to specify that the default python is python3
hisat2 -h
cd -
git clone https://github.com/gpertea/stringtie.git
cd stringtie
make release
stringtie -h
cd ${DIR}
```

## Download PRJNA558196 (Barakate et al 2021) anther RNAseq data
1. Download SRA sequences (Output: *.sra binary format)
```shell
time \
prefetch \
    PRJNA558196 \
    -O ${DIR}/SRA/
```

2. Extract SRA sequence identities (Output: new_names.tmp)
```shell
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
```shell
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
```

4. Clean-up
```shell
rm new_names.tmp
rm -R ${DIR}/SRA
rm -R ${DIR}/PRJNA558196
rm ${DIR_FASTQ}/m54053_*
```

## Trinity transcriptome assembly
```shell
cd $DIR
DIR_TRANSCRIPTOMES=${DIR}/TRANSCRIPTOMES
mkdir $DIR_TRANSCRIPTOMES

echo '#!/bin/bash
stage=$1
rep=$2
DIR_FASTQ=$3
DIR_TRANSCRIPTOMES=$4
# stage=M-PachDipl
# rep=1
zcat ${DIR_FASTQ}/${stage}*_${rep}.fastq.gz | gzip > ${DIR_TRANSCRIPTOMES}/${stage}_R${rep}.fastq.gz
' > concatenate_RNAseq_reps.sh
chmod +x concatenate_RNAseq_reps.sh
time \
parallel \
./concatenate_RNAseq_reps.sh \
    {1} \
    {2} \
    ${DIR_FASTQ} \
    ${DIR_TRANSCRIPTOMES} \
    ::: M-LepZyg M-PachDipl \
    ::: 1 2

time \
for stage in M-LepZyg M-PachDipl
do
# stage=M-PachDipl
Trinity \
    --seqType fq \
    --max_memory 200G \
    --trimmomatic \
    --left  ${DIR_TRANSCRIPTOMES}/${stage}_R1.fastq.gz \
    --right ${DIR_TRANSCRIPTOMES}/${stage}_R2.fastq.gz \
    --CPU 32 \
    --output ${DIR_TRANSCRIPTOMES}/trinity-${stage}
done
```


## Btr gene sequences
```{${DIR}/Btr_genes.fasta}
>HORVU.MOREX.r2.2HG0104610.1(Btr2-like-extra)
ATGGGCAAGCTCCTCTGCGACTCCTCCTCCGGCGGCGCCGCCGTCGCCGTTGCCGAGGCTCTCCTGCCCT
CCCCCGCCCCGCCCTCGCCGCCCAGCGCCGCCGCCTGCTCCGGATCTGGGCGCTCGACGACCAGCAGCGG
CTGGTCCGAGCGCCTGGCGCTCCACTGCTGCGGCCCCATCCCGTCGCTGCCGGTGGCCCGCATCGCCGAC
CTCCGCGACAACCACCACGACCGCCATGACGAGCGGCTCGCTCTGAACAGGCTCGAGGATGCGAGGGACT
TCGCCAAGCGCGCGCTCCGTGGGGTGGATGGAGCCCTCAAGCTCCTGGGCTCCGTCCAGTACATGCTTCA
CGACCTCGGCGCCGGCGCGGCCGGGGGTAGGCAGGCGATGGAAGAGCAGCTCCAGGCCGCCGCCCGCAAG
CTCCAGCTCGTGGCGGTCAGCACGTGCAACACGCGCTCGCTGGCCCGCATGGCCACCGAGCCTCCCATCG
GCAACCGCGTCCAGTGA
>HORVU.MOREX.r2.3HG0195160.1(Btr2-like-b1)
ATGGAGGAGTGGAGGAATCTGGCGTTGGCGGCGGCGGGTGAAAGCTTAGCCAACACCATGACCCGCCACA
TGGGTGTAGCGGAGGCCATCGCTCGCGCCGGCCAGCGGTACCGCCTGGCCGCCGAGGAGTGTCGCGGATT
CGGCCAGGGCGTGCACCCGACGCCCAACGCCGGCGAGCGCGCTTCAGCAGGCGGCGACTTCGTCGACCTC
GCCATCGACCGGATCAAGAGCATCAGCAGGTTCCACGCCGTGCGGGGCAGCGTCTTCTCCCTCTGCGTCC
GGCGTATCGGGCTCCAGGGCATCGGGCTCCACGGCGACGAGGCCCGCCACGCGGAGACGGCGATGCTGTC
GCTGCGTTCCGCCAAGTCGCACGGCCATGCGGCCCTCCGCGTCTTCGAGCGCATGCTCAGGCCGCCGTCG
CCGCAAGCGGTCGCCCGCGCGTGGGCGCCCGCGGCCGAGCAGCTCCTGCGCCGCGCGATCAACAATCTGG
ACATGGCGGCGGCCTCCGTGGGGCGGATACGCCCGGCCATCGTCGTCGAGTACAAGTACGCCCGGAAGCT
TCTGAATGCCTCCATATAG
>HORVU.MOREX.r2.3HG0195170.1(Btr1-like-b1)
ATGGCGCAGCCGGCAGGATGGAAGGCGATGTACCAGCAAGTGGTGATCGAGGCGGACGGCAGCTGCGCCG
ACGTCGAGCACAGAGTCGCCGCCGCGCGTACGGCGCTGGAGTCCCCGGAGGCGGTGCTGACCAGCCGCGA
CCCCACGGGGGTCTACACCTTGCTGAAGTCCGCGCTGGACGACGTCGAGCAAGCATCCGACTCCCTCTCC
GCCTTCATCATCCACGCGGTGGCGGCCGAGCGCCTGGCGCTCCACGGCTGCGGCGTCATCCCGTCGCAGC
CGGTGGCCCGCATCGCCGACCTCCGCGACGACCACCACGACCGCCACGACGAGCGGCTCGCTCTGAACAG
GCTCGAGGACGCCAGGGACTTCGCCAAGCGCGCGCTCCGCGGGGTGGATGGAGCCCTCAAGCTCCTGGGC
TCCGTCCAGTACATGCTTCGCGACCTCGGCGCCGGCGCGGCAGGGCGCAGGCAGGCGATGAAAGAGCAGC
TCCAGGCCGCCGCCCGCGAGCTCCAGCTCGTGGCGGTCAGCGTGTGCAACACGCGCTCGCTGGCCCGCAT
GGCCACCGAGCCTCCCATCGGCAACCGCGTCCAGTGA
>HORVU.MOREX.r2.3HG0195460.1(Btr1-like-a)
ATGGCGCAGCCGGCAGGATGGAAGGCGATGTACCAGCAAGTGGTGATCGAGGCGGACGGCAGCTGCGCCG
ACGTCGAGCACAGAGTCGCCGCCGCGCGTACGGCGCTGGAGTCCCCGGAGGCGGTGCTGACCTCCCGCGA
CCCCACGGGGGTCTACACCTTGCTGAAGTCCGCGCTGGACAACGTCGAGCAAGCATCCGACTCCCTCTCC
GCCTTCATCATCCACGCGGTGGCGGCCGAGCGCCTGGCGCTCCACGGCTGCGGCCCCATCCCGTCGCAGC
CGGTGGCCCGCATCGCCGACCTCCGCGACGACCACCACGACCGCCACGACGAGCGGCTCGCTCTGAACAG
GCTCGAGGACGCCAGGGACTTCGCCAAGCGCGCGCTCCGCGGGGTGGATGGAGCCCTCAAGCTCCTGGGC
TCCGTCCAGTACATGCTTCGCGACCTCGGCGCCGGCGCGGCAGGGCGCAGGCAGGCGATGAAAGAGCAGC
TCCAGGCCGCCGCCCGCGAGCTCCAGCTCGTGGCGGTCAGCGTGTGCAACACGCGCTCGCTGGCCCGCAT
GGCCACCGAGCCTCCCATCGGCCAGTGA
>HORVU.MOREX.r2.3HG0195470.1(Btr2-like-b2)
ATGGAGGAGTGGAGGAATCTGGCCTTGGCAGCGGCGGATGAAAGCTTCGCCAACACCGTGACCAATGGTG
TAGCGGAGACCATCGCTAGCGCCATCCAGCAGTACCGCCTGGCCGCCGAGGAGTGCCGCGGATTCGGCCA
GGGCGTGCACCCGACGCCCAACGCCGGCGAGCGCGCTTCAGCAGGCGGCGACTTCGTCGACCTCGCCATC
GACCGGATCAAGAGCATCAGCAGGTTCCACGCCGTGCGGGGCAGCGTCTTCTCGCTCTGCGTCCGGCGGA
TCGGGCTCCAGGGCAACGCGCTGTGGTACATGTGGCAGTTCTACCACGCCGACGAGGCCCGCCACGCGGA
GACGGCGATGCTGTCGCTGCGTTCCGCCAAGTCGCACGGCCATGCGGCCGTCCGCGTCTTCGAGCGCATG
CTCAGGCCGCCGTCGCCGCAAGCGGTCGCCCGCGCGTGGGCGCCCGCGGCCGTGCAGCTCCTGCGCCGCG
CAATCAAGAATCTGGCCATGGCGGAGGTCTCCGTGGGACAGATACGCCCGGCCATCGTCGTCGAGTACAA
CGACGCCCGGAGGCTTCTGCATGGCTGA
>HORVU.MOREX.r2.3HG0195480.1(Btr2-like-a)
ATGGCGGAGTGGATGAATCTGGCGTTGGCGGCGGCGTTTGACAGCTTCGCCTACACCGAGACCAATGGTG
TAGCGGAGGCCGTCGCTGGCGCCATCCAGCAGTACCGCCTGGCCGCCGAGGAGTGCCGCGGAATCGGCCA
GGGCGTGCACCCTACGCCCAACGCCGGCCAGGGCGCTTCAGCAGGCGGCGATTCCATCGACCTCGCCCTC
ACCCGGATCAAGAGCATCACCAGGTTCCACGCCGTGCGGGGTAGCGTCTTCTCCGTCTGCGTCCGCCGCA
TGGGGCTCCAGCCCGACACGCCGTGGCGGCTCCAGCACGCCACCGCGGCCCGCCACGCGGAGATGGCGAT
ACGGTGCCTGGGCACCGCCAAGTCGTACGGCCATGCGGCCCTCAGCGTCTTCCACCGCATGCTCAGGCCG
CCGTCGCCGCAAGCGGTCGCTCGCGCCTGGGCGCCCGCGGCCGAGCTGCTCCTGCGCCGCGCGATCGTCA
ATCTGGACATGGCGGAGGCCTCCGTGGGGAAGATACGCCCGGCCATCGGCGTCGAGTACAACGACGCCAG
GAGGCTTCTGCATGGCTGA
>HORVU.MOREX.r2.3HG0195510.1(Btr1)
ATGGCGCAGCCGCCGCAATGGAAGGCGATGTACCAGTATGTGACGCGACGGGCGCACGACGGCTGCGCCC
GCGTCGAGGAAAGCGTCGCCGCGGCGCGCGGAGCGCTGGCGACCCCGATGGTGCTGGACACCCGCGACGC
CGCGGGGCGGTGCACGTTGCTGCATTCCGCGGTGACCCACGTCGAGCACGCATCCGACTGCCTCTCCGGT
TTCATAGTCAGCGTGGTGGTGGCGGAGCTCCTGGTGCTCCATGGCTGCGGGGCCGTCCCGTCGAGGCCGG
TGGCCAGCATCGACGGCCTCCGCCGCAACCGCGACGACCACGACGAGTGGCTCGCTCTGAGCAGGCTCGA
GGCCGCCAGGGAGCACGGCCAGGACGCGCTCCGCGGAGTGGAGGGGGCCTTCACCCTCCTGGCCTCCGTC
CGGTTCATGCTTCGCAGCCGGACCCCCGACGCCGCCGGGCGCCGGCAAGCCATGGAAGAGCAGCTCCACG
CCGCCGCCGTCGAACTTCAGGCCGTGGTGGGCAGCGTGGCGAACATGTCCGCGCTGGCTTTCTTGGCCAC
TCAGCCTGCCATCCGCAACCGCATCCAGTGA
```

## Identify transcripts coding for the Btr-like genes
1. Find the most similar transcripts to the Btr-like genes
```shell
QUERY=${DIR}/Btr_genes.fasta

time \
for stage in M-LepZyg M-PachDipl
do
    DB=${DIR_TRANSCRIPTOMES}/trinity-${stage}.Trinity.fasta
    makeblastdb -in ${DB} \
                -dbtype nucl
    blastn -db ${DB} \
        -query ${QUERY} \
        -perc_identity 90 \
        -qcov_hsp_perc 0.90 \
        -out ${DIR_TRANSCRIPTOMES}/trinity-${stage}.Trinity.blastout
    blastn -db ${DB} \
        -query ${QUERY} \
        -perc_identity 90 \
        -qcov_hsp_perc 90 \
        -outfmt "6 qseqid staxids pident evalue qcovhsp bitscore stitle" \
        -out ${DIR_TRANSCRIPTOMES}/trinity-${stage}.Trinity.blastout.tbl
done
```

2. Generate a short list by identifying the transcripts containing >=99% of the gene sequences (Output: `${DIR}/Btr_genes_transcript_hits_per_stage.txt`)
```{R}
tryCatch(rm(OUT), error=function(e){print("First run")})
DIR_TRANSCRIPTOMES = "TRANSCRIPTOMES"
for (stage in c("M-LepZyg", "M-PachDipl")){
    # stage = "M-LepZyg"
    f=file.path(DIR_TRANSCRIPTOMES, paste0("trinity-", stage, ".Trinity.blastout.tbl"))
    dat = read.delim(f, header=FALSE)
    colnames(dat) = c("qseqid", "staxids", "pident", "evalue", "qcovhsp", "bitscore", "stitle")
    dat$transcript = as.factor(unlist(lapply(strsplit(as.character(dat$stitle), " "), FUN=function(x){x[1]})))
    dat$length = as.numeric(gsub("len=", "", unlist(lapply(strsplit(as.character(dat$stitle), " "), FUN=function(x){x[2]}))))
    gene_list = levels(dat$qseqid)
    gene_list = gene_list[order(unlist(lapply(strsplit(gene_list, "\\("), FUN=function(x){x[2]})), decreasing=FALSE)]
    for (gene in gene_list){
        # gene = levels(dat$qseqid)[2]
        subdf = dat[dat$qseqid == gene, ]
        subdf = subdf[subdf$qcovhsp>=99, ]
        subdf = subdf[order(subdf$bitscore, decreasing=TRUE), ]
        subdf = subdf[order(subdf$bitscore, decreasing=TRUE), ]
        outdf = data.frame(stage=rep(stage, nrow(subdf)),
                           gene=subdf$qseqid,
                           transcript=subdf$transcript,
                           score=subdf$bitscore,
                           length=subdf$length)
        if (exists("OUT")==FALSE){
            OUT = outdf
        } else {
            OUT = rbind(OUT, outdf)
        }
    }
}
write.table(OUT, file="Btr_genes_transcript_hits_per_stage.txt", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
```

3. Extract the Btr-like gene sequences and their corresponding top blast transcript hits (Output: `all_sequences_genes_and_transcripts.fasta`)
```shell
f=${DIR}/Btr_genes_transcript_hits_per_stage.txt
q=${DIR}/Btr_genes.fasta

cat $q > all_sequences_genes_and_transcripts.fasta

time \
for i in $(seq 1 $(cat $f | wc -l))
do
    # i=1
    stage=$(head -n${i} ${f} | tail -n1 | cut -f1); g=${DIR_TRANSCRIPTOMES}/trinity-${stage}.Trinity.fasta
    gene=$(head -n${i} ${f} | tail -n1 | cut -f2)
    transcript=$(head -n${i} ${f} | tail -n1 | cut -f3)
    echo ${gene} > ${gene}.tmp
    echo ${transcript} > ${transcript}.tmp
    seqtk subseq ${q} ${gene}.tmp > ${gene}.fa.tmp
    seqtk subseq ${g} ${transcript}.tmp > ${transcript}.fa.tmp
    needle \
        ${transcript}.fa.tmp \
        ${gene}.fa.tmp \
        -gapopen 100.0 \
        -gapextend 1.0 \
        -outfile ${stage}---${gene}---${transcript}.aln
    score=$(grep "# Score: " ${stage}---${gene}---${transcript}.aln | cut -d' ' -f3 | cut -d'.' -f1)
    if [ $score -lt 100 ]
    then
        rm ${stage}---${gene}---${transcript}.aln
        mv ${transcript}.fa.tmp ${transcript}.fa.tmp.bk
        seqtk seq -r ${transcript}.fa.tmp.bk > ${transcript}.fa.tmp
        rm ${transcript}.fa.tmp.bk
        needle \
            ${transcript}.fa.tmp \
            ${gene}.fa.tmp \
            -gapopen 100.0 \
            -gapextend 1.0 \
            -outfile ${stage}---${gene}---${transcript}.aln
    fi
    cat ${transcript}.fa.tmp >> all_sequences_genes_and_transcripts.fasta
done
rm *.tmp
```

4. Determine relationships between sequences by aligning them and building gene/transcript trees (Ouputs: `[Btr1, Btr2]_like_genes_and_transcripts.aln.cds`; and `[Btr1, Btr2]_like_genes_and_transcripts.aln.cds.iqtree`)
```shell
for btr in Btr1 Btr2 all
do
    # btr=Btr1
    if [ $btr == "all" ]
    then
        grep "^>" all_sequences_genes_and_transcripts.fasta | \
        sed 's/>//g' > seq_list.tmp
        btr=all_Btr
    else
        grep "$btr" ${DIR}/Btr_genes_transcript_hits_per_stage.txt | \
        cut -f2 | \
        sort | \
        uniq > seq_list.tmp
        grep "$btr" ${DIR}/Btr_genes_transcript_hits_per_stage.txt | \
        cut -f3 | \
        sort | \
        uniq >> seq_list.tmp
    fi
    seqtk subseq \
        all_sequences_genes_and_transcripts.fasta \
        seq_list.tmp \
        > ${btr}_like_genes_and_transcripts.fasta
    java -Xmx250G \
        -jar ${MACSE} \
        -prog alignSequences \
        -seq ${btr}_like_genes_and_transcripts.fasta \
        -out_NT ${btr}_like_genes_and_transcripts.aln.cds.tmp \
        -out_AA ${btr}_like_genes_and_transcripts.aln.pro.tmp
    # Convert stop codons and frameshifts as "---" for compatibility with downstream tools
    java -Xmx250G \
        -jar ${MACSE} \
        -prog exportAlignment \
        -align ${btr}_like_genes_and_transcripts.aln.cds.tmp \
        -codonForFinalStop --- \
        -codonForInternalStop NNN \
        -codonForExternalFS --- \
        -codonForInternalFS --- \
        -out_NT ${btr}_like_genes_and_transcripts.aln.cds \
        -out_AA ${btr}_like_genes_and_transcripts.aln.pro
    iqtree -s ${btr}_like_genes_and_transcripts.aln.cds
    cat ${btr}_like_genes_and_transcripts.aln.cds.iqtree
    rm *.tmp
done
```

## Assess expression levels of the transcripts
1. Build the indexes of the transcriptomes
```shell
echo '#!/bin/bash
f=$1
hisat2-build \
    ${f} \
    ${f%.fasta*}
' > build_index_with_hisat2.sh
chmod +x build_index_with_hisat2.sh
time \
parallel \
    ./build_index_with_hisat2.sh \
    {} \
    ::: $(find ${DIR_TRANSCRIPTOMES} -name '*.fasta')
```

2. Align in parallel
```shell
echo '#!/bin/bash
DIR_TRANSCRIPTOMES=$1
R1=$2
R2=$3
SAMOUT=${R1%.fastq*}.sam
BAMOUT=${R1%.fastq*}.bam
REF=${DIR_TRANSCRIPTOMES}/trinity-M-$(echo $R1 | grep -o "PachDipl\|LepZyg").Trinity
hisat2 \
    -x ${REF} \
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
    ${DIR_TRANSCRIPTOMES} \
    {1} \
    {2} \
    ::: $(find ${DIR_FASTQ} -name '*_1.fastq.gz' | grep "M-PachDipl\|M-LepZyg" | sort) \
    ::: $(find ${DIR_FASTQ} -name '*_2.fastq.gz' | grep "M-PachDipl\|M-LepZyg" | sort)
mv ${DIR_FASTQ}/*.bam ${DIR_BAM}
rm ${DIR_FASTQ}/*.sam
```

3. Create an annotation file using the transcriptome using itself as the psuedo-reference genome
```shell
for f in $(find ${DIR_TRANSCRIPTOMES} -name 'trinity-M-*.Trinity.fasta')
do
    # f=${DIR_TRANSCRIPTOMES}/trinity-M-PachDipl.Trinity.fasta
    # f=${DIR_TRANSCRIPTOMES}/trinity-M-LepZyg.Trinity.fasta
    grep "^>" $f | cut -d' ' -f1 | sed 's/>//g' > names.tmp
    grep "^>" $f | cut -d' ' -f2 | sed 's/len=//g' > lengths.tmp
    printf 'MOCK\n%.0s' $(seq 1 $(cat names.tmp | wc -l)) > mock_column2.tmp
    printf 'mRNA\n%.0s' $(seq 1 $(cat names.tmp | wc -l)) > mock_column3.tmp
    printf '1\n%.0s' $(seq 1 $(cat names.tmp | wc -l)) > start.tmp
    printf '.\n%.0s' $(seq 1 $(cat names.tmp | wc -l)) > dots.tmp
    printf '+\n%.0s' $(seq 1 $(cat names.tmp | wc -l)) > pluses.tmp
    sed 's/TRINITY_/ID=TRINITY_/g' names.tmp > ids.tmp
    paste names.tmp \
        mock_column2.tmp \
        mock_column3.tmp \
        start.tmp \
        lengths.tmp \
        dots.tmp \
        pluses.tmp \
        dots.tmp \
        ids.tmp > ${f}.mock.gff
done
rm *.tmp
```

4. Stringtie expression level assessment
```shell
echo '#!/bin/bash
BAM=$1
DIR_TRANSCRIPTOMES=$2
GTF=${DIR_TRANSCRIPTOMES}/trinity-M-$(echo $BAM | grep -o "PachDipl\|LepZyg").Trinity.fasta.mock.gff
stringtie \
    -G ${GTF} \
    -o ${BAM}.stringtie.gtf \
    ${BAM}
' > stringtie_run.sh
chmod +x stringtie_run.sh
time \
parallel \
./stringtie_run.sh \
    {} \
    ${DIR_TRANSCRIPTOMES} \
    ::: $(find ${DIR_BAM} -name '*.bam')
mv ${DIR_BAM}/*.stringtie.gtf ${DIR_GTF}
```

5. Extract expression levels of the short-listed putative Btr-like transcripts (Output: `${DIR}/Btr_genes_transcript_hits_per_stage_with_TPM.txt`)
```{R}
fname_summary = "Btr_genes_transcript_hits_per_stage.txt"
dirname_gtf = "GTF"
vec_fnames_gtf = list.files(dirname_gtf)[grepl(".stringtie.gtf$", list.files("GTF"))]

dat = read.delim(fname_summary, header=FALSE)
colnames(dat) = c("stage", "Btr_like_gene", "transcript", "bitscore", "transcript_length")
dat$transcript = as.character(dat$transcript)
dat$TPM_r1 = NA
dat$TPM_r2 = NA
dat$TPM_r3 = NA
dat$TPM_mu = NA

for (f in vec_fnames_gtf){
    # f = vec_fnames_gtf[1]
    fsplit = unlist(strsplit(unlist(strsplit(f, "_"))[1], ""))
    stage = paste(fsplit[1:(length(fsplit)-1)], collapse="")
    replicate = fsplit[length(fsplit)]
    
    d = read.delim(file.path(dirname_gtf, f), header=FALSE)
    d = d[grep("TPM", d$V9), ]
    d = d[(d$V1 %in% unique(dat$transcript)), ]
    transcripts = d$V1
    TPM = unlist(lapply(strsplit(as.character(d$V9), "; "), FUN=function(x){as.numeric(gsub(";", "", gsub("TPM ", "", x[length(x)])))}))
    d = aggregate(TPM~transcripts, data=data.frame(transcripts, TPM), FUN=sum)
    for (i in 1:nrow(d)){
        # i = 1
        ti = d$transcripts[i]
        ni = d$TPM[i]
        eval(parse(text=paste0("dat$TPM_r", replicate, "[(dat$stage == stage) & (dat$transcript == ti)] = ni")))
    }
}

dat$TPM_mu = (dat$TPM_r1 + dat$TPM_r2 + dat$TPM_r3) / 3

write.table(dat, file="Btr_genes_transcript_hits_per_stage_with_TPM.txt", sep="\t", quote=FALSE, row.names=FALSE)
```

## Map the BTR-like transcripts into the Golden Promise genome assembly
```shell
mkdir REF/
cd REF/
wget https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_902500625.1?download=true&gzip=true && \
rm wget-log && \
mv 'GCA_902500625.1?download=true' GoldenPromise.fasta

cd ..
mkdir BLAST_REF/
makeblastdb -in REF/GoldenPromise.fasta \
                -dbtype nucl

sed 's/-//g' all_Btr_like_genes_and_transcripts.aln.cds \
    > BLAST_REF/BLAST_Btr_like_genes_and_transcripts.cds

blastn -db REF/GoldenPromise.fasta \
        -query BLAST_REF/BLAST_Btr_like_genes_and_transcripts.cds \
        -perc_identity 80 \
        -qcov_hsp_perc 80 \
        -outfmt "6 qseqid sstart send pident evalue qcovhsp bitscore stitle" \
        -out BLAST_REF/BLAST_Btr_like_genes_and_transcripts.tsv

blastn -db REF/GoldenPromise.fasta \
        -query BLAST_REF/BLAST_Btr_like_genes_and_transcripts.cds \
        -perc_identity 80 \
        -qcov_hsp_perc 80 \
        -out BLAST_REF/BLAST_Btr_like_genes_and_transcripts.txt
```