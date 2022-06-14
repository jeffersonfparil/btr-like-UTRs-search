# btr-like-UTRs-search
Search for the exclusively-meiocyte-ProphaseI-expressed BTR1-like and BTR2-like genes in the "Golden Promise" barley cultivar

## Setup working directories
```{sh}
DIR=/data-weedomics-3/BTR-like-UTRs/search_test.md
DIR_SRA=${DIR}/SRA ### Assumes the reads are free of adapter sequences
DIR_FASTQ=${DIR}/FASTQ ### Assumes the reads are free of adapter sequences
DIR_TRANSCRIPTOMES=${DIR}/TRANSCRIPTOMES
PATH=${PATH}:${DIR}
MACSE=${DIR}/MACSE/macse_v2.06.jar
cd $DIR
```

## Install tools
```{sh}
cd $DIR
sudo apt install -y sra-toolkit \
                    ncbi-entrez-direct \
                    cmake \
                    bowtie2 \
                    salmon \
                    samtools \
                    ncbi-blast+ \
                    r-base-core \
                    seqtk \
                    emboss
wget https://github.com/trinityrnaseq/trinityrnaseq/releases/download/Trinity-v2.14.0/trinityrnaseq-v2.14.0.FULL.tar.gz
tar -xvzf trinityrnaseq-v2.14.0.FULL.tar.gz
cd trinityrnaseq-v2.14.0
make
PATH=${PATH}:$(pwd)
cd -
mkdir MACSE/
cd MACSE/
wget https://bioweb.supagro.inra.fr/macse/releases/macse_v2.06.jar
MACSE=${DIR}/MACSE/macse_v2.06.jar
java -Xmx250G -jar ${MACSE} -help
cd -
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

## Trinity transcriptome assembly
```{sh}
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
```{${DIR_REF}/GeMoMa_annotations/Btr_genes.fasta}
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
```{sh}
time \
for stage in M-LepZyg M-PachDipl
do
    DB=${DIR_TRANSCRIPTOMES}/trinity-${stage}.Trinity.fasta
    QUERY=${DIR_REF}/GeMoMa_annotations/Btr_genes.fasta
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
        -qcov_hsp_perc 0.90 \
        -outfmt "6 qseqid staxids pident evalue qcovhsp bitscore stitle" \
        -out ${DIR_TRANSCRIPTOMES}/trinity-${stage}.Trinity.blastout.tbl
done
```

2. Generate a short list by identifying the transcripts containing >=99% of the gene sequences (Output: ${DIR}/Btr_genes_transcript_hits_per_stage.txt)
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

3. Extract the Btr-like gene sequences and their corresponding top blast transcript hits (Output: all_sequences_genes_and_transcripts.fasta)
```{sh}
f=${DIR}/Btr_genes_transcript_hits_per_stage.txt
q=${DIR_REF}/GeMoMa_annotations/Btr_genes.fasta

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

4. Determine relationships between sequences by aligning them and building gene/transcript trees (Ouputs: [Btr1, Btr2]_like_genes_and_transcripts.aln.cds; and [Btr1, Btr2]_like_genes_and_transcripts.aln.cds.iqtree)
```{sh}
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
    java -Xmx8G \
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
