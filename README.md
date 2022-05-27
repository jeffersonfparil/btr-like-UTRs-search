# btr-like-UTRs-search
Search for UTRs of BTR-like genes

## Install gffread
```{sh}
git clone https://github.com/gpertea/gffread
cd gffread
make release
PATH=${PATH}:$(pwd)
cd -
```

## Note: no need to clean-up the Stringtie2-derived annotations by merging redundant transcripts into loci it's already clean after Stringtie2 it seems
```{sh}
# time \
# for SEQ in SRR9889998 SRR9890000
# do
# gffread \
#     --merge \
#     -T \
#     ${SEQ}.hisat.sorted.bam.gtf \
#     -o ${SEQ}.gffread.gtf
# done
```

## Extract transcript sequences from Stringtie2-derived annotations and the reference genome
```{sh}
time \
for SEQ in SRR9889998 SRR9890000
do
gffread \
    -g 180720_golden_promise_pseudomolecules_v1.fasta \
    ${SEQ}.hisat.sorted.bam.gtf \
    -w ${SEQ}.gffread.fasta
done
```

## Install NCBI-BLAST, EMBOSS (specifically needle for global pairwise sequence alignment) and seqtk on ubuntu
```{sh}
sudo apt -y install ncbi-blast+

git clone https://github.com/lh3/seqtk.git;
cd seqtk
make
PATH=${PATH}:$(pwd)
cd -

wget http://emboss.sourceforge.net/apps/release/6.6/emboss/apps/stretcher.html

```

## Set the extracted transcript sequences as BLAST databases
```{sh}
time \
for SEQ in SRR9889998 SRR9890000
do
makeblastdb -in ${SEQ}.gffread.fasta \
            -title "${SEQ}.gffread.blastdb" \
            -dbtype nucl
done
```

## Identify the BTR-like genes as the top hits in BLAST
```{sh}
time \
for SEQ in SRR9889998 SRR9890000
do
blastn -db ${SEQ}.gffread.fasta \
      -query Btr-like_ref_seq_.fa \
      -perc_identity 50 \
      -qcov_hsp_perc 0.90 \
      -max_hsps 5 \
      -outfmt "6 qseqid staxids pident evalue qcovhsp bitscore stitle" \
      -out ${SEQ}.blastout
done
```

## Extract transcript names of top hits for Btr-like genes
Note: Just naively taking the top blast hit. Will need to assess the other blast hits.
```{sh}
time \
for SEQ in SRR9889998 SRR9890000
do
    for BTR in "Btr1-like-a" "Btr2-like-a"
    do
        NAME=$(grep "${BTR}" ${SEQ}.blastout | head -n1 | cut -f7 )
        echo $NAME > ${SEQ}-${BTR}.tmp

        seqtk subseq ${SEQ}.gffread.fasta ${SEQ}-${BTR}.tmp > ${SEQ}-${BTR}.fasta.tmp
        sed -i "s/>/>${SEQ}-${BTR}-/g" ${SEQ}-${BTR}.fasta.tmp

        # extract BTR seq
        grep "^>" Btr-like_ref_seq_.fa | grep "${BTR}" | sed 's/>//g' > extract_BTR_seq.tmp
        seqtk subseq Btr-like_ref_seq_.fa extract_BTR_seq.tmp > ${BTR}.fasta.tmp

        # align BTR gene within the transcript seq
        needle \
            ${SEQ}-${BTR}.fasta.tmp \
            ${BTR}.fasta.tmp \
            -gapopen 10.0 \
            -gapextend 0.5 \
            -outfile ${SEQ}-${BTR}.aln.tmp
        # debugging
        echo ${SEQ}-${BTR}-${NAME}
        grep "${BTR}" ${SEQ}.blastout
    done
done



```

