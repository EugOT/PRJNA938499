#!/usr/bin/env bash

# esearch -db sra -query PRJNA938499 | efetch -format runinfo | cut -d ',' -f 1 | grep SRR | xargs ffq --ftp | jq -r '.[] | .url'
cd fastq && \
esearch -db sra -query PRJNA938499 | efetch -format runinfo | cut -d ',' -f 1 | grep SRR | xargs -n 1 -P 20 prefetch --max-size u && \
esearch -db sra -query PRJNA938499 | efetch -format runinfo | cut -d ',' -f 1 | grep SRR | xargs -n 1 -P 20 fasterq-dump -p -x --threads 10 --mem 20000M --outdir fastq --split-files --include-technical && pigz -p 20 fastq/SRR236150*.fastq && \
mv fastq/SRR23615084_1.fastq.gz fastq/SRR23615084_S1_L001_I1_001.fastq.gz && mv fastq/SRR23615084_2.fastq.gz fastq/SRR23615084_S1_L001_R1_001.fastq.gz && mv fastq/SRR23615084_3.fastq.gz fastq/SRR23615084_S1_L001_R2_001.fastq.gz && \
mv fastq/SRR23615083_1.fastq.gz fastq/SRR23615083_S1_L001_I1_001.fastq.gz && mv fastq/SRR23615083_2.fastq.gz fastq/SRR23615083_S1_L001_R1_001.fastq.gz && mv fastq/SRR23615083_3.fastq.gz fastq/SRR23615083_S1_L001_R2_001.fastq.gz && \
mv fastq/SRR23615082_1.fastq.gz fastq/SRR23615082_S1_L001_I1_001.fastq.gz && mv fastq/SRR23615082_2.fastq.gz fastq/SRR23615082_S1_L001_R1_001.fastq.gz && mv fastq/SRR23615082_3.fastq.gz fastq/SRR23615082_S1_L001_R2_001.fastq.gz && \
mv fastq/SRR23615081_1.fastq.gz fastq/SRR23615081_S1_L001_I1_001.fastq.gz && mv fastq/SRR23615081_2.fastq.gz fastq/SRR23615081_S1_L001_R1_001.fastq.gz && mv fastq/SRR23615081_3.fastq.gz fastq/SRR23615081_S1_L001_R2_001.fastq.gz && \
mv fastq/SRR23615080_1.fastq.gz fastq/SRR23615080_S1_L001_I1_001.fastq.gz && mv fastq/SRR23615080_2.fastq.gz fastq/SRR23615080_S1_L001_R1_001.fastq.gz && mv fastq/SRR23615080_3.fastq.gz fastq/SRR23615080_S1_L001_R2_001.fastq.gz && \
mv fastq/SRR23615079_1.fastq.gz fastq/SRR23615079_S1_L001_I1_001.fastq.gz && mv fastq/SRR23615079_2.fastq.gz fastq/SRR23615079_S1_L001_R1_001.fastq.gz && mv fastq/SRR23615079_3.fastq.gz fastq/SRR23615079_S1_L001_R2_001.fastq.gz && \
mv fastq/SRR23615078_1.fastq.gz fastq/SRR23615078_S1_L001_I1_001.fastq.gz && mv fastq/SRR23615078_2.fastq.gz fastq/SRR23615078_S1_L001_R1_001.fastq.gz && mv fastq/SRR23615078_3.fastq.gz fastq/SRR23615078_S1_L001_R2_001.fastq.gz && \
mv fastq/SRR23615077_1.fastq.gz fastq/SRR23615077_S1_L001_I1_001.fastq.gz && mv fastq/SRR23615077_2.fastq.gz fastq/SRR23615077_S1_L001_R1_001.fastq.gz && mv fastq/SRR23615077_3.fastq.gz fastq/SRR23615077_S1_L001_R2_001.fastq.gz