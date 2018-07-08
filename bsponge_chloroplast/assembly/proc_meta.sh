#!/bin/bash

for sample in 1747 1791
do
	sname=s$sample
	bowtie2 -x bchchlor -1 ${sname}_1.fq -2 ${sname}_2.fq -p 6 --no-unal -k 1| samtools view -Sbh -F u - >ch${sname}.bam
done 
for sample in 1747 1791 1820
do
	sname=d$sample
	bowtie2 -x bchchlor -1 ${sname}_1.fq -2 ${sname}_2.fq -p 6 --no-unal -k 1 | samtools view -Sbh -F u - >ch${sname}.bam
done 

for file in ch*.bam
do
	mfn=${file:2}
	sfn=${mfn%.bam}s
	samtools sort $file ${sfn}
	samtools index ${sfn}.bam
done

