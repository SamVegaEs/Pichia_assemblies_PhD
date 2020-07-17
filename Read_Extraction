Read extraction by position.

```bash
#Run at cd /projects/oldhome/groups/harrisonlab/project_files/Pichia/analysis/genome_alignment/minimap/589/vs_591. The -h option was added, because when sorting a message error appeared saying that the header was missing. 

samtools view 589_contigs_unmasked.fa_aligned_sorted.bam "contig_4:1880600-1880700" -h > 591_C3_subtel_region_100bp.bam

#Sorting of the bam file generated. 

samtools sort -n 591_C3_subtel_region_100bp.bam > 591_C3_subtel_region_100bp_sorted.bam

#Now that we have the reads in a bam file, we can extract the fastq reads using bedtools.

bedtools bamtofastq -i 591_C3_subtel_region_100bp_sorted.bam -fq 591_C3_subtel_region_100bp_sorted.bam.fq
```

Read extraction by sequence.

```bash
samtools view 589_contigs_unmasked.fa_aligned_sorted.bam | awk '$10 ~/CCGCAAGACGTGAAAAGATCCATA/' | awk '{print $1, $10}' > telseq.fa

cat telseq2.fa | grep "9aecc4bd" > 150.fa
```
