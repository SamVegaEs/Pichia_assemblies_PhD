```bash
  OutDir=genome_alignment/circos/Y-11454_v2_cov_plot
  mkdir -p $OutDir

  Ref_genome=$(ls assembly/misc_publications/P.stipitis/Y-11545_v2/pichia.allmasked)
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/fasta2circos.py --genome $Ref_genome --contig_prefix "Y-11545_" > $OutDir/ref_genome.txt
```


Calculate coverage over 10kb windows.
Parse read depth files:

```bash
OutDir=genome_alignment/circos/Y-11454_v2_cov_plot
for File in $(ls analysis/genome_alignment/bwa/P.stipitis/*/vs_assembly/misc_publications/P.stipitis/Y-11545_v2/pichia.AssembledScaffolds.fasta_unmasked/P.stipitis_P.stipitis_*_vs_genome_alignment_depth_10kb.tsv); do
  Strain=$(echo $File | cut -f5 -d '/')
  echo $Strain
  cat $File | awk '{print $1,$2-1000,$2,$3,$4}' OFS='\t' | cut -f1,2,3,4 | sed 's/chr_/Y-11545_chr_/g' > $OutDir/${Strain}_vs_ref_unmasked_scatterplot.tsv
done
for File in $(ls analysis/genome_alignment/bwa/P.stipitis/*/vs_assembly/misc_publications/P.stipitis/Y-11545_v2/pichia.AssembledScaffolds.fasta_unmasked/P.stipitis_P.stipitis_*_vs_genome_alignment_depth_10kb.tsv); do
  Strain=$(echo $File | cut -f5 -d '/')
  echo $Strain
  cat $File | awk '{print $1,$2-999,$2,$3,$4}' OFS='\t' | cut -f1,2,3,4 | sed 's/chr_/Y-11545_chr_/g' > $OutDir/${Strain}_vs_ref_masked_scatterplot.tsv
done
````

```bash
ProgDir=/home/armita/git_repos/emr_repos/scripts/pichia/genome_alignment/vs_ref/coverage_plot
circos -conf $ProgDir/Y-11454_v2_circos.conf -outputdir $OutDir
mv $OutDir/circos.png $OutDir/Y-11545_v2_cov_unmasked_circos.png
mv $OutDir/circos.svg $OutDir/Y-11545_v2_cov_unmasked_circos.svg
# mv $OutDir/circos.png $OutDir/Y-11545_v2_cov_masked_circos.png
# mv $OutDir/circos.svg $OutDir/Y-11545_v2_cov_masked_circos.svg
ls $PWD/$OutDir/Y-11545_v2*_circos.png
```


# Duplicated region?

Chromosome 5 of isolate 591 showed elevated read coverage at 5.0Mb.

The region with elevated coverage was identified:

```bash
File=$(ls analysis/genome_alignment/bwa/P.stipitis/591/vs_assembly/misc_publications/P.stipitis/Y-11545_v2/pichia.AssembledScaffolds.fasta_unmasked/P.stipitis_P.stipitis_591_vs_genome_alignment_depth_10kb.tsv)
cat $File | grep 'chr5.1' | less -S
```

```
chr_5.1 568000  61      P.stipitis_591
chr_5.1 569000  71      P.stipitis_591
chr_5.1 570000  322     P.stipitis_591
chr_5.1 571000  360     P.stipitis_591
chr_5.1 572000  372     P.stipitis_591
chr_5.1 573000  340     P.stipitis_591
chr_5.1 574000  347     P.stipitis_591
chr_5.1 575000  321     P.stipitis_591
chr_5.1 576000  346     P.stipitis_591
chr_5.1 577000  408     P.stipitis_591
chr_5.1 578000  325     P.stipitis_591
chr_5.1 579000  361     P.stipitis_591
chr_5.1 580000  342     P.stipitis_591
chr_5.1 581000  358     P.stipitis_591
chr_5.1 582000  347     P.stipitis_591
chr_5.1 583000  339     P.stipitis_591
chr_5.1 584000  354     P.stipitis_591
chr_5.1 585000  367     P.stipitis_591
chr_5.1 586000  366     P.stipitis_591
chr_5.1 587000  352     P.stipitis_591
chr_5.1 588000  401     P.stipitis_591
chr_5.1 589000  330     P.stipitis_591
chr_5.1 590000  384     P.stipitis_591
chr_5.1 591000  354     P.stipitis_591
chr_5.1 592000  388     P.stipitis_591
chr_5.1 593000  352     P.stipitis_591
chr_5.1 594000  329     P.stipitis_591
chr_5.1 595000  376     P.stipitis_591
chr_5.1 596000  394     P.stipitis_591
chr_5.1 597000  360     P.stipitis_591
chr_5.1 598000  407     P.stipitis_591
chr_5.1 599000  349     P.stipitis_591
chr_5.1 600000  352     P.stipitis_591
chr_5.1 601000  362     P.stipitis_591
chr_5.1 602000  394     P.stipitis_591
chr_5.1 603000  407     P.stipitis_591
chr_5.1 604000  412     P.stipitis_591
chr_5.1 605000  350     P.stipitis_591
chr_5.1 606000  352     P.stipitis_591
chr_5.1 607000  298     P.stipitis_591
chr_5.1 608000  301     P.stipitis_591
chr_5.1 609000  371     P.stipitis_591
chr_5.1 610000  354     P.stipitis_591
chr_5.1 611000  368     P.stipitis_591
chr_5.1 612000  387     P.stipitis_591
chr_5.1 613000  417     P.stipitis_591
chr_5.1 614000  352     P.stipitis_591
chr_5.1 615000  394     P.stipitis_591
chr_5.1 616000  379     P.stipitis_591
chr_5.1 617000  386     P.stipitis_591
chr_5.1 618000  353     P.stipitis_591
chr_5.1 619000  352     P.stipitis_591
chr_5.1 620000  329     P.stipitis_591
chr_5.1 621000  359     P.stipitis_591
chr_5.1 622000  349     P.stipitis_591
chr_5.1 623000  347     P.stipitis_591
chr_5.1 624000  372     P.stipitis_591
chr_5.1 625000  361     P.stipitis_591
chr_5.1 626000  364     P.stipitis_591
chr_5.1 627000  376     P.stipitis_591
chr_5.1 628000  366     P.stipitis_591
chr_5.1 629000  374     P.stipitis_591
chr_5.1 630000  371     P.stipitis_591
chr_5.1 631000  361     P.stipitis_591
chr_5.1 632000  392     P.stipitis_591
chr_5.1 633000  380     P.stipitis_591
chr_5.1 634000  365     P.stipitis_591
chr_5.1 635000  367     P.stipitis_591
chr_5.1 636000  400     P.stipitis_591
chr_5.1 637000  417     P.stipitis_591
chr_5.1 638000  420     P.stipitis_591
chr_5.1 639000  402     P.stipitis_591
chr_5.1 640000  382     P.stipitis_591
chr_5.1 641000  357     P.stipitis_591
chr_5.1 642000  357     P.stipitis_591
chr_5.1 643000  373     P.stipitis_591
chr_5.1 644000  422     P.stipitis_591
chr_5.1 645000  306     P.stipitis_591
chr_5.1 646000  411     P.stipitis_591
chr_5.1 647000  272     P.stipitis_591
chr_5.1 648000  121     P.stipitis_591
chr_5.1 649000  56      P.stipitis_591
chr_5.1 650000  38      P.stipitis_591
chr_5.1 651000  27      P.stipitis_591
chr_5.1 652000  54      P.stipitis_591
chr_5.1 653000  85      P.stipitis_591
chr_5.1 654000  78      P.stipitis_591
chr_5.1 655000  42      P.stipitis_591
chr_5.1 656000  88      P.stipitis_591
chr_5.1 657000  5       P.stipitis_591
chr_5.1 658000  6       P.stipitis_591
chr_5.1 659000  55      P.stipitis_591
chr_5.1 660000  76      P.stipitis_591
chr_5.1 661000  77      P.stipitis_591
```

Genes in this region (bp 569000-648000) were identified.

```bash
Gff=$(ls assembly/misc_publications/P.stipitis/Y-11545_v2/Pstipitisv2.FrozenGeneCatalog20080115.gff)
cat $Gff | grep 'chr_5.1' | grep -w -e "5......" -e "6....." | less
```

The first gene in this region was identified at bp 569801.
The final gene in this region stopped at 646810.

Genes in this region were identified:

```bash
# cat $Gff | grep 'chr_5.1' | grep -P -e "\t[56]\d\d\d\d\d\t" | grep -w 'exon' | grep -A 42 '569801'
cat $Gff | grep 'chr_5.1' | grep -P -e "\t[56]\d\d\d\d\d\t" | grep -w 'exon' | grep -A 42 '569801' | cut -f9 | sed "s/.*transcriptId //g" | uniq
cat $Gff | grep 'chr_5.1' | grep -P -e "\t[56]\d\d\d\d\d\t" | grep -w 'exon' | grep -A 42 '569801' | cut -f9 | sed "s/.*transcriptId //g" | uniq | wc -l
31
```

```
32054
89662
84281
90035
32059
61082
59814
89595
47494
67817
46396
83705
60417
22565
67819
32069
46113
65754
59994
67821
78208
32075
46516
46961
47531
83761
46807
32081
32082
32083
32085
```

```bash
Gff=$(ls assembly/misc_publications/P.stipitis/Y-11545_v2/Pstipitisv2.FrozenGeneCatalog20080115.gff)
OutDir=analysis/LS_regions/591_vs_ref
mkdir -p $OutDir
# cat $Gff | grep 'chr_5.1' | grep -P -e "\t[56]\d\d\d\d\d\t" | grep -w 'exon' | grep -A 42 '569801' | cut -f9 | sed "s/.*transcriptId //g" | uniq > tmp.txt
cat $Gff | grep 'chr_5.1' | grep -P -e "\t[56]\d\d\d\d\d\t" | grep -w 'exon' | grep -A 42 '569801' | cut -f9 | sed "s/.*transcriptId //g" | uniq > $OutDir/ref_genes_in_dup_591_region.txt

InterPro=$(ls assembly/misc_publications/P.stipitis/Y-11545_v2/Pstipitisv2.domaininfo_FrozenGeneCatalog20080115.tab)
# cat $InterPro | grep -w -f tmp.txt > ref_genes_in_dup_591_region_interpro.tsv
cat $InterPro | grep -w -f tmp.txt > $OutDir/ref_genes_in_dup_591_region_interpro.tsv

Kog=$(ls assembly/misc_publications/P.stipitis/Y-11545_v2/Pstipitisv2.koginfo_FrozenGeneCatalog20080115.tab)
# cat $Kog | grep -w -f tmp.txt > ref_genes_in_dup_591_region_kog.tsv
cat $Kog | grep -w -f tmp.txt > $OutDir/ref_genes_in_dup_591_region_kog.tsv

PathwayInfo=$(ls assembly/misc_publications/P.stipitis/Y-11545_v2/Pstipitisv2.ecpathwayinfo_FrozenGeneCatalog20080115.tab)
cat $PathwayInfo | grep -w -f tmp.txt > ref_genes_in_dup_591_region_pathwayinfo.tsv
cat $PathwayInfo | grep -w -f tmp.txt > $OutDir/ref_genes_in_dup_591_region_pathwayinfo.tsv
```
