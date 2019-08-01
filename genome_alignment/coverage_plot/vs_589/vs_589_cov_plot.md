  
```bash
 OutDir=vs_589_cov_plot
  mkdir -p $OutDir

  Ref_genome=$(ls repeat_masked/P.stipitis/589/filtered_contigs/*_contigs_hardmasked.fa)
  ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/fasta2circos.py --genome $Ref_genome --contig_prefix "589_" > $OutDir/ref_genome.txt
  ```

  # Calculate coverage over 10kb windows. Parse read depth files:

```bash
  OutDir=vs_589_cov_plot
for File in $(ls alignment/genome_alignment/vs_589/P.stipitis_59*/P.stipitis_59*_vs_589_depth_10kb.tsv); do
  Strain=$(echo $File | cut -f4 -d '/')
  echo $Strain
  cat $File | awk '{print $1,$2-1000,$2,$3,$4}' OFS='\t' | cut -f1,2,3,4 | sed 's/contig_/g' > $OutDir/${Strain}_vs_ref_hardmasked_scatterplot.tsv
done

```bash
ProgDir=/home/vegasa/git_repos/scripts/pichia/genome_alignment/coverage_plot/vs_589
circos -conf $ProgDir/vs_589_circos.conf -outputdir $OutDir
mv $OutDir/circos.png $OutDir/vs_589_cov_unmasked_circos.png
mv $OutDir/circos.svg $OutDir/vs_589_cov_unmasked_circos.svg
ls $PWD/$OutDir/vs_589*_circos.png
```
