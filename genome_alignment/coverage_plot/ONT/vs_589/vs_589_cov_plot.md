```bash
 OutDir=cov_chromosome_level/ONT
  mkdir -p $OutDir

  Ref_genome=$(ls chromosome_files/*/589/589_all.fa)
  ProgDir=~/git_repos/tools/seq_tools/circos
  python $ProgDir/fasta2circos.py --genome $Ref_genome --contig_prefix "589_" > $OutDir/query_genome.txt
```


Calculate coverage over 10kb windows for unmasked files. Parse read depth files:

```bash
 OutDir=cov_chromosome_level/ONT
for File in $(ls analysis/genome_alignment/minimap/chromosomes/*/*/vs_*/*_depth_10kb.tsv); do
  Strain=$(echo $File | cut -f7 -d '/')
  echo $Strain 
  cat $File | awk '{print $1,$2-1000,$2,$3,$4}' OFS='\t' | cut -f1,2,3,4 | sed 's/contig_/589_contig_/g' > $OutDir/${Strain}_vs_ref_hardmasked_scatterplot.tsv
done
```
```bash
OutDir=/projects/oldhome/groups/harrisonlab/project_files/Pichia/cov_chromosome_level/ONT
ProgDir=/projects/oldhome/vegasa/git_repos/scripts/pichia/genome_alignment/Chromosome_level_assembly/cov_plot/ONT/vs_589
circos -conf $ProgDir/vs_589_circos.conf -outputdir $OutDir
mv $OutDir/circos.png $OutDir/vs_589_cov_unmasked_circos.png
mv $OutDir/circos.svg $OutDir/vs_589_cov_unmasked_circos.svg
ls $PWD/$OutDir/vs_589*_circos.png
```
