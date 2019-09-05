```bash
OutDir=Chromosome_alignment_circos
```
```bash
ProgDir=~/git_repos/scripts/pichia/genome_alignment/Chromosome_alignment/vs_Parental/Parental_vs_594/
circos -conf $ProgDir/Parental_vs_594_circos.conf -outputdir $OutDir
mv $OutDir/circos.png $OutDir/Parental_vs_594_circos.png
mv $OutDir/circos.svg $OutDir/Parental_vs_594.svg
```
