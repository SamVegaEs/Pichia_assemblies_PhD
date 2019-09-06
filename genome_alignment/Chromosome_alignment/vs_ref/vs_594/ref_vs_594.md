
```bash
OutDir=Chromosome_alignment_circos
```

```bash
ProgDir=~/git_repos/scripts/pichia/genome_alignment/Chromosome_alignment/vs_ref/vs_594/
circos -conf $ProgDir/ref_vs_594_circos.conf -outputdir $OutDir
mv $OutDir/circos.png $OutDir/ref_vs_594_circos.png
mv $OutDir/circos.svg $OutDir/ref_vs_594_circos.svg
```
