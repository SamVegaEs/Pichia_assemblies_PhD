
OutDir=Chromosome_alignment_circos
mkdir -p $OutDir

ProgDir=~/git_repos/scripts/pichia/genome_alignment/Chromosome_alignment/vs_Parental/Parental_vs_591/
circos -conf $ProgDir/589_vs_591_circos.conf -outputdir $OutDir
mv $OutDir/circos.png $OutDir/Parental_vs_591_circos.png
mv $OutDir/circos.svg $OutDir/Parental_vs_591.svg
ls $PWD/$OutDir/589_vs_591_circos.png
