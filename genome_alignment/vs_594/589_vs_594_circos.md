#Create a folder for the genome alignment 594 vs 589

```bash
OutDir=589_vs_594_circos 
mkdir -p $OutDir
```

#Transform our .fa files with the genomes in .txt files (both the sample and the reference genomes.)

```bash
QueryGenome=$(ls repeat_masked/*/*/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep '594')
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
$ProgDir/fasta2circos.py --genome $QueryGenome --contig_prefix "594_" > $OutDir/query_genome.txt

Ref_genome=$(ls repeat_masked/*/*/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep '589')
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
$ProgDir/fasta2circos.py --genome $Ref_genome --contig_prefix "589_" > $OutDir/ref_genome.txt

cat $OutDir/query_genome.txt > $OutDir/query_ref_genome.txt
tac $OutDir/ref_genome.txt >> $OutDir/query_ref_genome.txt
```



```bash
OutDir=589_vs_594_circos
Coords=$(ls analysis/genome_alignment/mummer/P.stipitis/589/589_vs_594/589_vs_594_coords.tsv)
ProgDir=/home/armita/git_repos/emr_repos/scripts/alternaria/pathogen/genome_alignment
$ProgDir/nucmer_coords2circos.py --inp_coords $Coords --queery_id 589 --ref_id 594 > $OutDir/query_vs_ref_links.txt
cat $OutDir/query_vs_ref_links.txt > $OutDir/query_vs_ref_links_edited.txt
```

#A file showing contig orientations was made:

```bash
cat $OutDir/query_ref_genome.txt | cut -f3 -d ' ' | sed "s/$/; /g" | tr -d '\n' > $OutDir/query_contig_order.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/circos
$ProgDir/find_contig_orientation.py --links_file $OutDir/query_vs_ref_links_edited.txt > $OutDir/query_vs_ref_contig_orientation.txt
```

#Contig order was selected by taking the first line of that file and then also taking the reversed order of contigs using the command:

```bash
cat $OutDir/query_vs_ref_contig_orientation.txt | grep -A1 'Order of all seen contigs' | tail -n1 | sed "s/, /\n/g" > tmp.txt
cat $OutDir/query_ref_genome.txt | grep '589' | grep -w -v -f tmp.txt | cut -f3 -d ' '| tr -d '\n' | sed 's/589/, 589/g'
# cat $OutDir/Ag_As_genome.txt | grep 'Ag' | cut -f3 -d ' ' | tr -d '\n' | sed 's/Ag/, Ag/g' >> tmp.txt

echo "Order of unseen 594 contigs and remaining 594 contigs"
cat $OutDir/query_ref_genome.txt | grep -w -v -f tmp.txt | cut -f3 -d ' '| tr -d '\n' | sed 's/594/, 594/g' | sed 's/589/, 589/g'
```



# , 594_contig_1, 594_contig_2, 594_contig_3, 594_contig_4, 594_contig_5, 594_contig_6, 594_contig_7, 594_contig_8, 589_contig_11

```bash
ProgDir=~/git_repos/emr_repos/scripts/pichia/genome_alignment/vs_589
circos -conf $ProgDir/589_vs_594_circos.conf -outputdir $OutDir
mv $OutDir/circos.png $OutDir/589_vs_594_circos.png
mv $OutDir/circos.svg $OutDir/589_vs_594_circos.svg
ls $PWD/$OutDir/589_vs_594_circos.png
```
