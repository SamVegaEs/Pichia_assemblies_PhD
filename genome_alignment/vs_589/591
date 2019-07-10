#Create a folder for the genome alignment 591 vs 589

```bash
OutDir=genome_alignment/circos/589_vs_591_circos
mkdir -p $OutDir
```

#Transform our .fa files with the genomes in .txt files (both the sample and the reference genomes.)

```bash
QueryGenome=$(ls repeat_masked/*/*/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep '591')
ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
$ProgDir/fasta2circos.py --genome $QueryGenome --contig_prefix "591_" > $OutDir/query_genome.txt

Ref_genome=$(ls repeat_masked/*/*/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep '589')
ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
$ProgDir/fasta2circos.py --genome $Ref_genome --contig_prefix "589_" > $OutDir/ref_genome.txt

cat $OutDir/query_genome.txt > $OutDir/query_ref_genome.txt
tac $OutDir/ref_genome.txt >> $OutDir/query_ref_genome.txt
```




OutDir=genome_alignment/circos/589_vs_591_circos
Coords=$(ls analysis/genome_alignment/mummer/P.stipitis/589/589_vs_591/589_vs_591_coords.tsv)
ProgDir=/home/armita/git_repos/emr_repos/scripts/alternaria/pathogen/genome_alignment
$ProgDir/nucmer_coords2circos.py --inp_coords $Coords --queery_id 591 --ref_id 589 > $OutDir/query_vs_ref_links.txt
cat $OutDir/query_vs_ref_links.txt > $OutDir/query_vs_ref_links_edited.txt


#A file showing contig orientations was made:

cat $OutDir/query_ref_genome.txt | cut -f3 -d ' ' | sed "s/$/; /g" | tr -d '\n' > $OutDir/query_contig_order.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/circos
$ProgDir/find_contig_orientation.py --links_file $OutDir/query_vs_ref_links_edited.txt > $OutDir/query_vs_ref_contig_orientation.txt


#Contig order was selected by taking the first line of that file and then also taking the reversed order of contigs using the command:

cat $OutDir/query_vs_ref_contig_orientation.txt | grep -A1 'Order of all seen contigs' | tail -n1 | sed "s/, /\n/g" > tmp.txt
cat $OutDir/query_ref_genome.txt | grep '589' | grep -w -v -f tmp.txt | cut -f3 -d ' '| tr -d '\n' | sed 's/589/, 589/g'
# cat $OutDir/Ag_As_genome.txt | grep 'Ag' | cut -f3 -d ' ' | tr -d '\n' | sed 's/Ag/, Ag/g' >> tmp.txt

echo "Order of unseen 591 contigs and remaining 591 contigs"
cat $OutDir/query_ref_genome.txt | grep -w -v -f tmp.txt | cut -f3 -d ' '| tr -d '\n' | sed 's/591/, 591/g' | sed 's/589/, 589/g'



ProgDir=/home/armita/git_repos/emr_repos/scripts/pichia/genome_alignment/vs_589/591
circos -conf $ProgDir/589_vs_Y-11454_v2_circos.conf -outputdir $OutDir
mv $OutDir/circos.png $OutDir/589_vs_591_circos.png
mv $OutDir/circos.svg $OutDir/589_vs_591_circos.svg
ls $PWD/$OutDir/589_vs_591_circos.png
