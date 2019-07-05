```bash
  OutDir=genome_alignment/circos/589_vs_Y-11454_v2_circos
  mkdir -p $OutDir

  QueryGenome=$(ls repeat_masked/*/*/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep '589')
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/fasta2circos.py --genome $QueryGenome --contig_prefix "589_" > $OutDir/query_genome.txt

  Ref_genome=$(ls assembly/misc_publications/P.stipitis/Y-11545_v2/pichia.allmasked)
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/fasta2circos.py --genome $Ref_genome --contig_prefix "Y-11545_" > $OutDir/ref_genome.txt

  cat $OutDir/query_genome.txt > $OutDir/query_ref_genome.txt
  tac $OutDir/ref_genome.txt >> $OutDir/query_ref_genome.txt
```
<!--
Telomere locations on contigs:

```bash
cat analysis/telomere/A.gaisen/650/telomere_hits_circos.txt | sed 's/contig/Ag_contig/g' | sort -k3 -n -t'_' > $OutDir/Ag_vs_As_telomere_hits.txt
cat /home/groups/harrisonlab/project_files/alternaria/analysis/telomere/A.solani/altNL03003/telomere_hits_circos.txt  | sed 's/CP/As_CP/g' | sort -k3 -n -t'_' >> $OutDir/Ag_vs_As_telomere_hits.txt
```
-->

```bash
OutDir=genome_alignment/circos/589_vs_Y-11454_v2_circos
Coords=$(ls analysis/genome_alignment/mummer/P.stipitis/589/589_vs_Y-11454_v2/589_vs_Y-11454_v2_coords.tsv)
ProgDir=/home/armita/git_repos/emr_repos/scripts/alternaria/pathogen/genome_alignment
$ProgDir/nucmer_coords2circos.py --inp_coords $Coords --queery_id 589 --ref_id Y-11545 > $OutDir/query_vs_ref_links.txt
cat $OutDir/query_vs_ref_links.txt > $OutDir/query_vs_ref_links_edited.txt
```


A file showing contig orientations was made:
```bash
  cat $OutDir/query_ref_genome.txt | cut -f3 -d ' ' | sed "s/$/; /g" | tr -d '\n' > $OutDir/query_contig_order.txt
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/circos
  $ProgDir/find_contig_orientation.py --links_file $OutDir/query_vs_ref_links_edited.txt > $OutDir/query_vs_ref_contig_orientation.txt
```


Contig order was selected by taking the first line of that file and then also
taking the reversed order of contigs using the command:

```bash
cat $OutDir/query_vs_ref_contig_orientation.txt | grep -A1 'Order of all seen contigs' | tail -n1 | sed "s/, /\n/g" > tmp.txt
cat $OutDir/query_ref_genome.txt | grep 'Y-11545' | grep -w -v -f tmp.txt | cut -f3 -d ' '| tr -d '\n' | sed 's/Y-11545/, Y-11545/g'
# cat $OutDir/Ag_As_genome.txt | grep 'Ag' | cut -f3 -d ' ' | tr -d '\n' | sed 's/Ag/, Ag/g' >> tmp.txt

echo "Order of unseen 589 contigs and remaining 589 contigs"
cat $OutDir/query_ref_genome.txt | grep -w -v -f tmp.txt | cut -f3 -d ' '| tr -d '\n' | sed 's/589/, 589/g' | sed 's/Y-11545/, Y-11545/g'
```

```
, 589_contig_11, Y-11545_chr_8.1, Y-11545_chr_7.1, Y-11545_chr_6.1, Y-11545_chr_5.1, Y-11545_chr_4.1, Y-11545_chr_3.1, Y-11545_chr_2.1, Y-11545_chr_1.2, Y-11545_chr_1.1
```

```bash
ProgDir=/home/armita/git_repos/emr_repos/scripts/pichia/genome_alignment/vs_ref/589
circos -conf $ProgDir/589_vs_Y-11454_v2_circos.conf -outputdir $OutDir
mv $OutDir/circos.png $OutDir/589_vs_Y-11545_v2_circos.png
mv $OutDir/circos.svg $OutDir/589_vs_Y-11545_v2_circos.svg
ls $PWD/$OutDir/589_vs_Y-11545_v2_circos.png
```
<!--
```bash
OutDir=analysis/circos/Ag_vs_As_circos
cat $OutDir/At_Ag_genome_edited2.txt | grep -v '4287' > $OutDir/At_Ag_genome_final.txt
mkdir -p $OutDir/by_FoC_chr
for Num in $(seq 1 22); do
  Chr="contig_"$Num"_pilon"
  echo "$Chr"
  OrthologyTxt=analysis/orthology/orthomcl/FoC_vs_Fo_vs_FoL_publication_ncbi/FoC_vs_Fo_vs_FoL_publication_ncbi_orthogroups.txt
  ProgDir=~/git_repos/emr_repossh s/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/orthology2ribbons_internal.py \
  --chr1 $Chr \
  --orthology $OrthologyTxt \
  --name1 Fus2 \
  --gff1 gene_pred/final_genes/F.oxysporum_fsp_cepae/Fus2_canu_new/final/final_genes_appended.gff3 \
  | sort | uniq \
  > $OutDir/Ag_vs_As_links_edited.txt
  ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  circos -conf $ProgDir/Fus2/Fus2_FoL/Ag_vs_As_circos.conf -outputdir $OutDir
  mv $OutDir/circos.png $OutDir/by_FoC_chr/Ag_vs_As_LS_"$Chr"_circos.png
  mv $OutDir/circos.svg $OutDir/by_FoC_chr/Ag_vs_As_LS_"$Chr"_circos.svg
done
``` -->
<!--
The frequency of gene duplications within and between chromosomes was investigated:

```bash
OutDir=analysis/circos/Ag_vs_As_circos
for Num in $(seq 1 22); do
Chr="contig_"$Num"_pilon"
ChrList="$ChrList $Chr"
done
echo "$ChrList"
OrthologyTxt=analysis/orthology/orthomcl/FoC_vs_Fo_vs_FoL_publication_ncbi/FoC_vs_Fo_vs_FoL_publication_ncbi_orthogroups.txt
ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
$ProgDir/orthology2ribbons_internal.py \
--chr1 $ChrList \
--orthology $OrthologyTxt \
--name1 Fus2 \
--gff1 gene_pred/final_genes/F.oxysporum_fsp_cepae/Fus2_canu_new/final/final_genes_appended.gff3 \
| sort | uniq > $OutDir/Fus2_all_Fus2_links.txt
cat $OutDir/Fus2_all_Fus2_links.txt | cut -f1,4 | sort | uniq -c > $OutDir/Fus2_all_Fus2_links_occurence.txt

``` -->
