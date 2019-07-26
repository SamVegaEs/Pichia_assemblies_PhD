

OutDir=594_vs_Y-11454_v2_circos
mkdir -p $OutDir

QueryGenome=$(ls repeat_masked/*/*/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep '594')
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
$ProgDir/fasta2circos.py --genome $QueryGenome --contig_prefix "594_" > $OutDir/query_genome.txt

Ref_genome=$(ls assembly/misc_publications/P.stipitis/Y-11545_v2/pichia.allmasked)
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
$ProgDir/fasta2circos.py --genome $Ref_genome --contig_prefix "Y-11545_" > $OutDir/ref_genome.txt

cat $OutDir/query_genome.txt > $OutDir/query_ref_genome.txt
tac $OutDir/ref_genome.txt >> $OutDir/query_ref_genome.txt

OutDir=594_vs_Y-11454_v2_circos
Coords=$(ls analysis/genome_alignment/mummer/P.stipitis/594/594_vs_Y-11454_v2/594_vs_Y-11454_v2_coords.tsv)
ProgDir=/home/armita/git_repos/emr_repos/scripts/alternaria/pathogen/genome_alignment
$ProgDir/nucmer_coords2circos.py --inp_coords $Coords --queery_id 594 --ref_id Y-11545 > $OutDir/query_vs_ref_links.txt
cat $OutDir/query_vs_ref_links.txt > $OutDir/query_vs_ref_links_edited.txt



cat $OutDir/query_ref_genome.txt | cut -f3 -d ' ' | sed "s/$/; /g" | tr -d '\n' > $OutDir/query_contig_order.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/circos
$ProgDir/find_contig_orientation.py --links_file $OutDir/query_vs_ref_links_edited.txt > $OutDir/query_vs_ref_contig_orientation.txt

cat $OutDir/query_vs_ref_contig_orientation.txt | grep -A1 'Order of all seen contigs' | tail -n1 | sed "s/, /\n/g" > tmp.txt
cat $OutDir/query_ref_genome.txt | grep 'Y-11545' | grep -w -v -f tmp.txt | cut -f3 -d ' '| tr -d '\n' | sed 's/Y-11545/, Y-11545/g'
# cat $OutDir/Ag_As_genome.txt | grep 'Ag' | cut -f3 -d ' ' | tr -d '\n' | sed 's/Ag/, Ag/g' >> tmp.txt
echo "Order of unseen 594 contigs and remaining 594 contigs"
cat $OutDir/query_ref_genome.txt | grep -w -v -f tmp.txt | cut -f3 -d ' '| tr -d '\n' | sed 's/594/, 594/g' | sed 's/Y-11545/, Y-11545/g'


, 594_contig_1, 594_contig_2, 594_contig_3, 594_contig_4, 594_contig_5, 594_contig_6, 594_contig_7, 594_contig_8, Y-11545_chr_8.1, Y-11545_chr_7.1, Y-11545_chr_6.1, Y-11545_chr_5.1, Y-11545_chr_4.1, Y-11545_chr_3.1, Y-11545_chr_2.1, Y-11545_chr_1.2, Y-11545_chr_1.1


ProgDir=~/git_repos/scripts/pichia/genome_alignment/vs_ref/594
circos -conf $ProgDir/594_vs_Y-11454_v2_circos.conf -outputdir $OutDir
mv $OutDir/circos.png $OutDir/594_vs_Y-11545_v2_circos.png
mv $OutDir/circos.svg $OutDir/594_vs_Y-11545_v2_circos.svg
ls $PWD/$OutDir/594_vs_Y-11545_v2_circos.png
