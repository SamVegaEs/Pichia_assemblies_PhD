#SNP Calling: 589_vs_594

Parts of the commands had been already run by Andy for the SNP Calling of 591, so I took them for 594.

2.1. Pre SNP calling clean up 

```bash
Isolate="589"
Reference=$(ls analysis_aa/popgen/SNP_calling/589_contigs_unmasked.fa)
CurDir=$PWD
OutDir=analysis/popgen/SNP_calling/594_vs_589
mkdir -p $OutDir
cd $OutDir
ProgDir=~/git_repos/scripts/pichia/popgen
sbatch $ProgDir/slurm_SNP_calling_594_vs_589.sh $Reference $Isolate
cd $CurDir
```

#Filter SNPs based on this region being present in all isolates

Only retain biallelic high-quality SNPS with no missing data (for any individual) for genetic analyses below (in some cases, may allow some missing data in order to retain more SNPs, or first remove poorly sequenced individuals with too much missing data and then filter the SNPs)

```bash
srun --partition=long --pty bash
Vcf=$(ls analysis/popgen/SNP_calling/594_vs_589/vs_589/589_contigs_unmasked_temp.vcf)
vcftools=/projects/oldhome/sobczm/bin/vcftools/bin
vcflib=/projects/oldhome/sobczm/bin/vcflib/bin
mq=40
qual=30
dp=10
gq=30
na=1.00
removeindel=N
echo "count prefilter"
cat ${Vcf} | grep -v '#' | wc -l
$vcflib/vcffilter -f "QUAL > $qual & MQ > $mq" $Vcf \
| $vcflib/vcffilter -g "DP > $dp & GQ > $gq" > ${Vcf%.vcf}_qfiltered.vcf
echo "count qfilter"
cat ${Vcf%.vcf}_qfiltered.vcf | grep -v '#' | wc -l
$vcftools/vcftools --vcf ${Vcf%.vcf}_qfiltered.vcf --max-missing $na --remove-indels --recode --out ${Vcf%.vcf}_qfiltered_presence
```

```
After filtering, kept 2 out of 2 Individuals
Outputting VCF file...
After filtering, kept 30 out of a possible 44 Sites
Run Time = 0.00 seconds
```

#Collect VCF stats

General VCF stats (remember that vcftools needs to have the PERL library exported)

```bash
 Isolate="589"
  VcfTools=/projects/oldhome/sobczm/bin/vcftools/bin
  export PERL5LIB="$VcfTools:$PERL5LIB"
  VcfFiltered=$(ls analysis/popgen/SNP_calling/*/vs_${Isolate}/*_qfiltered_presence*.vcf)
  Stats=$(echo $VcfFiltered | sed 's/.vcf/.stat/g')
  perl $VcfTools/vcf-stats $VcfFiltered > $Stats
```
Calculate the index for percentage of shared SNP alleles between the individuals.

```bash
 Isolate="589"
  for Vcf in $(ls analysis/popgen/SNP_calling/*/vs_${Isolate}/*_qfiltered_presence*.vcf); do
      ProgDir=/projects/oldhome/armita/git_repos/emr_repos/scripts/popgen/snp
      $ProgDir/similarity_percentage.py $Vcf
  done
```

Identify SNPs in gene models:

Create custom SnpEff genome database

```
SnpEff=/projects/oldhome/sobczm/bin/snpEff
nano $SnpEff/snpEff.config
```

Add the following lines to the section with databases:

```
#---
# EMR Databases
#----
# Fus2 genome
Fus2v1.0.genome : Fus2
# Bc16 genome
Bc16v1.0.genome: BC-16
# P414 genome
P414v1.0.genome: 414
# 62471 genome
62471v1.0.genome: 62471
# R36_14 genome
R36_14v1.0.genome: R36_14
# SCRP371 genome
SCRP371v1.0.genome: SCRP371
# P. stipis
Ps589v1.0.genome: 589
```
Collect input files

```bash
Organism="P.stipitis"
Strain="589"
DbName="Ps589v1.0"
ProjDir=$PWD
Reference=$(ls $ProjDir/repeat_masked/${Organism}/${Strain}/filtered_contigs/${Strain}_contigs_unmasked.fa)
Gff=$(ls $ProjDir/gene_pred/fungap/${Organism}/${Strain}/fungap_out/fungap_out.gff3)
SnpEff=/projects/oldhome/sobczm/bin/snpEff
mkdir $SnpEff/data/${DbName}
cp $Reference $SnpEff/data/${DbName}/sequences.fa
cp $Gff $SnpEff/data/${DbName}/genes.gff

#Build database using GFF3 annotation
java -jar $SnpEff/snpEff.jar build -gff3 -v ${DbName}
```

Annotate VCF files

```bash
Organism="P.stipitis"
Isolate="589"
DbName="Ps589v1.0"
CurDir=/projects/oldhome/groups/harrisonlab/project_files/Pichia
cd $CurDir
  for Vcf in $(ls analysis/popgen/SNP_calling/*/vs_${Isolate}/*_qfiltered_presence*.recode.vcf); do
    echo $Vcf
    filename=$(basename "$Vcf")
    Prefix=${filename%.vcf}
    OutDir=$(dirname $Vcf)
    SnpEff=/projects/oldhome/sobczm/bin/snpEff
    java -Xmx4g -jar $SnpEff/snpEff.jar -v -ud 0 ${DbName} $Vcf > $OutDir/"$Prefix"_annotated.vcf
    mv snpEff_genes.txt $OutDir/snpEff_genes_"$Prefix".txt
    mv snpEff_summary.html $OutDir/snpEff_summary_"$Prefix".html
    # mv 414_v2_contigs_unmasked_filtered* $OutDir/.
    #-
    #Create subsamples of SNPs containing those in a given category
    #-
    #genic (includes 5', 3' UTRs)
    java -jar $SnpEff/SnpSift.jar filter "(ANN[*].EFFECT has 'missense_variant') || (ANN[*].EFFECT has 'nonsense_variant') || (ANN[*].EFFECT has 'synonymous_variant') || (ANN[*].EFFECT has 'intron_variant') || (ANN[*].EFFECT has '5_prime_UTR_variant') || (ANN[*].EFFECT has '3_prime_UTR_variant')" $OutDir/"$Prefix"_annotated.vcf > $OutDir/"$Prefix"_gene.vcf
    #coding
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant') || (ANN[0].EFFECT has 'synonymous_variant')" $OutDir/"$Prefix"_annotated.vcf > $OutDir/"$Prefix"_coding.vcf
    #non-synonymous
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant')" $OutDir/"$Prefix"_annotated.vcf > $OutDir/"$Prefix"_nonsyn.vcf
    #synonymous
    java -jar $SnpEff/SnpSift.jar filter "(ANN[0].EFFECT has 'synonymous_variant')" $OutDir/"$Prefix"_annotated.vcf > $OutDir/"$Prefix"_syn.vcf
    #Four-fold degenrate sites (output file suffix: 4fd)
    ProgDir=/projects/oldhome/sobczm/bin/popgen/summary_stats
    python $ProgDir/parse_snpeff_synonymous.py $OutDir/"$Prefix"_syn.vcf
    AllSnps=$(cat $OutDir/"$Prefix"_annotated.vcf | grep -v '#' | wc -l)
    GeneSnps=$(cat $OutDir/"$Prefix"_gene.vcf | grep -v '#' | wc -l)
    CdsSnps=$(cat $OutDir/"$Prefix"_coding.vcf | grep -v '#' | wc -l)
    NonsynSnps=$(cat $OutDir/"$Prefix"_nonsyn.vcf | grep -v '#' | wc -l)
    SynSnps=$(cat $OutDir/"$Prefix"_syn.vcf | grep -v '#' | wc -l)
    printf "$AllSnps\t$GeneSnps\t$CdsSnps\t$NonsynSnps\t$SynSnps\n"
done
```
