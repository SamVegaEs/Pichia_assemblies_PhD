
# SNP calling codes.

The alignment by BWA has ben conducted already during the genome assembly. So we will start from the Pre SNP calling clean up step.

# 1. vs published Ref

# 0.1 Aligning

All against reference genome.

```bash
cd /projects/oldhome/groups/harrisonlab/project_files/Pichia
for Reference in $(ls assembly/misc_publications/P.stipitis/CBS6054/GCF_000209165.1_ASM20916v1_genomic.fna); do
  Ref=$(echo $Reference | rev | cut -f2 -d '/' | rev)
for StrainPath in $(ls -d qc_dna/paired/P.stipitis/*); do
Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
F_Read=$(ls $StrainPath/F/*_trim.fq.gz)
R_Read=$(ls $StrainPath/R/*_trim.fq.gz)
echo $F_Read
echo $R_Read
Prefix="${Organism}_${Strain}"
OutDir=analysis_aa/genome_alignment/bwa/$Organism/$Strain/vs_${Ref}
ProgDir=/projects/oldhome/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/bwa
sbatch $ProgDir/slurm_bwa.sh $Prefix $Reference $F_Read $R_Read $OutDir
done
done
```

## 1.1. Pre SNP calling clean up

### 1.1.1  

Rename input mapping files in each folder by prefixing with the strain ID.

Alignment of the evolved strains were made against the parental strain.


```bash
  for File in $(ls analysis_aa/genome_alignment/bwa/P.stipitis/*/vs_CBS6054/*sorted.bam); do
    Strain=$(echo $File | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $File | rev | cut -f4 -d '/' | rev)
    Prefix="vs_CBS6054"
    echo "$Organism - $Strain"
    OutDir=analysis_aa/popgen/${Prefix}/$Organism/$Strain
    CurDir=$PWD
    echo $OutDir
    mkdir -p $OutDir
    cp $CurDir/$File $OutDir/${Strain}_${Prefix}_sorted.bam
  done
```

### 1.1.2
Remove multimapping reads, discordant reads. PCR and optical duplicates, and add read group and sample name to each mapped read (preferably, the shortest ID possible)

Convention used: qsub $ProgDir/sub_pre_snp_calling.sh <SAMPLE_ID>

```bash
conda activate gatk4
Reference=$(ls assembly/misc_publications/P.stipitis/CBS6054/GCF_000209165.1_ASM20916v1_genomic.fna)
for Sam in $(ls analysis_aa/popgen/vs_CBS6054/*/*/*CBS6054_sorted.bam); do
Strain=$(echo $Sam | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Sam | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/projects/oldhome/armita/git_repos/emr_repos/scripts/pichia/popgen
sbatch $ProgDir/slurm_pre_snp_calling.sh $Sam $Strain $Reference
done
```


Prepare genome reference indexes required by GATK. Prepare for 589, 591 and 594.

```bash
for Reference in $(ls assembly/misc_publications/P.stipitis/CBS6054/GCF_000209165.1_ASM20916v1_genomic.fna); do
OutName=$(echo $Reference | sed 's/.fna/.dict/g')
OutDir=$(dirname $Reference)
mkdir -p $OutDir
ProgDir=/projects/oldhome/sobczm/bin/picard-tools-2.5.0
java -jar $ProgDir/picard.jar CreateSequenceDictionary R=$Reference O=$OutName
samtools faidx $Reference
done
```

```bash
Isolate="CBS6054"
Reference=$(ls assembly/misc_publications/P.stipitis/CBS6054/GCF_000209165.1_ASM20916v1_genomic.fna)
CurDir=$PWD
OutDir=analysis_aa/popgen/SNP_calling
mkdir -p $OutDir
cd $OutDir
conda activate gatk4
ProgDir=/projects/oldhome/armita/git_repos/emr_repos/scripts/pichia/popgen
sbatch $ProgDir/slurm_snp_calling_ncbi.sh $Reference $Isolate
cd $CurDir
```


Visualising SNPs

https://knausb.github.io/vcfR_documentation/visualization_1.html
<!--
```
install.packages('vcfR')
library(vcfR)
vcf <- read.vcfR('/Users/armita/Downloads/589_contigs_softmasked_repeatmasker_TPSI_appended_temp.vcf', verbose = FALSE)

# Create a chromR object.
create.chromR(vcf, name = "CHROM", seq = NULL, ann = NULL, verbose = TRUE)
chromoqc(chrom, dp.alpha = 22)
plot(chrom)

```

```bash
Vcf=$(ls /Users/armita/Downloads/589_contigs_softmasked_repeatmasker_TPSI_appended_temp.vcf)
# cat $Vcf | grep -v '##' | cut -f6,9,10,11,12 | less
printf "Isolate\tCount\tDP\n" > /Users/armita/Downloads/vs_589_dp.tsv
cat $Vcf | grep -v '##' | cut -f10 | cut -f3 -d ':' | grep -v "^\.$" | sort -n | nl -d "\t" | sed "s/ //g" | sed $'s/^/589\t/g' >> /Users/armita/Downloads/vs_589_dp.tsv
cat $Vcf | grep -v '##' | cut -f11 | cut -f3 -d ':' | grep -v "^\.$" | sort -n | nl -d "\t" | sed "s/ //g" | sed $'s/^/591\t/g' >> /Users/armita/Downloads/vs_589_dp.tsv
cat $Vcf | grep -v '##' | cut -f12 | cut -f3 -d ':' | grep -v "^\.$" | sort -n | nl -d "\t" | sed "s/ //g" | sed $'s/^/594\t/g' >> /Users/armita/Downloads/vs_589_dp.tsv
printf "Isolate\tCount\tGQ\n" > /Users/armita/Downloads/vs_589_gq.tsv
cat $Vcf | grep -v '##' | cut -f10 | cut -f4 -d ':' | grep -v "^\.$" | sort -n | nl -d "\t" | sed "s/ //g" | sed $'s/^/589\t/g' >> /Users/armita/Downloads/vs_589_gq.tsv
cat $Vcf | grep -v '##' | cut -f11 | cut -f4 -d ':' | grep -v "^\.$" | sort -n | nl -d "\t" | sed "s/ //g" | sed $'s/^/591\t/g' >> /Users/armita/Downloads/vs_589_gq.tsv
cat $Vcf | grep -v '##' | cut -f12 | cut -f4 -d ':' | grep -v "^\.$" | sort -n | nl -d "\t" | sed "s/ //g" | sed $'s/^/594\t/g' >> /Users/armita/Downloads/vs_589_gq.tsv
```

```
library(ggplot2)
vs_589_dp <- read.delim("~/Downloads/vs_589_dp.tsv")
vs_589_dp$Isolate <- as.factor(vs_589_dp$Isolate)
ggplot(vs_589_dp, aes(x=Count, y=DP, color=Isolate)) + geom_line()
p1 <- ggplot(vs_589_dp, aes(x=Count, y=DP, color=Isolate)) + geom_line() + ylim(0,200)
ggsave("~/Downloads/vs_589_dp.pdf", p1)


vs_589_gq <- read.delim("~/Downloads/vs_589_gq.tsv")
vs_589_gq$Isolate <- as.factor(vs_589_gq$Isolate)
ggplot(vs_589_gq, aes(x=Count, y=GQ, color=Isolate)) + geom_line()
p2 <- ggplot(vs_589_gq, aes(x=Count, y=GQ, color=Isolate)) + geom_line() + ylim(0,100)
ggsave("~/Downloads/vs_589_gq.pdf", p2)
```
-->



## Filter SNPs based on this region being present in all isolates

Only retain biallelic high-quality SNPS with no missing data (for any individual) for genetic analyses below (in some cases, may allow some missing data in order to retain more SNPs, or first remove poorly sequenced individuals with too much missing data and then filter the SNPs).

```bash
srun --partition=long --pty bash
Vcf=$(ls analysis_aa/popgen/SNP_calling/vs_Y-11545/pichia.AssembledScaffolds_temp.vcf)
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
After filtering, kept 3 out of 3 Individuals
Outputting VCF file...
After filtering, kept 4002 out of a possible 4565 Sites
Run Time = 0.00 seconds
```

## Gene models Y-11545
```bash
cd assembly/misc_publications/P.stipitis/Y-11545_v2
wget https://genome.jgi.doe.gov/portal/Picst3/download/Pstipitis_FilteredModelsv2.0.gtf.gz
wget https://genome.jgi.doe.gov/portal/Picst3/download/Pstipitisv2.FilteredModels1.transcripts.gz
wget https://genome.jgi.doe.gov/portal/Picst3/download/Pstipitisv2.FilteredModels1.proteins.gz
wget https://genome.jgi.doe.gov/portal/Picst3/download/Pstipitisv2.ecpathwayinfo_FilteredModels1.tab.gz
wget https://genome.jgi.doe.gov/portal/Picst3/download/Pstipitisv2.goinfo_FilteredModels1.tab.gz
wget https://genome.jgi.doe.gov/portal/Picst3/download/Pstipitisv2.koginfo_FilteredModels1.tab.gz
wget https://genome.jgi.doe.gov/portal/Picst3/download/Pstipitisv2.domaininfo_FilteredModels1.tab.gz
wget https://genome.jgi.doe.gov/portal/Picst3/download/Pstipitisv2.allModels.proteins.gz
wget https://genome.jgi.doe.gov/portal/Picst3/download/Pstipitisv2.allModels.transcripts.gz
gunzip *.gz
cd /projects/oldhome/groups/harrisonlab/project_files/Pichia

mkdir -p assembly/misc_publications/P.stipitis/CBS6054
cd assembly/misc_publications/P.stipitis/CBS6054
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/209/165/GCF_000209165.1_ASM20916v1/GCF_000209165.1_ASM20916v1_cds_from_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/209/165/GCF_000209165.1_ASM20916v1/GCF_000209165.1_ASM20916v1_feature_table.txt.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/209/165/GCF_000209165.1_ASM20916v1/GCF_000209165.1_ASM20916v1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/209/165/GCF_000209165.1_ASM20916v1/GCF_000209165.1_ASM20916v1_genomic.gbff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/209/165/GCF_000209165.1_ASM20916v1/GCF_000209165.1_ASM20916v1_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/209/165/GCF_000209165.1_ASM20916v1/GCF_000209165.1_ASM20916v1_genomic.gtf.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/209/165/GCF_000209165.1_ASM20916v1/GCF_000209165.1_ASM20916v1_protein.faa.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/209/165/GCF_000209165.1_ASM20916v1/GCF_000209165.1_ASM20916v1_protein.gpff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/209/165/GCF_000209165.1_ASM20916v1/GCF_000209165.1_ASM20916v1_rm.out.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/209/165/GCF_000209165.1_ASM20916v1/GCF_000209165.1_ASM20916v1_rm.run
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/209/165/GCF_000209165.1_ASM20916v1/GCF_000209165.1_ASM20916v1_rna.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/209/165/GCF_000209165.1_ASM20916v1/GCF_000209165.1_ASM20916v1_rna.gbff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/209/165/GCF_000209165.1_ASM20916v1/GCF_000209165.1_ASM20916v1_rna_from_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/209/165/GCF_000209165.1_ASM20916v1/GCF_000209165.1_ASM20916v1_translated_cds.faa.gz
gunzip *.gz
```


# Collect VCF stats

General VCF stats (remember that vcftools needs to have the PERL library exported)

```bash
  Isolate="Y-11545"
  VcfTools=/projects/oldhome/sobczm/bin/vcftools/bin
  export PERL5LIB="$VcfTools:$PERL5LIB"
  VcfFiltered=$(ls analysis_aa/popgen/SNP_calling/vs_${Isolate}/*_qfiltered_presence*.vcf)
  Stats=$(echo $VcfFiltered | sed 's/.vcf/.stat/g')
  perl $VcfTools/vcf-stats $VcfFiltered > $Stats
```
Calculate the index for percentage of shared SNP alleles between the individuals.

```bash
  Isolate="Y-11545"
  for Vcf in $(ls analysis_aa/popgen/SNP_calling/vs_${Isolate}/*_qfiltered_presence*.vcf); do
      ProgDir=/projects/oldhome/armita/git_repos/emr_repos/scripts/popgen/snp
      $ProgDir/similarity_percentage.py $Vcf
  done
```

Visualise the output as heatmap and clustering dendrogram
```bash
Isolate="Y-11545"
for Log in $(ls analysis_aa/popgen/SNP_calling/vs_${Isolate}/*distance.log); do
  ProgDir=/projects/oldhome/armita/git_repos/emr_repos/scripts/popgen/snp
  Rscript --vanilla $ProgDir/distance_matrix.R $Log
  mv Rplots.pdf analysis_aa/popgen/SNP_calling/vs_${Isolate}/.
done
```

Identify SNPs in gene models:

Create custom SnpEff genome database

```bash
SnpEff=/projects/oldhome/sobczm/bin/snpEff
nano $SnpEff/snpEff.config
```
Add the following lines to the section with databases:

```bash
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
PsY-11545v1.0.genome: Y-11545
```

Collect input files

```bash
Organism="P.stipitis"
Strain="Y-11545"
DbName="PsY-11545v1.0"
ProjDir=$PWD
Reference=$(ls $ProjDir/assembly/misc_publications/P.stipitis/Y-11545_v2/pichia.AssembledScaffolds.fasta)
Gtf=$(ls $ProjDir/assembly/misc_publications/P.stipitis/Y-11545_v2/Pstipitisv2.FrozenGeneCatalog20080115.gff)
SnpEff=/projects/oldhome/sobczm/bin/snpEff
mkdir $SnpEff/data/${DbName}
cp $Reference $SnpEff/data/${DbName}/sequences.fa
cp $Gtf $SnpEff/data/${DbName}/genes.gtf

#Build database using GFF3 annotation
java -jar $SnpEff/snpEff.jar build -gtf22 -v ${DbName}

```

Annotate VCF files
```bash
Organism="P.stipitis"
Strain="589"
DbName="Ps589v1.0"
CurDir=/projects/oldhome/groups/harrisonlab/project_files/Pichia
cd $CurDir
  for Vcf in $(ls analysis_aa/popgen/SNP_calling/vs_${Isolate}/*_qfiltered_presence*.vcf); do
    echo $Vcf
    filename=$(basename "$Vcf")
    Prefix=${filename%.vcf}
    OutDir=$(ls -d analysis_aa/popgen/SNP_calling)
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
All Gene CDS nonsyn syn
27	9	9	7	2


# 2. Vs 589

# 2.0.1 Aligning

All against reference genome.

```bash
cd /projects/oldhome/groups/harrisonlab/project_files/Pichia
for Reference in $(ls repeat_masked/*/589/filtered_contigs/*_contigs_unmasked.fa); do
Ref=$(echo $Reference | rev | cut -f3 -d '/' | rev)
for StrainPath in $(ls -d qc_dna/paired/P.stipitis/*); do
Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
F_Read=$(ls $StrainPath/F/*_trim.fq.gz)
R_Read=$(ls $StrainPath/R/*_trim.fq.gz)
echo $F_Read
echo $R_Read
Prefix="${Organism}_${Strain}"
OutDir=analysis_aa/genome_alignment/bwa/$Organism/$Strain/vs_${Ref}
ProgDir=/projects/oldhome/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/bwa
sbatch $ProgDir/slurm_bwa.sh $Prefix $Reference $F_Read $R_Read $OutDir
done
done
```

## 2.1. Pre SNP calling clean up

### 2.1.1  Rename input mapping files in each folder by prefixing with the strain ID.

Alignment of the evolved strains were made against the parental strain.

```bash
  for File in $(ls analysis/genome_alignment/bwa/P.stipitis/*/vs_589_unmasked/*sorted.bam); do
    Strain=$(echo $File | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $File | rev | cut -f4 -d '/' | rev)
    Prefix=$(echo $File | rev | cut -f2 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=analysis_aa/popgen/$Organism/$Strain
    CurDir=$PWD
    echo $OutDir
    mkdir -p $OutDir
    cp $CurDir/$File $OutDir/${Strain}_${Prefix}_sorted.bam
  done
```

## 1.2 Remove multimapping reads, discordant reads. PCR and optical duplicates, and add read group and sample name to each mapped read (preferably, the shortest ID possible)

Convention used: qsub $ProgDir/sub_pre_snp_calling.sh <SAMPLE_ID>

```bash
conda activate gatk4
Reference=$(ls repeat_masked/*/589/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
for Sam in $(ls analysis_aa/popgen/*/*/*unmasked_sorted.bam); do
Strain=$(echo $Sam | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Sam | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/projects/oldhome/armita/git_repos/emr_repos/scripts/pichia/popgen
# sbatch $ProgDir/slurm_pre_snp_calling2.sh $Sam $Strain
sbatch $ProgDir/slurm_pre_snp_calling.sh $Sam $Strain $Reference
done
```


# 2. Run SNP calling

Prepare genome reference indexes required by GATK. Prepare for 589, 591 and 594.

```bash
for Reference in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
OutName=$(echo $Reference | sed 's/.fa/.dict/g')
OutDir=$(dirname $Reference)
mkdir -p $OutDir
ProgDir=/projects/oldhome/sobczm/bin/picard-tools-2.5.0
java -jar $ProgDir/picard.jar CreateSequenceDictionary R=$Reference O=$OutName
samtools faidx $Reference
done
```
Copy index file to same folder as BAM alignments

Move to the directory where the output of SNP calling should be placed. Then Start SNP calling with GATK. The submission script required need to be custom-prepared for each analysis, depending on what samples are being analysed. See inside the submission script below (GATK codes):

In order to run GATK I have used the next two set of commands, both of them start the run, but the outputs are empty. Folders are created normally but there is nothing inside.

```bash
Isolate="589"
Reference=$(ls repeat_masked/*/589/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
CurDir=$PWD
OutDir=analysis_aa/popgen/SNP_calling
mkdir -p $OutDir
cd $OutDir
conda activate gatk4
ProgDir=/projects/oldhome/armita/git_repos/emr_repos/scripts/pichia/popgen
sbatch $ProgDir/slurm_snp_calling.sh $Reference $Isolate
cd $CurDir
```


Visualising SNPs

https://knausb.github.io/vcfR_documentation/visualization_1.html

```
install.packages('vcfR')
library(vcfR)
vcf <- read.vcfR('/Users/armita/Downloads/589_contigs_softmasked_repeatmasker_TPSI_appended_temp.vcf', verbose = FALSE)

# Create a chromR object.
create.chromR(vcf, name = "CHROM", seq = NULL, ann = NULL, verbose = TRUE)
chromoqc(chrom, dp.alpha = 22)
plot(chrom)

```

```bash
Vcf=$(ls /Users/armita/Downloads/589_contigs_softmasked_repeatmasker_TPSI_appended_temp.vcf)
# cat $Vcf | grep -v '##' | cut -f6,9,10,11,12 | less
printf "Isolate\tCount\tDP\n" > /Users/armita/Downloads/vs_589_dp.tsv
cat $Vcf | grep -v '##' | cut -f10 | cut -f3 -d ':' | grep -v "^\.$" | sort -n | nl -d "\t" | sed "s/ //g" | sed $'s/^/589\t/g' >> /Users/armita/Downloads/vs_589_dp.tsv
cat $Vcf | grep -v '##' | cut -f11 | cut -f3 -d ':' | grep -v "^\.$" | sort -n | nl -d "\t" | sed "s/ //g" | sed $'s/^/591\t/g' >> /Users/armita/Downloads/vs_589_dp.tsv
cat $Vcf | grep -v '##' | cut -f12 | cut -f3 -d ':' | grep -v "^\.$" | sort -n | nl -d "\t" | sed "s/ //g" | sed $'s/^/594\t/g' >> /Users/armita/Downloads/vs_589_dp.tsv
printf "Isolate\tCount\tGQ\n" > /Users/armita/Downloads/vs_589_gq.tsv
cat $Vcf | grep -v '##' | cut -f10 | cut -f4 -d ':' | grep -v "^\.$" | sort -n | nl -d "\t" | sed "s/ //g" | sed $'s/^/589\t/g' >> /Users/armita/Downloads/vs_589_gq.tsv
cat $Vcf | grep -v '##' | cut -f11 | cut -f4 -d ':' | grep -v "^\.$" | sort -n | nl -d "\t" | sed "s/ //g" | sed $'s/^/591\t/g' >> /Users/armita/Downloads/vs_589_gq.tsv
cat $Vcf | grep -v '##' | cut -f12 | cut -f4 -d ':' | grep -v "^\.$" | sort -n | nl -d "\t" | sed "s/ //g" | sed $'s/^/594\t/g' >> /Users/armita/Downloads/vs_589_gq.tsv
```

```
library(ggplot2)
vs_589_dp <- read.delim("~/Downloads/vs_589_dp.tsv")
vs_589_dp$Isolate <- as.factor(vs_589_dp$Isolate)
ggplot(vs_589_dp, aes(x=Count, y=DP, color=Isolate)) + geom_line()
p1 <- ggplot(vs_589_dp, aes(x=Count, y=DP, color=Isolate)) + geom_line() + ylim(0,200)
ggsave("~/Downloads/vs_589_dp.pdf", p1)


vs_589_gq <- read.delim("~/Downloads/vs_589_gq.tsv")
vs_589_gq$Isolate <- as.factor(vs_589_gq$Isolate)
ggplot(vs_589_gq, aes(x=Count, y=GQ, color=Isolate)) + geom_line()
p2 <- ggplot(vs_589_gq, aes(x=Count, y=GQ, color=Isolate)) + geom_line() + ylim(0,100)
ggsave("~/Downloads/vs_589_gq.pdf", p2)
```



## Filter SNPs based on this region being present in all isolates

Only retain biallelic high-quality SNPS with no missing data (for any individual) for genetic analyses below (in some cases, may allow some missing data in order to retain more SNPs, or first remove poorly sequenced individuals with too much missing data and then filter the SNPs).

```bash
srun --partition=long --pty bash
Vcf=$(ls analysis_aa/popgen/SNP_calling/vs_589/589_contigs_softmasked_repeatmasker_TPSI_appended_temp.vcf)
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

# ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen/snp
# qsub $ProgDir/sub_vcf_parser.sh $Vcf 40 30 10 30 1 Y
```

```bash
# mv 62471_contigs_softmasked_repeatmasker_TPSI_appended_filtered.vcf analysis/popgen/SNP_calling/vs_62471/62471_contigs_softmasked_repeatmasker_TPSI_appended_filtered.vcf
mv R36_14_contigs_softmasked_repeatmasker_TPSI_appended_filtered.vcf analysis/popgen/SNP_calling/vs_R36_14/R36_14_contigs_softmasked_repeatmasker_TPSI_appended_filtered.vcf
```

```
After filtering, kept 20 out of 20 Individuals
Outputting VCF file...
After filtering, kept 65666 out of a possible 143123 Sites
Run Time = 9.00 seconds
```

## Gene models 589

Gene models were copied from the NRI cluster.

```bash
scp -r /Users/armita/OneDrive\ -\ University\ of\ Greenwich/Armitage/Projects/Pichia/fungap_out armita@149.155.34.73:/projects/oldhome/groups/harrisonlab/project_files/Pichia/gene_pred/.
```

```bash
cd /projects/oldhome/groups/harrisonlab/project_files/Pichia/
mkdir -p gene_pred/fungap/P.stipitis/589
mv gene_pred/fungap_out gene_pred/fungap/P.stipitis/589/.
```

# Collect VCF stats

General VCF stats (remember that vcftools needs to have the PERL library exported)

```bash
  Isolate="589"
  VcfTools=/projects/oldhome/sobczm/bin/vcftools/bin
  export PERL5LIB="$VcfTools:$PERL5LIB"
  VcfFiltered=$(ls analysis_aa/popgen/SNP_calling/vs_${Isolate}/*_qfiltered_presence*.vcf)
  Stats=$(echo $VcfFiltered | sed 's/.vcf/.stat/g')
  perl $VcfTools/vcf-stats $VcfFiltered > $Stats
```
Calculate the index for percentage of shared SNP alleles between the individuals.

```bash
  Isolate="589"
  for Vcf in $(ls analysis_aa/popgen/SNP_calling/vs_${Isolate}/*_qfiltered_presence*.vcf); do
      ProgDir=/projects/oldhome/armita/git_repos/emr_repos/scripts/popgen/snp
      $ProgDir/similarity_percentage.py $Vcf
  done
```

Visualise the output as heatmap and clustering dendrogram
```bash
Isolate="589"
for Log in $(ls analysis_aa/popgen/SNP_calling/vs_${Isolate}/*distance.log); do
  ProgDir=/projects/oldhome/armita/git_repos/emr_repos/scripts/popgen/snp
  Rscript --vanilla $ProgDir/distance_matrix.R $Log
  mv Rplots.pdf analysis_aa/popgen/SNP_calling/vs_${Isolate}/.
done
```

Identify SNPs in gene models:

Create custom SnpEff genome database

```bash
SnpEff=/projects/oldhome/sobczm/bin/snpEff
nano $SnpEff/snpEff.config
```
Add the following lines to the section with databases:

```bash
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
  for Vcf in $(ls analysis_aa/popgen/SNP_calling/vs_${Isolate}/*_qfiltered_presence*.recode.vcf); do
    echo $Vcf
    filename=$(basename "$Vcf")
    Prefix=${filename%.vcf}
    OutDir=$(ls -d analysis_aa/popgen/SNP_calling/vs_${Isolate})
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
All Gene CDS nonsyn syn
27	9	9	7	2
