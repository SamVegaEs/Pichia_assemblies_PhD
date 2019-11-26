```bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 16
#$ -l virtual_free=1G
#$ -l h=blacklace01.blacklace
# #$ -l h=blacklace01.blacklace|blacklace02.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace|blacklace12.blacklace

# Testing parallelisation of GATk HaplotypeCaller - may crash. (It did not! Resulted in 2x speedup)
# NOTE: this is a haploid organism. For diploid organism, change "ploidy" argument to 2.
# Changes required in the script:
# VARIABLES
# Reference - the genome reference used in read mapping.
# INSIDE THE GATK command:
# To specify which BAM mapping files (Out1 from pre_SNP_calling_cleanup.sh, RefName ending with "_rg" -> that is, with
# read group added) are to be used in SNP calling, use the -I argument with full path to each file following after that.
# Each new BAM file has to be specified after a separate -I

Reference=$1
Isolate=$2

# Project=/home/groups/harrisonlab/project_files/Pichia
Project=/home/groups/harrisonlab/project_files/Pichia
# OutDir=analysis/popgen/SNP_calling
OutDir=vs_${Isolate}
mkdir $OutDir
# Reference=$(ls /home/groups/harrisonlab/project_files/Pichia/repeat_masked/P.stipitis/589/filtered_contigs/589_contigs_unmasked.fa)
# Reference=$(ls /home/groups/harrisonlab/project_files/Pichia/repeat_masked/P.stipitis/589/filtereed_contigs/589_contigs_softmasked_repeatmasker_TPSI_appended.fa)


RefName=$(basename "$Reference")
Out1=$OutDir/"${RefName%.*}_temp.vcf"
Out2=$OutDir/"${RefName%.*}.vcf"

ProgDir=/home/sobczm/bin/GenomeAnalysisTK-3.6

java -jar $ProgDir/GenomeAnalysisTK.jar \
     -R $Reference \
     -T HaplotypeCaller \
     -ploidy 1 \
     -nct 15 \
     --allow_potentially_misencoded_quality_scores \
     -I $Project/analysis/popgen/P.stipitis/589/589_vs_589_unmasked_sorted_nomulti_proper_sorted_nodup_rg.bam \
     -I $Project/analysis/popgen/P.stipitis/589/589_vs_591_unmasked_sorted_nomulti_proper_sorted_nodup_rg.bam \
     -I $Project/analysis/popgen/P.stipitis/589/589_vs_594_unmasked_sorted_nomulti_proper_sorted_nodup_rg.bam
     -o $Out1

#Break down complex SNPs into primitive ones with VariantsToAllelicPrimitives
#This tool will take an MNP (e.g. ACCCA -> TCCCG) and break it up into separate records for each component part (A-T and A->G).

java -jar $ProgDir/GenomeAnalysisTK.jar \
   -T VariantsToAllelicPrimitives \
   -R $Reference \
   -V $Out1 \
   -o $Out2 \

#####################################
# Notes on GATK parallelisation
#####################################

# http://gatkforums.broadinstitute.org/gatk/discussion/1975/how-can-i-use-parallelism-to-make-gatk-tools-run-faster
(END)

```
