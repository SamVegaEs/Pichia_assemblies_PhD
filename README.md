# Pichia stipitis assembly and annotation.

commands used in the assembly and annotation of Pichia genomes

The program was run in:

cd /home/groups/harrisonlab/project_files/Pichia

The new folder after the change in the cluster:
/home/groups/harrisonlab/project_files/Pichia is now at cd /projects/oldhome/groups/harrisonlab/project_files/Pichia

My folder is now: /projects/oldhome/vegasa


# Demultiplexing reads using Albacore.

Data was basecalled again using Albacore 2.3.3 on the minion server:

```bash
screen -a  

#Sceen -a opens a new session.

ssh nanopore@nanopore

~/.local/bin/read_fast5_basecaller.py --version

mkdir Pichia_31_01_2019
cd Pichia_31_01_2019

# Oxford nanopore 07/03/17
Organism=P.stipitis
Date=31-01-19
FlowCell="FLO-MIN106"
Kit="SQK-LSK108"
RawDatDir=/data/seq_data/minion/2018/Pichia/GA20000/reads
OutDir=/data/scratch/nanopore_tmp_data/Pichia/albacore_v2.3.3
mkdir -p $OutDir

mkdir -p ~/Pichia_31_01_2019/$Date
cd ~/Pichia_31_01_2019/$Date
~/.local/bin/read_fast5_basecaller.py \
  --flowcell $FlowCell \
  --kit $Kit \
  --input $RawDatDir \
  --recursive \
  --worker_threads 24 \
  --save_path Pichia_albacore_v2.3.3_demultiplexed \
  --output_format fastq,fast5 \
  --reads_per_fastq_batch 4000 \
  --barcoding

#Run the commmand ls to check that the folders actually exist.

ls Pichia_albacore_v2.3.3_demultiplexed/workspace/pass/barcode01/*.fastq |wc -l

ls Pichia_albacore_v2.3.3_demultiplexed/workspace/pass/barcode02/*.fastq |wc -l

ls Pichia_albacore_v2.3.3_demultiplexed/workspace/pass/barcode03/*.fastq |wc -l

#Run each cat command individually and check that the output exists.

  cat Pichia_albacore_v2.3.3_demultiplexed/workspace/pass/barcode01/*.fastq | gzip -cf > Pichia_albacore_v2.3.3_barcode01.fastq.gz
  cat Pichia_albacore_v2.3.3_demultiplexed/workspace/pass/barcode02/*.fastq | gzip -cf > Pichia_albacore_v2.3.3_barcode02.fastq.gz
  cat Pichia_albacore_v2.3.3_demultiplexed/workspace/pass/barcode03/*.fastq | gzip -cf > Pichia_albacore_v2.3.3_barcode03.fastq.gz

#Set the space where the output data will be.

  OutDir=/data/scratch/nanopore_tmp_data/Pichia/albacore_v2.3.3
  mkdir -p $OutDir
  chmod +w $OutDir
  cp Pichia_albacore_v2.3.3_barcode0*.fastq.gz $OutDir/.
  chmod +rw $OutDir/Pichia_albacore_v2.3.3_barcode0*.fastq.gz
  # scp Alt_albacore_v2.10_barcode*.fastq.gz armita@192.168.1.200:$OutDir/.
  tar -cz -f Pichia_albacore_v2.3.3_demultiplexed.tar.gz Pichia_albacore_v2.3.3_demultiplexed
  OutDir=/data/scratch/nanopore_tmp_data/Pichia/albacore_v2.1.10
  mv Pichia_albacore_v2.3.3_demultiplexed.tar.gz $OutDir/.
  chmod +rw $OutDir/Pichia_albacore_v2.3.3_demultiplexed.tar.gz
```

# Building of directory structure

```bash
ProjDir=/home/groups/harrisonlab/project_files/Pichia
mkdir $ProjDir
cd $ProjDir

Organism=P.stipitis
Strain=589
OutDir=raw_dna/minion/$Organism/$Strain
mkdir -p $OutDir
RawDat=$(ls /data/scratch/nanopore_tmp_data/Pichia/albacore_v2.3.3/Pichia_albacore_v2.3.3_barcode01.fastq.gz)
cd $OutDir
cp -s $RawDat .
cd $ProjDir

Organism=P.stipitis
Strain=591
OutDir=raw_dna/minion/$Organism/$Strain
mkdir -p $OutDir
RawDat=$(ls /data/scratch/nanopore_tmp_data/Pichia/albacore_v2.3.3/Pichia_albacore_v2.3.3_barcode02.fastq.gz)
cd $OutDir
cp -s $RawDat .
cd $ProjDir

Organism=P.stipitis
Strain=594
OutDir=raw_dna/minion/$Organism/$Strain
mkdir -p $OutDir
RawDat=$(ls /data/scratch/nanopore_tmp_data/Pichia/albacore_v2.3.3/Pichia_albacore_v2.3.3_barcode03.fastq.gz)
cd $OutDir
cp -s $RawDat .
cd $ProjDir

```

Doing it in in this way was not possible since the message of no space left in the device appeared. So we will try another way: First running the OutDir variable and then running the tar command. So:

```bash
  OutDir=/data/scratch/nanopore_tmp_data/Alternaria/albacore_v2.1.10
  mv Pichia_albacore_v2.3.3_demultiplexed.tar.gz $OutDir/.
  chmod +rw $OutDir/Pichia_albacore_v2.3.3_demultiplexed.tar.gz

  tar -cz -f $OutDir/Pichia_albacore_v2.3.3_demultiplexed.tar.gz Pichia_albacore_v2.3.3_demultiplexed

```

# ASSEMBLY.

1. Removal of the adapters: Splitting reads and trimming adapters using porechop. We only need the FASTA files for this step, so we can do it while the "tar" command is runing. To do so, we have to create the directory in our project folder and not the nanopore node.


```bash
for RawReads in $(ls raw_dna/minion/*/*/*.fastq.gz); do
    Organism=$(echo $RawReads| rev | cut -f3 -d '/' | rev)
    Strain=$(echo $RawReads | rev | cut -f2 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=qc_dna/minion/$Organism/$Strain
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
    qsub $ProgDir/sub_porechop.sh $RawReads $OutDir
  done
```

2. Identify sequencing coverage

#For Minion data:

```bash
for RawData in $(ls qc_dna/minion/*/*/*q.gz); do
echo $RawData;
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
GenomeSz=16
OutDir=$(dirname $RawData)
mkdir -p $OutDir
qsub $ProgDir/sub_count_nuc.sh $GenomeSz $RawData $OutDir
done


  for StrainDir in $(ls -d qc_dna/minion/*/*); do
    Strain=$(basename $StrainDir)
    printf "$Strain\t"
    for File in $(ls $StrainDir/*.txt); do
      echo $(basename $File);
      cat $File | tail -n1 | rev | cut -f2 -d ' ' | rev;
    done | grep -v '.txt' | awk '{ SUM += $1} END { print SUM }'
  done
```

```
Coverage.

589     186.88
591     145.17
594     158.09
```


3. Read correction using Canu

```bash
for TrimReads in $(ls qc_dna/minion/*/*/*q.gz); do
Organism=$(echo $TrimReads | rev | cut -f3 -d '/' | rev)
Strain=$(echo $TrimReads | rev | cut -f2 -d '/' | rev)
OutDir=assembly/canu-1.8/$Organism/"$Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/canu
qsub $ProgDir/sub_canu_correction.sh $TrimReads 16m $Strain $OutDir
done
```

4. Assembbly using SMARTdenovo

```bash
for CorrectedReads in $(ls assembly/canu-1.8/*/*/*.trimmedReads.fasta.gz); do
Organism=$(echo $CorrectedReads | rev | cut -f3 -d '/' | rev)
Strain=$(echo $CorrectedReads | rev | cut -f2 -d '/' | rev)
Prefix="$Strain"_smartdenovo
OutDir=assembly/SMARTdenovo/$Organism/"$Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/SMARTdenovo
qsub $ProgDir/sub_SMARTdenovo.sh $CorrectedReads $Prefix $OutDir
done
```

5. Quast and busco were run to assess the effects of racon on assembly quality:

```bash

ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/SMARTdenovo/*/*/*.dmo.lay.utg); do
  Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/*.dmo.lay.utg); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/ascomycota_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

6. Error correction using racon:

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/*.dmo.lay.utg); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
ReadsFq=$(ls ../../../../../home/groups/harrisonlab/project_files/Pichia/qc_dna/minion/*/$Strain/*q.gz)
# ReadsFq=$(ls qc_dna/minion/*/$Strain/*q.gz)
Iterations=10
OutDir=$(dirname $Assembly)"/racon2_$Iterations"
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/racon
qsub $ProgDir/sub_racon.sh $Assembly $ReadsFq $Iterations $OutDir
done
```
```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/SMARTdenovo/P.*/*/racon2_10/*.fasta | grep 'round_10'); do
OutDir=$(dirname $Assembly)
echo "" > tmp.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
$ProgDir/remove_contaminants.py --keep_mitochondria --inp $Assembly --out $OutDir/racon_min_500bp_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
done
```

7. Quast and busco were run to assess the effects of racon on assembly quality:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/SMARTdenovo/P.*/*/racon2_10/racon_min_500bp_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```
```bash
# for Assembly in $(ls assembly/SMARTdenovo/F.*/*/racon*/*.fasta | grep 'FON_63' | grep 'racon_min_500bp_renamed'); do
for Assembly in $(ls assembly/SMARTdenovo/P.*/*/racon2_10/*.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/ascomycota_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly
# OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

```bash
printf "Filename\tComplete\tDuplicated\tFragmented\tMissing\tTotal\n"
for File in $(ls gene_pred/busco/P*/*/assembly/*/short_summary_*.txt); do
FileName=$(basename $File)
Complete=$(cat $File | grep "(C)" | cut -f2)
Duplicated=$(cat $File | grep "(D)" | cut -f2)
Fragmented=$(cat $File | grep "(F)" | cut -f2)
Missing=$(cat $File | grep "(M)" | cut -f2)
Total=$(cat $File | grep "Total" | cut -f2)
printf "$FileName\t$Complete\t$Duplicated\t$Fragmented\t$Missing\t$Total\n"
done
```

Busco Results.

```bash
short_summary_589_smartdenovo.dmo.lay.txt       531     13      306     478    1                                                                                        315
short_summary_589_smartdenovo_racon_round_10.txt        0       0       0      1                                                                                        315     1315
short_summary_589_smartdenovo_racon_round_1.txt 0       0       0       1315   1                                                                                        315
short_summary_589_smartdenovo_racon_round_2.txt 0       0       0       1315   1                                                                                        315
short_summary_589_smartdenovo_racon_round_3.txt 0       0       0       1315   1                                                                                        315
short_summary_589_smartdenovo_racon_round_6.txt 970     4       150     195    1                                                                                        315
short_summary_589_smartdenovo_racon_round_7.txt 0       0       0       1315   1                                                                                        315
short_summary_589_smartdenovo_racon_round_8.txt 0       0       0       1315   1                                                                                        315
short_summary_589_smartdenovo_racon_round_9.txt 0       0       0       1315   1                                                                                        315
short_summary_racon_min_500bp_renamed.txt       0       0       0       1315   1                                                                                        315
short_summary_591_smartdenovo.dmo.lay.txt       501     6       321     493    1                                                                                        315
short_summary_591_smartdenovo_racon_round_10.txt        0       0       0      1                                                                                        315     1315
short_summary_591_smartdenovo_racon_round_1.txt 0       0       0       1315   1                                                                                        315
short_summary_591_smartdenovo_racon_round_2.txt 0       0       0       1315   1                                                                                        315
short_summary_591_smartdenovo_racon_round_3.txt 0       0       0       1315   1                                                                                        315
short_summary_racon_min_500bp_renamed.txt       0       0       0       1315   1                                                                                        315
short_summary_594_smartdenovo.dmo.lay.txt       505     1       321     489    1                                                                                        315
short_summary_594_smartdenovo_racon_round_10.txt        0       0       0      1                                                                                        315     1315
short_summary_594_smartdenovo_racon_round_1.txt 0       0       0       1315   1                                                                                        315
short_summary_594_smartdenovo_racon_round_2.txt 0       0       0       1315   1                                                                                        315
short_summary_594_smartdenovo_racon_round_3.txt 0       0       0       1315   1                                                                                        315
short_summary_594_smartdenovo_racon_round_8.txt 0       0       0       1315   1                                                                                        315
short_summary_594_smartdenovo_racon_round_9.txt 0       0       0       1315   1                                                                                        315
short_summary_racon_min_500bp_renamed.txt       0       0       0       1315   1                                                                                        315

```

This results do not seem correct since some files have the value "0" for complete genes, that means that there are 0 genes identified. I checked the files after the run in racon, and they were empty, so I re-run the command. After running it I should check if there is information in them with ls -lh. After repeating it, the files had information on them, and BUSCO was run again. The results obtained after the second run are:

```bash
short_summary_589_smartdenovo.dmo.lay.txt         531     13      306     478    1315
short_summary_589_smartdenovo_racon_round_10.txt  975     5       158     182    1315
short_summary_589_smartdenovo_racon_round_1.txt   945     9       174     196    1315
short_summary_589_smartdenovo_racon_round_2.txt   949     13      169     197    1315
short_summary_589_smartdenovo_racon_round_3.txt   989     9       151     175    1315
short_summary_589_smartdenovo_racon_round_4.txt   977     10      142     196    1315
short_summary_589_smartdenovo_racon_round_5.txt   971     8       160     184    1315
short_summary_589_smartdenovo_racon_round_6.txt   970     4       150     195    1315
short_summary_589_smartdenovo_racon_round_7.txt   984     5       149     182    1315
short_summary_589_smartdenovo_racon_round_8.txt   970     5       145     200    1315
short_summary_589_smartdenovo_racon_round_9.txt   972     4       155     188    1315
short_summary_racon_min_500bp_renamed.txt         975     5       158     182    1315

short_summary_591_smartdenovo.dmo.lay.txt         501     6       321     493    1315
short_summary_591_smartdenovo_racon_round_10.txt  957     4       163     195    1315
short_summary_591_smartdenovo_racon_round_1.txt   944     3       159     212    1315
short_summary_591_smartdenovo_racon_round_2.txt   957     8       160     198    1315
short_summary_591_smartdenovo_racon_round_3.txt   960     7       170     185    1315
short_summary_591_smartdenovo_racon_round_4.txt   962     4       162     191    1315
short_summary_591_smartdenovo_racon_round_5.txt   973     4       161     181    1315
short_summary_591_smartdenovo_racon_round_6.txt   955     4       168     192    1315
short_summary_591_smartdenovo_racon_round_7.txt   962     4       160     193    1315
short_summary_591_smartdenovo_racon_round_8.txt   967     4       155     193    1315
short_summary_591_smartdenovo_racon_round_9.txt   973     4       149     193    1315
short_summary_racon_min_500bp_renamed.txt         957     4       163     195    1315


short_summary_594_smartdenovo.dmo.lay.txt         505     1       321     489    1315
short_summary_594_smartdenovo_racon_round_10.txt  950     1       170     195    1315
short_summary_594_smartdenovo_racon_round_1.txt   936     0       182     197    1315
short_summary_594_smartdenovo_racon_round_2.txt   949     1       166     200    1315
short_summary_594_smartdenovo_racon_round_3.txt   948     1       175     192    1315
short_summary_594_smartdenovo_racon_round_4.txt   951     0       178     186    1315
short_summary_594_smartdenovo_racon_round_5.txt   961     1       168     186    1315
short_summary_594_smartdenovo_racon_round_6.txt   949     2       173     193    1315
short_summary_594_smartdenovo_racon_round_7.txt   948     1       172     195    1315
short_summary_594_smartdenovo_racon_round_8.txt   973     0       165     177    1315
short_summary_594_smartdenovo_racon_round_9.txt   946     0       178     191    1315
short_summary_racon_min_500bp_renamed.txt         950     1       170     195    1315
```

# ASSEMBLY CORRECTION USING NANOPOLISH.

```bash
screen -a

ssh nanopore@nanopore

cd /home/nanopore/Pichia_31_01_2019/31-01-19/
mv Pichia_albacore_v2.3.3_demultiplexed.tar.gz $OutDir/.
chmod +rw $OutDir/Pichia_albacore_v2.3.3_demultiplexed.tar.gz

tar -cz -f $OutDir/Pichia_albacore_v2.3.3_demultiplexed.tar.gz Pichia_albacore_v2.3.3_demultiplexed
```
After problems running the tar command we ran it again from inside the nanopore node.

To check if the tar command is still running we have to resume the screen in which it was running.

Fast5 files are very large and need to be stored as gzipped tarballs. These needed temporarily unpacking but must be deleted after nanpolish has finished running.

Raw reads were moved onto the cluster scratch space for this step and unpacked:

```bash
screen -a

for Tar in $(ls /data/scratch/nanopore_tmp_data/Alternaria/albacore_v2.1.10/Pichia_albacore_v2.3.3_demultiplexed.tar.gz); do
  ScratchDir=/data2/scratch2/vegasa/pichia/albacore_v2.3.3
  mkdir -p $ScratchDir
  tar -zxvf $Tar -C $ScratchDir
done
```

When you try to run this part in a screen session it asks for an update in certain programs required to run nanopolish, so do not run it in screen.

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/589/racon2_10/racon_min_500bp_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
# Step 1 extract reads as a .fq file which contain info on the location of the fast5 files
# Note - the full path from home must be used
ReadDir=raw_dna/nanopolish/$Organism/$Strain
mkdir -p $ReadDir
ReadsFq=$(ls ../../../../../home/groups/harrisonlab/project_files/Pichia/raw_dna/minion/*/$Strain/*.fastq.gz)
ScratchDir=/data2/scratch2/vegasa/pichia/albacore_v2.3.3
Fast5Dir=$ScratchDir/Pichia_albacore_v2.3.3_demultiplexed/workspace/pass/barcode01
nanopolish index -d $Fast5Dir $ReadsFq
OutDir=$(dirname $Assembly)/nanopolish
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/nanopolish
# submit alignments for nanoppolish
# qsub $ProgDir/sub_minimap2_nanopolish.sh $Assembly $ReadsFq $OutDir/nanopolish
done

for Assembly in $(ls assembly/SMARTdenovo/*/591/racon2_10/racon_min_500bp_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
# Step 1 extract reads as a .fq file which contain info on the location of the fast5 files
# Note - the full path from home must be used
ReadDir=raw_dna/nanopolish/$Organism/$Strain
mkdir -p $ReadDir
ReadsFq=$(ls ../../../../../home/groups/harrisonlab/project_files/Pichia/raw_dna/minion/*/$Strain/*.fastq.gz)
ScratchDir=/data2/scratch2/vegasa/pichia/albacore_v2.3.3
Fast5Dir=$ScratchDir/Pichia_albacore_v2.3.3_demultiplexed/workspace/pass/barcode02
nanopolish index -d $Fast5Dir $ReadsFq
OutDir=$(dirname $Assembly)
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/nanopolish
# submit alignments for nanoppolish
# qsub $ProgDir/sub_minimap2_nanopolish.sh $Assembly $ReadsFq $OutDir/nanopolish
done

for Assembly in $(ls assembly/SMARTdenovo/*/594/racon2_10/racon_min_500bp_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
# Step 1 extract reads as a .fq file which contain info on the location of the fast5 files
# Note - the full path from home must be used
ReadDir=raw_dna/nanopolish/$Organism/$Strain
mkdir -p $ReadDir
ReadsFq=$(ls ../../../../../home/groups/harrisonlab/project_files/Pichia/raw_dna/minion/*/$Strain/*.fastq.gz)
ScratchDir=/data2/scratch2/vegasa/pichia/albacore_v2.3.3
Fast5Dir=$ScratchDir/Pichia_albacore_v2.3.3_demultiplexed/workspace/pass/barcode03
nanopolish index -d $Fast5Dir $ReadsFq
OutDir=$(dirname $Assembly)
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/nanopolish
# submit alignments for nanoppolish
# qsub $ProgDir/sub_minimap2_nanopolish.sh $Assembly $ReadsFq $OutDir/nanopolish
done
```

1. Split the assembly into 50Kb fragments an submit each to the cluster for nanopolish correction

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/racon2_10/racon_min_500bp_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=$(dirname $Assembly)/nanopolish
RawReads=$(ls ../../../../../home/groups/harrisonlab/project_files/Pichia/raw_dna/minion/*/$Strain/*.fastq.gz)
AlignedReads=$(ls $OutDir/reads.sorted.bam)

NanoPolishDir=/home/armita/prog/nanopolish/nanopolish/scripts
python $NanoPolishDir/nanopolish_makerange.py $Assembly --segment-length 50000 > $OutDir/nanopolish_range.txt

Ploidy=1
echo "nanopolish log:" > $OutDir/nanopolish_log.txt
for Region in $(cat $OutDir/nanopolish_range.txt); do
Jobs=$(qstat | grep 'sub_nanopo' | grep 'qw' | wc -l)
while [ $Jobs -gt 1 ]; do
sleep 1m
printf "."
Jobs=$(qstat | grep 'sub_nanopo' | grep 'qw' | wc -l)
done    
printf "\n"
echo $Region
echo $Region >> $OutDir/nanopolish_log.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/nanopolish
qsub $ProgDir/sub_nanopolish_variants.sh $Assembly $RawReads $AlignedReads $Ploidy $Region $OutDir/$Region
done
done
```

1.2 When the code for splitting was running we had to stop it. To take it from where it stopped and not from the beginning:

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/racon2_10/racon_min_500bp_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=$(dirname $Assembly)/nanopolish
RawReads=$(ls ../../../../../home/groups/harrisonlab/project_files/Pichia/raw_dna/minion/*/$Strain/*.fastq.gz)
AlignedReads=$(ls $OutDir/reads.sorted.bam)

NanoPolishDir=/home/armita/prog/nanopolish/nanopolish/scripts
python $NanoPolishDir/nanopolish_makerange.py $Assembly --segment-length 50000 > $OutDir/nanopolish_range.txt

Ploidy=1
echo "nanopolish log:" > $OutDir/nanopolish_log2.txt
ls -lh $OutDir/*/*.fa | grep -v ' 0 ' | cut -f8 -d '/' | sed 's/_consensus.fa//g' > $OutDir/files_present.txt
for Region in $(cat $OutDir/nanopolish_range.txt | grep -vwf "$OutDir/files_present.txt"); do
Jobs=$(qstat | grep 'sub_nanopo' | grep 'qw' | wc -l)
while [ $Jobs -gt 1 ]; do
sleep 1m
printf "."
Jobs=$(qstat | grep 'sub_nanopo' | grep 'qw' | wc -l)
done    
printf "\n"
echo $Region
echo $Region >> $OutDir/nanopolish_log2.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/nanopolish
qsub $ProgDir/sub_nanopolish_variants.sh $Assembly $RawReads $AlignedReads $Ploidy $Region $OutDir/$Region
done
done
```

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/racon2_10/racon_min_500bp_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
OutDir=assembly/SMARTdenovo/$Organism/$Strain/nanopolish
mkdir -p $OutDir
# cat "" > $OutDir/"$Strain"_nanoplish.fa
NanoPolishDir=/home/armita/prog/nanopolish/nanopolish/scripts
InDir=$(dirname $Assembly)
python $NanoPolishDir/nanopolish_merge.py $InDir/nanopolish/*/*.fa > $OutDir/"$Strain"_nanoplish.fa

echo "" > tmp.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
$ProgDir/remove_contaminants.py --keep_mitochondria --inp $OutDir/"$Strain"_nanoplish.fa --out $OutDir/"$Strain"_nanoplish_min_500bp_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
done
```

2. Quast and busco were run to assess the effects of nanopolish on assembly quality:

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/nanopolish/*_nanoplish_min_500bp_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
# Quast
OutDir=$(dirname $Assembly)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
# Busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/ascomycota_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```
3. For Pilon assembly correction we need the Illumina data, and we don't have them yet in our working directory. We need to copy them there.

```bash
RawDatDir=/data/seq_data/miseq/2018/RAW/181026_M04465_0088_000000000-C5M4V/Data/Intensities/BaseCalls

  ProjectDir=/home/groups/harrisonlab/project_files/Pichia

# Isolate 589

  OutDir=$ProjectDir/raw_dna/paired/P.stipitis/589

  mkdir -p $OutDir/F

  mkdir -p $OutDir/R

  cd $OutDir/F

  cp -s $RawDatDir/Pichia-589_S2_L001_R1_001.fastq.gz .

  cd $OutDir/R

  cp -s $RawDatDir/Pichia-589_S2_L001_R2_001.fastq.gz  .

  cd $ProjectDir

# Isolate 591

  OutDir=$ProjectDir/raw_dna/paired/P.stipitis/591

  mkdir -p $OutDir/F

  mkdir -p $OutDir/R

  cd $OutDir/F

  cp -s $RawDatDir/Pichia-591_S3_L001_R1_001.fastq.gz .

  cd $OutDir/R

  cp -s $RawDatDir/Pichia-591_S3_L001_R2_001.fastq.gz  .

  cd $ProjectDir

# Isolate 594

  OutDir=$ProjectDir/raw_dna/paired/P.stipitis/594

  mkdir -p $OutDir/F

  mkdir -p $OutDir/R

  cd $OutDir/F

  cp -s $RawDatDir/Pichia-594_S4_L001_R1_001.fastq.gz .

  cd $OutDir/R

  cp -s $RawDatDir/Pichia-594_S4_L001_R2_001.fastq.gz  .

  cd $ProjectDir
```

4. QC of MiSeq Data.

programs: fastqc fastq-mcf kmc

Data quality was visualised using fastqc:

```bash
for RawData in $(ls raw_dna/paired/*/*/*/*.fastq.gz); do
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
    echo $RawData;
    qsub $ProgDir/run_fastqc.sh $RawData
  done
```
```bash
for StrainPath in $(ls -d raw_dna/paired/*/*); do
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
    IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
    ReadsF=$(ls $StrainPath/F/*.fastq*)
    ReadsR=$(ls $StrainPath/R/*.fastq*)
    echo $ReadsF
    echo $ReadsR
    qsub $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters DNA
  done
```


Calculate sequnce coverage for illumina data:

```bash
for RawData in $(ls qc_dna/paired/*/*/*/*q.gz | grep -v 'appended'); do
echo $RawData;
GenomeSz=16
OutDir=$(dirname $RawData | sed 's/paired/tmp/g')
mkdir -p $OutDir
ProgDir=/projects/oldhome/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
sbatch $ProgDir/slurm_count_nuc.sh $GenomeSz $RawData $OutDir
done


  for StrainDir in $(ls -d qc_dna/tmp/*/*); do
    Strain=$(basename $StrainDir)
    printf "$Strain\t"
    for File in $(ls $StrainDir/*/*.txt); do
      echo $(basename $File);
      cat $File | tail -n1 | rev | cut -f2 -d ' ' | rev;
    done | grep -v '.txt' | awk '{ SUM += $1} END { print SUM }'
  done
```

```
589	102.76
591	76.43
594	106.91

```

5. Pilon assembly correction: Assemblies were polished using Pilon

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/nanopolish/*_nanoplish_min_500bp_renamed.fasta); do
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
IlluminaDir=$(ls -d ../../../../../home/groups/harrisonlab/project_files/Pichia/qc_dna/paired/*/$Strain)
TrimF1_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n1)
TrimR1_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n1)
OutDir=$(dirname $Assembly)/../pilon
Iterations=10
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
qsub $ProgDir/sub_pilon.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir $Iterations
done
```

6. Contigs were renamed.

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/pilon/*.fasta | grep 'pilon_10'); do
echo $Assembly
echo "" > tmp.txt
OutDir=$(dirname $Assembly)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
$ProgDir/remove_contaminants.py --keep_mitochondria --inp $Assembly --out $OutDir/pilon_min_500bp_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
done
```

7. Quast and busco were run to assess the effects of pilon on assembly quality:

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/pilon/*.fasta); do
  Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
  echo "$Organism - $Strain"
  # OutDir=$(dirname $Assembly)
  # ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  # qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  Jobs=$(qstat | grep 'sub_busco' | grep 'qw'| wc -l)
  while [ $Jobs -gt 1 ]; do
  sleep 1m
  printf "."
  Jobs=$(qstat | grep 'sub_busco' | grep 'qw'| wc -l)
  done
  printf "\n"
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
  BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/ascomycota_odb9)
  # OutDir=gene_pred/busco/$Organism/$Strain/assembly
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

```bash
printf "Filename\tComplete\tDuplicated\tFragmented\tMissing\tTotal\n"
for File in $(ls assembly/SMARTdenovo/*/*/pilon/*/short_summary_*.txt); do  
FileName=$(basename $File)
Complete=$(cat $File | grep "(C)" | cut -f2)
Duplicated=$(cat $File | grep "(D)" | cut -f2)
Fragmented=$(cat $File | grep "(F)" | cut -f2)
Missing=$(cat $File | grep "(M)" | cut -f2)
Total=$(cat $File | grep "Total" | cut -f2)
printf "$FileName\t$Complete\t$Duplicated\t$Fragmented\t$Missing\t$Total\n"
done
```

Summary results of Busco run:

```bash
short_summary_pilon_10.txt                      1290    11      5       20      1315
short_summary_pilon_1.txt                       1290    11      5       20      1315
short_summary_pilon_2.txt                       1290    11      5       20      1315
short_summary_pilon_3.txt                       1290    11      5       20      1315
short_summary_pilon_4.txt                       1290    11      5       20      1315
short_summary_pilon_5.txt                       1290    11      5       20      1315
short_summary_pilon_6.txt                       1290    11      5       20      1315
short_summary_pilon_7.txt                       1290    11      5       20      1315
short_summary_pilon_8.txt                       1290    11      5       20      1315
short_summary_pilon_9.txt                       1290    11      5       20      1315
short_summary_pilon_min_500bp_renamed.txt       1290    11      5       20      1315

short_summary_pilon_10.txt                      1284    3       7       24      1315
short_summary_pilon_1.txt                       1284    3       6       25      1315
short_summary_pilon_2.txt                       1285    3       6       24      1315
short_summary_pilon_3.txt                       1284    3       7       24      1315
short_summary_pilon_4.txt                       1284    3       7       24      1315
short_summary_pilon_5.txt                       1284    3       7       24      1315
short_summary_pilon_6.txt                       1284    3       7       24      1315
short_summary_pilon_7.txt                       1284    3       7       24      1315
short_summary_pilon_8.txt                       1284    3       7       24      1315
short_summary_pilon_9.txt                       1284    3       7       24      1315
short_summary_pilon_min_500bp_renamed.txt       1284    3       7       24      1315

short_summary_pilon_10.txt                      1288    1       5       22      1315
short_summary_pilon_1.txt                       1288    1       5       22      1315
short_summary_pilon_2.txt                       1288    1       5       22      1315
short_summary_pilon_3.txt                       1288    1       5       22      1315
short_summary_pilon_4.txt                       1288    1       5       22      1315
short_summary_pilon_5.txt                       1288    1       5       22      1315
short_summary_pilon_6.txt                       1288    1       5       22      1315
short_summary_pilon_7.txt                       1288    1       5       22      1315
short_summary_pilon_8.txt                       1288    1       5       22      1315
short_summary_pilon_9.txt                       1288    1       5       22      1315
short_summary_pilon_min_500bp_renamed.txt       1288    1       5       22      1315
```

To get the information about the total length of the genome, and the number of contigs we have we use this command:  

```bash
for File in $(ls assembly/SMARTdenovo/P.stipitis/*/racon2_10/report.txt); do
echo $File
cat $File
done

589

All statistics are based on contigs of size >= 500 bp, unless otherwise noted (e.g., "# contigs (>= 0 bp)" and "Total length (>= 0 bp)" include all contigs).

Assembly                   racon_min_500bp_renamed
# contigs (>= 0 bp)        11
# contigs (>= 1000 bp)     11
Total length (>= 0 bp)     15687345
Total length (>= 1000 bp)  15687345
# contigs                  11
Largest contig             2692878
Total length               15687345
GC (%)                     41.02
N50                        1881524
N75                        1677956
L50                        4
L75                        6
# N's per 100 kbp          0.00

591

Assembly                   racon_min_500bp_renamed
# contigs (>= 0 bp)        11
# contigs (>= 1000 bp)     11
Total length (>= 0 bp)     15639864
Total length (>= 1000 bp)  15639864
# contigs                  11
Largest contig             3983890
Total length               15639864
GC (%)                     41.11
N50                        1904357
N75                        1319906
L50                        3
L75                        6
# N's per 100 kbp          0.00

594

Assembly                   racon_min_500bp_renamed
# contigs (>= 0 bp)        9
# contigs (>= 1000 bp)     9
Total length (>= 0 bp)     15358923
Total length (>= 1000 bp)  15358923
# contigs                  9
Largest contig             3438607
Total length               15358923
GC (%)                     41.02
N50                        1898369
N75                        1767446
L50                        3
L75                        5
# N's per 100 kbp          0.00
```

# HYBRID ASSEMBLY.

1. Spades Assembly

```bash
for TrimReads in $(ls ../../../../../home/groups/harrisonlab/project_files/Pichia/raw_dna/minion/*/*/*.fastq.gz); do
Organism=$(echo $TrimReads | rev | cut -f3 -d '/' | rev)
Strain=$(echo $TrimReads | rev | cut -f2 -d '/' | rev)
IlluminaDir=$(ls -d ../../../../../home/groups/harrisonlab/project_files/Pichia/qc_dna/paired/*/$Strain)
TrimF1_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n1)
TrimR1_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n1)
OutDir=assembly/spades_minion/$Organism/"$Strain"
echo $TrimF1_Read
echo $TrimR1_Read
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades
qsub $ProgDir/sub_spades_minion.sh $TrimReads $TrimF1_Read $TrimR1_Read $OutDir
done
```

2. Contigs shorter than 500bp were removed from the assembly.

```bash
for Contigs in $(ls assembly/spades_minion/*/*/contigs.fasta); do
    AssemblyDir=$(dirname $Contigs)
    mkdir $AssemblyDir/filtered_contigs
    FilterDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/abyss
    $FilterDir/filter_abyss_contigs.py $Contigs 500 > $AssemblyDir/filtered_contigs/contigs_min_500bp.fasta
  done
```

3. Quast and BUSCO

```bash
for Assembly in $(ls assembly/spades_minion/*/*/filtered_contigs/contigs_min_500bp.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
OutDir=$(dirname $Assembly)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
# BuscoDB="Fungal"
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/ascomycota_odb9)
OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

To do yet:

```bash
 for File in $(ls assembly/spades_minion/*/*/*/polished/*/short_summary_*.txt ); do
  Strain=$(echo $File| rev | cut -d '/' -f5 | rev)
  Organism=$(echo $File | rev | cut -d '/' -f6 | rev)
  Prefix=$(basename $File)
  Complete=$(cat $File | grep "(C)" | cut -f2)
  Fragmented=$(cat $File | grep "(F)" | cut -f2)
  Missing=$(cat $File | grep "(M)" | cut -f2)
  Total=$(cat $File | grep "Total" | cut -f2)
  echo -e "$Organism\t$Strain\t$Prefix\t$Complete\t$Fragmented\t$Missing\t$Total"
  done
  ```

The folder polished is not present, so I have to create it before going on with the Merging of both assemblies.

NOTE: HYBRID ASSEMBBLY was stopped at this point. The rest is yet to do in case we are interested in the future.

# Re-run BUSCO

Repeat the run of Quast and busco were run to assess the effects of pilon on assembly quality, but in this case using the data from saccharomycetales_odb9, since it's Pichia's order.

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/pilon/*.fasta); do
  Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
  echo "$Organism - $Strain"
  # OutDir=$(dirname $Assembly)
  # ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  # qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  Jobs=$(qstat | grep 'sub_busco' | grep 'qw'| wc -l)
  while [ $Jobs -gt 1 ]; do
  sleep 1m
  printf "."
  Jobs=$(qstat | grep 'sub_busco' | grep 'qw'| wc -l)
  done
  printf "\n"
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
  BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/saccharomycetales_odb9)
  # OutDir=gene_pred/busco/$Organism/$Strain/assembly
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

```bash
printf "Filename\tComplete\tDuplicated\tFragmented\tMissing\tTotal\n"
for File in $(ls assembly/SMARTdenovo/*/*/pilon/*/short_summary_*.txt); do  
FileName=$(basename $File)
Complete=$(cat $File | grep "(C)" | cut -f2)
Duplicated=$(cat $File | grep "(D)" | cut -f2)
Fragmented=$(cat $File | grep "(F)" | cut -f2)
Missing=$(cat $File | grep "(M)" | cut -f2)
Total=$(cat $File | grep "Total" | cut -f2)
printf "$FileName\t$Complete\t$Duplicated\t$Fragmented\t$Missing\t$Total\n"
done
```

Summary results:

```bash
short_summary_pilon_10.txt                      1683    13      11      17      1711
short_summary_pilon_1.txt                       1682    13      12      17      1711
short_summary_pilon_2.txt                       1683    13      11      17      1711
short_summary_pilon_3.txt                       1683    13      11      17      1711
short_summary_pilon_4.txt                       1683    13      11      17      1711
short_summary_pilon_5.txt                       1683    13      11      17      1711
short_summary_pilon_6.txt                       1683    13      11      17      1711
short_summary_pilon_7.txt                       1683    13      11      17      1711
short_summary_pilon_8.txt                       1683    13      11      17      1711
short_summary_pilon_9.txt                       1683    13      11      17      1711
short_summary_pilon_min_500bp_renamed.txt       1683    13      11      17      1711

short_summary_pilon_10.txt                      1678    9       14      19      1711
short_summary_pilon_1.txt                       1676    9       15      20      1711
short_summary_pilon_2.txt                       1678    9       13      20      1711
short_summary_pilon_3.txt                       1678    9       14      19      1711
short_summary_pilon_4.txt                       1678    9       14      19      1711
short_summary_pilon_5.txt                       1678    9       14      19      1711
short_summary_pilon_6.txt                       1678    9       14      19      1711
short_summary_pilon_7.txt                       1678    9       14      19      1711
short_summary_pilon_8.txt                       1678    9       14      19      1711
short_summary_pilon_9.txt                       1678    9       14      19      1711
short_summary_pilon_min_500bp_renamed.txt       1678    9       14      19      1711

short_summary_pilon_10.txt                      1685    5       12      14      1711
short_summary_pilon_1.txt                       1684    5       13      14      1711
short_summary_pilon_2.txt                       1685    5       12      14      1711
short_summary_pilon_3.txt                       1685    5       12      14      1711
short_summary_pilon_4.txt                       1685    5       12      14      1711
short_summary_pilon_5.txt                       1685    5       12      14      1711
short_summary_pilon_6.txt                       1685    5       12      14      1711
short_summary_pilon_7.txt                       1685    5       12      14      1711
short_summary_pilon_8.txt                       1685    5       12      14      1711
short_summary_pilon_9.txt                       1685    5       12      14      1711
short_summary_pilon_min_500bp_renamed.txt       1685    5       12      14      1711
```

# Removal of mit. DNA

1. Using an exclusion database with deconseq:

```bash
 for Assembly in $(ls assembly/SMARTdenovo/*/*/pilon/pilon_min_500bp_renamed.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    for Exclude_db in "Calb_mtDNA"; do
      AssemblyDir=$(dirname $Assembly)
      OutDir=$AssemblyDir/../deconseq_$Exclude_db
      ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
      qsub $ProgDir/sub_deconseq_no_retain.sh $Assembly $Exclude_db $OutDir
    done
  done
```

Results were summarised using the commands:

```bash
for Exclude_db in "Calb_mtDNA"; do
echo $Exclude_db
for File in $(ls assembly/*/*/*/*/log.txt | grep "$Exclude_db"); do
Name=$(echo $File | rev | cut -f3 -d '/' | rev);
Good=$(cat $File |cut -f2 | head -n1 | tail -n1);
Bad=$(cat $File |cut -f2 | head -n3 | tail -n1);
printf "$Name\t$Good\t$Bad\n";
done
done
```

```bash
Calb_mtDNA
589     11      0
591     11      0
594     8       1
```

2. Quast was run on the removed mtDNA:

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/deconseq_Calb_mtDNA/*_cont.fa); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
OutDir=$(dirname $Assembly)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

# Repeat Masking.

Repeat masking was done for the Nanopore assembly polished with Illumina data (SMARTdenovo).

```bash
 # for Assembly in $(ls assembly/SMARTdenovo/*/*/pilon/pilon_min_500bp_renamed.fasta); do
  for Assembly in $(ls assembly/SMARTdenovo/*/*/deconseq_Calb_mtDNA/contigs_min_500bp_filtered_renamed.fasta); do
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=repeat_masked/$Organism/"$Strain"/filtered_contigs
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
    qsub $ProgDir/rep_modeling.sh $Assembly $OutDir
    qsub $ProgDir/transposonPSI.sh $Assembly $OutDir
  done
```
The TransposonPSI masked bases were used to mask additional bases from the repeatmasker / repeatmodeller softmasked and hardmasked files.

```bash
for File in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_softmasked.fa); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
echo "Number of masked bases:"
cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' '
done
# The number of N's in hardmasked sequence are not counted as some may be present within the assembly and were therefore not repeatmasked.
for File in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_hardmasked.fa); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_hardmasked.fa/_contigs_hardmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -fi $File -bed $TPSI -fo $OutFile
done
```

```bash
Output of previous run:

repeat_masked/P.stipitis/589/filtered_contigs/589_contigs_softmasked_repeatmasker_TPSI_appended.fa
Number of masked bases:
520465
repeat_masked/P.stipitis/591/filtered_contigs/591_contigs_softmasked_repeatmasker_TPSI_appended.fa
Number of masked bases:
591492
repeat_masked/P.stipitis/594/filtered_contigs/594_contigs_softmasked_repeatmasker_TPSI_appended.fa
Number of masked bases:
539679
```

```bash
for File in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_unmasked.fa.TPSI.allHits.chains.bestPerLocus.gff3); do
Strain=$(echo $File| rev | cut -d '/' -f3 | rev)
Organism=$(echo $File | rev | cut -d '/' -f4 | rev)
# echo "$Organism - $Strain"
DDE_1=$(cat $File | cut -f9 | cut -f2 -d';' | sort| uniq -c | grep 'DDE_1' | sed "s/^\s*//g" | cut -f1 -d ' ')
Gypsy=$(cat $File | cut -f9 | cut -f2 -d';' | sort| uniq -c | grep 'gypsy' | sed "s/^\s*//g" | cut -f1 -d ' ')
HAT=$(cat $File | cut -f9 | cut -f2 -d';' | sort| uniq -c | grep 'hAT' | sed "s/^\s*//g" | cut -f1 -d ' ')
TY1_Copia=$(cat $File | cut -f9 | cut -f2 -d';' | sort| uniq -c | grep 'TY1_Copia' | sed "s/^\s*//g" | cut -f1 -d ' ')
Mariner=$(cat $File | cut -f9 | cut -f2 -d';' | sort| uniq -c | grep -w 'mariner' | sed "s/^\s*//g" | cut -f1 -d ' ')
Cacta=$(cat $File | cut -f9 | cut -f2 -d';' | sort| uniq -c | grep 'cacta' | sed "s/^\s*//g" | cut -f1 -d ' ')
LINE=$(cat $File | cut -f9 | cut -f2 -d';' | sort| uniq -c | grep 'LINE' | sed "s/^\s*//g" | cut -f1 -d ' ')
MuDR_A_B=$(cat $File | cut -f9 | cut -f2 -d';' | sort| uniq -c | grep 'MuDR_A_B' | sed "s/^\s*//g" | cut -f1 -d ' ')
HelitronORF=$(cat $File | cut -f9 | cut -f2 -d';' | sort| uniq -c | grep 'helitronORF' | sed "s/^\s*//g" | cut -f1 -d ' ')
Mariner_ant1=$(cat $File | cut -f9 | cut -f2 -d';' | sort| uniq -c | grep 'mariner_ant1' | sed "s/^\s*//g" | cut -f1 -d ' ')
ISC1316=$(cat $File | cut -f9 | cut -f2 -d';' | sort| uniq -c | grep 'ISC1316' | sed "s/^\s*//g" | cut -f1 -d ' ')
Crypton=$(cat $File | cut -f9 | cut -f2 -d';' | sort| uniq -c | grep 'Crypton' | sed "s/^\s*//g" | cut -f1 -d ' ')
printf "$Organism\t$Strain\t$DDE_1\t$Gypsy\t$HAT\t$TY1_Copia\t$Mariner\t$Cacta\t$LINE\t$MuDR_A_B\t$HelitronORF\t$Mariner_ant1\t$ISC1316\t$Crypton\n"
done
```
```
Output of previous run:

P.stipitis      589             13      1       99      1               13     14       8
P.stipitis      591             12      1       112     1       1       14     14       8
P.stipitis      594             13      1       108     1               14     14       8
```

Quast and BUSCO

```bash
for Assembly in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_unmasked.fa); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
OutDir=$(dirname $Assembly)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/saccharomycetales_odb9)
OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

```bash
for File in $(ls repeat_masked/*/*/filtered_contigs/run_*_contigs_unmasked/short_summary_*.txt); do
  Strain=$(echo $File| rev | cut -d '/' -f4 | rev)
  Organism=$(echo $File | rev | cut -d '/' -f5 | rev)
  Complete=$(cat $File | grep "(C)" | cut -f2)
  Fragmented=$(cat $File | grep "(F)" | cut -f2)
  Missing=$(cat $File | grep "(M)" | cut -f2)
  Total=$(cat $File | grep "Total" | cut -f2)
  echo -e "$Organism\t$Strain\t$Complete\t$Fragmented\t$Missing\t$Total"
  done
```

Output of BUSCO:

```
P.stipitis      589     1683    11      17      1711
P.stipitis      591     1678    14      19      1711
P.stipitis      594     1685    12      14      1711
```

```bash
for File in $(ls repeat_masked/*/*/filtered_contigs/report.tsv); do
  Strain=$(echo $File| rev | cut -d '/' -f3 | rev)
  Organism=$(echo $File | rev | cut -d '/' -f4 | rev)
  Contigs=$(cat $File | grep "contigs (>= 0 bp)" | cut -f2)
  Length=$(cat $File | grep "Total length (>= 0 bp)" | cut -f2)
  Largest=$(cat $File | grep "Largest contig" | cut -f2)
  N50=$(cat $File | grep "N50" | cut -f2)
  echo -e "$Organism\t$Strain\t$Contigs\t$Length\t$Largest\t$N50"
  done
```

Output of general data:
```
P.stipitis      589     11      15681772        2691698 1880751
P.stipitis      591     11      15602566        3982952 1903767
P.stipitis      594     8       15273029        3437457 1897630
```

# Promer alignment of Assemblies.

Alignment of our assemblies against the reference S. stipitis genome (AB580, Y-11545_v2)

```bash
Reference=$(ls assembly/misc_publications/P.stipitis/Y-11545_v2/pichia.allmasked)
for Query in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_hardmasked.fa); do
Strain=$(echo $Query | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Query | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
Prefix="$Strain"_vs_Y-11454_v2
OutDir=analysis/genome_alignment/mummer/$Organism/$Strain/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/MUMmer
qsub $ProgDir/sub_nucmer.sh $Reference $Query $Prefix $OutDir
done
```

Alignment of all our assemblies against the assembly done for 589

```bash
Reference=$(ls repeat_masked/P.stipitis/589/filtered_contigs/589_contigs_hardmasked.fa)
for Query in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_hardmasked.fa); do
Strain=$(echo $Query | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Query | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
Prefix="$Strain"_vs_589
OutDir=analysis/genome_alignment/mummer/$Organism/$Strain/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/MUMmer
qsub $ProgDir/sub_nucmer.sh $Reference $Query $Prefix $OutDir
done
```

Alignment of all our assemblies against the assembly done for 591

```bash
Reference=$(ls repeat_masked/P.stipitis/591/*/591_contigs_hardmasked.fa)
for Query in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_hardmasked.fa); do
Strain=$(echo $Query | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Query | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
Prefix="$Strain"_vs_591
OutDir=analysis/genome_alignment/mummer/$Organism/$Strain/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/MUMmer
qsub $ProgDir/sub_nucmer.sh $Reference $Query $Prefix $OutDir
done
```

Alignment of all our assemblies against the assembly done for 594

```bash
Reference=$(ls repeat_masked/P.stipitis/594/*/594_contigs_hardmasked.fa)
for Query in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_hardmasked.fa); do
Strain=$(echo $Query | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Query | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
Prefix="$Strain"_vs_594
OutDir=analysis/genome_alignment/mummer/$Organism/$Strain/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/MUMmer
qsub $ProgDir/sub_nucmer.sh $Reference $Query $Prefix $OutDir
done
```

All against reference genome.

```bash
for Reference in $(ls assembly/misc_publications/P.stipitis/Y-11545_v2/pichia.allmasked); do
for StrainPath in $(ls -d qc_dna/paired/P.stipitis/*); do
Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
F_Read=$(ls $StrainPath/F/*_trim.fq.gz)
R_Read=$(ls $StrainPath/R/*_trim.fq.gz)
echo $F_Read
echo $R_Read
Prefix="${Organism}_${Strain}"
OutDir=analysis/genome_alignment/bwa/$Organism/$Strain/vs_${Reference}
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/bwa
qsub $ProgDir/sub_bwa.sh $Prefix $Reference $F_Read $R_Read $OutDir
done
done
```
# Read coverage

1. Identify read coverage over each bp for Illumina alignments.

```bash
  for Bam in $(ls analysis/genome_alignment/bwa/*/*/vs_*/*_sorted.bam); do
    Strain=$(echo $Bam | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Bam | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=$(dirname $Bam)
    samtools depth -aa $Bam > $OutDir/${Organism}_${Strain}_depth.tsv
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/coverage_analysis
    $ProgDir/cov_by_window.py --cov $OutDir/${Organism}_${Strain}_depth.tsv > $OutDir/${Organism}_${Strain}_depth_10kb.tsv
    sed -i "s/$/\t$Strain/g" $OutDir/${Organism}_${Strain}_depth_10kb.tsv
  done
```

2. Coverage of Nanopore reads over Assembly.

For the alignment of ONT reads versus Nanopore assembly use the program minimap. The code of the program is the following:

```bash
#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 16
#$ -l virtual_free=1G
#$ -l h=blacklace02.blacklace|blacklace03.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace

#Align raw reads to an assembly.

# ---------------
# Step 1
# Collect inputs
# ---------------

Assembly=$(basename $1)
Reads=$(basename $2)
OutDir=$3


CurDir=$PWD
echo  "Running Bowtie with the following inputs:"
echo "Assembly - $Assembly"
echo "Reads - $Reads"
echo "OutDir - $OutDir"

# ---------------
# Step 2
# Copy data
# ---------------

WorkDir=$TMPDIR/minimap2
mkdir -p $WorkDir
cd $WorkDir
cp $CurDir/$1 $Assembly
cp $CurDir/$2 $Reads


# ---------------
# Step 3
# Align seq reads
# ---------------
# Prepare the assembly for alignment
# Align reads against the assembly
# Convert the SAM file to BAM in preparation for sorting.
# Sort the BAM file, in preparation for SNP calling:
# Index the bam file


minimap2 -ax map-ont $Assembly $Reads > ${Assembly}_aligned.sam
samtools view --threads 16 -bS ${Assembly}_aligned.sam -o ${Assembly}_aligned.bam
samtools sort --threads 16 -o ${Assembly}_aligned_sorted.bam ${Assembly}_aligned.bam

rm $Assembly
rm $Reads
rm ${Assembly}_aligned.sam
mkdir -p $CurDir/$OutDir
cp -r $WorkDir/* $CurDir/$OutDir/.
```

This program is run with the following code:

```bash
for Assembly in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_unmasked.fa); do
Reference=$(ls repeat_masked/*/589/filtered_contigs/*_contigs_unmasked.fa)
Strain=$(echo $Assembly | rev | cut -f3 -d '/'| rev)
Organism=$(echo $Reference | rev | cut -f3 -d '/' | rev)
Reads=$(ls qc_dna/minion/*/$Strain/*_trim.fastq.gz)
Prefix="${Organism}_${Strain}"
OutDir=analysis/genome_alignment/minimap/$Organism/vs_${Strain}
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment
qsub $ProgDir/minimap/sub_minimap2.sh $Reference $Reads $OutDir
done
```

# Indexing files from Minimap.

To visualise in IGV the bam files we need the index of that bam file. Minimap did not generate the index files of the alignments, so I will try to index them using samtools.

```bash
for File in $(ls analysis/genome_alignment/minimap/*/*/*_contigs_unmasked.fa_aligned_sorted.bam); do
Strain=$(echo $File | rev | cut -f2 -d '/'| rev)
File=analysis/genome_alignment/minimap/589/$Strain/*_contigs_unmasked.fa_aligned_sorted.bam
samtools index <$File> -o $File
done
```

Once the alignment is done we have to plot the coverage, generating first the .tsv files:

```bash
for Sam in $(ls analysis/genome_alignment/minimap/*/vs_*/*_aligned_sorted.bam); do
  Target=$(echo $Sam | rev | cut -f2 -d '/' | rev)
  Strain=$(echo $Sam | rev | cut -f3 -d '/' | rev)
  echo "$Strain-$Target"
  OutDir=$(dirname $Sam)
  samtools depth -aa $Sam > $OutDir/${Strain}_${Target}_depth.tsv
done

for Strain in 589 591 594; do
  for Cov in $(ls analysis/genome_alignment/minimap/*/vs_*/*_depth.tsv); do
    echo ${Cov} | cut -f4,5,6 -d '/' --output-delimiter " - "
    cat $Cov | cut -f3 | sort -n | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }'
  done
done > analysis/genome_alignment/minimap/read_coverage.txt
```
#Plot the coverage using R. Andy run this codes since the high coverage was observable in the contig 3 of 591:

```
cat analysis/genome_alignment/minimap/589/vs_589/589_vs_589_depth.tsv  | grep 'contig_3' > tmp.tsv

tmp <- read.delim("~/Downloads/tmp.tsv", header=FALSE)
library(ggplot2)
p <- ggplot(data=tmp, aes(x=tmp$V2, y=tmp$V3)) + geom_line() + labs(x = "Position (bp)", y = "Coverage") + geom_vline(xintercept = 1071000, colour = 'red', linetype = "dashed") + geom_vline(xintercept = 1153000, colour = 'red', linetype = "dashed")
outfile= paste("589", "vs_589", "contig3", "minion.jpg", sep = "_")
ggsave(outfile , plot = p, width = 20, height = 5, units = 'in', limitsize = TRUE)


cat analysis/genome_alignment/minimap/589/vs_591/589_vs_591_depth.tsv  | grep 'contig_3' > tmp2.tsv
cat analysis/genome_alignment/minimap/589/vs_594/589_vs_594_depth.tsv  | grep 'contig_3' > tmp3.tsv
```
```
tmp2 <- read.delim("~/Downloads/tmp2.tsv", header=FALSE)
library(ggplot2)
p2 <- ggplot(data=tmp, aes(x=tmp2$V2, y=tmp2$V3)) + geom_line() + labs(x = "Position (bp)", y = "Coverage") + geom_vline(xintercept = 1071000, colour = 'red', linetype = "dashed") + geom_vline(xintercept = 1153000, colour = 'red', linetype = "dashed")
outfile= paste("591", "vs_589", "contig3", "minion.jpg", sep = "_")
ggsave(outfile , plot = p2, width = 20, height = 5, units = 'in', limitsize = TRUE)
tmp3 <- read.delim("~/Downloads/tmp3.tsv", header=FALSE)
library(ggplot2)
p4 <- ggplot(data=tmp, aes(x=tmp3$V2, y=tmp3$V3)) + geom_line() + labs(x = "Position (bp)", y = "Coverage") + geom_vline(xintercept = 1071000, colour = 'red', linetype = "dashed") + geom_vline(xintercept = 1153000, colour = 'red', linetype = "dashed")
outfile= paste("594", "vs_589", "contig3", "minion.jpg", sep = "_")
ggsave(outfile , plot = p4, width = 20, height = 5, units = 'in', limitsize = TRUE)
```

I will use the same to try to plot the whole genome.

```
cat analysis/genome_alignment/minimap/589/vs_589/589_vs_589_depth.tsv > tmp_589.tsv
cat analysis/genome_alignment/minimap/589/vs_591/589_vs_591_depth.tsv > tmp_591.tsv
cat analysis/genome_alignment/minimap/589/vs_594/589_vs_594_depth.tsv > tmp_594.tsv

tmp_589 <- read.delim("~/Alignments/Coverage_ONT/tmp_589.tsv", header=FALSE)
tmp_591 <- read.delim("~/Alignments/Coverage_ONT/tmp_591.tsv", header=FALSE)
tmp_594 <- read.delim("~/Alignments/Coverage_ONT/tmp_594.tsv", header=FALSE)

library(ggplot2)

p <- ggplot(data=tmp_589, aes(x=tmp_589$V2, y=tmp_589$V3)) + geom_line() + labs(x = "Position (bp)", y = "Coverage") + geom_vline(xintercept = 1071000, colour = 'red', linetype = "dashed") + geom_vline(xintercept = 1153000, colour = 'red', linetype = "dashed")

outfile= paste("589", "vs_589", "minion.jpg", sep = "_")
ggsave(outfile , plot = p, width = 20, height = 5, units = 'in', limitsize = TRUE)
```
I obtained the result, but the circos plots look better.


# Extracting reads from the bam files.

In order to extract the reads from the bam files we will use the bam files.

Run at cd /projects/oldhome/groups/harrisonlab/project_files/Pichia/analysis/genome_alignment/minimap/589/vs_591. The -h option was added, because when sorting a message error appeared saying that the header was missing.

```
samtools view 589_contigs_unmasked.fa_aligned_sorted.bam "contig_3:1082100-1082200" -h > 591_cen5_region_200Bp.bam
```
Sorting of the bam file generated.

```
samtools sort -n 591_cen5_region_200Bp.bam > 591_cen5_region_200bp_sorted.bam
```
Now that we have the reads in a bam file, we can extract the fastq reads using bedtools.

```bash
bedtools bamtofastq -i 591_cen5_region_200bp_sorted.bam -fq 591_cen5_region_200bp_sorted.bam.fq
```

# Investigate GC content in the genome:

```bash
 for Assembly in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_unmasked.fa); do
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    OutDir=analysis/GC-content/$Organism/$Strain
    mkdir -p $OutDir
    PlotProg=$(ls /home/armita/prog/occultercut/OcculterCut_v1.1/plot.plt)
    gnuplot $PlotProg
    mv plot.eps ${Strain}_GC-plot.eps
  done
```

# Identify Telomere repeats

Telomeric repeats were identified in assemblies:

```bash
for Assembly in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_unmasked.fa); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=analysis/telomere/$Organism/$Strain
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/telomeres
$ProgDir/annotate_telomeres.py --fasta $Assembly --out $OutDir/telomere_hits
# Motif="TTAGGG"
# for Strand in '+' '-'; do
# ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen/codon
# python $ProgDir/identify_telomere_repeats.py $Assembly $Motif $Strand $OutDir/${Strain}_telomere_${Motif}_${Strand}.bed
# $ProgDir/how_many_repeats_regions.py $OutDir/${Strain}_telomere_${Motif}_${Strand}.bed
# done
done
cat $OutDir/telomere_hits.txt | sort -nr -k5 | less
```
<!--
# SNP calling codes.

The alignment by BWA has ben conducted already during the genome assembly. So we will start from the Pre SNP calling clean up step.

1. Pre SNP calling clean up

1.1  Rename input mapping files in each folder by prefixing with the strain ID.

I will use the alignments of the Parental strain (589) versus the evolved strains (591 and 594).

```bash
  for File in $(ls analysis/genome_alignment/bwa/*/589/vs_5*_unmasked/*sorted.bam); do
    Strain=$(echo $File | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $File | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=analysis/popgen/$Organism/$Strain
    CurDir=$PWD
    echo $OutDir
    mkdir -p $OutDir
    cd $OutDir
    cp -s $CurDir/$File "$Strain"_$(echo $File | rev | cut -f2 -d '/' | rev)_sorted.bam
    cd $CurDir
  done
```
1.2 Remove multimapping reads, discordant reads. PCR and optical duplicates, and add read group and sample name to each mapped read (preferably, the shortest ID possible)

Convention used: qsub $ProgDir/sub_pre_snp_calling.sh <SAMPLE_ID>

```bash
for Sam in $(ls analysis/popgen/*/*/*_sorted.bam); do
Strain=$(echo $Sam | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Sam | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/Pcac_popgen
qsub $ProgDir/sub_pre_snp_calling.sh $Sam $Strain
done
```

2. Run SNP calling

Prepare genome reference indexes required by GATK. Prepare for 589, 591 and 594.

```bash
for Reference in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
OutName=$(echo $Reference | sed 's/.fa/.dict/g')
OutDir=$(dirname $Reference)
mkdir -p $OutDir
ProgDir=/home/sobczm/bin/picard-tools-2.5.0
java -jar $ProgDir/picard.jar CreateSequenceDictionary R=$Reference O=$OutName
samtools faidx $Reference
done
```
Copy index file to same folder as BAM alignments

Move to the directory where the output of SNP calling should be placed. Then Start SNP calling with GATK. The submission script required need to be custom-prepared for each analysis, depending on what samples are being analysed. See inside the submission script below (GATK codes):

In order to run GATK I have used the next two set of commands, both of them start the run, but the outputs are empty. Folders are created normally but there is nothing inside.

```bash
Isolate=589
Reference=$(ls repeat_masked/*/589/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa| grep -w "${Isolate}")
CurDir=$PWD
OutDir=analysis/popgen/SNP_calling
mkdir -p $OutDir
cd $OutDir
ProgDir=/home/vegasa/git_repos/scripts/pichia/popgen
qsub $ProgDir/sub_SNP_calling_multithreaded2.sh $Reference $Isolate
cd $CurDir
```
```bash
for File in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
  Reference=$(ls repeat_masked/P.stipitis/589/filtered_contigs/589_contigs_softmasked_repeatmasker_TPSI_appended.fa)
  Isolate=$(echo $File | rev | cut -f3 -d '/' | rev)
  CurDir=/home/groups/harrisonlab/project_files/Pichia
  OutDir=analysis/popgen/SNP_calling
  ProgDir=/home/vegasa/git_repos/scripts/pichia/popgen
  qsub $ProgDir/sub_SNP_calling_multithreaded2.sh $Reference $Isolate
  cd $CurDir
done
```
The program for SNP calling is GATK:/home/armita/git_repos/emr_repos/scripts/phytophthora/Pcac_popgen/sub_SNP_calling_multithreaded2.sh
The codes have been run with the program as:

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
     -I $Project/analysis/popgen/P.stipitis/589/589_vs_594_unmasked_sorted_nomulti_proper_sorted_nodup_rg.bam \
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

``` -->



# Gene prediction

Gne prediction was performed using fungap on the NRI cluster.

## Data download:


RNAseq reads were downloaded for the reference genome

```bash
conda create -n sra-tools
conda activate sra-tools
conda install -c bioconda sra-tools
```

```bash
ProjDir=
cd $ProjDir
OutDir=raw_rna/paired/P.stipitis/Y-7124
mkdir -p $OutDir
fastq-dump --split-files --gzip --outdir $OutDir SRR8420582
```


## Functional annotation
A) Interproscan

Interproscan was used to give gene models functional annotations. Annotation was run using the commands below:

Note: This is a long-running script. As such, these commands were run using 'screen' to allow jobs to be submitted and monitored in the background. This allows the session to be disconnected and reconnected over time.

Screen ouptut detailing the progress of submission of interporscan jobs was redirected to a temporary output file named interproscan_submission.log .

From the new cluster:
```bash
# screen -a
# cd /projects/oldhome/groups/harrisonlab/project_files/Pichia
#
# SplitfileDir=/projects/oldhome/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
# ProgDir=/projects/oldhome/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
# CurPath=$PWD
# for Proteome in $(ls gene_pred/fungap/P.stipitis/589/fungap_out/fungap_out_prot.faa); do
# Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
# Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
# SplitDir=gene_pred/fungap_split/$Organism/$Strain
# mkdir -p $SplitDir
# BaseName="$Organism""_$Strain"_fungap
# $SplitfileDir/splitfile_500.py --inp_fasta $Proteome --out_dir $SplitDir --out_base $BaseName
# done
#
# for Fasta in $(ls gene_pred/fungap/P.stipitis/589/fungap_out/fungap_out_prot.faa); do
# 	Jobs=$(squeue -n squeue -n slurm_interproscan.sh -t pd,s,cf,ca,f,to,pr,bf,nf,se | wc -l)
# 	while [ $Jobs -gt 3 ]; do
# 	sleep 10
# 	printf "."
# 	Jobs=$(squeue -n squeue -n slurm_interproscan.sh -t pd,s,cf,ca,f,to,pr,bf,nf,se | wc -l)
# 	done
# 	printf "\n"
# 	echo $File
# 	OutDir=$(dirname $Fasta)/raw
# 	ProgDir=/projects/oldhome/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
# 	sbatch $ProgDir/slurm_interproscan.sh $Fasta $OutDir
# done

Fasta=$(ls gene_pred/fungap/P.stipitis/589/fungap_out/fungap_out_prot.faa)
OutDir=$(dirname $Fasta | sed 's/fungap/interproscan/g')
ProgDir=/projects/oldhome/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
sbatch $ProgDir/slurm_interproscan.sh $Fasta $OutDir
```

# Assembly stats were collected for comparison to the reference assembly:


7. Quast and busco were run to assess the effects of pilon on assembly quality:

```bash
conda activate fungap
for Assembly in $(ls assembly/misc_publications/P.stipitis/*/* | grep -e 'pichia.AssembledScaffolds.fasta' -e 'GCF_000209165.1_ASM20916v1_genomic.fna' | grep -v '.fai'); do
  Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
  echo "$Organism - $Strain"
  OutDir=$(dirname $Assembly)
  ProgDir=/projects/oldhome/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  sbatch $ProgDir/slurm_quast.sh $Assembly $OutDir

  BuscoDB=$(ls -d /projects/oldhome/groups/harrisonlab/dbBusco/ascomycota_odb9)
  OutDir=$(dirname $Assembly)
  ProgDir=/projects/oldhome/armita/git_repos/emr_repos/tools/gene_prediction/busco
  # sbatch $ProgDir/slurm_busco_v3.sh $Assembly $BuscoDB $OutDir
done
```

```bash
printf "Filename\tComplete\tDuplicated\tFragmented\tMissing\tTotal\n"
for File in $(ls assembly/misc_publications/P.stipitis/*/short_summary_*.txt); do  
FileName=$(basename $File)
Complete=$(cat $File | grep "(C)" | cut -f2)
Duplicated=$(cat $File | grep "(D)" | cut -f2)
Fragmented=$(cat $File | grep "(F)" | cut -f2)
Missing=$(cat $File | grep "(M)" | cut -f2)
Total=$(cat $File | grep "Total" | cut -f2)
printf "$FileName\t$Complete\t$Duplicated\t$Fragmented\t$Missing\t$Total\n"
done
```
